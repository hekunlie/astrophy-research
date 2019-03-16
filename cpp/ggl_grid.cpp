#include<FQlib.h>
#include<mpi.h>

int main(int argc, char *argv[])
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	int i, j, k;
	int area_num = 4;
	int area_id;
	double st1, st2, st3, st4, st5, st6, st7,st8, sts, ste;

	int z_id = 0, ra_id = 1, dec_id = 2, mg1_id = 3, mg2_id = 4, mn_id = 5, mu_id = 6, mv_id = 7;
	int nib_id = 8, bs_id = 9, be_id = 10, bdy_id = 11, bdx_id = 12;
	// for the initial data
	int num_ini;
	int data_col = 8;	
	double *data_ini[13];
	// for others
	double *data_temp[13];
	int *mask, count;

	double ra, ra_min, ra_max;
	double dec, dec_min, dec_max;
	double *ra_bin, *dec_bin;
	double margin;
	int grid_num, grid_ny, grid_nx;
	int row, col, bin_label;
	int my_gal_s, my_gal_e;

	int *num_in_block, *block_start, *block_end;
	int *block_mask, *count_in_block;
	double *block_boundy, *block_boundx;

	// redshift range for foreground galaxy
	double low_z, high_z, block_scale;
	low_z = atof(argv[1]);
	high_z = atof(argv[2]);
	block_scale = atof(argv[3]);
	margin = 0.1*block_scale;

	char h5f_path_1[150], h5f_path_2[150], log_inform[200];
	char data_path[150], log_path[150];

	char  attrs_name[50];
	int shape[2];
	char set_name[50],*names[13];
	for (i = 0; i < 13; i++)
	{
		names[i] = new char[40];
	}
	sprintf(names[z_id], "Z");
	sprintf(names[ra_id], "RA");
	sprintf(names[dec_id], "DEC");
	sprintf(names[mg1_id], "G1");
	sprintf(names[mg2_id], "G2");
	sprintf(names[mn_id], "N");
	sprintf(names[mu_id], "U");
	sprintf(names[mv_id], "V");
	sprintf(names[nib_id], "num_in_block");
	sprintf(names[bs_id], "block_start");
	sprintf(names[be_id], "block_end");
	sprintf(names[bdy_id], "block_boundy");
	sprintf(names[bdx_id], "block_boundx");

	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/");
	sprintf(h5f_path_1, "%scata_result_ext_cut.hdf5", data_path);
	sprintf(h5f_path_2, "%scata_result_ext_grid.hdf5", data_path);
	sprintf(log_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/log/grid_log_%d.dat", rank);

	if (0 == rank)
	{
		sprintf(set_name, "/foreground/");
		create_h5_group(h5f_path_2, set_name, TRUE);
		for (i = 1; i < area_num + 1; i++)
		{
			sprintf(set_name, "/foreground/w_%d", i);
			create_h5_group(h5f_path_2, set_name, FALSE);
			sprintf(set_name, "/background/w_%d", i);
			create_h5_group(h5f_path_2, set_name, FALSE);
		}
	}
	//sprintf(log_inform, "Foreground: %d galaxies in W_%d", area_id, count + 1);
	//std::cout << log_inform << std::endl;
	
	for (area_id = 1; area_id < area_num+1; area_id++)
	{
		st1 = clock();

		sprintf(attrs_name, "shape");
		sprintf(set_name, "/w_%d/Z", area_id);
		read_h5_attrs(h5f_path_1, set_name, attrs_name, shape, "d");
		num_ini = shape[0];
		mask = new int[num_ini] {};
		for (i = 0; i < data_col; i++)
		{
			data_ini[i] = new double[num_ini];
			sprintf(set_name, "/w_%d/%s", area_id, names[i]);
			read_h5(h5f_path_1, set_name, data_ini[i]);
		}	
		
		// the foreground galaxy
		if (0 == rank)
		{	
			count = 0;
			for (i = 0; i < num_ini; i++)
			{
				// find the foreground galaxy
				if (low_z <= data_ini[z_id][i] and data_ini[z_id][i] <= high_z)
				{
					mask[i] = 1;
					count += 1;
				}
			}
			shape[0] = count;
			shape[1] = 1;
			
			for (i = 0; i < data_col; i++)
			{	
				// alloc the array for foreground galaxy
				data_temp[i] = new double[count];
			}
			count = 0;
			for (i = 0; i < num_ini; i++)
			{				
				if (1 == mask[i])
				{
					for (j = 0; j < data_col; j++)
					{
						data_temp[j][count] = data_ini[j][i];
					}
					count += 1;
				}				
			}
			
			sprintf(attrs_name, "shape");
			for (j = 0; j < data_col; j++)
			{
				// write to file
				sprintf(set_name, "/foreground/w_%d/%s", area_id, names[j]);
				write_h5(h5f_path_2, set_name, data_temp[j], count, 1, FALSE);
				write_h5_attrs(h5f_path_2, set_name, attrs_name, shape, 2, "d");
				// free the memory
				delete[] data_temp[j];
			}
		}
		st2 = clock();
		sprintf(log_inform, "Rank %d w_%d:  %d foreground galaxies in (%.2f sec)", rank, area_id, count, (st2 - st1) / CLOCKS_PER_SEC);
		write_log(log_path, log_inform);
		if (0 == rank)
		{			
			std::cout << log_inform << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);


		// build the grid
		ra_min = *std::min_element(data_ini[ra_id], data_ini[ra_id] + num_ini);
		ra_max = *std::max_element(data_ini[ra_id], data_ini[ra_id] + num_ini);
		ra_min = ra_min - margin;

		dec_min = *std::min_element(data_ini[dec_id], data_ini[dec_id] + num_ini);
		dec_max = *std::max_element(data_ini[dec_id], data_ini[dec_id] + num_ini);
		dec_min = dec_min - margin;

		grid_ny = (int) (dec_max - dec_min) / block_scale + 1;
		grid_nx = (int) (ra_max - ra_min) / block_scale + 1;
		grid_num = grid_nx * grid_ny;

		ra_bin = new double[grid_nx + 1]{};
		dec_bin = new double[grid_ny + 1]{};

		for (i = 0; i < grid_nx+1; i++)
		{
			ra_bin[i] = ra_min + i * block_scale;
		}
		for (i = 0; i < grid_ny + 1; i++)
		{
			dec_bin[i] = dec_min + i * block_scale;
		}

		// write the boundary of blocks to file
		if (0 == rank)
		{
			block_boundx = new double[grid_num * 4];
			block_boundy = new double[grid_num * 4];

			for (i = 0; i < grid_num; i++)
			{
				row = i / grid_nx;
				col = i % grid_nx;
				block_boundy[i * 4] = dec_bin[row];
				block_boundy[i * 4 + 1] = dec_bin[row];
				block_boundy[i * 4 + 2] = dec_bin[row + 1];
				block_boundy[i * 4 + 3] = dec_bin[row + 1];

				block_boundx[i * 4] = ra_bin[col];
				block_boundx[i * 4 + 1] = ra_bin[col + 1];
				block_boundx[i * 4 + 2] = ra_bin[col];
				block_boundx[i * 4 + 3] = ra_bin[col + 1];
			}

			shape[0] = grid_num;
			shape[1] = 4;
			sprintf(attrs_name, "shape");

			sprintf(set_name, "/background/w_%d/%s", area_id, names[bdx_id]);
			write_h5(h5f_path_2, set_name, block_boundx, grid_num, 4, FALSE);
			write_h5_attrs(h5f_path_2, set_name, attrs_name, shape, 2, "d");

			sprintf(set_name, "/background/w_%d/%s", area_id, names[bdy_id]);
			write_h5(h5f_path_2, set_name, block_boundy, grid_num, 4, FALSE);
			write_h5_attrs(h5f_path_2, set_name, attrs_name, shape, 2, "d");

			delete[] block_boundx;
			delete[] block_boundy;

		}		
		st3 = clock();
		sprintf(log_inform, "Rank %d w_%d:  Built RA, Dec and boundaries bin (%.2f sec)", rank, area_id, count, (st3 - st2) / CLOCKS_PER_SEC);
		write_log(log_path, log_inform);
		if (0 == rank)
		{			
			std::cout << log_inform << std::endl;
		}

		MPI_Win win;
		MPI_Aint size = num_ini * sizeof(int);
		if (0 == rank)
		{
			MPI_Win_allocate_shared(size, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &block_mask, &win);
		}
		else
		{
			int disp_unit;
			MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &block_mask, &win);
			MPI_Win_shared_query(win, 0, &size, &disp_unit, &block_mask);
		}	
		MPI_Barrier(MPI_COMM_WORLD);
		if (0 == rank)
		{
			initialize_arr(block_mask, num_ini, -1);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// distribute the tasks
		if (numprocs - 1 == rank)
		{
			my_gal_s = num_ini / numprocs * rank;
			my_gal_e = num_ini / numprocs * (rank+1) + num_ini % numprocs;
		}
		else
		{
			my_gal_s = num_ini / numprocs * rank;
			my_gal_e = num_ini / numprocs * (rank + 1);
		}
		sprintf(log_inform, "Rank %d w_%d: my job: %d ~ %d", rank, area_id, my_gal_s, my_gal_e);
		if (0 == rank)
		{
			std::cout << log_inform << std::endl;
		}
		write_log(log_path, log_inform);

		for (i = my_gal_s; i < my_gal_e; i++)
		{
			histogram2d_s(data_ini[dec_id][i], data_ini[ra_id][i], dec_bin, ra_bin, grid_ny, grid_nx, bin_label);
			block_mask[i] = bin_label;
		}
		MPI_Barrier(MPI_COMM_WORLD);

		st4 = clock();
		sprintf(log_inform, "Rank %d w_%d: Find all the blocks (%.2f sec)", rank, area_id, (st4-st3)/CLOCKS_PER_SEC);
		write_log(log_path, log_inform);
		if (0 == rank)
		{
			std::cout << log_inform << std::endl;
		}
		
		if (0 == rank)
		{
			for (i = 0; i < data_col; i++)
			{
				// alloc the array for background galaxy
				data_temp[i] = new double[num_ini];
				initialize_arr(data_temp[i], num_ini, 0);
			}

			num_in_block = new int[grid_num] {};
			block_start = new int[grid_num] {};
			block_end = new int[grid_num] {};
			count_in_block = new int[grid_num] {};

			initialize_arr(num_in_block, grid_num, 0);
			initialize_arr(block_start, grid_num, 0);
			initialize_arr(block_end, grid_num, 0);
			initialize_arr(count_in_block, grid_num, 0);

			// count
			for (i = 0; i < num_ini; i++)
			{
				for (j = 0; j < grid_num; j++)
				{
					if (block_mask[i] == j)
					{
						num_in_block[j] += 1;
					}
				}
			}
			st5 = clock();
			sprintf(log_inform, "Rank %d w_%d: Count the number in each block (%.2f sec)", rank, area_id, (st5 - st4) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			std::cout << log_inform << std::endl;

			// the start and end label of each block
			for (i = 0; i < grid_num; i++)
			{
				for (j = 0; j < i; j++)
				{
					block_start[i] += num_in_block[j];			
				}
				block_end[i] = block_start[i] + num_in_block[i];
			}
			st6 = clock();
			sprintf(log_inform, "Rank %d w_%d: Calculate the start- and end-point for each block (%.2f sec)", rank, area_id, (st6 - st5) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			std::cout << log_inform << std::endl;

			for (i = 0; i < num_ini; i++)
			{
				for (j = 0; j < grid_num; j++)
				{
					if (block_mask[i] == j)
					{
						row = count_in_block[j] + block_start[j];
						for (k = 0; k < data_col; k++)
						{
							data_temp[k][row] = data_ini[k][i];
						}
						count_in_block[j] += 1;
					}
				}
			}
			st7 = clock();
			sprintf(log_inform, "Rank %d w_%d: Reorganize the data (%.2f sec)", rank, area_id, (st7 - st6) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			std::cout << log_inform << std::endl;

			sprintf(attrs_name, "shape");
			shape[0] = num_ini;
			shape[1] = 1;
			for (j = 0; j < data_col; j++)
			{
				// write to file
				sprintf(set_name, "/background/w_%d/%s", area_id, names[j]);
				write_h5(h5f_path_2, set_name, data_temp[j], num_ini, 1, FALSE);
				write_h5_attrs(h5f_path_2, set_name, attrs_name, shape, 2, "d");
				// free the memory
				delete[] data_temp[i];
			}
			st8 = clock();
			sprintf(log_inform, "Rank %d w_%d: Write to file (%.2f sec)", rank, area_id, (st8 - st7) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			std::cout << log_inform << std::endl;

			delete[] num_in_block;
			delete[] block_start;
			delete[] block_end;
			delete[] count_in_block;
		}

		for (i = 0; i < data_col; i++)
		{
			delete[] data_ini[i];
		}
		delete[] mask;
		delete[] ra_bin;
		delete[] dec_bin;

		MPI_Win_free(&win);
		MPI_Finalize();
	}

	return 0;
}