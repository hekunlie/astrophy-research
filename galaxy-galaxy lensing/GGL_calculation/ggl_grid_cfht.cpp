#include<FQlib.h>
#include<mpi.h>

#define backgal_data_col 20
#define grid_data_col 5
#define max_data_col 30
#define mg_bin_num 8

int main(int argc, char *argv[])
{
	// argv[1]: block scale in unit of degree

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	int i, j, k;
	int area_num = 4, radius_num = 13;
	int area_id;
	double st1, st2, st3, st4, st5, st6, st7, st8, sts, ste;

	int nib_id = 0, bs_id = 1, be_id = 2, bdy_id = 3, bdx_id = 4;
	int z_id = 5, dist_id = 6, ra_id = 7, dec_id = 8, cos_dec_id = 9;
	int e1_id = 10, e2_id = 11, weight_id = 12, m_id = 13, c_id = 14;
	int starflag_id=15, zmin_lb = 16, zmax_lb = 17, odds_lb = 18, mag_lb = 19;

	// for the initial data
	int num_ini, data_len;
	double *data_ini[max_data_col];
	// for others
	double *data_temp[max_data_col], *data_buffer;
	int *mask, count;

	double ra, ra_min, ra_max;
	double dec, dec_min, dec_max;
	double *ra_bin, *dec_bin;
	double margin;
	int grid_num, grid_ny, grid_nx;
	int row, col, bin_label;
	int my_gal_s, my_gal_e, my_grid_s, my_grid_e;

	int *num_in_block, *block_start, *block_end;
	int *block_buffer, *count_in_block;
	double *block_boundy, *block_boundx;

	// redshift range for foreground galaxy
	double block_scale; // degree

	block_scale = atof(argv[1]);
	margin = 0.1*block_scale;

	char h5f_path_src[150], h5f_path_dst[150], log_inform[200];
	char data_path[150], log_path[150];

	char  attrs_name[50];
	int shape[2];
	double scale[1];
	double red_bin[2];
	char set_name[50], *names[max_data_col];
	for (i = 0; i < backgal_data_col; i++)
	{
		names[i] = new char[40];
	}

	sprintf(names[nib_id], "num_in_block");
	sprintf(names[bs_id], "block_start");
	sprintf(names[be_id], "block_end");
	sprintf(names[bdy_id], "block_boundy");
	sprintf(names[bdx_id], "block_boundx");

	sprintf(names[z_id], "Z");
	sprintf(names[dist_id], "DISTANCE");
	sprintf(names[ra_id], "RA");
	sprintf(names[dec_id], "DEC");
	sprintf(names[cos_dec_id], "COS_DEC");

	sprintf(names[e1_id], "E1");
	sprintf(names[e2_id], "E2");
	sprintf(names[weight_id], "WEIGHT");
	sprintf(names[m_id], "M");
	sprintf(names[c_id], "C");

	sprintf(names[starflag_id], "STARGLAG");
	sprintf(names[zmin_lb], "Z_MIN");
	sprintf(names[zmax_lb], "Z_MAX");
	sprintf(names[odds_lb], "ODDS");
	sprintf(names[mag_lb], "MAG");


	sprintf(data_path, "/mnt/perc/hklee/CFHT/gg_lensing/data/");
	sprintf(h5f_path_src, "%scfht_cata_cut.hdf5", data_path);
	sprintf(h5f_path_dst, "%scfht_cata_grid.hdf5", data_path);
	sprintf(log_path, "/mnt/perc/hklee/CFHT/gg_lensing/log/grid_log_%d.dat", rank);

	if (0 == rank)
	{
		sprintf(set_name, "/block_scale");
		scale[0] = block_scale;
		write_h5(h5f_path_dst, set_name, scale, 1, 1, TRUE);

		for (i = 1; i < area_num + 1; i++)
		{
			sprintf(set_name, "/w_%d", i);
			create_h5_group(h5f_path_dst, set_name, FALSE);
		}

		double *radius_bin = new double[radius_num + 1];
		sprintf(set_name, "/radius_bin");
		// the radius bin (log space)
		log_bin(0.04, 15, radius_num + 1, radius_bin);
		write_h5(h5f_path_dst, set_name, radius_bin, radius_num + 1, 1, FALSE);
	}

	for (area_id = 1; area_id < area_num + 1; area_id++)
	{
		st1 = clock();

		// read the former 'backgal_data_col' columns , Z, DISTANCE,... ~ V
		sprintf(set_name, "/w_%d/Z", area_id);
		read_h5_datasize(h5f_path_src, set_name, num_ini);
		mask = new int[num_ini] {};
		for (i = grid_data_col; i < backgal_data_col; i++)
		{
			data_ini[i] = new double[num_ini];
			sprintf(set_name, "/w_%d/%s", area_id, names[i]);
			read_h5(h5f_path_src, set_name, data_ini[i]);
		}
		
		st2 = clock();
		MPI_Barrier(MPI_COMM_WORLD);

		// build the grid
		ra_min = *std::min_element(data_ini[ra_id], data_ini[ra_id] + num_ini);
		ra_max = *std::max_element(data_ini[ra_id], data_ini[ra_id] + num_ini);
		ra_min = ra_min - 0.1*margin;

		dec_min = *std::min_element(data_ini[dec_id], data_ini[dec_id] + num_ini);
		dec_max = *std::max_element(data_ini[dec_id], data_ini[dec_id] + num_ini);
		dec_min = dec_min - 0.1*margin;

		grid_ny = (int)((dec_max - dec_min) / block_scale + 2);
		grid_nx = (int)((ra_max - ra_min) / block_scale + 2);
		grid_num = grid_nx * grid_ny;
		if (0 == rank)
		{
			std::cout << "Rank 0 w_" << area_id << ": ";
			std::cout << "DEC: " << dec_min << " " << dec_max << " (" << dec_max - dec_min << ") " << grid_ny << std::endl;
			std::cout << "Rank 0 w_" << area_id << ": ";
			std::cout << "RA: " << ra_min << " " << ra_max << " (" << ra_max - ra_min << ") " << grid_nx << std::endl;

			shape[0] = grid_ny;
			shape[1] = grid_nx;

			sprintf(set_name, "/w_%d/grid_shape", area_id);
			write_h5(h5f_path_dst, set_name, shape, 2, 1, FALSE);
		}

		ra_bin = new double[grid_nx + 1]{};
		dec_bin = new double[grid_ny + 1]{};

		for (i = 0; i < grid_nx + 1; i++)
		{
			ra_bin[i] = ra_min + i * block_scale;
		}
		for (i = 0; i < grid_ny + 1; i++)
		{
			dec_bin[i] = dec_min + i * block_scale;
		}
		if (ra_bin[grid_nx] <= ra_max or dec_bin[grid_ny] <= dec_max)
		{
			std::cout << "The RA & Dec bins don't cover all the region!!!";
			exit(0);
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

			sprintf(set_name, "/w_%d/%s", area_id, names[bdx_id]);
			write_h5(h5f_path_dst, set_name, block_boundx, grid_num, 4, FALSE);

			sprintf(set_name, "/w_%d/%s", area_id, names[bdy_id]);
			write_h5(h5f_path_dst, set_name, block_boundy, grid_num, 4, FALSE);

			shape[0] = grid_nx + 1;
			shape[1] = 1;
			sprintf(set_name, "/w_%d/RA_bin", area_id);
			write_h5(h5f_path_dst, set_name, ra_bin, shape[0], shape[1], FALSE);

			shape[0] = grid_ny + 1;
			shape[1] = 1;
			sprintf(set_name, "/w_%d/DEC_bin", area_id);
			write_h5(h5f_path_dst, set_name, dec_bin, shape[0], shape[1], FALSE);

			delete[] block_boundx;
			delete[] block_boundy;

		}
		st3 = clock();
		sprintf(log_inform, "Rank %d w_%d: Built RA, Dec and boundaries bin (%.2f sec)", rank, area_id, (st3 - st2) / CLOCKS_PER_SEC);
		write_log(log_path, log_inform);
		if (0 == rank)
		{
			std::cout << log_inform << std::endl;
		}

		// shared memory buffer  
		// "num_ini" elements (= lenof(data)) are the block labels of galaxies,
		// "grid_num" elements are the number counts for each block
		MPI_Win win, win_d;
		MPI_Aint size = (num_ini + grid_num) * sizeof(int), size_d = num_ini * sizeof(double);
		if (0 == rank)
		{
			MPI_Win_allocate_shared(size, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &block_buffer, &win);
			MPI_Win_allocate_shared(size_d, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &data_buffer, &win_d);
		}
		else
		{
			int disp_unit;
			MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &block_buffer, &win);
			MPI_Win_shared_query(win, 0, &size, &disp_unit, &block_buffer);

			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &data_buffer, &win_d);
			MPI_Win_shared_query(win_d, 0, &size_d, &disp_unit, &data_buffer);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (0 == rank)
		{
			initialize_arr(block_buffer, num_ini + grid_num, 0);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// distribute the tasks
		if (numprocs - 1 == rank)
		{
			my_gal_s = num_ini / numprocs * rank;
			my_gal_e = num_ini / numprocs * (rank + 1) + num_ini % numprocs;

			my_grid_s = grid_num / numprocs * rank;
			my_grid_e = grid_num / numprocs * (rank + 1) + grid_num % numprocs;
		}
		else
		{
			my_gal_s = num_ini / numprocs * rank;
			my_gal_e = num_ini / numprocs * (rank + 1);

			my_grid_s = grid_num / numprocs * rank;
			my_grid_e = grid_num / numprocs * (rank + 1);
		}
		sprintf(log_inform, "Rank %d w_%d: my gal: %d ~ %d (%d), my grid: %d ~ %d (%d x %d = %d)", rank, area_id, my_gal_s, my_gal_e, num_ini, my_grid_s, my_grid_e, grid_ny, grid_nx, grid_num);
		write_log(log_path, log_inform);
		if (0 == rank)
		{
			std::cout << log_inform << std::endl;
		}

		// find the block label for each galaxy
		for (i = my_gal_s; i < my_gal_e; i++)
		{
			histogram2d_s(data_ini[dec_id][i], data_ini[ra_id][i], dec_bin, ra_bin, grid_ny, grid_nx, bin_label);
			block_buffer[i] = bin_label;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		// count the number in each block
		for (i = 0; i < num_ini; i++)
		{
			for (j = my_grid_s; j < my_grid_e; j++)
			{
				if (block_buffer[i] == j)
				{
					block_buffer[num_ini + j] += 1;
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		st4 = clock();
		sprintf(log_inform, "Rank %d w_%d: Find all the blocks and count the number for all blocks (%.2f sec)", rank, area_id, (st4 - st3) / CLOCKS_PER_SEC);
		write_log(log_path, log_inform);
		if (0 == rank)
		{
			std::cout << log_inform << std::endl;
		}

		num_in_block = new int[grid_num] {};
		block_start = new int[grid_num] {};
		block_end = new int[grid_num] {};
		count_in_block = new int[grid_num] {};

		initialize_arr(num_in_block, grid_num, 0);
		initialize_arr(block_start, grid_num, 0);
		initialize_arr(block_end, grid_num, 0);
		initialize_arr(count_in_block, grid_num, 0);

		// the start and end label of each block
		for (i = 0; i < grid_num; i++)
		{
			num_in_block[i] = block_buffer[num_ini + i];
			for (j = 0; j < i; j++)
			{
				block_start[i] += num_in_block[j];
			}
			block_end[i] = block_start[i] + num_in_block[i];
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (k = grid_data_col; k < backgal_data_col; k++)
		{
			st5 = clock();
			initialize_arr(count_in_block, grid_num, 0);
			for (i = 0; i < num_ini; i++)
			{
				for (j = my_grid_s; j < my_grid_e; j++)
				{
					if (block_buffer[i] == j)
					{
						row = count_in_block[j] + block_start[j];
						data_buffer[row] = data_ini[k][i];
						count_in_block[j] += 1;
					}
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			st6 = clock();
			sprintf(log_inform, "Rank %d w_%d: Reorganize the data %s (%.2f sec)", rank, area_id, names[k], (st6 - st5) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);

			// write to file
			if (0 == rank)
			{
				std::cout << log_inform << std::endl;
				sprintf(set_name, "/w_%d/%s", area_id, names[k]);
				write_h5(h5f_path_dst, set_name, data_buffer, num_ini, 1, FALSE);

				if (k == grid_data_col)
				{
					// block_start
					sprintf(set_name, "/w_%d/%s", area_id, names[bs_id]);
					write_h5(h5f_path_dst, set_name, block_start, grid_num, 1, FALSE);

					// block_end
					sprintf(set_name, "/w_%d/%s", area_id, names[be_id]);
					write_h5(h5f_path_dst, set_name, block_end, grid_num, 1, FALSE);

					// num_in_block
					sprintf(set_name, "/w_%d/%s", area_id, names[nib_id]);
					write_h5(h5f_path_dst, set_name, num_in_block, grid_num, 1, FALSE);
				}
			}
			st7 = clock();
			sprintf(log_inform, "Rank %d w_%d: Write to file (%.2f sec)", rank, area_id, (st7 - st6) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			if (0 == rank)
			{
				std::cout << log_inform << std::endl;
			}

			MPI_Barrier(MPI_COMM_WORLD);
		}
		if (0 == rank)
		{
			std::cout << std::endl;
		}
		for (i = grid_data_col; i < backgal_data_col; i++)
		{
			delete[] data_ini[i];
		}

		delete[] mask;
		delete[] ra_bin;
		delete[] dec_bin;

		delete[] num_in_block;
		delete[] block_start;
		delete[] block_end;
		delete[] count_in_block;
		MPI_Win_free(&win);
		MPI_Win_free(&win_d);
	}

	MPI_Finalize();
	return 0;
}