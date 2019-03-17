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

	pts_info gal_info;

	int i, j, k;
	int area_id, area_num;
	area_num = 4;
	
	int foregal_num, gal_id;
	double *foregal_data[3];
	double z_f, ra_f, dec_f;


	int data_num;
	int backgal_num;
	int z_id = 0, ra_id = 1, dec_id = 2;
	int mg1_id = 3, mg2_id = 4, mn_id = 5, mu_id = 6, mv_id = 7;
	int nib_id = 8, bs_id = 9, be_id = 10, bdy_id = 11, bdx_id = 12;
	int ra_bin_id = 13, dec_bin_id = 14;
	double *backgal_data[13];
	double z_b, ra_b, dec_b;
	double *ra_bin, *dec_bin;
	int ra_bin_num, dec_bin_num;
	int my_gal_s, my_gal_e;


	int grid_num, grid_ny, grid_nx;
	int row, col, bin_label;
	double block_scale[1];
	int *block_mask;


	int radius_num, rad_id;
	double radius_s, radius_e;
	double *radius_bin;


	double *redshifts, *distances;
	int red_num;
	int shape[2];
	char data_path[150], log_path[150], h5f_path[150];
	char set_name[50], attrs_name[50];

	char *names[15];
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

	sprintf(names[ra_bin_id], "RA_bin");
	sprintf(names[dec_bin_id], "DEC_bin");


	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/");	
	sprintf(log_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/log/ggl_log_%d.dat", rank);

	// the Z-comoving distance have be calculated from z =0 ~ z =10
	sprintf(h5f_path, "%sredshift.hdf5", data_path);
	sprintf(set_name, "/redshift");
	sprintf(attrs_name, "shape");
	read_h5_attrs(h5f_path, set_name, attrs_name, shape, "d");
	red_num = shape[0];
	redshifts = new double[red_num] {};
	distances = new double[red_num] {};

	read_h5(h5f_path, set_name, redshifts);
	sprintf(set_name, "/distance");
	read_h5(h5f_path, set_name, distances);

	// read the search radius
	sprintf(h5f_path, "%scata_result_ext_grid.hdf5", data_path);
	sprintf(set_name, "/radius_bin");
	sprintf(attrs_name, "shape");
	read_h5_attrs(h5f_path, set_name, attrs_name, shape, "d");
	radius_num = shape[0];
	radius_bin = new double[radius_num] {};
	read_h5(h5f_path, set_name, radius_bin);
	radius_num = radius_num - 1;


	for (area_id = 1; area_id < area_num + 1; area_id++)
	{
		// read foreground information
		sprintf(attrs_name, "shape");
		// Z, RA, DEC
		for (i = 0; i < 3; i++)
		{
			sprintf(set_name, "/foreground/w_%d/%s", area_id, names[i]);
			read_h5_attrs(h5f_path, set_name, attrs_name, shape, "d");
			foregal_num = shape[0];

			foregal_data[i] = new double[foregal_num] {};			
			read_h5(h5f_path, set_name, foregal_data[i]);
		}


		// read background information
		sprintf(set_name, "/background/w_%d", area_id);
		sprintf(attrs_name, "grid_shape");
		read_h5_attrs(h5f_path, set_name, attrs_name, shape, "g");
		grid_ny = shape[0];
		grid_nx = shape[1];
		grid_num = grid_ny * grid_nx;
		sprintf(attrs_name, "block_scale");
		read_h5_attrs(h5f_path, set_name, attrs_name, block_scale, "g");
		
		block_mask = new int[grid_num];

		// Z, RA, DEC,  G1, G2, N, U, V,  num_in_block,  block_start, block_end, 
		// block_boundx, block_boundy
		for (i = 0; i < 15; i++)
		{			
			sprintf(set_name, "/background/w_%d/%s", area_id, names[i]);		
			sprintf(attrs_name, "shape");
			read_h5_attrs(h5f_path, set_name, attrs_name, shape, "d");
			if (0 == i)
			{
				backgal_num = shape[0];
			}
			if (13 == i)
			{
				ra_bin_num = shape[0] - 1;
			}
			if (14 == i)
			{
				dec_bin_num = shape[0] - 1;
			}
			backgal_data[i] = new double[shape[0] * shape[1]] {};
			read_h5(h5f_path, set_name, backgal_data[i]);
		}

		// task distribution
		if (numprocs - 1 == rank)
		{
			my_gal_s = foregal_num / numprocs * rank;
			my_gal_e = foregal_num / numprocs * (rank + 1) + foregal_num % numprocs;
		}
		else
		{
			my_gal_s = foregal_num / numprocs * rank;
			my_gal_e = foregal_num / numprocs * (rank + 1);
		}
	
		// loop the foreground galaxy
		for (gal_id = my_gal_s; gal_id < my_gal_e; gal_id++)
		{
			z_f = foregal_data[z_id][gal_id];
			ra_f = foregal_data[ra_id][gal_id];
			dec_f = foregal_data[dec_id][gal_id];

			histogram2d_s(dec_f, ra_f, backgal_data[dec_bin_id], backgal_data[ra_bin_id], dec_bin_num, ra_bin_num, bin_label);
			row = bin_label / grid_nx;
			col = bin_label % grid_nx;

			gal_info.idy = row;
			gal_info.idx = col;
			gal_info.y = dec_f;
			gal_info.x = ra_f;
			gal_info.ny = grid_ny;
			gal_info.nx = grid_nx;
			gal_info.blocks_num = grid_num;
			gal_info.scale = block_scale[0];

			initialize_arr(block_mask, grid_num, -1);

			// loop the search radius
			for (rad_id = 0; rad_id < radius_num; rad_id++)
			{				
				radius_s = radius_bin[rad_id];
				radius_e = radius_bin[rad_id + 1];

				// find the blocks needed
				find_block(&gal_info, radius_s, radius_e, backgal_data[bdy_id], backgal_data[bdx_id], block_mask);

				// loop the found blocks and calculate
				for (i = 0; i < grid_num; i++)
				{
					if (block_mask[i] > -1)
					{

					}
				}


			}
		}


		// free the memory
		for (i = 0; i < 3; i++)
		{
			delete[] foregal_data[i];
		}
		for (i = 0; i < 15; i++)
		{
			delete[] backgal_data[i];
		}
	}
	MPI_Finalize();
	return 0;
}