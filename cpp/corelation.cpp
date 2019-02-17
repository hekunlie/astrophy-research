#include<FQlib.h>
#include<mpi.h>

int main(int argc, char *argv[])
{
	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	// initialize GSL
	int seed;
	gsl_initialize(123);

	pts_info pts_data;
	
	int i, j, k;
	int area_id;
	int ix, iy, ik, ir, ir_1, ig, in, icg;
	int label, label_, label_g1[2], label_g2[2], tag, tag_, tag_g1[2], tag_g2[2];
	int gg_hist_label;

	// read the parameters of the data structure	

	int max_area[1]{}; // sky areas
	int radius_num[1]; // bin number 

	int grid_info[2]{}; // the grid shape to be read of each area
	int grid_ny, grid_nx, grid_num; // number of blocks of each area along each axis
	int max_block_length[2]{}; // [row, col] the maximum data number of all blocks in one sky area
	double block_scale[1]{}; // the length of the block side
	int mg_bin_num[1]{};

	int data_col, block_size; // columns of data & size of array of the biggest block

	// the column of each quantity in the data array
	int ra_idx = 0, dec_idx = 1;
	int mg1_idx = 2, mg2_idx = 3, mnpu_idx = 4, mnmu_idx = 5, mv_idx = 6;

	// for correlation <g1g1`> && <g2g2`>
	double g11, g12;
	double g21, g22;

	double *mg_bins;

	int chi_block_size;
	double *chi_1; 
	double *chi_2;

	int *mask;

	// the bounds of each block
	double *boundy;
	double *boundx;

	int *task_list;
	int *my_tasks;
	int *search_blocks;

	double covs[4]{};
	double mus[2]{};
	double corre_gs[2]{};
	int g_hat_num[1]{};

	double ra, de;
	double dy, dx, distance_sq;

	char h5f_path[150], set_name[50], attrs_name[50], cache_path[150], data_path[150];
	char log_inform[200];

	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/correlation/");
	sprintf(h5f_path, "%scata.hdf5", data_path);
	// read the number of the sky areas
	sprintf(set_name, "/grid");
	sprintf(attrs_name, "max_area");
	read_h5_attrs(h5f_path, set_name, attrs_name, max_area);
	// the bin number
	sprintf(set_name, "/radius");
	sprintf(attrs_name, "bin_num");
	read_h5_attrs(h5f_path, set_name, attrs_name, radius_num);
	// the bins of radian for correlation function
	double *radius = new double[radius_num[0] + 1]{};
	double *radius_sq = new double[radius_num[0] + 1]{};
	read_h5(h5f_path, set_name, radius);
	for (i = 0; i < radius_num[0] + 1; i++)
	{
		radius_sq[i] = radius[i] * radius[i];
	}

	// for the correlated shear 
	sprintf(set_name, "/g_hat_bins");
	sprintf(attrs_name, "bin_num");
	read_h5_attrs(h5f_path, set_name, attrs_name, g_hat_num);
	double *gg_hats = new double[g_hat_num[0]]{};
	double *gg_hats_cov = new double[g_hat_num[0]]{};
	read_h5(h5f_path, set_name, gg_hats);
	for (i = 0; i < g_hat_num[0]; i++)
	{
		gg_hats_cov[i] = fabs(gg_hats[i] * 2);
	}

	double *buffer_ptr;
	int *num_ptr;
	int *err_ptr;
	int arr_len;
	int disp_unit, disp_unit_num, disp_unit_err;
	MPI_Win win, win_num, win_err;
	MPI_Aint buffer_size, size_num, size_err;
	// loop the areas
	for (area_id = 0; area_id < max_area[0]; area_id++)
	{		
		sprintf(set_name, "/grid/w_%d", area_id);

		// read the shape of the grid "w_i"
		sprintf(attrs_name, "grid_info");
		read_h5_attrs(h5f_path, set_name, attrs_name, grid_info);
		grid_ny = grid_info[0];
		grid_nx = grid_info[1];
		grid_num = grid_ny * grid_nx;

		// of the single block
		sprintf(attrs_name, "max_block_length");
		read_h5_attrs(h5f_path, set_name, attrs_name, max_block_length);
		data_col = max_block_length[1];
		// the maximum block size. row * col
		block_size = max_block_length[0] * max_block_length[1];

		sprintf(attrs_name, "block_scale");
		read_h5_attrs(h5f_path, set_name, attrs_name, block_scale);

		sprintf(set_name, "/grid/w_%d/G_bins", area_id);
		sprintf(attrs_name, "bin_num");
		read_h5_attrs(h5f_path, set_name, attrs_name, mg_bin_num);
		mg_bins = new double[mg_bin_num[0] + 1]{};
		read_h5(h5f_path, set_name, mg_bins);
		
		chi_block_size = mg_bin_num[0] * mg_bin_num[0];
		chi_1 = new double[chi_block_size* g_hat_num[0]]{};
		chi_2 = new double[chi_block_size * g_hat_num[0]]{};

		// mask for non-repeating calculation
		mask = new int[max_block_length[0]]{};

		// the bounds of each block
		boundy = new double[grid_num * 4]{};
		boundx = new double[grid_num * 4]{};
		// the origin of the coordinates is the first cross of the grid (array)
		block_bound(block_scale[0], grid_ny, grid_nx, boundy, boundx);

		// each thread gets its own targets blocks in my_tasks
		task_list = new int[grid_num] {};
		my_tasks = new int[grid_num] {};
		search_blocks = new int[grid_num] {};

		// initialize the information of pts_info
		pts_data.scale = block_scale[0];
		pts_data.ny = grid_ny;
		pts_data.nx = grid_nx;
		pts_data.blocks_num = grid_num;

		// the shared buffer in the memory
		arr_len = block_size *grid_num;


		if (rank == 0)
		{
			buffer_size = arr_len * sizeof(double);
			MPI_Win_allocate_shared(buffer_size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &buffer_ptr, &win);
			size_num = grid_num * sizeof(int);
			MPI_Win_allocate_shared(size_num, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &num_ptr, &win_num);
			size_err = sizeof(int);
			MPI_Win_allocate_shared(size_err, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &err_ptr, &win_err);
		}
		else
		{		
			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &buffer_ptr, &win);
			MPI_Win_shared_query(win, 0, &buffer_size, &disp_unit, &buffer_ptr);

			MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &num_ptr, &win_num);
			MPI_Win_shared_query(win_num, 0, &size_num, &disp_unit_num, &num_ptr);

			MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &err_ptr, &win_err);
			MPI_Win_shared_query(win_err, 0, &size_err, &disp_unit_err, &err_ptr);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (0 == rank)
		{
			err_ptr[0] = 1;
			sprintf(log_inform, "Create the shared buffer for program. %d", err_ptr[0]);
			std::cout << log_inform << std::endl;
		}


		// initialize the data in the buffer for all threads
		if (0 == rank)
		{	
			int block_shape[2]{}, num_in_block;
			int block_id;
			double *temp = new double[block_size];

			sprintf(log_inform, "INITIALIZE THE BIG BUFFER FOR W_%d: ", area_id);
			std::cout << log_inform << std::endl;
			sprintf(log_inform, "the biggest blcok: %d x %d", max_block_length[0] , max_block_length[1]);
			std::cout << log_inform << std::endl;

			// read the number of points in each block
			initialize_arr(num_ptr, size_num, 0.);
			sprintf(set_name, "/grid/w_%d/num_in_block", area_id);
			read_h5(h5f_path, set_name, num_ptr);

			// check & cache
			sprintf(cache_path, "%scache/num.fits",data_path);
			write_fits(cache_path, num_ptr, grid_ny, grid_nx);

			initialize_arr(buffer_ptr, buffer_size, 0.);
			for (iy = 0; iy < grid_ny; iy++)
			{
				for (ix = 0; ix < grid_nx; ix++)
				{
					// initialize temp
					initialize_arr(temp, block_size,0.);

					sprintf(set_name, "/grid/w_%d/%d/%d", area_id, iy, ix);

					sprintf(attrs_name, "shape");
					read_h5_attrs(h5f_path, set_name, attrs_name, block_shape);
					num_in_block = block_shape[0] * block_shape[1];

					read_h5(h5f_path, set_name, temp);					

					sprintf(log_inform, "Initialzie the block [%d, %d]. The true shape (%d, %d)", iy, ix, block_shape[0], block_shape[1]);
					std::cout << log_inform << std::endl;

					block_id = iy * grid_nx + ix;
					for (k = 0; k < num_in_block; k++)
					{
						buffer_ptr[block_id*block_size + k] = temp[k];
					}
				}
			}
			delete[] temp;

			// check & cache
			sprintf(cache_path, "%scache/shared_buffer.hdf5", data_path);
			sprintf(set_name, "/data");
			write_h5(cache_path, set_name, buffer_ptr, grid_num*max_block_length[0], max_block_length[1]);

			sprintf(log_inform, "FINISH BUFFER INITIALZIATION");
			std::cout << log_inform << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);

			
		// tasks distribution
		for (k = 0; k < grid_num; k++)
		{
			task_list[k] = k;
		}
		// each thread gets its own task, the blocks to loop, in the "my_tasks".
		// -1 labels the end
		task_alloc(task_list, grid_num, numprocs, rank, my_tasks);

		// loop the target blocks 
		for (k = 0; k < grid_num; k++)
		{
			// initialize the mask for non-repeating calculation
			initialize_arr(mask, max_block_length[0], 1);

			if (my_tasks[k] > -1)
			{
				// the block (iy, ix) of area "w_i"
				iy = my_tasks[k] / grid_nx;
				ix = my_tasks[k] % grid_ny;
				label_ = my_tasks[k] * block_size;

				// loop the points in the block
				// the point number has been stored in num_ptr[my_tasks[k]]
				for (ik = 0; ik < num_ptr[my_tasks[k]]; ik++)
				{
					// used
					mask[ik] = 0;

					// remember 8 parameters in pts_data must be initialized
					// block [idy, idx]
					pts_data.idy = iy;
					pts_data.idx = ix;

					label = label_ + ik * data_col;
					// coordinate of point: y (dec), x (ra) 
					pts_data.y = buffer_ptr[label + dec_idx];
					pts_data.x = buffer_ptr[label + ra_idx];

					label_g1[0] = label + mg1_idx;
					label_g1[1] = label + mnmu_idx;
					label_g2[0] = label + mg2_idx;
					label_g2[1] = label + mnpu_idx;

					// loop the radius scales 
					for (ir = 0; ir < radius_num[0]; ir++)
					{
						ir_1 = ir + 1; // save computational cost
						// find the target blocks
						initialize_arr(search_blocks, grid_num, -1);
						find_block(&pts_data, radius[ir], radius[ir_1], boundy, boundx, search_blocks);
						
						// find the pairs in the searched blocks
						for (ig = 0; ig < grid_num; ig++)
						{
							if (search_blocks[ig] > -1)
							{
								tag_ = search_blocks[ig] * block_size;

								// if the target block is the block where the point belongs to
								if (search_blocks[ig] == my_tasks[k])
								{
									for (in = 0; in < num_ptr[search_blocks[ig]]; in++)
									{
										if (mask[in] == 1)
										{
											tag = tag_ + in * data_col;
											dy = buffer_ptr[tag + dec_idx] - pts_data.y;
											dx = buffer_ptr[tag + ra_idx] - pts_data.x;
											distance_sq = dy * dy + dx * dx;

											tag_g1[0] = tag + mg1_idx;
											tag_g1[1] = tag + mnmu_idx;
											tag_g2[0] = tag + mg2_idx;
											tag_g2[1] = tag + mnpu_idx;

											if (distance_sq >= radius_sq[ir] and distance_sq < radius_sq[ir_1])
											{

												for (icg = 0; icg < g_hat_num[0]; icg++)
												{
													covs[0] = gg_hats_cov[icg];
													covs[1] = gg_hats[icg];
													covs[2] = gg_hats[icg];
													covs[3] = gg_hats_cov[icg];
													rand_multi_gauss(covs, mus, 2, corre_gs);
													// rotation ...
													g11 = buffer_ptr[label_g1[0]] - corre_gs[0]* buffer_ptr[label_g1[1]];
													g12 = buffer_ptr[tag_g1[0]] - corre_gs[1] * buffer_ptr[tag_g1[1]];
													gg_hist_label = histogram2d_s(g11, g12, mg_bins, mg_bins, mg_bin_num[0], mg_bin_num[0]);
													chi_1[gg_hist_label + chi_block_size * icg] += 1;
													rand_multi_gauss(covs, mus, 2, corre_gs);
													// rotation ...
													g21 = buffer_ptr[label_g2[0]] - corre_gs[0] * buffer_ptr[label_g2[1]];
													g22 = buffer_ptr[tag_g2[0]] - corre_gs[1] * buffer_ptr[tag_g1[1]];
													gg_hist_label = histogram2d_s(g21, g22, mg_bins, mg_bins, mg_bin_num[0], mg_bin_num[0]);
													chi_2[gg_hist_label + chi_block_size * icg] += 1;
												}
											}
										}
									}
								}
								else
								{									
									for (in = 0; in < num_ptr[search_blocks[ig]]; in++)
									{
										tag = tag_ + in * data_col;
										dy = buffer_ptr[tag + dec_idx] - pts_data.y;
										dx = buffer_ptr[tag + ra_idx] - pts_data.x;
										distance_sq = dy * dy + dx * dx;

										tag_g1[0] = tag + mg1_idx;
										tag_g1[1] = tag + mnmu_idx;
										tag_g2[0] = tag + mg2_idx;
										tag_g2[1] = tag + mnpu_idx;

										if (distance_sq >= radius_sq[ir] and distance_sq < radius_sq[ir_1])
										{
											for (icg = 0; icg < g_hat_num[0]; icg++)
											{
												covs[0] = gg_hats_cov[icg];
												covs[1] = gg_hats[icg];
												covs[2] = gg_hats[icg];
												covs[3] = gg_hats_cov[icg];

												rand_multi_gauss(covs, mus, 2, corre_gs);
												g11 = buffer_ptr[label_g1[0]] - corre_gs[0] * buffer_ptr[label_g1[1]];
												g12 = buffer_ptr[tag_g1[0]] - corre_gs[1] * buffer_ptr[tag_g1[1]];
												gg_hist_label = histogram2d_s(g11, g12, mg_bins, mg_bins, mg_bin_num[0], mg_bin_num[0]);

												rand_multi_gauss(covs, mus, 2, corre_gs);
												// rotation ...
												g21 = buffer_ptr[label_g2[0]] - corre_gs[0] * buffer_ptr[label_g2[1]];
												g22 = buffer_ptr[tag_g2[0]] - corre_gs[1] * buffer_ptr[tag_g1[1]];
												gg_hist_label = histogram2d_s(g21, g22, mg_bins, mg_bins, mg_bin_num[0], mg_bin_num[0]);

											}
										}
									}
								}
							}
						}

					}

				}

			}
		}
		delete[] chi_1;
		delete[] chi_2;
		delete[] mask;
		delete[] boundy;
		delete[] boundx;
		delete[] task_list;
		delete[] my_tasks;
		delete[] search_blocks;
		delete[] mg_bins;


		MPI_Win_free(&win);
		MPI_Win_free(&win_num);
		MPI_Win_free(&win_err);

	}

	delete[] radius;
	delete[] radius_sq;
	delete[] gg_hats_cov;
	delete[] gg_hats;
	// free the GSL
	gsl_free();
	MPI_Finalize();
	return 0;
}