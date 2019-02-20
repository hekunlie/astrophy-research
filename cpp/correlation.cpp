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


	pts_info pts_data;
	
	char h5f_path[150], set_name[50], attrs_name[50], cache_path[150], data_path[150];
	char log_inform[200], chi_path[200], chi_set_name[30];

	int i, j, k;
	int task_num;
	int area_id;
	int ix, iy, ik, ir, ir_1, ig, in, icg;
	int pts_label, pts_label_, pts_s_tag, pts_s_tag_, pts_s_tags[5];
	int gg_hist_label;

	// the column of each quantity in the data array
	int ra_idx = 0, dec_idx = 1;
	int mg1_idx = 2, mg2_idx = 3, mn_idx = 4, mu_idx = 5, mv_idx = 6;
	double dy, dx, distance_sq;
	double sin2a, cos2a, sin4a, cos4a, gamma, gamma2;
	double mg1_r, mg2_r, mnu1_r, mnu2_r, mu_r;
	double mg1, mg2, mnu1, mnu2;

	// the parameters of the data structure	
	int max_area[1]{}; // sky areas
	int radius_num[1]; // bin number 
	double *radius, *radius_sq;

	// for the correlation bins
	int g_hat_num[1]{};
	double *gg_hats;
	double *gg_hats_cov;

	// the variables for shared memory buffer
	double *buffer_ptr;
	int *num_ptr, *err_ptr;
	int arr_len;
	int disp_unit, disp_unit_num, disp_unit_err;
	MPI_Win win, win_num, win_err;
	MPI_Aint buffer_size, size_num, size_err;


	int grid_info[2]{}; // the grid shape to be read of each area
	int grid_ny, grid_nx, grid_num; // number of blocks of each area along each axis
	
	int max_block_length[2]{}; // [row, col] the maximum data number of all blocks in one sky area
	int data_col, block_size; // columns of data & size of array of the biggest block
	
	double block_scale[1]{}; // the length of the block side
	
	// the bounds of each block
	double *boundy;
	double *boundx;

	int mg_bin_num[1]{};
	double *mg_bins;

	int chi_block_size, chi_block_size_ir, chi_label;
	long *chi_1;
	long *chi_2;

	// for correlation <g1g1`> && <g2g2`>
	double mg11, mg12;
	double mg21, mg22;
	double covs[4]{};
	double mus[2]{};
	double corre_gs[2]{};


	int *mask;
	int *task_list;
	int *my_tasks;
	int *search_blocks;


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
	radius = new double[radius_num[0] + 1]{};
	radius_sq = new double[radius_num[0] + 1]{};
	read_h5(h5f_path, set_name, radius);
	for (i = 0; i < radius_num[0] + 1; i++)
	{
		radius_sq[i] = radius[i] * radius[i];
	}

	// the correlation bins 
	sprintf(set_name, "/g_hat_bin");
	sprintf(attrs_name, "bin_num");
	read_h5_attrs(h5f_path, set_name, attrs_name, g_hat_num);
	gg_hats = new double[g_hat_num[0]]{};
	gg_hats_cov = new double[g_hat_num[0]]{};
	read_h5(h5f_path, set_name, gg_hats);
	for (i = 0; i < g_hat_num[0]; i++)
	{
		gg_hats_cov[i] = fabs(gg_hats[i] * 2);
	}


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
		chi_block_size_ir = chi_block_size * g_hat_num[0];
		chi_1 = new long[chi_block_size_ir* radius_num[0]]{};
		chi_2 = new long[chi_block_size_ir * radius_num[0]]{};

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

		// initialize the data in the buffer for all threads
		if (0 == rank)
		{	
			sprintf(log_inform, "Created the shared buffer for program.");
			std::cout << log_inform << std::endl;			

			sprintf(log_inform, "INITIALIZE THE BIG BUFFER FOR W_%d: ", area_id);
			std::cout << log_inform << std::endl;
			sprintf(log_inform, "the biggest blcok: %d x %d", max_block_length[0] , max_block_length[1]);
			std::cout << log_inform << std::endl;

			// read the number of points in each block
			initialize_arr(num_ptr, size_num, 0);
			sprintf(set_name, "/grid/w_%d/num_in_block", area_id);
			read_h5(h5f_path, set_name, num_ptr);
			// check & cache
			sprintf(cache_path, "%scache/num.fits",data_path);
			write_fits(cache_path, num_ptr, grid_ny, grid_nx);

			//read the data
			initialize_arr(buffer_ptr, buffer_size, 0.);
			sprintf(set_name, "/grid/w_%d/data", area_id);
			read_h5(h5f_path, set_name, buffer_ptr);
			// check & cache
			sprintf(cache_path, "%scache/shared_buffer.hdf5", data_path);
			sprintf(set_name, "/data");
			write_h5(cache_path, set_name, buffer_ptr, grid_num*max_block_length[0], max_block_length[1]);

			sprintf(log_inform, "FINISH BUFFER INITIALZIATION");
			std::cout << log_inform << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);

			
		// tasks distribution
		// each thread gets its own task, the blocks to loop, in the "my_tasks".
		task_num = 0;
		for (k = 0; k < grid_num; k++)
		{
			// skip the blocks contains nothing
			// because the observation areas are irregular
			// but our blocks are square
			if (num_ptr[k] > 0)
			{
				task_list[task_num] = k;
				task_num++;
			}
		}
		if (task_num == 0)
		{
			std::cout << "Something went wrong in task_num!!!" <<std::endl;
			exit(0);
		}
		// -1 denotes the end
		task_alloc(task_list, task_num, numprocs, rank, my_tasks);

		// loop the target blocks 
		for (k = 0; k < grid_num; k++)
		{
			if (my_tasks[k] > -1)
			{
				seed = (rank + 1)*1000000+ rank* my_tasks[k]*10;
				gsl_initialize(seed);

				// initialize the mask for non-repeating calculation
				initialize_arr(mask, max_block_length[0], 1);
				// the block (iy, ix) of area "w_i"
				iy = my_tasks[k] / grid_nx;
				ix = my_tasks[k] % grid_ny;
				pts_label_ = my_tasks[k] * block_size;

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

					pts_label = pts_label_ + ik * data_col;
					// coordinate of point: y (dec), x (ra) 
					pts_data.y = buffer_ptr[pts_label + dec_idx];
					pts_data.x = buffer_ptr[pts_label + ra_idx];

					mg1 = buffer_ptr[pts_label + mg1_idx];
					mg2 = buffer_ptr[pts_label + mg2_idx];
					// mnu1 = mn - mu, mnu2 = mn + mu for the quantity from the pipeline
					// while, the signs of mn, mu, which are different from that in paper, 
					//  have been corrected by Python.
					mnu1 = buffer_ptr[pts_label + mn_idx] + buffer_ptr[pts_label + mu_idx];
					mnu2 = buffer_ptr[pts_label + mn_idx] - buffer_ptr[pts_label + mu_idx];

					// loop the radius scales 
					for (ir = 0; ir < radius_num[0]; ir++)
					{
						// save computational cost
						ir_1 = ir + 1; 
						// find the target blocks
						initialize_arr(search_blocks, grid_num, -1);
						find_block(&pts_data, radius[ir], radius[ir_1], boundy, boundx, search_blocks);
						chi_label = chi_block_size_ir * ir;

						// find the pairs in the searched blocks
						for (ig = 0; ig < grid_num; ig++)
						{
							if (search_blocks[ig] > -1)
							{
								pts_s_tag_ = search_blocks[ig] * block_size;

								// if the target block is the block where the point belongs to
								if (search_blocks[ig] == my_tasks[k])
								{
									for (in = 0; in < num_ptr[search_blocks[ig]]; in++)
									{
										pts_s_tag = pts_s_tag_ + in * data_col;
										dy = buffer_ptr[pts_s_tag + dec_idx] - pts_data.y;
										dx = buffer_ptr[pts_s_tag + ra_idx] - pts_data.x;
										distance_sq = dy * dy + dx * dx;

										if (mask[in] == 1)
										{
											if (distance_sq >= radius_sq[ir] and distance_sq < radius_sq[ir_1])
											{
												// for the rotation
												gamma = dx / dy;
												gamma2 = gamma * gamma;
												// for mg1, mg2
												cos2a = (gamma2 - 1) / (gamma2 + 1);
												sin2a = 2 * gamma / (gamma2 + 1);
												// for mu, mv
												cos4a = (cos2a + sin2a)*(cos2a - sin2a);
												sin4a = 2 * cos2a*sin2a;

												// for saving computational cost
												pts_s_tags[0] = pts_s_tag + mg1_idx;
												pts_s_tags[1] = pts_s_tag + mg2_idx;
												pts_s_tags[2] = pts_s_tag + mn_idx;
												pts_s_tags[3] = pts_s_tag + mu_idx;
												pts_s_tags[4] = pts_s_tag + mv_idx;

												mg1_r = buffer_ptr[pts_s_tags[0]] * cos2a - buffer_ptr[pts_s_tags[1]] * sin2a;
												mg2_r = buffer_ptr[pts_s_tags[0]] * sin2a + buffer_ptr[pts_s_tags[1]] * cos2a;
												mu_r = buffer_ptr[pts_s_tags[3]] * cos4a - buffer_ptr[pts_s_tags[4]] * sin4a;
												mnu1_r = buffer_ptr[pts_s_tags[2]] + mu_r;
												mnu2_r = buffer_ptr[pts_s_tags[2]] - mu_r;
				
												for (icg = 0; icg < g_hat_num[0]; icg++)
												{
													covs[0] = gg_hats_cov[icg];
													covs[1] = gg_hats[icg];
													covs[2] = gg_hats[icg];
													covs[3] = gg_hats_cov[icg];
													rand_multi_gauss(covs, mus, 2, corre_gs);

													mg11 = mg1 - corre_gs[0]* mnu1;
													mg12 = mg1_r - corre_gs[1] * mnu1_r;

													gg_hist_label = histogram2d_s(mg11, mg12, mg_bins, mg_bins, mg_bin_num[0], mg_bin_num[0]);
													chi_1[chi_label + chi_block_size * icg + gg_hist_label] += 1;

													//rand_multi_gauss(covs, mus, 2, corre_gs);

													mg21 = mg2 - corre_gs[0] *mnu2;
													mg22 = mg2_r - corre_gs[1] * mnu2_r;
													gg_hist_label = histogram2d_s(mg21, mg22, mg_bins, mg_bins, mg_bin_num[0], mg_bin_num[0]);
													chi_2[chi_label + chi_block_size * icg + gg_hist_label] += 1;
												}
											}
										}
									}
								}
								else
								{									
									for (in = 0; in < num_ptr[search_blocks[ig]]; in++)
									{
										pts_s_tag = pts_s_tag_ + in * data_col;
										dy = buffer_ptr[pts_s_tag + dec_idx] - pts_data.y;
										dx = buffer_ptr[pts_s_tag + ra_idx] - pts_data.x;
										distance_sq = dy * dy + dx * dx;

										if (distance_sq >= radius_sq[ir] and distance_sq < radius_sq[ir_1])
										{

											// for the rotation
											gamma = dx / dy;
											gamma2 = gamma * gamma;
											// for mg1, mg2
											cos2a = (gamma2 - 1) / (gamma2 + 1);
											sin2a = 2 * gamma / (gamma2 + 1);
											// for mu, mv
											cos4a = (cos2a + sin2a)*(cos2a - sin2a);
											sin4a = 2 * cos2a*sin2a;

											// for saving computational cost
											pts_s_tags[0] = pts_s_tag + mg1_idx;
											pts_s_tags[1] = pts_s_tag + mg2_idx;
											pts_s_tags[2] = pts_s_tag + mn_idx;
											pts_s_tags[3] = pts_s_tag + mu_idx;
											pts_s_tags[4] = pts_s_tag + mv_idx;

											mg1_r = buffer_ptr[pts_s_tags[0]] * cos2a - buffer_ptr[pts_s_tags[1]] * sin2a;
											mg2_r = buffer_ptr[pts_s_tags[0]] * sin2a + buffer_ptr[pts_s_tags[1]] * cos2a;
											mu_r = buffer_ptr[pts_s_tags[3]] * cos4a - buffer_ptr[pts_s_tags[4]] * sin4a;
											mnu1_r = buffer_ptr[pts_s_tags[2]] + mu_r;
											mnu2_r = buffer_ptr[pts_s_tags[2]] - mu_r;

											for (icg = 0; icg < g_hat_num[0]; icg++)
											{
												covs[0] = gg_hats_cov[icg];
												covs[1] = gg_hats[icg];
												covs[2] = gg_hats[icg];
												covs[3] = gg_hats_cov[icg];
												rand_multi_gauss(covs, mus, 2, corre_gs);

												mg11 = mg1 - corre_gs[0] * mnu1;
												mg12 = mg1_r - corre_gs[1] * mnu1_r;

												gg_hist_label = histogram2d_s(mg11, mg12, mg_bins, mg_bins, mg_bin_num[0], mg_bin_num[0]);
												chi_1[chi_label + chi_block_size * icg + gg_hist_label] += 1;

												//rand_multi_gauss(covs, mus, 2, corre_gs);

												mg21 = mg2 - corre_gs[0] * mnu2;
												mg22 = mg2_r - corre_gs[1] * mnu2_r;
												gg_hist_label = histogram2d_s(mg21, mg22, mg_bins, mg_bins, mg_bin_num[0], mg_bin_num[0]);
												chi_2[chi_label + chi_block_size * icg + gg_hist_label] += 1;
											}
										}
									}
								}
							}
						}

					}

				}
			// free the GSL
			gsl_free();
			}

		}
		// save the chi array
		sprintf(chi_path, "%scache/chi_1_%d.hdf5", area_id);
		sprintf(chi_set_name, "/data");
		write_h5(chi_path, chi_set_name, chi_1, mg_bin_num[0], mg_bin_num[0] * g_hat_num[0]*radius_num[0]);

		sprintf(chi_path, "%scache/chi_2_%d.hdf5", area_id);
		write_h5(chi_path, chi_set_name, chi_2, mg_bin_num[0], mg_bin_num[0] * g_hat_num[0]*radius_num[0]);

		delete[] boundy;
		delete[] boundx;
		delete[] mg_bins;
		delete[] chi_1;
		delete[] chi_2;
		delete[] mask;
		delete[] task_list;
		delete[] my_tasks;
		delete[] search_blocks;

		MPI_Win_free(&win);
		MPI_Win_free(&win_num);
		MPI_Win_free(&win_err);

	}

	delete[] radius;
	delete[] radius_sq;
	delete[] gg_hats;
	delete[] gg_hats_cov;

	MPI_Finalize();
	return 0;
}