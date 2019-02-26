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

	char h5f_path[200], set_name[50], attrs_name[50], cache_path[200], data_path[200];
	char log_inform[200], log_path[200],chi_path[200], chi_set_name[30];
	char temp[200];

	int i, j, k, i_t;
	
	int area_id, block_id, block_id_s;
	int ix, iy, ik, ir, ir_1, ig, in, icg;	
	double st1, st2, st3, st4, st5, t_block_s, t_block_e, t_area_s, t_area_e;
	int int_attrs[2]{};
	double double_attrs[2]{};

	st1 = clock();
	// the column of each quantity in the data array
	int ra_idx = 0, dec_idx = 1;
	int mg1_idx = 2, mg2_idx = 3, mn_idx = 4, mu_idx = 5, mv_idx = 6;
	double dy, dx, distance_sq;
	double sin2a, cos2a, sin4a, cos4a, gamma, gamma2;
	double mg1_2r, mg2_2r, mnu1_2r, mnu2_2r, mu2_2r;
	double mg1_1r, mg2_1r, mnu1_1r, mnu2_1r, mu2_1r;
	double mg1_1, mg2_1, mn1_1, mu1_1, mv1_1;

	// the parameters of the data structure
	int max_area; // sky areas
	int radius_bin_num;// bin number = radius_num - 1
	double *radius, *radius_sq;

	// for the correlation bins
	int g_hat_num;
	int gg_hist_label;
	double *gg_hats;
	double *gg_hats_cov;

	// the variables for shared memory buffer
	//MPI_Win win_mg1, win_mg2, win_mn, win_mu, win_mv, win_ra, win_dec;
	//int disp_mg1, disp_mg2, disp_mn, disp_mu, disp_mv, disp_ra, disp_dec;
		//MPI_Aint size_pts;
	double *MG1, *MG2, *MN, *MU, *MV, *RA, *DEC;
	int data_num;


	int *num_in_block, *block_start, *block_end;
	int task_num;
	int grid_ny, grid_nx, grid_num; // number of blocks of each area along each axis

	int block_row, block_col, block_size; // columns of data & size of array of the biggest block

	double block_scale; // the length of the block side

	// the bounds of each block
	double *boundy;
	double *boundx;

	int mg_bin_len[1]{};// number of the boundary of bins, bin_num+1
	int mg_bin_num; // bin number
	double *mg_bins; //  the boundary of each bin

	int chi_block_size, chi_block_size_ir, chi_label;
	long *chi_1, *chi_2;
	long *chi_1_total, *chi_2_total;
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

	sprintf(log_path, "/mnt/ddnfs/data_users/hkli/CFHT/correlation/logs/log_%d.dat", rank);
	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/data/");

	//sprintf(log_path, "/mnt/perc/hklee/CFHT/correlation/logs/log_%d.dat", rank);
	//sprintf(data_path, "/mnt/perc/hklee/CFHT/data/");

	sprintf(h5f_path, "%scf_cata_result_int.hdf5", data_path);
	// read the number of the sky areas
	sprintf(set_name, "/grid");
	sprintf(attrs_name, "max_area");
	read_h5_attrs(h5f_path, set_name, attrs_name, int_attrs,"g");
	max_area = int_attrs[0];

	// the bin number
	sprintf(set_name, "/radius_bin");
	sprintf(attrs_name, "shape");
	read_h5_attrs(h5f_path, set_name, attrs_name, int_attrs,"d");
	radius_bin_num = int_attrs[0] - 1;
	// the bins of radian for correlation function
	radius = new double[int_attrs[0]]{};
	radius_sq = new double[int_attrs[0]]{};
	read_h5(h5f_path, set_name, radius);
	for (i = 0; i < int_attrs[0]; i++)
	{
		radius_sq[i] = radius[i] * radius[i];
	}

	// the correlation bins
	sprintf(set_name, "/g_hat_bin");
	sprintf(attrs_name, "shape");
	read_h5_attrs(h5f_path, set_name, attrs_name, int_attrs, "d");
	g_hat_num = int_attrs[0];
	gg_hats = new double[g_hat_num]{};
	gg_hats_cov = new double[g_hat_num]{};
	read_h5(h5f_path, set_name, gg_hats);
	for (i = 0; i < g_hat_num; i++)
	{
		gg_hats_cov[i] = fabs(gg_hats[i] * 2);
	}

	sprintf(log_inform, "RANK: %d begins to loop areas", rank);
	write_log(log_path, log_inform);
	if (rank == 0)
	{
		std::cout << log_inform << std::endl;
	}


	// loop the areas
	for (area_id = 1; area_id < max_area+1; area_id++)
	{
		t_area_s = clock();

		sprintf(log_inform, "RANK: %d read the parameters in area %d", rank, area_id);
		write_log(log_path, log_inform);

		sprintf(set_name, "/grid/w_%d", area_id);

		// read the shape of the grid "w_i"
		sprintf(attrs_name, "grid_shape");
		read_h5_attrs(h5f_path, set_name, attrs_name, int_attrs,"g");
		grid_ny = int_attrs[0];
		grid_nx = int_attrs[1];
		grid_num = grid_ny * grid_nx;

		// of the single block
		sprintf(attrs_name, "max_block_size");
		read_h5_attrs(h5f_path, set_name, attrs_name, int_attrs,"g");
		block_row = int_attrs[0];
		block_col = int_attrs[1];
		// the maximum block size. row * col
		block_size = int_attrs[0] * int_attrs[1];
		// mask for non-repeating calculation
		mask = new int[block_row]{};

		sprintf(attrs_name, "block_scale");
		read_h5_attrs(h5f_path, set_name, attrs_name, double_attrs,"g");
		block_scale = double_attrs[0];

		sprintf(set_name, "/grid/w_%d/G1_bin", area_id);
		sprintf(attrs_name, "shape");
		read_h5_attrs(h5f_path, set_name, attrs_name, int_attrs, "d");
		mg_bins = new double[int_attrs[0]]{};
		mg_bin_num = int_attrs[0] - 1;
		read_h5(h5f_path, set_name, mg_bins);

		chi_block_size = mg_bin_num* mg_bin_num;
		chi_block_size_ir = chi_block_size * g_hat_num;
		chi_1 = new long[chi_block_size_ir* radius_bin_num]{};
		chi_2 = new long[chi_block_size_ir * radius_bin_num]{};
		chi_1_total = new long[chi_block_size_ir* radius_bin_num*numprocs]{};
		chi_2_total = new long[chi_block_size_ir* radius_bin_num*numprocs]{};

		// the bounds of each block
		boundy = new double[grid_num * 4]{};
		boundx = new double[grid_num * 4]{};
		// the origin of the coordinates is the first cross of the grid (array)
		sprintf(set_name, "/grid/w_%d/boundy", area_id);
		read_h5(h5f_path, set_name, boundy);
		sprintf(set_name, "/grid/w_%d/boundx", area_id);
		read_h5(h5f_path, set_name, boundx);
		//block_bound(block_scale[0], grid_ny, grid_nx, boundy, boundx);

		// each thread gets its own targets blocks in my_tasks
		task_list = new int[grid_num] {};
		my_tasks = new int[grid_num] {};
		search_blocks = new int[grid_num] {};

		// initialize the information of pts_info
		pts_data.scale = block_scale;
		pts_data.ny = grid_ny;
		pts_data.nx = grid_nx;
		pts_data.blocks_num = grid_num;

		// the shared buffer in the memory
		sprintf(set_name, "/grid/w_%d/data/G1", area_id);
		sprintf(attrs_name, "shape");
		read_h5_attrs(h5f_path, set_name, attrs_name, int_attrs, "d");
		data_num = int_attrs[0];

		RA = new double[data_num] {};
		DEC = new double[data_num] {};

		MG1 = new double[data_num] {};
		MG2 = new double[data_num] {};
		MN = new double[data_num] {};
		MU = new double[data_num] {};
		MV = new double[data_num] {};

		// read the number of points in each block
		num_in_block = new int[grid_num] {};
		block_start = new int[grid_num] {};
		block_end = new int[grid_num] {};
		initialize_arr(num_in_block, grid_num, 0);
		sprintf(set_name, "/grid/w_%d/num_in_block", area_id);
		read_h5(h5f_path, set_name, num_in_block);

		initialize_arr(block_start, grid_num, 0);
		sprintf(set_name, "/grid/w_%d/block_start", area_id);
		read_h5(h5f_path, set_name, block_start);

		initialize_arr(block_end, grid_num, 0);
		sprintf(set_name, "/grid/w_%d/block_end", area_id);
		read_h5(h5f_path, set_name, block_end);

		for (k = 0; k < numprocs; k++)
		{
			if (k == rank)
			{	
				sprintf(log_inform, "rank: %d, area: %d, radius_num: %d, g_num: %d, grid: %d x %d\(%d)", rank, max_area, radius_bin_num, g_hat_num, grid_ny, grid_nx, grid_num);
				std::cout <<log_inform<< std::endl;
				sprintf(log_inform, "rank: %d, block scale: %.1f, G_bins: %d, data_num: %d ", rank, block_scale, mg_bin_num, data_num);
				std::cout << log_inform << std::endl;
				std::cout << std::endl;
			}			
			MPI_Barrier(MPI_COMM_WORLD);
		}

	/*	if (rank == 0)
		{
			size_pts = data_num * sizeof(double);
			MPI_Win_allocate_shared(size_pts, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &MG1, &win_mg1);
			MPI_Win_allocate_shared(size_pts, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &MG2, &win_mg2);
			MPI_Win_allocate_shared(size_pts, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &MN, &win_mn);
			MPI_Win_allocate_shared(size_pts, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &MU, &win_mu);
			MPI_Win_allocate_shared(size_pts, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &MV, &win_mv);
			MPI_Win_allocate_shared(size_pts, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &RA, &win_ra);
			MPI_Win_allocate_shared(size_pts, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &DEC, &win_dec);

		}
		else
		{
			MPI_Win_allocate_shared(size_pts, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &MG1, &win_mg1);
			MPI_Win_shared_query(win_mg1, 0, &size_pts, &disp_mg1, &MG1);

			MPI_Win_allocate_shared(size_pts, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &MG2, &win_mg2);
			MPI_Win_shared_query(win_mg2, 0, &size_pts, &disp_mg2, &MG2);

			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &MN, &win_mn);
			MPI_Win_shared_query(win_mn, 0, &size_pts, &disp_mn, &MN);

			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &MU, &win_mu);
			MPI_Win_shared_query(win_mu, 0, &size_pts, &disp_mu, &MU);

			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &MV, &win_mv);
			MPI_Win_shared_query(win_mv, 0, &size_pts, &disp_mv, &MV);

			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &RA, &win_ra);
			MPI_Win_shared_query(win_ra, 0, &size_pts, &disp_ra, &RA);

			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &DEC, &win_dec);
			MPI_Win_shared_query(win_dec, 0, &size_pts, &disp_dec, &DEC);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// initialize the data in the buffer for all threads
		if (0 == rank)
		{
			sprintf(log_inform, "Created the shared buffer for program.");
			std::cout << log_inform << std::endl;

			sprintf(log_inform, "INITIALIZE THE BIG BUFFER FOR W_%d: ", area_id);
			std::cout << log_inform << std::endl;
			sprintf(log_inform, "the biggest blcok: %d x %d", max_block_size[0] , max_block_size[1]);
			std::cout << log_inform << std::endl;

			// check & cache
			sprintf(cache_path, "%scache/num_%d.fits",data_path,area_id);
			write_fits(cache_path, num_in_block, grid_ny, grid_nx);
			sprintf(cache_path, "%scache/block_start_%d.fits", data_path, area_id);
			write_fits(cache_path, block_start, grid_ny, grid_nx);
			sprintf(cache_path, "%scache/block_end_%d.fits", data_path, area_id);
			write_fits(cache_path, block_end, grid_ny, grid_nx);

			//read the data
			initialize_arr(MG1, data_num, 0.);
			sprintf(set_name, "/grid/w_%d/data/G1", area_id);
			read_h5(h5f_path, set_name, MG1);

			initialize_arr(MG2, data_num, 0.);
			sprintf(set_name, "/grid/w_%d/data/G2", area_id);
			read_h5(h5f_path, set_name, MG2);

			initialize_arr(MN, data_num, 0.);
			sprintf(set_name, "/grid/w_%d/data/N", area_id);
			read_h5(h5f_path, set_name, MN);

			initialize_arr(MU, data_num, 0.);
			sprintf(set_name, "/grid/w_%d/data/U", area_id);
			read_h5(h5f_path, set_name, MU);

			initialize_arr(MV, data_num, 0.);
			sprintf(set_name, "/grid/w_%d/data/U", area_id);
			read_h5(h5f_path, set_name, MV);

			initialize_arr(RA, data_num, 0.);
			sprintf(set_name, "/grid/w_%d/data/RA", area_id);
			read_h5(h5f_path, set_name, RA);

			initialize_arr(DEC, data_num, 0.);
			sprintf(set_name, "/grid/w_%d/data/DEC", area_id);
			read_h5(h5f_path, set_name, DEC);

			sprintf(log_inform, "FINISH BUFFER INITIALZIATION");
			std::cout << log_inform << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		sprintf(log_inform, "RANK: %d shared buffer created in area %d", rank, area_id);
		write_log(log_path, log_inform);*/



		sprintf(log_inform, "INITIALIZE THE BIG BUFFER FOR W_%d: the biggest blcok: %d x %d", area_id, block_row, block_col);
		write_log(log_path, log_inform);
		if (rank == 0)
		{
			std::cout << log_inform << std::endl;
			// check & cache
			sprintf(cache_path, "!%scache/num_%d.fits", data_path, area_id);
			write_fits(cache_path, num_in_block, grid_ny, grid_nx);
			sprintf(cache_path, "!%scache/block_start_%d.fits", data_path, area_id);
			write_fits(cache_path, block_start, grid_ny, grid_nx);
			sprintf(cache_path, "!%scache/block_end_%d.fits", data_path, area_id);
			write_fits(cache_path, block_end, grid_ny, grid_nx);
		}

		//read the data
		initialize_arr(MG1, data_num, 0.);
		sprintf(set_name, "/grid/w_%d/data/G1", area_id);
		read_h5(h5f_path, set_name, MG1);

		initialize_arr(MG2, data_num, 0.);
		sprintf(set_name, "/grid/w_%d/data/G2", area_id);
		read_h5(h5f_path, set_name, MG2);

		initialize_arr(MN, data_num, 0.);
		sprintf(set_name, "/grid/w_%d/data/N", area_id);
		read_h5(h5f_path, set_name, MN);

		initialize_arr(MU, data_num, 0.);
		sprintf(set_name, "/grid/w_%d/data/U", area_id);
		read_h5(h5f_path, set_name, MU);

		initialize_arr(MV, data_num, 0.);
		sprintf(set_name, "/grid/w_%d/data/V", area_id);
		read_h5(h5f_path, set_name, MV);

		initialize_arr(RA, data_num, 0.);
		sprintf(set_name, "/grid/w_%d/data/RA", area_id); 
		read_h5(h5f_path, set_name, RA);

		initialize_arr(DEC, data_num, 0.);
		sprintf(set_name, "/grid/w_%d/data/DEC", area_id);
		read_h5(h5f_path, set_name, DEC);

		sprintf(log_inform, "FINISH BUFFER INITIALZIATION");
		write_log(log_path, log_inform);
		if (rank == 0)
		{
			std::cout << log_inform << std::endl;
		}

		// tasks distribution
		// each thread gets its own task, the blocks to loop, in the "my_tasks".
		task_num = 0;
		for (k = 0; k < grid_num; k++)
		{
			// skip the blocks contains nothing
			// because the observation areas are irregular
			// but our blocks are square
			if (num_in_block[k] > 0)
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
        for (k=0; k<grid_num; k++)
        {
            my_tasks[k] = -1;
        }
		task_alloc(task_list, task_num, numprocs, rank, my_tasks);
		//for (k = 0; k < numprocs; k++)
		//{
		//	if (k == rank)
		//	{
		//		std::cout << rank << " My tasks:  "<<task_num<<"  ";
		//		for (ik = 0; ik < grid_num; ik++)
		//		{
		//			if (my_tasks[ik] > -1)
		//			{
		//				std::cout << my_tasks[ik] << ", ";
		//			}
		//		}
		//		std::cout << std::endl;
		//	}
		//	MPI_Barrier(MPI_COMM_WORLD);
		//}
	
		sprintf(log_inform, "RANK: %d begin to loop the task blocks in area %d", rank, area_id);
		write_log(log_path, log_inform);
	
		// loop the target blocks
		for (k = 0; k < grid_num; k++)
		{
			block_id = my_tasks[k];

			if (block_id > -1)
			{
				sprintf(log_inform, "RANK: %d begin to loop the %d 'th block in area %d. ", rank, block_id, area_id);
				write_log(log_path, log_inform);

				seed = (rank + 1)*100000+ rank* block_id*2;
				gsl_initialize(seed);

				// initialize the mask for non-repeating calculation
				//initialize_arr(mask, max_block_size[0], 1);
				// the block (iy, ix) of area "w_i"
				iy = block_id / grid_nx;
				ix = block_id % grid_nx;

				t_block_s = clock();
				// loop the points in the block
				// the point number has been stored in num_ptr[my_tasks[k]]
				for (ik = block_start[block_id]; ik < block_end[block_id]; ik++)
				{	
					// used
					mask[ik - block_start[block_id]] = 0;

					// remember 8 parameters in pts_data must be initialized
					// block [idy, idx]
					pts_data.idy = iy;
					pts_data.idx = ix;

					// coordinate of point: y (dec), x (ra)
					pts_data.y = DEC[ik];
					pts_data.x = RA[ik];

					mg1_1 = MG1[ik];
					mg2_1 = MG2[ik];
					// mnu1 = mn - mu, mnu2 = mn + mu for the quantity from the pipeline
					// while, the signs of mn, mu, which are different from that in paper,
					//  have been corrected by Python.
					mn1_1 = MN[ik];
					mu1_1 = MU[ik];
					mv1_1 = MV[ik];

					// loop the radius scales
					for (ir = 0; ir < radius_bin_num; ir++)
					{
						// save computational cost
						ir_1 = ir + 1;
						// find the target blocks
						initialize_arr(search_blocks, grid_num, -1);
						find_block(&pts_data, radius[ir], radius[ir_1], boundy, boundx, search_blocks);
						chi_label = chi_block_size_ir * ir;
						
						for (i_t = 0; i_t < grid_num; i_t++)
						{
							if (search_blocks[i_t] >= grid_num)
							{
								std::cout << "CROSS: " << i_t<<", "<< search_blocks[i_t]<<", "<<rank << std::endl;
								exit(0);
							}
						}

						// find the pairs in the searched blocks
						for (ig = 0; ig < grid_num; ig++)
						{
							block_id_s = search_blocks[ig];
							if (block_id_s > -1)
							{
								// if the target block is the block where the point belongs to
								if (block_id_s == block_id)
								{
									for (in = block_start[block_id_s]; in < block_end[block_id_s]; in++)
									{
										dy = DEC[in] - pts_data.y;
										dx = RA[in] - pts_data.x;
										distance_sq = dy * dy + dx * dx;
										
										if (mask[in- block_start[block_id_s]] == 1 and distance_sq >= radius_sq[ir] and distance_sq < radius_sq[ir_1])
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
						
											mg1_1r = mg1_1 * cos2a - mg2_1 * sin2a;
											mg2_1r = mg1_1 * sin2a + mg2_1 * cos2a;
											mu2_1r = mu1_1 * cos4a - mv1_1 * sin4a;
											mnu1_1r = mn1_1 + mu1_1;
											mnu2_1r = mn1_1 - mu1_1;

											mg1_2r = MG1[in] * cos2a - MG2[in] * sin2a;
											mg2_2r = MG1[in] * sin2a + MG2[in] * cos2a;
											mu2_2r = MU[in] * cos4a - MV[in] * sin4a;
											mnu1_2r = MN[in] + mu2_2r;
											mnu2_2r = MN[in] - mu2_2r;

											for (icg = 0; icg < g_hat_num; icg++)
											{
												covs[0] = gg_hats_cov[icg];
												covs[1] = gg_hats[icg];
												covs[2] = gg_hats[icg];
												covs[3] = gg_hats_cov[icg];
												rand_multi_gauss(covs, mus, 2, corre_gs);

												mg11 = mg1_1r - corre_gs[0]* mnu1_1r;
												mg12 = mg1_2r - corre_gs[1] * mnu1_2r;

												gg_hist_label = histogram2d_s(mg11, mg12, mg_bins, mg_bins, mg_bin_num, mg_bin_num);
												//if (gg_hist_label >= chi_block_size)
												//{
												//	std::cout << "CROSS " << gg_hist_label << std::endl;
												//}
												chi_1[chi_label + chi_block_size * icg + gg_hist_label] += 1;

												//rand_multi_gauss(covs, mus, 2, corre_gs);

												mg21 = mg2_1r - corre_gs[0] * mnu2_1r;
												mg22 = mg2_2r - corre_gs[1] * mnu2_2r;
												gg_hist_label = histogram2d_s(mg21, mg22, mg_bins, mg_bins, mg_bin_num, mg_bin_num);
												//if (gg_hist_label >= chi_block_size)
												//{
												//	std::cout << "CROSS " << gg_hist_label << std::endl;
												//}
												chi_2[chi_label + chi_block_size * icg + gg_hist_label] += 1;

											}

										}
									}
								}
								else
								{
									for (in = block_start[block_id_s]; in < block_end[block_id_s]; in++)
									{
										dy = DEC[in] - pts_data.y;
										dx = RA[in] - pts_data.x;
										distance_sq = dy * dy + dx * dx;
										
										if (mask[in - block_start[block_id_s]] == 1 and distance_sq >= radius_sq[ir] and distance_sq < radius_sq[ir_1])
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
						
											mg1_1r = mg1_1 * cos2a - mg2_1 * sin2a;
											mg2_1r = mg1_1 * sin2a + mg2_1 * cos2a;
											mu2_1r = mu1_1 * cos4a - mv1_1 * sin4a;
											mnu1_1r = mn1_1 + mu1_1;
											mnu2_1r = mn1_1 - mu1_1;

											mg1_2r = MG1[in] * cos2a - MG2[in] * sin2a;
											mg2_2r = MG1[in] * sin2a + MG2[in] * cos2a;
											mu2_2r = MU[in] * cos4a - MV[in] * sin4a;
											mnu1_2r = MN[in] + mu2_2r;
											mnu2_2r = MN[in] - mu2_2r;

											for (icg = 0; icg < g_hat_num; icg++)
											{
												covs[0] = gg_hats_cov[icg];
												covs[1] = gg_hats[icg];
												covs[2] = gg_hats[icg];
												covs[3] = gg_hats_cov[icg];
												rand_multi_gauss(covs, mus, 2, corre_gs);

												mg11 = mg1_1r - corre_gs[0] * mnu1_1r;
												mg12 = mg1_2r - corre_gs[1] * mnu1_2r;

												gg_hist_label = histogram2d_s(mg11, mg12, mg_bins, mg_bins, mg_bin_num, mg_bin_num);
												//if (gg_hist_label >= chi_block_size)
												//{
												//	std::cout << "CROSS " << gg_hist_label << std::endl;
												//}
												chi_1[chi_label + chi_block_size * icg + gg_hist_label] += 1;

												//rand_multi_gauss(covs, mus, 2, corre_gs);

												mg21 = mg2_1r - corre_gs[0] * mnu2_1r;
												mg22 = mg2_2r - corre_gs[1] * mnu2_2r;
												gg_hist_label = histogram2d_s(mg21, mg22, mg_bins, mg_bins, mg_bin_num, mg_bin_num);
												//if (gg_hist_label >= chi_block_size)
												//{
												//	std::cout << "CROSS " << gg_hist_label << std::endl;
												//}
												chi_2[chi_label + chi_block_size * icg + gg_hist_label] += 1;
											}

										}
									}
								}
							}
						}

					}

				}

				t_block_e = clock();
				sprintf(log_inform, "RANK: %d finish %d 'th block in area %d in %.2f sec", rank, block_id, area_id, (t_block_e-t_block_s)/CLOCKS_PER_SEC);
				write_log(log_path, log_inform);
				if (0 == rank)
				{
					std::cout << log_inform << std::endl;
				}
				// free the GSL
				gsl_free();
			}

		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(chi_1, chi_block_size_ir * radius_bin_num, MPI_LONG, chi_1_total, chi_block_size_ir * radius_bin_num, MPI_LONG, 0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(chi_2, chi_block_size_ir * radius_bin_num, MPI_LONG, chi_2_total, chi_block_size_ir * radius_bin_num, MPI_LONG, 0, MPI_COMM_WORLD);
		for (k = 0; k < numprocs; k++)
		{
			if (k == rank)
			{
				std::cout << rank << " SYNC finish" << std::endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		sprintf(chi_path, "%scache/w_%d_chi_1_%d.hdf5", data_path, area_id, rank);
		sprintf(chi_set_name, "/data");
		//creat_h5_group(chi_path, set_name, TRUE);
		write_h5(chi_path, chi_set_name,  chi_1, mg_bin_num, mg_bin_num * g_hat_num * radius_bin_num);

		for (k = 0; k < numprocs; k++)
		{
			if (k == rank)
			{
				std::cout << rank << " chi_1 finish" << std::endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		sprintf(chi_set_name, "/data");		
		//creat_h5_group(chi_path, set_name, TRUE);
		sprintf(chi_path, "%scache/w_%d_chi_2_%d.hdf5", data_path, area_id, rank);
		write_h5(chi_path, chi_set_name,  chi_2, mg_bin_num, mg_bin_num * g_hat_num * radius_bin_num);

		for (k = 0; k < numprocs; k++)
		{
			if (k == rank)
			{
				std::cout << rank << " chi_2 finish" << std::endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		t_area_e = clock();
		sprintf(log_inform, "RANK: %d  finish the aread %d in %.2f sec", rank, area_id,(t_area_e- t_area_s)/CLOCKS_PER_SEC);
		write_log(log_path, log_inform);

		if (0 == rank)
		{
			std::cout << log_inform << std::endl;
			sprintf(chi_path, "%scache/total_chi_1_%d.hdf5", data_path, area_id);
			write_h5(chi_path, chi_set_name,  chi_1_total, mg_bin_num, mg_bin_num * g_hat_num * radius_bin_num*numprocs);

			sprintf(chi_path, "%scache/total_chi_2_%d.hdf5", data_path, area_id);
			write_h5(chi_path, chi_set_name,  chi_2_total, mg_bin_num, mg_bin_num * g_hat_num * radius_bin_num*numprocs);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		delete[] boundy;
		delete[] boundx;
		delete[] mg_bins;
		delete[] chi_1;
		delete[] chi_2;
		delete[] chi_1_total;
		delete[] chi_2_total;
		delete[] mask;
		delete[] task_list;
		delete[] my_tasks;
		delete[] search_blocks;
		delete[] num_in_block;
		delete[] block_start;
		delete[] block_end;

		delete[] RA ;
		delete[] DEC ;

		delete[] MG1;
		delete[] MG2;
		delete[] MN ;
		delete[] MU;
		delete[] MV ;
		//MPI_Win_free(&win_mg1);
		//MPI_Win_free(&win_mg2);
		//MPI_Win_free(&win_mn);
		//MPI_Win_free(&win_mu);
		//MPI_Win_free(&win_mv);
		//MPI_Win_free(&win_ra);
		//MPI_Win_free(&win_dec);
	}

	delete[] radius;
	delete[] radius_sq;
	delete[] gg_hats;
	delete[] gg_hats_cov;

	MPI_Finalize();
	return 0;
}
