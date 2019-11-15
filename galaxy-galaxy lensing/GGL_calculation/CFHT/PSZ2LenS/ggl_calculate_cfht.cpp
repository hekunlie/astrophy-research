#include<FQlib.h>
#include<mpi.h>
#include<vector>
#include<hk_iolib.h>

#define MAX_DATA_COL 40
#define FOREGAL_DATA_COL 6
#define LENS_FIT_COL 5
#define BACKGAL_DATA_COL 14
#define VEC_DATA_COL 8
#define MAX_PAIR 20000000
#define SMALL_CATA

#ifdef SMALL_CATA
#define MY_INT int
#else
#define MY_INT long
#endif // SMALL_CATA

int main(int argc, char *argv[])
{
	/* the macros should be adjusted in each case depends on how many columns will be read */
	/* The input paremeters:																									*/
	/*		1. the sky area label																								*/
	/*		2. the radius bin label, the search radius																			*/
	/*     3. the name of the foreground data set																			    */

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	// if the memory of the system could contain all the pairs
	// it will gather all the data and estimate the signal with SYM-PDF method
#if defined (SMALL_CATA)
	std::vector<double> data_cache;
	int data_col[1];
#endif

	int i, j, k, temp;
	char parent_path[250], data_path[250], result_path[250], log_path[250], h5f_path_grid[250], h5f_path_fore[250], h5f_res_path[250], temp_path[300];
	char set_name[50], set_name_2[50], attrs_name[80], log_infom[300];
	char foreground_name[50];

	// the controllers
	int area_id = atoi(argv[1]);
	int radius_label = atoi(argv[2]);
	strcpy(foreground_name, argv[3]);


	sprintf(parent_path, "/mnt/perc/hklee/CFHT/gg_lensing/");
	sprintf(data_path, "%sdata/", parent_path);
	sprintf(result_path, "%sresult/", parent_path);

	sprintf(h5f_path_grid, "%scfht_cata/cfht_cata_cut.hdf5", data_path);

	sprintf(h5f_res_path, "%s%s/cfht/w_%d/radius_%d.hdf5", result_path, foreground_name, area_id, radius_label);

	sprintf(log_path, "%slog/ggl_log_%d.dat", parent_path, rank); 

	sprintf(h5f_path_fore, "%sforeground/%s/w_%d_3.hdf5", data_path, foreground_name, area_id);
	
	//sprintf(parent_path, "/mnt/perc/hklee/CFHT/gg_lensing/");
	//sprintf(data_path, "%sdata/", parent_path);
	//sprintf(h5f_res_path, "%sresult/%s/fourier/w_%d/radius_%d.hdf5", parent_path, foreground_name, area_id, radius_label);
	//sprintf(log_path, "%slog/ggl_log_%d.dat", parent_path, rank);
	//sprintf(h5f_path_fore, "%sdata/foreground/%s/w_%d.hdf5", parent_path, foreground_name, area_id);

	pts_info gal_info;

	double st_start, st_end, st1, st2, st3, st4, stgal;
	int process_per;
	int per_n;

	int foregal_num;
	int my_gal_s, my_gal_e, gal_id;
	double *foregal_data[MAX_DATA_COL];
	long my_pair_count;
	int total_pair_count;
	double z_f, ra_f, dec_f;
	double dist_len, dist_integ_len, dist_source, dist_source_integ, dist_len_coeff, dist_integ_len_coeff;
	double rad_to_deg, deg_to_rad, coeff_rad_dist;
	double crit_surf_density_com, crit_surf_density_com_integ;


	int backgal_num;
	int *num_in_block, *block_start, *gal_in_block;
	double *ra_bin, *dec_bin, *block_boundx, *block_boundy;
	double *backgal_data[MAX_DATA_COL];
	double backgal_cos_2phi, backgal_sin_2phi, backgal_cos_4phi, backgal_sin_4phi;
	double backgal_e_t, backgal_e_x, backgal_m, backgal_c2;
	double z_b, z_thresh, z_b_sig95, z_b_odds;
	double ra_b, dec_b;
	double diff_ra, diff_dec, diff_r, diff_theta, diff_theta_sq, diff_z_thresh;

	int subset_num[1];
	double *my_data_buf, *final_buf;
	// the buf
	double my_signal_buf[2];
	// the signal from the all areas
	double *delta_sigma;
	int *pair_count_shared;



	int grid_num, grid_ny, grid_nx;
	int row, col, bin_label;
	int block_id, block_s, block_e, ib, ibb, ibg;
	double block_scale[1];
	int *block_mask;
	int ra_bin_num, dec_bin_num;


	int radius_num;
	double radius_s, radius_e, radius_e_sq;
	double *radius_bin;

	// source Z > len Z + diff_z_thresh
	diff_z_thresh = 0.05;


	///////////////////////////////////// radius bin ////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	radius_num = 24;
	radius_bin = new double[radius_num + 1]{};
	log_bin(0.1, 25.12, radius_num + 1, radius_bin);
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////

	int e1_id = 0, e2_id = 1, weight_id = 2, m_id = 3, c_id = 4;
	int z_id = 5, dist_id = 6, dist_integ_id=7, ra_id = 8, dec_id = 9, cos_dec_id = 10;
	int zmin_lb = 11, zmax_lb = 12, odds_lb = 13;

	int shape[2];

	char *names[MAX_DATA_COL];
	for (i = 0; i < BACKGAL_DATA_COL; i++)
	{
		names[i] = new char[40];
	}

	sprintf(names[e1_id], "E1");
	sprintf(names[e2_id], "E2");
	sprintf(names[weight_id], "WEIGHT");
	sprintf(names[m_id], "M");
	sprintf(names[c_id], "C");

	sprintf(names[z_id], "Z");
	sprintf(names[dist_id], "DISTANCE");
	sprintf(names[dist_integ_id], "DISTANCE_INTEG");
	sprintf(names[ra_id], "RA");
	sprintf(names[dec_id], "DEC");
	sprintf(names[cos_dec_id], "COS_DEC");

	sprintf(names[zmin_lb], "Z_MIN");
	sprintf(names[zmax_lb], "Z_MAX");
	sprintf(names[odds_lb], "ODDS");


	sprintf(log_infom, "RANK: %d. Start area: w_%d, radius: %d", rank, area_id, radius_label);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// the shared buffer for the total \Delta\Sigma_t , \Delta\Sigma_x, and the errors
#if ! defined(SMALL_CATA)

	MPI_Win win_delta_sigma;
	MPI_Aint size_delta_sigma;
	// the signals from each field could be added into the final one to get the signal of the whole area
	// so, the denominator and numerator are stored separately
	//[\Sum weight*\Sigma_t,  \Sum weight*\Sigma_t_err,  
	//  \Sum weight*\Sigma_x, \Sum weight*\Sigma_x_err,  \Sum weight]
	size_delta_sigma = radius_num*5;

#endif
	MPI_Win win_pair_count;
	MPI_Aint size_pair_count;

	size_pair_count = numprocs;

	if (0 == rank)
	{
#if ! defined( SMALL_CATA)
		// for the chi square of tangential shear
		MPI_Win_allocate_shared(size_delta_sigma * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &delta_sigma, &win_delta_sigma);
#endif
		MPI_Win_allocate_shared(size_pair_count * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &pair_count_shared, &win_pair_count);
	}
	else
	{
		int dispu_total;
#if ! defined( SMALL_CATA)
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &delta_sigma, &win_delta_sigma);
		MPI_Win_shared_query(win_delta_sigma, 0, &size_delta_sigma, &dispu_total, &delta_sigma);

#endif
		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &pair_count_shared, &win_pair_count);
		MPI_Win_shared_query(win_pair_count, 0, &size_pair_count, &dispu_total, &pair_count_shared);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// initialization
	if (0 == rank)
	{
#if ! defined( SMALL_CATA)
		initialize_arr(delta_sigma, radius_num * 5, 0);
#endif
		initialize_arr(pair_count_shared, numprocs, 0);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// read foreground information
	// Z, DISTANCE, DISTANCE_INTEG, RA, DEC, COS_DEC
	for (i = LENS_FIT_COL; i < FOREGAL_DATA_COL + LENS_FIT_COL; i++)
	{
		sprintf(set_name, "/%s", names[i]);
		read_h5_datasize(h5f_path_fore, set_name, foregal_num);

		foregal_data[i] = new double[foregal_num];
		read_h5(h5f_path_fore, set_name, foregal_data[i]);
	}
	sprintf(log_infom, "RANK: %d. w_%d. Read foreground data. %d galaxies", rank, area_id, foregal_num);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////


	////////////////////////////  read background information ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// E1, E2, WEIGHT, M, C,  Z,  DISTANCE, DISTANCE_INTEG,  RA, DEC,  COS_DEC,  Z_MIN, Z_MAX, ODDS, 
	for (i = 0; i < BACKGAL_DATA_COL; i++)
	{
		sprintf(set_name, "/w_%d/%s", area_id, names[i]);
		read_h5_datasize(h5f_path_grid, set_name, backgal_num);
		backgal_data[i] = new double[backgal_num] {};
		read_h5(h5f_path_grid, set_name, backgal_data[i]);
	}

	sprintf(set_name, "/w_%d/grid_shape", area_id);
	read_h5(h5f_path_grid, set_name, shape);
	grid_ny = shape[0];
	grid_nx = shape[1];
	grid_num = grid_ny * grid_nx;
	sprintf(set_name, "/block_scale");
	read_h5(h5f_path_grid, set_name, block_scale); // degree

	block_mask = new int[grid_num] {};

	num_in_block = new int[grid_num] {};
	block_start = new int[grid_num] {};
	gal_in_block = new int[backgal_num] {};

	block_boundx = new double[grid_num * 4];
	block_boundy = new double[grid_num * 4];

	sprintf(set_name, "/w_%d/num_in_block", area_id);
	read_h5(h5f_path_grid, set_name, num_in_block);

	sprintf(set_name, "/w_%d/block_start", area_id);
	read_h5(h5f_path_grid, set_name, block_start);

	sprintf(set_name, "/w_%d/gal_in_block", area_id);
	read_h5(h5f_path_grid, set_name, gal_in_block);

	sprintf(set_name, "/w_%d/RA_bin", area_id);
	read_h5_datasize(h5f_path_grid, set_name, ra_bin_num);
	ra_bin = new double[ra_bin_num];
	ra_bin_num = ra_bin_num - 1;
	read_h5(h5f_path_grid, set_name, ra_bin);

	sprintf(set_name, "/w_%d/DEC_bin", area_id);
	read_h5_datasize(h5f_path_grid, set_name, dec_bin_num);
	dec_bin = new double[dec_bin_num];
	dec_bin_num = dec_bin_num - 1;
	read_h5(h5f_path_grid, set_name, dec_bin);


	sprintf(set_name, "/w_%d/block_boundx", area_id);
	read_h5(h5f_path_grid, set_name, block_boundx);

	sprintf(set_name, "/w_%d/block_boundy", area_id);
	read_h5(h5f_path_grid, set_name, block_boundy);

	sprintf(log_infom, "RANK: %d. w_%d. Read background data. %d galaxies. %d grids (%d x %d)", rank, area_id, backgal_num, grid_num, grid_ny, grid_nx);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// FOR DEBUG THE DETECTION LABEL
	// int *detect_mask = new int[backgal_num];
	// char mask_path[200];
	// for(i=0;i<backgal_num;i++)
	// {
	// 	detect_mask[i] = -1;
	// }
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	st1 = clock();
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

	sprintf(log_infom, "RANK: %d. w_%d. Task distribution: my gal: %d ~ %d. ", rank, area_id, my_gal_s, my_gal_e);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
	}

	rad_to_deg = 180. / Pi;
	deg_to_rad = Pi / 180;
	my_pair_count = 0;

	st1 = clock();
	for (gal_id = my_gal_s; gal_id < my_gal_e; gal_id++)
	{
		z_f = foregal_data[z_id][gal_id];
		// the source must be at z = z_f + diff_z_thresh
		z_thresh = z_f + diff_z_thresh;

		// comoving distance
		dist_len = foregal_data[dist_id][gal_id];
		dist_len_coeff = 1. / dist_len / (1 + z_f);
		// the integrate part of the comoving distance
		dist_integ_len = foregal_data[dist_integ_id][gal_id];
		dist_integ_len_coeff = 1. / dist_integ_len / (1 + z_f);

		// the max searching radius depend on the redshift of lens  //	
		radius_e = radius_bin[radius_label + 1] *rad_to_deg/ dist_len / foregal_data[cos_dec_id][gal_id] * 2; // degree

		// degree
		ra_f = foregal_data[ra_id][gal_id];
		dec_f = foregal_data[dec_id][gal_id];
		//std::cout<<ra_f<<" "<<dec_f<<std::endl;

		// all the data has been arranged into blocks, find the block label of this foreground galaxy
		histogram2d_s(dec_f, ra_f, dec_bin, ra_bin, dec_bin_num, ra_bin_num, bin_label);
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

		// find the blocks needed //
		initialize_arr(block_mask, grid_num, -1);
		find_block(&gal_info, radius_e, block_boundy, block_boundx, block_mask);

		
		for (ib = 0; ib < grid_num; ib++)
		{
			// the block id [i,j] = i*nx + j 
			block_id = block_mask[ib];
			// if(block_mask[ib] > -1)
			// {	
			// 	std::cout<<ib<<" "<<block_id<<" "<<block_id/22<<" "<<block_id%22<<std::endl;
				
			// }
			if (block_id > -1 and num_in_block[block_id] > 0)
			{
				// the start and end point of the block //
				// the start- & end-point					      //
				block_s = block_start[block_id];

				//std::cout << block_mask[block_id] << " " << block_s  << std::endl;
				for (ibb = 0; ibb < num_in_block[block_id]; ibb++)
				{	
					// the real data label in the array
					ibg = gal_in_block[block_s + ibb];
					z_b_sig95 = z_f + (backgal_data[zmax_lb][ibg] - backgal_data[zmin_lb][ibg]) / 2;
					z_b_odds = backgal_data[odds_lb][ibg];
					z_b = backgal_data[z_id][ibg];

					if (z_b >= z_thresh and z_b >= z_b_sig95)
					{
						ra_b = backgal_data[ra_id][ibg];
						dec_b = backgal_data[dec_id][ibg];

						// comoving distance
						dist_source = backgal_data[dist_id][ibg];
						// the integrate part of the comoving distance
						dist_source_integ = backgal_data[dist_integ_id][ibg];

						// times cos(dec) due to the different length the arc corresponding to the same delta R.A. at different Dec
						diff_ra = (ra_b - ra_f)*foregal_data[cos_dec_id][gal_id];
						diff_dec = dec_b - dec_f;
						diff_theta_sq = diff_ra * diff_ra + diff_dec * diff_dec; // degree^2

						// the seperation in comoving coordinate, 
						//diff_r = dist_len * sqrt(diff_theta_sq)*coeff_inv;

						separation(ra_b, dec_b, ra_f, dec_f, diff_theta);
						diff_r = dist_len * diff_theta;
						//std::cout << radius_bin[radius_label] << " " << diff_r <<" "<< diff_theta<< " "<< sqrt(diff_theta_sq)<<" " << radius_bin[radius_label + 1] << std::endl;

						if (radius_bin[radius_label] <= diff_r and diff_r < radius_bin[radius_label + 1])
						{	
							my_pair_count += 1;
							crit_surf_density_com = dist_source / (dist_source - dist_len) *dist_len_coeff;
							crit_surf_density_com_integ = dist_source_integ / (dist_source_integ - dist_integ_len) *dist_integ_len_coeff;

							// rotation for shear calculation, see the NOTE of gg_lensing for the details 
							backgal_sin_2phi = 2 * diff_ra*diff_dec / diff_theta_sq;
							backgal_cos_2phi = (diff_dec - diff_ra)*(diff_ra + diff_dec) / diff_theta_sq;

#if defined (SMALL_CATA)
							//// tangential component, e_t =  e_1 *cos2\phi - e_2*sin2\phi
							backgal_e_t = backgal_data[e1_id][ibg] * backgal_cos_2phi - backgal_data[e2_id][ibg] * backgal_sin_2phi;
							//// the cross component, e_x =  e_1 *sin2\phi + e_2*cos2\phi
							backgal_e_x = backgal_data[e1_id][ibg] * backgal_sin_2phi + backgal_data[e2_id][ibg] * backgal_cos_2phi;


							data_cache.push_back(backgal_e_t);
							data_cache.push_back(backgal_e_x);

							// multiplicative bias
							data_cache.push_back(backgal_data[m_id][ibg]);
							data_cache.push_back(backgal_data[weight_id][ibg]);

							data_cache.push_back(crit_surf_density_com_integ);
							data_cache.push_back(crit_surf_density_com);

							data_cache.push_back(diff_r);
							data_cache.push_back(z_b);
							// FOR DEBUG THE DETECTION LABEL
							// detect_mask[ibg] = radius_label;

#else

							// calculate the PDF of the estimator for shear
							for (ig = 0; ig < gh_num; ig++)
							{
								mg_t = backgal_mg_tan - gh[ig] * backgal_mnu1_tan_c;
								histogram_s(mg_t, mg1_bin, mg_bin_num, chi_bin_label);
								my_chi_tan[ig*mg_bin_num + chi_bin_label] += 1;

								mg_x = backgal_mg_cross - gh[ig] * backgal_mnu2_tan_c;
								histogram_s(mg_x, mg2_bin, mg_bin_num, chi_bin_label);
								my_chi_cross[ig*mg_bin_num + chi_bin_label] += 1;
							}

							// calculate the PDF of the estimator for 'shear*critical_surface_density'
							for (ig = 0; ig < gh_crit_num; ig++)
							{
								mg_t = backgal_mg_tan - gh_crit[ig] * backgal_mnu1_tan;
								histogram_s(mg_t, mg1_bin, mg_bin_num, chi_bin_label);
								my_chi_crit_tan[ig*mg_bin_num + chi_bin_label] += 1;

								mg_x = backgal_mg_cross - gh_crit[ig] * backgal_mnu2_tan;
								histogram_s(mg_x, mg2_bin, mg_bin_num, chi_bin_label);
								my_chi_crit_cross[ig*mg_bin_num + chi_bin_label] += 1;
							}
#endif
						}

					}
				}
			}
		}
		stgal = clock();
		if (my_gal_e - my_gal_s > 10)
		{
			process_per = (gal_id - my_gal_s) % (int((my_gal_e - my_gal_s) *0.1));
			per_n = int((gal_id - my_gal_s) *1.0 / (my_gal_e - my_gal_s) * 100);
		}
		else
		{
			per_n = gal_id;
		}
		if (0 == process_per)
		{
			sprintf(log_infom, "RANK: %d. w_%d. Looping foreground %d%%. Time since begin: %.2f sec.", rank, area_id, per_n, (stgal - st1) / CLOCKS_PER_SEC);
			if (0 == rank)
			{
				std::cout << log_infom << std::endl;
			}
		}
	}
	pair_count_shared[rank] = my_pair_count;
	MPI_Barrier(MPI_COMM_WORLD);

	// FOR DEBUG THE DETECTION LABEL
	// sprintf(mask_path,"/home/hklee/transfer/mask.hdf5");
	// sprintf(set_name, "/mask_%d",radius_label);
	//  if(radius_label == 0)
	//  {
	//  	write_h5(mask_path, set_name, detect_mask, backgal_num, 1, TRUE);
	//  }
	//  else
	//  {
	//  	write_h5(mask_path, set_name, detect_mask, backgal_num, 1, FALSE);
	//  }
	// std::cout<<backgal_num<<" galaxies"<<std::endl;

	sum_arr(pair_count_shared, numprocs, 0, numprocs, total_pair_count);
	if (1.0*my_pair_count*VEC_DATA_COL >= INT_MAX or 1.0*total_pair_count * VEC_DATA_COL >= INT_MAX)
	{
		std::cout << "INT overflow. Too many pairs." << std::endl;
		exit(0);
	}
	MPI_Barrier(MPI_COMM_WORLD);

#if defined(SMALL_CATA)

	if (rank == 0)
	{
		sprintf(set_name, "/pair_count");
		write_h5(h5f_res_path, set_name, pair_count_shared, numprocs, 1, TRUE);
		sprintf(set_name, "/radius_bin");
		write_h5(h5f_res_path, set_name, radius_bin, radius_num + 1, 1, FALSE);

		data_col[0] = VEC_DATA_COL;
		sprintf(set_name, "/data_col");
		write_h5(h5f_res_path, set_name, data_col, 1, 1, FALSE);

		if (total_pair_count < 1)
		{
			subset_num[0] = 0;
			sprintf(set_name, "/subset_num");
			write_h5(h5f_res_path, set_name, subset_num, 1, 1, FALSE);
		}
		sprintf(log_infom, "RANK: %d. w_%d. %d galaxies have been found in Radius [%.4f, %.4f].", rank, area_id, total_pair_count, radius_bin[radius_label], radius_bin[radius_label + 1]);
		std::cout << log_infom << std::endl;
	}

	if (total_pair_count > 1)
	{
		// final_buf will store the data of all the pairs
		my_data_buf = new double[pair_count_shared[rank] * VEC_DATA_COL]{};
		// copy the data in the vector into the buffer 
		if (!data_cache.empty())
		{
			memcpy(my_data_buf, &data_cache[0], data_cache.size() * sizeof(double));
		}

		if (numprocs > 1)
		{
			if (total_pair_count <= MAX_PAIR)
			{
				// calculate the entry of each rank in the big buffer
				int *displ = new int[numprocs] {};
				int *num_of_thread = new int[numprocs] {};

				for (i = 0; i < numprocs; i++)
				{
					num_of_thread[i] = pair_count_shared[i] * VEC_DATA_COL;
					for (j = 0; j < i; j++)
					{
						displ[i] += num_of_thread[j];
					}
				}
				if (rank == 0)
				{
					final_buf = new double[total_pair_count * VEC_DATA_COL];
					//show_arr(displ, 1, numprocs);
					//show_arr(num_of_thread, 1, numprocs);
				}
				MPI_Barrier(MPI_COMM_WORLD);

				MPI_Gatherv(my_data_buf, num_of_thread[rank], MPI_DOUBLE, final_buf, num_of_thread, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				MPI_Barrier(MPI_COMM_WORLD);

				if (rank == 0)
				{
					// only one subset
					subset_num[0] = 1;
					sprintf(set_name, "/subset_num");
					write_h5(h5f_res_path, set_name, subset_num, 1, 1, FALSE);

					sprintf(set_name, "/pair_data_0");
					write_h5(h5f_res_path, set_name, final_buf, total_pair_count, VEC_DATA_COL, FALSE);
					delete[] final_buf;
				}
			}
			else
			{
				if (rank == 0)
				{
					subset_num[0] = numprocs;
					sprintf(set_name, "/subset_num");
					write_h5(h5f_res_path, set_name, subset_num, 1, 1, FALSE);
				}
				for (i = 0; i < numprocs; i++)
				{
					if (i == rank)
					{
						if (pair_count_shared[rank] > 0)
						{
							sprintf(set_name, "/pair_data_%d", rank);
							write_h5(h5f_res_path, set_name, my_data_buf, pair_count_shared[rank], VEC_DATA_COL, FALSE);
						}
					}
					MPI_Barrier(MPI_COMM_WORLD);
				}
			}
		}
		else
		{
			subset_num[0] = 1;
			sprintf(set_name, "/subset_num");
			write_h5(h5f_res_path, set_name, subset_num, 1, 1, FALSE);

			sprintf(set_name, "/pair_data_0");
			write_h5(h5f_res_path, set_name, my_data_buf, total_pair_count, VEC_DATA_COL, FALSE);
		}
		delete[] my_data_buf;
		MPI_Barrier(MPI_COMM_WORLD);
	}
#else
	if (rank == 0)
	{
		sprintf(log_infom, "RANK: %d. w_%d. %ld galaxies have been found in Radius [%.4f, %.4f].", rank, area_id, pair_count, radius_bin[radius_label], radius_bin[radius_label + 1]);
		std::cout << log_infom << std::endl;
	}
	// merge the PDF of signal estimator into the shared buffer
	for (i = 0; i < numprocs; i++)
	{
		if (rank == i)
		{
			for (j = 0; j < mg_bin_num*gh_num; j++)
			{
				chi_tan_shared[j] += my_chi_tan[j];
				chi_cross_shared[j] += my_chi_cross[j];
			}
			for (j = 0; j < mg_bin_num*gh_crit_num; j++)
			{
				chi_crit_tan_shared[j] += my_chi_crit_tan[j];
				chi_crit_cross_shared[j] += my_chi_crit_cross[j];
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	

	// free the memory	
	MPI_Win_free(&win_delta_sigma);
#endif

	MPI_Win_free(&win_pair_count);

	sprintf(log_infom, "RANK: %d. Write file of Radius bin [%.4f,  %.4f] ---- end. ", rank, radius_bin[radius_label], radius_bin[radius_label + 1]);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		char times[50];
		get_time(times, 50);
		std::cout << log_infom << std::endl;
		std::cout << h5f_res_path << std::endl;
		std::cout << times << std::endl;
	}


	for (i = 0; i < BACKGAL_DATA_COL; i++)
	{
		delete[] backgal_data[i];
	}
	for (i = LENS_FIT_COL; i < FOREGAL_DATA_COL + LENS_FIT_COL; i++)
	{
		delete[] foregal_data[i];
	}

	delete[] num_in_block;
	delete[] block_boundx;
	delete[] block_boundy;
	delete[] gal_in_block;
	delete[] block_start;

	delete[] ra_bin;
	delete[] dec_bin;
	delete[] block_mask;

	delete[] radius_bin;


	MPI_Finalize();
	return 0;
}
