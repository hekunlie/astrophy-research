#include<FQlib.h>
#include<mpi.h>
#include<vector>

#define max_data_col 40
#define foregal_data_col 5
#define grid_data_col 5
#define backgal_data_col 22

#define SMALL_CATA

#ifdef SMALL_CATA
#define MY_INT int
#else
#define MY_INT long
#endif // SMALL_CATA

int main(int argc, char *argv[])
{
	/* the macros should be adjusted in each case depends on how many columns will be read */
	/* The input paremeters:																												*/
	/*		1. the sky area label																												*/
	/*		2. the radius bin label, the search raidus																			*/
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
	int vec_data_col = 9;
#endif

	int i, j, k, temp;
	char parent_path[250], data_path[250], log_path[250], h5f_path_grid[250], h5f_path_fore[250], h5f_res_path[250], temp_path[300];
	char set_name[50], set_name_2[50], attrs_name[80], log_infom[300];
	char foreground_name[50];

	// the controllers
	int radius_label = atoi(argv[1]);
	int area_id = 1;

	sprintf(parent_path, "/home/hklee/work/cpp/code_test/separation_angle_test/");
	sprintf(data_path, "%sdata/", parent_path);
	sprintf(h5f_res_path, "%sradius_%d.hdf5", data_path, radius_label);
	sprintf(log_path, "%slog/ggl_log_%d.dat", data_path, rank);
	sprintf(h5f_path_fore, "%s/cmass_w1_sub.hdf5", data_path);

	sprintf(h5f_path_grid, "/mnt/perc/hklee/CFHT/gg_lensing/data/cfht_cata_grid.hdf5");

	pts_info gal_info;

	double st_start, st_end, st1, st2, st3, st4, stgal;
	int process_per;
	int per_n;

	int foregal_num;
	int my_gal_s, my_gal_e, gal_id;
	double *foregal_data[max_data_col];
	MY_INT pair_count;// be carefull, the pair number may be too many, long or double 
	double z_f, ra_f, dec_f;
	double dist_len, dist_source, dist_len_coeff;
	double coeff, coeff_inv, coeff_rad_dist;
	double crit_surf_density_com;


	int backgal_num;
	double *backgal_data[max_data_col];
	double backgal_cos_2phi, backgal_sin_2phi, backgal_cos_4phi, backgal_sin_4phi;
	double backgal_e_t, backgal_e_x, backgal_m, backgal_c2;
	double z_b, z_thresh, z_b_sig95, z_b_odds;
	double ra_b, dec_b;
	double diff_ra, diff_dec, diff_r, diff_theta, diff_theta_sq, diff_z_thresh = 0.1;

	double *my_data_buf, *final_buf;
	// the buf
	double my_signal_buf[2];
	// the signal from the all areas
	double *delta_sigma;
	MY_INT *pair_count_shared;


	int grid_num, grid_ny, grid_nx;
	int row, col, bin_label;
	int block_id, block_s, block_e, ib;
	double block_scale[1];
	int *block_mask;
	int ra_bin_num, dec_bin_num;


	int radius_num;
	double radius_s, radius_e, radius_e_sq;
	double *radius_bin;
	// radius bin
	radius_num = 20;
	radius_bin = new double[radius_num + 1]{};
	log_bin(0.01, 12, radius_num + 1, radius_bin);


	int nib_id = 0, bs_id = 1, be_id = 2, bdy_id = 3, bdx_id = 4;
	int z_id = 5, dist_id = 6, ra_id = 7, dec_id = 8, cos_dec_id = 9;
	int e1_id = 10, e2_id = 11, weight_id = 12, m_id = 13, c_id = 14;
	int starflag_id = 15, zmin_lb = 16, zmax_lb = 17, odds_lb = 18, mag_lb = 19;
	int ra_bin_id = 20, dec_bin_id = 21;
	int block_scale_id = 22, grid_shape_id = 23;

	int shape[2];


	char *names[max_data_col];//backgal_data_col
	for (i = 0; i < backgal_data_col + 2; i++)
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

	sprintf(names[ra_bin_id], "RA_bin");
	sprintf(names[dec_bin_id], "DEC_bin");

	sprintf(names[block_scale_id], "block_scale");
	sprintf(names[grid_shape_id], "grid_shape");

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
	size_delta_sigma = radius_num * 5;

#endif
	MPI_Win win_pair_count;
	MPI_Aint  size_pair_count;

	size_pair_count = numprocs;

	if (0 == rank)
	{
#if ! defined( SMALL_CATA)
		// for the chi square of tangential shear
		MPI_Win_allocate_shared(size_delta_sigma * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &delta_sigma, &win_delta_sigma);
#endif
		MPI_Win_allocate_shared(size_pair_count * sizeof(MY_INT), sizeof(MY_INT), MPI_INFO_NULL, MPI_COMM_WORLD, &pair_count_shared, &win_pair_count);
	}
	else
	{
		int dispu_total;
#if ! defined( SMALL_CATA)
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &delta_sigma, &win_delta_sigma);
		MPI_Win_shared_query(win_delta_sigma, 0, &size_delta_sigma, &dispu_total, &delta_sigma);

#endif
		MPI_Win_allocate_shared(0, sizeof(MY_INT), MPI_INFO_NULL, MPI_COMM_WORLD, &pair_count_shared, &win_pair_count);
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


	// read the search radius
	//sprintf(set_name, "/radius_bin");
	//read_h5_datasize(h5f_path_grid, set_name, radius_num);
	//radius_bin = new double[radius_num] {};
	//read_h5(h5f_path_grid, set_name, radius_bin);
	//radius_num = radius_num - 1;


	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// read foreground information
	// Z, DISTANCE, RA, DEC, COS_DEC
	for (i = grid_data_col; i < foregal_data_col + grid_data_col; i++)
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


	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// read background information
	sprintf(set_name, "/w_%d/grid_shape", area_id);
	read_h5(h5f_path_grid, set_name, shape);
	grid_ny = shape[0];
	grid_nx = shape[1];
	grid_num = grid_ny * grid_nx;
	sprintf(set_name, "/block_scale");
	read_h5(h5f_path_grid, set_name, block_scale); // degree

	block_mask = new int[grid_num] {};

	// Z, Z_MIN, Z_MAX, ODDS,  DISTANCE, RA, DEC,  E1, E2, WEIGHT,  num_in_block,  block_start, block_end, 
	// block_boundx, block_boundy
	for (i = 0; i < backgal_data_col; i++)
	{
		sprintf(set_name, "/w_%d/%s", area_id, names[i]);
		read_h5_datasize(h5f_path_grid, set_name, temp);

		if (z_id == i)
		{
			backgal_num = temp;
		}
		if (ra_bin_id == i)
		{
			ra_bin_num = temp - 1;
		}
		if (dec_bin_id == i)
		{
			dec_bin_num = temp - 1;
		}
		backgal_data[i] = new double[temp] {};
		read_h5(h5f_path_grid, set_name, backgal_data[i]);
	}

	sprintf(log_infom, "RANK: %d. w_%d. Read background data. %d galaxies. %d grids (%d x %d)", rank, area_id, backgal_num, grid_num, grid_ny, grid_nx);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////


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

	coeff = 0.18 / C_0_hat / Pi;
	coeff_inv = C_0_hat * Pi / 0.18;
	coeff_rad_dist = C_0_hat * 1000;

	radius_e = radius_bin[radius_label + 1] * coeff;

	st1 = clock();
	for (gal_id = my_gal_s; gal_id < my_gal_e; gal_id++)
	{
		z_f = foregal_data[z_id][gal_id];
		// the source must be at z = z_f + diff_z_thresh
		z_thresh = z_f + diff_z_thresh;

		//find_near(redshifts, z_f, red_num, tag_f);
		//dist_len = distances[tag_f];
		dist_len = foregal_data[dist_id][gal_id];
		dist_len_coeff = 1. / dist_len / (1 + z_f);

		// the max searching radius depend on the redshift of lens  //	
		// the max seperation of the source at the z = z_f  //
		radius_e = radius_bin[radius_label + 1] * coeff / dist_len / foregal_data[cos_dec_id][gal_id] * 2; // degree
		radius_e_sq = radius_e * radius_e; // degree^2

		// degree
		ra_f = foregal_data[ra_id][gal_id];
		dec_f = foregal_data[dec_id][gal_id];

		// all the data has been aranged into blocks, find the block label of this foreground galaxy
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

		// find the blocks needed //
		initialize_arr(block_mask, grid_num, -1);
		find_block(&gal_info, radius_e, backgal_data[bdy_id], backgal_data[bdx_id], block_mask);

		for (block_id = 0; block_id < grid_num; block_id++)
		{
			if (block_mask[block_id] > -1)
			{
				// the start and end point of the block //
				// the start- & end-point					      //
				block_s = backgal_data[bs_id][block_mask[block_id]];
				block_e = backgal_data[be_id][block_mask[block_id]];

				//std::cout << block_mask[block_id] << " " << block_s << " " << block_e << std::endl;
				for (ib = block_s; ib < block_e; ib++)
				{
					z_b_sig95 = z_f + (backgal_data[zmax_lb][ib] - backgal_data[zmin_lb][ib]) / 2;
					z_b_odds = backgal_data[odds_lb][ib];

					//if (backgal_data[z_id][ib] >= z_thresh)
					if (backgal_data[z_id][ib] >= z_thresh and backgal_data[z_id][ib] > z_b_sig95)
					{
						ra_b = backgal_data[ra_id][ib];
						dec_b = backgal_data[dec_id][ib];

						dist_source = backgal_data[dist_id][ib];

						// times cos(dec) due to the different length the arc corresponding to the same delta R.A. at different Dec
						diff_ra = (ra_b - ra_f)*foregal_data[cos_dec_id][gal_id];
						diff_dec = dec_b - dec_f;
						diff_theta_sq = diff_ra * diff_ra + diff_dec * diff_dec; // degree^2

						// the seperation in comving coordinate, 
						//diff_r = dist_source * sqrt(diff_theta_sq)*coeff_inv;

						separation(ra_b, dec_b, ra_f, dec_f, diff_theta);
						diff_r = dist_source * sin(diff_theta) * coeff_rad_dist;

						if (diff_r >= radius_bin[radius_label] and diff_r < radius_bin[radius_label + 1])
						{
							pair_count_shared[rank] += 1;
							crit_surf_density_com = dist_source / (dist_source - dist_len) *dist_len_coeff;

							// rotation for shear calculation, see the NOTE of gg_lensing for the detials 
							backgal_sin_2phi = 2 * diff_ra*diff_dec / diff_theta_sq;
							backgal_cos_2phi = (diff_dec - diff_ra)*(diff_ra + diff_dec) / diff_theta_sq;

#if defined (SMALL_CATA)

							//// the direction of R.A. is oppsite, actually,  e2 = -e2
							//// tangential component, e_t =  e_1 *cos2\phi - e_2*sin2\phi
							//backgal_e_t = backgal_data[e1_id][ib] * backgal_cos_2phi + backgal_data[e2_id][ib] * backgal_sin_2phi;
							//data_cache.push_back(backgal_e_t);

							//// the cross component, e_x =  e_1 *sin2\phi + e_2*cos2\phi
							//backgal_e_x = backgal_data[e1_id][ib] * backgal_sin_2phi - backgal_data[e2_id][ib] * backgal_cos_2phi;
							//data_cache.push_back(backgal_e_x);

							data_cache.push_back(backgal_data[e1_id][ib]);
							data_cache.push_back(backgal_data[e2_id][ib]);
							data_cache.push_back(backgal_cos_2phi);
							data_cache.push_back(backgal_sin_2phi);
							// mutiplicative bias
							data_cache.push_back(backgal_data[m_id][ib]);
							// additive bias
							data_cache.push_back(backgal_data[c_id][ib]);
							data_cache.push_back(backgal_data[weight_id][ib]);
							data_cache.push_back(crit_surf_density_com);
							data_cache.push_back(diff_r);

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
		process_per = (gal_id - my_gal_s) % (int((my_gal_e - my_gal_s) *0.1));
		per_n = int((gal_id - my_gal_s) *1.0 / (my_gal_e - my_gal_s) * 100);

		if (0 == process_per)
		{
			sprintf(log_infom, "RANK: %d. w_%d. Looping foreground %d%%. Time since begin: %.2f sec.", rank, area_id, per_n, (stgal - st1) / CLOCKS_PER_SEC);
			if (0 == rank)
			{
				std::cout << log_infom << std::endl;
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	sum_arr(pair_count_shared, numprocs, 0, numprocs, pair_count);

	MPI_Barrier(MPI_COMM_WORLD);

#if defined(SMALL_CATA)
	if (rank == 0)
	{
		sprintf(log_infom, "RANK: %d. w_%d. %d galaxies have been found in Radius [%.4f, %.4f].", rank, area_id, pair_count, radius_bin[radius_label], radius_bin[radius_label + 1]);
		std::cout << log_infom << std::endl;
	}
	if (pair_count > 1)
	{
		// final_buf will store the data of all the pairs
		my_data_buf = new double[pair_count_shared[rank] * vec_data_col]{};
		// copy the data in the vector into the buffer 
		if (!data_cache.empty())
		{
			memcpy(my_data_buf, &data_cache[0], data_cache.size() * sizeof(double));
		}

		// if more than 2 cpus
		if (numprocs > 1)
		{
			// calculate the entry of each rank in the big buffer
			MY_INT *displ = new MY_INT[numprocs]{};
			MY_INT *num_of_thread = new MY_INT[numprocs]{};

			for (i = 0; i < numprocs; i++)
			{
				num_of_thread[i] = pair_count_shared[i] * vec_data_col;
				for (j = 0; j < i; j++)
				{
					displ[i] += num_of_thread[j];
				}
			}
			if (rank == 0)
			{
				final_buf = new double[pair_count * vec_data_col];
				//show_arr(displ, 1, numprocs);
				//show_arr(num_of_thread, 1, numprocs);
			}
			MPI_Barrier(MPI_COMM_WORLD);

			//char test_path[200];
			//sprintf(test_path, "/home/hkli/work/test/%d.hdf5", rank);
			//sprintf(set_name, "/pair_data");
			//write_h5(test_path, set_name, my_data_buf, pair_count_shared[rank], 5, TRUE);
			// gather the data from each thread, empty data from some threads are allowed

			MPI_Gatherv(my_data_buf, num_of_thread[rank], MPI_DOUBLE, final_buf, num_of_thread, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);

			if (rank == 0)
			{
				//show_arr(final_buf, pair_count, vec_data_col);
				sprintf(set_name, "/pair_data");
				write_h5(h5f_res_path, set_name, final_buf, pair_count, vec_data_col, TRUE);

				sprintf(temp_path, "%sresult/%s/cfht/w_%d/radius_bin.hdf5", parent_path, foreground_name, area_id);
				sprintf(set_name, "/radius_bin");
				write_h5(temp_path, set_name, radius_bin, radius_num + 1, 1, TRUE);

				delete[] final_buf;
			}
		}
		else
		{
			// only one cpu
			sprintf(set_name, "/pair_data");
			write_h5(h5f_res_path, set_name, my_data_buf, pair_count, vec_data_col, TRUE);

			sprintf(temp_path, "%sresult/%s/cfht/w_%d/radius_bin.hdf5", parent_path, foreground_name, area_id);
			sprintf(set_name, "/radius_bin");
			write_h5(temp_path, set_name, radius_bin, radius_num + 1, 1, TRUE);
		}
		delete[] my_data_buf;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else
	{
		if (rank == 0)
		{
			my_data_buf = new double[vec_data_col];
			initialize_arr(my_data_buf, vec_data_col, -1);
			sprintf(set_name, "/pair_data");
			write_h5(h5f_res_path, set_name, my_data_buf, 1, vec_data_col, TRUE);

			sprintf(temp_path, "%sresult/%s/cfht/w_%d/radius_bin.hdf5", parent_path, foreground_name, area_id);
			sprintf(set_name, "/radius_bin");
			write_h5(temp_path, set_name, radius_bin, radius_num + 1, 1, TRUE);

			delete[] my_data_buf;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
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

}

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
		std::cout << log_infom << times << std::endl;
	}


	for (i = 0; i < backgal_data_col; i++)
	{
		delete[] backgal_data[i];
	}
	for (i = grid_data_col; i < foregal_data_col + grid_data_col; i++)
	{
		delete[] foregal_data[i];
	}
	delete[] block_mask;
	delete[] radius_bin;

	MPI_Finalize();
	return 0;
}
