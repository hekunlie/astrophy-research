#include<FQlib.h>
#include<mpi.h>
#include<vector>

#define max_data_col 40
#define foregal_data_col 5
#define grid_data_col 5
#define backgal_data_col 23
#define mg_bin_num 12

int main(int argc, char *argv[])
{
	/* the macros should be adjusted in each case depends on how many columns will be read */
	/* The input paremeters:																												*/
	/*		1. the sky area label																												*/
	/*		2. the radius bin label, the search raidus																			*/
	/*		3. the start of the guess signal, PDF-SYM method																*/
	/*		4. the end of the guess signal	, PDF-SYM method																*/
	/*		5. the number of points between the end and start of guess signal, PDF-SYM method	*/

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);


	int i, j, k;
	char data_path[250], log_path[250], h5f_path_grid[250], h5f_path_fore[250], h5f_res_path[250], temp_path[300];
	char set_name[50], set_name_2[50], attrs_name[80], log_infom[300];
	char foreground_name[50];

	// the controllers
	int area_id = atoi(argv[1]);
	int radius_label = atoi(argv[2]);
	double gh_start = atof(argv[3]);
	double gh_end = atof(argv[4]);
	strcpy(foreground_name, argv[5]);

	   
	//sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/");
	//sprintf(log_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/log/ggl_log_%d.dat", rank);
	//sprintf(h5f_path_fore, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/foreground/%s", foreground_name);

	sprintf(data_path, "/mnt/perc/hklee/CFHT/gg_lensing/data/");
	sprintf(log_path, "/mnt/perc/hklee/CFHT/gg_lensing/log/ggl_log_%d.dat", rank);
	sprintf(h5f_path_fore, "/mnt/perc/hklee/CFHT/gg_lensing/data/foreground/%s", foreground_name);

	sprintf(h5f_res_path, "%s/w_%d/radius_%d.hdf5", data_path, area_id, radius_label);
	
	double g_step;
	double gh_left = -130 + 10 * rank;
	double gh_right = fabs(gh_left);
	int gh_num = int(gh_right * 2);
	g_step = (gh_end - gh_start) / (gh_num - 1);

	double *gh = new double[gh_num];
	for (i = 0; i < gh_num; i++)
	{
		gh[i] = gh_left + i;
	}

	pts_info gal_info;

	double st_start, st_end, st1, st2, st3, st4, stgal;
	int process_per;
	int per_n;


	int foregal_num;
	int my_gal_s, my_gal_e, gal_id;
	double *foregal_data[max_data_col];
	long pair_count;// be carefull, the pair number may be too many, long or double 
	double z_f, ra_f, dec_f;
	double dist_len, dist_source, dist_len_coeff;
	double coeff, coeff_inv;
	double crit_surf_density_com;


	int backgal_num;
	double *backgal_data[max_data_col]; //backgal_data_col = 17
	double backgal_cos_2phi, backgal_sin_2phi, backgal_cos_4phi, backgal_sin_4phi;
	double backgal_mg_tan, backgal_mg_cross, backgal_mn_tan, backgal_mu_tan;
	double z_b, z_thresh, ra_b, dec_b;
	double diff_ra, diff_dec, diff_r, diff_theta, diff_theta_sq, diff_z_thresh = 0.1;
	   
	// the chi and the shear guess
	int ig, ic, ig_label;
	int mg_bin_num2 = mg_bin_num / 2;

	// the bin of G1(2) for shear estimation
	double *mg_bin = new double[mg_bin_num + 1];

	// chi of the signal from the all areas
	long *chi_tan_shared, *chi_cross_shared, *pair_count_shared;
	// chi square of the signal of each thread in each areas
	long *my_chi_tan = new long[mg_bin_num*gh_num];
	long *my_chi_cross = new long[mg_bin_num*gh_num];
	double mg_t, mg_x;
	int chi_bin_label;
	
	int grid_num, grid_ny, grid_nx;
	int row, col, bin_label;
	int block_id, block_s, block_e, ib;
	double block_scale[1];
	int *block_mask;
	int ra_bin_num, dec_bin_num;


	int radius_num;
	double radius_s, radius_e, radius_e_sq;
	double *radius_bin;

	int nib_id = 0, bs_id = 1, be_id = 2, bdy_id = 3, bdx_id = 4;
	int z_id = 5, dist_id = 6, ra_id = 7, dec_id = 8, cos_dec_id = 9;
	int mg1_id = 10, mg2_id = 11, mn_id = 12, mu_id = 13, mv_id = 14;
	int zmin_lb = 15, zmax_lb = 16, odds_lb = 17, mag_lb = 18;
	int ra_bin_id = 19, dec_bin_id = 20, block_scale_id = 21, grid_shape_id = 22;

	int shape[2];

	char *names[max_data_col];//backgal_data_col
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

	sprintf(names[mg1_id], "G1");
	sprintf(names[mg2_id], "G2");
	sprintf(names[mn_id], "N");
	sprintf(names[mu_id], "U");
	sprintf(names[mv_id], "V");

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

	// read the search radius
	sprintf(h5f_path_grid, "%scata_result_ext_grid.hdf5", data_path);
	sprintf(set_name, "/radius_bin");
	read_h5_datasize(h5f_path_grid, set_name, radius_num);
	radius_bin = new double[radius_num] {};
	read_h5(h5f_path_grid, set_name, radius_bin);
	radius_num = radius_num - 1;

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// the shared buffer for the total chi square of the signal
	MPI_Win win_chi_tan_total, win_chi_cross_total, win_pair_count;
	MPI_Aint size_chi_tan, size_chi_cross, size_pair_count;

	// [chi_tan, chi_cross]
	size_chi_tan =  mg_bin_num*gh_num;
	size_chi_cross = mg_bin_num*gh_num;
	size_pair_count = numprocs;
	if (0 == rank)
	{
		MPI_Win_allocate_shared(size_chi_tan * sizeof(long), sizeof(long), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_tan_shared, &win_chi_tan_total);
		MPI_Win_allocate_shared(size_chi_cross * sizeof(long), sizeof(long), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_cross_shared, &win_chi_cross_total);
		MPI_Win_allocate_shared(size_pair_count * sizeof(long), sizeof(long), MPI_INFO_NULL, MPI_COMM_WORLD, &pair_count_shared, &win_pair_count);
	}
	else
	{
		int dispu_total;
		MPI_Win_allocate_shared(0, sizeof(long), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_tan_shared, &win_chi_tan_total);
		MPI_Win_shared_query(win_chi_tan_total, 0, &size_chi_tan, &dispu_total, &chi_tan_shared);

		MPI_Win_allocate_shared(0, sizeof(long), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_cross_shared, &win_chi_cross_total);
		MPI_Win_shared_query(win_chi_cross_total, 0, &size_chi_cross, &dispu_total, &chi_cross_shared);

		MPI_Win_allocate_shared(0, sizeof(long), MPI_INFO_NULL, MPI_COMM_WORLD, &pair_count_shared, &win_pair_count);
		MPI_Win_shared_query(win_pair_count, 0, &size_pair_count, &dispu_total, &pair_count_shared);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// initialization
	if (0 == rank)
	{	
		initialize_arr(chi_tan_shared, mg_bin_num*gh_num, 0);
		initialize_arr(chi_cross_shared, mg_bin_num*gh_num, 0);
		initialize_arr(pair_count_shared, numprocs, 0);
	}
	initialize_arr(my_chi_tan, mg_bin_num*gh_num, 0);
	initialize_arr(my_chi_cross, mg_bin_num*gh_num, 0);
	MPI_Barrier(MPI_COMM_WORLD);
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// read foreground information
	// Z, DISTANCE, RA, DEC, COS_DEC
	for (i = 0; i < foregal_data_col; i++)
	{
		sprintf(set_name, "/w_%d/%s", area_id, names[i + grid_data_col]);
		read_h5_datasize(h5f_path_grid, set_name, foregal_num);

		foregal_data[i] = new double[foregal_num];
		read_h5(h5f_path_grid, set_name, foregal_data[i]);
	}
	sprintf(set_name, "/w_%d/mg_bin", area_id);
	read_h5_attrs(h5f_path_grid, set_name, attrs_name, shape, "d");
	if (shape[0] != mg_bin_num + 1)
	{	//check
		std::cout << "The shapes of the G-bin in date and code dosen't match each other. (" << shape[0] << ", " << mg_bin_num + 1 << ").";
		exit(0);
	}
	read_h5(h5f_path_grid, set_name, mg_bin);
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
	sprintf(set_name, "/background/w_%d", area_id);
	sprintf(attrs_name, "grid_shape");
	read_h5_attrs(h5f_path_grid, set_name, attrs_name, shape, "g");
	grid_ny = shape[0];
	grid_nx = shape[1];
	grid_num = grid_ny * grid_nx;
	sprintf(attrs_name, "block_scale");
	read_h5_attrs(h5f_path_grid, set_name, attrs_name, block_scale, "g");

	// Z, RA, DEC,  G1, G2, N, U, V,  num_in_block,  block_start, block_end, 
	// block_boundx, block_boundy
	for (i = 0; i < backgal_data_col; i++)
	{
		sprintf(set_name, "/background/w_%d/%s", area_id, names[i]);
		sprintf(attrs_name, "shape");
		read_h5_attrs(h5f_path_grid, set_name, attrs_name, shape, "d");

		if (z_id == i)
		{
			backgal_num = shape[0];
		}
		if (ra_bin_id == i)
		{
			ra_bin_num = shape[0] - 1;
		}
		if (dec_bin_id == i)
		{
			dec_bin_num = shape[0] - 1;
		}
		backgal_data[i] = new double[shape[0] * shape[1]]{};
		read_h5(h5f_path_grid, set_name, backgal_data[i]);
	}

	block_mask = new int[grid_num];

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
		radius_e = radius_bin[radius_label + 1] * coeff / dist_len / foregal_data[cos_dec_id][gal_id] * 1.4; // degree
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
		//std::cout << row << " " << col << " " << block_scale[0] << " " << radius_e << " " << dist_len <<" "<< coeff << std::endl;

		// find the blocks needed //
		initialize_arr(block_mask, grid_num, -1);
		find_block(&gal_info, radius_e, backgal_data[bdy_id], backgal_data[bdx_id], block_mask);
		//for (i = 0; i < grid_num; i++)
		//{
		//	if (block_mask[i] > -1)
		//	{
		//		std::cout << block_mask[i] << " " << block_mask[i] / grid_nx << "  " << block_mask[i] % grid_nx << " Row col " << row << " " << col<<" Ra Dec "<<ra_f<<" "<< dec_f << std::endl;
		//	}
		//}
		//sprintf(temp_path, "!/home/hklee/work/mask.fits");
		//write_fits(temp_path, block_mask, grid_ny, grid_nx);
		//std::cout << radius_e << " Row col " << row << " " << col << " Ra Dec " << ra_f << " " << dec_f <<std::endl;

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
					if (backgal_data[z_id][ib] >= z_thresh)
					{
						ra_b = backgal_data[ra_id][ib];
						dec_b = backgal_data[dec_id][ib];

						dist_source = backgal_data[dist_id][ib];

						// times cos(dec) due to the different length the arc corresponding to the same delta R.A. at different Dec
						diff_ra = (ra_b - ra_f)*foregal_data[cos_dec_id][gal_id];
						diff_dec = dec_b - dec_f;
						diff_theta_sq = diff_ra * diff_ra + diff_dec * diff_dec; // degree^2

						// the seperation in comving coordinate, 
						diff_r = dist_source * sqrt(diff_theta_sq)*coeff_inv;

						if (diff_r >= radius_bin[radius_label] and diff_r < radius_bin[radius_label + 1])
						{
							
							crit_surf_density_com = dist_source / (dist_source - dist_len) *dist_len_coeff;

							// rotation for shear calculation
							backgal_sin_2phi = 2 * diff_ra*diff_dec / diff_theta_sq;
							backgal_cos_2phi = (diff_dec - diff_ra)*(diff_ra + diff_dec) / diff_theta_sq;

							backgal_sin_4phi = 2 * backgal_sin_2phi * backgal_cos_2phi;
							backgal_cos_4phi = (backgal_cos_2phi + backgal_sin_2phi)*(backgal_cos_2phi - backgal_sin_2phi);

							// G_t = (G_1 + i*G_2)*EXP(2i\phi) =  G_1 *cos2\phi - G_2*sin2\phi
							// the direction of R.A. is oppsite, actually,  G_t =  G_1 *cos2\phi + G_2*sin2\phi
							// \Sigma_crit *G_t(x)
							backgal_mg_tan = crit_surf_density_com*( backgal_data[mg1_id][ib] * backgal_cos_2phi + backgal_data[mg2_id][ib] * backgal_sin_2phi);
							// the cross components
							backgal_mg_cross = crit_surf_density_com * (backgal_data[mg1_id][ib] * backgal_sin_2phi - backgal_data[mg2_id][ib] * backgal_cos_2phi);
							// scalar
							backgal_mn_tan = backgal_data[mn_id][ib];
							// U_t = Re[(U+i*V)*EXP(-4i\phi)] = U*cos4\phi + V*sin\4phi
							backgal_mu_tan = backgal_data[mu_id][ib] * backgal_cos_4phi - backgal_data[mv_id][ib] * backgal_sin_4phi;
							
							// calculate the PDF of the estimator
							for (ig = 0; ig < gh_num; ig++)
							{
								mg_t = backgal_mg_tan - gh[ig] * (backgal_mn_tan + backgal_mu_tan);
								histogram_s(mg_t, mg_bin, mg_bin_num, chi_bin_label);
								my_chi_tan[ig*mg_bin_num + chi_bin_label];

								mg_x = backgal_mg_tan - gh[ig] * (backgal_mn_tan - backgal_mu_tan);
								histogram_s(mg_x, mg_bin, mg_bin_num, chi_bin_label);
								my_chi_cross[ig*mg_bin_num + chi_bin_label];
							}
							pair_count_shared[rank] += 1;

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

	// merge the PDF of signal estimator into the shared buffer
	for(i=0; i<numprocs; i++)
	{
		if (rank == i)
		{	
			if (rank == 0)
			{
				sum_arr(pair_count_shared, numprocs, 0, numprocs, pair_count);
				sprintf(log_infom, "RANK: %d. w_%d. %ld galaxies have been found in Radius [%.4f, %.4f].", rank, area_id, pair_count, radius_bin[radius_label], radius_bin[radius_label+1]);
				std::cout << log_infom << std::endl;
			}
			for (j = 0; j < mg_bin_num*gh_num; j++)
			{
				chi_tan_shared[j] += my_chi_tan[j];
				chi_cross_shared[j] += my_chi_cross[j];
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// calculate the result and write to file
	if (0 == rank)
	{
		// the chi square for fitting shear
		long *chi_block = new long[mg_bin_num];
		double *chisq_tan = new double[gh_num];
		double *chisq_cross = new double[gh_num];
		double chi_temp, crit, crit_sig;
		double crit_result[4];

		for (i = 0; i < gh_num; i++)
		{	// tangential 
			for (j = 0; j < mg_bin_num; j++)
			{
				chi_block[j] = chi_tan_shared[i*mg_bin_num + j];
			}
			cal_chisq_1d(chi_block, mg_bin_num, chi_temp);
			chisq_tan[i] = chi_temp;
			// cross
			for (j = 0; j < mg_bin_num; j++)
			{
				chi_block[j] = chi_cross_shared[i*mg_bin_num + j];
			}
			cal_chisq_1d(chi_block, mg_bin_num, chi_temp);
			chisq_cross[i] = chi_temp;
		}

		// write the result to file
		sprintf(set_name, "/chi_tan");
		write_h5(h5f_res_path, set_name, chi_tan_shared, gh_num, mg_bin_num, TRUE);
		sprintf(set_name, "/chi_cross");
		write_h5(h5f_res_path, set_name, chi_cross_shared, gh_num, mg_bin_num, FALSE);

		sprintf(set_name, "/chi_square_tan");
		write_h5(h5f_res_path, set_name, chisq_tan, gh_num, 1, FALSE);
		sprintf(set_name, "/chi_sqaure_cross");
		write_h5(h5f_res_path, set_name, chisq_cross, gh_num, 1, FALSE);

		// tangential
		fit_shear(gh, chisq_tan, gh_num, crit, crit_sig, 50);
		crit_result[0] = crit; 
		crit_result[1] = crit_sig;
		// cross
		fit_shear(gh, chisq_cross, gh_num, crit, crit_sig, 50);
		crit_result[2] = crit;
		crit_result[3] = crit_sig;
		sprintf(set_name, "/result");
		write_h5(h5f_res_path, set_name, crit_result, 4, 1, FALSE);

		delete[] chi_block;
		delete[] chisq_cross;
		delete[] chisq_tan;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// free the memory	
	MPI_Win_free(&win_chi_tan_total);
	MPI_Win_free(&win_chi_cross_total);
	MPI_Win_free(&win_pair_count);


	for (i = 0; i < backgal_data_col; i++)
	{
		delete[] backgal_data[i];
	}
	for (i = 0; i < foregal_data_col; i++)
	{
		delete[] foregal_data[i];
	}
	delete[] block_mask;
	delete[] mg_bin;
	delete[] gh;
	delete[] radius_bin;
	delete[] my_chi_tan;
	delete[] my_chi_cross;

	sprintf(log_infom, "RANK: %d. Write file of Radius bin [%.4f,  %.4f] ---- end. ", rank, radius_bin[radius_label], radius_bin[radius_label + 1]);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
	}

	MPI_Finalize();
	return 0;
}
