#include<FQlib.h>
#include<mpi.h>

# define backgal_data_col 15

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
	double time_label[10], st_start, st_end;
	int process_per;
	double per_n;
	int area_id, area_num;

	area_num = 4;
	int tag;
	int foregal_num, gal_id;
	int my_gal_s, my_gal_e;
	double *foregal_data[3];
	double *foregal_data_shared;
	int foregal_data_shared_col, shared_col_in_radi, shared_elem_id;
	double pair_count;
	double shear_tan, shear_tan_sig, shear_cros;
	double z_f, ra_f, dec_f;
	double dist_len, dist_source, dist_len_source;
	double coeff, da_coeff;
	double crit_surf_density_comoving;

	int data_num;
	int backgal_num;
	double *backgal_data[backgal_data_col];
	int *backgal_mask;
	double *backgal_cos_2phi, *backgal_sin_2phi;
	double z_b, z_thresh, ra_b, dec_b;
	double diff_ra, diff_dec, diff_r;
	double back_mg1, back_mg2, back_mnu1, back_mnu2;
	int back_tag;


	double *ra_bin, *dec_bin;
	int ra_bin_num, dec_bin_num;
	

	// the chi and the shear guess
	int g_num = 100, ig, ic, ig_label;
	int mg_bin_num = 12, mg_bin_num2 = mg_bin_num / 2;
	double *mg1_bin = new double[mg_bin_num + 1];
	double *mg2_bin = new double[mg_bin_num + 1];
	double *chi_1 = new double[mg_bin_num*g_num];
	double *chi_2 = new double[mg_bin_num*g_num];
	double *chisq_1 = new double[g_num];
	double *chisq_2 = new double[g_num];
	double *gh = new double[g_num];
	double g_step;
	double chi_temp1, chi_temp2;
	double gh1, gh1_sig, gh2, gh2_sig;
	g_step = 0.17 / g_num;
	for (i = 0; i < g_num; i++)
	{
		gh[i] = -0.085 + i * g_step;
	}

	int grid_num, grid_ny, grid_nx;
	int row, col, bin_label;
	int block_id, block_s, block_e, ib;
	double block_scale[1];
	int *block_mask;


	int radius_num, radi_id;
	double radius_s, radius_e;
	double radius_s_sq, radius_e_sq;
	double *radius_bin;

	int z_id = 0, ra_id = 1, dec_id = 2;
	int mg1_id = 3, mg2_id = 4, mn_id = 5, mu_id = 6, mv_id = 7;
	int nib_id = 8, bs_id = 9, be_id = 10, bdy_id = 11, bdx_id = 12;
	int ra_bin_id = 13, dec_bin_id = 14;

	double *redshifts, *distances;
	int red_num;
	int shape[2];
	char data_path[200], log_path[200], h5f_path[200], h5f_res_path[200];
	char set_name[80], attrs_name[80], log_infom[200];

	char *names[backgal_data_col];
	for (i = 0; i < backgal_data_col; i++)
	{
		names[i] = new char[25];
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

	sprintf(names[ra_bin_id], "RA_bin");
	sprintf(names[dec_bin_id], "DEC_bin");

	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/");
	sprintf(log_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/log/ggl_log_%d.dat", rank);
	sprintf(h5f_res_path, "%sggl_result.hdf5", data_path);

	sprintf(log_infom, "RANK: %d. Start ...", rank);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
	}

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

	sprintf(log_infom, "RANK: %d. Read redshift.hdf5", rank);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
	}

	// read the search radius
	sprintf(h5f_path, "%scata_result_ext_grid.hdf5", data_path);
	sprintf(set_name, "/radius_bin");
	sprintf(attrs_name, "shape");
	read_h5_attrs(h5f_path, set_name, attrs_name, shape, "d");
	radius_num = shape[0];
	radius_bin = new double[radius_num] {};
	read_h5(h5f_path, set_name, radius_bin);
	radius_num = radius_num - 1;

	sprintf(log_infom, "RANK: %d. Read radius_bin in cata_result_ext_grid.hdf5", rank);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
		std::cout << red_num << " " << radius_num << std::endl;
	}

	MPI_Win win_1;
	MPI_Aint size_1;
	for (area_id = 1; area_id < area_num + 1; area_id++)
	{
		st_start = clock();
		time_label[0] = clock();
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

		sprintf(log_infom, "RANK: %d. w_%d. Read foreground data. %d galaxies", rank, area_id, foregal_num);
		write_log(log_path, log_infom);
		if (0 == rank)
		{
			std::cout << log_infom << std::endl;
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
		for (i = 0; i < backgal_data_col; i++)
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
			backgal_data[i] = new double[shape[0] * shape[1]]{};
			read_h5(h5f_path, set_name, backgal_data[i]);
		}


		backgal_mask = new int[backgal_num];
		backgal_cos_2phi = new double[backgal_num];
		backgal_sin_2phi = new double[backgal_num];

		// set the bins for g1 & g2 estimation
		set_bin(backgal_data[mg1_id], backgal_num, mg1_bin, mg_bin_num, 100);
		set_bin(backgal_data[mg2_id], backgal_num, mg2_bin, mg_bin_num, 100);
		if (0 == rank)
		{
			show_arr(mg1_bin, 1, mg_bin_num + 1);
			show_arr(mg2_bin, 1, mg_bin_num + 1);
		}
		sprintf(log_infom, "RANK: %d. w_%d. Read background data. %d galaxies. %d grids (%d x %d)", rank, area_id, backgal_num, grid_num, grid_ny, grid_nx);
		write_log(log_path, log_infom);
		if (0 == rank)
		{
			std::cout << log_infom << std::endl;
		}


		shared_col_in_radi = 7;
		foregal_data_shared_col = radius_num * shared_col_in_radi;
		size_1 = foregal_num * foregal_data_shared_col * sizeof(double);
		if (0 == rank)
		{
			MPI_Win_allocate_shared(size_1, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &foregal_data_shared, &win_1);
		}
		else
		{
			int disp_unit;
			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &foregal_data_shared, &win_1);
			MPI_Win_shared_query(win_1, 0, &size_1, &disp_unit, &foregal_data_shared);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (0 == rank)
		{
			for (i = 0; i < foregal_num * foregal_data_shared_col; i++)
			{
				foregal_data_shared[i] = 0;
			}
		}
		sprintf(log_infom, "RANK: %d. w_%d. Create & initialize the shared buffer [%d,  %d]", rank, area_id, foregal_num, foregal_data_shared_col);
		write_log(log_path, log_infom);
		if (0 == rank)
		{
			std::cout << log_infom << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);

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
		time_label[1] = clock();
		sprintf(log_infom, "RANK: %d. w_%d. Task distribution: my gal: %d ~ %d. Time since begin: %.2f sec", rank, area_id, my_gal_s, my_gal_e, (time_label[1] - time_label[0]) / CLOCKS_PER_SEC);
		write_log(log_path, log_infom);
		if (0 == rank)
		{
			std::cout << log_infom << std::endl;
		}

		// loop the foreground galaxy
		coeff = 0.18 / C_0_hat / Pi;
		da_coeff = C_0_hat / H_0_hat;

		for (gal_id = my_gal_s; gal_id < my_gal_e; gal_id++)
		{
			time_label[2] = clock();

			//process_per = (gal_id - my_gal_s) % (int((my_gal_e - my_gal_s) *0.1));
			//per_n = double((gal_id - my_gal_s) / (my_gal_e - my_gal_s) * 100);
			//
			//if (0 == process_per)
			//{
			//	sprintf(log_infom, "RANK: %d. w_%d. Looping foreground %.2f%%. Time since begin: %.2f sec", rank, area_id, per_n, (time_label[2] - time_label[1]) / CLOCKS_PER_SEC);
			//	write_log(log_path, log_infom);
			//	if (0 == rank)
			//	{
			//		std::cout << log_infom << std::endl;
			//	}
			//}

			z_f = foregal_data[z_id][gal_id];
			z_thresh = z_f + 0.3;
			find_near(redshifts, z_f, red_num, tag);
			dist_len = distances[tag];
			// degree
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
			
			// loop the search radius//
			for (radi_id = 0; radi_id < radius_num; radi_id++)
			{
				// the searching radius depend on the redshift of lens
				radius_s = radius_bin[radi_id] * coeff / dist_len * (1 + z_f); // the physical distance
				radius_e = radius_bin[radi_id + 1] * coeff / dist_len * (1 + z_f);
				radius_s_sq = radius_s * radius_s;
				radius_e_sq = radius_e * radius_e;
				//std::cout << radius_s << " " << radius_s_sq << " " << radius_e << " " << radius_e_sq << std::endl;
				// find the blocks needed
				initialize_arr(block_mask, grid_num, -1);
				initialize_arr(backgal_mask, backgal_num, 0);
				find_block(&gal_info, radius_s, radius_e, backgal_data[bdy_id], backgal_data[bdx_id], block_mask);
				// the background galaxy may belong to many foreground galaxies in the same (or not) radius scale
				// depending on the foreground galaxy density
				// it must be initialized for each foreground galaxy each radius bin
				initialize_arr(backgal_sin_2phi, backgal_num, 0);
				initialize_arr(backgal_cos_2phi, backgal_num, 0);

				initialize_arr(chi_1, mg_bin_num*g_num, 0);
				initialize_arr(chi_2, mg_bin_num*g_num, 0);

				pair_count = 0;
				// loop the found blocks and calculate//
				for (block_id = 0; block_id < grid_num; block_id++)
				{
					if (block_mask[block_id] > -1)
					{
						// the start and end point of the block//
						// the start- & end-point
						block_s = backgal_data[bs_id][block_id];
						block_e = backgal_data[be_id][block_id];
						for (ib = block_s; ib < block_e; ib++)
						{
							ra_b = backgal_data[ra_id][ib];
							dec_b = backgal_data[dec_id][ib];
							diff_ra = ra_b - ra_f;
							diff_dec = dec_b - dec_f;
							diff_r = diff_ra * diff_ra + diff_dec * diff_dec;

							if (diff_r >= radius_s_sq and diff_r < radius_e_sq and backgal_data[z_id][ib]>z_thresh)
							{
								backgal_mask[ib] = 1;
								pair_count += 1;
								// the position angle of background galaxy respect to the foreground
								// the cos(-2\phi) & sin(-2\phi)
								backgal_cos_2phi[ib] = (diff_ra*diff_ra - diff_dec * diff_dec) / diff_r;
								backgal_sin_2phi[ib] = -2 * diff_ra*diff_dec / diff_r;

								back_mnu1 = backgal_data[mn_id][ib] + backgal_data[mu_id][ib];
								back_mnu2 = backgal_data[mn_id][ib] - backgal_data[mu_id][ib];

								for (ig = 0; ig < g_num; ig++)
								{
									ig_label = ig * mg_bin_num;
									back_mg1 = backgal_data[mg1_id][ib] - gh[ig] * back_mnu1;
									back_mg2 = backgal_data[mg2_id][ib] - gh[ig] * back_mnu2;

									histogram_s(back_mg1, mg1_bin, mg_bin_num, back_tag);
									chi_1[ig_label + back_tag] += 1;
									//std::cout << back_mg1 << " " << back_tag << std::endl;
									histogram_s(back_mg2, mg2_bin, mg_bin_num, back_tag);
									chi_2[ig_label + back_tag] += 1;
									//std::cout << back_mg2 << " " << back_tag << std::endl;
								}
							}

						}
					}
				}

				//show_arr(chi_1, mg_bin_num, g_num);
				// esitmate the shear in [raidus_s, radius_e]//
				shared_elem_id = gal_id * foregal_data_shared_col + radi_id * shared_col_in_radi;
				if (pair_count > 10000)
				{
					// skip if too little pair 
					initialize_arr(chisq_1, g_num, 0);
					initialize_arr(chisq_2, g_num, 0);
					for (ig = 0; ig < g_num; ig++)
					{
						ig_label = ig * mg_bin_num;

						for (ic = 0; ic < mg_bin_num2; ic++)
						{
							chi_temp1 = chi_1[ig_label + ic] - chi_1[ig_label + ic + mg_bin_num2];
							chi_temp2 = chi_1[ig_label + ic] + chi_1[ig_label + ic + mg_bin_num2];
							chisq_1[ig] += chi_temp1 * chi_temp1 / chi_temp2;

							chi_temp1 = chi_2[ig_label + ic] - chi_2[ig_label + ic + mg_bin_num2];
							chi_temp2 = chi_2[ig_label + ic] + chi_2[ig_label + ic + mg_bin_num2];
							chisq_2[ig] += chi_temp1 * chi_temp1 / chi_temp2;
						}
						chisq_1[ig] = chisq_1[ig] * 0.5;
						chisq_2[ig] = chisq_2[ig] * 0.5;
					}
					std::cout << rank << " " << radi_id<<" "<<radius_s<<radius_e<<std::endl;
					show_arr(chisq_1, 1, g_num);
					show_arr(chisq_2, 1, g_num);

					fit_shear(gh, chisq_1, g_num, gh1, gh1_sig);
					fit_shear(gh, chisq_2, g_num, gh2, gh2_sig);

					foregal_data_shared[shared_elem_id + 0] = gh1;
					// because the direction of RA is opposite to that of the real
					foregal_data_shared[shared_elem_id + 1] = -gh2;

					foregal_data_shared[shared_elem_id + 2] = gh1_sig;
					foregal_data_shared[shared_elem_id + 3] = gh2_sig;

					// calculate the \Sum g_t \Sigma_crit
					for (ib = 0; ib < backgal_num; ib++)
					{
						if (backgal_mask[ib] == 1)
						{
							z_b = backgal_data[z_id][ib];
							find_near(redshifts, z_b, red_num, tag);
							dist_source = distances[tag];
							dist_len_source = dist_source - dist_len;

							// only the comoving distance part of the real critical surface density
							// it will be multiplied by the factor
							crit_surf_density_comoving = dist_source / dist_len_source / dist_len / (1 + z_f);
							shear_tan = -gh1 * backgal_cos_2phi[ib] - gh2 * backgal_sin_2phi[ib];
							shear_tan_sig = -gh1_sig * backgal_cos_2phi[ib] - gh2_sig * backgal_sin_2phi[ib];
							//shear_cros = gh1*backgal_sin_2phi[ib] - gh2*backgal_cos_2phi[ib];

							shared_elem_id = ib * foregal_data_shared_col + radi_id * shared_col_in_radi;
							foregal_data_shared[shared_elem_id + 4] += shear_tan * crit_surf_density_comoving;
							foregal_data_shared[shared_elem_id + 5] += shear_tan_sig * crit_surf_density_comoving;
						}
					}
					foregal_data_shared[gal_id * foregal_data_shared_col + radi_id * shared_col_in_radi + 6] = pair_count;
				}
			}
			time_label[3] = clock();
			if (0 == process_per)
			{
				sprintf(log_infom, "RANK: %d. w_%d. Looping foreground %.2f%%. Time used: %.2f sec", rank, area_id, per_n, (time_label[3] - time_label[2]) / CLOCKS_PER_SEC);
				write_log(log_path, log_infom);
				if (0 == rank)
				{
					std::cout << log_infom << std::endl;
				}
			}
		

		}

		// finish the area_i
		//MPI_Barrier(MPI_COMM_WORLD);
		//if(0 == rank)
		//{
		//	sprintf(set_name, "/w_%d", area_id);
		//	if (1 == area_id)
		//	{
		//		write_h5(h5f_res_path, set_name, foregal_data_shared, foregal_data_shared_col, foregal_num, TRUE);
		//	}
		//	else
		//	{
		//		write_h5(h5f_res_path, set_name, foregal_data_shared, foregal_data_shared_col, foregal_num, FALSE);
		//	}
		//}
		//MPI_Barrier(MPI_COMM_WORLD);
		//st_end = clock();
		//sprintf(log_infom, "RANK: %d. w_%d. Finish aread %d. Time: %.2f sec", rank, area_id, area_id, (st_end - st_start) / CLOCKS_PER_SEC);
		//write_log(log_path, log_infom);
		//if (0 == rank)
		//{
		//	std::cout << log_infom << std::endl;
		//}

		// free the memory	
		MPI_Win_free(&win_1);

		delete[] backgal_sin_2phi;
		delete[] backgal_cos_2phi;
		delete[] backgal_mask;

		for (i = 0; i < backgal_data_col; i++)
		{
			delete[] backgal_data[i];
		}

		delete[] block_mask;

		for (i = 0; i < 3; i++)
		{
			delete[] foregal_data[i];
		}

		sprintf(log_infom, "RANK: %d. w_%d. Free the memory", rank, area_id);
		write_log(log_path, log_infom);
		if (0 == rank)
		{
		std::cout << log_infom << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	delete[] mg1_bin;
	delete[] mg2_bin;
	delete[] chi_1;
	delete[] chi_2;
	delete[] chisq_1;
	delete[] chisq_2;
	delete[] gh;
	delete[] redshifts;
	delete[] distances;
	delete[] radius_bin;

	sprintf(log_infom, "RANK: %d. w_%d. Free the memory the last one", rank, area_id);
	write_log(log_path, log_infom);
	std::cout << log_infom << std::endl;

	MPI_Finalize();
	return 0;
}