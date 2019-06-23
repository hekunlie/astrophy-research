#include<FQlib.h>
#include<vector>

#define max_data_col 40
#define foregal_data_col 5
#define grid_data_col 5
#define backgal_data_col 21
#define mg_bin_num 8

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



	// if the memory of the system could contain all the pairs
	// it will gather all the data and estimate the signal with SYM-PDF method
#if defined (SMALL_CATA)
	std::vector<double> data_cache;
#endif

	int i, j, k, temp;
	char data_path[250], log_path[250], h5f_path_grid[250], h5f_path_fore[250], h5f_res_path[250], temp_path[300];
	char set_name[50], set_name_2[50], attrs_name[80], log_infom[300];
	char foreground_name[50];
	
	pts_info gal_info;

	double st_start, st_end, st1, st2, st3, st4, stgal;
	int process_per;
	int per_n;

	int foregal_num;
	int my_gal_s, my_gal_e, gal_id;
	double *foregal_data[max_data_col];
	MY_INT pair_count = 0;// be carefull, the pair number may be too many, long or double 
	double z_f, ra_f, dec_f;
	double dist_len, dist_source, dist_len_coeff;
	double coeff, coeff_inv, coeff_rad_dist;
	double crit_surf_density_com;


	int backgal_num;
	double *backgal_data[max_data_col];
	double backgal_cos_2phi, backgal_sin_2phi, backgal_cos_4phi, backgal_sin_4phi;
	double backgal_mg_tan, backgal_mg_cross, backgal_mn_tan, backgal_mu_tan;
	double backgal_mnu1_tan, backgal_mnu2_tan, backgal_mnu1_tan_c, backgal_mnu2_tan_c;
	double z_b, z_thresh, z_b_sig95, z_b_odds;
	double ra_b, dec_b;
	double diff_ra, diff_dec, diff_r, diff_theta, diff_theta_sq, diff_z_thresh = 0.1;

	double *my_data_buf, *final_buf;
	// the chi and the shear guess
	int ig, ic, ig_label;
	

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
	int ra_bin_id = 19, dec_bin_id = 20;
	int block_scale_id = 21, grid_shape_id = 22;

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

	// the controllers
	int area_id = atoi(argv[1]);
	int radius_label = atoi(argv[2]);
	strcpy(foreground_name, argv[3]);
	my_gal_s = atoi(argv[4]);
	my_gal_e = my_gal_s + 1;


	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/");
	sprintf(h5f_res_path, "/home/hkli/work/test/radius_%d.hdf5", foreground_name, area_id, radius_label);

	sprintf(h5f_path_fore, "/home/hkli/work/test/cluster_w_%d.hdf5",area_id);

	sprintf(h5f_path_grid, "%scata_result_ext_grid.hdf5", data_path);

	// read the search radius
	sprintf(set_name, "/radius_bin");
	read_h5_datasize(h5f_path_grid, set_name, radius_num);
	radius_bin = new double[radius_num] {};
	read_h5(h5f_path_grid, set_name, radius_bin);
	radius_num = radius_num - 1;

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
	sprintf(log_infom, " w_%d. Read foreground data. %d galaxies", area_id, foregal_num);

	std::cout << log_infom << std::endl;

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

	// Z, Z_MIN, Z_MAX, ODDS,  DISTANCE, RA, DEC,  G1, G2, N, U, V,  num_in_block,  block_start, block_end, 
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

	sprintf(log_infom, "w_%d. Read background data. %d galaxies. %d grids (%d x %d)",  area_id, backgal_num, grid_num, grid_ny, grid_nx);
	std::cout << log_infom << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	st1 = clock();
	

	sprintf(log_infom, "w_%d. Task distribution: my gal: %d ~ %d. ",  area_id, my_gal_s, my_gal_e);
	write_log(log_path, log_infom);

	std::cout << log_infom << std::endl;

	coeff = 0.18 / C_0_hat / Pi;
	coeff_inv = C_0_hat * Pi / 0.18;
	coeff_rad_dist = C_0_hat * 1000;

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
					z_b_sig95 = z_f + (backgal_data[zmin_lb][ib] + backgal_data[zmax_lb][ib]) / 2;
					z_b_odds = backgal_data[odds_lb][ib];

					//if (backgal_data[z_id][ib] >= z_thresh and backgal_data[z_id][ib] > z_b_sig95 and z_b_odds > 0.5)	
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
						//diff_r = dist_source * sqrt(diff_theta_sq)*coeff_inv;
						separation(ra_b, dec_b, ra_f, dec_f, diff_theta);

						diff_r = dist_source * sin(diff_theta) * coeff_rad_dist;
						if (diff_r >= radius_bin[radius_label] and diff_r < radius_bin[radius_label + 1])
						{
							pair_count += 1;
							crit_surf_density_com = dist_source / (dist_source - dist_len) *dist_len_coeff;

							// rotation for shear calculation, see the NOTE of gg_lensing for the detials 
							backgal_sin_2phi = 2 * diff_ra*diff_dec / diff_theta_sq;
							backgal_cos_2phi = (diff_dec - diff_ra)*(diff_ra + diff_dec) / diff_theta_sq;

							backgal_sin_4phi = 2 * backgal_sin_2phi * backgal_cos_2phi;
							backgal_cos_4phi = (backgal_cos_2phi + backgal_sin_2phi)*(backgal_cos_2phi - backgal_sin_2phi);

							// G_t = (G_1 + i*G_2)*EXP(2i\phi) =  G_1 *cos2\phi - G_2*sin2\phi
							// the direction of R.A. is oppsite, actually,  G_t =  G_1 *cos2\phi + G_2*sin2\phi
							// \Sigma_crit *G_t(x)
							backgal_mg_tan = backgal_data[mg1_id][ib] * backgal_cos_2phi + backgal_data[mg2_id][ib] * backgal_sin_2phi;
							// the cross components
							backgal_mg_cross = backgal_data[mg1_id][ib] * backgal_sin_2phi - backgal_data[mg2_id][ib] * backgal_cos_2phi;
							// scalar
							backgal_mn_tan = backgal_data[mn_id][ib];
							// U_t = Re[(U+i*V)*EXP(-4i\phi)] = U*cos4\phi + V*sin\4phi
							backgal_mu_tan = backgal_data[mu_id][ib] * backgal_cos_4phi - backgal_data[mv_id][ib] * backgal_sin_4phi;
							//if (fabs(backgal_mg_tan) < 1.e-4)
							//{
							//	std::cout << "Rank " << rank << " " << backgal_mg_tan <<" "<< backgal_mg_cross <<" "<< backgal_mn_tan << std::endl;
							//}
							data_cache.push_back(backgal_mg_tan);
							data_cache.push_back(backgal_mg_cross);
							data_cache.push_back(backgal_mn_tan);
							data_cache.push_back(backgal_mu_tan);
							data_cache.push_back(crit_surf_density_com);
							
						}

					}
				}
			}
		}
		stgal = clock();	
	}
	
	my_data_buf = new double[pair_count * 5]{};
	// copy the data in the vector into the buffer 
	if (!data_cache.empty())
	{
		memcpy(my_data_buf, &data_cache[0], data_cache.size() * sizeof(double));
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

	return 0;
}
