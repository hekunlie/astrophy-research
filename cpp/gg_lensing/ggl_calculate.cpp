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
	double time_label[10], st_start, st_end;
	int process_per;
	double per_n;
	int area_id, area_num;

	area_num = 4;
	int tag_f, tag_b;
	int foregal_num, gal_id;
	int my_gal_s, my_gal_e;
	double *foregal_data[3];
	double *foregal_data_shared;
	int foregal_data_shared_col, shared_col_in_radi, shared_elem_id;
	int pair_count;// be carefull, the pair number may be too many, long or double 
	double shear_tan, shear_tan_sig, shear_cros,shear_cros_sig;
	double z_f, ra_f, dec_f;
	double dist_len, dist_source, dist_len_source;
	double coeff;
	double crit_surf_density_com, crit_surf_density_com_sig;
	double *delta_crit;

	shared_col_in_radi = 5;


	int data_num, backgal_data_col=15;
	int backgal_num;
	double *backgal_data[15]; //backgal_data_col
	int *backgal_mask, *my_backgal_mask;
	double backgal_cos_2phi, backgal_sin_2phi, dx_over_dy, dx_over_dy_sq;
	double z_b, z_thresh, ra_b, dec_b;
	double diff_ra, diff_dec, diff_r, diff_theta;
	double back_mg1, back_mg2, back_mnu1, back_mnu2;
	int back_tag;

	
	double *mg1, *mg2, *mnu1, *mnu2;

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
	double *shear_in_radius;


	int z_id = 0, ra_id = 1, dec_id = 2;
	int mg1_id = 3, mg2_id = 4, mn_id = 5, mu_id = 6, mv_id = 7;
	int nib_id = 8, bs_id = 9, be_id = 10, bdy_id = 11, bdx_id = 12;
	int ra_bin_id = 13, dec_bin_id = 14;

	double *redshifts, *distances;
	int red_num;
	int shape[2];
	char data_path[200], log_path[200], h5f_path[200], h5f_res_path[200];
	char set_name[50], set_name_2[50], attrs_name[80], log_infom[200];

	char *names[15];//backgal_data_col
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
	// to store the \Delta \Sigma(R) of each radius bin
	delta_crit = new double[radius_num*2]{};

	sprintf(log_infom, "RANK: %d. Read radius_bin in cata_result_ext_grid.hdf5", rank);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
		std::cout << red_num << " " << radius_num << std::endl;
	}


	for (area_id = 1; area_id < area_num + 1; area_id++)
	{
		st_start = clock();
		time_label[0] = clock();
		// read foreground information
		// Z, RA, DEC
		sprintf(attrs_name, "shape");
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

		my_backgal_mask = new int[backgal_num] {};

		// set the bins for g1 & g2 estimation
		set_bin(backgal_data[mg1_id], backgal_num, mg1_bin, mg_bin_num, 100, 50000);
		set_bin(backgal_data[mg2_id], backgal_num, mg2_bin, mg_bin_num, 100, 50000);
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

		MPI_Win win_foregal, win_mask, win_shear;
		MPI_Aint size_fore, size_mask, size_shear;

		foregal_data_shared_col = radius_num * shared_col_in_radi;
		size_fore = foregal_num * foregal_data_shared_col * sizeof(double);
		size_mask = backgal_num * sizeof(int);
		size_shear = 4*radius_num * sizeof(double);
		if (0 == rank)
		{
			MPI_Win_allocate_shared(size_fore, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &foregal_data_shared, &win_foregal);

			MPI_Win_allocate_shared(size_mask, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &backgal_mask, &win_mask);

			MPI_Win_allocate_shared(size_shear, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shear_in_radius, &win_shear);
		}
		else
		{
			int disp_unit, disp_unit_mask, disp_unit_shear;
			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &foregal_data_shared, &win_foregal);
			MPI_Win_shared_query(win_foregal, 0, &size_fore, &disp_unit, &foregal_data_shared);

			MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &backgal_mask, &win_mask);
			MPI_Win_shared_query(win_mask, 0, &size_mask, &disp_unit_mask, &backgal_mask);

			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shear_in_radius, &win_shear);
			MPI_Win_shared_query(win_shear, 0, &size_shear, &disp_unit_shear, &shear_in_radius);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (0 == rank)
		{
			for (i = 0; i < foregal_num * foregal_data_shared_col; i++)
			{
				foregal_data_shared[i] = 0;
			}
			for (i = 0; i < 4*radius_num; i++)
			{
				shear_in_radius[i] = 0;
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
		
		coeff = 0.18 / C_0_hat / Pi;

		// loop the search radius//
		for (radi_id = 0; radi_id < radius_num; radi_id++)
		{	
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////			1). loop the foreground galaxy to find the background galaxies and label them					   //////
			////			2). rank 0 calculate the shear in this radius bin,  [radius_s, radius_e].										  //////
			////			1). loop the foreground galaxy to calculate the tangential shear and the \delta \Sigma(R)  //////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			//		1). loop the foreground galaxyto find the background galaxies in [radius_s, radius_e].	//
			///////////////////////////////////////////////////////////////////////////////////////////////////////
			if (0 == rank)
			{
				initialize_arr(backgal_mask, backgal_num, 0);
			}
			initialize_arr(my_backgal_mask, backgal_num, 0);
			MPI_Barrier(MPI_COMM_WORLD);

			radius_s = radius_bin[radi_id] * coeff ; 
			radius_e = radius_bin[radi_id + 1] * coeff ;

			for (gal_id = my_gal_s; gal_id < my_gal_e; gal_id++)
			{
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
				find_near(redshifts, z_f, red_num, tag_f);
				dist_len = distances[tag_f];
				// the searching radius depend on the redshift of lens  //
				radius_s = radius_s / dist_len; //comving distance
				radius_e = radius_e / dist_len ;
				radius_s_sq = radius_s * radius_s;
				radius_e_sq = radius_e * radius_e;

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
				find_block(&gal_info, radius_s, radius_e, backgal_data[bdy_id], backgal_data[bdx_id], block_mask);

				for (block_id = 0; block_id < grid_num; block_id++)
				{
					if (block_mask[block_id] > -1)
					{
						// the start and end point of the block //
						// the start- & end-point					      //
						block_s = backgal_data[bs_id][block_id];
						block_e = backgal_data[be_id][block_id];
						for (ib = block_s; ib < block_e; ib++)
						{
							ra_b = backgal_data[ra_id][ib];
							dec_b = backgal_data[dec_id][ib];
							diff_ra = ra_b - ra_f;
							diff_dec = dec_b - dec_f;
							diff_theta = sqrt(diff_ra * diff_ra + diff_dec * diff_dec);

							z_b = backgal_data[z_id][ib];
							find_near(redshifts, z_b, red_num, tag_b);
							diff_r = distances[tag_b]*diff_theta;

							if (backgal_data[z_id][ib] > z_thresh and diff_r < radius_e and diff_r >= radius_s)
							{

								diff_theta = r
								if(dist_len)
								my_backgal_mask[ib] = +1;
							}
						}
					}
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			sprintf(log_infom, "RANK: %d. All the source galaxies have been found in radius bin [%.4f,  %.4f]", rank, radius_bin[radi_id], radius_bin[radi_id+1]);
			if (0 == rank)
			{
				std::cout << log_infom << std::endl;
			}

			for (i = 0; i < numprocs; i++)
			{
				if (i == rank)
				{
					for (j = 0; j < backgal_num; j++)
					{	
						//  if the galaxy is identified as a background galaxy of one foreground galaxy  // 
						//   then it will be labeled. One background galaxy may be the background of   // 
						//   foreground galaxy, but it will be used only once in the shear calculation     // 
						if (my_backgal_mask[j] > 0)
						{
							//  only be used once, otherwise, backgal_mask[j] += 1, // 
							// 	and fix the following process to use it n times            // 
							//backgal_mask[j] = 1;
							backgal_mask[j] += my_backgal_mask[j];
						}
					}
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}		
			sprintf(log_infom, "RANK: %d. Assign the mask to the total mask in radius bin [%.4f,  %.4f]", rank, radius_bin[radi_id], radius_bin[radi_id + 1]);
			if (0 == rank)
			{
				std::cout << log_infom << std::endl;
			}


			/////////////////////////////////////////////////////////////
			//		2). calculate the shear in [radius_s, radius_e]		//
			///////////////////////////////////////////////////////////
			MPI_Barrier(MPI_COMM_WORLD);
			pair_count = 0;
			if (0 == rank)
			{
				// the number of source galaxies
				for (gal_id = 0; gal_id < backgal_num; gal_id++)
				{
					if (backgal_mask[gal_id] > 0)
					{
						pair_count += 1;
					}
				}

				mg1 = new double[pair_count];
				mg2 = new double[pair_count];
				mnu1 = new double[pair_count];
				mnu2 = new double[pair_count];

				pair_count = 0;
				for (gal_id = 0; gal_id < backgal_num; gal_id++)
				{
					// each source galaxy is used just once
					if (backgal_mask[gal_id] > 0)
					{
						mg1[pair_count] = backgal_data[mg1_id][gal_id];
						mg2[pair_count] = backgal_data[mg2_id][gal_id];
						mnu1[pair_count] = backgal_data[mn_id][gal_id] + backgal_data[mu_id][gal_id];
						mnu2[pair_count] = backgal_data[mn_id][gal_id] - backgal_data[mu_id][gal_id];
						pair_count += 1;
					}
				}
				find_shear(mg1, mnu1, pair_count, 12, gh1, gh1_sig);
				find_shear(mg2, mnu2, pair_count, 12, gh2, gh2_sig);

				shear_in_radius[radi_id] = gh1;
				shear_in_radius[radi_id + radius_num] = gh1_sig;
				shear_in_radius[radi_id + radius_num * 2] = gh2;
				shear_in_radius[radi_id + radius_num * 3] = gh2_sig;

				delete[] mg1;
				delete[] mg2;
				delete[] mnu1;
				delete[] mnu2;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			sprintf(log_infom, "RANK: %d. Calculate the shear in radius bin [%.4f,  %.4f]", rank, radius_bin[radi_id], radius_bin[radi_id + 1]);
			if (0 == rank)
			{
				std::cout << log_infom << std::endl;
			}


			///////////////////////////////////////////////
			//		3). caculate the tangential shear		//
			/////////////////////////////////////////////
			radius_s = radius_bin[radi_id] * coeff;
			radius_e = radius_bin[radi_id + 1] * coeff;

			gh1 = shear_in_radius[radi_id];
			gh1_sig = shear_in_radius[radi_id + radius_num];
			gh2 = shear_in_radius[radi_id + radius_num * 2];
			gh2_sig = shear_in_radius[radi_id + radius_num * 3];

			for (gal_id = my_gal_s; gal_id < my_gal_e; gal_id++)
			{
				z_f = foregal_data[z_id][gal_id];
				z_thresh = z_f + 0.3;
				find_near(redshifts, z_f, red_num, tag_f);
				dist_len = distances[tag_f];

				// the searching radius depend on the redshift of lens
				radius_s = radius_s / dist_len * (1 + z_f); // the physical distance
				radius_e = radius_e / dist_len * (1 + z_f);
				radius_s_sq = radius_s * radius_s;
				radius_e_sq = radius_e * radius_e;

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

				// find the blocks needed
				initialize_arr(block_mask, grid_num, -1);
				find_block(&gal_info, radius_s, radius_e, backgal_data[bdy_id], backgal_data[bdx_id], block_mask);

				for (block_id = 0; block_id < grid_num; block_id++)
				{
					if (block_mask[block_id] > -1)
					{
						// the start and end point of the block//
						block_s = backgal_data[bs_id][block_id];
						block_e = backgal_data[be_id][block_id];
						for (ib = block_s; ib < block_e; ib++)
						{
							ra_b = backgal_data[ra_id][ib];
							dec_b = backgal_data[dec_id][ib];
							diff_ra = ra_b - ra_f;
							diff_dec = dec_b - dec_f;
							diff_r = diff_ra * diff_ra + diff_dec * diff_dec;

							z_b = backgal_data[z_id][ib];

							if (diff_r >= radius_s_sq and diff_r < radius_e_sq and z_b>z_thresh)
							{	
								find_near(redshifts, z_b, red_num, tag_b);
								dist_source = distances[tag_b];

								backgal_sin_2phi = 2 * diff_ra*diff_dec / diff_r;
								backgal_cos_2phi = (diff_ra*diff_ra - diff_dec * diff_dec) / diff_r;

								shear_tan = -gh1 * backgal_cos_2phi - gh2 * backgal_sin_2phi;
								shear_tan_sig = -gh1_sig * backgal_cos_2phi - gh2_sig * backgal_sin_2phi;

								//shear_cros = gh1 * backgal_sin_2phi - gh2 * backgal_cos_2phi;
								//shear_cros_sig = gh1_sig * backgal_sin_2phi - gh2_sig * backgal_cos_2phi;

								crit_surf_density_com = dist_source / dist_len / (dist_source - dist_len) / (1 + z_f);

								foregal_data_shared[gal_id*foregal_data_shared_col + radi_id * radius_num] += 1;
								foregal_data_shared[gal_id*foregal_data_shared_col + radi_id * radius_num + 1] += crit_surf_density_com * shear_tan;
								foregal_data_shared[gal_id*foregal_data_shared_col + radi_id * radius_num + 2] += crit_surf_density_com * shear_tan_sig;

								//foregal_data_shared[gal_id*foregal_data_shared_col+2] += crit_surf_density_com * shear_cros;
								//foregal_data_shared[gal_id*foregal_data_shared_col+3] += crit_surf_density_com * shear_cros_sig;
							}
						}
					}
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);

			if (0 == rank)
			{
				crit_surf_density_com = 0;
				crit_surf_density_com_sig = 0;
				pair_count = 0;

				for (gal_id = 0; gal_id < backgal_num; gal_id++)
				{
					pair_count += 1;
				}

				for (gal_id = 0; gal_id < backgal_num; gal_id++)
				{
					crit_surf_density_com += foregal_data_shared[gal_id*foregal_data_shared_col + radi_id * radius_num + 1];
					crit_surf_density_com_sig += foregal_data_shared[gal_id*foregal_data_shared_col + radi_id * radius_num + 2];
				}

				delta_crit[radi_id] = crit_surf_density_com/pair_count;
				delta_crit[radi_id + radius_num] = crit_surf_density_com_sig/ pair_count;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		// finish the area_i
		MPI_Barrier(MPI_COMM_WORLD);
		if(0 == rank)
		{	
			sprintf(set_name, "/w_%d/result", area_id);
			if (1 == area_id)
			{
				write_h5(h5f_res_path, set_name, foregal_data_shared, foregal_num, foregal_data_shared_col, TRUE);
			}
			else
			{
				write_h5(h5f_res_path, set_name, foregal_data_shared, foregal_num, foregal_data_shared_col, FALSE);
			}

			sprintf(set_name, "/w_%d/mask", area_id);
			write_h5(h5f_res_path, set_name, backgal_mask, backgal_num, 1, FALSE);

			sprintf(set_name, "/w_%d/crit_surf_dens", area_id);
			write_h5(h5f_res_path, set_name, delta_crit, 2, radius_num, FALSE);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		st_end = clock();
		sprintf(log_infom, "RANK: %d. w_%d. Finish aread %d. Time: %.2f sec", rank, area_id, area_id, (st_end - st_start) / CLOCKS_PER_SEC);
		write_log(log_path, log_infom);
		if (0 == rank)
		{
			std::cout << log_infom << std::endl;
		}

		// free the memory	
		MPI_Win_free(&win_foregal);
		MPI_Win_free(&win_mask);
		MPI_Win_free(&win_shear);

		delete[] my_backgal_mask;

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

	sprintf(log_infom, "RANK: %d. Finish. Free the memory at last", rank);
	write_log(log_path, log_infom);
	std::cout << log_infom << std::endl;

	MPI_Finalize();
	return 0;
}
/*
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
	std::cout << rank << " " << radi_id << " " << radius_s << radius_e << std::endl;
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

time_label[3] = clock();
if (0 == process_per)
{
	sprintf(log_infom, "RANK: %d. w_%d. Looping foreground %.2f%%. Time used: %.2f sec", rank, area_id, per_n, (time_label[3] - time_label[2]) / CLOCKS_PER_SEC);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
	}
}*/