#include<FQlib.h>
#include<mpi.h>
#include<vector>

#define max_data_col 40
#define foregal_data_col 6
#define grid_data_col 5
#define backgal_data_col 21
#define mg_bin_num 8
#define vec_data_col 8
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
#endif

	int i, j, k, temp;
	char parent_path[250], data_path[250],result_path[250], log_path[250], h5f_path_grid[250], h5f_path_fore[250], h5f_res_path[250], temp_path[300];
	char set_name[50], set_name_2[50], attrs_name[80], log_infom[300];
	char foreground_name[50];
	char cata_name[30];

	int radius_num;
	double radius_s, radius_e, radius_e_sq;
	double *radius_bin;
	
	// the controllers
	int area_id = atoi(argv[1]);
	int radius_label = atoi(argv[2]);
	strcpy(foreground_name, argv[3]);

	// radius bin
	radius_num = 13;
	radius_bin = new double[radius_num + 1]{};
	log_bin(0.04, 15, radius_num + 1, radius_bin);

	   
	//sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/");
	//sprintf(h5f_res_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/result/%s/fourier/w_%d/radius_%d.hdf5", foreground_name, area_id, radius_label);
	//sprintf(log_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/log/ggl_log_%d.dat", rank);
	//sprintf(h5f_path_fore, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/foreground/%s/w_%d.hdf5", foreground_name, area_id);

	sprintf(parent_path, "/mnt/perc/hklee/CFHT/gg_lensing/");
	sprintf(data_path, "%sdata/", parent_path);

	sprintf(cata_name, "fourier_cata_old");
	
	sprintf(result_path, "%sresult/%s/%s/", parent_path, foreground_name, cata_name);
	sprintf(h5f_path_grid, "%s%s/fourier_cata_grid.hdf5", data_path, cata_name);

	sprintf(h5f_res_path, "%sw_%d/radius_%d.hdf5", result_path, area_id, radius_label);
	sprintf(h5f_path_fore, "%sforeground/%s/w_%d.hdf5", data_path, foreground_name, area_id);

	sprintf(log_path, "%slog/ggl_log_%d.dat", parent_path, rank);
	
	// be careful with the boundary of the guess of critical density and shear 
	double gh_crit_step = 0.001;
	double gh_crit_left = -0.32 + 0.02 * radius_label;
	double gh_crit_right = fabs(gh_crit_left);
	int gh_crit_num = int(gh_crit_right * 2 / gh_crit_step) + 1;
	double *gh_crit = new double[gh_crit_num];
	for (i = 0; i < gh_crit_num; i++)
	{
		gh_crit[i] = gh_crit_left + gh_crit_step * i;
	}

	double gh_step = 0.0001;
	double gh_left = -0.04 + 0.003 * radius_label;
	double gh_right = fabs(gh_left);
	int gh_num = int(gh_right * 2 / gh_step) + 1;
	double *gh = new double[gh_num];
	for (i = 0; i < gh_num; i++)
	{
		gh[i] = gh_left + gh_step * i;
	}

	//double gh_crit_step[13]{ 0.005, 0.005, 0.004, 0.002, 0.002, 0.001, 0.001, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002 };
	//double gh_crit_left[13]{-0.4,  -0.25, -0.2, -0.1, -0.05, -0.05,  -0.03, -0.01, -0.01, -0.01, -0.01,  -0.01,  -0.01};
	//double gh_crit_right[13]{0.4,  0.25,  0.2,  0.15,    0.1,     0.1,  0.08,   0.05,  0.05,   0.02,  0.02,  0.012,  0.012};
	//int gh_crit_num = int((gh_crit_right[radius_label] -gh_crit_left[radius_label])/ gh_crit_step[radius_label])+1;
	//double *gh_crit = new double[gh_crit_num];
	//for (i = 0; i < gh_crit_num; i++)
	//{
	//	gh_crit[i] = gh_crit_left[radius_label] + gh_crit_step[radius_label] * i;
	//}

	//double gh_step[13]{ 0.0005, 0.0005, 0.0004, 0.0002, 0.0002, 0.0001, 0.0001, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002 };
	//double gh_left[13]{ -0.04,  -0.025, -0.02, -0.01, -0.005, -0.005,  -0.003, -0.001, -0.001, -0.001, -0.001,  -0.001,  -0.001 };
	//double gh_right[13]{ 0.04,  0.025,  0.02,  0.015,    0.01,     0.01,  0.008,   0.005,  0.005,   0.002,  0.002,  0.0012,  0.0012 };
	//int gh_num = int((gh_right[radius_label] - gh_left[radius_label]) / gh_step[radius_label]) + 1;
	//double *gh = new double[gh_num];
	//for (i = 0; i < gh_crit_num; i++)
	//{
	//	gh[i] = gh_left[radius_label] + gh_step[radius_label] * i;
	//}

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
	double dist_len_integ, dist_source_integ, dist_len_integ_coeff;
	double coeff, coeff_inv, coeff_rad_dist;
	double crit_surf_density_com, crit_surf_density_com_integ;


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
	int mg_bin_num2 = mg_bin_num / 2;

	// the bin of G1(2) for shear estimation
	double *mg1_bin = new double[mg_bin_num + 1];
	double *mg2_bin = new double[mg_bin_num + 1];

	// chi of the signal from the all areas
	MY_INT *chi_tan_shared, *chi_cross_shared, *chi_crit_tan_shared, *chi_crit_cross_shared, *pair_count_shared;
	// chi square of the signal of each thread in each areas
	double *my_chi_tan = new double[mg_bin_num*gh_num];
	double *my_chi_cross = new double[mg_bin_num*gh_num];
	double *my_chi_crit_tan = new double[mg_bin_num*gh_crit_num];
	double *my_chi_crit_cross = new double[mg_bin_num*gh_crit_num];
	double mg_t, mg_x;
	int chi_bin_label;
	
	int grid_num, grid_ny, grid_nx;
	int row, col, bin_label;
	int block_id, block_s, block_e, ib;
	double block_scale[1];
	int *block_mask;
	int ra_bin_num, dec_bin_num;


	int nib_id = 0, bs_id = 1, be_id = 2, bdy_id = 3, bdx_id = 4;
	int z_id = 5, dist_id = 6, dist_integ_id = 7;
	int ra_id = 8, dec_id = 9, cos_dec_id = 10;
	int mg1_id = 11, mg2_id = 12, mn_id = 13, mu_id = 14, mv_id = 15;
	int zmin_lb = 16, zmax_lb = 17, odds_lb = 18;
	int ra_bin_id = 19, dec_bin_id = 20;


	int shape[2];

	char *names[max_data_col];
	//backgal_data_col includes the data to "DEC_BIN",
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
	sprintf(names[dist_integ_id], "DISTANCE_INTEG");
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
	
	sprintf(names[ra_bin_id], "RA_bin");
	sprintf(names[dec_bin_id], "DEC_bin");



	sprintf(log_infom, "RANK: %d. Start area: w_%d, radius: %d", rank, area_id, radius_label);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		std::cout << log_infom << std::endl;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// the shared buffer for the total chi square of the signal
#if ! defined(SMALL_CATA)

	MPI_Win win_chi_tan_total, win_chi_cross_total, win_chi_crit_tan_total, win_chi_crit_cross_total;
	MPI_Aint size_chi_tan, size_chi_cross, size_chi_crit_tan, size_chi_crit_cross;

	// [chi_tan, chi_cross]
	size_chi_tan =  mg_bin_num*gh_num;
	size_chi_cross = mg_bin_num*gh_num;
	size_chi_crit_tan = mg_bin_num * gh_crit_num;
	size_chi_crit_cross = mg_bin_num * gh_crit_num;
#endif
	MPI_Win win_pair_count;
	MPI_Aint  size_pair_count;

	size_pair_count = numprocs;

	if (0 == rank)
	{
#if ! defined( SMALL_CATA)
		// for the chi square of tangential shear
		MPI_Win_allocate_shared(size_chi_tan * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_tan_shared, &win_chi_tan_total);
		MPI_Win_allocate_shared(size_chi_cross * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_cross_shared, &win_chi_cross_total);
		// for the chi square of shear*critical_surface_density
		MPI_Win_allocate_shared(size_chi_crit_tan * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_crit_tan_shared, &win_chi_crit_tan_total);
		MPI_Win_allocate_shared(size_chi_crit_cross * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_crit_cross_shared, &win_chi_crit_cross_total);
#endif
		MPI_Win_allocate_shared(size_pair_count * sizeof(MY_INT), sizeof(MY_INT), MPI_INFO_NULL, MPI_COMM_WORLD, &pair_count_shared, &win_pair_count);
	}
	else
	{
		int dispu_total;
#if ! defined( SMALL_CATA)
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_tan_shared, &win_chi_tan_total);
		MPI_Win_shared_query(win_chi_tan_total, 0, &size_chi_tan, &dispu_total, &chi_tan_shared);

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_cross_shared, &win_chi_cross_total);
		MPI_Win_shared_query(win_chi_cross_total, 0, &size_chi_cross, &dispu_total, &chi_cross_shared);

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_crit_tan_shared, &win_chi_crit_tan_total);
		MPI_Win_shared_query(win_chi_crit_tan_total, 0, &size_chi_crit_tan, &dispu_total, &chi_crit_tan_shared);

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chi_crit_cross_shared, &win_chi_crit_cross_total);
		MPI_Win_shared_query(win_chi_crit_cross_total, 0, &size_chi_crit_cross, &dispu_total, &chi_crit_cross_shared);
#endif
		MPI_Win_allocate_shared(0, sizeof(MY_INT), MPI_INFO_NULL, MPI_COMM_WORLD, &pair_count_shared, &win_pair_count);
		MPI_Win_shared_query(win_pair_count, 0, &size_pair_count, &dispu_total, &pair_count_shared);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// initialization
	if (0 == rank)
	{	
#if ! defined( SMALL_CATA)
		initialize_arr(chi_tan_shared, mg_bin_num*gh_num, 0);
		initialize_arr(chi_cross_shared, mg_bin_num*gh_num, 0);
		initialize_arr(chi_crit_tan_shared, mg_bin_num*gh_crit_num, 0);
		initialize_arr(chi_crit_cross_shared, mg_bin_num*gh_crit_num, 0);
#endif
		initialize_arr(pair_count_shared, numprocs, 0);
	}
	initialize_arr(my_chi_tan, mg_bin_num*gh_num, 0);
	initialize_arr(my_chi_cross, mg_bin_num*gh_num, 0);
	initialize_arr(my_chi_crit_tan, mg_bin_num*gh_crit_num, 0);
	initialize_arr(my_chi_crit_cross, mg_bin_num*gh_crit_num, 0);
	MPI_Barrier(MPI_COMM_WORLD);
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////



	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// read foreground information
	// Z, DISTANCE, DISTANCE_INTEG, RA, DEC, COS_DEC
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

	// Z, Z_MIN, Z_MAX, ODDS,  DISTANCE, DISTANCE_INTEG,  RA, DEC,  G1, G2, N, U, V,  num_in_block,  block_start, block_end, 
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
		backgal_data[i] = new double[temp]{};
		read_h5(h5f_path_grid, set_name, backgal_data[i]);
	}
	//set_bin(backgal_data[mg1_id], backgal_num, mg_bin, mg_bin_num, 1000, 50000);
	sprintf(set_name, "/w_1/mg1_bin");
	read_h5(h5f_path_grid, set_name, mg1_bin);
	sprintf(set_name, "/w_1/mg2_bin");
	read_h5(h5f_path_grid, set_name, mg2_bin);

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

	coeff = 180 / Pi;
	coeff_inv = Pi / 180;

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
		dist_len_integ = foregal_data[dist_integ_id][gal_id];
		dist_len_integ_coeff = 1. / dist_len_integ / (1 + z_f);

		// the max searching radius depend on the redshift of lens 
		radius_e = radius_bin[radius_label + 1] * coeff / dist_len / foregal_data[cos_dec_id][gal_id] * 1.5; // degree

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
					z_b = backgal_data[z_id][ib];
					//if (backgal_data[z_id][ib] >= z_thresh)
					if (z_b >= z_thresh and z_b > z_b_sig95 and z_b_odds> 0.5)
					{
						ra_b = backgal_data[ra_id][ib];
						dec_b = backgal_data[dec_id][ib];

						// comoving distance
						dist_source = backgal_data[dist_id][ib];
						// the integrate part of the comoving distance
						dist_source_integ = backgal_data[dist_integ_id][ib];

						// times cos(dec) due to the different length the arc corresponding to the same delta R.A. at different Dec
						diff_ra = (ra_b - ra_f)*foregal_data[cos_dec_id][gal_id];
						diff_dec = dec_b - dec_f;
						diff_theta_sq = diff_ra * diff_ra + diff_dec * diff_dec; // degree^2

						// the seperation in comving coordinate, 
						//diff_r = dist_len * sqrt(diff_theta_sq)*coeff_inv;
						
						separation(ra_b, dec_b, ra_f, dec_f, diff_theta);
						diff_r = dist_len * diff_theta;
						//diff_r = dist_len * sin(diff_theta);
						//std::cout << radius_bin[radius_label] << " " << diff_r << " " << radius_bin[radius_label + 1] << " " << diff_theta << std::endl;
						if (radius_bin[radius_label] <= diff_r and diff_r < radius_bin[radius_label + 1])
						{
							pair_count_shared[rank] += 1;

							crit_surf_density_com = dist_source / (dist_source - dist_len) *dist_len_coeff;
							crit_surf_density_com_integ = dist_source_integ / (dist_source_integ - dist_len_integ) *dist_len_integ_coeff;

							// rotation for shear calculation, see the NOTE of gg_lensing for the detials 
							backgal_sin_2phi = 2 * diff_ra*diff_dec / diff_theta_sq;
							backgal_cos_2phi = (diff_dec - diff_ra)*(diff_ra + diff_dec) / diff_theta_sq;

							backgal_sin_4phi = 2 * backgal_sin_2phi * backgal_cos_2phi;
							backgal_cos_4phi = (backgal_cos_2phi + backgal_sin_2phi)*(backgal_cos_2phi - backgal_sin_2phi);
							
#if defined (SMALL_CATA)
							// G_t = (G_1 + i*G_2)*EXP(2i\phi) =  G_1 *cos2\phi - G_2*sin2\phi
							// the direction of R.A. is oppsite, actually,  G_t =  G_1 *cos2\phi + G_2*sin2\phi
							backgal_mg_tan = backgal_data[mg1_id][ib] * backgal_cos_2phi + backgal_data[mg2_id][ib] * backgal_sin_2phi;
							// the cross components
							backgal_mg_cross = backgal_data[mg1_id][ib] * backgal_sin_2phi - backgal_data[mg2_id][ib] * backgal_cos_2phi;
							// scalar
							backgal_mn_tan = backgal_data[mn_id][ib];
							// U_t = (U+i*V)*EXP(4i\phi)] = U*cos4\phi - V*sin\4phi
							backgal_mu_tan = backgal_data[mu_id][ib] * backgal_cos_4phi - backgal_data[mv_id][ib] * backgal_sin_4phi;

							data_cache.push_back(backgal_mg_tan);
							data_cache.push_back(backgal_mg_cross);
							data_cache.push_back(backgal_mn_tan);
							data_cache.push_back(backgal_mu_tan);
							data_cache.push_back(crit_surf_density_com);
							data_cache.push_back(crit_surf_density_com_integ);
							data_cache.push_back(diff_r);
							data_cache.push_back(z_b);

#else
							// G_t = (G_1 + i*G_2)*EXP(2i\phi) =  G_1 *cos2\phi - G_2*sin2\phi
							// the direction of R.A. is oppsite, actually,  G_t =  G_1 *cos2\phi + G_2*sin2\phi
							// \Sigma_crit *G_t(x)
							backgal_mg_tan = crit_surf_density_com * (backgal_data[mg1_id][ib] * backgal_cos_2phi + backgal_data[mg2_id][ib] * backgal_sin_2phi);
							// the cross components
							backgal_mg_cross = crit_surf_density_com * (backgal_data[mg1_id][ib] * backgal_sin_2phi - backgal_data[mg2_id][ib] * backgal_cos_2phi);
							// scalar
							backgal_mn_tan = backgal_data[mn_id][ib];
							// U_t = Re[(U+i*V)*EXP(-4i\phi)] = U*cos4\phi + V*sin\4phi
							backgal_mu_tan = backgal_data[mu_id][ib] * backgal_cos_4phi - backgal_data[mv_id][ib] * backgal_sin_4phi;
							backgal_mnu1_tan_c = crit_surf_density_com * (backgal_mn_tan + backgal_mu_tan);
							backgal_mnu2_tan_c = crit_surf_density_com * (backgal_mn_tan - backgal_mu_tan);
							backgal_mnu1_tan = backgal_mn_tan + backgal_mu_tan;
							backgal_mnu2_tan = backgal_mn_tan - backgal_mu_tan;

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
								my_chi_crit_tan[ig*mg_bin_num + chi_bin_label] +=1;

								mg_x = backgal_mg_cross - gh_crit[ig] * backgal_mnu2_tan;
								histogram_s(mg_x, mg2_bin, mg_bin_num, chi_bin_label);
								my_chi_crit_cross[ig*mg_bin_num + chi_bin_label]+=1;
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

			MPI_Gatherv(my_data_buf, num_of_thread[rank], MPI_DOUBLE, final_buf, num_of_thread, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);

			if (rank == 0)
			{
				// if mask == 1, the catalog can be used
				sprintf(set_name, "/mask");
				int *mask = new int[1];
				mask[0] = 1;
				write_h5(h5f_res_path, set_name, mask, 1, 1, TRUE);

				sprintf(set_name, "/pair_data");
				write_h5(h5f_res_path, set_name, final_buf, pair_count, vec_data_col, FALSE);


				sprintf(temp_path, "%sw_%d/radius_bin.hdf5", result_path, area_id);
				sprintf(set_name, "/radius_bin");
				write_h5(temp_path, set_name, radius_bin, radius_num + 1, 1, TRUE);

				delete[] final_buf;
			}
		}
		else
		{
			// if mask == 1, the catalog can be used
			sprintf(set_name, "/mask");
			int *mask = new int[1];
			mask[0] = 1;
			write_h5(h5f_res_path, set_name, mask, 1, 1, TRUE);

			sprintf(set_name, "/pair_data");
			write_h5(h5f_res_path, set_name, my_data_buf, pair_count, vec_data_col, FALSE);

			sprintf(temp_path, "%sw_%d/radius_bin.hdf5", result_path, area_id);
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
			// if mask == 1, the catalog can be used
			sprintf(set_name, "/mask");
			int *mask = new int[1];
			mask[0] = 0;
			write_h5(h5f_res_path, set_name, mask, 1, 1, TRUE);

			sprintf(temp_path, "%sw_%d/radius_bin.hdf5", result_path, area_id);
			sprintf(set_name, "/radius_bin");
			write_h5(temp_path, set_name, radius_bin, radius_num + 1, 1, TRUE);


		}
	}
#else
	if (rank == 0)
	{
		sprintf(log_infom, "RANK: %d. w_%d. %ld galaxies have been found in Radius [%.4f, %.4f].", rank, area_id, pair_count, radius_bin[radius_label], radius_bin[radius_label + 1]);
		std::cout << log_infom << std::endl;
	}
	// merge the PDF of signal estimator into the shared buffer
	for(i=0; i<numprocs; i++)
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

	// calculate the result and write to file
	if (0 == rank)
	{
		// the chi square for fitting shear
		long *chi_block = new long[mg_bin_num];
		double *chisq_tan = new double[gh_num];
		double *chisq_cross = new double[gh_num];
		double *chisq_crit_tan = new double[gh_crit_num];
		double *chisq_crit_cross = new double[gh_crit_num];
		double chi_temp, crit, crit_sig;
		double crit_result[4];
		// chi square of shear
		for (i = 0; i < gh_num; i++)
		{	
			// tangential 
			for (j = 0; j < mg_bin_num; j++)
			{
				chi_block[j] = chi_tan_shared[i*mg_bin_num + j];
			}
			try
			{
				cal_chisq_1d(chi_block, mg_bin_num, chi_temp);
			}
			catch (const char* msg)
			{
				std::cerr << "Rank: " << rank << ".  Tangential g_guess: " << i <<". "<< msg << std::endl;
				show_arr(chi_block, 1, mg_bin_num);
				chi_temp = 0;
			}
			chisq_tan[i] = chi_temp;

			// cross
			for (j = 0; j < mg_bin_num; j++)
			{
				chi_block[j] = chi_cross_shared[i*mg_bin_num + j];
			}
			try
			{
				cal_chisq_1d(chi_block, mg_bin_num, chi_temp);
			}
			catch (const char* msg)
			{
				std::cerr << "Rank: " << rank << ".  Cross g_guess: " << i << ". " << msg << std::endl;
				show_arr(chi_block, 1, mg_bin_num);
				chi_temp = 0;
			}
			chisq_cross[i] = chi_temp;
		}
		// chi square of shear*critical_surface_density
		for (i = 0; i < gh_crit_num; i++)
		{	// tangential 
			for (j = 0; j < mg_bin_num; j++)
			{
				chi_block[j] = chi_crit_tan_shared[i*mg_bin_num + j];
			}
			try
			{
				cal_chisq_1d(chi_block, mg_bin_num, chi_temp);
			}
			catch (const char* msg)
			{
				std::cerr << "Rank: " << rank << ".  Tangential g_crit_guess: " << i << ". " << msg << std::endl;
				show_arr(chi_block, 1, mg_bin_num);
				chi_temp = 0;
			}
			chisq_crit_tan[i] = chi_temp;
			// cross
			for (j = 0; j < mg_bin_num; j++)
			{
				chi_block[j] = chi_crit_cross_shared[i*mg_bin_num + j];
			}
			try
			{
				cal_chisq_1d(chi_block, mg_bin_num, chi_temp);
			}
			catch (const char* msg)
			{
				std::cerr << "Rank: " << rank << ".  Cross g_crit_guess: " << i << ". " << msg << std::endl;
				show_arr(chi_block, 1, mg_bin_num);
				chi_temp = 0;
			}
			chisq_crit_cross[i] = chi_temp;
		}

		// write the result to file
		sprintf(set_name, "/chi_tan");
		write_h5(h5f_res_path, set_name, chi_tan_shared, gh_num, mg_bin_num, TRUE);
		sprintf(set_name, "/chi_cross");
		write_h5(h5f_res_path, set_name, chi_cross_shared, gh_num, mg_bin_num, FALSE);

		sprintf(set_name, "/chi_crit_tan");
		write_h5(h5f_res_path, set_name, chi_crit_tan_shared, gh_crit_num, mg_bin_num, FALSE);
		sprintf(set_name, "/chi_crit_cross");
		write_h5(h5f_res_path, set_name, chi_crit_cross_shared, gh_crit_num, mg_bin_num, FALSE);

		sprintf(set_name, "/chisq_tan");
		write_h5(h5f_res_path, set_name, chisq_tan, gh_num, 1, FALSE);
		sprintf(set_name, "/chisq_cross");
		write_h5(h5f_res_path, set_name, chisq_cross, gh_num, 1, FALSE);

		sprintf(set_name, "/chisq_crit_tan");
		write_h5(h5f_res_path, set_name, chisq_crit_tan, gh_crit_num, 1, FALSE);
		sprintf(set_name, "/chisq_crit_cross");
		write_h5(h5f_res_path, set_name, chisq_crit_cross, gh_crit_num, 1, FALSE);

		sprintf(set_name, "/gh");
		write_h5(h5f_res_path, set_name, gh, gh_num, 1, FALSE);
		sprintf(set_name, "/gh_crit");
		write_h5(h5f_res_path, set_name, gh_crit, gh_crit_num, 1, FALSE);

		//// tangential
		//fit_shear(gh, chisq_tan, gh_num, crit, crit_sig, 50);
		//crit_result[0] = crit; 
		//crit_result[1] = crit_sig;
		//// cross
		//fit_shear(gh, chisq_cross, gh_num, crit, crit_sig, 50);
		//crit_result[2] = crit;
		//crit_result[3] = crit_sig;
		//sprintf(set_name, "/result");
		//write_h5(h5f_res_path, set_name, crit_result, 4, 1, FALSE);

		delete[] chi_block;
		delete[] chisq_cross;
		delete[] chisq_tan;
		delete[] chisq_crit_cross;
		delete[] chisq_crit_tan;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// free the memory	
	MPI_Win_free(&win_chi_tan_total);
	MPI_Win_free(&win_chi_cross_total);
	MPI_Win_free(&win_chi_crit_tan_total);
	MPI_Win_free(&win_chi_crit_cross_total);
#endif

	MPI_Win_free(&win_pair_count);
	
	sprintf(log_infom, "RANK: %d. Write file of Radius bin [%.4f,  %.4f] ---- end.", rank, radius_bin[radius_label], radius_bin[radius_label + 1]);
	write_log(log_path, log_infom);
	if (0 == rank)
	{
		char times[50];
		get_time(times, 50);
		std::cout << log_infom << std::endl;
		std::cout << h5f_res_path << std::endl;
		std::cout<< times << std::endl;
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
	delete[] mg1_bin;
	delete[] mg2_bin;
	delete[] gh;
	delete[] radius_bin;
	delete[] my_chi_tan;
	delete[] my_chi_cross;

	MPI_Finalize();
	return 0;
}
