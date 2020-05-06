#include<FQlib.h>
#include<hk_mpi.h>
#include<hk_iolib.h>
#define IMG_CHECK_LABEL 3
#define FLUX_PDF_UNI
#define EPSF


void arr_sqrt(double *arr_in, double *arr_out, const int length)
{
    for(int i=0;i<length;i++)
    {
        arr_out[i] = sqrt(arr_in[i]);
    }
}


void arr_add(double *arr1, const double*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] += arr2[i];
    }
}

void arr_add(double *arr1, const double*arr2,const double*arr3,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] = arr2[i] + arr3[i];
    }
}

void arr_deduct(double *result_buff, const double *arr1, const double*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        result_buff[i] = arr1[i] - arr2[i];
    }
}
void arr_deduct(double *result_buff, const double *arr1, const double*arr2, const double *arr3, const int length)
{
    for(int i=0;i<length;i++)
    {
        result_buff[i] = arr1[i] - arr2[i] - arr3[i];
    }
}
void arr_copy(double *arr1, const double*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] = arr2[i];
    }
}

void noise_subtraction_new(double *pow_arr1, const double*pow_arr2,const double *noise_pow_arr, const int length)
{
    for(int i=0;i<length;i++)
    {
        pow_arr1[i] = 2*pow_arr1[i] - pow_arr2[i];// + 2*noise_pow_arr[i];
    }
}
int main(int argc, char*argv[])
{	
	/* it is designed to run on the cluster (should save the memory) */
	/* the total chips should be divisible by the numprocs !!!! */
	/* loop the shear points, the chips will be distributed to threads */

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	fq_paras all_paras;

	char parent_path[300], chip_path[300], para_path[300], shear_path[300], result_path[300], log_path[300];
	char buffer[300], log_inform[250], set_name[50];

	
	
	int i, j, k, ib;
	int sss1, sss2, seed_pts, seed_n1, seed_n2, seed_step;
	int rotation;

	int num_p, size, shear_pairs, img_len;
    double max_radius;
	int total_chips, sub_chip_num, sub_data_row, total_data_row;
	int stamp_num, stamp_nx, shear_data_cols;		
	int row, chip_st, chip_ed, shear_id, psf_type, temp_s, detect_label;
	int seed_ini;

	double psf_scale, psf_thresh_scale, sig_level, psf_noise_sig, gal_noise_sig, flux_i, mag_i;
	double g1, g2, ts, te, t1, t2;
	double *g1t, *g2t;
	double psf_ellip, ellip_theta;
	double img_cent;
	double pts_step;
	double gal_fluxs[5]{3000, 6000, 9000, 16000, 800000};
	int flux_tag;

	double theta;

	strcpy(parent_path, argv[1]);
	seed_ini = atoi(argv[2]);
	flux_tag = atoi(argv[3]);
	total_chips = atoi(argv[4]);
	size = atoi(argv[5]);
	shear_pairs = atoi(argv[6]);
	pts_step = 1;//atof(argv[1]);

	rotation = 0;
	num_p = 30;
	max_radius= 7;
	stamp_num = 10000;
	shear_data_cols = 5;

	psf_type = 2;
	psf_scale = 4;//psf_scales[psf_scale_id];
	psf_ellip = 0.15;

	ellip_theta = 0;
	psf_thresh_scale = 2.;
	temp_s = rank;
	sig_level = 1.5;
	psf_noise_sig = 0;
    gal_noise_sig = 60;

    img_len = size*size;
	img_cent = size*0.5 - 0.5;
    stamp_nx = 100;

	all_paras.stamp_size = size;
	all_paras.img_x = size;
	all_paras.img_y = size;
	// because the max half light radius of the galsim source is 5.5 pixels
	all_paras.max_distance = max_radius; 

	sprintf(log_path, "%s/logs/%02d.dat", parent_path, rank);

#ifdef IMG_CHECK_LABEL
    double *big_img_noise_free = new double[stamp_nx*stamp_nx*img_len]{};
    double *big_img_noise_residual = new double[stamp_nx*stamp_nx*img_len]{};
    double *big_img_cross_term = new double[stamp_nx*stamp_nx*img_len]{};
	double *big_img_noisy = new double[stamp_nx*stamp_nx*img_len]{};
	
#endif

	double *point = new double[2 * num_p]{};
	double *point_r = new double[2 * num_p]{};


	double *gal = new double[img_len]{};
	double *pgal = new double[img_len]{};
	double *pgal_dn = new double[img_len]{};
	
	double *gal_noisy = new double[img_len]{};
	double *pgal_noisy = new double[img_len]{};

	double *gal_noisy_a = new double[img_len]{};
	double *pgal_noisy_a = new double[img_len]{};

	double *pgal_cross_term = new double[img_len]{};
	double *pgal_cross_term_est = new double[img_len]{};

	double *psf = new double[img_len]{};
	double *ppsf = new double[img_len]{};
	double *ppsf_sqrt = new double[img_len]{};

	double *noise_1 = new double[img_len]{};
	double *pnoise_1 = new double[img_len]{};

    double *noise_2 = new double[img_len]{};
	double *pnoise_2 = new double[img_len]{};
    double *noise_pow_diff = new double[img_len]{};

	double *noise_3 = new double[img_len]{};
	double *pnoise_3 = new double[img_len]{};
	double *pnoise_cross = new double[img_len]{};

	// the shear estimators data matrix  
	double *total_data;	
	double *sub_noise_free_data;
    double *sub_noise_residual_data;	
	double *sub_cross_term_data;
	double *sub_cross_term_sqrt_data;
    double *sub_noisy_data;
	double *sub_cross_term_estimate;
	double *sub_noise_cross_term_estimate;
	double *sub_noise_residual_term_estimate;	


	double *total_flux, *sub_flux;
	int *scatter_count,*gather_count;

	// for scatterring the flux to each thread
	scatter_count = new int[numprocs]{};
	// for the gatherv when finish in each shear point
	gather_count = new int[numprocs]{};

	///////////////////// task distribution /////////////////////////////////////
	sub_chip_num = total_chips / numprocs;
	j = total_chips%numprocs;
	for(i=0;i<numprocs;i++)
	{
		scatter_count[i] = sub_chip_num;
	}
	for(i=0;i<j;i++)
	{	
		// the chip number
		scatter_count[i] += 1;
	}
	sub_chip_num = scatter_count[rank];
	// the start- & end-label of chip of each thread
	chip_st = 0;
	for(i=0;i<rank;i++)
	{
		chip_st += scatter_count[i];
	}
	chip_ed = chip_st+scatter_count[rank];
	// the final data from all the source in one shear point
	total_data_row = total_chips * stamp_num;
	// the sub-data from the source processed by each thread
	sub_data_row = sub_chip_num * stamp_num;
	for(i=0;i<numprocs;i++)
	{
		// the real count of galaxies for each thread
		scatter_count[i] *= stamp_num;
		// the real amount of data of each thread
		gather_count[i] = scatter_count[i]*shear_data_cols;
	}
	
	sub_noise_free_data = new double[gather_count[rank]]{};
	sub_noise_residual_data = new double[gather_count[rank]]{};
	sub_cross_term_data = new double[gather_count[rank]]{};
	sub_cross_term_sqrt_data = new double[gather_count[rank]]{};
	sub_noisy_data = new double[gather_count[rank]]{};
	sub_cross_term_estimate = new double[gather_count[rank]]{};
	sub_noise_cross_term_estimate = new double[gather_count[rank]]{};
	sub_noise_residual_term_estimate = new double[gather_count[rank]]{};

	double *total_theta;
	double *sub_theta = new double[scatter_count[rank]]{};


	// seed distribution, different thread gets different seed
	seed_step = 1;
	sss1 = 2*seed_step*shear_pairs*total_chips;
	seed_pts = sss1*rank + 1 + seed_ini;//35000;
	seed_n1 = sss1*rank + 1 + seed_ini*2;// 4001*(rotation+1);
	seed_n2 = sss1*rank + 1 + seed_ini*4;//2300*(rotation+1);

	if (0 == rank)
	{
		total_theta = new double[total_data_row]{};
		total_data = new double[total_data_row*shear_data_cols]{};
	}

	// read shear
	g1t = new double[shear_pairs]{};
	g2t = new double[shear_pairs]{};

	sprintf(shear_path,"%s/parameters/shear.hdf5", parent_path);
	sprintf(set_name,"/g1");
	read_h5(shear_path, set_name, g1t);
	sprintf(set_name,"/g2");
	read_h5(shear_path, set_name, g2t);

    
#ifdef EPSF
	create_psf_e(psf, psf_scale, size, img_cent, psf_ellip, ellip_theta, psf_type);
#else
	// create PSF
	create_psf(psf, psf_scale, size, img_cent, psf_type);
#endif
	pow_spec(psf, ppsf, size, size);
	get_psf_radius(ppsf, &all_paras, psf_thresh_scale);
	arr_sqrt(ppsf, ppsf_sqrt, img_len);
	// std::cout<<rank<<std::endl;
	// show_arr(scatter_count,1,numprocs);
	// std::cout<<"---------------------------------------------------------------------------"<<std::endl;
	
	if (0 == rank)
	{	
		std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		std::cout << parent_path << std::endl;
		std::cout<<"Shear num: "<<shear_pairs<<std::endl;
        std::cout<<"Rotation: "<<rotation<<std::endl;
		std::cout << "Total chip: " << total_chips<< ", Stamp size: " << size  << std::endl;
		std::cout << "Total cpus: " << numprocs << std::endl;
		std::cout <<"PSF Scale: "<<psf_scale<< "PSF THRESH: " << all_paras.psf_pow_thresh <<" PSF HLR: " << all_paras.psf_hlr << std::endl;
		std::cout <<"MAX RADIUS: "<< max_radius <<" ,Step: "<<pts_step<< ", SIG_LEVEL: " << sig_level <<"sigma"<< std::endl;
#ifdef EPSF
		sprintf(buffer, "!%s/epsf_%.2f.fits", parent_path,psf_scale);
#else
		sprintf(buffer, "!%s/psf_%.2f.fits", parent_path,psf_scale);
#endif
		write_fits(buffer,psf, size, size);
		std::cout<<"Gal Num of each thread: ";
		show_arr(scatter_count,1,numprocs);
		std::cout<<"---------------------------------------------------------------------------"<<std::endl<<std::endl;
	}
	for(i=0;i<10;i++)
	{
		if(i == rank)
		{
			std::cout<<rank<<" "<<chip_st<<" "<<chip_ed<<" "<<seed_pts<<" "<<seed_n1<<" "<<seed_n2<<" "<<scatter_count[i]<<" "<<gather_count[i]<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	
	// loop the shear points
	for (shear_id = 0; shear_id < shear_pairs; shear_id++)
	{
		ts = clock();

		g1 = g1t[shear_id];
		g2 = g2t[shear_id];

		sprintf(log_inform, "size: %d, total chips: %d (%d cpus),  point num: %d , noise sigma: %.2f ", size, total_chips, numprocs, num_p, gal_noise_sig);
		write_log(log_path, log_inform);
		sprintf(log_inform, "PSF scale: %.2f, max radius: %.2f, Step: %.4f", psf_scale, max_radius, pts_step);
		write_log(log_path, log_inform);
		sprintf(log_inform, "RANK: %03d, SHEAR %02d: my chips: %d - %d", rank, shear_id, chip_st, chip_ed);
		write_log(log_path, log_inform);

		// rank 0 reads the total flux array and scatters to each thread
		if(rank == 0)
		{
			std::cout<<"---------------------------------------------------------------------------"<<std::endl;
			std::cout<<"g1: "<<g1<<" g2: "<<g2<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);


		// initialize_arr(point, num_p * 2, 0);
		// gsl_initialize(shear_id+seed_ini,0);
		// create_points(point, num_p, max_radius, pts_step, rng0);
		// for(i=0;i<5;i++)
		// {
		// 	if(rank == i)
		// 	{std::cout<<"------------------"<<std::endl;show_arr(point,2,num_p);std::cout<<"------------------"<<std::endl;}
		// 	MPI_Barrier(MPI_COMM_WORLD);
		// }
		// gsl_free(0);
		

		// loop the chips
		for (i = chip_st; i < chip_ed; i++)
		{
			t1 = clock();

			sprintf(log_inform, "RANK: %03d, SHEAR %02d:, chip: %05d, start. seed: %d, %d, %d", rank,shear_id, i, seed_pts, seed_n1, seed_n2);
			write_log(log_path, log_inform);
			if (rank == 0)
			{
				std::cout << log_inform << std::endl;
			}

			// initialize GSL
			gsl_initialize(seed_pts,0);
			gsl_initialize(seed_n1,1);
			gsl_initialize(seed_n2,2);
			seed_pts += seed_step;
			seed_n1 += seed_step;
			seed_n2 += seed_step;

			row = (i-chip_st)*stamp_num*shear_data_cols;

			// loop the stamps
			for (j = 0; j < stamp_num; j++)
			{	

				flux_i = gal_fluxs[flux_tag]/ num_p;				
				
				initialize_arr(pgal, size*size, 0);
				initialize_arr(point, num_p * 2, 0);
				initialize_arr(noise_1, img_len, 0);
				initialize_arr(pnoise_1, img_len, 0);
				initialize_arr(noise_2, img_len, 0);
				initialize_arr(pnoise_2, img_len, 0);

				initialize_arr(pgal_cross_term, img_len, 0);
				initialize_arr(pgal_noisy, img_len, 0);

				addnoise(noise_1, img_len, gal_noise_sig, rng1);
				pow_spec(noise_1, pnoise_1, size, size);
				addnoise(noise_2, img_len, gal_noise_sig, rng1);
				pow_spec(noise_2, pnoise_2, size, size);

				// rand_uniform(0, 2*Pi, theta, rng2);
				// coord_rotation(point, num_p, theta, point_r);
				// sub_theta[row/shear_data_cols+j] = theta;

				create_points(point, num_p, max_radius, pts_step, rng1);

#ifdef EPSF
				convolve_e(gal, point, flux_i, size, img_cent, num_p, rotation, psf_scale, g1, g2, psf_type, psf_ellip, ellip_theta);
#else
				// convolve(gal, point, flux_i, size, img_cent, num_p, rotation, psf_scale, g1, g2, psf_type);
				convolve(gal, point_r, flux_i, size, img_cent, num_p, rotation, psf_scale, g1, g2, psf_type);
#endif
				// noise free
				pow_spec(gal, pgal, size, size);
				// noisy image
				arr_add(gal_noisy, gal, noise_1,img_len);
				pow_spec(gal_noisy, pgal_noisy, size, size);
				arr_deduct(pgal_dn, pgal_noisy, pnoise_2, img_len);
				// noise residual
				arr_deduct(noise_pow_diff, pnoise_1, pnoise_2, img_len);
				// cross term
                arr_deduct(pgal_cross_term, pgal_noisy, pgal, pnoise_1, img_len);
				// cross term estimate
				arr_add(gal_noisy_a, gal_noisy, noise_2, img_len);
				pow_spec(gal_noisy_a, pgal_noisy_a, size, size);
				arr_deduct(pgal_cross_term_est, pgal_noisy_a, pgal_noisy, pnoise_2, img_len);


				///////////////// Noise free /////////////////////////
				shear_est(pgal, ppsf, &all_paras);
				sub_noise_free_data[row + j * shear_data_cols] = all_paras.n1;
				sub_noise_free_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_noise_free_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_noise_free_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_noise_free_data[row + j * shear_data_cols + 4] = all_paras.dv;

				/////////////////////////// noisy image /////////////////////////////////
 				shear_est(pgal_dn, ppsf, &all_paras);
				sub_noisy_data[row + j * shear_data_cols] = all_paras.n1;
				sub_noisy_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_noisy_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_noisy_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_noisy_data[row + j * shear_data_cols + 4] = all_paras.dv;

 				////////////////// noise-residual image //////////////////////
				shear_est(noise_pow_diff, ppsf, &all_paras);
				sub_noise_residual_data[row + j * shear_data_cols] = all_paras.n1;
				sub_noise_residual_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_noise_residual_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_noise_residual_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_noise_residual_data[row + j * shear_data_cols + 4] = all_paras.dv;
				
				////////////////// cross-term image //////////////////////
				shear_est(pgal_cross_term, ppsf, &all_paras);
				sub_cross_term_data[row + j * shear_data_cols] = all_paras.n1;
				sub_cross_term_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_cross_term_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_cross_term_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_cross_term_data[row + j * shear_data_cols + 4] = all_paras.dv;

				////////////////// cross-term divided by |P_psf| image //////////////////////
				shear_est(pgal_cross_term, ppsf_sqrt, &all_paras);
				sub_cross_term_sqrt_data[row + j * shear_data_cols] = all_paras.n1;
				sub_cross_term_sqrt_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_cross_term_sqrt_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_cross_term_sqrt_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_cross_term_sqrt_data[row + j * shear_data_cols + 4] = all_paras.dv;

				////////////////// cross-term-est image //////////////////////
				shear_est(pgal_cross_term_est, ppsf, &all_paras);				
				sub_cross_term_estimate[row + j * shear_data_cols] = all_paras.n1;
				sub_cross_term_estimate[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_cross_term_estimate[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_cross_term_estimate[row + j * shear_data_cols + 3] = all_paras.du;
				sub_cross_term_estimate[row + j * shear_data_cols + 4] = all_paras.dv;

               

#ifdef IMG_CHECK_LABEL
				if(i <= IMG_CHECK_LABEL)
				{
					stack(big_img_noise_free, gal, j, size, stamp_nx, stamp_nx);
					stack(big_img_noise_residual, noise_pow_diff, j, size, stamp_nx, stamp_nx);
					stack(big_img_cross_term, pgal_cross_term, j, size, stamp_nx, stamp_nx);
					stack(big_img_noisy, gal_noisy, j, size, stamp_nx, stamp_nx);
				}
#endif
		
			}
			gsl_free(0);
			gsl_free(1);
			gsl_free(2);


#ifdef IMG_CHECK_LABEL
			if(i <= IMG_CHECK_LABEL)
			{
				sprintf(chip_path, "!%s/%d/gal_chip_%05d_noise_free.fits", parent_path, shear_id, i);
				write_fits(chip_path, big_img_noise_free, stamp_nx*size, stamp_nx*size);
				sprintf(chip_path, "!%s/%d/gal_chip_%05d_noise_pow_residual.fits", parent_path, shear_id, i);
				write_fits(chip_path, big_img_noise_residual, stamp_nx*size, stamp_nx*size);
                sprintf(chip_path, "!%s/%d/gal_chip_%05d_cross_term.fits", parent_path, shear_id, i);
				write_fits(chip_path, big_img_cross_term, stamp_nx*size, stamp_nx*size);
				sprintf(chip_path, "!%s/%d/gal_chip_%05d_noisy.fits", parent_path, shear_id, i);
				write_fits(chip_path, big_img_noisy, stamp_nx*size, stamp_nx*size);
			}
#endif
			t2 = clock();
			sprintf(log_inform, "RANK: %03d, SHEAR %02d: chip: %05d, done in %.2f s.", rank, shear_id, i, (t2 - t1) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			if (rank == 0)
			{
				std::cout << log_inform << std::endl;
			}
		}	
		// finish the chip loop
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_noise_free_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
			sprintf(set_name, "/data");
#ifdef EPSF
			sprintf(result_path, "%s/data/data_noise_free_epsf_%d.hdf5", parent_path,shear_id);
#else
			sprintf(result_path, "%s/data/data_noise_free_%d.hdf5", parent_path, shear_id);
#endif
			write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols, true);
		}
		MPI_Barrier(MPI_COMM_WORLD);


		my_Gatherv(sub_noisy_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
#ifdef EPSF
			sprintf(result_path, "%s/data/data_noisy_cpp_epsf_%d.hdf5", parent_path, shear_id);
#else
			sprintf(result_path, "%s/data/data_noisy_cpp_%d.hdf5", parent_path,shear_id);
#endif
			write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols, true);
			std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_noise_residual_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
#ifdef EPSF
			sprintf(result_path, "%s/data/data_noise_residual_epsf_%d.hdf5", parent_path, shear_id);
#else
			sprintf(result_path, "%s/data/data_noise_residual_%d.hdf5", parent_path,shear_id);
#endif
			write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols, true);
		}
		MPI_Barrier(MPI_COMM_WORLD);


		my_Gatherv(sub_cross_term_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
#ifdef EPSF
			sprintf(result_path, "%s/data/data_cross_term_epsf_%d.hdf5", parent_path, shear_id);
#else
			sprintf(result_path, "%s/data/data_cross_term_%d.hdf5", parent_path,shear_id);
#endif
			write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols, true);
		}
		MPI_Barrier(MPI_COMM_WORLD);


		my_Gatherv(sub_cross_term_sqrt_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
#ifdef EPSF
			sprintf(result_path, "%s/data/data_cross_term_sqrt_epsf_%d.hdf5", parent_path, shear_id);
#else
			sprintf(result_path, "%s/data/data_cross_term_sqrt_%d.hdf5", parent_path,shear_id);
#endif
			write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols, true);
		}
		MPI_Barrier(MPI_COMM_WORLD);


		my_Gatherv(sub_cross_term_estimate, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
#ifdef EPSF
			sprintf(result_path, "%s/data/data_cross_term_est_epsf_%d.hdf5", parent_path, shear_id);
#else
			sprintf(result_path, "%s/data/data_cross_term_est_%d.hdf5", parent_path,shear_id);
#endif
			write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols, true);
			std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);


		// my_Gatherv(sub_theta, scatter_count, total_theta, numprocs, rank);
		// if (0 == rank)
		// {
		// 	sprintf(result_path, "%s/data/theta_%d.hdf5", parent_path,shear_id);
		// 	write_h5(result_path, set_name, total_theta, total_data_row, 1, true);
		// 	std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		// }
		// MPI_Barrier(MPI_COMM_WORLD);

		te = clock();
		if (rank == 0)
		{
			sprintf(buffer, "rank %d: done in %g %s/%d\n", rank, (te - ts) / CLOCKS_PER_SEC, parent_path, shear_id);
			std::cout << buffer<<std::endl;
		}
	}


	if (0 == rank)
	{	
		std::cout<<parent_path<<std::endl;
		std::cout << "FINISH ALL JOBS" << std::endl;
		delete[] total_data;
#ifdef FLUX_PDF
		delete[] total_flux;
#endif
	}

#ifdef IMG_CHECK_LABEL
	delete[] big_img_noise_free;
	delete[] big_img_noisy;
	delete[] big_img_cross_term;
	delete[] big_img_noise_residual;
#endif

	delete[] point;
	delete[] gal;
	delete[] pgal;
	delete[] psf;
	delete[] ppsf;
	delete[] noise_1;
	delete[] pnoise_1;
    delete[] noise_2;
	delete[] pnoise_2;
	delete[] sub_noise_free_data;
	delete[] sub_noise_residual_data;
	delete[] sub_cross_term_data;
	delete[] sub_noisy_data;
#ifdef FLUX_PDF
	delete[] sub_flux;
#endif
	delete[] scatter_count;
	delete[] gather_count;
	
	MPI_Finalize();
	return 0;
}
