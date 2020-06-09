#include<FQlib.h>
#include<hk_mpi.h>
#include<hk_iolib.h>
#define IMG_CHECK_LABEL 3
#define FLUX_PDF_UNI
#define EPSF

#define SAVE_MEM

//#define DATA_SEP

#ifdef SAVE_MEM
#define MY_FLOAT float
#else
#define MY_FLOAT double
#endif


void arr_sqrt(MY_FLOAT *arr_in, MY_FLOAT *arr_out, const int length)
{
    for(int i=0;i<length;i++)
    {
        arr_out[i] = sqrt(arr_in[i]);
    }
}


void arr_add(MY_FLOAT *arr1, const MY_FLOAT*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] += arr2[i];
    }
}

void arr_add(MY_FLOAT *arr1, const MY_FLOAT*arr2,const MY_FLOAT*arr3,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] = arr2[i] + arr3[i];
    }
}

void arr_deduct(MY_FLOAT *result_buff, const MY_FLOAT *arr1, const MY_FLOAT*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        result_buff[i] = arr1[i] - arr2[i];
    }
}
void arr_deduct(MY_FLOAT *result_buff, const MY_FLOAT *arr1, const MY_FLOAT*arr2, const MY_FLOAT *arr3, const int length)
{
    for(int i=0;i<length;i++)
    {
        result_buff[i] = arr1[i] - arr2[i] - arr3[i];
    }
}
void arr_copy(MY_FLOAT *arr1, const MY_FLOAT*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] = arr2[i];
    }
}

void noise_subtraction_new(MY_FLOAT *pow_arr1, const MY_FLOAT*pow_arr2,const MY_FLOAT *noise_pow_arr, const int length)
{
    for(int i=0;i<length;i++)
    {
        pow_arr1[i] = 2*pow_arr1[i] - pow_arr2[i];// + 2*noise_pow_arr[i];
    }
}
void data_sep(const MY_FLOAT *data_arr, MY_FLOAT **sub_data, const int data_row, const int data_col)
{
	int i, j;
	for(i=0;i<data_row;i++)
	{
		for(j=0;j<data_col;j++)
		{
			sub_data[j][i] = data_arr[i*data_col+j];
		}
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


#ifdef SAVE_MEM
	fq_paras_float all_paras;
#else
	fq_paras all_paras;
#endif

	char parent_path[300], chip_path[300], para_path[300], shear_path[300], result_path[300], log_path[300];
	char buffer[300], log_inform[250], set_name[50];

	
	
	int i, j, k, ib;
	int sss1, sss2, seed_pts, seed_n1, seed_n2, seed_step;
	int rotation;

	int num_p, size, shear_pairs, img_len;
    MY_FLOAT max_radius;
	int total_chips, sub_chip_num, sub_data_row, total_data_row;
	int stamp_num, stamp_nx, shear_data_cols;		
	int row, chip_st, chip_ed, shear_id, psf_type, temp_s, detect_label;
	int seed_ini;

	MY_FLOAT psf_scale, psf_thresh_scale, sig_level, psf_noise_sig, gal_noise_sig, flux_i, mag_i;
	MY_FLOAT g1, g2, ts, te, t1, t2;
	MY_FLOAT *g1t, *g2t;
	MY_FLOAT psf_ellip, ellip_theta;
	MY_FLOAT img_cent;
	MY_FLOAT pts_step;
	MY_FLOAT gal_fluxs[8]{4000, 6000, 9000, 13500, 20250, 30375, 64000, 96000};
	int flux_tag;

	MY_FLOAT theta;

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
	psf_ellip = 0.10;

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
    MY_FLOAT *big_img_noise_free = new MY_FLOAT[stamp_nx*stamp_nx*img_len]{};
    MY_FLOAT *big_img_noise_residual = new MY_FLOAT[stamp_nx*stamp_nx*img_len]{};
    MY_FLOAT *big_img_cross_term = new MY_FLOAT[stamp_nx*stamp_nx*img_len]{};
	MY_FLOAT *big_img_noisy = new MY_FLOAT[stamp_nx*stamp_nx*img_len]{};
	
#endif

	MY_FLOAT *point = new MY_FLOAT[2 * num_p]{};
	MY_FLOAT *point_r = new MY_FLOAT[2 * num_p]{};


	MY_FLOAT *gal = new MY_FLOAT[img_len]{};
	MY_FLOAT *pgal = new MY_FLOAT[img_len]{};
	MY_FLOAT *pgal_sqrt = new MY_FLOAT[img_len]{};
	MY_FLOAT *pgal_dn = new MY_FLOAT[img_len]{};
	
	MY_FLOAT *gal_noisy = new MY_FLOAT[img_len]{};
	MY_FLOAT *pgal_noisy = new MY_FLOAT[img_len]{};

	MY_FLOAT *gal_noisy_a = new MY_FLOAT[img_len]{};
	MY_FLOAT *pgal_noisy_a = new MY_FLOAT[img_len]{};

	MY_FLOAT *pgal_cross_term = new MY_FLOAT[img_len]{};
	MY_FLOAT *pgal_cross_term_est = new MY_FLOAT[img_len]{};
	MY_FLOAT *pgal_pure_cross_term_est = new MY_FLOAT[img_len]{};


	MY_FLOAT *psf = new MY_FLOAT[img_len]{};
	MY_FLOAT *ppsf = new MY_FLOAT[img_len]{};
	MY_FLOAT *ppsf_sqrt = new MY_FLOAT[img_len]{};

	MY_FLOAT *noise_1 = new MY_FLOAT[img_len]{};
	MY_FLOAT *pnoise_1 = new MY_FLOAT[img_len]{};

    MY_FLOAT *noise_2 = new MY_FLOAT[img_len]{};
	MY_FLOAT *pnoise_2 = new MY_FLOAT[img_len]{};

	MY_FLOAT *noise_3 = new MY_FLOAT[img_len]{};
	MY_FLOAT *pnoise_3 = new MY_FLOAT[img_len]{};

	MY_FLOAT *noise_4 = new MY_FLOAT[img_len]{};
	MY_FLOAT *pnoise_4 = new MY_FLOAT[img_len]{};
	MY_FLOAT *temp = new MY_FLOAT[img_len]{};
	MY_FLOAT *noise_pow_diff = new MY_FLOAT[img_len]{};
	MY_FLOAT *noise_cross = new MY_FLOAT[img_len]{};
	MY_FLOAT *pnoise_cross = new MY_FLOAT[img_len]{};
	MY_FLOAT *pnoise_cross_est = new MY_FLOAT[img_len]{};


	// the shear estimators data matrix  
	MY_FLOAT *total_data;	
	MY_FLOAT *sub_noise_free_data;
	MY_FLOAT *sub_noise_free_sqrt_data;
    MY_FLOAT *sub_noise_residual_data;	
	MY_FLOAT *sub_cross_term_data;
	MY_FLOAT *sub_pure_cross_term_estimate;
	MY_FLOAT *sub_cross_term_sqrt_data;
    MY_FLOAT *sub_noisy_data;
	MY_FLOAT *sub_cross_term_estimate;
	MY_FLOAT *sub_noise_cross_term_estimate;
	MY_FLOAT *sub_noise_cross_term;	
	MY_FLOAT *mg_data[5];
	char *mg_name[5];

	MY_FLOAT *total_flux, *sub_flux;
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
	
	sub_noise_free_data = new MY_FLOAT[gather_count[rank]]{};
	// sub_noise_free_sqrt_data = new MY_FLOAT[gather_count[rank]]{};
	sub_noise_residual_data = new MY_FLOAT[gather_count[rank]]{};
	sub_cross_term_data = new MY_FLOAT[gather_count[rank]]{};
	// sub_cross_term_sqrt_data = new MY_FLOAT[gather_count[rank]]{};
	sub_noisy_data = new MY_FLOAT[gather_count[rank]]{};
	sub_cross_term_estimate = new MY_FLOAT[gather_count[rank]]{};
	sub_pure_cross_term_estimate = new MY_FLOAT[gather_count[rank]]{};
	sub_noise_cross_term = new MY_FLOAT[gather_count[rank]]{};
	sub_noise_cross_term_estimate = new MY_FLOAT[gather_count[rank]]{};



	// seed distribution, different thread gets different seed
	seed_step = 1;
	sss1 = 2*seed_step*shear_pairs*total_chips;
	seed_pts = sss1*rank + 1 + seed_ini;//35000;
	seed_n1 = sss1*rank + 1 + seed_ini*2;// 4001*(rotation+1);
	seed_n2 = sss1*rank + 1 + seed_ini*4;//2300*(rotation+1);

	if (0 == rank)
	{
		total_data = new MY_FLOAT[total_data_row*shear_data_cols]{};
		
		for(i=0; i<shear_data_cols; i++)
		{
			mg_name[i] = new char[50];
			mg_data[i] = new MY_FLOAT[total_data_row];
		}
		sprintf(mg_name[0], "/mg1");
		sprintf(mg_name[1], "/mg2");
		sprintf(mg_name[2], "/mn");
		sprintf(mg_name[3], "/mu");
		sprintf(mg_name[4], "/mv");
	}

	// read shear
	g1t = new MY_FLOAT[shear_pairs]{};
	g2t = new MY_FLOAT[shear_pairs]{};

	// sprintf(shear_path,"%s/parameters/shear.hdf5", parent_path);
	// sprintf(set_name,"/g1");
	// read_h5(shear_path, set_name, g1t);
	// sprintf(set_name,"/g2");
	// read_h5(shear_path, set_name, g2t);

    
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

		g1 = -0.04;//g1t[shear_id];
		g2 = 0.033;//g2t[shear_id];

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
				initialize_arr(noise_3, img_len, 0);
				initialize_arr(pnoise_3, img_len, 0);
				initialize_arr(noise_4, img_len, 0);
				initialize_arr(pnoise_4, img_len, 0);

				initialize_arr(pgal_cross_term, img_len, 0);
				initialize_arr(pgal_noisy, img_len, 0);

				addnoise(noise_1, img_len, gal_noise_sig, rng1);
				pow_spec(noise_1, pnoise_1, size, size);
				addnoise(noise_2, img_len, gal_noise_sig, rng1);
				pow_spec(noise_2, pnoise_2, size, size);

				addnoise(noise_3, img_len, gal_noise_sig, rng2);
				pow_spec(noise_3, pnoise_3, size, size);
				addnoise(noise_4, img_len, gal_noise_sig, rng2);
				pow_spec(noise_4, pnoise_4, size, size);

				// rand_uniform(0, 2*Pi, theta, rng2);
				// coord_rotation(point, num_p, theta, point_r);
				// sub_theta[row/shear_data_cols+j] = theta;

				create_points(point, num_p, max_radius, pts_step, rng0);

#ifdef EPSF
				convolve_e(point,num_p,flux_i, g1, g2, gal, size, img_cent, psf_scale,psf_type,psf_ellip, ellip_theta);
#else
				convolve(point,num_p,flux_i, g1, g2, gal, size, img_cent, psf_scale,psf_type);

#endif
				// noise free
				pow_spec(gal, pgal, size, size);
				arr_sqrt(pgal, pgal_sqrt, img_len);
				// noisy image
				arr_add(gal_noisy, gal, noise_1,img_len);
				pow_spec(gal_noisy, pgal_noisy, size, size);
				arr_deduct(pgal_dn, pgal_noisy, pnoise_2, img_len);
				// noise residual
				arr_deduct(noise_pow_diff, pnoise_1, pnoise_2, img_len);
				// noise cross term
				arr_add(noise_cross,noise_1, noise_2, img_len);
				pow_spec(noise_cross, temp, size, size);
				arr_deduct(pnoise_cross, temp, pnoise_1, pnoise_2,  img_len);
				// noise cross term estimate
				arr_add(noise_cross, noise_3, noise_4, img_len);
				pow_spec(noise_cross, temp, size, size);
				arr_deduct(pnoise_cross_est, temp, pnoise_3, pnoise_4, img_len);
				// galaxy cross term
                arr_deduct(pgal_cross_term, pgal_noisy, pgal, pnoise_1, img_len);
				// galaxy cross term estimate
				arr_add(gal_noisy_a, gal_noisy, noise_2, img_len);
				pow_spec(gal_noisy_a, pgal_noisy_a, size, size);
				arr_deduct(pgal_cross_term_est, pgal_noisy_a, pgal_noisy, pnoise_2, img_len);
				// pure galaxy cross term estimate
				arr_deduct(pgal_pure_cross_term_est, pgal_cross_term_est, pnoise_cross, img_len);



				///////////// Noise free /////////////////////////
				shear_est(pgal, ppsf, &all_paras);
				sub_noise_free_data[row + j * shear_data_cols] = all_paras.n1;
				sub_noise_free_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_noise_free_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_noise_free_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_noise_free_data[row + j * shear_data_cols + 4] = all_paras.dv;

				// ///////////////// Noise free /////////////////////////
				// shear_est(pgal_sqrt, ppsf, &all_paras);
				// sub_noise_free_sqrt_data[row + j * shear_data_cols] = all_paras.n1;
				// sub_noise_free_sqrt_data[row + j * shear_data_cols + 1] = all_paras.n2;
				// sub_noise_free_sqrt_data[row + j * shear_data_cols + 2] = all_paras.dn;
				// sub_noise_free_sqrt_data[row + j * shear_data_cols + 3] = all_paras.du;
				// sub_noise_free_sqrt_data[row + j * shear_data_cols + 4] = all_paras.dv;

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
				
				////////////////// true cross-term image //////////////////////
				shear_est(pgal_cross_term, ppsf, &all_paras);
				sub_cross_term_data[row + j * shear_data_cols] = all_paras.n1;
				sub_cross_term_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_cross_term_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_cross_term_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_cross_term_data[row + j * shear_data_cols + 4] = all_paras.dv;

				// ////////////////// cross-term divided by |P_psf| image //////////////////////
				// shear_est(pgal_cross_term, ppsf_sqrt, &all_paras);
				// sub_cross_term_sqrt_data[row + j * shear_data_cols] = all_paras.n1;
				// sub_cross_term_sqrt_data[row + j * shear_data_cols + 1] = all_paras.n2;
				// sub_cross_term_sqrt_data[row + j * shear_data_cols + 2] = all_paras.dn;
				// sub_cross_term_sqrt_data[row + j * shear_data_cols + 3] = all_paras.du;
				// sub_cross_term_sqrt_data[row + j * shear_data_cols + 4] = all_paras.dv;

				////////////////// galaxy-noise cross-term-est image //////////////////////
				shear_est(pgal_cross_term_est, ppsf, &all_paras);				
				sub_cross_term_estimate[row + j * shear_data_cols] = all_paras.n1;
				sub_cross_term_estimate[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_cross_term_estimate[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_cross_term_estimate[row + j * shear_data_cols + 3] = all_paras.du;
				sub_cross_term_estimate[row + j * shear_data_cols + 4] = all_paras.dv;

               	////////////////// pure galaxy-noise cross-term-est image //////////////////////
				shear_est(pgal_pure_cross_term_est, ppsf, &all_paras);				
				sub_pure_cross_term_estimate[row + j * shear_data_cols] = all_paras.n1;
				sub_pure_cross_term_estimate[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_pure_cross_term_estimate[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_pure_cross_term_estimate[row + j * shear_data_cols + 3] = all_paras.du;
				sub_pure_cross_term_estimate[row + j * shear_data_cols + 4] = all_paras.dv;

               	////////////////// true noise-noise cross-term image //////////////////////
				shear_est(pnoise_cross, ppsf, &all_paras);				
				sub_noise_cross_term[row + j * shear_data_cols] = all_paras.n1;
				sub_noise_cross_term[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_noise_cross_term[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_noise_cross_term[row + j * shear_data_cols + 3] = all_paras.du;
				sub_noise_cross_term[row + j * shear_data_cols + 4] = all_paras.dv;

               	////////////////// noise-noise cross-term-est image //////////////////////
				shear_est(pnoise_cross_est, ppsf, &all_paras);				
				sub_noise_cross_term_estimate[row + j * shear_data_cols] = all_paras.n1;
				sub_noise_cross_term_estimate[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_noise_cross_term_estimate[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_noise_cross_term_estimate[row + j * shear_data_cols + 3] = all_paras.du;
				sub_noise_cross_term_estimate[row + j * shear_data_cols + 4] = all_paras.dv;


#ifdef IMG_CHECK_LABEL
				if(i <= IMG_CHECK_LABEL and shear_id<=IMG_CHECK_LABEL)
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
			if(i <= IMG_CHECK_LABEL and shear_id<=IMG_CHECK_LABEL)
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
		sprintf(set_name, "/data");

		my_Gatherv(sub_noise_free_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
#ifdef EPSF
			sprintf(result_path, "%s/data/data_noise_free_epsf_%d.hdf5", parent_path,shear_id);
#else
			sprintf(result_path, "%s/data/data_noise_free_%d.hdf5", parent_path, shear_id);
#endif	

#ifdef DATA_SEP
			data_sep(total_data, mg_data, total_data_row, shear_data_cols);
			for(k=0; k<shear_data_cols;k++)
			{	
				if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
				else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
			}
#else
			write_h5(result_path, set_name, total_data, total_data_row,shear_data_cols,true);
#endif
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
#ifdef DATA_SEP
			data_sep(total_data, mg_data, total_data_row, shear_data_cols);
			for(k=0; k<shear_data_cols;k++)
			{	
				if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
				else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
			}
#else
			write_h5(result_path, set_name, total_data, total_data_row,shear_data_cols,true);
#endif
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
#ifdef DATA_SEP
			data_sep(total_data, mg_data, total_data_row, shear_data_cols);
			for(k=0; k<shear_data_cols;k++)
			{	
				if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
				else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
			}
#else
			write_h5(result_path, set_name, total_data, total_data_row,shear_data_cols,true);
#endif
		}
		MPI_Barrier(MPI_COMM_WORLD);


		my_Gatherv(sub_cross_term_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
#ifdef EPSF
			sprintf(result_path, "%s/data/data_gal_noise_cross_term_epsf_%d.hdf5", parent_path, shear_id);
#else
			sprintf(result_path, "%s/data/data_gal_noise_cross_term_%d.hdf5", parent_path,shear_id);
#endif

#ifdef DATA_SEP
			data_sep(total_data, mg_data, total_data_row, shear_data_cols);
			for(k=0; k<shear_data_cols;k++)
			{	
				if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
				else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
			}
#else
			write_h5(result_path, set_name, total_data, total_data_row,shear_data_cols,true);
#endif
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_cross_term_estimate, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
#ifdef EPSF
			sprintf(result_path, "%s/data/data_gal_noise_cross_term_est_epsf_%d.hdf5", parent_path, shear_id);
#else
			sprintf(result_path, "%s/data/data_gal_noise_cross_term_est_%d.hdf5", parent_path,shear_id);
#endif

#ifdef DATA_SEP
			data_sep(total_data, mg_data, total_data_row, shear_data_cols);
			for(k=0; k<shear_data_cols;k++)
			{	
				if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
				else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
			}
#else
			write_h5(result_path, set_name, total_data, total_data_row,shear_data_cols,true);
#endif
			std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);


		my_Gatherv(sub_pure_cross_term_estimate, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
#ifdef EPSF
			sprintf(result_path, "%s/data/data_pure_gal_noise_cross_term_est_epsf_%d.hdf5", parent_path, shear_id);
#else
			sprintf(result_path, "%s/data/data_pure_gal_noise_cross_term_est_%d.hdf5", parent_path,shear_id);
#endif

#ifdef DATA_SEP
			data_sep(total_data, mg_data, total_data_row, shear_data_cols);
			for(k=0; k<shear_data_cols;k++)
			{	
				if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
				else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
			}
#else
			write_h5(result_path, set_name, total_data, total_data_row,shear_data_cols,true);
#endif
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_noise_cross_term, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
#ifdef EPSF
			sprintf(result_path, "%s/data/data_noise_noise_cross_term_epsf_%d.hdf5", parent_path, shear_id);
#else
			sprintf(result_path, "%s/data/data_noise_noise_cross_term_%d.hdf5", parent_path,shear_id);
#endif

#ifdef DATA_SEP
			data_sep(total_data, mg_data, total_data_row, shear_data_cols);
			for(k=0; k<shear_data_cols;k++)
			{	
				if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
				else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
			}
#else
			write_h5(result_path, set_name, total_data, total_data_row,shear_data_cols,true);
#endif
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_noise_cross_term_estimate, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
#ifdef EPSF
			sprintf(result_path, "%s/data/data_noise_noise_cross_term_est_epsf_%d.hdf5", parent_path, shear_id);
#else
			sprintf(result_path, "%s/data/data_noise_noise_cross_term_est_%d.hdf5", parent_path,shear_id);
#endif

#ifdef DATA_SEP
			data_sep(total_data, mg_data, total_data_row, shear_data_cols);
			for(k=0; k<shear_data_cols;k++)
			{	
				if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
				else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
			}
#else
			write_h5(result_path, set_name, total_data, total_data_row,shear_data_cols,true);
#endif
		}
		MPI_Barrier(MPI_COMM_WORLD);

// 		my_Gatherv(sub_noise_free_sqrt_data, gather_count, total_data, numprocs, rank);
// 		if (0 == rank)
// 		{
// #ifdef EPSF
// 			sprintf(result_path, "%s/data/data_noise_free_sqrt_epsf_%d.hdf5", parent_path,shear_id);
// #else
// 			sprintf(result_path, "%s/data/data_noise_free_sqrt_%d.hdf5", parent_path, shear_id);
// #endif
// 			data_sep(total_data, mg_data, total_data_row, shear_data_cols);
// 			for(k=0; k<shear_data_cols;k++)
// 			{	
// 				if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
// 				else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
// 			}
// 		}
// 		MPI_Barrier(MPI_COMM_WORLD);


// 		my_Gatherv(sub_cross_term_sqrt_data, gather_count, total_data, numprocs, rank);
// 		if (0 == rank)
// 		{
// #ifdef EPSF
// 			sprintf(result_path, "%s/data/data_cross_term_sqrt_epsf_%d.hdf5", parent_path, shear_id);
// #else
// 			sprintf(result_path, "%s/data/data_cross_term_sqrt_%d.hdf5", parent_path,shear_id);
// #endif
// 			data_sep(total_data, mg_data, total_data_row, shear_data_cols);
// 			for(k=0; k<shear_data_cols;k++)
// 			{	
// 				if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
// 				else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
// 			}
// 		}
// 		MPI_Barrier(MPI_COMM_WORLD);


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

#ifdef FLUX_PDF
	delete[] sub_flux;
#endif
	delete[] scatter_count;
	delete[] gather_count;
	
	MPI_Finalize();
	return 0;
}
