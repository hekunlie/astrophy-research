#include<FQlib.h>
#include<hk_mpi.h>
#include<hk_iolib.h>
#define IMG_CHECK_LABEL 3
#define FLUX_PDF_UNI
#define EPSF

#define SAVE_MEM

#define DATA_SEP

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

void arr_scale(MY_FLOAT *arr1, const MY_FLOAT scale,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] += arr1[i]*scale;
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

void write_data(char *result_path, MY_FLOAT *total_data, char **mg_name, MY_FLOAT **mg_data, int total_data_row, int shear_data_cols)
{	
	int k;
	data_sep(total_data, mg_data, total_data_row, shear_data_cols);
	for(k=0; k<shear_data_cols;k++)
	{	
		if(0 == k){write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, true);}
		else{write_h5(result_path, mg_name[k], mg_data[k], total_data_row, 1, false);}
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

	char parent_path[500], chip_path[500], para_path[500], shear_path[500], result_path[500], log_path[500];
	char buffer[500], log_inform[500], set_name[50];

	
	
	int i, j, k, ib, m,n;
	int sss1, sss2, seed_pts, seed_n1, seed_n2, seed_step;

	int num_p, size, shear_pairs, img_len;
    MY_FLOAT max_radius;
	int total_chips, sub_chip_num, sub_data_row, total_data_row;
	int stamp_num, stamp_nx, shear_data_cols;		
	int row, chip_st, chip_ed, shear_id, psf_type;
	int seed_ini;

	int detect_label;
	std::string detect_info;
	double temp_val;

	MY_FLOAT psf_scale, psf_thresh_scale, sig_level, psf_noise_sig, gal_noise_sig, flux_i, mag_i;
	MY_FLOAT g1, g2, ts, te, t1, t2;
	MY_FLOAT *g1t, *g2t;
	MY_FLOAT psf_ellip, ellip_theta;
	MY_FLOAT img_cent;
	MY_FLOAT pts_step;
	MY_FLOAT gal_fluxs[7]{4000, 6000, 9000, 13500, 20250, 30375, 45562},flux_scale;
	// MY_FLOAT gal_noise_sigma[5]{48,54,60,66,72};
	MY_FLOAT gal_noise_sigma[5]{60,120,240,150,240};
	MY_FLOAT pts_steps[4]{0.5, 1, 1.5, 2};

	int flux_tag;
	int pts_tag;
	int flux_num;
	int kernel_size1, kernel_size2;

	strcpy(parent_path, argv[1]);
	seed_ini = atoi(argv[2]);
	flux_tag = atoi(argv[3]);
	total_chips = atoi(argv[4]);
	size = atoi(argv[5]);
	shear_pairs = atoi(argv[6]);
	// pts_tag = atof(argv[3]);
	pts_step = 1;//pts_steps[pts_tag];


	num_p = 30;
	max_radius= 6;
	stamp_num = 10000;
	shear_data_cols = 5;

	psf_type = 2;
	psf_scale = 4;//psf_scales[psf_scale_id];
	psf_ellip = 0.10;
	ellip_theta = 0;
	psf_thresh_scale = 2.;

	sig_level = 1.5;
	psf_noise_sig = 0;
    gal_noise_sig = 60;
	
    img_len = size*size;
	img_cent = size*0.5 - 0.5;
    stamp_nx = 100;

	
	flux_scale = gal_fluxs[flux_tag]/gal_fluxs[0];
	flux_i = gal_fluxs[flux_tag]/num_p;
	
	flux_num = 1000000;

	kernel_size1 = 3;
	kernel_size2 = 5;

	all_paras.stamp_size = size;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.gal_noise_sig = gal_noise_sig;
	
	sprintf(log_path, "%s/logs/%002d.dat", parent_path, rank);

	///////////////////// task distribution /////////////////////////////////////
	int *scatter_count,*gather_count;
	// for scatterring the flux to each thread
	scatter_count = new int[numprocs]{};
	// for the gatherv when finish in each shear point
	gather_count = new int[numprocs]{};

	task_alloc(total_chips, numprocs, rank, chip_st, chip_ed, scatter_count);
	sub_chip_num = scatter_count[rank];
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


#ifdef IMG_CHECK_LABEL
	MY_FLOAT *big_img_check[10];
	for(i=0;i<1;i++)
	{
		big_img_check[i] = new MY_FLOAT[stamp_nx*stamp_nx*img_len]{};
	}
	
#endif

	MY_FLOAT *point = new MY_FLOAT[2 * num_p]{};
	MY_FLOAT *point_r = new MY_FLOAT[2 * num_p]{};

	MY_FLOAT *stamp_img[20];
	MY_FLOAT *stamp_pow_img[20];
	MY_FLOAT *psf_img[10];
	MY_FLOAT *noise_img[10];
	MY_FLOAT *noise_pow_img[10];

	
	// the shear estimators data matrix  
	MY_FLOAT *total_data;	

	MY_FLOAT *sub_noise_free_data;
    MY_FLOAT *sub_noisy_data;
    MY_FLOAT *sub_noise_residual_data;
	MY_FLOAT *sub_cross_term_data;
	MY_FLOAT *sub_cross_term_sqrt_data;


	MY_FLOAT *sub_cross_term_est_data[5];
	MY_FLOAT *sub_cross_term_est_data_r[5];

	MY_FLOAT *galaxy_flux;
	MY_FLOAT *mg_data[5];
	char *mg_name[5];

	MY_FLOAT *total_flux, *sub_flux;
	
	for(i=0;i<10;i++)
	{
		stamp_img[i] = new MY_FLOAT[img_len];
		stamp_pow_img[i] = new MY_FLOAT[img_len];
		noise_img[i] = new MY_FLOAT[img_len];
		noise_pow_img[i] = new MY_FLOAT[img_len];
		if(i<6){psf_img[i] = new MY_FLOAT[img_len]{};}
	}

	sub_noise_free_data = new MY_FLOAT[gather_count[rank]]{};
	sub_noisy_data = new MY_FLOAT[gather_count[rank]]{};
	sub_cross_term_data = new MY_FLOAT[gather_count[rank]]{};
	sub_cross_term_sqrt_data = new MY_FLOAT[gather_count[rank]]{};
	sub_noise_residual_data = new MY_FLOAT[gather_count[rank]]{};

	for(i=0;i<3;i++)
	{
		sub_cross_term_est_data[i] = new MY_FLOAT[gather_count[rank]]{};
		sub_cross_term_est_data_r[i] = new MY_FLOAT[gather_count[rank]]{};
	}

	// galaxy_flux = new MY_FLOAT[flux_num]{};

	// seed distribution, different thread gets different seed
	seed_step = 1;
	sss1 = 2*seed_step*shear_pairs;
	seed_pts = sss1*rank + 1 + seed_ini;//35000;
	seed_n1 = sss1*rank + 1 + seed_ini;// 4001*(rotation+1);
	seed_n2 = sss1*rank + 1 + seed_ini;//2300*(rotation+1);

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

	sprintf(shear_path,"%s/parameters/shear.hdf5", parent_path);
	sprintf(set_name,"/g1");
	read_h5(shear_path, set_name, g1t);
	sprintf(set_name,"/g2");
	read_h5(shear_path, set_name, g2t);


#ifdef EPSF
	create_psf_e(psf_img[0], psf_scale, size, img_cent, psf_ellip, ellip_theta, psf_type);
#else
	// create PSF
	create_psf(psf_img[0], psf_scale, size, img_cent, psf_type);
#endif
	pow_spec(psf_img[0], psf_img[1], size, size);
	arr_sqrt(psf_img[1], psf_img[4], img_len);
	get_psf_radius(psf_img[1], &all_paras, psf_thresh_scale);


	image_rotation(psf_img[0], psf_img[2], size);
	pow_spec(psf_img[2], psf_img[3], size, size);

	
	if (0 == rank)
	{	
		std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		std::cout << parent_path << std::endl;
		std::cout<<"Shear num: "<<shear_pairs<<std::endl;
		std::cout << "Total chip: " << total_chips<< ", Stamp size: " << size  << std::endl;
		std::cout << "Total cpus: " << numprocs << std::endl;
		std::cout <<"PSF Scale: "<<psf_scale<< " PSF THRESH: " << all_paras.psf_pow_thresh <<" PSF HLR: " << all_paras.psf_hlr << std::endl;
		std::cout <<"MAX RADIUS: "<< max_radius <<" , Step: "<<pts_step<< ", SIG_LEVEL: " << sig_level <<"sigma"<< std::endl;
#ifdef EPSF
		sprintf(buffer, "!%s/imgs/epsf_%.2f.fits", parent_path,psf_scale);
		sprintf(chip_path, "!%s/imgs/epsf_r_%.2f.fits", parent_path,psf_scale);

#else
		sprintf(buffer, "!%s/psf_%.2f.fits", parent_path,psf_scale);
		sprintf(chip_path, "!%s/psf_r_%.2f.fits", parent_path,psf_scale);
#endif
		write_fits(buffer,psf_img[0], size, size);
		write_fits(chip_path,psf_img[2], size, size);

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

		// // read flux
		// sprintf(shear_path,"%s/parameters/para_%d.hdf5", parent_path, shear_id);
		// sprintf(set_name,"/flux");
		// read_h5(shear_path, set_name, galaxy_flux);
		
		// rank 0 reads the total flux array and scatters to each thread
		if(rank == 0)
		{
			std::cout<<"---------------------------------------------------------------------------"<<std::endl;
			std::cout<<"g1: "<<g1<<" g2: "<<g2<<" FLux:"<<flux_i<<" Flux Scale:"<<flux_scale<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);


		// initialize GSL
		gsl_initialize(seed_pts,0);
		gsl_initialize(seed_n1,1);
		// gsl_initialize(seed_n2,2);
		seed_pts += seed_step;
		seed_n1 += seed_step;
		seed_n2 += seed_step;
		// loop the chips
		for (i = chip_st; i < chip_ed; i++)
		{
			t1 = clock();



			sprintf(log_inform, "RANK: %03d, SHEAR %02d:, chip: %05d, start. seed: %d, %d", rank,shear_id, i, seed_pts, seed_n1);
			write_log(log_path, log_inform);
			if (rank == 0)
			{
				std::cout << log_inform << std::endl;
			}

			row = (i-chip_st)*stamp_num*shear_data_cols;

			sprintf(chip_path, "%s/imgs/%d/gal_chip_%05d_noise_free.fits", parent_path, shear_id, i);
			read_fits(chip_path, big_img_check[0]);

			sprintf(log_inform, "RANK: %03d, SHEAR %02d:, chip: %05d read", rank,shear_id, i);
			write_log(log_path, log_inform);
			if (rank == 0)
			{
				std::cout << log_inform << std::endl;
			}
			
			// loop the stamps
			for (j = 0; j < stamp_num; j++)
			{	
				if(j == 0 or j==stamp_num/5 or j == stamp_num/5*2 or j == stamp_num/5*3 or j == stamp_num/5*4)
				{
					sprintf(log_inform, "RANK: %03d, SHEAR %02d:, chip: %05d stamp: %d", rank,shear_id, i, j);
					write_log(log_path, log_inform);
				}

				// n = (i*stamp_num + j)%flux_num;
				// flux_i = galaxy_flux[n]/num_p;

				initialize_arr(point, num_p * 2, 0);				
				for(k=0;k<8;k++)
				{
					initialize_arr(stamp_img[k], size*size, 0);
					initialize_arr(noise_img[k], size*size, 0);
				}
				for(k=0;k<3;k++)
				{					
					addnoise(noise_img[k], img_len, gal_noise_sig, rng1);
					pow_spec(noise_img[k], noise_pow_img[k], size, size);
				}

				segment(big_img_check[0], stamp_img[0], j, size, stamp_nx, stamp_nx);
				arr_scale(stamp_img[0], flux_scale, img_len);

// 				create_points(point, num_p, max_radius, pts_step, rng0);

// #ifdef EPSF
// 				convolve_e(point,num_p,flux_i, g1, g2, stamp_img[0], size, img_cent, psf_scale,psf_type,psf_ellip, ellip_theta);
// #else
// 				convolve(point,num_p,flux_i, g1, g2, stamp_img[0], size, img_cent, psf_scale, psf_type);

// #endif
				//if(rank == 0 and shear_id==0 and i < 3)
				// {
					// stack(big_img_check[0], stamp_img[0], j, size, stamp_nx, stamp_nx);
				// }
				// noise free
				pow_spec(stamp_img[0], stamp_pow_img[0], size, size);
				// noisy image
				arr_add(stamp_img[1], stamp_img[0], noise_img[0], img_len);
				pow_spec(stamp_img[1], stamp_pow_img[1], size, size);
				arr_deduct(stamp_pow_img[2], stamp_pow_img[1], noise_pow_img[1], img_len);
				// true galaxy-noise cross term
                arr_deduct(stamp_pow_img[3], stamp_pow_img[1], stamp_pow_img[0], noise_pow_img[0], img_len);
				// noise power residual
				arr_deduct(noise_pow_img[3], noise_pow_img[0], noise_pow_img[1], img_len);

				// estimated galaxy-noise cross term (include noise-noise cross term)
				m = 2;
				arr_add(stamp_img[4], stamp_img[1], noise_img[m], img_len);
				pow_spec(stamp_img[4], stamp_pow_img[4], size, size);
				arr_deduct(stamp_pow_img[5], stamp_pow_img[4], stamp_pow_img[1], noise_pow_img[m], img_len);
				// noise-noise cross term in estimated galaxy-noise cross term  
				arr_add(noise_img[4], noise_img[0], noise_img[m], img_len);
				pow_spec(noise_img[4], noise_pow_img[4], size, size);
				arr_deduct(noise_pow_img[5], noise_pow_img[4], noise_pow_img[0], noise_pow_img[m], img_len);
				// pure galaxy noise cross term
				arr_deduct(stamp_pow_img[6], stamp_pow_img[5], noise_pow_img[5], img_len);


				/////////////////////// Noise free /////////////////////////
				shear_est(stamp_pow_img[0], psf_img[1], &all_paras);
				sub_noise_free_data[row + j * shear_data_cols] = all_paras.n1;
				sub_noise_free_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_noise_free_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_noise_free_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_noise_free_data[row + j * shear_data_cols + 4] = all_paras.dv;

				/////////////////////////// noisy image /////////////////////////////////
 				shear_est(stamp_pow_img[2], psf_img[1], &all_paras);
				sub_noisy_data[row + j * shear_data_cols] = all_paras.n1;
				sub_noisy_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_noisy_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_noisy_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_noisy_data[row + j * shear_data_cols + 4] = all_paras.dv;
				
				////////////////// true cross-term image //////////////////////
				shear_est(stamp_pow_img[3], psf_img[1], &all_paras);
				sub_cross_term_data[row + j * shear_data_cols] = all_paras.n1;
				sub_cross_term_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_cross_term_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_cross_term_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_cross_term_data[row + j * shear_data_cols + 4] = all_paras.dv;

				shear_est(stamp_pow_img[3], psf_img[4], &all_paras);
				sub_cross_term_sqrt_data[row + j * shear_data_cols] = all_paras.n1;
				sub_cross_term_sqrt_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_cross_term_sqrt_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_cross_term_sqrt_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_cross_term_sqrt_data[row + j * shear_data_cols + 4] = all_paras.dv;

				////////////////// noise power residual image //////////////////////
				shear_est(noise_pow_img[3], psf_img[1], &all_paras);
				sub_noise_residual_data[row + j * shear_data_cols] = all_paras.n1;
				sub_noise_residual_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_noise_residual_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_noise_residual_data[row + j * shear_data_cols + 3] = all_paras.du;
				sub_noise_residual_data[row + j * shear_data_cols + 4] = all_paras.dv;




				////////////////// galaxy-noise cross-term-est image ////////////////
				shear_est(stamp_pow_img[5], psf_img[1], &all_paras);				
				sub_cross_term_est_data[0][row + j * shear_data_cols] = all_paras.n1;
				sub_cross_term_est_data[0][row + j * shear_data_cols + 1] = all_paras.n2;
				sub_cross_term_est_data[0][row + j * shear_data_cols + 2] = all_paras.dn;
				sub_cross_term_est_data[0][row + j * shear_data_cols + 3] = all_paras.du;
				sub_cross_term_est_data[0][row + j * shear_data_cols + 4] = all_paras.dv;

				////////////////// pure galaxy-noise cross-term-est image ////////////////
				shear_est(stamp_pow_img[6], psf_img[1], &all_paras);				
				sub_cross_term_est_data[1][row + j * shear_data_cols] = all_paras.n1;
				sub_cross_term_est_data[1][row + j * shear_data_cols + 1] = all_paras.n2;
				sub_cross_term_est_data[1][row + j * shear_data_cols + 2] = all_paras.dn;
				sub_cross_term_est_data[1][row + j * shear_data_cols + 3] = all_paras.du;
				sub_cross_term_est_data[1][row + j * shear_data_cols + 4] = all_paras.dv;

				////////////////// noise-noise cross-term in galaxy-noise cross-term est image ///////////
				shear_est(noise_pow_img[5], psf_img[1], &all_paras);				
				sub_cross_term_est_data[2][row + j * shear_data_cols] = all_paras.n1;
				sub_cross_term_est_data[2][row + j * shear_data_cols + 1] = all_paras.n2;
				sub_cross_term_est_data[2][row + j * shear_data_cols + 2] = all_paras.dn;
				sub_cross_term_est_data[2][row + j * shear_data_cols + 3] = all_paras.du;
				sub_cross_term_est_data[2][row + j * shear_data_cols + 4] = all_paras.dv;

				// shear_est(stamp_pow_img[6], psf_img[4], &all_paras);				
				// sub_cross_term_est_data[2][row + j * shear_data_cols] = all_paras.n1;
				// sub_cross_term_est_data[2][row + j * shear_data_cols + 1] = all_paras.n2;
				// sub_cross_term_est_data[2][row + j * shear_data_cols + 2] = all_paras.dn;
				// sub_cross_term_est_data[2][row + j * shear_data_cols + 3] = all_paras.du;
				// sub_cross_term_est_data[2][row + j * shear_data_cols + 4] = all_paras.dv;



				//////////////// galaxy-noise cross-term-est image psf_rotation ////////////////
				shear_est(stamp_pow_img[5], psf_img[3], &all_paras);				
				sub_cross_term_est_data_r[0][row + j * shear_data_cols] = all_paras.n1;
				sub_cross_term_est_data_r[0][row + j * shear_data_cols + 1] = all_paras.n2;
				sub_cross_term_est_data_r[0][row + j * shear_data_cols + 2] = all_paras.dn;
				sub_cross_term_est_data_r[0][row + j * shear_data_cols + 3] = all_paras.du;
				sub_cross_term_est_data_r[0][row + j * shear_data_cols + 4] = all_paras.dv;

				// ////////////////// pure galaxy-noise cross-term-est image ////////////////
				// shear_est(stamp_pow_img[6], psf_img[3], &all_paras);				
				// sub_cross_term_est_data_r[1][row + j * shear_data_cols] = all_paras.n1;
				// sub_cross_term_est_data_r[1][row + j * shear_data_cols + 1] = all_paras.n2;
				// sub_cross_term_est_data_r[1][row + j * shear_data_cols + 2] = all_paras.dn;
				// sub_cross_term_est_data_r[1][row + j * shear_data_cols + 3] = all_paras.du;
				// sub_cross_term_est_data_r[1][row + j * shear_data_cols + 4] = all_paras.dv;

				// ////////////////// noise-noise cross-term in galaxy-noise cross-term est image ////////////////
				// shear_est(noise_pow_img[5], psf_img[3], &all_paras);				
				// sub_cross_term_est_data_r[2][row + j * shear_data_cols] = all_paras.n1;
				// sub_cross_term_est_data_r[2][row + j * shear_data_cols + 1] = all_paras.n2;
				// sub_cross_term_est_data_r[2][row + j * shear_data_cols + 2] = all_paras.dn;
				// sub_cross_term_est_data_r[2][row + j * shear_data_cols + 3] = all_paras.du;
				// sub_cross_term_est_data_r[2][row + j * shear_data_cols + 4] = all_paras.dv;

				// ////////////////// noise-noise cross-term in galaxy-noise cross-term est image ////////////////
				// shear_est(stamp_pow_img[6], psf_img[5], &all_paras);				
				// sub_cross_term_est_data_r[2][row + j * shear_data_cols] = all_paras.n1;
				// sub_cross_term_est_data_r[2][row + j * shear_data_cols + 1] = all_paras.n2;
				// sub_cross_term_est_data_r[2][row + j * shear_data_cols + 2] = all_paras.dn;
				// sub_cross_term_est_data_r[2][row + j * shear_data_cols + 3] = all_paras.du;
				// sub_cross_term_est_data_r[2][row + j * shear_data_cols + 4] = all_paras.dv;

			}

			// if(rank == 0 and shear_id==0 and i < 3)
			// {
				// sprintf(chip_path, "!%s/imgs/%d/gal_chip_%05d_noise_free.fits", parent_path, shear_id, i);
			 	// write_fits(chip_path, big_img_check[0], stamp_nx*size, stamp_nx*size);
			// }
			t2 = clock();
			sprintf(log_inform, "RANK: %03d, SHEAR %02d: chip: %05d, done in %.2f s.", rank, shear_id, i, (t2 - t1) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			if (rank == 0)
			{
				std::cout << log_inform << std::endl;
			}
		}
		gsl_free(0);
		gsl_free(1);
		// gsl_free(2);


		// finish the chip loop
		MPI_Barrier(MPI_COMM_WORLD);
		
		my_Gatherv(sub_noise_free_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data_%d/data_noise_free_%d.hdf5", parent_path,flux_tag,shear_id);
			write_data(result_path, total_data, mg_name, mg_data, total_data_row, shear_data_cols);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_noisy_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data_%d/data_noisy_cpp_%d.hdf5", parent_path, flux_tag,shear_id);
			write_data(result_path, total_data, mg_name, mg_data, total_data_row, shear_data_cols);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_cross_term_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data_%d/data_gal_noise_cross_term_%d.hdf5", parent_path,flux_tag,  shear_id);
			write_data(result_path, total_data, mg_name, mg_data, total_data_row, shear_data_cols);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_cross_term_sqrt_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data_%d/data_gal_noise_cross_term_sqrt_%d.hdf5", parent_path,flux_tag,  shear_id);
			write_data(result_path, total_data, mg_name, mg_data, total_data_row, shear_data_cols);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_noise_residual_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data_%d/data_noise_residual_%d.hdf5", parent_path,flux_tag,  shear_id);
			write_data(result_path, total_data, mg_name, mg_data, total_data_row, shear_data_cols);
		}
		MPI_Barrier(MPI_COMM_WORLD);



		my_Gatherv(sub_cross_term_est_data[0], gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data_%d/data_gal_noise_cross_term_est_%d.hdf5", parent_path,flux_tag,shear_id);
			write_data(result_path, total_data, mg_name, mg_data, total_data_row, shear_data_cols);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_cross_term_est_data[1], gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data_%d/data_pure_gal_noise_cross_term_est_%d.hdf5", parent_path,flux_tag, shear_id);
			write_data(result_path, total_data, mg_name, mg_data, total_data_row, shear_data_cols);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_cross_term_est_data[2], gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data_%d/data_noise_noise_cross_term_est_%d.hdf5", parent_path,flux_tag, shear_id);
			write_data(result_path, total_data, mg_name, mg_data, total_data_row, shear_data_cols);
		}
		MPI_Barrier(MPI_COMM_WORLD);



		my_Gatherv(sub_cross_term_est_data_r[0], gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data_%d/data_gal_noise_cross_term_est_r_%d.hdf5", parent_path, flux_tag, shear_id);
			write_data(result_path, total_data, mg_name, mg_data, total_data_row, shear_data_cols);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// my_Gatherv(sub_cross_term_est_data_r[1], gather_count, total_data, numprocs, rank);
		// if (0 == rank)
		// {
		// 	sprintf(result_path, "%s/data/data_pure_gal_noise_cross_term_est_r_%d.hdf5", parent_path,shear_id);
		// 	write_data(result_path, total_data, mg_name, mg_data, total_data_row, shear_data_cols);
		// }
		// MPI_Barrier(MPI_COMM_WORLD);

		// my_Gatherv(sub_cross_term_est_data_r[2], gather_count, total_data, numprocs, rank);
		// if (0 == rank)
		// {
		// 	sprintf(result_path, "%s/data/data_pure_gal_noise_cross_term_est_r_sqrt_%d.hdf5", parent_path,shear_id);
		// 	write_data(result_path, total_data, mg_name, mg_data, total_data_row, shear_data_cols);
		// }
		// MPI_Barrier(MPI_COMM_WORLD);


		te = clock();
		if (rank == 0)
		{
			sprintf(buffer, "rank %d: done in %.2f sec %s/imgs/%d\n", rank, (te - ts) / CLOCKS_PER_SEC, parent_path, shear_id);
			std::cout << buffer<<std::endl;
			std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		}
	}

	if (0 == rank)
	{	
		std::cout<<parent_path<<std::endl;
		std::cout << "FINISH ALL JOBS" << std::endl;
	}
	
	MPI_Finalize();
	return 0;
}
