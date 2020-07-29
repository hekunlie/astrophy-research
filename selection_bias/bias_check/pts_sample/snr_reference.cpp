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

	
    double t1, t2, te, ts;
	int i, j, k, ib, m,n;
	int sss1, sss2, seed_pts, seed_n1, seed_n2, seed_step;

	int num_p, size, img_len;
    MY_FLOAT max_radius;
	int total_chips, data_row, data_col;
	int stamp_num, stamp_nx;
	int seed_ini;

	int detect_label;
	std::string detect_info;
	double temp_val;
 

	MY_FLOAT psf_scale, psf_type, psf_thresh_scale, sig_level, psf_noise_sig, gal_noise_sig, flux_i;
	MY_FLOAT psf_ellip, ellip_theta;
	MY_FLOAT img_cent;
	MY_FLOAT pts_step;
	MY_FLOAT gal_flux[8], noise_sigma[5];
    MY_FLOAT flux_st;
    int ix, iy, tag;


	strcpy(parent_path, argv[1]);
	seed_ini = atoi(argv[2]);
	total_chips = atoi(argv[3]);
	size = atoi(argv[4]);
    flux_st = atof(argv[5]);

    for(i=0;i<8;i++){ gal_flux[i] = (i+flux_st)*1000; }
    for(i=0;i<5;i++){ noise_sigma[i] = i*10+30; }


    iy = rank/8;
    ix = rank%8;

	pts_step = 1;//atof(argv[1]);

	num_p = 30;
	max_radius = 7;
	stamp_num = 10000;
    img_len = size*size;
	img_cent = size*0.5 - 0.5;
    stamp_nx = 100;

    data_row = total_chips*stamp_num;
	data_col = 4;


	sig_level = 1.5;
	psf_noise_sig = 0;

    gal_noise_sig = noise_sigma[iy];
    flux_i = gal_flux[ix]/num_p;



	psf_type = 2;
	psf_scale = 4;//psf_scales[psf_scale_id];
	psf_ellip = 0.10;

	ellip_theta = 0;
	psf_thresh_scale = 2.;

	all_paras.gal_noise_sig = gal_noise_sig;
	all_paras.psf_noise_sig = 0;
	all_paras.stamp_size = size;
	all_paras.max_source = 30;
	all_paras.area_thresh = 5;
	all_paras.detect_thresh = gal_noise_sig*sig_level;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.max_distance = 6; 


	MY_FLOAT *big_img_check[10];
	int *big_mask = new int[stamp_num*img_len];

	for(i=0;i<2;i++)
	{
		big_img_check[i] = new MY_FLOAT[stamp_num*img_len]{};
	}

    MY_FLOAT *snr_data = new MY_FLOAT[data_col*data_row];

	MY_FLOAT *point = new MY_FLOAT[2 * num_p]{};

    int *mask = new int[img_len]{};
	MY_FLOAT *gal = new MY_FLOAT[img_len]{};
	MY_FLOAT *gal_noisy = new MY_FLOAT[img_len]{};
	MY_FLOAT *noise_1 = new MY_FLOAT[img_len]{};


	// seed distribution, different thread gets different seed
	seed_step = 1;
	sss1 = 2*seed_step*total_chips;
	seed_pts = sss1*rank + 1 + seed_ini;//35000;
	seed_n1 = sss1*rank + 1 + seed_ini*2;// 4001*(rotation+1);

	for(i=0;i<numprocs;i++)
    {
        if(i == rank){std::cout<<rank<<" "<<iy<<" "<<ix<<" "<<noise_sigma[iy]<<" "<<gal_flux[ix]<<" "<<std::endl;}
        MPI_Barrier(MPI_COMM_WORLD);
    }

    ts = clock();

    // loop the chips
    for (i = 0; i < total_chips; i++)
    {
        t1 = clock();

        // initialize GSL
        gsl_initialize(seed_pts,0);
        gsl_initialize(seed_n1,1);
        seed_pts += seed_step;
        seed_n1 += seed_step;

        k = i*stamp_num*data_col;
        // loop the stamps
        for (j = 0; j < stamp_num; j++)
        {	
            initialize_arr(point, num_p * 2, 0);

            initialize_arr(noise_1, img_len, 0);
            addnoise(noise_1, img_len, gal_noise_sig, rng1);

            create_points(point, num_p, max_radius, pts_step, rng0);

#ifdef EPSF
            convolve_e(point,num_p, flux_i, 0, 0, gal, size, img_cent, psf_scale, psf_type, psf_ellip, ellip_theta);
#else
            convolve(point, num_p, flux_i, 0, 0, gal, size, img_cent, psf_scale,psf_type);

#endif
            arr_add(gal_noisy, gal, noise_1, img_len);

            galaxy_finder(gal_noisy, mask, &all_paras, false, detect_label, detect_info);

            snr_data[k + j * data_col] = detect_label;
            snr_data[k + j * data_col + 1] = all_paras.gal_osnr;
            snr_data[k + j * data_col + 2] = all_paras.gal_size;
            snr_data[k + j * data_col + 3] = all_paras.gal_hflux/sqrt(all_paras.gal_hsize)/gal_noise_sig;

            stack(big_img_check[0], gal, j, size, stamp_nx, stamp_nx);
            stack(big_img_check[1], gal_noisy, j, size, stamp_nx, stamp_nx);
            stack(big_mask, mask, j, size, stamp_nx, stamp_nx);
            // if(j < 10)
            // {std::cout<<detect_label<<" "<<all_paras.gal_osnr<<" "<<all_paras.gal_size<<std::endl;}
        }
        gsl_free(0);
        gsl_free(1);

        // show_arr(snr_data,10,3);
        
        sprintf(chip_path, "!%s/img/gal_chip_%d_%d_noise_free.fits", parent_path, iy, ix);
        write_fits(chip_path, big_img_check[0], stamp_nx*size, stamp_nx*size);
        sprintf(chip_path, "!%s/img/gal_chip_%d_%d_noisy.fits", parent_path, iy, ix);
        write_fits(chip_path, big_img_check[1], stamp_nx*size, stamp_nx*size);

        sprintf(chip_path, "!%s/img/gal_chip_%d_%d_mask.fits", parent_path, iy, ix);
        write_fits(chip_path, big_mask, stamp_nx*size, stamp_nx*size);      

        t2 = clock();
        sprintf(log_inform, "RANK: %03d, chip: %d, done in %.2f s.", rank, i, (t2 - t1) / CLOCKS_PER_SEC);
        // write_log(log_path, log_inform);
        if (rank == 0)
        {
            std::cout << log_inform << std::endl;
        }	

    }

    sprintf(result_path,"%s/data/snr_%d_%d.hdf5", parent_path, iy, ix);
    sprintf(set_name,"/data");
    write_h5(result_path, set_name, snr_data, data_row, data_col, true);
    MPI_Barrier(MPI_COMM_WORLD);

    te = clock();


	if (0 == rank)
	{	
		std::cout<<parent_path<<std::endl;
		std::cout << "FINISH ALL JOBS" << std::endl;
	}
	
	MPI_Finalize();
	return 0;
}
