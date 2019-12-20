#include<FQlib.h>
#include<hk_mpi.h>
#include<hk_iolib.h>
#define IMG_CHECK_LABEL 2
#define TIME_CAL

void arr_add(double *arr1, const double*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] += arr2[i];
    }
}
int main(int argc, char*argv[])
{	
	fq_paras all_paras;

	char parent_path[100], chip_path[150], para_path[150], shear_path[150], result_path[150], log_path[150];
    char time_now[50];

	int i, j, k, ib;
	int sss1, sss2, seed_pts, seed_n1, seed_n2, seed_step;

	int num_p, size, shear_pairs;
    double max_radius;
	int total_chips, sub_chip_num, sub_data_row, total_data_row;
	int stamp_num, stamp_nx, shear_data_cols;		
	int row, chip_st, chip_ed, shear_id, psf_type, temp_s, detect_label;
	int img_len;

	double psf_scale, psf_thresh_scale, sig_level, psf_noise_sig, gal_noise_sig;
	double g1, g2, ts, te, t1, t2;
    double flux_i;


	num_p = 50;
	max_radius= 8;
	stamp_num = 3000;
	shear_data_cols = 4;

    shear_pairs = 10;

	psf_type=2;
	psf_scale = 4;
	psf_thresh_scale = 2.;
	sig_level = 1.5;
	psf_noise_sig = 0;
    gal_noise_sig = 60;

	size = 64;
    img_len = size*size;
	total_chips = 1000;
    stamp_nx = 100;

	all_paras.stamp_size = size;
	all_paras.img_x = size;
	all_paras.img_y = size;
	// because the max half light radius of the galsim source is 5.5 pixels
	all_paras.max_distance = max_radius; 

	double *point = new double[2 * num_p]{};

	double *gal = new double[size*size]{};
	double *pgal = new double[size*size]{};

	double *psf = new double[size*size]{};
	double *ppsf = new double[size*size]{};

	double *noise_1 = new double[size*size]{};
	double *pnoise_1 = new double[size*size]{};

    double *noise_2 = new double[size*size]{};
	double *pnoise_2 = new double[size*size]{};

	// the shear estimators data matrix  
	double *total_data = new double[shear_data_cols*stamp_num]{};
    double *time_label = new double[10]{};
    double *time_collection = new double[10]{};
	// create PSF
	create_psf(psf, psf_scale, size, psf_type);
	pow_spec(psf, ppsf, size, size);
	get_psf_radius(ppsf, &all_paras, psf_thresh_scale);

    t1 = clock();

    // initialize GSL
    gsl_initialize(1,0);
    gsl_initialize(2,1);
    gsl_initialize(3,2);
    seed_pts += seed_step;
    seed_n1 += seed_step;
    seed_n2 += seed_step;

    get_time(time_now, 50);
    std::cout<<"Start..."<<time_now<<std::endl;
    // loop the stamps
    for (j = 0; j < stamp_num; j++)
    {	
        
#ifdef TIME_CAL
        time_label[0] = clock();
#endif
        flux_i = 9000./ num_p;
        
        initialize_arr(gal, img_len, 0);
        initialize_arr(pgal, img_len, 0);
        initialize_arr(point, num_p * 2, 0);
        initialize_arr(noise_1, img_len, 0);
        initialize_arr(pnoise_1, img_len, 0);

#ifdef TIME_CAL
        time_label[1] = clock();
#endif
        create_points(point, num_p, max_radius, rng0);

#ifdef TIME_CAL
        time_label[2] = clock();
#endif
        convolve(gal, point, flux_i, size, num_p, 0, psf_scale, g1, g2, psf_type, 1, &all_paras);	

#ifdef TIME_CAL
        time_label[3] = clock();
#endif
        addnoise(gal, img_len, gal_noise_sig, rng1);
        addnoise(noise_1, img_len, gal_noise_sig, rng1);

#ifdef TIME_CAL
        time_label[4] = clock();
#endif
        pow_spec(noise_1, pnoise_1, size, size);
        pow_spec(gal, pgal, size, size);

#ifdef TIME_CAL
        time_label[5] = clock();
#endif
        noise_subtraction(pgal, pnoise_1, &all_paras, 1, 1);

#ifdef TIME_CAL
        time_label[6] = clock();
#endif
        shear_est(pgal, ppsf, &all_paras);

#ifdef TIME_CAL
        time_label[7] = clock();
#endif
        total_data[j * shear_data_cols] = all_paras.n1;
        total_data[j * shear_data_cols + 1] = all_paras.n2;
        total_data[j * shear_data_cols + 2] = all_paras.dn;
        total_data[j * shear_data_cols + 3] = all_paras.du;
        //{std::cout<<j<<std::endl;}
#ifdef TIME_CAL
        time_label[8] = clock();
        for(i=0;i<10;i++)
        {
            time_collection[i] += time_label[i];
        }
#endif

    }

    gsl_free(0);
    gsl_free(1);
    gsl_free(2);
#ifdef TIME_CAL
    for(i=1;i<9;i++)
    {
        std::cout<<(time_collection[i]-time_collection[i-1])/CLOCKS_PER_SEC<<std::endl;
    }
#endif
    t2 = clock();
    get_time(time_now, 50);
    std::cout<<time_now<<" Total: "<<(t2-t1)/CLOCKS_PER_SEC<<std::endl;



	delete[] point;
	delete[] gal;
	delete[] pgal;
	delete[] psf;
	delete[] ppsf;
	delete[] noise_1;
	delete[] pnoise_1;
    delete[] noise_2;
	delete[] pnoise_2;
    delete[] total_data;


	return 0;
}
