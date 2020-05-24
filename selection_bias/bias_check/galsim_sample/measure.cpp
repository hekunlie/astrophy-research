#include<stdlib.h>
#include<hk_mpi.h>
#include<hk_iolib.h>
#include<FQlib.h>

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

void arr_deduct(double *arr1, const double*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] -= arr2[i];
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

void arr_sep(double *data, float *mg1,float *mg2, float *mn, float *mu,float *mv, int data_num, int data_col)
{
	for(int i=0; i<data_num; i++)
	{
		mg1[i] = data[i*data_col];
		mg2[i] = data[i*data_col + 1];
		mn[i] = data[i*data_col + 2];
		mu[i] = data[i*data_col + 3];
		mv[i] = data[i*data_col + 4];
	}
}
int main(int argc, char*argv[])
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	fq_paras all_paras;

	char parent_path[300], data_path[300], result_path[300], chip_path[300], set_name[50];
	char log_path[300], log_inform[300];

	int i, j, k, row;
	int seed_ini, seed, seed_step;
	double ts, te, t1, t2;

	int size, img_len;
	int total_chips, chip_num, chip_id_s, chip_id_e;
	int total_data_row, total_data_cols, sub_data_row;
	int stamp_num, stamp_nx;
	
	int shear_id, shear_pairs;
	double psf_thresh_scale, gal_noise_sig;

	float *mg[5];
	char *mg_name[5];

	strcpy(parent_path,argv[1]);
	total_chips = atoi(argv[2]);
	shear_pairs = atoi(argv[3]);
	size = atoi(argv[4]);
	seed_ini = atoi(argv[5]);

	total_data_cols = 5;
	stamp_nx = 100;
	stamp_num = stamp_nx*stamp_nx;

	psf_thresh_scale = 2.;
	gal_noise_sig = 60;

	img_len = size*size;
	all_paras.stamp_size = size;
	all_paras.img_x = size;
	all_paras.img_y = size;

	chip_num = total_chips / numprocs;
	total_data_row = total_chips * stamp_num;
	sub_data_row = chip_num * stamp_num;

	chip_id_s = chip_num * rank;
	chip_id_e = chip_num * (rank + 1);

	sprintf(log_path,"%s/logs/mlog_%d.dat", parent_path, rank);
	if (0 == rank)
	{	
		std::cout<<parent_path<<std::endl;
		std::cout << "Total chip: " << total_chips<<", Total cpus: "<<numprocs <<", Stamp size: "<<size <<std::endl;
		sprintf(log_inform, "RANK: %03d,  thread: %d, total cpus: %d, individual chip: %d , size：%d, stamp_col: %d", rank, numprocs, total_chips, chip_num, size, stamp_nx);
		std::cout<<log_inform<<std::endl;
		write_log(log_path, log_inform);
	}


	double *psf = new double[img_len]{};
	double *ppsf = new double[img_len]{};
	double *ppsf_sqrt = new double[img_len]{};

	double *gal = new double[img_len]{};
	double *pgal = new double[img_len]{};
	double *pgal_dn = new double[img_len]{};
	
	double *gal_noisy = new double[img_len]{};
	double *pgal_noisy = new double[img_len]{};

	double *gal_noisy_a = new double[img_len]{};
	double *pgal_noisy_a = new double[img_len]{};

	double *pgal_cross_term = new double[img_len]{};
	double *pgal_cross_term_est = new double[img_len]{};

	
	double *noise_1 = new double[img_len]{};
	double *pnoise_1 = new double[img_len]{};

    double *noise_2 = new double[img_len]{};
	double *pnoise_2 = new double[img_len]{};

	double *noise_3 = new double[img_len]{};
	double *pnoise_3 = new double[img_len]{};

	double *pnoise_cross = new double[img_len]{};
    double *noise_pow_diff_1 = new double[img_len]{};
	double *noise_pow_diff_2 = new double[img_len]{};

	double *big_img = new double[img_len*stamp_num]{};

	double *recvbuf;
	
	double *sub_noisy_data = new double[sub_data_row*total_data_cols]{};
	double *sub_noise_free_data = new double[sub_data_row*total_data_cols]{};
	double *sub_noise_residual_data = new double[sub_data_row*total_data_cols]{};
	double *sub_noise_residual_estimate = new double[sub_data_row*total_data_cols]{};

	double *sub_cross_term_data = new double[sub_data_row*total_data_cols]{};
	double *sub_cross_term_estimate = new double[sub_data_row*total_data_cols]{};
	double *sub_cross_term_sqrt_data = new double[sub_data_row*total_data_cols]{};

	int *gather_count = new int[numprocs]{};
	for(i=0;i<numprocs;i++){gather_count[i] = sub_data_row*total_data_cols;}

	if (0 == rank)
	{	
		// for(i=0;i<5;i++)
		// {
		// 	mg[i] = new float[total_chips*total_data_row];
		// 	mg_name[i] = new char[50];
		// }
		// sprintf(mg_name[0], "/mg1");
		// sprintf(mg_name[1], "/mg2");
		// sprintf(mg_name[2], "/mn");
		// sprintf(mg_name[3], "/mu");
		// sprintf(mg_name[4], "/mv");
		recvbuf = new double[total_data_row*total_data_cols];
	}
	sprintf(data_path, "%s/psf.fits", parent_path);
	read_fits(data_path, psf);
	pow_spec(psf, ppsf, size, size);
	arr_sqrt(ppsf, ppsf_sqrt, img_len);
	get_psf_radius(ppsf, &all_paras, psf_thresh_scale);

	seed_step = 2;
	i = 3*seed_step*shear_pairs*(total_chips/numprocs+1);
	seed = i*rank + 1 + seed_ini;

	if (0 == rank)
	{
		std::cout << "PSF THRESH: " << all_paras.psf_pow_thresh <<" "<< all_paras.psf_hlr << std::endl;
	}

	for (shear_id = 0; shear_id < shear_pairs; shear_id++)
	{
		ts = clock();
		sprintf(log_inform, "RANK: %03d, SHEAR %02d: my chips: %d - %d, total chips: %d (%d cpus)", rank, shear_id, chip_id_s, chip_id_e, total_chips, numprocs);
		write_log(log_path, log_inform);
		
		if (0 == rank)
		{
			std::cout<<log_inform<<std::endl;
			initialize_arr(recvbuf, total_data_row*total_data_cols, 0);
		}

		for (i = chip_id_s; i < chip_id_e; i++)
		{
			t1 = clock();
			sprintf(log_inform, "RANK: %03d, SHEAR %02d: %04d 's chip start. Seed: %d", rank, shear_id, i, seed);
			write_log(log_path, log_inform);
			if (0 == rank)
			{
				std::cout << log_inform << std::endl;
			}
			gsl_initialize(seed, 1);
			seed += seed_step;

			sprintf(chip_path, "/imgs/%s/%d/chip_%04d.fits", parent_path, shear_id, i);
			initialize_arr(big_img, stamp_num*img_len, 0);
			
			read_fits(chip_path, big_img);

			row = (i - chip_id_s) *stamp_num*total_data_cols;
			
			for (j = 0; j < stamp_num; j++)			
			{	
				initialize_arr(gal, img_len, 0);
				segment(big_img, gal, j, size, stamp_nx, stamp_nx);	

				initialize_arr(pgal, size*size, 0);
				
				initialize_arr(noise_1, img_len, 0);
				initialize_arr(pnoise_1, img_len, 0);
				initialize_arr(noise_2, img_len, 0);
				initialize_arr(pnoise_2, img_len, 0);
				initialize_arr(noise_3, img_len, 0);
				initialize_arr(pnoise_3, img_len, 0);

				initialize_arr(pgal_cross_term, img_len, 0);
				initialize_arr(pgal_noisy, img_len, 0);

				addnoise(noise_1, img_len, gal_noise_sig, rng1);
				pow_spec(noise_1, pnoise_1, size, size);
				addnoise(noise_2, img_len, gal_noise_sig, rng1);
				pow_spec(noise_2, pnoise_2, size, size);
				addnoise(noise_3, img_len, gal_noise_sig, rng1);
				pow_spec(noise_3, pnoise_3, size, size);

				// noise free
				pow_spec(gal, pgal, size, size);
				// noisy image
				arr_add(gal_noisy, gal, noise_1,img_len);
				pow_spec(gal_noisy, pgal_noisy, size, size);
				arr_deduct(pgal_dn, pgal_noisy, pnoise_2, img_len);
				// noise residual
				arr_deduct(noise_pow_diff_1, pnoise_1, pnoise_2, img_len);
				arr_deduct(noise_pow_diff_2, pnoise_2, pnoise_3, img_len);
				// cross term
                arr_deduct(pgal_cross_term, pgal_noisy, pgal, pnoise_1, img_len);
				// cross term estimate
				arr_add(gal_noisy_a, gal_noisy, noise_2, img_len);
				pow_spec(gal_noisy_a, pgal_noisy_a, size, size);
				arr_deduct(pgal_cross_term_est, pgal_noisy_a, pgal_noisy, pnoise_2, img_len);

				///////////////// Noise free /////////////////////////
				shear_est(pgal, ppsf, &all_paras);
				sub_noise_free_data[row + j * total_data_cols] = all_paras.n1;
				sub_noise_free_data[row + j * total_data_cols + 1] = all_paras.n2;
				sub_noise_free_data[row + j * total_data_cols + 2] = all_paras.dn;
				sub_noise_free_data[row + j * total_data_cols + 3] = all_paras.du;
				sub_noise_free_data[row + j * total_data_cols + 4] = all_paras.dv;

				/////////////////////////// noisy image /////////////////////////////////
 				shear_est(pgal_dn, ppsf, &all_paras);
				sub_noisy_data[row + j * total_data_cols] = all_paras.n1;
				sub_noisy_data[row + j * total_data_cols + 1] = all_paras.n2;
				sub_noisy_data[row + j * total_data_cols + 2] = all_paras.dn;
				sub_noisy_data[row + j * total_data_cols + 3] = all_paras.du;
				sub_noisy_data[row + j * total_data_cols + 4] = all_paras.dv;

 				////////////////// noise-residual image //////////////////////
				shear_est(noise_pow_diff_1, ppsf, &all_paras);
				sub_noise_residual_data[row + j * total_data_cols] = all_paras.n1;
				sub_noise_residual_data[row + j * total_data_cols + 1] = all_paras.n2;
				sub_noise_residual_data[row + j * total_data_cols + 2] = all_paras.dn;
				sub_noise_residual_data[row + j * total_data_cols + 3] = all_paras.du;
				sub_noise_residual_data[row + j * total_data_cols + 4] = all_paras.dv;

 				////////////////// noise-residual image //////////////////////
				shear_est(noise_pow_diff_2, ppsf, &all_paras);
				sub_noise_residual_estimate[row + j * total_data_cols] = all_paras.n1;
				sub_noise_residual_estimate[row + j * total_data_cols + 1] = all_paras.n2;
				sub_noise_residual_estimate[row + j * total_data_cols + 2] = all_paras.dn;
				sub_noise_residual_estimate[row + j * total_data_cols + 3] = all_paras.du;
				sub_noise_residual_estimate[row + j * total_data_cols + 4] = all_paras.dv;

				////////////////// cross-term image //////////////////////
				shear_est(pgal_cross_term, ppsf, &all_paras);
				sub_cross_term_data[row + j * total_data_cols] = all_paras.n1;
				sub_cross_term_data[row + j * total_data_cols + 1] = all_paras.n2;
				sub_cross_term_data[row + j * total_data_cols + 2] = all_paras.dn;
				sub_cross_term_data[row + j * total_data_cols + 3] = all_paras.du;
				sub_cross_term_data[row + j * total_data_cols + 4] = all_paras.dv;

				////////////////// cross-term divided by |P_psf| image //////////////
				shear_est(pgal_cross_term, ppsf_sqrt, &all_paras);
				sub_cross_term_sqrt_data[row + j * total_data_cols] = all_paras.n1;
				sub_cross_term_sqrt_data[row + j * total_data_cols + 1] = all_paras.n2;
				sub_cross_term_sqrt_data[row + j * total_data_cols + 2] = all_paras.dn;
				sub_cross_term_sqrt_data[row + j * total_data_cols + 3] = all_paras.du;
				sub_cross_term_sqrt_data[row + j * total_data_cols + 4] = all_paras.dv;

				////////////////// cross-term-est image //////////////////////
				shear_est(pgal_cross_term_est, ppsf, &all_paras);				
				sub_cross_term_estimate[row + j * total_data_cols] = all_paras.n1;
				sub_cross_term_estimate[row + j * total_data_cols + 1] = all_paras.n2;
				sub_cross_term_estimate[row + j * total_data_cols + 2] = all_paras.dn;
				sub_cross_term_estimate[row + j * total_data_cols + 3] = all_paras.du;
				sub_cross_term_estimate[row + j * total_data_cols + 4] = all_paras.dv;
			}		

			t2 = clock();
			sprintf(log_inform, "RANK: %03d, SHEAR %02d: %04d 's chip finish in %.2f sec", rank, shear_id, i, (t2 - t1) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			if (0 == rank)
			{
				std::cout << log_inform << std::endl;
			}
			gsl_free(1);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		sprintf(set_name, "/data");
		my_Gatherv(sub_noise_free_data, gather_count, recvbuf, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data/data_noise_free_%d.hdf5", parent_path, shear_id);
			write_h5(result_path, set_name, recvbuf, total_data_row, total_data_cols, TRUE);
			// arr_sep(recvbuf, mg[0], mg[1], mg[2], mg[3], mg[4], total_data_row, total_data_cols)
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_noisy_data, gather_count, recvbuf, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data/data_noisy_cpp_%d.hdf5", parent_path, shear_id);
			write_h5(result_path, set_name, recvbuf, total_data_row, total_data_cols, TRUE);
		}

		MPI_Barrier(MPI_COMM_WORLD);


		my_Gatherv(sub_noise_residual_data, gather_count, recvbuf, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data/data_noise_residual_%d.hdf5", parent_path, shear_id);
			write_h5(result_path, set_name, recvbuf, total_data_row, total_data_cols, TRUE);
		}
		MPI_Barrier(MPI_COMM_WORLD);


		my_Gatherv(sub_noise_residual_estimate, gather_count, recvbuf, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data/data_noise_residual_est_%d.hdf5", parent_path, shear_id);
			write_h5(result_path, set_name, recvbuf, total_data_row, total_data_cols, TRUE);
		}
		MPI_Barrier(MPI_COMM_WORLD);


		my_Gatherv(sub_cross_term_data, gather_count, recvbuf, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data/data_cross_term_%d.hdf5", parent_path, shear_id);
			write_h5(result_path, set_name, recvbuf, total_data_row, total_data_cols, TRUE);
		}
		MPI_Barrier(MPI_COMM_WORLD);	

		my_Gatherv(sub_cross_term_estimate, gather_count, recvbuf, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data/data_cross_term_est_%d.hdf5", parent_path, shear_id);
			write_h5(result_path, set_name, recvbuf, total_data_row, total_data_cols, TRUE);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_cross_term_sqrt_data, gather_count, recvbuf, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/data/data_cross_term_sqrt_%d.hdf5", parent_path, shear_id);
			write_h5(result_path, set_name, recvbuf, total_data_row, total_data_cols, TRUE);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		te = clock();
		sprintf(log_inform, "RANK: %03d, SHEAR %02d: write data to file and finish jobs in %.2f sec", rank, shear_id, (te - ts) / CLOCKS_PER_SEC);
		write_log(log_path, log_inform);
		if (0 == rank)
		{
			std::cout << log_inform << std::endl;
			std::cout<<parent_path<<std::endl;
		}
	}

	if (0 == rank)
	{
		delete[] recvbuf;
	}		
	delete[] psf;
	delete[] ppsf;
	delete[] big_img;
	delete[] gal;
	delete[] pgal;

	MPI_Finalize();
	return 0;
}
