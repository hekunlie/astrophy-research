#include<FQlib.h>
#include<hk_mpi.h>
#include<hk_iolib.h>
#define IMG_CHECK_LABEL 2
//#define FLUX_PDF

void arr_add(double *arr1, const double*arr2,const int length)
{
    for(int i=0;i<length;i++)
    {
        arr1[i] += arr2[i];
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

	char parent_path[100], chip_path[150], para_path[150], shear_path[150], result_path[150], log_path[150];
	char buffer[200], log_inform[250], set_name[50];

	//sprintf(parent_path, "/mnt/perc/hklee/bias_check");
	//sprintf(parent_path, "/mnt/ddnfs/data_users/hkli/bias_check/data_from_pi/step_test");
	sprintf(parent_path, "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/selection_bias/bias_check/step_test");

	//strcpy(parent_path, argv[1]);
	std::string str_data_path,str_shear_path;

	char_to_str(parent_path, str_data_path);
	str_shear_path = str_data_path + "/parameters/shear.dat";
	
	int i, j, k, ib;
	int sss1, sss2, seed_pts, seed_n1, seed_n2, seed_step;
	int rotation;

	int num_p, size, shear_pairs;
    double max_radius;
	int total_chips, sub_chip_num, sub_data_row, total_data_row;
	int stamp_num, stamp_nx, shear_data_cols;		
	int row, chip_st, chip_ed, shear_id, psf_type, temp_s, detect_label;

	double psf_scale, psf_thresh_scale, sig_level, psf_noise_sig, gal_noise_sig, flux_i, mag_i;
	double g1, g2, ts, te, t1, t2;
	double psf_ellip, psf_ang, psf_norm_factor;

	double pts_step;

	pts_step = atof(argv[1]);
	//rotation = atoi(argv[1]);
	rotation = 0;
	num_p = 60;
	max_radius= 10;
	stamp_num = 10000;
	shear_data_cols = 4;

    shear_pairs = 10;

	psf_type=2;
	psf_scale = 4;
	psf_thresh_scale = 2.;
	temp_s = rank;
	sig_level = 1.5;
	psf_noise_sig = 0;
    gal_noise_sig = 60;

	size = 48;
    total_chips = 1000;
    stamp_nx = 100;

	all_paras.stamp_size = size;
	all_paras.img_x = size;
	all_paras.img_y = size;
	// because the max half light radius of the galsim source is 5.5 pixels
	all_paras.max_distance = max_radius; 

	sprintf(log_path, "%s/logs/%02d_R%d_S%.2f.dat", parent_path, rank, rotation, pts_step);

#ifdef IMG_CHECK_LABEL	
	double *big_img_noisy = new double[stamp_nx*stamp_nx*size*size]{};
	double *big_img_noise_free = new double[stamp_nx*stamp_nx*size*size]{};

#endif

	double *point = new double[2 * num_p]{};

	double *gal = new double[size*size]{};
	double *pgal = new double[size*size]{};

	double *psf = new double[size*size]{};
	double *ppsf = new double[size*size]{};

	double *noise_1 = new double[size*size]{};
	double *pnoise_1 = new double[size*size]{};

    double *noise_2 = new double[size*size]{};
	double *pnoise_2 = new double[size*size]{};

	double *shear = new double[2 * shear_pairs]{}; // [...g1,...,..g2,...]

	// the shear estimators data matrix  
	double *total_data,*sub_noisy_data;
	double *sub_noise_free_data;
	double *sub_noise_residual_data;

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
	sub_noisy_data = new double[gather_count[rank]]{};
	sub_noise_free_data = new double[gather_count[rank]]{};
	sub_noise_residual_data = new double[gather_count[rank]]{};

	// seed distribution, different thread gets different seed
	seed_step = 1;
	sss1 = 1*seed_step*shear_pairs*scatter_count[0]/stamp_num;
	seed_pts = sss1*rank + 1 + 350000;
	seed_n1 = sss1*rank + 1 + 400100*(rotation+1);
	seed_n2 = sss1*rank + 1 + 2300100*(rotation+1);

	if (0 == rank)
	{
#ifdef FLUX_PDF
		total_flux = new double[total_data_row]{};
#endif
		total_data = new double[total_data_row*shear_data_cols]{};
	}

#ifdef FLUX_PDF
	sub_flux = new double[scatter_count[rank]]{};
#endif
	// read shear
	read_text(str_shear_path, shear, 2*shear_pairs);
    
	// create PSF
	create_psf(psf, psf_scale, size, psf_type);
	pow_spec(psf, ppsf, size, size);
	get_psf_radius(ppsf, &all_paras, psf_thresh_scale);

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
		std::cout << "PSF THRESH: " << all_paras.psf_pow_thresh << " PSF HLR: " << all_paras.psf_hlr << std::endl;
		std::cout <<"MAX RADIUS: "<< max_radius <<" ,Step: "<<pts_step<< ", SIG_LEVEL: " << sig_level <<"sigma"<< std::endl;
		sprintf(buffer, "!%s/psf.fits", parent_path);
		write_fits(buffer,psf, size, size);
		std::cout<<"Chip Num of each thread: ";
		show_arr(scatter_count,1,numprocs);
		std::cout<<"---------------------------------------------------------------------------"<<std::endl<<std::endl;
	}
	for(i=0;i<numprocs;i++)
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

		g1 = shear[shear_id];
		g2 = shear[shear_id + shear_pairs];

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
#ifdef FLUX_PDF
			sprintf(para_path, "%s/parameters/para_%d.hdf5", parent_path, shear_id);
			sprintf(set_name, "/flux");
			read_h5(para_path, set_name, total_flux);
#endif
		}
#ifdef FLUX_PDF
		my_Scatterv(total_flux, scatter_count, sub_flux, numprocs, rank);
#endif
		MPI_Barrier(MPI_COMM_WORLD);

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


#ifdef IMG_CHECK_LABEL
			initialize_arr(big_img_noisy, stamp_nx*stamp_nx*size*size, 0);
			initialize_arr(big_img_noise_free, stamp_nx*stamp_nx*size*size, 0);
#endif

			row = (i-chip_st)*stamp_num*shear_data_cols;

			// loop the stamps
			for (j = 0; j < stamp_num; j++)
			{	
#ifdef FLUX_PDF
				flux_i = sub_flux[(i-chip_st)*stamp_num + j] / num_p;
#else	
				flux_i = 9000./ num_p;
#endif
				// ////////////////// pure noise image //////////////////////
				// initialize_arr(noise_1, size*size, 0);
				// initialize_arr(pnoise_1, size*size, 0);
                // initialize_arr(noise_2, size*size, 0);
				// initialize_arr(pnoise_2, size*size, 0);

				// addnoise(noise_1, size*size, gal_noise_sig, rng1);
				// addnoise(noise_2, size*size, gal_noise_sig, rng2);

                // pow_spec(noise_1, pnoise_1, size, size);
                // pow_spec(noise_2, pnoise_2, size, size);


                // // measure on the noise power spectrum residual
                // noise_subtraction(pnoise_1, pnoise_2, &all_paras, 1, 0);
                // shear_est(pnoise_1, ppsf, &all_paras);

                // sub_noise_residual_data[row + j * shear_data_cols] = all_paras.n1;
				// sub_noise_residual_data[row + j * shear_data_cols + 1] = all_paras.n2;
				// sub_noise_residual_data[row + j * shear_data_cols + 2] = all_paras.dn;
				// sub_noise_residual_data[row + j * shear_data_cols + 3] = all_paras.du;


				///////////////// Noise free /////////////////////////		
				initialize_arr(point, num_p * 2, 0);
				
				initialize_arr(gal, size*size, 0);
				initialize_arr(pgal, size*size, 0);
				
				create_points(point, num_p, max_radius, pts_step, rng0);
				convolve(gal, point, flux_i, size, num_p, rotation, psf_scale, g1, g2, psf_type, 1, &all_paras);
#ifdef IMG_CHECK_LABEL
				stack(big_img_noise_free, gal, j, size, stamp_nx, stamp_nx);
				//std::cout<<flux_i<<" "<<sub_flux[(i-chip_st)*stamp_num + j]<<std::endl;
#endif			
				pow_spec(gal, pgal, size, size);
				shear_est(pgal, ppsf, &all_paras);

				sub_noise_free_data[row + j * shear_data_cols] = all_paras.n1;
				sub_noise_free_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_noise_free_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_noise_free_data[row + j * shear_data_cols + 3] = all_paras.du;


                ////////////////// noisy image //////////////////////
                initialize_arr(pgal, size*size, 0);
				initialize_arr(noise_2, size*size, 0);
				initialize_arr(pnoise_2, size*size, 0);

				addnoise(gal, size*size, gal_noise_sig, rng1);
				pow_spec(gal, pgal, size, size);

				addnoise(noise_2, size*size, gal_noise_sig, rng1);
                pow_spec(noise_2, pnoise_2, size, size);

				noise_subtraction(pgal, pnoise_2, &all_paras, 1, 1);

#ifdef IMG_CHECK_LABEL
				stack(big_img_noisy, gal, j, size, stamp_nx, stamp_nx);
				//std::cout<<flux_i<<" "<<sub_flux[(i-chip_st)*stamp_num + j]<<std::endl;
#endif

				shear_est(pgal, ppsf, &all_paras);

				sub_noisy_data[row + j * shear_data_cols] = all_paras.n1;
				sub_noisy_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_noisy_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_noisy_data[row + j * shear_data_cols + 3] = all_paras.du;

			}
			gsl_free(0);
			gsl_free(1);
			gsl_free(2);


#ifdef IMG_CHECK_LABEL
			if(i == IMG_CHECK_LABEL)
			{
				sprintf(chip_path, "!%s/%d/gal_chip_%05d_r%d_step_%.2f_noise_free.fits", parent_path, shear_id, i, rotation, pts_step);
				write_fits(chip_path, big_img_noise_free, stamp_nx*size, stamp_nx*size);
				sprintf(chip_path, "!%s/%d/gal_chip_%05d_r%d_step_%.2f_noisy.fits", parent_path, shear_id, i, rotation, pts_step);
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
			sprintf(result_path, "%s/result/data_%.2f/data_%d_noise_free.hdf5", parent_path, pts_step, shear_id);
			write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols, true);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// my_Gatherv(sub_noise_residual_data, gather_count, total_data, numprocs, rank);
		// if (0 == rank)
		// {
		// 	sprintf(result_path, "%s/result/data_%.2f/data_%d_noise_residual.hdf5", parent_path, pts_step, shear_id);
		// 	write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols, true);
		// }
		// MPI_Barrier(MPI_COMM_WORLD);

		my_Gatherv(sub_noisy_data, gather_count, total_data, numprocs, rank);
		if (0 == rank)
		{
			sprintf(result_path, "%s/result/data_%.2f/data_%d_noisy_cpp.hdf5", parent_path, pts_step, shear_id);
			write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols, true);
			std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		te = clock();
		if (rank == 0)
		{
			sprintf(buffer, "rank %d: done in %g \n", rank, (te - ts) / CLOCKS_PER_SEC);
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
	delete[] sub_noisy_data;
	delete[] shear;
#ifdef FLUX_PDF
	delete[] sub_flux;
#endif
	delete[] scatter_count;
	delete[] gather_count;
	
	MPI_Finalize();
	return 0;
}
