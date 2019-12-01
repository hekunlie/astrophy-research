#include<FQlib.h>
#include<hk_iolib.h>
#include<hk_mpi.h>

int main(int argc, char**argv)
{
    /* it is designed to run on the cluster (should save the memory) */
	/* measure the FQ shear estimators on the noise powerspectrum residual */
    /* two PSF's will be used, Gaussian & Moffat */
	/* loop the shear points, the chips will be distributed to threads */

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	fq_paras all_paras_1;
	fq_paras all_paras_2;

	char parent_path[100], chip_path[150], para_path[150], shear_path[150], result_path[150], log_path[150];
	char buffer[200], log_inform[250], set_name[50];

	sprintf(parent_path, "/mnt/ddnfs/data_users/hkli/noise_img");

    int i, j, k, ib;
	double ts,te, t1, t2;
	int shear_pairs;

	int total_chips, sub_chip_num, sub_data_row, total_data_row;
	int stamp_size, stamp_num, stamp_nx, shear_data_cols, img_len;		
	int row, chip_st, chip_ed, shear_id, temp_s;
    double psf_scale, psf_thresh_scale;
    int psf_type;

    psf_type = 2;
    psf_scale = 4;
    psf_thresh_scale = 2;

    total_chips = 1000;
    stamp_size = 64;
	img_len = stamp_size*stamp_size;
    stamp_nx = 100;
    stamp_num = stamp_nx*stamp_nx;
	shear_data_cols = 24;
    total_data_row = total_chips*stamp_num;

    all_paras_1.stamp_size = stamp_size;
	all_paras_1.img_x = stamp_size;
	all_paras_1.img_y = stamp_size;

	all_paras_2.stamp_size = stamp_size;
	all_paras_2.img_x = stamp_size;
	all_paras_2.img_y = stamp_size;

    double *gauss_psf_img = new double[img_len]{};
	double *gauss_psf_img_pow = new double[img_len]{};    
    double *moffat_psf_img = new double[img_len]{};
	double *moffat_psf_img_pow = new double[img_len]{};    

	double *noise_img1 = new double[img_len]{};
	double *noise_img2 = new double[img_len]{};
	double *noise_img1_pow = new double[img_len]{};
	double *noise_img2_pow = new double[img_len]{};
    double *noise_pow_residual = new double[img_len]{};

    double *big_img1 = new double[stamp_nx*stamp_nx*img_len]{};
    double *big_img2 = new double[stamp_nx*stamp_nx*img_len]{};

    double *total_data,*sub_data;
	double *total_flux, *sub_flux;
	int *gather_count;

	// for the gatherv when finish in each shear point
	gather_count = new int[numprocs]{};

	///////////////////// task distribution /////////////////////////////////////
	sub_chip_num = total_chips / numprocs;
	j = total_chips%numprocs;
	for(i=0;i<numprocs;i++)
	{
		gather_count[i] = sub_chip_num;
	}
	for(i=0;i<j;i++)
	{	
		// the chip number
		gather_count[i] += 1;
	}
	sub_chip_num = gather_count[rank];
	// the start- & end-label of chip of each thread
	chip_st = 0;
	for(i=0;i<rank;i++)
	{
		chip_st += gather_count[i];
	}
	chip_ed = chip_st+gather_count[rank];
	// the final data from all the source in one shear point
	total_data_row = total_chips * stamp_num;
	// the sub-data from the source processed by each thread
	sub_data_row = sub_chip_num * stamp_num;
	for(i=0;i<numprocs;i++)
	{
		// the real amount of data of each thread
		gather_count[i] = gather_count[i]*shear_data_cols*stamp_num;
	}
	sub_data = new double[gather_count[rank]]{};
	if (0 == rank)
	{
		total_data = new double[total_data_row*shear_data_cols]{};
	}


    // create Gaussian PSF
	create_psf(gauss_psf_img, psf_scale, stamp_size, 0);
	pow_spec(gauss_psf_img, gauss_psf_img_pow, stamp_size, stamp_size);
	get_psf_radius(gauss_psf_img_pow, &all_paras_1, psf_thresh_scale);
    // Moffat PSF
	create_psf(moffat_psf_img, psf_scale, stamp_size, 1);
	pow_spec(moffat_psf_img, moffat_psf_img_pow, stamp_size, stamp_size);
	get_psf_radius(moffat_psf_img_pow, &all_paras_2, psf_thresh_scale);

	if(rank==0)show_arr(gather_count,1,numprocs);
	for(i=0;i<numprocs;i++)
	{
		if(i == rank)
		{
			std::cout<<rank<<" "<<chip_st<<" "<<chip_ed<<" "<<gather_count[i]<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
    for (shear_id = 0; shear_id < shear_pairs; shear_id++)
	{
		ts = clock();
        sprintf(log_inform, "RANK: %03d, SHEAR %02d: my chips: %d - %d", rank, shear_id, chip_st, chip_ed);

		if(rank == 0)
		{
			std::cout<<"---------------------------------------------------------------------------"<<std::endl;
			std::cout<<log_inform<<std::endl;
		}
        // loop chips
        for (i = chip_st; i < chip_ed; i++)
		{
			t1 = clock();

			sprintf(log_inform, "RANK: %03d, SHEAR %02d: chip: %05d, start.", rank,shear_id, i);
			if (rank == 0)
			{
				std::cout << log_inform << std::endl;
			}

            // read noise images
            initialize_arr(big_img1, stamp_nx*stamp_nx*img_len,0);
            initialize_arr(big_img2, stamp_nx*stamp_nx*img_len,0);

            sprintf(chip_path,"%s/%d/chip_%d.hdf5",parent_path, shear_id, i);
            sprintf(set_name,"/noise_1");
            read_h5(chip_path, set_name, big_img1);
            sprintf(set_name,"/noise_2");
            read_h5(chip_path, set_name, big_img2);


            row = (i-chip_st)*stamp_num*shear_data_cols;

            for (j = 0; j < stamp_num; j++)
			{   
                initialize_arr(noise_img1,img_len,0);
                initialize_arr(noise_img2,img_len,0);
                initialize_arr(noise_img1_pow,img_len,0);
                initialize_arr(noise_img2_pow,img_len,0);

                segment(big_img1, noise_img1, j, stamp_size, stamp_nx, stamp_nx);
                segment(big_img2, noise_img2, j, stamp_size, stamp_nx, stamp_nx);

                pow_spec(noise_img1, noise_img1_pow, stamp_size, stamp_size);
                pow_spec(noise_img2, noise_img2_pow, stamp_size, stamp_size);
				for(k=0;k<img_len;k++)
				{
					noise_pow_residual[k] = noise_img1_pow[k] - noise_img2_pow[k];
				}

				// Gaussian PSF
				shear_est(noise_img1_pow, gauss_psf_img_pow, &all_paras_1);
				sub_data[row + j*shear_data_cols] = all_paras_1.n1;
				sub_data[row + j*shear_data_cols + 1] = all_paras_1.n2;
				sub_data[row + j*shear_data_cols + 2] = all_paras_1.dn;
				sub_data[row + j*shear_data_cols + 3] = all_paras_1.du;

				shear_est(noise_img2_pow, gauss_psf_img_pow, &all_paras_1);
				sub_data[row + j*shear_data_cols + 4] = all_paras_1.n1;
				sub_data[row + j*shear_data_cols + 5] = all_paras_1.n2;
				sub_data[row + j*shear_data_cols + 6] = all_paras_1.dn;
				sub_data[row + j*shear_data_cols + 7] = all_paras_1.du;

                shear_est(noise_pow_residual, gauss_psf_img_pow, &all_paras_1);
				sub_data[row + j*shear_data_cols + 8] = all_paras_1.n1;
				sub_data[row + j*shear_data_cols + 9] = all_paras_1.n2;
				sub_data[row + j*shear_data_cols + 10] = all_paras_1.dn;
				sub_data[row + j*shear_data_cols + 11] = all_paras_1.du;

				// Moffat PSF
				shear_est(noise_img1_pow, moffat_psf_img_pow, &all_paras_2);
				sub_data[row + j*shear_data_cols + 12] = all_paras_2.n1;
				sub_data[row + j*shear_data_cols + 13] = all_paras_2.n2;
				sub_data[row + j*shear_data_cols + 14] = all_paras_2.dn;
				sub_data[row + j*shear_data_cols + 15] = all_paras_2.du;

				shear_est(noise_img2_pow, moffat_psf_img_pow, &all_paras_2);
				sub_data[row + j*shear_data_cols + 16] = all_paras_2.n1;
				sub_data[row + j*shear_data_cols + 17] = all_paras_2.n2;
				sub_data[row + j*shear_data_cols + 18] = all_paras_2.dn;
				sub_data[row + j*shear_data_cols + 19] = all_paras_2.du;

                shear_est(noise_pow_residual, moffat_psf_img_pow, &all_paras_2);
				sub_data[row + j*shear_data_cols + 20] = all_paras_2.n1;
				sub_data[row + j*shear_data_cols + 21] = all_paras_2.n2;
				sub_data[row + j*shear_data_cols + 22] = all_paras_2.dn;
				sub_data[row + j*shear_data_cols + 23] = all_paras_2.du;
            }

			t2 = clock();
			sprintf(log_inform, "RANK: %03d, SHEAR %02d: chip: %05d, done in %.2f s.", rank, shear_id, i, (t2 - t1) / CLOCKS_PER_SEC);
			if (rank == numprocs-1)
			{
				std::cout << log_inform << std::endl;
			}

        }

		te = clock();
		if (rank == numprocs-1)
		{
			sprintf(log_inform, "rank %d: done in %g \n", rank, (te - ts) / CLOCKS_PER_SEC);
			std::cout << log_inform<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		my_Gatherv(sub_data, gather_count, total_data,numprocs, rank);
		if (0 == rank)
		{
			sprintf(set_name, "/data");
			sprintf(result_path, "%s/result/data/data_%d.hdf5", parent_path, shear_id);
			write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols, true);
			std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		}

    }
    if(rank == 0)
    {
        delete[] total_data;
    }
    delete[] gauss_psf_img;
    delete[] gauss_psf_img_pow;
	delete[] moffat_psf_img;
    delete[] moffat_psf_img_pow;
    delete[] noise_img1;
    delete[] noise_img1_pow;
    delete[] noise_img2;
    delete[] noise_img2_pow;
    delete[] noise_pow_residual;
    delete[] big_img1;
    delete[] big_img2;
    delete[] sub_data;
   	MPI_Finalize();
	return 0;
}