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

	para all_paras;

	char parent_path[100], chip_path[150], para_path[150], shear_path[150], result_path[150], log_path[150];
	char buffer[200], log_inform[250], set_name[50];

	sprintf(parent_path, "/mnt/ddnfs/data_users/hkli/bias_check");

    int i, j, k, ib;

	int shear_pairs;

	int total_chips, sub_chip_num, sub_data_row, total_data_row;
	int stamp_size,stamp_num, stamp_nx, shear_data_cols;		
	int row, chip_st, chip_ed, shear_id, temp_s;
    double psf_scale, psf_thresh_scale;
    int psf_type;

    psf_type = 2;
    psf_scale = 4;
    psf_thresh_scale = 2;

    total_chips = 1000;
    stamp_size = 64;
    stamp_nx = 100;
    stamp_num = stamp_nx*stamp_nx;
	shear_data_cols = 12;
    total_data_row = total_chips*stamp_num;

    all_paras.stamp_size = stamp_size;
	all_paras.img_x = stamp_size;
	all_paras.img_y = stamp_size;

    double *psf_img = new double[stamp_size*stamp_size]{};
	double *psf_img_pow = new double[stamp_size*stamp_size]{};    

	double *noise_img1 = new double[stamp_size*stamp_size]{};
	double *noise_img2 = new double[stamp_size*stamp_size]{};
	double *noise_img1_pow = new double[stamp_size*stamp_size]{};
	double *noise_img2_pow = new double[stamp_size*stamp_size]{};
    
    double *big_img1 = new double[stamp_nx*stamp_nx*stamp_size*stamp_size]{};
    double *big_img2 = new double[stamp_nx*stamp_nx*stamp_size*stamp_size]{};

    double *total_data,*sub_data;
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
	sub_data = new double[gather_count[rank]]{};
	if (0 == rank)
	{
		total_data = new double[total_data_row*shear_data_cols]{};
	}


    // create PSF
	create_psf(psf_img, psf_scale, stamp_size, psf_type);
	pow_spec(psf_img, psf_img_pow, stamp_size, stamp_size);
	get_psf_radius(psf_img_pow, &all_paras, psf_thresh_scale);
    
	for(i=0;i<numprocs;i++)
	{
		if(i == rank)
		{
			std::cout<<rank<<" "<<chip_st<<" "<<chip_ed<<" "<<scatter_count[i]<<" "<<gather_count[i]<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

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

			sprintf(log_inform, "RANK: %03d, SHEAR %02d:, chip: %05d, start. seed:%d", rank,shear_id, i, seed);
			if (rank == 0)
			{
				std::cout << log_inform << std::endl;
			}

            // read noise images
            initialize_arr(big_img1, stamp_nx*stamp_nx*stamp_size*stamp_size);
            initialize_arr(big_img2, stamp_nx*stamp_nx*stamp_size*stamp_size);

            sprintf(chip_path,"%s/%d/chip_%d.hdf5",parent_path, shear_id, i);
            sprintf(set_name,"/noise_1");
            read_h5(chip_path, set_name, big_img1);
            sprintf(set_name,"/noise_2");
            read_h5(chip_path, set_name, big_img2);


            row = (i-chip_st)*stamp_num*shear_data_cols;

            for (j = 0; j < stamp_num; j++)
			{   
                initialize_arr(noise_img1,stamp_size*stamp_size,0);
                initialize_arr(noise_img2,stamp_size*stamp_size,0);
                initialize_arr(noise_img1_pow,stamp_size*stamp_size,0);
                initialize_arr(noise_img2_pow,stamp_size*stamp_size,0);

                segment(big_img1, noise_img1, j, stamp_size, stamp_nx, stamp_nx);
                segment(big_img2, noise_img2, j, stamp_size, stamp_nx, stamp_nx);

                pow_spec(noise_img1, noise_img1_pow, stamp_size, stamp_size);
                pow_spec(noise_img2, noise_img2_pow, stamp_size, stamp_size);
                
            }

        }
    

    }
    if(rank == 0)
    {
        show_arr(scatter_count,1,numprocs);
        show_arr(gather_count,1,numprocs);

        delete[] total_data;
    }
    delete[] psf_img;
    delete[] psf_img_pow;
    delete[] noise_img1;
    delete[] noise_img1_pow;
    delete[] noise_img2;
    delete[] noise_img2_pow;
    delete[] big_img1;
    delete[] big_img2;
    delete[] sub_data;
   	MPI_Finalize();
	return 0;
}