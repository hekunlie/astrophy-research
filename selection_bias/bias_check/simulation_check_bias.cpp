#include <stdlib.h>
#include <hk_mpi.h>
#include "FQlib.h"
#include<cstdio>
#include<hk_iolib.h>

#define IMG_CHECK_LABEL 2

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

	para all_paras;

	char parent_path[100], chip_path[150], para_path[150], shear_path[150], result_path[150], log_path[150];
	char buffer[200], log_inform[250], set_name[50];

	sprintf(parent_path, "/lustre/home/acct-phyzj/phyzj/hklee/work/selection_bias/bias_check/data");
	//strcpy(parent_path, argv[1]);
	std::ifstream fin;
	std::string str_stampsize = "stamp_size", str_total_num = "total_num", str_noise = "noise_sig", str_shear_num = "shear_num", str_nx = "stamp_col";
	std::string dect_info;
	std::string str_data_path,str_para_path,str_shear_path;


	int i, j, k, ib;
	int sss1, sss2, seed;

	int num_p, max_radius, size, shear_pairs;
	int total_chips, sub_chip_num, sub_data_row, total_data_row;
	int stamp_num, stamp_nx, shear_data_cols;		
	int row, chip_st, chip_ed, shear_id, psf_type, temp_s, detect_label;

	double psf_scale, psf_thresh_scale, sig_level, psf_noise_sig, gal_noise_sig, flux_i, mag_i;
	double g1, g2, ts, te, t1, t2;
	double psf_ellip, psf_ang, psf_norm_factor;

	num_p = 100;
	max_radius=9;
	stamp_num = 10000;
	shear_data_cols = 8;

	psf_type=2;
	psf_scale = 4;
	psf_thresh_scale = 2.;
	temp_s = rank;
	sig_level = 1.5;
	psf_noise_sig = 0;

	char_to_str(parent_path, str_data_path);
	str_para_path = str_data_path + "/parameters/para.ini";
	str_shear_path = str_data_path + "/parameters/shear.dat";
	sprintf(log_path, "%s/logs/%02d.dat", parent_path, rank);

	read_para(str_para_path, str_stampsize, size);
	read_para(str_para_path, str_total_num, total_chips);
	read_para(str_para_path, str_noise, gal_noise_sig);
	read_para(str_para_path, str_shear_num, shear_pairs);
	read_para(str_para_path, str_nx, stamp_nx);

	all_paras.stamp_size = size;
	all_paras.img_x = size;
	all_paras.img_y = size;
	// because the max half light radius of the galsim source is 5.5 pixels
	all_paras.max_distance = max_radius; 

#ifdef IMG_CHECK_LABEL	
	double *big_img = new double[stamp_nx*stamp_nx*size*size]();
#endif

	int img_len;
	img_len = size*size;
	double *point = new double[2 * num_p]();
	double *gal = new double[img_len]();
	double *pgal = new double[img_len]();
	double *psf = new double[img_len]();
	double *ppsf = new double[img_len]();
	double *noise = new double[img_len]();
	double *pnoise = new double[img_len]();
	double *shear = new double[2 * shear_pairs](); // [...g1,...,..g2,...]

	float *fft_img_in = new float[img_len]{};
	float *fft_img_out = new float[img_len]{};

	// the shear estimators data matrix  
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
		total_flux = new double[total_data_row]{};
		total_data = new double[total_data_row*shear_data_cols]{};
	}
	sub_flux = new double[scatter_count[rank]]{};

	// read shear
	read_text(str_shear_path, shear, 2*shear_pairs);

	// create PSF
	create_psf(psf, psf_scale, size, psf_type);
	for(i=0;i<img_len;i++)
	{
		fft_img_in[i] = psf[i];
	}
	//pow_spec(psf, ppsf, size, size);
	pow_spec(fft_img_in,fft_img_out,size,size);
	for(i=0;i<img_len;i++)
	{
		ppsf[i] = fft_img_out[i];
	}
	get_psf_radius(ppsf, &all_paras, psf_thresh_scale);

	// std::cout<<rank<<std::endl;
	// show_arr(scatter_count,1,numprocs);
	// std::cout<<"---------------------------------------------------------------------------"<<std::endl;
	
	if (0 == rank)
	{	
		std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		std::cout << parent_path << std::endl;
		std::cout<<"Shear num: "<<shear_pairs<<std::endl;
		std::cout << "Total chip: " << total_chips<< ", Stamp size: " << size  << std::endl;
		std::cout << "Total cpus: " << numprocs << std::endl;
		std::cout << "PSF THRESH: " << all_paras.psf_pow_thresh << " PSF HLR: " << all_paras.psf_hlr << std::endl;
		std::cout <<"MAX RADIUS: "<< max_radius << ", SIG_LEVEL: " << sig_level <<"sigma"<< std::endl;
		sprintf(buffer, "%s/psf.hdf5", parent_path);
		sprintf(set_name, "/data");
		write_h5(buffer,set_name,psf,size, size, true);
		//write_fits(buffer,psf, size, size);
		std::cout<<"Chip Num of each thread: ";
		show_arr(scatter_count,1,numprocs);
		std::cout<<"---------------------------------------------------------------------------"<<std::endl<<std::endl;
	}
	for(i=0;i<numprocs;i++)
	{
		if(i == rank)
		{
			std::cout<<rank<<" "<<chip_st<<" "<<chip_ed<<" "<<scatter_count[i]<<" "<<gather_count[i]<<std::endl;
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
		sprintf(log_inform, "PSF scale: %.2f, max radius: %.2f", psf_scale, max_radius);
		write_log(log_path, log_inform);
		sprintf(log_inform, "RANK: %03d, SHEAR %02d: my chips: %d - %d", rank, shear_id, chip_st, chip_ed);
		write_log(log_path, log_inform);

		// rank 0 reads the total flux array and scatters to each thread
		if(rank == 0)
		{
			std::cout<<"---------------------------------------------------------------------------"<<std::endl;
			std::cout<<"g1: "<<g1<<" g2: "<<g2<<std::endl;
			sprintf(para_path, "%s/parameters/para_%d.hdf5", parent_path, shear_id);
			sprintf(set_name, "/flux");
			read_h5(para_path, set_name, total_flux);
		}
		my_Scatterv(total_flux, scatter_count, sub_flux, numprocs, rank);

		// sprintf(para_path, "%s/parameters/para_%d_check_%d.hdf5", parent_path, shear_id,rank);
		// sprintf(set_name, "/data");
		// write_h5(para_path, set_name, sub_flux, sub_data_row,1,true);


		MPI_Barrier(MPI_COMM_WORLD);

		// loop the chips
		for (i = chip_st; i < chip_ed; i++)
		{
			t1 = clock();

			sprintf(log_inform, "RANK: %03d, SHEAR %02d:, chip: %05d, start.", rank,shear_id, i);
			write_log(log_path, log_inform);
			if (0 == rank)
			{
				std::cout << log_inform << std::endl;
			}

			// initialize GSL
			seed = shear_id*500000 + rank*50000 + i*5;
			gsl_initialize(seed);

#ifdef IMG_CHECK_LABEL
			initialize_arr(big_img, stamp_nx*stamp_nx*size*size, 0);
#endif

			row = (i-chip_st)*stamp_num*shear_data_cols;

			// loop the stamps
			for (j = 0; j < stamp_num; j++)
			{	

				flux_i = sub_flux[(i-chip_st)*stamp_num + j] / num_p;
				
				initialize_arr(point, num_p * 2, 0);
				create_points(point, num_p, max_radius);
				

				///////////////// Noise free /////////////////////////
				all_paras.n1 = 0;
				all_paras.n2 = 0;
				all_paras.dn = 0;
				all_paras.du = 0;

				initialize_arr(gal, size*size, 0);
				initialize_arr(pgal, size*size, 0);

				convolve(gal, point, flux_i, size, num_p, 0, psf_scale, g1, g2, psf_type, 1, &all_paras);
				for(i=0;i<img_len;i++)
				{
					fft_img_in[i] = gal[i];
				}
				pow_spec(fft_img_in,fft_img_out,size,size);
				for(i=0;i<img_len;i++)
				{
					pgal[i] = fft_img_out[i];
				}			
				//pow_spec(gal, pgal, size, size);
				shear_est(pgal, ppsf, &all_paras);

				sub_data[row + j * shear_data_cols + 0] = all_paras.n1;
				sub_data[row + j * shear_data_cols + 1] = all_paras.n2;
				sub_data[row + j * shear_data_cols + 2] = all_paras.dn;
				sub_data[row + j * shear_data_cols + 3] = all_paras.du;


				////////////////// noisy image //////////////////////
				all_paras.n1 = 0;
				all_paras.n2 = 0;
				all_paras.dn = 0;
				all_paras.du = 0;

				initialize_arr(noise, size*size, 0);
				initialize_arr(pnoise, size*size, 0);
				initialize_arr(pgal, size*size, 0);

				addnoise(gal, size*size, gal_noise_sig);

#ifdef IMG_CHECK_LABEL
				stack(big_img, gal, j, size, stamp_nx, stamp_nx);
				//std::cout<<flux_i<<" "<<sub_flux[(i-chip_st)*stamp_num + j]<<std::endl;
#endif
				for(i=0;i<img_len;i++)
				{
					fft_img_in[i] = gal[i];
				}
				pow_spec(fft_img_in,fft_img_out,size,size);
				for(i=0;i<img_len;i++)
				{
					pgal[i] = fft_img_out[i];
				}
				//pow_spec(gal, pgal, size, size);

				addnoise(noise, size*size, gal_noise_sig);
				for(i=0;i<img_len;i++)
				{
					fft_img_in[i] = noise[i];
				}
				pow_spec(fft_img_in,fft_img_out,size,size);
				for(i=0;i<img_len;i++)
				{
					pnoise[i] = fft_img_out[i];
				}	
				//pow_spec(noise, pnoise, size, size);

				noise_subtraction(pgal, pnoise, &all_paras, 1, 1);
				shear_est(pgal, ppsf, &all_paras);

				sub_data[row + j * shear_data_cols + 4] = all_paras.n1;
				sub_data[row + j * shear_data_cols + 5] = all_paras.n2;
				sub_data[row + j * shear_data_cols + 6] = all_paras.dn;
				sub_data[row + j * shear_data_cols + 7] = all_paras.du;

			}
			gsl_free();

#ifdef IMG_CHECK_LABEL
			if(i == IMG_CHECK_LABEL)
			{
				sprintf(chip_path, "%s/%d/gal_chip_%05d.hdf5", parent_path, shear_id, i);
				sprintf(set_name,"/data");
				write_h5(chip_path, set_name,big_img, stamp_nx*size, stamp_nx*size,true);
				//write_fits(chip_path, big_img, stamp_nx*size, stamp_nx*size);
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

		te = clock();
		if (rank == 0)
		{
			sprintf(buffer, "rank %d: done in %g \n", rank, (te - ts) / CLOCKS_PER_SEC);
			std::cout << buffer<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// sprintf(result_path, "%s/result/data/data_%d_%d_check.hdf5", parent_path, shear_id,rank);
		// sprintf(set_name, "/data");
		// write_h5(result_path, set_name, sub_data, sub_data_row,shear_data_cols,true);

		sprintf(set_name, "/data");
		//MPI_Gather(sub_data, sub_data_row*shear_data_cols, MPI_DOUBLE, total_data, sub_data_row*shear_data_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		my_Gatherv(sub_data, gather_count, total_data, numprocs, rank);

		if (0 == rank)
		{
			sprintf(result_path, "%s/result/data/data_%d.hdf5", parent_path, shear_id);
			write_h5(result_path, set_name, total_data, total_data_row, shear_data_cols, true);
			std::cout<<"---------------------------------------------------------------------------"<<std::endl;
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (0 == rank)
	{
		std::cout << "FINISH ALL JOBS" << std::endl;
		delete[] total_data;
		delete[] total_flux;
	}
#ifdef IMG_CHECK_LABEL
	delete[] big_img;
#endif
	delete[] point;
	delete[] gal;
	delete[] pgal;
	delete[] psf;
	delete[] ppsf;
	delete[] noise;
	delete[] pnoise;
	delete[] sub_data;
	delete[] shear;
	delete[] sub_flux;
	delete[] scatter_count;
	delete[] gather_count;
	
	MPI_Finalize();
	return 0;
}
