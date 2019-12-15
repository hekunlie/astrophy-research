#include <cmath>
#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include <hk_mpi.h>
#include "FQlib.h"
#include<cstdio>
#include<string>
#include<hk_iolib.h>

int main(int argc, char*argv[])
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	fq_paras all_paras;

	std::ifstream fin;
	std::string s, str_stampsize = "stamp_size", str_total_num = "total_num", str_noise = "noise_sig", str_shear_num = "shear_num", str_nx = "stamp_col";
	std::string dect_info;

	char data_path[100], chip_path[150], snr_h5_path[150], para_path[150], shear_path[150],h5_path[150], log_path[150];
	char buffer[200], log_inform[250], set_1[50], set_2[50], finish_path[150];
	
	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/simu_test_1");
	//sprintf(data_path, "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/selection_bias/simu/simu_test_1");

	std::string str_data_path;
	char_to_str(data_path, str_data_path);
	sprintf(shear_path,"%s/parameters/shear.hdf5",data_path);	
	std::string str_para_path = str_data_path + "/parameters/para.ini";
	sprintf(log_path, "%s/logs/m_%02d.dat", data_path, rank);
	
	int i, j, k, ib;
	int sss1, sss2;
	int seed, seed_step;

	int total_chips, chip_num, shear_pairs, data_row, total_data_row;
	int stamp_num, stamp_nx, size, num_p;		
	int row, row_s, shear_esti_data_cols, snr_para_data_cols;
	int chip_id_s, chip_id_e, shear_id, psf_type, temp_s, detect_label;
	double max_radius, psf_scale, psf_thresh_scale, sig_level, psf_noise_sig, gal_noise_sig, psf_peak, flux_i, mag_i;

	double g1, g2, ts, te, t1, t2;
	double psf_ellip, psf_ang, psf_norm_factor;

	num_p = 50;
	stamp_num = 10000;
	shear_esti_data_cols = 7;
	snr_para_data_cols = 10;

	max_radius=8;
	psf_scale=4.;
	psf_type = 2;
	psf_thresh_scale = 2.;
	sig_level = 1.5;
	psf_noise_sig = 0;
	psf_peak = 0;
	// if (0 == rank)
	// {
	// 	for (i = 0; i < numprocs; i++)
	// 	{
	// 		sprintf(finish_path, "/home/hkli/work/test/job/debug/finish_%d.dat", i);
	// 		if (remove(finish_path))
	// 		{
	// 			std::cout << "REMOVE: " << finish_path << std::endl;
	// 		}
	// 		else
	// 		{
	// 			std::cout << "FAILURE REMOVE: " << finish_path << std::endl;
	// 		}
	// 	}
	// }
	// MPI_Barrier(MPI_COMM_WORLD);

	read_para(str_para_path, str_stampsize, size);
	read_para(str_para_path, str_total_num, total_chips);
	read_para(str_para_path, str_noise, gal_noise_sig);
	read_para(str_para_path, str_shear_num, shear_pairs);
	read_para(str_para_path, str_nx, stamp_nx);

	chip_num = total_chips / numprocs;
	total_data_row = total_chips * stamp_num;
	data_row = chip_num * stamp_num;

	chip_id_s = chip_num * rank;
	chip_id_e = chip_num * (rank + 1);

	all_paras.gal_noise_sig = gal_noise_sig;
	all_paras.psf_noise_sig = psf_noise_sig;
	all_paras.stamp_size = size;
	all_paras.max_source = 30;
	all_paras.area_thresh = 5;
	all_paras.detect_thresh = gal_noise_sig * sig_level;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.max_distance = max_radius; // because the max half light radius of the galsim source is 5.5 pixels

	initialize_para(&all_paras);
	
	double *big_img = new double[stamp_nx*stamp_nx*size*size]{};
	float *big_img_buffer = new float[stamp_nx*stamp_nx*size*size]{};
	double *point = new double[2 * num_p]{};
	double *gal = new double[size*size]{};
	double *pgal = new double[size*size]{};
	double *psf = new double[size*size]{};
	double *ppsf = new double[size*size]{};
	double *noise = new double[size*size]{};
	double *pnoise = new double[size*size]{};
	double *shear = new double[2 * shear_pairs]{}; // [...g1,...,..g2,...]
	int *mask = new int[size*size]{};

	// the shear estimators data matrix  
	double *data = new double[data_row*shear_esti_data_cols]{};
	// the snr parameters data matrix 
	double *data_s = new double[data_row*snr_para_data_cols]{};
	double *recvbuf, *recvbuf_s;
	double *flux = new double[total_data_row];
	double *mag = new double[total_data_row];
	if (0 == rank)
	{
		recvbuf = new double[total_chips*stamp_num*shear_esti_data_cols];
		recvbuf_s = new double[total_chips*stamp_num*snr_para_data_cols];
	}

	double *g1_t = new double[shear_pairs]{};
	double *g2_t = new double[shear_pairs]{};
	sprintf(set_1,"/g1");
	read_h5(shear_path, set_1, g1_t);
	sprintf(set_1,"/g2");
	read_h5(shear_path, set_1, g2_t);

	// circle PSF
	//create_psf(psf, psf_scale, size, psf_type);

	// elliptical PSF, e = 0.05, position angle = Pi/4
	psf_ellip = 0.05;
	psf_ang = Pi / 4;
	psf_norm_factor = 19.0643; // by numercal integrate
	create_psf(psf, psf_scale, size, psf_ellip, psf_ang, psf_norm_factor, psf_type);

	pow_spec(psf, ppsf, size, size);
	get_psf_radius(ppsf, &all_paras, psf_thresh_scale);
	if (0 == rank)
	{
		std::cout << data_path << std::endl;
		std::cout << "PSF THRES: " << all_paras.psf_pow_thresh << "PSF HLR: " << all_paras.psf_hlr << std::endl;
		std::cout <<"MAX RADIUS: "<< max_radius << ", SIG_LEVEL: " << sig_level <<"sigma"<< std::endl;
		std::cout << "Total chip: " << total_chips << ", Total cpus: " << numprocs << ", Stamp size: " << size << std::endl;
		sprintf(buffer, "!%spsf.fits", data_path);
		write_fits(buffer,psf, size, size);
	}

	seed_step = 2;
	sss1 = 4*seed_step*shear_pairs*total_chips;
	seed = sss1*rank + 1 + 5000;

	for (shear_id = 0; shear_id < shear_pairs; shear_id++)
	{
		ts = clock();

		sprintf(log_inform, "size: %d, total chips: %d (%d cpus),  point num: %d , noise sigma: %.2f", size, total_chips, numprocs, num_p, gal_noise_sig);
		write_log(log_path, log_inform);
		sprintf(log_inform, "PSF scale: %.2f, max radius: %.2f", psf_scale, max_radius);
		write_log(log_path, log_inform);
		sprintf(log_inform, "RANK: %03d, SHEAR %02d: my chips: %d - %d", rank, shear_id, chip_id_s, chip_id_e);
		write_log(log_path, log_inform);

		sprintf(para_path, "%s/parameters/para_%d.hdf5", data_path, shear_id);
		sprintf(set_1, "/flux");
		sprintf(set_2, "/mag");
		read_h5(para_path, set_1, flux);
		read_h5(para_path, set_2, mag);

		g1 = g1_t[shear_id];
		g2 = g2_t[shear_id];

		for (i = chip_id_s; i < chip_id_e; i++)
		{
			t1 = clock();
			
			gsl_initialize(seed, 0);
			gsl_initialize(seed, 1);

			initialize_arr(big_img, stamp_nx*stamp_nx*size*size, 0);

			sprintf(log_inform, "RANK: %03d, SHEAR %02d: (%.4f, %.4f), chip: %04d, seed: %d. start.", rank, shear_id, g1, g2, i, seed);
			write_log(log_path, log_inform);
			if (0 == rank)
			{
				std::cout << log_inform << std::endl;
			}

			seed += seed_step;
			row = (i - chip_id_s) *stamp_num*shear_esti_data_cols;
			row_s = (i - chip_id_s) *stamp_num*snr_para_data_cols;

			for (j = 0; j < stamp_num; j++)
			{
				initialize_arr(gal, size*size, 0);
				initialize_arr(pgal, size*size, 0);
				initialize_arr(point, num_p * 2, 0);
				initialize_arr(noise, size*size, 0);
				initialize_arr(pnoise, size*size, 0);
				initialize_para(&all_paras);

				create_points(point, num_p, max_radius, rng0);
				flux_i = flux[i*stamp_num + j] / num_p;
				// for measuring the intrinsic ellipticity
				// circle PSF
				//convolve(gal, point, flux_i, size, num_p, 0, psf_scale, 0, 0, psf_type, 0, &all_paras);
				// elliptical PSF
				//convolve(gal, point, flux_i, size, num_p, 0, psf_scale, 0, 0, psf_type, 0, psf_ellip, psf_ang, psf_norm_factor, &all_paras);

				initialize_arr(gal, size*size, 0);
				// circle PSF
				//convolve(gal, point, flux_i, size, num_p, 0, psf_scale, g1, g2, psf_type, 1, &all_paras);					
				// elliptical PSF
				convolve(gal, point, flux_i, size, num_p, 0, psf_scale, g1, g2, psf_type, 1, psf_ellip, psf_ang, psf_norm_factor, &all_paras);

				addnoise(gal, size*size, gal_noise_sig, rng1);

				stack(big_img, gal, j, size, stamp_nx, stamp_nx);

				galaxy_finder(gal, mask, &all_paras, false, detect_label, dect_info);

				pow_spec(gal, pgal, size, size);

				snr_est(pgal, &all_paras, 2);

				addnoise(noise, size*size, gal_noise_sig, rng1);
				pow_spec(noise, pnoise, size, size);

				noise_subtraction(pgal, pnoise, &all_paras, 1, 1);
				shear_est(pgal, ppsf, &all_paras);

				data[row + j * shear_esti_data_cols + 0] = 0;
				data[row + j * shear_esti_data_cols + 1] = 0;
				data[row + j * shear_esti_data_cols + 2] = all_paras.n1;
				data[row + j * shear_esti_data_cols + 3] = all_paras.n2;
				data[row + j * shear_esti_data_cols + 4] = all_paras.dn;
				data[row + j * shear_esti_data_cols + 5] = all_paras.du;
				data[row + j * shear_esti_data_cols + 6] = all_paras.dv;

				data_s[row_s + j * snr_para_data_cols + 0] = all_paras.gal_flux2;
				data_s[row_s + j * snr_para_data_cols + 1] = all_paras.gal_flux_alt;
				data_s[row_s + j * snr_para_data_cols + 2] = all_paras.gal_flux;
				data_s[row_s + j * snr_para_data_cols + 3] = all_paras.gal_osnr;

				data_s[row_s + j * snr_para_data_cols + 4] = all_paras.gal_flux2_ext[0];//P(k=0)
				data_s[row_s + j * snr_para_data_cols + 5] = all_paras.gal_flux2_ext[1];//P(k=0)_fit
				data_s[row_s + j * snr_para_data_cols + 6] = all_paras.gal_flux2_ext[2];//MAX(P(k=0), P(k=0)_fit)
				data_s[row_s + j * snr_para_data_cols + 7] = all_paras.gal_flux2_ext[3];
				data_s[row_s + j * snr_para_data_cols + 8] = -mag[i*stamp_num + j];
				data_s[row_s + j * snr_para_data_cols + 9] = detect_label;

			}
			gsl_free(0);
			gsl_free(1);

			// float array for saving disk volume
			for (ib = 0; ib < stamp_nx*stamp_nx*size*size; ib++)
			{
				big_img_buffer[ib] = big_img[ib];
			}
			sprintf(chip_path, "!%s/%d/gal_chip_%04d.fits", data_path, shear_id, i);
			write_fits(chip_path, big_img_buffer, stamp_nx*size, stamp_nx*size);

			t2 = clock();
			sprintf(log_inform, "RANK: %03d, SHEAR %02d: chip: %04d, done in %.2f s.", rank, shear_id, i, (t2 - t1) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			if (rank == 0)
			{
				std::cout << log_inform << std::endl;
			}
		}	

		te = clock();
		if (rank == 0)
		{
			sprintf(buffer, "rank %d:  done in %g \n", rank, (te - ts) / CLOCKS_PER_SEC);
			std::cout << buffer<<std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		sprintf(set_1, "/data");
		MPI_Gather(data, data_row*shear_esti_data_cols, MPI_DOUBLE, recvbuf, data_row*shear_esti_data_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (0 == rank)
		{
			sprintf(h5_path, "%s/result/data/data_%d.hdf5", data_path, shear_id);
			write_h5(h5_path, set_1, recvbuf, total_data_row, shear_esti_data_cols, TRUE);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		// the SNR.. parameters
		MPI_Gather(data_s, data_row*snr_para_data_cols, MPI_DOUBLE, recvbuf_s, data_row*snr_para_data_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (0 == rank)
		{
			sprintf(snr_h5_path, "%s/result/data/data_%.1fsig/data_%d.hdf5", data_path, sig_level, shear_id);
			write_h5(snr_h5_path, set_1, recvbuf_s, total_data_row, snr_para_data_cols, TRUE);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	// sprintf(log_path, "/home/hkli/work/test/job/debug/finish_%d.dat", rank);
	// write_log(log_path, log_path);

	if (0 == rank)
	{
		std::cout << data_path << std::endl;
		delete[] recvbuf;
		delete[] recvbuf_s;
	}

	delete[] big_img;
	delete[] big_img_buffer;
	delete[] point;
	delete[] gal;
	delete[] pgal;
	delete[] mask;
	delete[] psf;
	delete[] ppsf;
	delete[] noise;
	delete[] pnoise;
	delete[] data;
	delete[] data_s;
	delete[] shear;
	delete[] flux;
	delete[] mag;		
	MPI_Finalize();
	return 0;
}
