#include <stdlib.h>
#include "mpi.h"
#include "FQlib.h"
#include<cstdio>
#include<hk_iolib.h>

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

	char data_path[100], chip_path[150], snr_h5_path[150], para_path[150], shear_path[150],h5_path[150], log_path[150];
	char buffer[200], log_inform[250], set_1[50], set_2[50], finish_path[150];

	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/simu_test");

	std::ifstream fin;
	std::string str_stampsize = "stamp_size", str_total_num = "total_num", str_noise = "noise_sig", str_shear_num = "shear_num", str_nx = "stamp_col";
	std::string dect_info;
	std::string str_data_path,str_para_path,str_shear_path;


	int i, j, k, ib;
	int sss1, sss2, seed;

	int num_p, max_radius, size, total_chips, chip_num, shear_pairs, data_row, total_data_row;
	int stamp_num, stamp_nx, shear_data_cols;		
	int row, chip_st, chip_ed, shear_id, psf_type, temp_s, detect_label;

	double psf_scale, psf_thresh_scale, sig_level, psf_noise_sig, gal_noise_sig, flux_i, mag_i;
	double g1, g2, ts, te, t1, t2;
	double psf_ellip, psf_ang, psf_norm_factor;

	num_p = 100;
	max_radius=9;
	stamp_num = 10000;
	shear_data_cols = 7;

	psf_type=2;
	psf_scale = 4;
	psf_thresh_scale = 2.;
	temp_s = rank;
	sig_level = 1.5;
	psf_noise_sig = 0;

	char_to_str(data_path, str_data_path);
	str_para_path = str_data_path + "/parameters/para.ini";
	str_shear_path = str_data_path + "/parameters/shear.dat";
	sprintf(log_path, "%s/logs/%02d.dat", data_path, rank);

	read_para(str_para_path, str_stampsize, size);
	read_para(str_para_path, str_total_num, total_chips);
	read_para(str_para_path, str_noise, gal_noise_sig);
	read_para(str_para_path, str_shear_num, shear_pairs);
	read_para(str_para_path, str_nx, stamp_nx);

	/* the total chips should be divisible by the numprocs !!!! */
	chip_num = total_chips / numprocs;
	chip_st = chip_num * rank;
	chip_ed = chip_num * (rank + 1);

	total_data_row = total_chips * stamp_num;
	data_row = chip_num * stamp_num;

	all_paras.stamp_size = size;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.max_distance = max_radius; // because the max half light radius of the galsim source is 5.5 pixels

	initialize_para(&all_paras);
	
	double *big_img = new double[stamp_nx*stamp_nx*size*size]();
	double *point = new double[2 * num_p]();
	double *gal = new double[size*size]();
	double *pgal = new double[size*size]();
	double *psf = new double[size*size]();
	double *ppsf = new double[size*size]();
	double *noise = new double[size*size]();
	double *pnoise = new double[size*size]();
	double *shear = new double[2 * shear_pairs](); // [...g1,...,..g2,...]
	int *mask = new int[size*size]();

	// the shear estimators data matrix  
	double *data = new double[data_row*shear_data_cols]();
	double *recvbuf;
	double *flux = new double[total_data_row];
	double *mag = new double[total_data_row];
	if (0 == rank)
	{
		recvbuf = new double[total_chips*stamp_num*shear_data_cols];
	}

	read_text(str_shear_path, shear, 2*shear_pairs);

	// circle PSF
	create_psf(psf, psf_scale, size, psf_type);

	// elliptical PSF, e = 0.05, position angle = Pi/4
	//psf_ellip = 0.05;
	//psf_ang = Pi / 4;
	//psf_norm_factor = 19.0643; // by numercal integrate
	//create_psf(psf, psf_scale, size, psf_ellip, psf_ang, psf_norm_factor, psf_type);

	pow_spec(psf, ppsf, size, size);
	get_psf_radius(ppsf, &all_paras, psf_thresh_scale);
	if (0 == rank)
	{
		std::cout << data_path << std::endl;
		std::cout << "PSF THRESH: " << all_paras.psf_pow_thresh << "PSF HLR: " << all_paras.psf_hlr << std::endl;
		std::cout <<"MAX RADIUS: "<< max_radius << ", SIG_LEVEL: " << sig_level <<"sigma"<< std::endl;
		std::cout << "Total chip: " << total_chips << ", Total cpus: " << numprocs << ", Stamp size: " << size << std::endl;
		sprintf(buffer, "!%s/psf.fits", data_path);
		write_fits(buffer,psf, size, size);
	}

	for (shear_id = 0; shear_id < shear_pairs; shear_id++)
	{
		ts = clock();

		sprintf(log_inform, "size: %d, total chips: %d (%d cpus),  point num: %d , noise sigma: %.2f ", size, total_chips, numprocs, num_p, gal_noise_sig);
		write_log(log_path, log_inform);
		sprintf(log_inform, "PSF scale: %.2f, max radius: %.2f", psf_scale, max_radius);
		write_log(log_path, log_inform);
		sprintf(log_inform, "RANK: %03d, SHEAR %02d: my chips: %d - %d", rank, shear_id, chip_st, chip_ed);
		write_log(log_path, log_inform);

		sprintf(para_path, "%s/parameters/para_%d.hdf5", data_path, shear_id);
		sprintf(set_1, "/flux");
		sprintf(set_2, "/mag");
		read_h5(para_path, set_1, flux);
		read_h5(para_path, set_2, mag);

		g1 = shear[shear_id];
		g2 = shear[shear_id + shear_pairs];

		for (i = chip_st; i < chip_ed; i++)
		{
			t1 = clock();

			seed = rank * i + shear_id + 15115+i+temp_s;
			temp_s++;
			gsl_initialize(seed + i);

			initialize_arr(big_img, stamp_nx*stamp_nx*size*size, 0);

			sprintf(log_inform, "RANK: %03d, SHEAR %02d:, chip: %04d, start.", rank,shear_id, i);
			write_log(log_path, log_inform);
			if (0 == rank)
			{
				std::cout << log_inform << std::endl;
			}

			row = (i - chip_st) *stamp_num*shear_data_cols;

			for (j = 0; j < stamp_num; j++)
			{
				initialize_arr(gal, size*size, 0);
				initialize_arr(pgal, size*size, 0);
				initialize_arr(point, num_p * 2, 0);
				initialize_arr(noise, size*size, 0);
				initialize_arr(pnoise, size*size, 0);
				initialize_para(&all_paras);

				create_points(point, num_p, max_radius);
				flux_i = flux[i*stamp_num + j] / num_p;
				// for measuring the intrinsic ellipticity
				// circle PSF
				//convolve(gal, point, flux_i, size, num_p, 0, psf_scale, 0, 0, psf_type, 0, &all_paras);
				// elliptical PSF
				//convolve(gal, point, flux_i, size, num_p, 0, psf_scale, 0, 0, psf_type, 0, psf_ellip, psf_ang, psf_norm_factor, &all_paras);

				// circle PSF
				//initialize_arr(gal, size*size, 0);
				convolve(gal, point, flux_i, size, num_p, 0, psf_scale, g1, g2, psf_type, 1, &all_paras);					
				// elliptical PSF
				//convolve(gal, point, flux_i, size, num_p, 0, psf_scale, g1, g2, psf_type, 1, psf_ellip, psf_ang, psf_norm_factor, &all_paras);

				addnoise(gal, size*size, gal_noise_sig);

				stack(big_img, gal, j, size, stamp_nx, stamp_nx);

				// galaxy_finder(gal, mask, &all_paras, false, detect_label, dect_info);

				pow_spec(gal, pgal, size, size);

				addnoise(noise, size*size, gal_noise_sig);
				pow_spec(noise, pnoise, size, size);

				noise_subtraction(pgal, pnoise, &all_paras, 1, 1);
				shear_est(pgal, ppsf, &all_paras);

				data[row + j * shear_data_cols + 0] = all_paras.gal_e1;
				data[row + j * shear_data_cols + 1] = all_paras.gal_e2;
				data[row + j * shear_data_cols + 2] = all_paras.n1;
				data[row + j * shear_data_cols + 3] = all_paras.n2;
				data[row + j * shear_data_cols + 4] = all_paras.dn;
				data[row + j * shear_data_cols + 5] = all_paras.du;
				data[row + j * shear_data_cols + 6] = all_paras.dv;
			}
			gsl_free();

			sprintf(chip_path, "!%s%d/gal_chip_%04d.fits", data_path, shear_id, i);
			write_fits(chip_path, big_img, stamp_nx*size, stamp_nx*size);

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
		MPI_Gather(data, data_row*shear_data_cols, MPI_DOUBLE, recvbuf, data_row*shear_data_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (0 == rank)
		{
			sprintf(h5_path, "%sresult/data/data_%d.hdf5", data_path, shear_id);
			write_h5(h5_path, set_1, recvbuf, total_data_row, shear_data_cols, TRUE);
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}
	sprintf(log_path, "/home/hkli/work/test/job/debug/finish_%d.dat", rank);
	write_log(log_path, log_path);

	if (0 == rank)
	{
		std::cout << data_path << std::endl;
		delete[] recvbuf;
	}

	delete[] big_img;
	delete[] point;
	delete[] gal;
	delete[] pgal;
	delete[] mask;
	delete[] psf;
	delete[] ppsf;
	delete[] noise;
	delete[] pnoise;
	delete[] data;
	delete[] shear;
	delete[] flux;
	delete[] mag;		
	MPI_Finalize();
	return 0;
}
