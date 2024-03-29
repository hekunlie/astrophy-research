#include <cmath>
#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include "mpi.h"
#include "FQlib.h"
#include<hdf5.h>
#include<cstdio>
#include<string>

int main(int argc, char*argv[])
{
	int myid, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	para all_paras;

	std::ifstream fin;
	std::string dect_info;

	char data_path[100], chip_path[150], snr_h5_path[150], para_path[150], shear_path[150], h5_path[150], log_path[150];
	char buffer[200], log_inform[250], set_1[50], set_2[50], finish_path[150];

	sprintf(data_path, "/mnt/perc/hklee/multi_shear/data/");
	sprintf(log_path, "%slogs/m_%02d.dat", data_path, myid);

	int num_p = 100, size, total_chips, chip_num, shear_pairs, data_row, total_data_row;
	int stamp_num = 10000, stamp_nx, shear_esti_data_cols = 7, snr_para_data_cols = 10;
	int row, row_s, seed, chip_id_s, chip_id_e, shear_id, psf_type = 2, temp_s = myid, detect_label;
	double max_radius = 9, psf_scale = 4., psf_thres_scale = 2., sig_level = 1.5, psf_noise_sig = 0, gal_noise_sig, psf_peak = 0, flux_i, mag_i;
	int i, j, k, ib;
	int sss1, sss2;
	double g1, g2, ts, te, t1, t2;
	double psf_ellip, psf_ang, psf_norm_factor;

	size = 50;
	total_chips = 20;
	gal_noise_sig = 10;
	shear_pairs = 10;
	stamp_nx = 100;

	chip_num = total_chips / numprocs;
	total_data_row = total_chips * stamp_num;
	data_row = chip_num * stamp_num;

	chip_id_s = chip_num * myid;
	chip_id_e = chip_num * (myid + 1);

	all_paras.gal_noise_sig = gal_noise_sig;
	all_paras.psf_noise_sig = psf_noise_sig;
	all_paras.stamp_size = size;
	all_paras.max_source = 30;
	all_paras.area_thres = 5;
	all_paras.detect_thres = gal_noise_sig * sig_level;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.max_distance = max_radius; // because the max half light radius of the galsim source is 5.5 pixels

	initialize_para(&all_paras);

	double *big_img = new double[stamp_nx*stamp_nx*size*size]();
	float *big_img_buffer = new float[stamp_nx*stamp_nx*size*size]();
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
	double *data = new double[data_row*shear_esti_data_cols]();
	// the snr parameters data matrix 
	double *data_s = new double[data_row*snr_para_data_cols]();
	double *recvbuf, *recvbuf_s;
	double *flux = new double[total_data_row];
	double *mag = new double[total_data_row];
	if (0 == myid)
	{
		recvbuf = new double[total_chips*stamp_num*shear_esti_data_cols];
		recvbuf_s = new double[total_chips*stamp_num*snr_para_data_cols];
	}

	for (i = 0; i < shear_pairs; i++)
	{
		shear[i] = -0.05 + i * 0.01;
		shear[i + shear_pairs] = 0.05 - i * 0.01;
	}

	// circle PSF
	create_psf(psf, psf_scale, size, psf_type);

	pow_spec(psf, ppsf, size, size);
	get_psf_radius(ppsf, &all_paras, psf_thres_scale);
	if (0 == myid)
	{
		std::cout << data_path << std::endl;
		std::cout << "PSF THRES: " << all_paras.psf_pow_thres << "PSF HLR: " << all_paras.psf_hlr << std::endl;
		std::cout << "MAX RADIUS: " << max_radius << ", SIG_LEVEL: " << sig_level << "sigma" << std::endl;
		std::cout << "Total chip: " << total_chips << ", Total cpus: " << numprocs << ", Stamp size: " << size << std::endl;
		sprintf(buffer, "!%spsf.fits", data_path);
		write_fits(buffer, psf, size, size);
	}

	for (shear_id = 0; shear_id < shear_pairs; shear_id++)
	{
		ts = clock();

		sprintf(para_path, "%sparameters/para_%d.hdf5", data_path, shear_id);
		sprintf(set_1, "/flux");
		sprintf(set_2, "/mag");
		read_h5(para_path, set_1, flux);
		read_h5(para_path, set_2, mag);

		g1 = shear[shear_id];
		g2 = shear[shear_id + shear_pairs];

		for (i = chip_id_s; i < chip_id_e; i++)
		{
			t1 = clock();

			seed = myid * i + shear_id + 115 + i + temp_s;
			temp_s++;
			gsl_initialize(seed + i);

			initialize_arr(big_img, stamp_nx*stamp_nx*size*size, 0);

			sprintf(log_inform, "RANK: %03d, SHEAR %02d:, chip: %04d, start.", myid, shear_id, i);
			write_log(log_path, log_inform);
			if (0 == myid)
			{
				std::cout << log_inform << std::endl;
			}

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

				create_points(point, num_p, max_radius);
				flux_i = flux[i*stamp_num + j] / num_p;

				initialize_arr(gal, size*size, 0);
				convolve(gal, point, flux_i, size, num_p, 0, psf_scale, g1, g2, psf_type, 1, &all_paras);

				addnoise(gal, size*size, gal_noise_sig);

				stack(big_img, gal, j, size, stamp_nx, stamp_nx);

				galaxy_finder(gal, mask, &all_paras, false, detect_label, dect_info);

				pow_spec(gal, pgal, size, size);

				snr_est(pgal, &all_paras, 2);

				addnoise(noise, size*size, gal_noise_sig);
				pow_spec(noise, pnoise, size, size);

				noise_subtraction(pgal, pnoise, &all_paras, 1, 1);
				shear_est(pgal, ppsf, &all_paras);

				data[row + j * shear_esti_data_cols + 0] = all_paras.gal_e1;
				data[row + j * shear_esti_data_cols + 1] = all_paras.gal_e2;
				data[row + j * shear_esti_data_cols + 2] = all_paras.n1;
				data[row + j * shear_esti_data_cols + 3] = all_paras.n2;
				data[row + j * shear_esti_data_cols + 4] = all_paras.dn;
				data[row + j * shear_esti_data_cols + 5] = all_paras.du;
				data[row + j * shear_esti_data_cols + 6] = all_paras.dv;

			}
			gsl_free();

			// float array for saving disk volume
			for (ib = 0; ib < stamp_nx*stamp_nx*size*size; ib++)
			{
				big_img_buffer[ib] = big_img[ib];
			}
			sprintf(chip_path, "!%s%d/gal_chip_%04d.fits", data_path, shear_id, i);
			write_fits(chip_path, big_img_buffer, stamp_nx*size, stamp_nx*size);

			t2 = clock();
			sprintf(log_inform, "RANK: %03d, SHEAR %02d: chip: %04d, done in %.2f s.", myid, shear_id, i, (t2 - t1) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			if (myid == 0)
			{
				std::cout << log_inform << std::endl;
			}
		}

		te = clock();
		if (myid == 0)
		{
			sprintf(buffer, "myid %d:  done in %g \n", myid, (te - ts) / CLOCKS_PER_SEC);
			std::cout << buffer << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		sprintf(set_1, "/data");
		MPI_Gather(data, data_row*shear_esti_data_cols, MPI_DOUBLE, recvbuf, data_row*shear_esti_data_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (0 == myid)
		{
			sprintf(h5_path, "%sresult/data/data_%d.hdf5", data_path, shear_id);
			write_h5(h5_path, set_1, recvbuf, total_data_row, shear_esti_data_cols, TRUE);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		// the SNR.. parameters
		MPI_Gather(data_s, data_row*snr_para_data_cols, MPI_DOUBLE, recvbuf_s, data_row*snr_para_data_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (0 == myid)
		{
			sprintf(snr_h5_path, "%sresult/data/data_%.1fsig/data_%d.hdf5", data_path, sig_level, shear_id);
			write_h5(snr_h5_path, set_1, recvbuf_s, total_data_row, snr_para_data_cols, TRUE);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	sprintf(log_path, "/home/hkli/work/test/job/debug/finish_%d.dat", myid);
	write_log(log_path, log_path);

	if (0 == myid)
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
