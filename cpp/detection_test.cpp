﻿#include <cmath>
#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include "mpi.h"
#include "FQlib.h"
#include<hdf5.h>
#include<stdio.h>
#include<string>


using namespace std;

int main(int argc, char*argv[])
{
	int myid, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	para all_paras;
	ifstream fin;
	string s;

	/* 14 (g1,g2) points and each pairs contain 500 chips which cotians 10000 gals */
	int total_chip_num = 20, chip_num, stamp_num = 10000, shear_pairs = 1;
	/* remember to change the data_cols when you change the number of estimators recorded */
	int i, j, seed, data_rows, shear_esti_data_cols = 7, snr_para_data_cols = 7, chip_id, shear_id, detect_label;
	int size = 60, num_p = 40, stamp_nx = 100, psf_type = 2;
	double psf_scale = 4., max_radius = 8., st, ed, s1, s2;
	double g1 = 0., g2 = 0.;
	double gal_noise_sig = 60, psf_noise_sig = 0., scale = 2.;
	int total_num = total_chip_num * stamp_num;

	chip_num = 1;
	data_rows = chip_num*stamp_num;
	chip_id = myid;

	all_paras.stamp_size = size;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.detect_thres = gal_noise_sig * 2;
	all_paras.gal_noise_sig = gal_noise_sig;
	all_paras.psf_noise_sig = psf_noise_sig;
	all_paras.max_distance = 8;
	initialize_para(&all_paras);

	double *big_img = new double[stamp_nx*stamp_nx*size*size]();
	double *point = new double[2 * num_p]();
	double *gal = new double[size*size]();
	double *gpow = new double[size*size]();
	double *psf = new double[size*size]();
	double *ppow = new double[size*size]();
	double *noise = new double[size*size]();
	double *pnoise = new double[size*size]();
	char para_inform[500], data_path[60], shear_path[120], para_path[120], log_path[120], log_inform[120], chip_path[120], buffer[200];
	char h5_path[120], snr_h5_path[120];
	char set_name1[20], set_name2[20];
	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/test/");// the total path of all the data

#ifndef PRECISION
	DATA_TYPE *cp = new DATA_TYPE[stamp_nx*stamp_nx*size*size]();
#endif

	/* the shear estimators data matrix data[i][j] */
	double *data_m = new double[data_rows*shear_esti_data_cols]();
	double **data = new double*[data_rows];
	for (i = 0; i < data_rows; i++)
	{
		data[i] = data_m + i*shear_esti_data_cols;
	}

	/* the snr parameters data matrix data_snr[i][j] */
	double *data_s = new double[data_rows*snr_para_data_cols]();
	double **data_snr = new double*[data_rows];
	for (i = 0; i < data_rows; i++)
	{
		data_snr[i] = data_s + i*snr_para_data_cols;
	}

	// initialize gsl
	int sss1, sss2;
	sss1 = 2430;
	sss2 = 130;
	//sss1 = 5811430;
	//sss2 = 7161130;
	seed = myid *sss1 + sss2;
	//seed = myid *380 + 1401;// no bias
	gsl_rng_initialize(seed);


	// read parameters
	double *flux = new double[total_num];
	double *mag = new double[total_num];

	sprintf(para_path, "%spara.hdf5", data_path);
	sprintf(set_name1, "/flux");
	sprintf(set_name2, "/mag");
	read_h5(para_path, set_name1, flux, set_name2, mag, NULL, NULL);


	//PSF
	create_psf(psf, psf_scale, size, psf_type);

	pow_spec(psf, ppow, size, size);

	get_radius(ppow, &all_paras, scale, 1, psf_noise_sig);

	st = clock();

	for (i = 0; i < chip_num; i++)
	{
		s1 = clock();

		sprintf(chip_path, "!%s/gal_chip_%04d.fits", data_path, chip_id + i);
		initialize_arr(big_img, stamp_nx*stamp_nx*size*size);

		for (j = 0; j < stamp_num; j++)
		{
			initialize_arr(gal, size*size);
			initialize_arr(gpow, size*size);
			initialize_arr(point, num_p * 2);
			initialize_arr(noise, size*size);
			initialize_arr(pnoise, size*size);
			initialize_para(&all_paras);

			create_points(point, num_p, max_radius);

			//create_epoints(point, num_p, ellip[i*stamp_num + j]);

			convolve(gal, point, flux[(chip_id + i)*stamp_num + j] / num_p, size, num_p, 0, psf_scale, g1, g2, psf_type);

			addnoise(gal, size*size, gal_noise_sig);

			stack(big_img, gal, j, size, stamp_nx, stamp_nx);

			//get_radius(gal, &all_paras, 9999999999.*thres, 2, gal_noise_sig);
			detect_label = galaxy_finder(gal, &all_paras, false);
			if (detect_label > -1)
			{
				pow_spec(gal, gpow, size, size);
				f_snr(gpow, &all_paras, 1);
			}

			addnoise(noise, size*size, gal_noise_sig);


			data_snr[i*stamp_num + j][0] = all_paras.gal_osnr;
			data_snr[i*stamp_num + j][1] = all_paras.gal_flux;
			data_snr[i*stamp_num + j][2] = all_paras.gal_flux_alt;
			data_snr[i*stamp_num + j][3] = all_paras.gal_flux2;
			data_snr[i*stamp_num + j][4] = all_paras.gal_snr;
			data_snr[i*stamp_num + j][5] = all_paras.gal_size;
			data_snr[i*stamp_num + j][6] = mag[i*stamp_num + j];

		}

#ifdef PRECISION
		write_img(big_img, size*stamp_nx, size*stamp_nx, chip_path);
#else
		copy(big_img, big_img + stamp_nx*stamp_nx*size*size, cp);
		write_img(cp, size*stamp_nx, size*stamp_nx, chip_path);
#endif	

		s2 = clock();
		sprintf(log_inform, "Thread: %d, chip: %d, done in %.2f s.", myid, i, (s2 - s1) / CLOCKS_PER_SEC);
		write_log(log_path, log_inform);
		if (myid == 0)
		{
			cout << log_inform << endl;
		}
	}

	sprintf(set_name1, "/data");

	sprintf(snr_h5_path, "%sdata_%d.hdf5", data_path, myid);
	write_h5(snr_h5_path, set_name1, data_rows, snr_para_data_cols, data_snr[0], NULL);


	ed = clock();
	if (myid == 0)
	{
		sprintf(buffer, "myid %d:  done in %g \n", myid, (ed - st) / CLOCKS_PER_SEC);
		cout << buffer;
	}

	delete[] flux;
	delete[] big_img;
	delete[] point;
	delete[] gal;
	delete[] gpow;
	delete[] psf;
	delete[] ppow;
	delete[] noise;
	delete[] pnoise;
	delete[] data_m;
	delete[] data;
	delete[] data_s;
	delete[] data_snr;
	delete[] mag;
	gsl_rng_free();
	MPI_Finalize();
	return 0;
}