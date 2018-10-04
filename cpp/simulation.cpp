#include <cmath>
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
		int total_chip_num = 500, chip_num, stamp_num = 10000, shear_pairs = 14;
		/* remember to change the data_cols when you change the number of estimators recorded */
		int i, j, seed, data_rows, data_cols=15, chip_id, shear_id, detect_label;
		int size =60, num_p = 40, stamp_nx = 100,  psf_type = 2;
		double psf_scale = 4., max_radius = 9., st, ed, s1, s2;
		double g1=0., g2=0.;
		double gal_noise_sig = 380.86, psf_noise_sig = 0., scale = 2.;
		int total_num = total_chip_num * stamp_num;

		chip_num = total_chip_num / (numprocs/shear_pairs);
		data_rows = chip_num*stamp_num;
		shear_id = myid - myid / shear_pairs*shear_pairs;
		chip_id = myid / shear_pairs*chip_num;

		all_paras.stamp_size = size;
		all_paras.img_x = size;
		all_paras.img_y = size;
		all_paras.detect_thres = gal_noise_sig*1.5;
		all_paras.gal_noise_sig = gal_noise_sig;
		all_paras.psf_noise_sig = psf_noise_sig;
		initialize_para(&all_paras);

		double *big_img = new double[stamp_nx*stamp_nx*size*size]();
		double *point = new double[2 * num_p]();
		double *gal = new double[size*size]();
		double *gpow = new double[size*size]();
		double *psf = new double[size*size]();
		double *ppow = new double[size*size]();
		double *noise = new double[size*size]();
		double *pnoise = new double[size*size]();
		char para_inform[500], log_inform[100], log_path[100], chip_path[200], para_path[200], buffer[200], h5_path[100];
		char set_name1[20], set_name2[20];

		#ifndef PRECISION
		DATA_TYPE *cp = new DATA_TYPE[stamp_nx*stamp_nx*size*size]();
		#endif

		/* the data matrix data[i][j] */
		double *data_m = new double[data_rows*data_cols]();
		double **data = new double*[data_rows];
		for (i = 0; i < data_rows; i++)
		{
			data[i] = data_m + i*data_cols;
		}

		// initialize gsl
		int sss1, sss2;
		sss1 = 430;
		sss2 = 130;
		//sss1 = 5811430;
		//sss2 = 7161130;
		seed = myid *sss1 +sss2;
		//seed = myid *380 + 1401;// no bias
		gsl_rng_initialize(seed);

		
		const char *str;
		// read shear, the input signal
		double *shear = new double[2 * shear_pairs](); // [...g1,...,..g2,...]
		fin.open("/mnt/ddnfs/data_users/hkli/selection_bias_64_dimmer/parameters/shear.dat");
		i = 0;
		while (!fin.eof())
		{
			getline(fin, s);
			str = s.c_str();
			shear[i] = atof(str);
			i++;
		}
		fin.close();
		
		// read parameters
		double *flux = new double[total_num];
		double *mag = new double[total_num];

		sprintf(para_path,"/mnt/ddnfs/data_users/hkli/selection_bias_64_dimmer/parameters/para_%d.hdf5", shear_id);
		sprintf(set_name1, "/flux");
		sprintf(set_name2, "/mag");
		read_h5(para_path, set_name1, flux, set_name2, mag, NULL, NULL);
		
		
		 //PSF
		create_psf(psf, psf_scale, size, psf_type);

		pow_spec(psf, ppow, size, size);

		get_radius(ppow, &all_paras, scale, 1, psf_noise_sig);

		st = clock();

		g1 = shear[shear_id];
		g2 = shear[shear_id + shear_pairs];

		sprintf(log_path, "/mnt/ddnfs/data_users/hkli/selection_bias_64_dimmer/logs/log_%d_%d.dat", shear_id, myid / shear_pairs);

		sprintf(para_inform, "myid: %03d, chip_id: %d, shear_id: %d, size: %d, chip_num: %d * %d\n \
seed: myid*%d + %d, noise sigma: %.2f, point num: %d \n \
PSF scale: %.2f, max radius: %.2f", myid, chip_id, shear_id, size, chip_num, (numprocs / shear_pairs), sss1, sss2, gal_noise_sig, num_p, psf_scale, max_radius);
		write_log(log_path, para_inform);

		for (i = 0; i < chip_num;  i++)
		{
			s1 = clock();

			sprintf(log_inform, "Thread: %d, chip: %d, start.", myid, i);
			write_log(log_path, log_inform);
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

				convolve(gal, point, flux[(chip_id+i)*stamp_num + j]/num_p, size, num_p, 0, psf_scale, g1, g2, psf_type);

				addnoise(gal, size*size,  gal_noise_sig);

				stack(big_img, gal, j, size, stamp_nx, stamp_nx);

				//get_radius(gal, &all_paras, 9999999999.*thres, 2, gal_noise_sig);
				detect_label = galaxy_finder(gal, &all_paras, false);

				pow_spec(gal, gpow, size, size);
				f_snr(gpow, &all_paras, 1);

				addnoise(noise, size*size, gal_noise_sig);

				pow_spec(noise, pnoise, size, size);

				shear_est(gpow, ppow, pnoise, &all_paras);

				data[i*stamp_num + j][0] = g1;
				data[i*stamp_num + j][1] = g2;
				data[i*stamp_num + j][2] = all_paras.n1;
				data[i*stamp_num + j][3] = all_paras.n2;
				data[i*stamp_num + j][4] = all_paras.dn;
				data[i*stamp_num + j][5] = all_paras.du;
				data[i*stamp_num + j][6] = all_paras.dv;
				data[i*stamp_num + j][7] = all_paras.gal_osnr;
				data[i*stamp_num + j][8] = all_paras.gal_flux;
				data[i*stamp_num + j][9] = all_paras.gal_flux_alt;
				data[i*stamp_num + j][10] = all_paras.gal_flux2;
				data[i*stamp_num + j][11] = all_paras.gal_snr;
				data[i*stamp_num + j][12] = all_paras.gal_size;
				data[i*stamp_num + j][13] = mag[i*stamp_num + j];
				data[i*stamp_num + j][14] = 0;


			}
			sprintf(chip_path, "!/mnt/ddnfs/data_users/hkli/selection_bias_64_dimmer/%d/gal_chip_%04d.fits", shear_id, chip_id + i);

			#ifdef PRECISION
			write_img(big_img, size*stamp_nx, size*stamp_nx, chip_path);			
			#else
			copy(big_img, big_img + stamp_nx*stamp_nx*size*size, cp);
			write_img(cp, size*stamp_nx, size*stamp_nx, chip_path);
			#endif

			initialize_arr(big_img, stamp_nx*stamp_nx*size*size);		

			s2 = clock();
			sprintf(log_inform, "Thread: %d, chip: %d, done in %.2f s.", myid, i, (s2 - s1) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			if (myid == 0)
			{
				cout << log_inform<<endl;
			}
		}

		sprintf(h5_path, "/mnt/ddnfs/data_users/hkli/selection_bias_64_dimmer/result/data/data_%d_%d.hdf5", shear_id, myid / shear_pairs);
		sprintf(set_name1, "/data");
		write_h5(h5_path, set_name1, data_rows, data_cols, data[0],NULL);
		ed = clock();
		if (myid == 0)
		{
			sprintf(buffer, "myid %d:  done in %g \n", myid,  (ed - st) / CLOCKS_PER_SEC);
			cout << buffer;
		}
		delete[] shear;
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
		delete[] mag;
		gsl_rng_free();
		MPI_Finalize();
		return 0;
}
