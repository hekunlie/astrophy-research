#include <cmath>
#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include "mpi.h"
#include "FQlib.h"
#include<hdf5.h>
#include<stdio.h>

//#define TRANS_S_STD 0.5
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

	int size = 90, shear_pairs = 14, chip_num, stamp_num=10000, stamp_nx =100;
	chip_num = 500 /(numprocs / 14);
	int data_rows = chip_num*stamp_num, data_cols = 17;
	int i, j, k, seed, chip_id, shear_id;
	double thres = 2.,  psf_scale = 4., psf_noise_sig = 0, gal_noise_sig = 380.86, ts, te, t1, t2;

	all_paras.gal_noise_sig = gal_noise_sig;
	all_paras.psf_noise_sig = psf_noise_sig;
	all_paras.img_size = size;

	shear_id = myid - myid / shear_pairs*shear_pairs;
	chip_id = myid / shear_pairs*chip_num;

	ts = clock();
	seed = myid * 15322 + 43132;
	gsl_rng_initialize(seed);

	double *psf = new double[size*size]();
	double *ppsf = new double[size*size]();
	double *big_img = new double[size*size*stamp_num]();
	double *gal = new double[size*size]();
	double *pgal = new double[size*size]();
	double *noise = new double[size*size]();
	double *pnoise = new double[size*size]();
	double *matrix = new double[data_rows*data_cols]();
	double **data = new double*[data_rows];
	double*temp_g = new double[size*size]();
	double*gal_fit_img = new double[size*size]();
	double *noise_fit_img = new double[size*size]();

	double coeffs[150 * 25];

	for (i = 0; i < data_rows; i++)
	{
		data[i] = matrix + i*data_cols;
	}

	char chip_path[150], buffer[200], h5_path[150], set_name[50], log_path[150], log_inform[150];
	sprintf(buffer, "/home/hkli/work/c/coeffs.hdf5");
	sprintf(set_name, "/data");
	read_h5(buffer, set_name, coeffs, NULL, NULL, NULL, NULL);

	sprintf(log_path, "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/logs/m_%d_log.dat", myid);

	sprintf(chip_path, "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/psf.fits");
	read_img(psf, chip_path);
	//create_psf(psf, psf_scale, size, 2);
	pow_spec(psf, ppsf, size, size);
	get_radius(ppsf, &all_paras, thres, 1, psf_noise_sig);
	get_psf_thres(ppsf, &all_paras);
	if (0 == myid)
	{
		cout << "PSF THRES: " << all_paras.psf_pow_thres << endl;
	}
	for (i = 0; i < chip_num; i++)
	{	
		t1 = clock();
		sprintf(log_inform, "%03d 's chip start...", i);
		write_log(log_path, log_inform);
		if (0 == myid)
		{
			cout << log_inform << endl;
		}

		sprintf(chip_path, "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/%d/gal_chip_%04d.fits", shear_id, i+chip_id);
		read_img(big_img, chip_path);

		for (j = 0; j < stamp_num; j++)
		{
			addnoise(noise, size*size, gal_noise_sig);
			pow_spec(noise, pnoise, size, size);
			smooth(pnoise, noise_fit_img, ppsf, coeffs, &all_paras);

			segment(big_img, gal, j, size, stamp_nx, stamp_nx);
			get_radius(gal, &all_paras, 99999999 * thres, 2, gal_noise_sig);
			pow_spec(gal, pgal, size, size);

			f_snr(pgal, &all_paras);
			data[i*stamp_num + j][10] = all_paras.gal_fsnr_c;// original fsnr

			smooth(pgal, gal_fit_img, ppsf, coeffs, &all_paras);
			f_snr(gal_fit_img, &all_paras);
			data[i*stamp_num + j][13] = all_paras.gal_fsnr_c;//the fsnr on the fitting image

			shear_est(gal_fit_img, ppsf, noise_fit_img, &all_paras);

			/* subtact the noise background */
			//for (k = 0; k < size*size; k++)
			//{
			//	temp_g[k] = pgal[k] - pnoise[k];
			//}
			//f_snr(temp_g, &all_paras, size);
			//data[i*stamp_num + j][14] = all_paras.gal_fsnr_c;// fsnr on the image of which the background is subtracted

			//initialize(fit_img, size*size);
			//paraboloid_fit(temp_g, fit_img, &all_paras, size);
			//f_snr(fit_img, &all_paras, size);
			//data[i*stamp_num + j][15] = all_paras.gal_fsnr_c;
			

			data[i*stamp_num + j][0] = 0;
			data[i*stamp_num + j][1] = 0;
			data[i*stamp_num + j][2] = all_paras.n1;
			data[i*stamp_num + j][3] = all_paras.n2;
			data[i*stamp_num + j][4] = all_paras.dn;
			data[i*stamp_num + j][5] = all_paras.du;
			data[i*stamp_num + j][6] = all_paras.dv;
			data[i*stamp_num + j][7] = all_paras.gal_osnr;
			data[i*stamp_num + j][8] = all_paras.gal_flux;
			data[i*stamp_num + j][9] = all_paras.gal_peak;
			
			data[i*stamp_num + j][11] = all_paras.gal_snr;
			data[i*stamp_num + j][12] = all_paras.gal_size;
			
			data[i*stamp_num + j][14] = 0.;
			data[i*stamp_num + j][15] = 0.;
			data[i*stamp_num + j][16] = 0.;

			initialize(noise);
			initialize(pnoise);
			initialize(gal);
			initialize(pgal);
			initialize(temp_g);
			initialize(gal_fit_img);
			initialize(noise_fit_img);
		}
		initialize(big_img);

		t2 = clock();
		sprintf(log_inform, "%03d 's chip finish in %.2f sec", i, (t2-t1)/CLOCKS_PER_SEC);
		write_log(log_path, log_inform);
		if (0 == myid)
		{
			cout << log_inform << endl;
		}
	}

	sprintf(h5_path, "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/result/data/data_f_%d_%d.hdf5", shear_id, myid / shear_pairs);
	sprintf(set_name, "/data");
	write_h5(h5_path, set_name, data_rows, data_cols, data[0], NULL);

	te = clock();
	sprintf(log_inform, "write data to file and finish jobs in %.2f sec", (te-ts)/CLOCKS_PER_SEC);
	write_log(log_path, log_inform);
	if (0 == myid)
	{
		cout << log_inform << endl;
	}
	delete[] psf;
	delete[] ppsf;
	delete[] big_img;
	delete[] gal;
	delete[] pgal;
	delete[] noise;
	delete[] pnoise;
	delete[] data;
	delete[] matrix;
	delete[] gal_fit_img;
	delete[] noise_fit_img;
	delete[] temp_g;
	gsl_rng_free();
	MPI_Finalize();
	return 0;
}