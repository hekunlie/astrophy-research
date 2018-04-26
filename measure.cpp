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
	
	int size = 64, shear_pairs = 14, chip_num, stamp_num = 10000, stamp_nx = 100;
	chip_num = 500 / (numprocs / 14);
	int data_rows = chip_num*stamp_num, data_cols = 17;
	int i, j, k, seed, chip_id, shear_id;
	double thres = 2., psf_scale = 5., psf_noise_sig = 0, gal_noise_sig = 380.86, ts, te, t1, t2;
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
	double*fit_img = new double[size*size]();

	for (i = 0; i < data_rows; i++)
	{
		data[i] = matrix + i*data_cols;
	}
	char chip_path[150], buffer[200], h5_path[150], set_name[50], log_path[150], log_inform[150];

	sprintf(log_path, "/lmc/selection_bias_0/logs/m_%d_log.dat", myid);

	sprintf(chip_path, "/lmc/selection_bias_0/psf.fits");
	//read_img(psf, chip_path);
	create_psf(psf, psf_scale, size, 2);
	pow_spec(psf, ppsf, size, size);
	get_radius(ppsf, &all_paras, thres, 1,psf_noise_sig);
	get_psf_thres(ppsf, &all_paras);

	if (myid == 0)
	{
		cout << all_paras.psf_pow_thres << endl;
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

		sprintf(chip_path, "/lmc/selection_bias_0/%d/gal_chip_%04d.fits", shear_id, i + chip_id);
		read_img(big_img, chip_path);

		for (j = 0; j < stamp_num; j++)
		{
			addnoise(noise, size*size, gal_noise_sig);
			pow_spec(noise, pnoise, size, size);

			segment(big_img, gal, j, size, stamp_nx, stamp_nx);
			get_radius(gal, &all_paras, 2., 2, gal_noise_sig);
			pow_spec(gal, pgal, size, size);

			f_snr(pgal, &all_paras);
			data[i*stamp_num + j][10] = all_paras.gal_fsnr_c;// original fsnr

			//paraboloid_fit(pgal, fit_img, &all_paras, size);
//			f_snr(fit_img, &all_paras);
	//		data[i*stamp_num + j][13] = all_paras.gal_fsnr_c;//the fsnr on the fitting image

			shear_est(pgal, ppsf, pnoise, &all_paras);



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
			//initialize(fit_img);
		}
		initialize(big_img);

		t2 = clock();
		sprintf(log_inform, "%03d 's chip finish in %.2f sec", i, (t2 - t1) / CLOCKS_PER_SEC);
		write_log(log_path, log_inform);
		if (0 == myid)
		{
			cout << log_inform << endl;
		}
	}

	sprintf(h5_path, "/lmc/selection_bias_0/result/data/data_f_%d_%d.hdf5", shear_id, myid / shear_pairs);
	sprintf(set_name, "/data");
	cout << myid << " " << h5_path << endl;
	write_h5(h5_path, set_name, data_rows, data_cols, data[0], NULL);

	te = clock();
	sprintf(log_inform, "write data to file and finish jobs in %.2f sec", (te - ts) / CLOCKS_PER_SEC);
	write_log(log_path, log_inform);
	if (0 == myid)
	{
		cout << log_inform << endl;
	}

	gsl_rng_free();
	delete[] psf;
	delete[] ppsf;
	delete[] big_img;
	delete[] gal;
	delete[] pgal;
	delete[] noise;
	delete[] pnoise;
	delete[] data;
	delete[] matrix;
	delete[] fit_img;
	
	MPI_Finalize();
	return 0;
}