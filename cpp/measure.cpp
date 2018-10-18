#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include "mpi.h"
#include "FQlib.h"
#include<hdf5.h>
#include<stdio.h>
#include<string>

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


	int size = 60, shear_pairs = 14, chip_num, stamp_num=10000, stamp_nx =100;
	chip_num = 500 / (numprocs / 14);
	int data_rows = chip_num*stamp_num, shear_esti_data_cols = 7, snr_para_data_cols=7;
	int i, j, k, tag, seed, chip_id, shear_id;
	double thres = 2.,  psf_scale = 4., psf_noise_sig = 0, gal_noise_sig = 380.86, ts, te, t1, t2;
	int cmd = 1;

	all_paras.gal_noise_sig = gal_noise_sig;
	all_paras.psf_noise_sig = psf_noise_sig;
	all_paras.stamp_size = size;
	all_paras.max_source = 30;
	all_paras.area_thres = 6;
	all_paras.detect_thres = gal_noise_sig*2;
	all_paras.img_x = size;
	all_paras.img_y = size;
	all_paras.max_distance = 8; /* because the max half light radius of the galsim source is 5.5 pixels */
	
	shear_id = myid - myid / shear_pairs*shear_pairs;
	chip_id = myid / shear_pairs*chip_num;

	ts = clock();
	seed = myid * 15322 + 4332;
	gsl_rng_initialize(seed);

	double *psf = new double[size*size]();
	double *ppsf = new double[size*size]();
	double *big_img = new double[size*size*stamp_num]();
	double *gal = new double[size*size]();
	double *pgal = new double[size*size]();
	double *noise = new double[size*size]();
	double *pnoise = new double[size*size]();
	double *mag = new double[data_rows]();// the data_row could arise bug if runs on more than 14 threads,  data_row = 500/(numprocs/14)*stamp_num
	//double*psf_fit_img = new double[size*size]();
	//double*gal_fit_img = new double[size*size]();
	//double *noise_fit_img = new double[size*size]();
	
	double coeffs[150 * 25];
	double *matrix = new double[data_rows*shear_esti_data_cols]();
	double **data = new double*[data_rows];
	for (i = 0; i < data_rows; i++)
	{
		data[i] = matrix + i*shear_esti_data_cols;
	}

	/* the snr parameters data matrix data_snr[i][j] */
	double *data_s = new double[data_rows*snr_para_data_cols]();
	double **data_snr = new double*[data_rows];
	for (i = 0; i < data_rows; i++)
	{
		data_snr[i] = data_s + i*snr_para_data_cols;
	}

	char data_path[100],chip_path[150], snr_h5_path[150], para_path[150], buffer[200], h5_path[150], set_name[50], log_path[150], log_inform[150];
	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/selection_bias_64_bright/");

	sprintf(buffer, "/home/hkli/work/c/coeffs.hdf5");
	sprintf(set_name, "/data");
	read_h5(buffer, set_name, coeffs, NULL, NULL, NULL, NULL);

	sprintf(para_path, "%sparameters/para_%d.hdf5", data_path, shear_id);
	sprintf(set_name, "/mag");

	read_h5(para_path, set_name, mag, NULL, NULL, NULL, NULL);

	sprintf(log_path, "%slogs/m2_%d_log.dat", data_path, myid);

	//sprintf(chip_path, "%spsf.fits", data_path);
	//read_img(psf, chip_path);
	create_psf(psf, psf_scale, size, 2);

	pow_spec(psf, ppsf, size, size);
	get_radius(ppsf, &all_paras, thres, 1, psf_noise_sig);
	
	get_psf_thres(ppsf, &all_paras);

	if (0 == myid)
	{
		cout << "PSF THRES: " << all_paras.psf_pow_thres << endl<< all_paras.psf_hlr << endl;
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

		sprintf(chip_path, "%s%d/gal_chip_%04d.fits", data_path,shear_id, i+chip_id);
		initialize_arr(big_img, stamp_nx*stamp_nx*size*size);
		read_img(big_img, chip_path);		

		for (j = 0; j < stamp_num; j++)
		{	
			initialize_arr(noise, size*size);
			initialize_arr(pnoise, size*size);
			initialize_arr(gal, size*size);
			initialize_arr(pgal, size*size);
			
			segment(big_img, gal, j, size, stamp_nx, stamp_nx);
			pow_spec(gal, pgal, size, size);
			//smooth(pgal, gal_fit_img, ppsf, coeffs, &all_paras);
			galaxy_finder(gal, &all_paras, false);

			f_snr(pgal, &all_paras, 2);

			if (cmd == 0)
			{
				addnoise(noise, size*size, gal_noise_sig);
				pow_spec(noise, pnoise, size, size);
				//smooth(pnoise, noise_fit_img, psf_fit_img, coeffs, &all_paras);

				shear_est(pgal, ppsf, pnoise, &all_paras);

				data[i*stamp_num + j][0] = 0;
				data[i*stamp_num + j][1] = 0;
				data[i*stamp_num + j][2] = all_paras.n1;
				data[i*stamp_num + j][3] = all_paras.n2;
				data[i*stamp_num + j][4] = all_paras.dn;
				data[i*stamp_num + j][5] = all_paras.du;
				data[i*stamp_num + j][6] = all_paras.dv;
			}

			data_snr[i*stamp_num + j][0] = all_paras.gal_osnr;
			data_snr[i*stamp_num + j][1] = all_paras.gal_flux;
			data_snr[i*stamp_num + j][2] = all_paras.gal_flux_alt;
			data_snr[i*stamp_num + j][3] = all_paras.gal_flux2;
			data_snr[i*stamp_num + j][4] = all_paras.gal_snr;
			data_snr[i*stamp_num + j][5] = all_paras.gal_size;
			data_snr[i*stamp_num + j][6] = mag[i*stamp_num + j];

			//data[i*stamp_num + j][7] = all_paras.gal_flux;
			//data[i*stamp_num + j][8] = all_paras.gal_hflux;
			//data[i*stamp_num + j][9] = all_paras.gal_peak;
			//data[i*stamp_num + j][10] = all_paras.gal_size;
			//data[i*stamp_num + j][11] = all_paras.gal_hsize;
			//data[i*stamp_num + j][12] = all_paras.gal_snr;
			//data[i*stamp_num + j][13] = all_paras.gal_flux2;// original fsnr
			//data[i*stamp_num + j][14] = all_paras.gal_flux_alt;// fitted fsnr
			//data[i*stamp_num + j][15] = 0.;
			//data[i*stamp_num + j][16] = 0.;

		}
		
		t2 = clock();
		sprintf(log_inform, "%03d 's chip finish in %.2f sec", i, (t2-t1)/CLOCKS_PER_SEC);
		write_log(log_path, log_inform);
		if (0 == myid)
		{
			cout << log_inform << endl;
		}
	}

	sprintf(set_name, "/data");
	if (cmd == 0)
	{
		sprintf(h5_path, "%sresult/data/data_%d_%d.hdf5", data_path, shear_id, myid / shear_pairs);		
		write_h5(h5_path, set_name, data_rows, shear_esti_data_cols, data[0], NULL);
	}

	sprintf(snr_h5_path, "%sresult/data/data_2sig/data_%d_%d.hdf5", data_path, shear_id, myid / shear_pairs);
	write_h5(snr_h5_path, set_name, data_rows, snr_para_data_cols, data_snr[0], NULL);

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
	delete[] data_s;
	delete[] data_snr;
	delete[] mag;
	//delete[] gal_fit_img;
	//delete[] noise_fit_img;
	//delete[] psf_fit_img;
	gsl_rng_free();
	MPI_Finalize();
	return 0;
}