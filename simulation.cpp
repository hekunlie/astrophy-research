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

		/* 10 (g1,g2) points and each pairs contain 100 chips which cotians 10000 gals */
		int chip_num = 600, stamp_num = 10000, shear_pairs = 14;
		/* remember to change the data_cols when you change the number of estimators recorded */
		int i, j, seed, data_rows, data_cols=20;
		int size =60, num_p = 45, stamp_nx = 100,  psf_type = 2;
		double psf_scale = 5., max_radius = 8, st, ed, s1, s2;
		double g1=0., g2=0.;
		double gal_noise_sig = 380.64, psf_noise_sig = 0., thres = 2.;
		data_rows = chip_num*stamp_num;

		all_paras.gal_noise_sig = gal_noise_sig;
		all_paras.psf_noise_sig = psf_noise_sig;

		double *big_img = new double[stamp_nx*stamp_nx*size*size]();
		double *point = new double[2 * num_p]();
		double *gal = new double[size*size]();
		double *gpow = new double[size*size]();
		double *psf = new double[size*size]();
		double *ppow = new double[size*size]();
		double *noise = new double[size*size]();
		double *pnoise = new double[size*size]();
		char log_inform[100], log_path[100], chip_path[200], data_path[200], buffer[200], h5_path[100], set_name[20];

		/* the data matrix data[i][j] */
		double *data_m = new double[data_rows*data_cols]();
		double **data = new double*[data_rows];
		for (i = 0; i < data_rows; i++)
		{
			data[i] = data_m + i*data_cols;
		}

		// initialize gsl
		//seed = myid * 84324 + 46331;
		seed = myid *3800 + 1401;
		gsl_rng_initialize(seed);
		
		string s;
		const char *str;
		// read shear, the input signal
		double *shear = new double[2 * shear_pairs](); // [...g1,...,..g2,...]
		fin.open("/lmc/selection_bias/parameters/shear.dat");
		i = 0;
		while (!fin.eof())
		{
			getline(fin, s);
			str = s.c_str();
			shear[i] = atof(str);
			i++;
		}
		fin.close();

		// read flux of galaxies
		double *flux = new double[chip_num*stamp_num];
		fin.open("/lmc/selection_bias/parameters/lsstmagsims.dat");
		i = 0;
		while (!fin.eof())
		{
			getline(fin, s);
			str = s.c_str();
			flux[i] = atof(str)/num_p;	
			i++;
			if (i == chip_num*stamp_num) // in case of the length of files is longger than the array
			{
				break;
			}
		}
		fin.close();

		// read ellipticity of galaxies
		double *ellip = new double[chip_num*stamp_num];
		sprintf(buffer, "/lmc/selection_bias/parameters/ellip_%d.dat", myid);
		fin.open(buffer);
		i = 0;
		while (!fin.eof())
		{
			getline(fin, s);
			str = s.c_str();
			ellip[i] = atof(str);
			i++;
			if (i == chip_num*stamp_num) // in case of the length of files is longger than the array
			{	
				break;
			}
		}
		fin.close();

		 //PSF
		create_psf(psf, psf_scale, size, psf_type);
		pow_spec(psf, ppow, size, size);		
		get_radius(ppow, &all_paras, thres, size, 1, psf_noise_sig);

		st = clock();

		g1 = shear[myid];
		g2 = shear[myid + shear_pairs];

		sprintf(h5_path, "/lmc/selection_bias/result/data/data_%d.hdf5", myid);
		sprintf(set_name, "/data");
		sprintf(log_path, "/lmc/selection_bias/logs/log_%d.dat", myid);

		for (i = 0; i < chip_num;  i++)
		{	
			s1 = clock();

			sprintf(log_inform, "Thread: %d, chip: %d, start.", myid, i);
			write_log(log_path, log_inform);

			for (j = 0; j < stamp_num; j++)
			{
				//create_points(point, num_p, max_radius);

				create_epoints(point, num_p, ellip[i*stamp_num + j]);

				convolve(gal, point, flux[i*stamp_num + j], size, num_p, 0, psf_scale, g1, g2, psf_type);

				addnoise(gal, size*size,  gal_noise_sig);

				stack(big_img, gal, j, size, stamp_nx, stamp_nx);

				get_radius(gal, &all_paras, 9999999999.*thres, size, 2, gal_noise_sig);

				pow_spec(gal, gpow, size, size);
				f_snr(gpow, &all_paras, size, 1);

				addnoise(noise, size*size, gal_noise_sig);				
				
				pow_spec(noise, pnoise, size, size);
				
				shear_est(gpow, ppow, pnoise, &all_paras, size);
				
				initialize(gal, size*size);
				initialize(gpow, size*size);
				initialize(point, num_p * 2);
				initialize(noise, size*size);
				initialize(pnoise, size*size);

				data[i*stamp_num + j][0] = g1;
				data[i*stamp_num + j][1] = g2;
				data[i*stamp_num + j][2] = all_paras.n1;
				data[i*stamp_num + j][3] = all_paras.n2;
				data[i*stamp_num + j][4] = all_paras.dn;
				data[i*stamp_num + j][5] = all_paras.du;
				data[i*stamp_num + j][6] = all_paras.dv;
				data[i*stamp_num + j][7] = all_paras.gal_osnr;
				data[i*stamp_num + j][8] = all_paras.gal_flux;
				data[i*stamp_num + j][9] = all_paras.gal_peak;
				data[i*stamp_num + j][10] = all_paras.gal_fsnr;
				data[i*stamp_num + j][11] = all_paras.gal_fsnr1;
				data[i*stamp_num + j][12] = all_paras.gal_fsnr4;
				data[i*stamp_num + j][13] = all_paras.gal_fsnr9;
				data[i*stamp_num + j][14] = all_paras.gal_snr;
				data[i*stamp_num + j][15] = flux[i*stamp_num + j];
				data[i*stamp_num + j][16] = 0.;
				data[i*stamp_num + j][17] = (double)(myid);
				data[i*stamp_num + j][18] = (double)(i);
				data[i*stamp_num + j][19] = (double)(j);
			}
			sprintf(chip_path, "!/lmc/selection_bias/%d/gal_chip_%04d.fits", myid, i);
			write_img(big_img, stamp_nx*size, stamp_nx*size, chip_path, myid,i);		
			
			initialize(big_img, stamp_nx*stamp_nx*size*size);	

			s2 = clock();
			sprintf(log_inform, "Thread: %d, chip: %d, done in %.2f s.", myid, i, (s2 - s1) / CLOCKS_PER_SEC);
			write_log(log_path, log_inform);
			if (myid == 0)
			{
				sprintf(buffer, "myid %d: chip %d done in %g \n", myid, i, (s2 - s1) / CLOCKS_PER_SEC);
				cout << buffer;
			}
		}
		write_h5(h5_path, set_name, data_rows, data_cols, data[0],NULL);
		ed = clock();
		if (myid == 0)
		{
			sprintf(buffer, "myid %d:  done in %g \n", myid, i, (ed - st) / CLOCKS_PER_SEC);
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
		delete[] ellip;
		gsl_rng_free();				
		MPI_Finalize();
		return 0;
}
