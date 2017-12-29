// exercise.cpp : 定义控制台应用程序的入口点。

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include "mpi.h"
#include "FQlib.h"

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
		ofstream fout;
		ifstream fin;

		/* 10 (g1,g2) points and each pairs contain 100 chips which cotians 10000 gals */
		int chip_num = 1000, stamp_num = 10000, shear_pairs = 10;

		int i, j, seed;
		int size = 52, num_p = 45, stamp_nx = 100,  psf_type = 2;
		double psf_scale = 4., max_radius = 9., st, ed, s1, s2;
		double g1=0., g2=0.;
		double gal_noise_sig = 380.64, psf_noise_sig = 0., thres = 2.;
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
		char chip_path[200], data_path[200], buffer[300], buffer1[200] ;

		// initialize gsl
		//seed = myid * 84324 + 46331;
		seed = myid * 2255234 + 3455331;
		gsl_rng_initialize(seed);
		
		string s;
		const char *str;
		// read shear, the input signal
		double *shear = new double[2 * shear_pairs](); // [...g1,...,..g2,...]
		fin.open("/home/hklee/work/selection_bias/parameters/shear.dat");
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
		fin.open("/home/hklee/work/selection_bias/parameters/lsstmagsims.dat");
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

		 //PSF
		create_psf(psf, psf_scale, size, psf_type);
		pow_spec(psf, ppow, size, size);		
		get_radius(ppow, &all_paras, thres, size, 1, psf_noise_sig);
	
		st = clock();

		g1 = shear[myid];
		g2 = shear[myid + shear_pairs];
		sprintf(buffer1, "myid %d  g1: %4f  g2: %4f \n", myid, g1, g2);
		cout << buffer1;

		for (i = 0; i < chip_num;  i++)
		{	
			s1 = clock();

			sprintf(chip_path, "/lmc/selection_bias/%d/gal_chip_%04d.fits", myid, i);
			sprintf(data_path, "/lmc/selection_bias/result/data/%d_gal_chip_%04d.dat", myid, i);

			fout.open(data_path);

			for (j = 0; j < stamp_num; j++)
			{
				create_points(point, num_p, max_radius);
				convolve(gal, point, flux[i*stamp_num+j], size, num_p, 0, psf_scale, g1, g2, psf_type);

				addnoise(gal, size*size, &all_paras, gal_noise_sig);

				stack(big_img, gal, j, size, stamp_nx, stamp_nx);

				get_radius(gal, &all_paras, 9999999999.*thres, size, 2, gal_noise_sig);

				pow_spec(gal, gpow, size, size);
				f_snr(gpow, &all_paras, size, 1);

				addnoise(noise, size*size, &all_paras, gal_noise_sig);
				pow_spec(noise, pnoise, size, size);
				
				shear_est(gpow, ppow, pnoise, &all_paras, size);
				
				initialize(gal, size*size);
				initialize(gpow, size*size);
				initialize(point, num_p * 2);
				initialize(noise, size*size);
				initialize(pnoise, size*size);
				sprintf(buffer, "%g %g %g %.6f %.6f %g %g %g %.6f %.6f %.6f %.6f %.6f %g %g %g %d %.3f %.3f %.3f %.3f %.3f %d %d \n", 0., 0., 0., all_paras.n1, g1, 0.,0.,0., all_paras.n2, g2, all_paras.dn, all_paras.du, all_paras.dv,
				0., 0., 0., all_paras.gal_size, all_paras.gal_flux , all_paras.gal_peak , all_paras.gal_snr, all_paras.gal_fsnr, all_paras.gal_osnr ,myid, i);
				fout << buffer;
			}

			write_img(big_img, stamp_nx*size, stamp_nx*size, chip_path);		
			
			initialize(big_img, stamp_nx*stamp_nx*size*size);	
			fout.close();
			s2 = clock();

			sprintf(buffer1, "myid %d: chip %d done in %g \n", myid, i, (s2 - s1) / CLOCKS_PER_SEC);
			cout << buffer1;
		}
		
		ed = clock();

		sprintf(buffer1, "myid %d:  done in %g \n", myid, i, (ed - st) / CLOCKS_PER_SEC);
		cout << buffer1;

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
		gsl_rng_free();		
		
		MPI_Finalize();
		return 0;
}
