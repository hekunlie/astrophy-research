#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include "mpi.h"
#include "FQlib.h"
#include<stdio.h>

//#define TRANS_S_STD 0.5
using namespace std;
para all_paras;


int main(int argc, char*argv[])
{
	ofstream fout;
	ifstream fin;

	all_paras.t1 = 0;
	all_paras.t2 = 0;
	/* 10 (g1,g2) points and each pairs contain 100 chips which cotians 10000 gals */
	int chip_num = 1, stamp_num = 10000, shear_pairs = 10;
	int myid = 0;
	int i, j, seed;
	int size = 2043, y_size=4606, num_p = 100, stamp_nx = 100, psf_type = 2;
	double psf_scale = 5., max_radius = 11., st, ed;
	double g1 = 0., g2 = 0.;
	double gal_noise_sig = 380.64, psf_noise_sig = 0., thres = 2.;
	all_paras.gal_noise_sig = gal_noise_sig;
	all_paras.psf_noise_sig = psf_noise_sig;

	//double *big_img = new double[stamp_nx*stamp_nx*size*size]();
	double *point = new double[2 * num_p]();
	double *points_r = new double[2 * num_p]();
	double *gal = new double[y_size*size]();
	double *gpow = new double[size*size]();
	double *psf = new double[size*size]();
	double *ppow = new double[size*size]();
	double *noise = new double[size*size]();
	double *pnoise = new double[size*size]();
	char chip_path[200], data_path[200], buffer[300], buffer1[200];
	int *chain = new int[2*size * y_size + 1]();

	seed = 155301;
	gsl_rng_initialize(seed);

	//PSF
	//create_psf(psf, psf_scale, size, psf_type);
	//pow_spec(psf, ppow, size, size);
	//get_radius(ppow, &all_paras, thres, size, 1, psf_noise_sig);

	st = clock();

	g1 = 0.03;
	g2 = 0.04;

	double s1, s2, s3, s4, s5, s6, s7, s8, flux_m = 200.;
	double d1 = 0., d2 = 0., d3 = 0., d4 = 0, d5 = 0, r1, rd, rs;
	int k, ii, jj;
	int rows = 1000;
	int columns = 27;
	//for (i = 0; i < 42; i++)
	//{
	//	cout <<"I: "<<i<<"\t"<<i/14<< endl;
	//}

	for (i = 0; i < 1; i++)
	{
		s1 = clock();
		sprintf(chip_path, "/home/hklee/831549p_1OFCBC.fits");
		cout << chip_path << endl;
		sprintf(data_path, "/home/hklee/data.hdf5");
		cout << chip_path << endl;
		read_img(gal,chip_path);
		cout << "read img" << endl;
		detector(gal, chain, 100, y_size, size);
		cout << "detecting img" << endl;
		sprintf(chip_path, "/data");
		s2 = clock();
		write_h5(data_path, chip_path, 1, 2 * y_size*size + 1, NULL,chain);

		sprintf(buffer1, "myid %d:  %g %g %g %g \n", myid, g1, d1/d3, g2, (s2-s1) / CLOCKS_PER_SEC);
		cout << buffer1;
	}

	ed = clock();

	sprintf(buffer1, "myid %d:  done in %g \n", myid, i, (ed - st) / CLOCKS_PER_SEC);
	cout << buffer1;

	//delete[] big_img;
	delete[] point;
	delete[] gal;
	delete[] gpow;
	delete[] psf;
	delete[] ppow;
	delete[] noise;
	delete[] pnoise;
	delete[] points_r;
	gsl_rng_free();
	delete[] chain;
	return 0;
}
