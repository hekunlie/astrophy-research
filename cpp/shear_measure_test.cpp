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
	int total_chip_num = 50, chip_num, stamp_num = 10000, shear_pairs = 14;
	/* remember to change the data_cols when you change the number of estimators recorded */
	int i, j, seed, data_rows, data_cols = 7, chip_id, shear_id, detect_label;
	int size = 52, num_p = 40, stamp_nx = 100, psf_type = 2;
	double psf_scale = 4., max_radius = 8., st, ed, s1, s2, s3, s4, s5, s6;
	double g1 = 0., g2 = 0.;
	double gal_noise_sig = 0, psf_noise_sig = 0., scale = 2.;
	int total_num = total_chip_num * stamp_num;

	chip_num = total_chip_num / (numprocs / shear_pairs);
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


	/* the data matrix data[i][j] */
	double *data_m = new double[data_rows*data_cols]();
	double **data = new double*[data_rows];
	for (i = 0; i < data_rows; i++)
	{
		data[i] = data_m + i*data_cols;
	}

	// initialize gsl
	int sss1, sss2;
	sss1 = 581140;
	sss2 = 716110;
	seed = myid *sss1 + sss2;
	//seed = myid *380 + 1401;// no bias
	gsl_rng_initialize(seed);


	const char *str;
	// read shear, the input signal
	double *shear = new double[2 * shear_pairs](); // [...g1,...,..g2,...]
	fin.open("/home/hklee/work/test/shear.dat");
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
	double *flux = new double[5000000];
	double *mag = new double[5000000];

	//sprintf(para_path, "/mnt/ddnfs/data_users/hkli/selection_bias_64/parameters/para_%d.hdf5", shear_id);
	//sprintf(set_name1, "/flux");
	//sprintf(set_name2, "/mag");
	//read_h5(para_path, set_name1, flux, set_name2, mag, NULL, NULL);


	//PSF
	create_psf(psf, psf_scale, size, psf_type);

	pow_spec(psf, ppow, size, size);

	get_radius(ppow, &all_paras, scale, 1, psf_noise_sig);

	st = clock();

	g1 = shear[shear_id];
	g2 = shear[shear_id + shear_pairs];
	double t1=0, t2=0, t3=0, t4=0, t5=0;
	for (i = 0; i < chip_num; i++)
	{
		

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

			convolve(gal, point, 100, size, num_p, 0, psf_scale, g1, g2, psf_type);

			//addnoise(gal, size*size, gal_noise_sig);
			if(j<stamp_nx*stamp_nx)
			stack(big_img, gal, j, size, stamp_nx, stamp_nx);			

			//get_radius(gal, &all_paras, 9999999999.*thres, 2, gal_noise_sig);
			detect_label = galaxy_finder(gal, &all_paras, false);

			pow_spec(gal, gpow, size, size);
			f_snr(gpow, &all_paras, 1);

			//addnoise(noise, size*size, gal_noise_sig);

			pow_spec(noise, pnoise, size, size);
			
			shear_est(gpow, ppow, pnoise, &all_paras);
			
			data[i*stamp_num + j][0] = g1;
			data[i*stamp_num + j][1] = g2;
			data[i*stamp_num + j][2] = all_paras.n1;
			data[i*stamp_num + j][3] = all_paras.n2;
			data[i*stamp_num + j][4] = all_paras.dn;
			data[i*stamp_num + j][5] = all_paras.du;
			data[i*stamp_num + j][6] = all_paras.dv;

		}
		if (myid == 0 && i == 0)
		{
			sprintf(chip_path, "!/home/hklee/work/test/gal_chip_%04d.fits", chip_id + i);
			write_img(big_img, size*stamp_nx, size*stamp_nx, chip_path);
		}
		
		initialize_arr(big_img, stamp_nx*stamp_nx*size*size);	
		
	}

	sprintf(h5_path, "/home/hklee/work/test/data/data_%d_%d.hdf5", shear_id, myid / shear_pairs);
	sprintf(set_name1, "/data");
	write_h5(h5_path, set_name1, data_rows, data_cols, data[0], NULL);
	ed = clock();
	if (myid == 0)
	{
		sprintf(buffer, "myid %d:  done in %g \n", myid, (ed - st) / CLOCKS_PER_SEC);
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


/*
	int size = 40, s_num, area_s=0, area_e;
	para paras;
	paras.img_x = size;
	paras.img_y = size;
	paras.detect_thres = 4.5;
	paras.stamp_size = size;

	double t1, t2;
	double *img = new double[size*size]{};
	double* img_t = new double[size*size]{};
	int *s_x = new int[size*size]{};
	int *s_y = new int[size*size]{};
	double *s_c = new double[8*paras.max_source]{};
	int detect;
	char buffer[60];
	cout << "starting..." << endl;
	sprintf(buffer, "/home/hkli/temp/test.fits");
	read_img(img, buffer);
	cout << "read img..." << endl;
	s_num = source_detector(img, s_x, s_y, s_c, &paras, false);
	detect = galaxy_finder(img, &paras, false);
	cout << paras.gal_py << " " << paras.gal_px << endl;
	cout << detect<<"finished, find: "<< s_num << endl;
	for (int m = 0; m < s_num; m++)
	{
		cout << "AREA: "<<s_c[8 * m] << endl;
		area_e = area_s + s_c[m * 8];

		for (int i = area_s; i < area_e; i++)
		{
			img_t[s_x[i] + s_y[i] * size] = 1;
		}
		img_t[(int)s_c[8 * m + 1] * size + (int)s_c[8 * m + 2]] += 2;
		area_s = area_e;
	}
	img_t[paras.gal_py*size + paras.gal_px] += 1;
	sprintf(buffer, "!/home/hkli/temp/t_img.fits");
	write_img(img_t, size, size, buffer);
	
	/*double t1, t2;
	t1 = clock();
	f_snr(img, &paras);
	cout << paras.gal_flux2 << " " << paras.gal_flux_alt << endl;
	int x[20]{ -1,  0,  1, -2, -1,  0,  1,  2, -2, -1,  1,  2, -2, -1,  0,  1,  2, -1,  0,  1 };
	int y[20]{ -2, -2, -2, -1, -1, -1, -1, -1,  0,  0,  0,  0,  1,  1,  1,  1,  1,  2, 2, 2 };
	double fz[20], fit_paras[6];
	int xc = size / 2, i;
	cout << xc << endl;
	for (i = 0; i < 20; i++)
	{
		fz[i] = img[(xc + y[i])*size + xc + x[i]];
	}

	hyperfit_5(fz, fit_paras,&paras);
	for (i = 0; i < 6; i++)
	{
		cout << fit_paras[i] << " ";
	}
	cout << endl;
	cout << pow(10, fit_paras[0]) <<" "<<img[xc*size+xc]<< endl;
	/*
	}
	t2 = clock();
	cout << (t2 - t1) / CLOCKS_PER_SEC;
	MPI_Finalize();
	return 0;
}*/
