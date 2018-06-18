#ifndef FQLIB_H
#define FQLIB_H

#pragma once
#include <iostream>
#include<string.h>
#include<fstream>
#include <iomanip>
#include <cmath>
#include "fitsio.h"
#include<stdio.h>
#include <fftw3.h>
#include <ctime>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include<gsl/gsl_cblas.h>
#include<hdf5.h>
//#include<mkl.h>

#define PRECISION
#ifdef  PRECISION
typedef double DATA_TYPE;
#define  IMG_PRE DOUBLE_IMG
#define  T_IMG TDOUBLE
#else
typedef float DATA_TYPE;
#define IMG_PRE FLOAT_IMG
#define T_IMG TFLOAT
#endif 
struct para
{
	int psf_size, psf_px, psf_py;
	double psf_peak, psf_hlr, psf_flux, psf_fluxsq, psf_noise_sig, psf_pow_thres = 0.0001;

	int gal_size, gal_hsize, gal_px, gal_py;
	double gal_peak, gal_hlr, gal_flux, gal_hflux, gal_fluxsq, gal_flux2,gal_flux_alt, gal_snr, gal_osnr, gal_noise_sig;

	double n1, n2, dn, du, dv, dp1, dp2;
	double t1, t2, t3, t4;


	/*parameters for detection*/
	int img_size;
	int img_x, img_y; /* for the 'detection()'*/
	int area_thres;
	double noise_sig;
	double detect_thres; /* the threshold value of source */
	int max_source = 100; /* the maximum of sources in each chip*/


	/* hyper_fit_5 matrix elements of order 2 of xy polynomials */
	/* this is the final matrix and the data value is the only things needed */
	double fit_matrix[6][20] =
	{
		{ -0.0530303,0.0113636,-0.0530303,-0.0530303,0.1401515,0.2045455,0.1401515,-0.0530303,0.0113636,0.2045455,0.2045455,0.0113636,-0.0530303,0.1401515,0.2045455,0.1401515,-0.0530303,-0.0530303,0.0113636,-0.0530303 },
		{ -0.0294118,0.0,0.0294118,-0.0588235,-0.0294118,0.0,0.0294118,0.0588235,-0.0588235,-0.0294118,0.0294118,0.0588235,-0.0588235,-0.0294118,0.0,0.0294118,0.0588235,-0.0294118,0.0,0.0294118 },
		{ -0.0588235,-0.0588235,-0.0588235,-0.0294118,-0.0294118,-0.0294118,-0.0294118,-0.0294118,0.0,0.0,0.0,0.0,0.0294118,0.0294118,0.0294118,0.0294118,0.0294118,0.0588235,0.0588235,0.0588235 },
		{ 0.0088745,-0.0172078,0.0088745,0.0517316,-0.0265152,-0.0525974,-0.0265152,0.0517316,0.0399351,-0.0383117,-0.0383117,0.0399351,0.0517316,-0.0265152,-0.0525974,-0.0265152,0.0517316,0.0088745,-0.0172078,0.0088745 },
		{ 0.0555556,0.0,-0.0555556,0.0555556,0.0277778,0.0,-0.0277778,-0.0555556,0.0,0.0,0.0,0.0,-0.0555556,-0.0277778,0.0,0.0277778,0.0555556,-0.0555556,0.0,0.0555556 },
		{ 0.0517316,0.0399351,0.0517316,0.0088745,-0.0265152,-0.0383117,-0.0265152,0.0088745,-0.0172078,-0.0525974,-0.0525974,-0.0172078,0.0088745,-0.0265152,-0.0383117,-0.0265152,0.0088745,0.0517316,0.0399351,0.0517316 },
	};

};


using namespace std;
const double Pi = 3.1415926535897932384626433832795;
extern const gsl_rng_type *T;
extern gsl_rng *rng;
extern ofstream loggers;

void write_log(char *filename, char *inform);
void read_h5(char *filename, char *set_name1, double *matrix1, char*set_name2, double *matrix2, char*set_name3, double*matrix3);
void write_h5(char *filename, char *set_name, int row, int column, double*d_matrix, int *i_matrix);
void read_img(DATA_TYPE *arr, char *path);
void write_img(DATA_TYPE *img, int ysize, int xsize, char *filename);
void pow_spec(double *in_img, double *out_img, int column, int row);
void get_radius(double *in_img, para *paras, double scale, int type, double sig_level);
void detector(double *source_img, int *soucrce_x, int*source_y, double *source_paras, para* paras, bool cross);
void convolve(double *in_img, double * points, double flux, int size, int num_p, int rotate, double scale, double g1, double g2, int psf);
void shear_est(double *gal_img, double *psf_img, double *noise_img, para *paras);
void create_points(double *point, int num_p, double radius);
void create_epoints(double *point, int num_p, double ellip);
void create_psf(double*in_img, double scale, int size, int psf);
void get_psf_thres(double *ppsf, para*paras);
void initialize(double *array, int size);
void stack(double *container, double *stamp, int tag, int size, int row, int col);
void segment(double *chip, double *stamp, int tag, int size, int row, int col);
void addnoise(double *image, int pixel_num,  double sigma);
void f_snr(double *image, para *paras);
void gsl_rng_initialize(int seed);
void gsl_rng_free();
void smooth(double *image,double*fit_image, double *psf_pow, double *coeffs, para *paras);
void hyperfit_5(double *data,double*fit_para, para *paras);
#endif // !FQLIB_H

