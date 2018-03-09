#ifndef FQLIB_H
#define FQLIB_H

#pragma once
#include <iostream>
#include<fstream>
#include <iomanip>
#include <cmath>
#include "fitsio.h"
#include<stdio.h>
#include <fftw3.h>
#include <ctime>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include<hdf5.h>

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
	int psf_szie, psf_px, psf_py;
	double psf_peak, psf_hlr, psf_flux, psf_fluxsq, psf_noise_sig;

	int gal_size, gal_px, gal_py;
	double gal_peak, gal_hlr, gal_flux, gal_fluxsq, gal_fsnr_c, gal_snr, gal_osnr, gal_noise_sig;

	double n1, n2, dn, du, dv, dp1, dp2;
	double t1, t2;
	double noise_sig;
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
void get_radius(double *in_img, para *paras, double scale, int size, int type, double sig_level);
void detector(double *source_img, int *soucrce_chain, double thres, int y_size, int x_size);
void convolve(double *in_img, double * points, double flux, int size, int num_p, int rotate, double scale, double g1, double g2, int psf);
void shear_est(double *gal_img, double *psf_img, double *noise_img, para *paras, int size);
void create_points(double *point, int num_p, double radius);
void create_epoints(double *point, int num_p, double ellip);
void create_psf(double*in_img, double scale, int size, int psf);
void initialize(double *array, int size );
void stack(double *container, double *stamp, int tag, int size, int row, int col);
void segment(double *chip, double *stamp, int tag, int size, int row, int col);
void addnoise(double *image, int pixel_num,  double sigma);
void f_snr(double *image, para *paras, int size);
void gsl_rng_initialize(int seed);
void gsl_rng_free();

#endif // !FQLIB_H

