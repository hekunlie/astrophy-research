#ifndef FQLIB_H
#define FQLIB_H

#pragma once
#include <iostream>
#include<fstream>
#include<string>
#include<string.h>
#include<sstream>
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
#include<stdlib.h>
#include<algorithm> // sort()
#include<functional> // std::less, std::greater..

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
	double gal_peak, gal_hlr, gal_flux, gal_hflux, gal_fluxsq;
	double gal_flux2, gal_flux_alt, gal_snr, gal_osnr, gal_noise_sig;

	double n1, n2, dn, du, dv, dp1, dp2;
	double t1, t2, t3, t4;


	/*parameters for detection which should be initialized before */
	int stamp_size; /* the stamp size for get_radius() */
	int img_x, img_y; /* the size of chip image for the 'source_detector()' and galaxy_finder()*/
	int area_thres=6; /* the minimun pixels for a detection */
	double detect_thres; /* the threshold of pixel value of source */
	double noise_sig;
	int max_source = 1000; /* the maximum of sources allowed in each chip, changeable */
	double max_distance= 8.;/* the max distance of peak away from the center of the source candidate */


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

	/* for calculating the exponential exp() quickly */
	double exp_val[15] = { exp(1.), exp(2.), exp(3.), exp(4.), exp(5.), exp(6.), exp(7.), exp(8.), exp(9.),  exp(10.), exp(11.), exp(12.), exp(13.), exp(14.), exp(15.), };
};


//using namespace std;
const double Pi = 3.1415926535897932384626433832795;
extern const gsl_rng_type *T;
extern gsl_rng *rng;
extern std::ofstream loggers;

/********************************************************************************************************************************************/
/* file reading and writting*/
/********************************************************************************************************************************************/
void write_log(char *filename, char *inform); 
/* write char to log file */

void read_para(const std::string path, const std::string name, int &para);
void read_para(const std::string path, const std::string name, double &para);
void read_para(const std::string path, const std::string name, float &para);

void read_h5(char *filename, char *set_name1, double *matrix1, char*set_name2, double *matrix2, char*set_name3, double*matrix3);
/* read hdf5 file 
	set_name and matrix should be used in pair (developing)
*/

void write_h5(char *filename, char *set_name, int row, int column, double*d_matrix, int *i_matrix);
/* write to hdf5 file with double or integer matrix 
	only one of "d_matrix" and "i_matrix" should be inputted each time
*/

void read_img(DATA_TYPE *arr, char *path);
/* read fits file, the preciseion dependences on the DATA_TYPE */

void write_img(DATA_TYPE *img, int ysize, int xsize, char *filename);
/* write the array to  fits file, the preciseion dependences on the DATA_TYPE 
	the size of each axises should be inputted，ARR(y, x)
*/

void stack(double *big_arr, double *stamp, int tag, int size, int row, int col); 
/* from stamps to a integrate image 
	big_arr:   the big image that will contain all the stamps
	stamps: the small array that will be stacked into the big image
	tag: the location label of stamps in tha big image, 
				| 0,1,2,..., i |
				| i+1,........  |
				....
	size: the pixel of the side of each square stamp
	row and col:  how many stamps in row and column
*/

void segment(double *big_arr, double *stamp, int tag, int size, int row, int col);
/* to cut the specific square area in the big_arr
	see the annotation of stack()	
*/

/********************************************************************************************************************************************/
/* operations on the image */
/********************************************************************************************************************************************/
void create_points(double *point, int num_p, double radius);
void create_epoints(double *point, int num_p, double ellip);
void create_psf(double*in_img, double scale, int size, int psf);
void convolve(double *in_img, double * points, double flux, int size, int num_p, int rotate, double scale, double g1, double g2, int psf);
void pow_spec(double *in_img, double *out_img, int column, int row);

void get_radius(double *in_img, para *paras, double scale, int type, double sig_level);
void get_psf_radius(const double *psf_pow, para*para, const double scale);
/*measure the size of psf power spectrum for the \beta parameter in the measurement.
	power of k=0 may be not the maximun, be careful!!!! */

int source_detector(double *source_img, int *soucrce_x, int*source_y, double *source_paras,para* paras, bool cross);
/* operates on the copy,
	if the method finds too many sources ( > para.max_source), the overflows will be ignored.
	source_img: the inputted array where to find the source galaxies
	source_x, _y:  the array to store the coordinates of sources detected
	source_paras: the array to store the parameters of sources detected,
						   8 elemets for each source, [....,area, peak_y, peak_x, peak_val, half_light_area, total_flux, half_light_flux, flux_sq,...]
	cross: boolean, True for detection on the nearest four pixels, "+", upper, lower, left, right
							False for detecion on the nearest eight pixels, "x" and "+"  
	return : int, the total number of detection */


int galaxy_finder(double *stamp_arr, para *paras, bool cross);
/* to indentify the galaxy on each stamp basing on source_detector(), because of many detections on it
	the biggest source which peaks in the central circle with a radius of 6 pixels.	
	return: int, "-1" means no detection
	*/

void addnoise(double *image, int pixel_num, double sigma); 
/* add Gaussian noise to an array */

void initialize_arr(double *array, int size);
/* set every elements to zero*/

/********************************************************************************************************************************************/
/* Fourier Quad */
/********************************************************************************************************************************************/
void possion_subtraction(double *image_pow, para *paras, int edge);
void noise_subtraction(double *image_pow, double *noise_pow, para *paras, const int edge, const int possion);
void shear_est(double *gal_img, double *psf_img, para *paras);
void snr_est(const double *image, para *paras, int fit);
/* if fit=2 for both flux2 and flux_alt estimations 
	else just estimate the flux2 
*/

/********************************************************************************************************************************************/
/* fitting */
/********************************************************************************************************************************************/
void smooth(double *image, const double *coeffs, para *paras);
/* smooth all the region */
void smooth(double *image, const double *psf_pow, const double *coeffs, para *paras);
/* the image will be repalced by the smoothed one. 
	the psf_pow and the paras->psf_thres_pow are the mask and  threshold to label the region where to be smoothed
	to fit the curve: a1 + a2*x +a3*y + a4*x^2 +a5*x*y + a6*y^2  */
void hyperfit_5(double *data, double*fit_para, para *paras);


/********************************************************************************************************************************************/
/* general methods */
/********************************************************************************************************************************************/
void initialize_para(para *paras);
/* set the "gal_: parameters zero */

void set_bin(const double *data, const int data_num, double * bins, const int bin_num);
void set_bin(const float *data, const int data_num, float * bins, const int bin_num);
void set_bin(const int *data, const int data_num, int * bins, const int bin_num);

void histogram(const double *data, const double *bins, int *num, const int data_num, const int bin_num);
void histogram(const float *data, const float *bins, int *num, const int data_num, const int bin_num);
void histogram(const int *data, const  int *bins, int *num, const  int data_num, const  int bin_num);

void histogram2d(const double *data_y, const double*data_x, const double *bin_y, const double *bin_x, int *num, const int data_num, const int ybin_num, const  int xbin_num);
void histogram2d(const float *data_y, const float*data_x, const float *bin_y, const float *bin_x, int *num, const int data_num, const  int ybin_num, const int xbin_num);
void histogram2d(const int *data_y, const int*data_x, const int *bin_y, const int *bin_x, int *num, const int data_num, const int ybin_num, const int xbin_num);

void sort_double(double *arr, int size, int order);
/* sort the double array according to the order, order =1 for ascend, else for descend*/
void sort_float(float *arr, int size, int order);
void sort_int(int *arr, int size, int order);

int com_double_ascend(const double a, const double b);
/* the compare function for the qsort() method */
int com_float_ascend(const float a, const float b);
int com_int_ascend(const int a, const int b);
int com_double_descend(const  double a, const double b);
int com_float_descend(const float a, const float b);
int com_int_descend(const int a, const int b);

void get_time(char *str_time, int length);
/* get the current time.
	the length of str_time should be larger than 40 */

/*double qucik_exp(double x, double precision, para* paras);*/

/********************************************************************************************************************************************/
/* GSL library */
/********************************************************************************************************************************************/
void gsl_rng_initialize(int seed);
void gsl_rng_free();

#endif // !FQLIB_H

