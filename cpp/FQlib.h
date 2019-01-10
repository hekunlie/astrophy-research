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
#include<algorithm> // sort(), std::max()
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
	double gal_peak, gal_hlr, gal_flux, gal_hflux, gal_fluxsq, gal_total_flux;
	double gal_flux2, gal_flux_alt, gal_snr, gal_osnr, gal_noise_sig, gal_flux2_new;
	double gal_size_ext[5];
	double gal_flux_ext[5];
	double gal_flux2_ext[5];

	double n1, n2, dn, du, dv, dp1, dp2;
	double t1, t2, t3, t4;


	/*parameters for detection which should be initialized before */
	int stamp_size; /* the stamp size for get_radius() */
	int img_x, img_y; /* the size of chip image for the 'source_detector()' and galaxy_finder()*/
	int area_thres=6; /* the minimun pixels for a detection */
	double detect_thres; /* the threshold of pixel value of source */
	double noise_sig;
	int max_source = 20; /* the maximum of sources allowed in each chip, changeable */
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
/* read the parameters ("name") value from parameter file */

void read_text(const std::string path, double *arr, const int read_lines);
void read_text(const std::string path, float *arr, const int read_lines);
void read_text(const std::string path, int *arr, const int read_lines);
/* read data from txt file which should be just one column.
	the read_lines limits the maximum lines to read
*/

void read_h5(const char *filename, const char *set_name, double *arr);
void read_h5(const char *filename, const char *set_name, float *arr);
void read_h5(const char *filename, const char *set_name, int *arr);
void write_h5(const char *filename, const char *set_name, const double *arr, const int row, const int column);
void write_h5(const char *filename, const char *set_name, const float *arr, const int row, const int column);
void write_h5(const char *filename, const char *set_name, const int *arr, const int row, const int column);
/* read and write the hdf5 file */

void read_fits(const char *filename, double *arr);
void read_fits(const char *filename, float *arr);
void read_fits(const char *filename, int *arr);
void write_fits(const char *filename, double *img, const int ysize, const int xsize);
void write_fits(const char *filename, float *img, const int ysize, const int xsize);
void write_fits(const char *filename, int *img, const int ysize, const int xsize);
/* read and write the array to  fits file, 
	be careful with the datetype "TINT" and "LONG_IMG"!!! 
	the length of INT may be different in different platform,
	however, the "TINT32BIT" doesn't work with "LONG_IMG".

	the size of each axises should be inputted，ARRAY(y, x)
*/

void stack(double *big_arr, const double *stamp, const int tag, const int size, const int row, const int col);
void stack(float *big_arr, const float *stamp, const int tag, const int size, const int row, const int col);
void stack(int *big_arr, const int *stamp, const int tag, const int size, const int row, const int col);
/* 
    from stamps to a integrate image 
	big_arr:   the big image that will contain all the stamps
	stamps: the small array that will be stacked into the big image
	tag: the location label of stamps in tha big image, 
				| 0,1,2,..., i |
				| i+1,........  |
				....
	size: the pixel of the side of each square stamp
	row and col:  how many stamps in row and column
*/

void segment(const double *big_arr, double *stamp, const int tag, const int size, const int row, const int col);
void segment(const float *big_arr, float *stamp, const int tag, const int size, const int row, const int col);
void segment(const int *big_arr, int *stamp, const int tag, const int size, const int row, const int col);
/* to cut the specific square area in the big_arr
	see the annotation of stack()	
*/

/********************************************************************************************************************************************/
/* operations on the image */
/********************************************************************************************************************************************/
void create_points(double *point, const int num_p, const double radius);
void create_epoints(double *point, const int num_p, const double ellip);
void create_psf(double*in_img, const double scale, const int size, const int psf);
void convolve(double *in_img, double * points, double flux, int size, int num_p, int rotate, double scale, double g1, double g2, int psf);

void pow_spec(const double *in_img, double *out_img, const int column, const int row);
void pow_spec(const float *in_img, float *out_img, const int column, const int row);

void get_radius(double *in_img, para *paras, double scale, int type, double sig_level);
void get_psf_radius(const double *psf_pow, para*para, const double scale);
void get_psf_radius(const float *psf_pow, para*para, const float scale);
/*measure the size of psf power spectrum for the \beta parameter in the measurement.
	power of k=0 may be not the maximun, be careful!!!! */

int source_detector(const double *source_img, int *soucrce_x, int*source_y, double *source_paras,para* paras, bool cross);
int source_detector(const float *source_img, int *soucrce_x, int*source_y, float *source_paras, para* paras, bool cross);
/* operates on the copy,
	if the method finds too many sources ( > para.max_source), the overflows will be ignored.
	source_img: the inputted array in where to find the source galaxies
	source_x, _y:  the array to store the coordinates of sources detected
	source_paras: the array to store the parameters of sources detected,
						   8 elemets for each source, [....,area, peak_y, peak_x, peak_val, half_light_area, total_flux, half_light_flux, flux_sq,...]
	cross: boolean, True for detection on the nearest four pixels, "+", upper, lower, left, right
							False for detecion on the nearest eight pixels, "x" and "+"  
	return : int, the total number of detection 
*/


int galaxy_finder(const double *stamp_arr, int *check_mask, para *paras, bool cross);
int galaxy_finder(const float *stamp_arr, int *check_mask, para *paras, bool cross);
/* to indentify the galaxy on each stamp basing on source_detector(), because of many detections on it
	the biggest source which peaks in the central circle with a radius of 6 pixels.	
	return: int, "-1" means no detection
*/

int edge_extend(int *mask, const int *source_y, const int* source_x, const int source_id, const int source_len, para *paras, const int iters);
/*	"stamp_size" in the structure paras will be used !!! 

	extend the border of a source galaxy by 1 pixel each time
	the coordinates source_x(y), will be copied to a new array.

	mask: array, on which the source and extendion edge will be labeled, 1 means source pixel
	source_y(x): array of source coordinates from source_detector(), may be more than one source
	source_len : the length of the target source
	source_id : the start of the source coordinates in the source_y(x), because there may be more than one source.
	iters:  iterations, each time will extend the edge by one pixel.
	return: the area of the extended source
*/

void addnoise(double *image, const int pixel_num, const double sigma);
void addnoise(float *image, const int pixel_num, const float sigma);
/* add Gaussian noise to an array */

void initialize_arr(double *array, const int size);
void initialize_arr(float *array, const int size);
void initialize_arr(int *array, const int size);
/* set every elements to zero*/

void normalize_arr(double *arr, const int size);
void normalize_arr(float *arr, const int size);
/* normalize the PSF power spectrum,
	divide each pixel by the peak 
*/

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
	to fit the curve: a1 + a2*x +a3*y + a4*x^2 +a5*x*y + a6*y^2  
*/

void smooth_real(double*image, const double *coeffs, para *paras);
/* smooth the image by fitting a polynomial */

void hyperfit_5(const double *data, double*fit_para, para *paras);


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
void sort_float(float *arr, int size, int order);
void sort_int(int *arr, int size, int order);
/* sort the double array according to the order, order =1 for ascend, else for descend*/


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

