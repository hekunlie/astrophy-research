// please compile it with C++11 standard(2011)

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
#include<gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include<hdf5.h>
#include<stdlib.h>
#include<algorithm> // sort(), std::max()
#include<functional> // std::less, std::greater..
#include<ciso646> // for "and, not, or, ..."

// something relates to the shear measurement and source dection
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

	double gal_e1, gal_e2;
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

// something relates to the correlation function calculation
struct pts_info
{
	int idy, idx; // block id of the point
	double y, x; // the coordinates of the point

	double scale; // the length of the side of the square blocks
	int ny, nx; // the number of blocks along each axis 
	int blocks_num; // the numbers of total blocks
};


const double Pi = 3.1415926535897932384626433832795;
extern const gsl_rng_type *T;
extern gsl_rng *rng;
extern std::ofstream loggers;

/********************************************************************************************************************************************/
/* array operation by Fortran*/
/********************************************************************************************************************************************/
extern "C"
{	
	// the addresses of the variables, array or int or float..., should be passed to the subroutines in Fortran
	void arr_pow_d_(double *arr_in,double *arr_out, int *size, double *alpha, double *beta, double *pows);
	void arr_pow_f_(float *arr_in, float *arr_out, int *size, float *alpha, float *beta, float *pows);
	void arr_base_d_(double *arr_in, double *arr_out, int *size, double *alpha, double *beta, double *base);
	void arr_base_f_(float *arr_in, float *arr_out, int *size, float *alpha, float *beta, float *base);
	/* 
		calculate (alpha*arr_in +beta)^pows and base^(alpha*arr_in + beta)
		results will be stored in arr_out
		except for the array, all the variables are passed by address, "&var"
		arr_in: array
		arr_out: array, the result
		size: int, the length of array
		alpha, beta, pows: double or float, 
	*/

	void arr_exp_d_(double *arr_in, double *arr_out, int *size, double *alpha, double *beta);
	void arr_exp_f_(float *arr_in, float *arr_out, int *size, float *alpha, float *beta);
	/*
		calculate EXP(alpha*arr_in +beta),
		results will be stored in arr_out
		except for the array, all the variables are passed by address, "&var"
		arr_in: array
		arr_out: array, the result
		size: int, the length of array
		alpha, beta: double or float,
	*/

	void arr_sum_d_(double *arr_in, double *total, int *size, double *alpha, double *beta);
	void arr_sum_f_(float *arr_in,  float *total, int *size, float *alpha, float *beta);
	/*
		calculate SUM(alpha*arr_in +beta),
		results will be stored in result
		except for the array, all the variables are passed by address, "&var"
		arr_in: array
		total: float or double, the sum
		size: int, the length of array
		alpha, beta: double or float,
	*/
}

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
void convolve(double *in_img, const double * points, const double flux, const int size, const int num_p, const int rotate, const double scale, const double g1, const double g2, const int psf, const int flag, para *paras);

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
void snr_est(const double *image, para *paras, int fit);
/* if fit=2 for both flux2 and flux_alt estimations
	else just estimate the flux2
*/
void possion_subtraction(double *image_pow, para *paras, int edge);
void noise_subtraction(double *image_pow, double *noise_pow, para *paras, const int edge, const int possion);
void shear_est(double *gal_img, double *psf_img, para *paras);
void ellip_est(const double *gal_img, const int size, para*paras);

//void find_block(const int ny, const int nx, const double y, const double x, const double scale, const int num_y , const int num_x,)

/********************************************************************************************************************************************/
/* random */
/********************************************************************************************************************************************/

double rand_gauss(double sigma, double mean);
/* return a double from the normal distribution with sigma and mean. 
*/

double rand_uniform(double start, double end);
/* return a double, [ start, end ), with a unifrom distribution.
*/

void rand_shuffle(double *seq, int length);
void rand_shuffle(float *seq, int length);
void rand_shuffle(int *seq, int length);
/* 	shuffle the value of the elements for choosing elements without repeating.
	then one can choose arbitrary elements from the disordered array.
	
	seq: array
	length: the length of the array
	
	i.e.  seq contains value in [a,b], then one can shuffle it with this method,
		  and choose n-elements from it as the random-chosen labels of other array.
		  
		  rand_shuffle(seq, length);
		  for(i=0;i<..;i++)
		  {	....
			some_method(array[seq[i]]);
			....
		  }
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

void poly_fit_1d(const double *x, const double *fx, const double *fx_err, const int data_num, double *coeffs, int weight);
/* polynomial fitting by GSL
	FX = c + m*X
	x: array, coordinates
	fx: array, measured values
	fx_err: array, uncertainty (\sigma), of which the reciprocal will be the weights of data points.
	data_num: the number of data points
	coeffs: array[4], [c, sig_c, m, sig_m]. 
	weight: 1 for weighted version, others for fitting directly without weigths
*/

void poly_fit_2d(const double *x, const double *y, const double *fxy, const int data_num, const int order, double *coeffs);
/* involves the method cov_matrix_2d().
	fit 2d polynomial f(x,y) to n order for background removing.
	
	f(x,y) = a1 + a2*x + a3*y + a4*x^2 + a5*x*y + a6*y^2 .....
	
	!!! one should be very careful with problem of the precision that comes from the larger power gap between
	!!! the maximum and minimum in the matrix A (the square matrix from least square). So, one should shift and rescale the "x" and "y" before fiiting.
	!!! if the coordinates have been shifted and scaled, the orignal f(x,y) should be calculated at new (x,y)

	x(y) : array, coordinates
	fxy: array, the measured value at (x,y)
	data_num : int, the length of data
	order : the highest order of the target polynomial
	coeffs: array, the results, stores the target parameters, with length (order + 1)*(order + 2) / 2
*/

double fval_at_xy(const double x, const double y, const int order, const double *coeffs );
/* calculate the f(x,y) when given the "order" and "coefficients".
	the polynomial f(x,y) = a1 + a2*x + a3*y + a4*x^2 + a5*x*y + a6*y^2 ..... 
	coeffs:array, [a1,a2,a3...]
*/

void background_remove(const double *arr, const int size_x, const int size_y);
/* fit the stand deviation of background noise of the a chip
*/

void cov_martix_2d(const double *x, const double *y, const double *fxy, const int data_num, const int order, double *cov_matrix, double *f_vertor);
/*	 solve the matrix equation A*P = F which comes from the least square method ( not the eaquations A*X = Y)
	A is the target covariance matrix which will be obtained by this method.
	P is the vector of coefficients ( this method has nothing to do with it)
	F is the target vector which contains the terms like f(x,y)*x^n*y^m and will be obtain by this method.
	
	!!! one should be very careful with problem of the precision that comes from the larger power gap between
	!!! the maximum and minimum in the matrix A. So, one should shift and rescale the "x" and "y" before fiiting.

	f(x,y) = a1 + a2*x + a3*y + a4*x^2 + a5*x*y + a6*y^2 .....

	x(y) : array, coordinates
	fxy: array, the measured value at (x,y)
	data_num : int, the length of data
	order : the highest order of the target polynomial
	cov_matrix : array, with length of  ((order + 1)*(order + 2) / 2)^2, to store the matrix "A" of A*P = F
	f_vector: array, the right of the A*P = F	
*/

void sum_arr(const double *arr, const int size, const int start_t, const int end_t, double &total);
/* sum the array from "start_t" to "end_t"(excluded) and assign to the "total"
	"total" will be set to be zero in the begin of the method for safety.
*/

void arr_pow(const double *arr, double *arr_out, const int size, const int alpha, const int beta, const double power);
/* arr_out = (alpha*arr+beta )^power
*/

void arr_rescale(double *x, const double dx, const double scale, const int num);
/* rescale the array. x= scale*(x+dx)
	"num" is the length of array
*/

void matrix_product(const double*arr_left, const int size_1, const int size_2, const int size_3, const double *arr_right, double *result);
/* C = A*B, calculated by GSL
	arr_left: array, size_1 x size_2
	arr_right: array, size_2 x size_3
	result: array, size_1 x size_3, the result
*/

void matrix_inv(const double *arr, const int size, double *arr_inv);
/* square matrix
*/

/********************************************************************************************************************************************/
/* general methods */
/********************************************************************************************************************************************/

void show_arr(const double*arr, const int size_1, const int size_2);
/* print the elements on the screen
*/
void initialize_para(para *paras);
/* set the "gal_" parameters zero */

void set_bin(const double *data, const int data_num, double * bins, const int bin_num, const double max_scale); //checked
void set_bin(const float *data, const int data_num, float * bins, const int bin_num, const float max_scale);//checked
void set_bin(const int *data, const int data_num, int * bins, const int bin_num, const int max_scale);//checked
/* operate on the copy of data, involving sort_arr().
	the length of bins is bin_num+1
	data: array
	data_num: length of data
	bins: array with length = bin_num+1
	max_scale: the scale (times the outer boundary of bins) of the bins boundary,
					  to make it big enough to contain all the data, even when the data
					  have been shifted.
*/

void histogram(const double *data, const double *bins, int *num, const int data_num, const int bin_num);//checked
void histogram(const float *data, const float *bins, int *num, const int data_num, const int bin_num);//checked
void histogram(const int *data, const  int *bins, int *num, const  int data_num, const  int bin_num);//checked
void histogram2d(const double *data_y, const double*data_x, const double *bin_y, const double *bin_x, int *num, const int data_num, const int ybin_num, const  int xbin_num);
void histogram2d(const float *data_y, const float*data_x, const float *bin_y, const float *bin_x, int *num, const int data_num, const  int ybin_num, const int xbin_num);
void histogram2d(const int *data_y, const int*data_x, const int *bin_y, const int *bin_x, int *num, const int data_num, const int ybin_num, const int xbin_num);
/* data (data_x, data_y): array
	bins (bin_x, bin_y): array, boundary of bins with length=bin_num+1
	num: 1d-array, the counts of the number of the data fall into each bin,
			 for the 2d histogram, the layout of the num-array is the same as the usual 2d array
			 |(0,0)     (x,0)       |
			 |(y,0)					|
			 |							|
	bin_num, y(x)bin_num: array
*/

void sort_arr(double *arr, int size, int order);//checked
void sort_arr(float *arr, int size, int order);//checked
void sort_arr(int *arr, int size, int order);//checked
/* sort the double array according to the "order", 
    order =1 for ascend, else for descend
*/

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

