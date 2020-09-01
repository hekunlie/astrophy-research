// please compile it with C++11 standard(2011)

#ifndef FQLIB_H
#define FQLIB_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<stdio.h>
#include<string>
#include<cstring>
#include<ctime>
#include<cmath>

#include<stdlib.h>

#include<algorithm> // sort(), std::max()
#include<functional> // std::less, std::greater..
#include<ciso646> // for "and, not, or, ..."


#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include<gsl/gsl_cblas.h>
#include<gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_histogram.h>



// something relates to the shear measurement and source detection
struct fq_paras
{
	int psf_size, psf_px, psf_py;
	double psf_peak, psf_hlr, psf_flux, psf_fluxsq, psf_noise_sig, psf_pow_thresh = 0.0001;

	int gal_size, gal_hsize, gal_px, gal_py;
	double gal_peak, gal_hlr, gal_flux, gal_hflux, gal_fluxsq, gal_total_flux,gal_effective_radius;
	double gal_flux2, gal_flux_alt, gal_snr, gal_osnr, gal_noise_sig, gal_flux2_new;
	double gal_size_ext[5];
	double gal_flux_ext[5];
	double gal_flux2_ext[5];

	double gal_e1, gal_e2;
	double n1, n2, dn, du, dv, dp1, dp2;
	double t1, t2, t3, t4;


	/*parameters for detection which should be initialized before */
	int stamp_size; /* the stamp size,  for get_radius() */
	double stamp_cent;
	int img_x, img_y; /* the size of chip image for the 'source_detector()' and galaxy_finder()*/
	int area_thresh=6; /* the minimun pixels for a detection */
	double detect_thresh; /* the threshold of pixel value of source */
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
};


struct fq_paras_float
{
	int psf_size, psf_px, psf_py;
	float psf_peak, psf_hlr, psf_flux, psf_fluxsq, psf_noise_sig, psf_pow_thresh = 0.0001;

	int gal_size, gal_hsize, gal_px, gal_py;
	float gal_peak, gal_hlr, gal_flux, gal_hflux, gal_fluxsq, gal_total_flux,gal_effective_radius;
	float gal_flux2, gal_flux_alt, gal_snr, gal_osnr, gal_noise_sig, gal_flux2_new;
	float gal_size_ext[5];
	float gal_flux_ext[5];
	float gal_flux2_ext[5];

	float gal_e1, gal_e2;
	float n1, n2, dn, du, dv, dp1, dp2;
	float t1, t2, t3, t4;


	/*parameters for detection which should be initialized before */
	int stamp_size; /* the stamp size,  for get_radius() */
	float stamp_cent;
	int img_x, img_y; /* the size of chip image for the 'source_detector()' and galaxy_finder()*/
	int area_thresh=6; /* the minimun pixels for a detection */
	float detect_thresh; /* the threshold of pixel value of source */
	float noise_sig;
	int max_source = 20; /* the maximum of sources allowed in each chip, changeable */
	float max_distance= 8.;/* the max distance of peak away from the center of the source candidate */


	/* hyper_fit_5 matrix elements of order 2 of xy polynomials */
	/* this is the final matrix and the data value is the only things needed */
	float fit_matrix[6][20] =
	{
		{ -0.0530303,0.0113636,-0.0530303,-0.0530303,0.1401515,0.2045455,0.1401515,-0.0530303,0.0113636,0.2045455,0.2045455,0.0113636,-0.0530303,0.1401515,0.2045455,0.1401515,-0.0530303,-0.0530303,0.0113636,-0.0530303 },
		{ -0.0294118,0.0,0.0294118,-0.0588235,-0.0294118,0.0,0.0294118,0.0588235,-0.0588235,-0.0294118,0.0294118,0.0588235,-0.0588235,-0.0294118,0.0,0.0294118,0.0588235,-0.0294118,0.0,0.0294118 },
		{ -0.0588235,-0.0588235,-0.0588235,-0.0294118,-0.0294118,-0.0294118,-0.0294118,-0.0294118,0.0,0.0,0.0,0.0,0.0294118,0.0294118,0.0294118,0.0294118,0.0294118,0.0588235,0.0588235,0.0588235 },
		{ 0.0088745,-0.0172078,0.0088745,0.0517316,-0.0265152,-0.0525974,-0.0265152,0.0517316,0.0399351,-0.0383117,-0.0383117,0.0399351,0.0517316,-0.0265152,-0.0525974,-0.0265152,0.0517316,0.0088745,-0.0172078,0.0088745 },
		{ 0.0555556,0.0,-0.0555556,0.0555556,0.0277778,0.0,-0.0277778,-0.0555556,0.0,0.0,0.0,0.0,-0.0555556,-0.0277778,0.0,0.0277778,0.0555556,-0.0555556,0.0,0.0555556 },
		{ 0.0517316,0.0399351,0.0517316,0.0088745,-0.0265152,-0.0383117,-0.0265152,0.0088745,-0.0172078,-0.0525974,-0.0525974,-0.0172078,0.0088745,-0.0265152,-0.0383117,-0.0265152,0.0088745,0.0517316,0.0399351,0.0517316 },
	};

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
const double DEG2RAD = 1./180*Pi;
const double C_0_hat = 2.99792458; // 10^8
const double C_0 = 2.99792458*1.e8;// speed of light m/s

const double H_0_hat = 0.70;
const double H_0 = 70; //  Km/s/Mpc

const double G_0_hat= 6.6740831; //  10^{-11}
const double G_0 = 6.6740831*1.e-11; //  m^3 s^{-2} Kg^{-1}

const double One_Light_Year_hat = 9.4607304725808;// 10^15
const double One_Light_Year = 9.4607304725808*1.e15;// meter

const double One_Mpc_hat = 3.085677581; // 10^22
const double One_Mpc = 3.085677581*1.e22;// meter

const double M_sun_hat = 1.9885;
const double M_sun = 1.9885*1.e30;//Kg

extern const gsl_rng_type *T0,*T1, *T2, *T3;
extern gsl_rng *rng0, *rng1, *rng2, *rng3; 
extern int GSL_SETUP_LABEL;
extern std::ofstream loggers;



void stack(double *big_arr, const double *stamp, const int tag, const int size, const int row, const int col);//checked
void stack(float *big_arr, const float *stamp, const int tag, const int size, const int row, const int col);//checked
void stack(int *big_arr, const int *stamp, const int tag, const int size, const int row, const int col);//checked
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

void segment(const double *big_arr, double *stamp, const int tag, const int size, const int row, const int col);//checked
void segment(const float *big_arr, float *stamp, const int tag, const int size, const int row, const int col);//checked
void segment(const int *big_arr, int *stamp, const int tag, const int size, const int row, const int col);//checked
/* to cut the specific square area in the big_arr
	see the annotation of stack()	
*/

/********************************************************************************************************************************************/
/* operations on the image */
/********************************************************************************************************************************************/
void create_points(double *point, const int num_p, const double radius, const double step, const gsl_rng *gsl_rand_rng);
void create_points(float *point, const int num_p, const float radius, const float step, const gsl_rng *gsl_rand_rng);

void create_epoints(double *point, const int num_p, const double ellip, const gsl_rng *gsl_rand_rng);
void create_epoints(float *point, const int num_p, const float ellip, const gsl_rng *gsl_rand_rng);


void coord_rotation(const double*xy, const int pts_num, const double theta, double *xy_r);
void coord_rotation(const float*xy, const int pts_num, const float theta, float *xy_r);


void create_psf(double *in_img, const double psf_scale, const int size, const double img_cent, const int psf);
void create_psf(float *in_img, const float psf_scale, const int size, const float img_cent, const int psf);


void create_psf_e(double*in_img, const double psf_scale, const int size, const double img_cent, const double ellip, const double theta, const int psf_type);
void create_psf_e(float*in_img, const float psf_scale, const int size, const float img_cent, const float ellip, const float theta, const int psf_type);


void convolve(const double * points, const int num_p, const double flux_per_pts, const double g1, const double g2,
 				double *in_img, const int size, const double img_cent, const double psf_scale, const int psf_type);
void convolve(const float * points, const int num_p, const float flux_per_pts, const float g1, const float g2, 
				float *in_img, const int size, const float img_cent, const float psf_scale, const int psf_type);


void convolve_e(const double * points, const int num_p, const double flux_per_pts, const double g1, const double g2,
				double *in_img, const int size, const double img_cent, const double psf_scale, const int psf_type, const double ellip, const double theta);
void convolve_e(const float * points, const int num_p, const float flux_per_pts, const float g1, const float g2, 
				float *in_img, const int size, const float img_cent, const float psf_scale, const int psf_type, const float ellip, const float theta);


void deconvolution(const double *gal_pow, const double *psf_pow, double *out_img, const int column, const int row);
void deconvolution(const double *gal_pow, const double*gal_pow_real, const double*gal_pow_imag, const double *psf_pow, const double*psf_pow_real, const double*psf_pow_imag, 
     double *img_out, const int column, const int row);


void pow_spec(const double *in_img, double *out_img_pow, const int column, const int row);
void pow_spec(const float *in_img, float *out_img, const int column, const int row);

void pow_spec(const double *in_img, double *out_img_pow, double *out_img_real, double *out_img_imag, const int column, const int row);
void pow_spec(const float *in_img, float *out_img, float *out_img_real, float *out_img_imag,const int column, const int row);

void pow_spec(const double *in_img, double *out_img_pow, double *out_img_real, double *out_img_imag, double*phase_arg, const int column, const int row);
void pow_spec(const float *in_img, float *out_img_pow, float *out_img_real, float *out_img_imag, float*phase_arg, const int column, const int row);



void get_radius(double *in_img, fq_paras *paras, double scale, int type, double sig_level);
void get_psf_radius(const double *psf_pow, fq_paras*para, const double scale);
void get_psf_radius(const float *psf_pow, fq_paras_float*para, const float scale);
/*measure the size of psf power spectrum for the \beta parameter in the measurement.
	power of k=0 may be not the maximun, be careful!!!! */

void get_quad(const double *img, const int img_size, const double img_cent, const double weight_sigma_sq, double &quad_size);
void get_quad(const float *img, const int img_size, const float img_cent, const float weight_sigma_sq, float &quad_size);
/* calculate the gaussian weighted quadrupole of the image   */
/* weight_sigma_sq: the squared sigma of the gaussian weight */
/* quad size := \Sum {weight*r^2*img} / \Sum {weight*img}      */

void source_detector(const double *source_img, int *source_x, int*source_y, double *source_paras,fq_paras* paras, bool cross, int &detection, std::string &info);
void source_detector(const float *source_img, int *source_x, int*source_y, float *source_paras, fq_paras_float* paras, bool cross, int &detection, std::string &info);
/* operates on the copy,
	if the method finds too many sources ( > para.max_source), the overflows will be ignored.
	source_img: the inputted array in where to find the source galaxies
	source_x, _y:  the array to store the coordinates of sources detected
	source_paras: the array to store the parameters of sources detected,
						   8 elemets for each source, [....,area, peak_y, peak_x, peak_val, half_light_area, total_flux, half_light_flux, flux_sq,...]
	cross: boolean, True for detection on the nearest four pixels, "+", upper, lower, left, right
							False for detection on the nearest eight pixels, "x" and "+"  
	detection: int, the total number of detection 
	info: the information of detection
*/


void galaxy_finder(const double *stamp_arr, int *check_mask, fq_paras *paras, bool cross, int &detect_label, std::string &info);
void galaxy_finder(const float *stamp_arr, int *check_mask, fq_paras_float *paras, bool cross, int &detect_label, std::string &info);
/* to identify the galaxy on each stamp basing on source_detector(), because of many detections on it
	the biggest source which peaks in the central circle with a radius of 6 pixels.	
	return: int, "-1" means no detection

	needs: "stamp_size","max_distance","img_y", "img_x", "detect_thresh", "area_thresh", "max_source"
*/

int edge_extend(int *mask, const int *source_y, const int* source_x, const int source_id, const int source_len, fq_paras *paras, const int iters);
/*	"stamp_size" in the structure fq_paras will be used !!! 

	extend the border of a source galaxy by 1 pixel each time
	the coordinates source_x(y), will be copied to a new array.

	mask: array, on which the source and extendion edge will be labeled, 1 means source pixel
	source_y(x): array of source coordinates from source_detector(), may be more than one source
	source_len : the length of the target source
	source_id : the start of the source coordinates in the source_y(x), because there may be more than one source.
	iters:  iterations, each time will extend the edge by one pixel.
	return: the area of the extended source
*/

void initialize_arr(long *arr, const int length, const long x);
void initialize_arr(double *arr, const int length, const double x);//checked
void initialize_arr(float *arr, const int length, const float x);//checked
void initialize_arr(int *arr, const int length, const int x);//checked
/* set every elements to x*/

void normalize_arr(double *arr, const int size);//checked
void normalize_arr(float *arr, const int size);//checked
/* normalize the PSF power spectrum,
	divide each pixel by the peak 
*/


/********************************************************************************************************************************************/
/* Fourier Quad */
/********************************************************************************************************************************************/
void snr_est(const double *image, fq_paras *paras, int fit);//checked
void snr_est(const float *image, fq_paras_float *paras, int fit);
/* if fit=2 for both flux2 and flux_alt estimations
	else just estimate the flux2
*/

void possion_subtraction(double *image_pow, fq_paras *paras, int edge);//checked
void possion_subtraction(float *image_pow, fq_paras_float *paras, int edge);

void noise_subtraction(double *image_pow, double *noise_pow, fq_paras *paras, const int edge, const int possion);//checked
void noise_subtraction(float *image_pow, float *noise_pow, fq_paras_float *paras, const int edge, const int possion);


void shear_est(double *gal_pow, double *psf_pow, fq_paras *paras);//checked
void shear_est(float *gal_pow, float *psf_pow, fq_paras_float *paras);

void shear_est(const double *gal_pow, const double *gal_pow_real, const double *gal_pow_imag,
               double *psf_pow,const double *psf_pow_real,const double *psf_pow_imag, fq_paras *paras);//checked


void ellip_est(const double *gal_img, const int size, fq_paras*paras);
void ellip_est(const float *gal_img, const int size, fq_paras_float*paras);


void find_block(const pts_info *infos, const double radius_s, const double radius_e, const double *bound_y, const double *bound_x, int *block_mask);//checked
void find_block(const pts_info *infos, const double radius, const double *bound_y, const double *bound_x, int *block_mask);//checked
/* find the target blocks for a specific point for the calculation of correlation function.
	
*/

void block_bound(const double scale, const int ny, const int nx, double *bound_y, double *bound_x);//checked


void chisq_Gbin_1d(const double *mg, const double *mnu, const int data_num, const double *bins, const int bin_num, const double gh, double &result);
void chisq_Gbin_1d(const float *mg, const float *mnu, const int data_num, const float *bins, const int bin_num, const float gh, float &result);
/* calculate the chi square with the input G1, N, U and the guess of shear, gh 
	Fourier Quad shear estimators: G1, G2, N, U, V
	mg: array, G1(2), for g1(2)
	mnu: array, N+U for g1, N-U for g2
	gh: the guess of shear
	bins: the bin for G1(2) , which must be set up (call set_bin()) before , for chi square calculation
	chisq: the chi square with shear guess, gh
*/


void chisq_Gbin_1d(const double *mg,const int data_num, const double *bins, const int bin_num, double &result);
/* 	new for PDF_SYM
*/


void cal_chisq_2d(const double *hist_arr, const int bin_num, double &result);//checked
void cal_chisq_2d(const long *hist_arr, const int bin_num, double &result);//checked
void cal_chisq_2d(const int *hist_arr, const int bin_num, double &result);//checked
/* 2d chi square for correlation calculation
	hist_arr: the 2d histogram of G1, G1~
*/
void cal_chisq_1d(const double *hist_num,  const int bin_num, double &result);
void cal_chisq_1d(const long *hist_num, const int bin_num, double &result);
void cal_chisq_1d(const int *hist_num,  const int bin_num,  double &result);//checked
void cal_chisq_1d(const int *hist_num,  const int bin_num,  float &result);
void cal_chisq_1d(const int *hist_num,  const int bin_num, const int num, double &result);//checked
/* calculate the 1d chi square using the histogram of G1(2) for the shear estimation 
	hist_num: the histogrom of G1(2), the count of G1(2) 
*/


void find_shear_mean(const double *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col, 
						double &gh, double &gh_sig, const int sub_block_num, const double scale=1);
void find_shear_mean(const float *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col, 
						float &gh, float &gh_sig, const int sub_block_num, const float scale=1);

void find_shear_mean(const double *mg, const double *mn, const int data_row, double &gh, double &gh_sig, const int sub_block_num, const double scale=1);
void find_shear_mean(const float *mg, const float *mn, const int data_row, float &gh, float &gh_sig, const int sub_block_num, const float scale=1);
/* 	if the data array is very large, say > 10^7\, then summing it directly may cause numerical problem
	it's better to sum the sub-block and then add the these quantities together.
*/


void find_shear(const double *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col, const int mu_col, const int bin_num, int g_label, 
				double &gh, double &gh_sig, double & chisq_min_fit, double *chi_check, const int chi_fit_num = 20, const int choice=0, 
				const double max_scale=100., const double ini_left = -0.1, const double ini_right = 0.1, const double chi_gap = 40);

void find_shear(const float *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col, const int mu_col, const int bin_num, int g_label, 
				float &gh, float &gh_sig, float & chisq_min_fit, float *chi_check, const int chi_fit_num = 20, const int choice=0, 
				const float max_scale=100., const float ini_left = -0.1, const float ini_right = 0.1, const float chi_gap = 40);
// these two for the data that contain all the estimators in one array
// g_label: 1 for g1, others for g2
void find_shear_iter(const double *mg, const double *mn, const double *mu, const int data_num, const int bin_num, const double *bins, int g_label, const int iters, double *result, 
						const int chi_fit_num=20, const double ini_left= -0.1, const double ini_right= 0.1, const double chi_gap= 40);
void find_shear_iter(const float *mg, const float *mn, const float *mu, const int data_num, const int bin_num, const float *bins, int g_label, const int iters, float *result, 
						const int chi_fit_num=20, const float ini_left= -0.1, const float ini_right= 0.1, const float chi_gap=40);

void find_shear(const double *mg, const double *mn, const double *mu, const int data_num, const int bin_num, const double*bins, const int g_label, 
				double &gh, double &gh_sig, double & chisq_min_fit, double *chi_check, const int chi_fit_num = 20, const double ini_left = -0.1, 
				const double ini_right = 0.1, const double chi_gap = 40);

void find_shear(const float *mg, const float *mn, const float *mu, const int data_num, const int bin_num, const float*bins, const int g_label, 
				float &gh, float &gh_sig, float & chisq_min_fit, float *chi_check, const int chi_fit_num = 20, const float ini_left = -0.1, 
				const float ini_right = 0.1, const float chi_gap = 40);
// checked
/* estimate shear and sigma using dichotomy 
	Fourier Quad shear estimators: G1, G2, N, U, V

	find the minimum range of g, then call "fit_shear()" to fit a quadratic function

	mg: array, G1(2), for g1(2)
	mnu: array, N+U for g1, N-U for g2
	data_num: data number
	bin_num: must be even number, bin number, >= 4
	gh (gh_sig): the result, g and sigma of g
	chi_check: store the chi squares for checking, length = 2*chi_fit_num.
	chi_fit_num: the number of point in the final interval of signal for fitting.
	choice: if > 0, "randomly" choose a sub-sample to set up the bin for shear estimation to save time
	max_scale: by which the boundary of the bins (the first and last elements) will be multiplied to include all data points.
	ini_left: the initial guess of shear of the left end
	ini_right: the initial guess of shear of the right end
	chi_gap: the difference between left- (right-) chi square and  middle chi square,  >= 40 recommended
*/

void fit_shear(const double *shear, const double *chisq, const int num, double &gh, double &gh_sig, double &chisq_min_fit, const double chi_gap = 40);// checked
void fit_shear(const float *shear, const float *chisq, const int num, float &gh, float &gh_sig, float &chisq_min_fit, const float chi_gap = 40);
/* fitting a quadratic function to estimate shear 
	
	shear: array, the shears [start, end] for fitting, the X
	chisq: array, chi square corresponding to the points in "shear" array, the Y
	num: the number of point
	gh (gh_sig): the result, g and sigma of g
	chi_gap: the difference between left- (right-) chi square and  middle chi square,  >= 40 recommended
				  if < 0, all the points will be used for fitting.
*/

void find_shear_NU(const double *mg, const double *mn, const double *mu, const int data_num, const int bin_num, const double *bins, 
					int NU_label, double fix_g, double &gh, double &gh_sig, double &chisq_min_fit, double *chi_check, 
					const int chi_fit_num, const double ini_left, const double ini_right, const double chi_gap);

void find_shear_NU(const float *mg, const float *mn, const float *mu, const int data_num, const int bin_num, const float *bins,
					 int NU_label, float fix_g, float &gh, float &gh_sig, float &chisq_min_fit, float *chi_check, 
					const int chi_fit_num, const float ini_left, const float ini_right, const float chi_gap);
		
void fourier_hist_NU(const double *mg, const double *mn,const double *mu, const int data_row, const double gN, const double gU, 
					const double *bins, int *num_in_bin, const int bin_num);
void fourier_hist_NU(const float *mg, const float *mn,const float *mu, const int data_row, const float gN, const float gU, 
					const float *bins, int *num_in_bin, const int bin_num);

void fourier_hist(const double *mg, const double *mn,const double *mu, const int data_row, const double gh, const int g_label, const double *bins, int *num_in_bin, const int bin_num);
void fourier_hist(const float *mg, const float *mn,const float *mu, const int data_row, const float gh, const int g_label, const float *bins, int *num_in_bin, const int bin_num);

void fourier_hist(const double *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col, const int mu_col, const double gh, 
					const int g_label, const double *bins, int *num_in_bin, const int bin_num);
void fourier_hist(const float *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col, const int mu_col, const float gh, 
					const int g_label, const float *bins, int *num_in_bin, const int bin_num);

void estimator_rotation(const double theta,const double mg1, const double mg2, const double mn, const double mu, const double mv, double *output);
void estimator_rotation(const float theta,const float mg1, const float mg2, const float mn, const float mu, const float mv, float *output);

void image_rotation(const double *image_in, double *image_out, const int size);
void image_rotation(const float *image_in, float *image_out, const int size);
// 90 degrees rotation

/********************************************************************************************************************************************/
/* cosmology */
/********************************************************************************************************************************************/
void com_distance(const double low_z, const double high_z, const double omg_m, const double omg_lam, double &result, const double precision_thresh = 1.e-4, const bool integ_only=false);//checked
/* calculate the comoving distance [ Mpc/h ].																																	*/
/*There are only matter and dark energy in the universe																													 */
/* the precision_thresh is set to be 1.e-8 by default to reduce the difference with astropy.cosmology module (<0.001 Mpc)			*/
/* integ_only: false calculate the comoving distance,																														*/
/*			   true calculate the integrate part of the comoving distance																						*/												

void log_bin(const double start, const double end, const int num, double * bins);
/* the logrithmical bin, including the start and end point
	num: the number of bin borders
*/
void linspace(const double start, const double end, const int num, double *bins);
/* numpy.linspace(st, ed, num), include the end point*/



/********************************************************************************************************************************************/
/* fitting */
/********************************************************************************************************************************************/
void image_convole(double *image_in, double *image_out, const int img_size, const double *kernel, const int kernel_size);
void image_convole(float *image_in, float*image_out, const int img_size, const float *kernel, const int kernel_size);

void smooth(double *image, const double *coeffs, fq_paras *paras);//checked
/* smooth all the region */
void smooth(double *image, const double *psf_pow, const double *coeffs, fq_paras *paras);//checked
/* the image will be replaced by the smoothed one. 
	the psf_pow and the fq_paras->psf_thresh_pow are the mask and  threshold to label the region where to be smoothed
	to fit the curve: a1 + a2*x +a3*y + a4*x^2 +a5*x*y + a6*y^2  
*/

void smooth_real(double*image, const double *coeffs, fq_paras *paras);
/* smooth the image by fitting a polynomial */

void hyperfit_5(const double *data, double*fit_para, fq_paras *paras);//checked
void hyperfit_5(const float *data, float*fit_para, fq_paras_float *paras);


void poly_fit_1d(const double *x, const double *fx, const int data_num, const int order, double *coeffs);//checked
void poly_fit_1d(const float *x, const float *fx, const int data_num, const int order, float *coeffs);
/* fit y = a1 +a2*x +a3*x^2 + a4*x^3 ... */

void poly_fit_1d(const double *x, const double *fx, const double *fx_err, const int data_num, double *coeffs, int weight);// checked
void poly_fit_1d(const float *x, const float *fx, const float *fx_err, const int data_num, float *coeffs, int weight);
/* polynomial fitting by GSL
	FX = c + m*X
	x: array, coordinates
	fx: array, measured values
	fx_err: array, uncertainty (\sigma), of which the reciprocal will be the weights of data points.
	data_num: the number of data points
	coeffs: array[4], [c, sig_c, m, sig_m]. 
	weight: 1 for weighted version, others for fitting directly without weights
*/

void poly_fit_2d(const double *x, const double *y, const double *fxy, const int data_num, const int order, double *coeffs);//checked
/* involves the method cov_matrix_2d().
	fit 2d polynomial f(x,y) to n order for background removing.
	
	f(x,y) = a1 + a2*x + a3*y + a4*x^2 + a5*x*y + a6*y^2 .....
	
	!!! one should be very careful with problem of the precision that comes from the larger power gap between
	!!! the maximum and minimum in the matrix A (the square matrix from least square). So, one should shift and rescale the "x" and "y" before fitting.
	!!! if the coordinates have been shifted and scaled, the original f(x,y) should be calculated at new (x,y)

	x(y) : array, coordinates
	fxy: array, the measured value at (x,y)
	data_num : int, the length of data
	order : the highest order of the target polynomial
	coeffs: array, the results, stores the target parameters, with length (order + 1)*(order + 2) / 2
*/

double fval_at_xy(const double x, const double y, const int order, const double *coeffs );//checked
/* calculate the f(x,y) when given the "order" and "coefficients".
	the polynomial f(x,y) = a1 + a2*x + a3*y + a4*x^2 + a5*x*y + a6*y^2 ..... 
	coeffs:array, [a1,a2,a3...]
*/

//void background_remove(const double *arr, const int size_x, const int size_y);
/* fit the stand deviation of background noise of the a chip
*/

void cov_matrix_1d(const double *x, const double *fx, const int data_num, const int order, double *cov_matrix, double *f_vector);//checked
void cov_matrix_1d(const float *x, const float *fx, const int data_num, const int order, float *cov_matrix, float *f_vector);
/* calculate the matrix on the left and right of the equation from the least square method to "order" */

void cov_matrix_2d(const double *x, const double *y, const double *fxy, const int data_num, const int order, double *cov_matrix, double *f_vertor);//checked
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

void sum_arr(const double *arr, const int size, const int start_t, const int end_t, double &total);//checked
void sum_arr(const long *arr, const int size, const int start_t, const int end_t, long &total);//checked
void sum_arr(const int *arr, const int size, const int start_t, const int end_t, int &total);//checked
/* sum the array from "start_t" to "end_t"(excluded) and assign to the "total"
	"total" will be set to be zero in the begin of the method for safety.
*/

void arr_pow(const double *arr, double *arr_out, const int size, const int alpha, const int beta, const double power);//checked
/* arr_out = (alpha*arr+beta )^power
*/

void arr_rescale(double *x, const double dx, const double scale, const int num);//checked
/* rescale the array. x= scale*(x+dx)
	"num" is the length of array
*/

void matrix_product(const double*arr_left, const int size_1, const int size_2, const int size_3, const double *arr_right, double *result);//checked
/* C = A*B, calculated by GSL
	arr_left: array, size_1 x size_2
	arr_right: array, size_2 x size_3
	result: array, size_1 x size_3, the result
*/

void matrix_inv(const double *arr, const int size, double *arr_inv);//checked
/* square matrix
*/

/********************************************************************************************************************************************/
/* general methods */
/********************************************************************************************************************************************/
void separation(const double RA1, const double DEC1, const double RA2, const double DEC2, double &sep_radian);// checked
void separation(const float RA1, const float DEC1, const float RA2, const float DEC2, float &sep_radian);
/* calculate the separation angle between two points on the sphere */
/* RA & DEC in unit of degree */
/* https://en.wikipedia.org/wiki/Great-circle_distance */

void find_near(const double *arr, const double tar_val, const int arr_len, int & label);//checked
/* find the label of the element in "arr" which is the nearest to "tar_val" */

void check_buffer(double *target_arr, double *buffer, const int start_t, const int buffer_size, int & count, int count_line);
/* if the count > count_line, the data in buffer will be added to the target_arr and be cleared then.
	call it anywhere you want, if count>count_line, it will work.

	it is designed for the large number.
	it may be useless if one add a small numer to a very large number because of the finite significant digits 

	target_arr: array, a big container
	buffer: array, a temp container
	start_t: int, the start element in the targer_arr where the buffer will be added to
	buffer_size: int, the size of buffer
	count: int, the count, it will be 
	count_line: int, the upper bound for count
*/

void task_alloc(const int *label_list, const int total_task_num, const int portion,  const int portion_label, int *allocated_list );//checked
/* 
	distribute the tasks
	-1 is used to label the end of the tasks list
	please keep the lengths of "label_list" and "allocated_list" the same

	label_list: array, contains the labels of each task, 1,2,3... or something, non-negative,
	total_task_num: int, number of total tasks
	portion: int, how many portions the total task will be divided into
	portion_label: int, "0" means the first part, which portion to be returned
	allocated_list: array, the labels of the returned tasks	
*/

void task_alloc(const int total_task_num, const int portion, const int my_part_id, int &my_start, int &my_end, int *task_count);
void task_alloc(const int total_task_num, const int portion, const int my_part_id, int &my_start, int &my_end, int *task_count, int *entry_for_gather);
/* distribute the tasks */

void show_arr(const double*arr, const int rows, const int cols);//checked
void show_arr(const long*arr, const int rows, const int cols);//checked
void show_arr(const int*arr, const int rows, const int cols);//checked
void show_arr(const float*arr, const int rows, const int cols);//checked
/* print the elements on the screen
*/
void initialize_para(fq_paras *paras);
void initialize_para(fq_paras_float *paras);

/* set the "gal_" parameters zero */


void set_bin(double *data, const int data_num, double * bins, const int bin_num, const double max_scale);
void set_bin(float *data, const int data_num, float * bins, const int bin_num, const float max_scale);
void set_bin(const double *data, const int data_num, const int bin_num, double * bins, const double max_scale, int choice); //checked
void set_bin(const float *data, const int data_num, const int bin_num, float * bins, const float max_scale, int choice);//checked
void set_bin(const int *data, const int data_num, int * bins, const int bin_num, const int max_scale, int choice);//checked
/* operate on the copy of data, involving sort_arr().
	the length of bins is bin_num+1
	data: array
	data_num: length of data
	bins: array with length = bin_num+1
	max_scale: the scale (times the outer boundary of bins) of the bins boundary,
					  to make it big enough to contain all the data, even when the data
					  have been shifted.
	choice: if > 0, "randomly" choose a sub-sample to set up the bin for shear estimation to save time
*/

void histogram(const double *data, const double *bins, int *num_in_bin, const int data_num, const int bin_num);//checked
void histogram(const double *data, const double *bins, long *num_in_bin, const int data_num, const int bin_num);
void histogram(const float *data, const float *bins, int *num_in_bin, const int data_num, const int bin_num);//checked
void histogram(const float *data, const float *bins, long *num_in_bin, const int data_num, const int bin_num);
void histogram(const int *data, const  int *bins, int *num_in_bin, const  int data_num, const  int bin_num);//checked
void histogram(const int *data, const  int *bins, long *num_in_bin, const  int data_num, const  int bin_num);

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
void histogram_s(const double data, const double *data_bins, int bin_num, int &bin_label);
void histogram2d_s(const double data_y, const double data_x, const double *bin_y, const double *bin_x, const int ybin_num, const  int xbin_num, int &bin_label);


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


/************************************* random *********************************************************/

/**************************************** GSL *********************************************************/
// if use GSL, you must setup the gsl_env.
// then initialize the gsl_rng, gsl_initialize(gsl_rng),
// it will also call gsl_setup() if the gsl_env has not been setup
// finally, you must free the gsl_rng by calling gsl_free(gsl_rng)
void gsl_setup();
void gsl_initialize(int seed, const int rng_label);
// three gsl_rng_type, T1, T2 and T3, and three gsl_rng, rng1, rng2 and rng3
// have been created.
// call gsl_initialize() to initialize they.
// type: 0 T0 & rng0, 1 T1 & rng1, 2 T2 & rng2, 3 T3 & rng3
void gsl_free(const int rng_label);

void addnoise(double *image, const int pixel_num, const double sigma, const gsl_rng *gsl_rand_rng);//checked
void addnoise(float *image, const int pixel_num, const float sigma, const gsl_rng *gsl_rand_rng);//checked
/* add Gaussian noise to an array */
void noise2pow(double *pow_img, const int size, const double sigma, gsl_rng *gsl_rand_rng);

void rand_gauss(const double sigma, const double mean, double &rand_n, const gsl_rng *gsl_rand_rng);//checked
/* return a double from the normal distribution with sigma and mean. 
*/

void rand_multi_gauss(const double*cov, const double *mu, const int num, double *result, const gsl_rng *gsl_rand_rng);//checked
/* calling the gsl_ran_multivariate_gaussian() to generate the k-dimensional multivariate Gaussian 

	cov: array, the covariance matrix
	mu: array, the means
	num: int, the dimensions k
	result: array, the k numbers
*/

void rand_uniform(const double start, const double end, double &rand_n, const gsl_rng *gsl_rand_rng);
/* return a double, [ start, end ), with a unifrom distribution.
*/

void rand_shuffle(double *seq, int length, const gsl_rng *gsl_rand_rng);//checked
void rand_shuffle(float *seq, int length, const gsl_rng *gsl_rand_rng);//checked
void rand_shuffle(int *seq, int length, const gsl_rng *gsl_rand_rng);//checked
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


#endif // !FQLIB_H

