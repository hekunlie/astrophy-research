﻿#include "FQlib.h"

const gsl_rng_type *T;
gsl_rng *rng;
std::ofstream loggers;

char buffer[1000], exception_name[50];
/********************************************************************************************************************************************/
/* file reading and writting*/
/********************************************************************************************************************************************/
void write_log(char*filename, char *inform)
{
	char time_now[40];
	get_time(time_now, 40);
	loggers.open(filename, std::ios::out | std::ios::app);
	loggers << time_now << " ---- " << inform << std::endl;
	loggers.close();
}


void read_para(const std::string path, const std::string name, int &para)
{
	std::ifstream infile;
	std::string str, str1, str2, str3;
	std::stringstream strs;
	infile.open(path);
	int f = 0;
	while (!infile.eof())
	{
		str.clear();
		str1.clear();
		str2.clear();
		str3.clear();
		strs.clear();

		getline(infile, str);
		strs << str;
		strs >> str1 >> str2 >> str3;

		if (str1 == name)
		{
			para = std::stoi(str3);
			f = 1;
			break;
		}
	}
	infile.close();
	if (f == 0)
	{	
		str.clear();
		str = name + "can not be found!!";
		std::cout << str << std::endl;
		exit(0);
	}
}

void read_para(const std::string path, const std::string name, double &para)
{
	std::ifstream infile;
	std::string str, str1, str2, str3;
	std::stringstream strs;
	int f = 0;
	infile.open(path);
	while (!infile.eof())
	{
		str.clear();
		str1.clear();
		str2.clear();
		str3.clear();
		strs.clear();

		getline(infile, str);
		strs << str;
		strs >> str1 >> str2 >> str3;

		if (str1 == name)
		{
			para = std::stod(str3);
			f = 1;
			break;
		}
	}
	infile.close();
	if (f == 0)
	{
		str.clear();
		str = name + "can not be found!!";
		std::cout << str << std::endl;
		exit(0);
	}
}

void read_para(const std::string path, const std::string name, float &para)
{
	std::ifstream infile;
	std::string str, str1, str2, str3;
	std::stringstream strs;
	int f = 0;
	infile.open(path);
	while (!infile.eof())
	{
		str.clear();
		str1.clear();
		str2.clear();
		str3.clear();
		strs.clear();

		getline(infile, str);
		strs << str;
		strs >> str1 >> str2 >> str3;

		if (str1 == name)
		{
			para = std::stof(str3);
			f = 1;
			break;
		}
	}
	infile.close();
	if (f == 0)
	{
		str.clear();
		str = name + "can not be found!!";
		std::cout << str << std::endl;
		exit(0);
	}
}


void read_text(const std::string path, double *arr, const int read_lines)
{
	std::ifstream infile;
	int i = 0;
	infile.open(path);
	while (i<read_lines && infile>>arr[i])
	{
		i++;
	}	
	infile.close();
}

void read_text(const std::string path, float *arr, const int read_lines)
{
	std::ifstream infile;
	int i = 0;
	infile.open(path);
	while (i<read_lines && infile >> arr[i])
	{
		i++;
	}
	infile.close();
}

void read_text(const std::string path, int *arr, const int read_lines)
{
	std::ifstream infile;
	int i = 0;
	infile.open(path);
	while (i<read_lines && infile >> arr[i])
	{
		i++;
	}
	infile.close();
}


void read_h5(const char *filename, const char *set_name, double *arr)
{
	hid_t file_id;
	herr_t status;
	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	hid_t dataset_id;
	dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
}

void read_h5(const char *filename, const char *set_name, float *arr)
{
	hid_t file_id;
	herr_t status;
	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	hid_t dataset_id;
	dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
}

void read_h5(const char *filename, const char *set_name, int *arr)
{
	hid_t file_id;
	herr_t status;
	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	hid_t dataset_id;
	dataset_id = H5Dopen(file_id, set_name, H5P_DEFAULT);
	status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
}

void write_h5(const char *filename, const char *set_name, const double*arr, const int row, const int column)
{
	hid_t file_id;
	herr_t status;
	remove(filename);
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);	
	unsigned rank = 2;
	hsize_t dims[2];
	dims[0] = row;
	dims[1] = column;
	hid_t dataspace_id;  
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	hid_t dataset_id;    
	dataset_id = H5Dcreate(file_id, set_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Fclose(file_id);
 }

void write_h5(const char *filename, const char *set_name, const float*arr, const int row, const int column)
{
	hid_t file_id;
	herr_t status;
	remove(filename);
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	unsigned rank = 2;
	hsize_t dims[2];
	dims[0] = row;
	dims[1] = column;
	hid_t dataspace_id;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	hid_t dataset_id;
	dataset_id = H5Dcreate(file_id, set_name, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Fclose(file_id);
}

void write_h5(const char *filename, const char *set_name, const int*arr, const int row, const int column)
{
	hid_t file_id;
	herr_t status;
	remove(filename);
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	unsigned rank = 2;
	hsize_t dims[2];
	dims[0] = row;
	dims[1] = column;
	hid_t dataspace_id;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	hid_t dataset_id;
	dataset_id = H5Dcreate(file_id, set_name, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Fclose(file_id);
}


void read_fits(const char *filename, double *arr)
{
	fitsfile *fptr;															/* pointer to the FITS file, defined in fitsio.h */
	int status = 0, nfound, anynull;
	long naxes[2], fpixel = 1, nbuffer, npixels, ii;
	double datamin, datamax, nullval = 0;

	fits_open_file(&fptr, filename, READONLY, &status);
	fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);		/* read the NAXIS1 and NAXIS2 keyword to get image size */
	npixels = naxes[0] * naxes[1];	/* number of pixels in the image, python_arr[naxes[1], naxes[0]] */
	double *buffer = new double[npixels];										/* create a new array */
	fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval, buffer, &anynull, &status);
	for (ii = 0; ii < npixels; ii++) arr[ii] = buffer[ii];
	delete[] buffer;														/* (have to) delete the array */
	fits_close_file(fptr, &status);
}

void read_fits(const char *filename, float *arr)
{
	fitsfile *fptr;			/* pointer to the FITS file, defined in fitsio.h */
	int status = 0, nfound, anynull;
	long naxes[2], fpixel = 1, nbuffer, npixels, ii;
	double datamin, datamax, nullval = 0;

	fits_open_file(&fptr, filename, READONLY, &status);
	fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);		/* read the NAXIS1 and NAXIS2 keyword to get image size */
	npixels = naxes[0] * naxes[1];	/* number of pixels in the image, python_arr[naxes[1], naxes[0]] */
	float *buffer = new float[npixels];	/* create a new array */
	fits_read_img(fptr, TFLOAT, fpixel, npixels, &nullval, buffer, &anynull, &status);
	for (ii = 0; ii < npixels; ii++) arr[ii] = buffer[ii];
	delete[] buffer;		/* (have to) delete the array */
	fits_close_file(fptr, &status);
}

void read_fits(const char *filename, int *arr)
{
	fitsfile *fptr;			/* pointer to the FITS file, defined in fitsio.h */
	int status = 0, nfound, anynull;
	long naxes[2], fpixel = 1, nbuffer, npixels, ii;
	double datamin, datamax, nullval = 0;

	fits_open_file(&fptr, filename, READONLY, &status);
	fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);		/* read the NAXIS1 and NAXIS2 keyword to get image size */
	npixels = naxes[0] * naxes[1];	/* number of pixels in the image, python_arr[naxes[1], naxes[0]] */
	int *buffer = new int[npixels];	/* create a new array */
	fits_read_img(fptr, TINT, fpixel, npixels, &nullval, buffer, &anynull, &status);
	for (ii = 0; ii < npixels; ii++) arr[ii] = buffer[ii];
	delete[] buffer;		/* (have to) delete the array */
	fits_close_file(fptr, &status);
}

void write_fits(const char *filename, double *img, const int ysize, const int xsize)
{
	fitsfile *fptr;		/* pointer to the FITS file; defined in fitsio.h */
	int status, ii, jj;
	long fpixel = 1, naxis = 2, nelements, exposure;
	long naxes[2] ;	    /* x, y */
	naxes[0] = xsize;
	naxes[1] = ysize;	 /* the same as numpy array */	
	status = 0;			/* initialize status before calling fitsio routines */
	fits_create_file(&fptr, filename, &status);		/* create new file */
	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);		 /* Create the primary array image (16-bit short integer pixels */
	//fits_update_key(fptr, TDOUBLE, "EXPOSURE", &exposure, "Total Exposure Time", &status);  /* Write a keyword; must pass the ADDRESS of the value */	
	nelements = xsize*ysize;         /* number of pixels to write */
	fits_write_img(fptr, TDOUBLE, fpixel, nelements, img, &status);     /* Write the array of integers to the image */
	fits_close_file(fptr, &status);              /* close the file */
	fits_report_error(stderr, status);      /* print out any error messages */	
}

void write_fits(const char *filename, float *img, const int ysize, const int xsize)
{
	fitsfile *fptr;				/* pointer to the FITS file; defined in fitsio.h */
	int status, ii, jj;
	long fpixel = 1, naxis = 2, nelements, exposure;
	long naxes[2];			/* x, y */
	naxes[0] = xsize;
	naxes[1] = ysize;   	/* the same as numpy array */
	status = 0;		/* initialize status before calling fitsio routines */
	fits_create_file(&fptr, filename, &status);		/* create new file */
	fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);	 /* Create the primary array image (16-bit short integer pixels */
	nelements = xsize * ysize;         /* number of pixels to write */
	fits_write_img(fptr, TFLOAT, fpixel, nelements, img, &status);     /* Write the array of integers to the image */
	fits_close_file(fptr, &status);              /* close the file */
	fits_report_error(stderr, status);      /* print out any error messages */
}

void write_fits(const char *filename, int *img, const int ysize, const int xsize)
{
	fitsfile *fptr;			/* pointer to the FITS file; defined in fitsio.h */
	int status, ii, jj;
	long fpixel = 1, naxis = 2, nelements;
	long naxes[2];		/* x, y */
	naxes[0] = xsize;
	naxes[1] = ysize; 	/* the same as numpy array */
	status = 0;	/* initialize status before calling fitsio routines */
	fits_create_file(&fptr, filename, &status);		/* create new file */
	fits_create_img(fptr, LONG_IMG, naxis, naxes, &status);	 /* Create the primary array image (16-bit short integer pixels */
	nelements = xsize * ysize;         /* number of pixels to write */
	fits_write_img(fptr, TINT, fpixel, nelements, img, &status);     /* Write the array of integers to the image */
	fits_close_file(fptr, &status);              /* close the file */
	fits_report_error(stderr, status);      /* print out any error messages */
}


void stack(double *container, const double *stamp, const int tag, const int size, const int row, const int col)
{
	int i, j, m;
	/* 'container' is the container of 'row'*'col' stamps which is array with 'size'*'size' pixels.
	'tag' is the tag of the stamp in all stamps.
	the (i,j)'th pixel of stamp will loacte at '(tag - tag%col)/col*size*col*size + i*n*size + tag%col8size + j' */
	m = (tag - tag%col)*size*size + tag%col*size;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			container[m + i*col*size + j] = stamp[i*size + j];
		}
	}
}

void stack(float *container, const float *stamp, const int tag, const int size, const int row, const int col)
{
	int i, j, m;
	/* 'container' is the container of 'row'*'col' stamps which is array with 'size'*'size' pixels.
	'tag' is the tag of the stamp in all stamps.
	the (i,j)'th pixel of stamp will loacte at '(tag - tag%col)/col*size*col*size + i*n*size + tag%col8size + j' */
	m = (tag - tag % col)*size*size + tag % col*size;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			container[m + i * col*size + j] = stamp[i*size + j];
		}
	}
}

void stack(int *container, const int *stamp, const int tag, const int size, const int row, const int col)
{
	int i, j, m;
	/* 'container' is the container of 'row'*'col' stamps which is array with 'size'*'size' pixels.
	'tag' is the tag of the stamp in all stamps.
	the (i,j)'th pixel of stamp will loacte at '(tag - tag%col)/col*size*col*size + i*n*size + tag%col8size + j' */
	m = (tag - tag % col)*size*size + tag % col*size;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			container[m + i * col*size + j] = stamp[i*size + j];
		}
	}
}

void segment(const double *chip, double *stamp, const int tag, const int size, const int row, const int col)
{
	int i, j, m;
	m = (tag - tag%col)*size*size + tag%col*size;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			stamp[i*size + j] = chip[m + i*col*size + j];
		}
	}
}

void segment(const float *chip, float *stamp, const int tag, const int size, const int row, const int col)
{
	int i, j, m;
	m = (tag - tag % col)*size*size + tag % col*size;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			stamp[i*size + j] = chip[m + i * col*size + j];
		}
	}
}

void segment(const int *chip, int *stamp, const int tag, const int size, const int row, const int col)
{
	int i, j, m;
	m = (tag - tag % col)*size*size + tag % col*size;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			stamp[i*size + j] = chip[m + i * col*size + j];
		}
	}
}


/********************************************************************************************************************************************/
/* operations on the image */
/********************************************************************************************************************************************/
void create_points(double *point, int num_p, double radius)
{	/* point is the container of the coordinates of the points created which should has 2*num_p elements ,[x..., y....]*/
	int i;
	double x = 0., y = 0., xm = 0., ym = 0., theta, step;

	for (i = 0; i < num_p; i++)
	{
		theta = 2.*Pi*gsl_rng_uniform(rng);
		x += cos(theta);
		y += sin(theta);
		if (x*x + y*y > radius*radius)
		{
			x = cos(theta);
			y = sin(theta);
		}
		point[i] = x;
		point[i + num_p] = y;
		xm += x;
		ym += y;
	}
	//re-move the mass center of the points cluster to (0,0)
	for (i = 0; i < num_p; i++)
	{
		point[i] = point[i] - xm / num_p;
		point[i + num_p] = point[i + num_p] - ym / num_p;
	}
}

void create_epoints(double *point, int num_p, double ellip)
{
	int i = 0;
	double theta, beta, radius, r1, r2, b, y, x;
	radius = 5 + 4 * gsl_rng_uniform(rng);
	b = sqrt((1 - ellip) / (1 + ellip))*radius;
	beta = 2 * Pi*gsl_rng_uniform(rng);
	r1 = cos(beta);
	r2 = sin(beta);
	while (i < num_p)
	{
		x = gsl_ran_gaussian(rng, radius);
		y = gsl_ran_gaussian(rng, b);
		if (x*x / (radius*radius) + y*y / (b*b) <= 1)
		{
			point[i] = x*r1 - y*r2;
			point[i + num_p] = x*r2 + r1*y;
			i++;
		}
	}
}

void create_psf(double*in_img, double scale, int size, int psf)
{
	int i, j;
	double rs, r1, val, flux_g, flux_m, rd;

	flux_g = 1. / (2 * Pi *scale*scale);     /* 1 / sqrt(2*Pi*sig_x^2)/sqrt(2*Pi*sig_x^2) */
	flux_m = 1. / (Pi*scale*scale*(1. - pow(10, -2.5))*0.4); /* 1 / ( Pi*scale^2*( (1 + alpha^2)^(1-beta) - 1) /(1-beta)), where alpha = 3, beta = 3.5 */

	rd = 1. / scale / scale;
	for (i = 0; i < size; i++)
	{
		r1 = (i - size / 2.)*(i - size / 2.)*rd;
		for (j = 0; j < size; j++)
		{
			rs = r1 + (j - size / 2.)*(j - size / 2.)*rd;
			if (psf == 1)  // Gaussian PSF
			{
				if (rs <= 9) in_img[i*size + j] += flux_g*exp(-rs*0.5);
			}
			else              // Moffat PSF
			{
				if (rs <= 9.) in_img[i*size + j] += flux_m*pow(1. + rs, -3.5);
			}
		}
	}
}

void convolve(double *in_img, double * points, double flux, int size, int num_p, int rotate, double scale, double g1, double g2, int psf)
{	 /* will not change the inputted array */
	 /* in_img is the container of the final image,
	 points is the array of points' coordinates,
	 rotate is the radian in units of pi/4,
	 scale is the scale length of PSF,
	 psf=1 means the Gaussian PSF, psf=2 means the Moffat PSF	*/
	int i, j, k, m;
	double r1, r2, n, flux_g, flux_m;
	// |rot1, - rot2 |
	// |rot2,  rot1  |
	double rot1 = cos(rotate*Pi / 4.), rot2 = sin(rotate*Pi / 4.), val, rs, rd;
	rd = 1. / scale / scale;  // scale of PSF	

	double *points_r = new double[num_p * 2];
	/* rotate and shear */
	if (rotate != 0)
	{
		for (i = 0; i < num_p; i++)
		{
			points_r[i] = (1. + g1)*(rot1 * points[i] - rot2*points[i + num_p]) + g2*(rot2 * points[i] + rot1*points[i + num_p]) + size / 2.;
			points_r[i + num_p] = g2*(rot1 * points[i] - rot2*points[i + num_p]) + (1. - g1)*(rot2 * points[i] + rot1*points[i + num_p]) + size / 2.;
		}
	}
	else
	{
		for (i = 0; i < num_p; i++)
		{
			points_r[i] = (1. + g1)* points[i] + g2*points[i + num_p] + size / 2.;
			points_r[i + num_p] = g2*points[i] + (1. - g1)*points[i + num_p] + size / 2.;
		}
	}

	/*  convolve PSF and draw the image */
	flux_g = flux / (2 * Pi *scale*scale);     /* 1 / sqrt(2*Pi*sig_x^2)/sqrt(2*Pi*sig_x^2) */
	flux_m = flux / (Pi*scale*scale*(1. - pow(10, -2.5))*0.4);   /* 1 / ( Pi*scale^2*( (1 + alpha^2)^(1-beta) - 1) /(1-beta)), where alpha = 3, beta = 3.5 */
	for (k = 0; k < num_p; k++)
	{
		for (i = 0; i < size; i++)  /* y coordinate */
		{
			r1 = (i - points_r[k + num_p])*(i - points_r[k + num_p])*rd;
			for (j = 0; j < size; j++) /* x coordinate */
			{
				rs = r1 + (j - points_r[k])*(j - points_r[k])*rd;
				if (psf == 1)  // Gaussian PSF
				{
					if (rs <= 9) in_img[i*size + j] += flux_g*exp(-rs*0.5);
				}
				else				// Moffat PSF
				{
					if (rs <= 9) in_img[i*size + j] += flux_m*pow(1. + rs, -3.5);
				}
			}
		}
	}
	delete[] points_r;
}


void pow_spec(const double *in_img, double *out_img, const int column, const int row)
{   /* will not change the inputted array */
	/* in_img is the inputted array and the out_img is the container of the outputted image */
	fftwl_complex *in, *out;
	fftwl_plan p;
	long int i, j, m, n;

	in = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex) *(row*column));
	for (i = 0; i < (row*column); i++)
	{
		in[i][1] = 0;
		in[i][0] = in_img[i];
	}

	out = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex) *(row*column));
	p = fftwl_plan_dft_2d(row, column, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwl_execute(p);															 /* repeat as needed */

	 /* shift the signal to the center */
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < column; j++)
		{
			m = i*column + j;
			if (i<row / 2)
			{
				if (j < column / 2)
					n = (i + row / 2)*column + j + column / 2;
				else
					n = (i + row / 2)*column + j - column / 2;
			}
			else
			{
				if (j<column / 2)
					n = (i - row / 2)*column + j + column / 2;
				else
					n = (i - row / 2)*column + j - column / 2;
			}

			out_img[n] = out[m][0] * out[m][0] + out[m][1] * out[m][1];
		}
	}
	/* shift the signal to the center */

	fftwl_destroy_plan(p);
	fftwl_free(in);
	fftwl_free(out);
}

void pow_spec(const float *in_img, float *out_img, const int column, const int row)
{   /* will not change the inputted array */
	/* in_img is the inputted array and the out_img is the container of the outputted image */
	fftwf_complex *in, *out;
	fftwf_plan p;
	long int i, j, m, n;

	in = (fftwf_complex*)fftwl_malloc(sizeof(fftwf_complex) *(row*column));
	for (i = 0; i < (row*column); i++)
	{
		in[i][1] = 0;
		in[i][0] = in_img[i];
	}

	out = (fftwf_complex*)fftwl_malloc(sizeof(fftwf_complex) *(row*column));
	p = fftwf_plan_dft_2d(row, column, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(p);															 /* repeat as needed */

	 /* shift the signal to the center */
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < column; j++)
		{
			m = i * column + j;
			if (i < row / 2)
			{
				if (j < column / 2)
					n = (i + row / 2)*column + j + column / 2;
				else
					n = (i + row / 2)*column + j - column / 2;
			}
			else
			{
				if (j < column / 2)
					n = (i - row / 2)*column + j + column / 2;
				else
					n = (i - row / 2)*column + j - column / 2;
			}

			out_img[n] = out[m][0] * out[m][0] + out[m][1] * out[m][1];
		}
	}
	/* shift the signal to the center */

	fftwf_destroy_plan(p);
	fftwf_free(in);
	fftwf_free(out);
}


void get_radius(double *in_img, para *paras, double scale, int type, double sig_level)
{	 /* will not change the inputted array */
	/* the image should be larger than 12*12 */
	/* setting scale = infinity ,  one can obtain the area of the signal */

	int x, y, xp = 0, yp = 0, num0 = 0, num = 1, nump, p = 1, size = paras->stamp_size;
	double max = 0, flux = 0, flux_sq=0.;	
	double *cp_img = new double[size*size]();
	int *col = new int[size*size];
	int *row = new int[size*size];
	/* find the maximun */
	for (y = size / 2 - 5; y < size / 2 +5; y++)
	{
		for (x = size / 2 - 5; x < size / 2 + 5; x++)
		{
			if (in_img[y*size + x] > max)
			{
				max = in_img[y*size + x];
				xp = x;
				yp = y;
			}
		}
	}
	/* copy the image and wrap out the value smaller than the specific one */
	for (x = 0; x < size*size; x++)
	{	
		if (in_img[x] > max / scale && in_img[x] > 1.5 * sig_level)
		{
			cp_img[x] = in_img[x];
		}
	}
	row[0] = yp;
	col[0] = xp;
	flux = max;
	flux_sq = max*max;
	/* the peak coordinates have been put in the coordinate lists */
	cp_img[ yp*size + xp] = 0.;

	while (num0 != num)
	{
		nump = num - num0;
		num0 = p;
		for (x = num0 - nump; x < num0; x++)
		{
			if ((row[x] - 1) > -1 && fabs(cp_img[col[x] + (row[x] - 1) * size]) > 0)
			{
				flux += cp_img[col[x] + (row[x] - 1) * size];
				flux_sq += cp_img[col[x] + (row[x] - 1) * size]* cp_img[col[x] + (row[x] - 1) * size];
				cp_img[col[x] + (row[x] - 1) * size] = 0.;
				row[p] = row[x] - 1;
				col[p] = col[x];
				p++;
			}
			if ((row[x] + 1) < size && fabs(cp_img[col[x] + (row[x] + 1) * size]) > 0 )
			{
				flux += cp_img[col[x] + (row[x] + 1) * size];
				flux_sq += cp_img[col[x] + (row[x] + 1) * size]* cp_img[col[x] + (row[x] + 1) * size];
				cp_img[col[x] + (row[x] + 1) * size] = 0.;
				row[p] = row[x] + 1;
				col[p] = col[x];
				p++;
			}
			if ((col[x] - 1) > -1 && fabs(cp_img[col[x] - 1 + row[x] * size]) > 0 )
			{
				flux += cp_img[col[x] - 1 + row[x] * size];
				flux_sq += cp_img[col[x] - 1 + row[x] * size]* cp_img[col[x] - 1 + row[x] * size];
				cp_img[col[x] - 1 + row[x] * size] = 0.;
				row[p] = row[x];
				col[p] = col[x] - 1;
				p++;
			}
			if ((col[x] + 1) < size && fabs(cp_img[col[x] + 1 + row[x] * size]) > 0 )
			{
				flux += cp_img[col[x] + 1 + row[x] * size];
				flux_sq += cp_img[col[x] + 1 + row[x] * size] * cp_img[col[x] + 1 + row[x] * size];
				cp_img[col[x] + 1 + row[x] * size] = 0.;
				row[p] = row[x];
				col[p] = col[x] + 1;
				p++;
			}
		}
		num = p;
	}
	if (type == 1)
	{	
		paras->psf_peak = max;
		paras->psf_py = yp;
		paras->psf_px = xp;
		paras->psf_hlr = sqrt(p / Pi);
		paras->psf_flux = flux;
	}
	else
	{
		if (max > 1.5 * sig_level && p>4) 
		{
			paras->gal_peak = max;
			paras->gal_py = yp;
			paras->gal_px = xp;

			paras->gal_hlr = sqrt(p / Pi);
			paras->gal_flux = flux;
			paras->gal_size = p;
			paras->gal_fluxsq = flux_sq;
			paras->gal_snr = sqrt(flux_sq) / paras->gal_noise_sig;
			paras->gal_osnr = flux / sqrt(p) / paras->gal_noise_sig;
		}
		else
		{	
			paras->gal_peak = 0;
			paras->gal_py = 0;
			paras->gal_px = 0;

			paras->gal_hlr =0;
			paras->gal_flux = 0;
			paras->gal_size = 0;
			paras->gal_fluxsq = 0;
			paras->gal_snr = 0;
			paras->gal_osnr = 0;
		}

	}
	delete[] col;
	delete[] row;
	delete[] cp_img;
}


void get_psf_radius(const double *psf_pow, para*paras, const double scale)
{
	int x, y, xp = 0, yp = 0, num0 = 0, num = 1, nump, p, size = paras->stamp_size;
	double max = 0;
	double *cp_img = new double[size*size]();
	int *col = new int[size*size];
	int *row = new int[size*size];
	/* find the maximun, but power of k=0 may be not the maximun, be careful!!!! */
	for (y = size / 2 - 5; y < size / 2 + 5; y++)
	{
		for (x = size / 2 - 5; x < size / 2 + 5; x++)
		{
			if (psf_pow[y*size + x] > max)
			{
				max = psf_pow[y*size + x];
				xp = x;
				yp = y;
			}
		}
	}
	paras->psf_pow_thres = max / 10000.;
	/* copy the image and wrap out the value smaller than the specific one */	
	for (x = 0; x < size*size; x++)
	{
		if (psf_pow[x] > max / scale)
		{
			cp_img[x] = psf_pow[x];
		}
	}
	row[0] = yp;
	col[0] = xp;
	p = 1;
	/* the peak coordinates have been put in the coordinate lists */
	cp_img[yp*size + xp] = 0.;

	while (num0 != num)
	{
		nump = num - num0;
		num0 = p;
		for (x = num0 - nump; x < num0; x++)
		{
			if ((row[x] - 1) > -1 && fabs(cp_img[col[x] + (row[x] - 1) * size]) > 0)
			{
				cp_img[col[x] + (row[x] - 1) * size] = 0.;
				row[p] = row[x] - 1;
				col[p] = col[x];
				p++;
			}
			if ((row[x] + 1) < size && fabs(cp_img[col[x] + (row[x] + 1) * size]) > 0)
			{
				cp_img[col[x] + (row[x] + 1) * size] = 0.;
				row[p] = row[x] + 1;
				col[p] = col[x];
				p++;
			}
			if ((col[x] - 1) > -1 && fabs(cp_img[col[x] - 1 + row[x] * size]) > 0)
			{
				cp_img[col[x] - 1 + row[x] * size] = 0.;
				row[p] = row[x];
				col[p] = col[x] - 1;
				p++;
			}
			if ((col[x] + 1) < size && fabs(cp_img[col[x] + 1 + row[x] * size]) > 0)
			{
				cp_img[col[x] + 1 + row[x] * size] = 0.;
				row[p] = row[x];
				col[p] = col[x] + 1;
				p++;
			}
		}
		num = p;
	}
	paras->psf_hlr = sqrt(p / Pi);
	delete[] cp_img;
	delete[] col;
	delete[] row;
}

void get_psf_radius(const float *psf_pow, para*paras, const float scale)
{
	int x, y, xp = 0, yp = 0, num0 = 0, num = 1, nump, p, size = paras->stamp_size;
	float max = 0;
	float *cp_img = new float[size*size]();
	int *col = new int[size*size];
	int *row = new int[size*size];
	/* find the maximun, but power of k=0 may be not the maximun, be careful!!!! */
	for (y = size / 2 - 5; y < size / 2 + 5; y++)
	{
		for (x = size / 2 - 5; x < size / 2 + 5; x++)
		{
			if (psf_pow[y*size + x] > max)
			{
				max = psf_pow[y*size + x];
				xp = x;
				yp = y;
			}
		}
	}
	paras->psf_pow_thres = max / 10000.;
	/* copy the image and wrap out the value smaller than the specific one */	
	for (x = 0; x < size*size; x++)
	{
		if (psf_pow[x] > max / scale)
		{
			cp_img[x] = psf_pow[x];
		}
	}

	row[0] = yp;
	col[0] = xp;
	p = 1;
	/* the peak coordinates have been put in the coordinate lists */
	cp_img[yp*size + xp] = 0.;

	while (num0 != num)
	{
		nump = num - num0;
		num0 = p;
		for (x = num0 - nump; x < num0; x++)
		{
			if ((row[x] - 1) > -1 && fabs(cp_img[col[x] + (row[x] - 1) * size]) > 0)
			{
				cp_img[col[x] + (row[x] - 1) * size] = 0.;
				row[p] = row[x] - 1;
				col[p] = col[x];
				p++;
			}
			if ((row[x] + 1) < size && fabs(cp_img[col[x] + (row[x] + 1) * size]) > 0)
			{
				cp_img[col[x] + (row[x] + 1) * size] = 0.;
				row[p] = row[x] + 1;
				col[p] = col[x];
				p++;
			}
			if ((col[x] - 1) > -1 && fabs(cp_img[col[x] - 1 + row[x] * size]) > 0)
			{
				cp_img[col[x] - 1 + row[x] * size] = 0.;
				row[p] = row[x];
				col[p] = col[x] - 1;
				p++;
			}
			if ((col[x] + 1) < size && fabs(cp_img[col[x] + 1 + row[x] * size]) > 0)
			{
				cp_img[col[x] + 1 + row[x] * size] = 0.;
				row[p] = row[x];
				col[p] = col[x] + 1;
				p++;
			}
		}
		num = p;
	}
	paras->psf_hlr = sqrt(p / Pi);
	delete[] cp_img;
	delete[] col;
	delete[] row;
}


int source_detector(const double *source_img, int *source_x, int*source_y, double*source_paras, para*paras, bool cross)
{	/* remember to add the peak detection! to locate the source */
	/* it will not change the inputted array */
	int i, j, k, m, c, ix, iy, tx, ty, x, y, y_size = paras->img_y, x_size = paras->img_x;
	int s = y_size*x_size;
	int peak = 0, yp, xp, area = 0, half_light_area = 0;
	double flux = 0, flux_sq = 0, half_light_flux = 0;
	int  len0 = 0, len=0, s_num = 0, num0, num, num_new;
	double detect_thres = paras->detect_thres;
	double *cp_img = new double[s] {};
	int *temp_x = new int[s] {};
	int *temp_y = new int[s] {};
	int cross_x[4]{ -1,1,0,0 };
	int cross_y[4]{ 0,0,-1,1 };
	
	/* copy and mask the candidates pixels */
	for (i = 0; i < s; i++)
	{
			if (source_img[i] >= detect_thres)
			{
				cp_img[i] = source_img[i];
			}
	}
	
	/* search the source by FoF */	
	for (i = 0; i < y_size; i++)
	{
		for (j = 0; j < x_size; j++)
		{
			if (cp_img[i*x_size + j] > 0)
			{	
				peak = cp_img[i*x_size + j];
				half_light_area = 0;
				half_light_flux = 0;
				flux = cp_img[i*x_size + j];
				flux_sq += flux*flux;

				len = 0;
				num0 = 0;
				num = 1;
				/* the temp_x, _y are the temporal array to store the coordinates of the source detected*/
				temp_x[len] = j;
				temp_y[len] = i;
				len = 1;
				cp_img[i*x_size + j] = 0;
				while (num0 != num)
				{
					num_new = num - num0;
					num0 = len;
					for (k = num0 - num_new; k < num0; k++)
					{	
						if (cross)/* search the nearest four pixels */
						{	
							for (c = 0; c < 4; c++)
							{
								tx = temp_x[k] + cross_x[c];
								ty = temp_y[k] + cross_y[c];
								if (ty<y_size&&ty >= 0 && tx<x_size&&tx >= 0 && cp_img[ty*x_size + tx]>0)
								{
									temp_x[len] = tx;
									temp_y[len] = ty;
									if (cp_img[ty*x_size + tx] > peak) /* to find the peak of the source */
									{
										peak = cp_img[ty*x_size + tx];
										yp = ty;
										xp = tx;
									}
									flux = flux + cp_img[ty*x_size + tx];
									flux_sq += cp_img[ty*x_size + tx] * cp_img[ty*x_size + tx];
									cp_img[ty*x_size + tx] = 0; /* mask each pixel after detection */
									
									len++;
								}
							}							
						}

						else /* search the nearest nine pixels */
						{
							for (iy = -1; iy < 2; iy++)
							{
								ty = temp_y[k] + iy;
								for (ix = -1; ix < 2; ix++)
								{
									tx = temp_x[k] + ix;									
									if (ty<y_size&&ty >= 0 && tx<x_size&&tx >= 0 && cp_img[ty*x_size + tx]>0)
									{
										temp_x[len] = tx;
										temp_y[len] = ty;
										if (cp_img[ty*x_size + tx] > peak) /* to find the peak of the source */
										{
											peak = cp_img[ty*x_size + tx];
											yp = ty;
											xp = tx;
										}
										flux = flux + cp_img[ty*x_size + tx];
										flux_sq += cp_img[ty*x_size + tx] * cp_img[ty*x_size + tx];
										cp_img[ty*x_size + tx] = 0;
										len++;
									}
								}
							}
						}
					}
					num = len;
				}

				if (len >= paras->area_thres)
				{	
					if (s_num >= paras->max_source)
					{
						std::cout << "Too many source!" << std::endl;
						break;
					}
					for (m = 0; m < len; m++)
					{

						if (source_img[temp_y[m] * y_size + temp_x[m]] >= peak*0.5) /* to calculate the half light radius */							
						{
							half_light_area++;
							half_light_flux = half_light_flux + source_img[temp_y[m] * y_size + temp_x[m]];
						}
						source_x[len0 + m] = temp_x[m];
						source_y[len0 + m] = temp_y[m];						
					}

					source_paras[8 * s_num] = len; /* length of source */
					source_paras[8 * s_num + 1] = yp; /* 'y' of peak of source */
					source_paras[8 * s_num + 2] = xp; /* 'x' of peak of source */
					source_paras[8 * s_num + 3] = peak ; /* peak value of source */
					source_paras[8 * s_num + 4] = half_light_area; /* length of source that the pixel value bigger than 0.5*peak */
					source_paras[8 * s_num + 5] = flux; /* total flux of source */
					source_paras[8 * s_num + 6] = half_light_flux; /* half light flux */
					source_paras[8 * s_num + 7] = flux_sq; /* sum of square of flux of each source pixel */
					len0 += len;
					s_num++;
				}
			}
		}
	}
	delete[] cp_img;
	delete[] temp_x;
	delete[] temp_y;
	return s_num;
}

int source_detector(const float *source_img, int *source_x, int*source_y, float*source_paras, para*paras, bool cross)
{	/* remember to add the peak detection! to locate the source */
	/* it will not change the inputted array */
	int i, j, k, m, c, ix, iy, tx, ty, x, y, y_size = paras->img_y, x_size = paras->img_x;
	int s = y_size * x_size;
	int peak = 0, yp, xp, area = 0, half_light_area = 0;
	float flux = 0, flux_sq = 0, half_light_flux = 0;
	int  len0 = 0, len = 0, s_num = 0, num0, num, num_new;
	float detect_thres = paras->detect_thres;
	float *cp_img = new float[s] {};
	int *temp_x = new int[s] {};
	int *temp_y = new int[s] {};
	int cross_x[4]{ -1,1,0,0 };
	int cross_y[4]{ 0,0,-1,1 };

	/* copy and mask the candidates pixels */
	for (i = 0; i < s; i++)
	{
		if (source_img[i] >= detect_thres)
		{
			cp_img[i] = source_img[i];
		}
	}

	/* search the source by FoF */
	for (i = 0; i < y_size; i++)
	{
		for (j = 0; j < x_size; j++)
		{
			if (cp_img[i*x_size + j] > 0)
			{
				peak = cp_img[i*x_size + j];
				half_light_area = 0;
				half_light_flux = 0;
				flux = cp_img[i*x_size + j];
				flux_sq += flux * flux;

				len = 0;
				num0 = 0;
				num = 1;
				/* the temp_x, _y are the temporal array to store the coordinates of the source detected*/
				temp_x[len] = j;
				temp_y[len] = i;
				len = 1;
				cp_img[i*x_size + j] = 0;
				while (num0 != num)
				{
					num_new = num - num0;
					num0 = len;
					for (k = num0 - num_new; k < num0; k++)
					{
						if (cross)/* search the nearest four pixels */
						{
							for (c = 0; c < 4; c++)
							{
								tx = temp_x[k] + cross_x[c];
								ty = temp_y[k] + cross_y[c];
								if (ty < y_size&&ty >= 0 && tx < x_size&&tx >= 0 && cp_img[ty*x_size + tx]>0)
								{
									temp_x[len] = tx;
									temp_y[len] = ty;
									if (cp_img[ty*x_size + tx] > peak) /* to find the peak of the source */
									{
										peak = cp_img[ty*x_size + tx];
										yp = ty;
										xp = tx;
									}
									flux = flux + cp_img[ty*x_size + tx];
									flux_sq += cp_img[ty*x_size + tx] * cp_img[ty*x_size + tx];
									cp_img[ty*x_size + tx] = 0; /* mask each pixel after detection */

									len++;
								}
							}
						}

						else /* search the nearest nine pixels */
						{
							for (iy = -1; iy < 2; iy++)
							{
								ty = temp_y[k] + iy;
								for (ix = -1; ix < 2; ix++)
								{
									tx = temp_x[k] + ix;
									if (ty < y_size&&ty >= 0 && tx < x_size&&tx >= 0 && cp_img[ty*x_size + tx]>0)
									{
										temp_x[len] = tx;
										temp_y[len] = ty;
										if (cp_img[ty*x_size + tx] > peak) /* to find the peak of the source */
										{
											peak = cp_img[ty*x_size + tx];
											yp = ty;
											xp = tx;
										}
										flux = flux + cp_img[ty*x_size + tx];
										flux_sq += cp_img[ty*x_size + tx] * cp_img[ty*x_size + tx];
										cp_img[ty*x_size + tx] = 0;
										len++;
									}
								}
							}
						}
					}
					num = len;
				}

				if (len >= paras->area_thres)
				{
					if (s_num >= paras->max_source)
					{
						std::cout << "Too many source!" << std::endl;
						break;
					}
					for (m = 0; m < len; m++)
					{

						if (source_img[temp_y[m] * y_size + temp_x[m]] >= peak * 0.5) /* to calculate the half light radius */
						{
							half_light_area++;
							half_light_flux = half_light_flux + source_img[temp_y[m] * y_size + temp_x[m]];
						}
						source_x[len0 + m] = temp_x[m];
						source_y[len0 + m] = temp_y[m];
					}

					source_paras[8 * s_num] = len; /* length of source */
					source_paras[8 * s_num + 1] = yp; /* 'y' of peak of source */
					source_paras[8 * s_num + 2] = xp; /* 'x' of peak of source */
					source_paras[8 * s_num + 3] = peak; /* peak value of source */
					source_paras[8 * s_num + 4] = half_light_area; /* length of source that the pixel value bigger than 0.5*peak */
					source_paras[8 * s_num + 5] = flux; /* total flux of source */
					source_paras[8 * s_num + 6] = half_light_flux; /* half light flux */
					source_paras[8 * s_num + 7] = flux_sq; /* sum of square of flux of each source pixel */
					len0 += len;
					s_num++;
				}
			}
		}
	}
	delete[] cp_img;
	delete[] temp_x;
	delete[] temp_y;
	return s_num;
}


int galaxy_finder(const double *stamp_arr, int *check_mask, para *paras, bool cross)
{	
	int elem_unit = 8; // the number of parameters for each source detected 
	int source_num, area=0, hlr_area, yp, xp, pix_num = paras->stamp_size*paras->stamp_size;
	int detect=-1, xc = paras->stamp_size / 2, yc = paras->stamp_size / 2;
	int tag_s = 0, tag_e, i, j;
	double hlr, flux, radius, max_distance=paras->max_distance*paras->max_distance;
	int *source_x = new int[pix_num] {};
	int *source_y = new int[pix_num] {};
	int *mask = new int[pix_num] {};
	double *source_para = new double[elem_unit*paras->max_source]{}; // determined by 'max_sources' in paras.
	source_num = source_detector(stamp_arr, source_x, source_y, source_para, paras, cross);
	//std::cout << source_num << std::endl;
	for ( i = 0; i < source_num; i++)
	{	
		// start point of source_y(x) of i'th source
		if (i > 0)
		{
			tag_s += source_para[(i-1)*elem_unit];
		}
		else
		{
			tag_s = 0;
		}
		for (j = tag_s; j < tag_s + source_para[i * elem_unit]; j++)
		{
			if (source_y[j] == yc &&  source_x[j] == xc)
			{
				if (source_para[i * elem_unit] > area)
				{
					area = source_para[i * elem_unit];
					detect = i;
					//std::cout << detect << std::endl;
				}
			}
		}
		radius = (source_para[i * elem_unit + 1] - xc)*(source_para[i * elem_unit + 1] - xc) + (source_para[i * elem_unit + 2] - yc)*(source_para[i * elem_unit + 2] - yc);
		//std::cout << "PARA  "<<source_para[i * elem_unit + 1]<<" "<< source_para[i * elem_unit + 2]<<" "<< source_para[i*elem_unit] <<" "<<max_distance<< std::endl;
		if (radius <= max_distance) // if it peaks within 5.5 pixels from the stamp center, it will be indentified as a galaxy
		{
			if (source_para[i * elem_unit] > area)
			{
				area = source_para[i * elem_unit];
				detect = i;
				//std::cout << i << " " << source_para[i * elem_unit] << " " << detect << " " << i << std::endl;
			}
		}
	}
	tag_s = 0;
	for (i = 0; i < detect; i++)
	{
		tag_s += source_para[i*elem_unit];
	}
	//std::cout << detect << "  "<<tag_s<<" "<<area<<std::endl;
	if (detect > -1)
	{
		double temp_flux;
		int area_ext;

		for (i = 0;  i < 5; i++)
		{
			temp_flux = 0;
			initialize_arr(mask, pix_num);
			area_ext = edge_extend(mask, source_y, source_x, tag_s, area, paras, 2*i+1);
			for (j = 0; j < pix_num; j++)
			{
				if (mask[j] > 0)
				{
					temp_flux += stamp_arr[j];
				}
				if (i == 4)
				{
					check_mask[j] = mask[j];
				}
			}
			paras->gal_flux_ext[i] = fabs(temp_flux);
			paras->gal_size_ext[i] = area_ext;
		}
		
		
		paras->gal_size = area;
		paras->gal_py = source_para[detect * elem_unit + 1];
		paras->gal_px = source_para[detect * elem_unit + 2];
		paras->gal_peak = source_para[detect * elem_unit + 3];
		paras->gal_hsize = source_para[detect * elem_unit + 4];
		paras->gal_flux = source_para[detect * elem_unit + 5];
		paras->gal_hflux = source_para[detect * elem_unit + 6];

		paras->gal_hlr = sqrt(area / Pi);
		paras->gal_snr = sqrt(source_para[detect * elem_unit + 7]) / paras->gal_noise_sig;
		paras->gal_osnr = source_para[detect * elem_unit + 5] / sqrt(area) / paras->gal_noise_sig;
	}
	else
	{
		for (i = 0; i < 5; i++)
		{
			paras->gal_flux_ext[i] = 0;
			paras->gal_size_ext[i] = 0;
			paras->gal_flux2_ext[i] = 0;
		}
		initialize_para(paras); // set the relative parameters to zero
	}
	delete[] source_x;
	delete[] source_y;
	delete[] source_para;
	delete[] mask;
	return detect;
}

int galaxy_finder(const float *stamp_arr, int *check_mask, para *paras, bool cross)
{
	int elem_unit = 8; // the number of parameters for each source detected 
	int source_num, area = 0, hlr_area, yp, xp, pix_num = paras->stamp_size*paras->stamp_size;
	int detect = -1, xc = paras->stamp_size / 2, yc = paras->stamp_size / 2;
	int tag_s = 0, tag_e, i, j;
	float hlr, flux, radius, max_distance = paras->max_distance*paras->max_distance;
	int *source_x = new int[pix_num] {};
	int *source_y = new int[pix_num] {};
	int *mask = new int[pix_num] {};
	float *source_para = new float[elem_unit*paras->max_source]{}; // determined by 'max_sources' in paras.
	source_num = source_detector(stamp_arr, source_x, source_y, source_para, paras, cross);
	//std::cout << source_num << std::endl;
	for (i = 0; i < source_num; i++)
	{
		// start point of source_y(x) of i'th source
		if (i > 0)
		{
			tag_s += source_para[(i - 1)*elem_unit];
		}
		else
		{
			tag_s = 0;
		}
		for (j = tag_s; j < tag_s + source_para[i * elem_unit]; j++)
		{
			if (source_y[j] == yc && source_x[j] == xc)
			{
				if (source_para[i * elem_unit] > area)
				{
					area = source_para[i * elem_unit];
					detect = i;
					//std::cout << detect << std::endl;
				}
			}
		}
		radius = (source_para[i * elem_unit + 1] - xc)*(source_para[i * elem_unit + 1] - xc) + (source_para[i * elem_unit + 2] - yc)*(source_para[i * elem_unit + 2] - yc);
		//std::cout << "PARA  "<<source_para[i * elem_unit + 1]<<" "<< source_para[i * elem_unit + 2]<<" "<< source_para[i*elem_unit] <<" "<<max_distance<< std::endl;
		if (radius <= max_distance) // if it peaks within 5.5 pixels from the stamp center, it will be indentified as a galaxy
		{
			if (source_para[i * elem_unit] > area)
			{
				area = source_para[i * elem_unit];
				detect = i;
				//std::cout << i << " " << source_para[i * elem_unit] << " " << detect << " " << i << std::endl;
			}
		}
	}
	tag_s = 0;
	for (i = 0; i < detect; i++)
	{
		tag_s += source_para[i*elem_unit];
	}
	//std::cout << detect << "  "<<tag_s<<" "<<area<<std::endl;
	if (detect > -1)
	{
		float temp_flux;
		int area_ext;

		for (i = 0; i < 5; i++)
		{
			temp_flux = 0;
			initialize_arr(mask, pix_num);
			area_ext = edge_extend(mask, source_y, source_x, tag_s, area, paras, 2 * i + 1);
			for (j = 0; j < pix_num; j++)
			{
				if (mask[j] > 0)
				{
					temp_flux += stamp_arr[j];
				}
				if (i == 4)
				{
					check_mask[j] = mask[j];
				}
			}
			paras->gal_flux_ext[i] = fabs(temp_flux);
			paras->gal_size_ext[i] = area_ext;
		}


		paras->gal_size = area;
		paras->gal_py = source_para[detect * elem_unit + 1];
		paras->gal_px = source_para[detect * elem_unit + 2];
		paras->gal_peak = source_para[detect * elem_unit + 3];
		paras->gal_hsize = source_para[detect * elem_unit + 4];
		paras->gal_flux = source_para[detect * elem_unit + 5];
		paras->gal_hflux = source_para[detect * elem_unit + 6];

		paras->gal_hlr = sqrt(area / Pi);
		paras->gal_snr = sqrt(source_para[detect * elem_unit + 7]) / paras->gal_noise_sig;
		paras->gal_osnr = source_para[detect * elem_unit + 5] / sqrt(area) / paras->gal_noise_sig;
	}
	else
	{
		for (i = 0; i < 5; i++)
		{
			paras->gal_flux_ext[i] = 0;
			paras->gal_size_ext[i] = 0;
			paras->gal_flux2_ext[i] = 0;
		}
		initialize_para(paras); // set the relative parameters to zero
	}
	delete[] source_x;
	delete[] source_y;
	delete[] source_para;
	delete[] mask;
	return detect;
}

int edge_extend(int *mask, const int *source_y, const int* source_x, const int source_id, const int source_len, para *paras, const int iters)
{
	int size = paras->stamp_size, pix_len=0, pix_len_0, pix_new,ix, iy, i, j, m,n,sub=2;
	int *cp_y = new int[size*size]{};
	int *cp_x = new int[size*size]{};
	char buffer[150];
	for (i = source_id; i < source_len+ source_id; i++)
	{
		// copy of source coordinates 
		cp_y[i - source_id] = source_y[i];
		cp_x[i - source_id] = source_x[i];
		// label the source pixel on the mask
		mask[source_y[i] * size + source_x[i]] = 1;
	}
	
	pix_len_0 = 0;
	pix_len = source_len;

	for (i = 0; i < iters; i++)
	{	
		// count the newly added pixels
		pix_new = pix_len - pix_len_0;
		// label the source length at the biginning
		pix_len_0 = pix_len;
		// each time search around the newly added pixels
		for (j = pix_len_0 - pix_new; j < pix_len_0; j++)
		{
			for (m = -1; m < 2; m++)
			{
				iy = cp_y[j] + m;
				if (iy >= 0 && iy < size)
				{
					for (n = -1; n < 2; n++)
					{
						ix = cp_x[j] + n;
						if (ix >= 0 && ix < size && mask[iy*size + ix] == 0)
						{	
							mask[iy*size + ix] = sub;
							// add the newly found pixels to source list
							// and increase the pixel counts
							cp_y[pix_len] = iy;
							cp_x[pix_len] = ix;
							//sprintf(buffer, "%d, (%d, %d) --> (%d, %d), %d", j, cp_y[j], cp_x[j], iy, ix, pix_len);
							//std::cout << buffer << std::endl;
							pix_len++;
						}
					}
				}
			}
		}
		sub++;
	}
	
	delete[] cp_y;
	delete[] cp_x;
	return pix_len;
}


void addnoise(double *image, const int pixel_num, const double sigma)
{
	for (int i = 0; i < pixel_num; i++)
	{
		image[i] = image[i] + gsl_ran_gaussian(rng, sigma);
	}
}

void addnoise(float *image, const int pixel_num, const float sigma)
{
	for (int i = 0; i < pixel_num; i++)
	{
		image[i] = image[i] + gsl_ran_gaussian(rng, sigma);
	}
}


void initialize_arr(double *in_img, const int length)
{/* will set all the elements to zero */
	for (int i = 0; i < length; i++)
	{
		in_img[i] = 0.;
	}
}

void initialize_arr(float *in_img, const int length)
{/* will set all the elements to zero */
	for (int i = 0; i < length; i++)
	{
		in_img[i] = 0.;
	}
}

void initialize_arr(int *in_img, const int length)
{/* will set all the elements to zero */
	for (int i = 0; i < length; i++)
	{
		in_img[i] = 0;
	}
}


void normalize_arr(double * arr,const int size)
{
	int i;
	double peak=0;
	for (i = 0; i < size*size; i++)
	{
		if (arr[i] > peak)
		{
			peak = arr[i];
		}
	}
	if (peak != 0)
	{
		for (i = 0; i < size*size; i++)
		{
			arr[i] = arr[i] / peak;
		}
	}
	else
	{
		std::cout << "Divided by ZERO!!!" << std::endl;
		exit(0);
	}
}

void normalize_arr(float * arr, const int size)
{
	int i;
	float peak = 0;
	for (i = 0; i < size*size; i++)
	{
		if (arr[i] > peak)
		{
			peak = arr[i];
		}
	}
	if (peak != 0)
	{
		for (i = 0; i < size*size; i++)
		{
			arr[i] = arr[i] / peak;
		}
	}
	else
	{
		std::cout << "Divided by ZERO!!!" << std::endl;
		exit(0);
	}
}


/********************************************************************************************************************************************/
/* Fourier Quad */
/********************************************************************************************************************************************/

void snr_est(const double *image, para *paras, int fit)
{	/* will not change the inputted array */
	/* estimate the snr in Fourier space */
	double n = 0, noise;
	int size = paras->stamp_size;
	int i, k, edge = 1, xc = size / 2;
	int x[20]{ -1,  0,  1, -2, -1,  0,  1,  2, -2, -1,  1,  2, -2, -1,  0,  1,  2, -1,  0,  1 };
	int y[20]{ -2, -2, -2, -1, -1, -1, -1, -1,  0,  0,  0,  0,  1,  1,  1,  1,  1,  2, 2, 2 };
	double fz[20], fit_paras[6];

	for (i = 0; i < size; i++) //y coordinates
	{
		for (k = 0; k < size; k++) // x coordinates
		{
			if (i< edge || i> size - edge - 1)
			{
				n += image[i*size + k];
			}
			else if (k < edge || k>size - edge - 1)
			{
				n += image[i*size + k];
			}
		}
	}
	noise = sqrt(n*0.25 / ((size - edge)*edge));
	paras->gal_flux2 = sqrt(image[xc*size + xc]) / noise;
	for (i = 0; i < 5; i++)
	{
		paras->gal_flux2_ext[i] = paras->gal_flux_ext[i] / noise;
	}

	if (fit == 2)
	{
		for (i = 0; i < 20; i++)
		{
			fz[i] = image[(xc + y[i])*size + xc + x[i]];
		}
		hyperfit_5(fz, fit_paras, paras);
		paras->gal_flux_alt = sqrt(pow(10, fit_paras[0]))/ noise;
	}
}

void possion_subtraction(double *arr, para *paras, int edge)
{
	int i, j,size = paras->stamp_size;
	double noise = 0;
	for (i = 0; i < size; i++) //y
	{
		for (j = 0; j < size; j++) // x
		{
			if (i< edge || i> size - edge - 1)
			{
				noise += arr[i*size + j];
			}
			else if (j < edge || j>size - edge - 1)
			{
				noise += arr[i*size + j];
			}
		}
	}
	noise = noise * 0.25 / ((size - edge)*edge);
	for (i = 0; i < size*size; i++) 
	{
		arr[i] = arr[i] - noise;
	}
}

void noise_subtraction(double *image_pow, double *noise_pow, para *paras, const int edge, const int possion)
{
	int i, j, size = paras->stamp_size;
	double inoise = 0, pnoise=0;

	if (0 == possion)
	{
		for (i = 0; i < size*size; i++)
		{
			image_pow[i] = image_pow[i] - noise_pow[i];
		}
	}
	else
	{
		for (i = 0; i < size; i++) //y
		{
			for (j = 0; j < size; j++) // x
			{
				if (i< edge || i> size - edge - 1)
				{
					inoise += image_pow[i*size + j];
					pnoise += noise_pow[i*size + j];
				}
				else if (j < edge || j>size - edge - 1)
				{
					inoise += image_pow[i*size + j];
					pnoise += noise_pow[i*size + j];
				}
			}
		}
		inoise = inoise * 0.25 / ((size - edge)*edge);
		pnoise = pnoise * 0.25 / ((size - edge)*edge);

		for (i = 0; i < size*size; i++)
		{
			image_pow[i] = image_pow[i] - inoise - noise_pow[i] + pnoise;
		}
	}
}

void shear_est(double *gal_img, double *psf_img, para *paras)
{	 /* will not change the inputted array */
	 /* all the inputted images are the powerspectrums */
	/* if there's no backgroud noise, a array of '0' should be inputted */
	double mg1 = 0., mg2 = 0., mn = 0., mu = 0., mv = 0., beta, thres, alpha, kx, kx2, ky2, ky, tk, k2, k4;
	double mp1=0., mp2=0.;
	int i, j, k, size;
	size = paras->stamp_size;

	alpha = 16*Pi*Pi*Pi*Pi/ (size*size*size*size);
	/* beta is the beta_square in the estimators */
	beta = 1./ paras->psf_hlr / paras->psf_hlr;
	
	//find the maximum of psf power spectrum and set the threshold of max/10000 above which the pixel value will be taken into account
	thres = paras->psf_pow_thres;

	for (i = 0; i < size; i++)//y coordinates
	{
		ky = i - size*0.5;
		for (j = 0; j < size; j++) // x coordinates
		{
			kx = j - size*0.5;
			if (psf_img[i*size + j] >= thres)
			{	
				k2 = kx*kx + ky*ky;
				k4 = k2*k2;
				kx2 = kx*kx;
				ky2 = ky*ky;

				tk = exp( - k2 * beta ) / psf_img[i*size + j] * gal_img[i*size + j] * alpha;
				//tk = exp(-k2 * beta) / psf_img[i*size + j] * (gal_img[i*size + j] - noise_img[i*size + j]) * alpha;
				mg1 += -0.5 * ( kx2 - ky2 ) * tk;
				mg2 += -kx*ky*tk;
				mn += ( k2 - 0.5*beta*k4 ) * tk;
				mu += (k4  - 8 *kx2*ky2)*tk * (-0.5*beta);
				mv += kx*ky*(kx2-ky2)* tk * (-2.* beta);
				//mp1 += (-4.*(kx*kx - ky*ky) + 8.*beta*( pow(kx, 4) - pow(ky, 4) ) - 2.*beta*beta*( pow(kx, 6) + pow(kx, 4)*ky*ky - kx*kx*pow(ky, 4) - pow(ky, 6) ) )*tk;
				//mp2 += ( -8.*kx*ky + 16.*beta*( kx*kx*kx*ky + kx*ky*ky*ky ) - 4*beta*beta*( pow(kx, 5)*ky + 2*kx*kx*kx*ky*ky*ky + kx*pow(ky, 5) ) )*tk;				
			}
		}
	}

	paras->n1 = mg1;
	paras->n2 = mg2;
	paras->dn = mn;
	paras->du = mu;
	paras->dv = mv;
	//paras->dp1 = mp1;
	//paras->dp2 = mp2;

}


/********************************************************************************************************************************************/
/* fitting */
/********************************************************************************************************************************************/
void smooth(double *image,  const double *coeffs, para*paras)//be careful of the memset()
{
	/*  to fit the curve: a1 + a2*x +a3*y + a4*x^2 +a5*x*y + a6*y^2  */
	int i, j, m, n, q, p, pk = 0, tag, cen, coe, jx, iy, size = paras->stamp_size;
	double fz[6]{}, z[25]{}, fit_para_6, max = 0., thres;
	double*temp = new double[size*size]{};
	double fit_temp[6 * 25]{};
	double ones[6]{ 1.,1.,1.,1.,1.,1. };

	cen = (size*size + size)*0.5; // the position of the k0 of the power spectrum
	for (i = 0; i < size*size; i++)
	{
		temp[i] = log10(image[i]);
	}
	for (i = 0; i < size; i++) //y
	{
		for (j = 0; j < size; j++)//x
		{
			tag = 0;
			pk = 0;
			for (m = -2; m < 3; m++)
			{
				if (2 < i && i < size - 2)
				{
					p = (i + m)*size;
				}
				else
				{
					p = (i + m + size) % size*size;// the periodic boundry condition
				}
				for (n = -2; n < 3; n++)
				{
					if (2 < j && j < size - 2)
					{
						q = j + n;
					}
					else
					{
						q = (j + n + size) % size; // the periodic boundry condition
					}

					if (tag != 0 && tag != 4 && tag != 20 && tag != 24)
						// exclude the corner and center of each 5*5 block and the k0 of the power spectrum
					{
						if (p + q != cen)
						{
							z[tag] = temp[p + q];
						}
						else
						{
							pk = tag; // tag dicides the row of invers coefficients matrix to be used
						}
					}
					tag++;
				}
			}
			coe = pk * 150;
			for (q = 0; q < 150; q++)
			{
				fit_temp[q] = coeffs[coe + q];
			}
			cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 25, 1, fit_temp, 25, z, 1, 0, fz, 1);
			fit_para_6 = cblas_ddot(6, fz, 1, ones, 1);
			image[i*size + j] = pow(10., fit_para_6);
			memset(fz, 0, sizeof(fz));
			memset(z, 0, sizeof(z));			
		}
	}
	delete[] temp;
}

void smooth(double *image, const double* psf_pow, const double *coeffs, para*paras)//be careful of the memset()
{
	/*  to fit the curve: a1 + a2*x +a3*y + a4*x^2 +a5*x*y + a6*y^2  */
	int i, j, m, n, q, p, pk = 0, tag, cen, coe, jx, iy, size = paras->stamp_size;
	double fz[6]{}, z[25]{}, fit_para_6, max = 0., thres = paras->psf_pow_thres/1.2; // the lower thres (divided by 1.2) for safety!!
	double*temp = new double[size*size]{};
	double fit_temp[6 * 25]{};
	double ones[6]{ 1.,1.,1.,1.,1.,1. };

	cen = (size*size + size)*0.5; // the position of the k0 of the power spectrum
	for (i = 0; i < size*size; i++)
	{
		temp[i] = log10(image[i]);
	}
	
	for (i = 0; i < size; i++) //y
	{
		for (j = 0; j < size; j++)//x
		{
			if (psf_pow[i*size + j] >= thres)
			{				
				tag = 0;
				pk = 0;
				for (m = -2; m < 3; m++)
				{
					if (2 < i && i < size - 2)
					{
						p = (i + m)*size;
					}
					else
					{
						p = (i + m + size) % size*size;// the periodic boundry condition
					}
					for (n = -2; n < 3; n++)
					{
						if (2 < j && j < size - 2)
						{
							q = j + n;
						}
						else
						{
							q = (j + n + size) % size; // the periodic boundry condition
						}

						if (tag != 0 && tag != 4 && tag != 20 && tag != 24)
							// exclude the corner and center of each 5*5 block and the k0 of the power spectrum
						{
							if (p + q != cen)
							{
								z[tag] = temp[p + q];
							}
							else
							{
								pk = tag; // tag dicides the row of invers coefficients matrix to be used
							}
						}
						tag++;
					}
				}
				coe = pk * 150;
				for (q = 0; q < 150; q++)
				{
					fit_temp[q] = coeffs[coe + q];
				}

				cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 25, 1, fit_temp, 25, z, 1, 0, fz, 1);
				fit_para_6 = cblas_ddot(6, fz, 1, ones, 1);
				image[i*size + j] = pow(10., fit_para_6);

				memset(fz, 0, sizeof(fz));
				memset(z, 0, sizeof(z));
			}
		}
	}
	delete[] temp;
}

void smooth_real(double *image, const double *coeffs, para*paras)//be careful of the memset()
{
	/*  to fit the curve: a1 + a2*x +a3*y + a4*x^2 +a5*x*y + a6*y^2  */
	int i, j, m, n, q, p, pk = 0, tag, cen, coe, jx, iy, size = paras->stamp_size;
	double fz[6]{}, z[25]{}, fit_para_6, flux=0;
	double*temp = new double[size*size]{};
	double fit_temp[6 * 25]{};
	double ones[6]{ 1.,1.,1.,1.,1.,1. };

	for (q = 0; q < 150; q++)
	{
		fit_temp[q] = coeffs[q];
	}

	cen = (size*size + size)*0.5; // the position of the k0 of the power spectrum
	for (i = 0; i < size*size; i++)
	{
		temp[i] = image[i];
	}
	for (i = 0; i < size; i++) //y
	{
		for (j = 0; j < size; j++)//x
		{
			tag = 0;

			for (m = -2; m < 3; m++)
			{
				if (2 < i && i < size - 2)
				{
					p = (i + m)*size;
				}
				else
				{
					p = (i + m + size) % size*size;// the periodic boundry condition
				}
				for (n = -2; n < 3; n++)
				{
					if (2 < j && j < size - 2)
					{
						q = j + n;
					}
					else
					{
						q = (j + n + size) % size; // the periodic boundry condition
					}

					if (tag != 0 && tag != 4 && tag != 20 && tag != 24)
						// exclude the corner and center of each 5*5 block 
					{
						z[tag] = temp[p + q];
					}
					tag++;
				}
			}

			cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 25, 1, fit_temp, 25, z, 1, 0, fz, 1);
			fit_para_6 = cblas_ddot(6, fz, 1, ones, 1);
			//image[i*size + j] = fit_para_6;
			flux += fit_para_6;
			memset(fz, 0, sizeof(fz));
			memset(z, 0, sizeof(z));
		}
	}
	paras->gal_total_flux = fabs(flux);
	delete[] temp;
}

void hyperfit_5(double *data, double *fit_paras, para *paras)
{
	double temp = 0;

	for (int i = 0; i < 6; i++)
	{
		temp = 0;
		for (int j = 0; j < 20; j++)
		{
			temp += paras->fit_matrix[i][j] * log10(data[j]);
		}
		fit_paras[i] = temp;
	}
}


/********************************************************************************************************************************************/
/* general methods */
/********************************************************************************************************************************************/
void initialize_para(para *paras)
{	
	paras->gal_size = 0;
	paras->gal_hsize = 0;
	paras->gal_px = 0;
	paras->gal_py = 0;
	paras->gal_peak = 0;
	paras->gal_hlr = 0;
	paras->gal_flux = 0;
	paras->gal_hflux = 0;
	paras->gal_fluxsq = 0;
	paras->gal_flux2 = 0;
	paras->gal_flux_alt = 0;
	paras->gal_snr = 0;
	paras->gal_osnr = 0;
}

void set_bin(const double *data, const int data_num, double * bins, const int bin_num)
{
	double *data_cp = new double[data_num];
	int i, mid = bin_num / 2, step = data_num / bin_num * 2;
	for (i = 0; i < data_num; i++)
	{
		data_cp[i] = fabs(data[i]);
	}
	sort_double(data_cp, data_num, 1);
	bins[0] = -data_cp[data_num - 1] * 100;
	bins[bin_num] = data_cp[data_num - 1] * 100;
	bins[mid] = 0;
	for (i = 1; i < bin_num / 2; i++)
	{
		bins[mid + i] = data_cp[step*i];
		bins[mid - i] = -data_cp[step*i];
	}
	delete[] data_cp;
}

void set_bin(const float *data, const int data_num, float * bins, const int bin_num)
{
	float *data_cp = new float[data_num];
	int i, mid = bin_num / 2, step = data_num / bin_num * 2;
	for (i = 0; i < data_num; i++)
	{
		data_cp[i] = fabs(data[i]);
	}
	sort_float(data_cp, data_num, 1);
	bins[0] = -data_cp[data_num - 1] * 100.;
	bins[bin_num] = data_cp[data_num - 1] * 100.;
	bins[mid] = 0;
	for (i = 1; i < bin_num / 2; i++)
	{
		bins[mid + i] = data_cp[step*i];
		bins[mid - i] = -data_cp[step*i];
	}
	delete[] data_cp;
}

void set_bin(const int *data, const int data_num, int * bins, const int bin_num)
{
	int *data_cp = new int[data_num];
	int i, mid = bin_num / 2, step = data_num / bin_num * 2;
	for (i = 0; i < data_num; i++)
	{
		data_cp[i] = fabs(data[i]);
	}
	sort_int(data_cp, data_num, 1);
	bins[0] = -data_cp[data_num - 1] * 100;
	bins[bin_num] = data_cp[data_num - 1] * 100;
	bins[mid] = 0;
	for (i = 1; i < bin_num / 2; i++)
	{
		bins[mid + i] = data_cp[step*i];
		bins[mid - i] = -data_cp[step*i];
	}
	delete[] data_cp;
}

void histogram(const double *data, const double *bins, int *num, const int data_num, const int bin_num)
{
	for (int i = 0; i < data_num; i++)
	{
		for ( int j = 0; j < bin_num; j++)
		{
			if (data[i] < bins[j + 1] && data[i] >= bins[j])
			{
				num[j] += 1;
				break;
			}
		}
	}
}

void histogram(const float *data, const float *bins, int *num, const int data_num, const int bin_num)
{	
	for (int i = 0; i < data_num; i++)
	{
		for (int j = 0; j < bin_num; j++)
		{
			if (data[i] < bins[j + 1] && data[i] >= bins[j])
			{
				num[j] += 1;
				break;
			}
		}
	}
}

void histogram(const int *data, const  int *bins, int *num, const  int data_num, const  int bin_num)
{
	for (int i = 0; i < data_num; i++)
	{
		for (int j = 0; j < bin_num; j++)
		{
			if (data[i] < bins[j + 1] && data[i] >= bins[j])
			{
				num[j] += 1;
				break;
			}
		}
	}
}

void histogram2d(const double *data_y, const double*data_x, const double *bin_y, const double *bin_x, int *num, const int data_num, const int ybin_num, const int xbin_num)
{
	// bin_num = len(bins) - 1
	for (int k = 0; k < data_num; k++)
	{
		// loop y-bins
		for (int i = 0; i < ybin_num; i++)
		{
			if (data_y[k] < bin_y[i + 1] && data_y[k] >= bin_y[i])
			{
				// loop x-bins
				for (int j = 0; j < xbin_num; j++)
				{
					if (data_x[k] < bin_x[j + 1] && data_x[k] >= bin_x[j])
					{
						num[i * xbin_num + j] += 1;
						break;
					}
				}
				break;
			}
		}
	}
}

void histogram2d(const float *data_y, const float*data_x, const float *bin_y, const float *bin_x, int *num, const int data_num, const  int ybin_num, const int xbin_num)
{
	// bin_num = len(bins) - 1
	for (int k = 0; k < data_num; k++)
	{
		// loop y-bins
		for (int i = 0; i < ybin_num; i++)
		{
			if (data_y[k] < bin_y[i + 1] && data_y[k] >= bin_y[i])
			{
				// loop x-bins
				for (int j = 0; j < xbin_num; j++)
				{
					if (data_x[k] < bin_x[j + 1] && data_x[k] >= bin_x[j])
					{
						num[i * xbin_num + j] += 1;
						break;
					}
				}
				break;
			}
		}
	}
}

void histogram2d(const int *data_y, const int*data_x, const int *bin_y, const int *bin_x, int *num, const int data_num, const int ybin_num, const int xbin_num)
{
	// bin_num = len(bins) - 1
	for (int k = 0; k < data_num; k++)
	{
		// loop y-bins
		for (int i = 0; i < ybin_num; i++)
		{
			if (data_y[k] < bin_y[i + 1] && data_y[k] >= bin_y[i])
			{
				// loop x-bins
				for (int j = 0; j < xbin_num; j++)
				{
					if (data_x[k] < bin_x[j + 1] && data_x[k] >= bin_x[j])
					{
						num[i * xbin_num + j] += 1;
						break;
					}
				}
				break;
			}
		}
	}
}

void sort_double(double* arr, int size, int order=1)
{
	if (order == 1)
	{
		std::sort(arr, arr+size, std::less<double>());
	}
	else
	{
		std::sort(arr, arr+size,std::greater<double>());
	}
}

void sort_float(float *arr, int size, int order=1)
{
	if (order == 1)
	{
		std::sort(arr, arr+size,  std::less<float>());
	}
	else
	{
		std::sort(arr, arr+size, std::greater<float>());
	}
}

void sort_int(int *arr, int size, int order=1)
{
	if (order == 1)
	{
		std::sort(arr, arr+size, std::less<int>());
	}
	else
	{
		std::sort(arr, arr+size, std::greater<int>());
	}
}

int com_double_ascend(const double a, const double b)
{
	return a < b;
}

int com_float_ascend(const float a, const float b)
{
	return a<b;
}

int com_int_ascend(const int a, const int b)
{
	return a<b;
}

int com_double_descend(const double a, const double b)
{
	return a>b;
}

int com_float_descend(const float a, const float b)
{
	return a>b;
}

int com_int_descend(const int a, const int b)
{
	return a>b;
}

void get_time(char *str_time, int length)
{	
	/* the length of str_time should be larger than 40 */
	time_t timer;
	time(&timer);
	strftime(str_time, length, "%Y-%m-%d %H:%M:%S", localtime(&timer));
}



/*double quick_exp(double x, double precision, para* paras)
{
	int i, j, m, mm;
	double temp=1, expo = 0, n;
	if (x >1 or x<-1)
	{
		m = (int)x;
		n = x - m;
	}
	else
	{
		m = 0;
		n = x;
	}
	for (i = 0; temp >= precision; i++)
	{
		temp *= n / i;
		expo += temp;
	}

	if(m>0)
	{
		expo *= paras->exp_val[m-1];
	}
	else if (m < 0)
	{	
		m = abs(m);

		if (m >= -15)
		{
			expo /= paras->exp_val[abs(m + 1)];
		}
		else
		{

		}
	}

	return expo;
}*/

/********************************************************************************************************************************************/
/* GSL library */
/********************************************************************************************************************************************/
void gsl_rng_initialize(int seed)
{
	gsl_rng_env_setup();
	//gsl_rng_default_seed = seed;
	//T = gsl_rng_mt19937;
	T = gsl_rng_ranlxs0;
	rng = gsl_rng_alloc(T);
	gsl_rng_set(rng, seed);
}

void gsl_rng_free()
{
	gsl_rng_free(rng);
}

