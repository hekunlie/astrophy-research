
//#include "stdafx.h"
#include "FQlib.h"

using namespace std;
const gsl_rng_type *T;
gsl_rng *rng;

void write_h5(char *filename, char *set_name, int row, int column, double *matrix)
{
	hid_t file_id;
	herr_t status;
	remove(filename);
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	// �������ݼ���metadata�е�dataspace��Ϣ��Ŀ
	unsigned rank = 2;
	hsize_t dims[2];
	dims[0] = row;
	dims[1] = column;
	hid_t dataspace_id;  // ���ݼ�metadata��dataspace��id
						 // dataspace_id = H5Screate_simple(int rank, �ռ�ά��
						 //              const hsize_t* current_dims, ÿ��ά��Ԫ�ظ���
						 //                    - ����Ϊ0����ʱ�޷�д������
						 //                  const hsize_t* max_dims, ÿ��ά��Ԫ�ظ�������
						 //                    - ��ΪNULLָ�룬���current_dim��ͬ��
						 //                    - ��ΪH5S_UNLIMITED���������ޣ���datasetһ���Ƿֿ��(chunked).
	dataspace_id = H5Screate_simple(rank, dims, NULL);

	// �������ݼ��е����ݱ���
	hid_t dataset_id;    // ���ݼ������id
						 // dataset_id = H5Dcreate(loc_id, λ��id
						 //              const char *name, ���ݼ���
						 //                hid_t dtype_id, ��������
						 //                hid_t space_id, dataspace��id
						 //             ����(link)��������,
						 //                 ���ݼ���������,
						 //                 ���ݼ���������)
	dataset_id = H5Dcreate(file_id, set_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// ������д�����ݼ�
	// herr_t д��״̬ = H5Dwrite(д��Ŀ�����ݼ�id,
	//                               �ڴ����ݸ�ʽ,
	//                       memory_dataspace_id, �����ڴ�dataspace�����е�ѡ��
	//                          - H5S_ALL: �ļ���dataspace�����ڴ�dataspace��file_dataspace_id�е�ѡ����Ϊ�ڴ�dataspace��ѡ��
	//                         file_dataspace_id, �����ļ���dataspace��ѡ��
	//                          - H5S_ALL: �ļ���datasapce��ȫ��������Ϊ���ݼ���dataspace�����ȫ��ά������
	//                        ����IO������ת������,
	//                          const void * buf, �ڴ������ݵ�λ��
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix);

	// �ر�dataset��ض���
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);

	// �ر��ļ�����
	status = H5Fclose(file_id);
}

void read_img(double *arr, char *path)
{
	fitsfile *fptr;															/* pointer to the FITS file, defined in fitsio.h */
	int status = 0, nfound, anynull;
	long naxes[2], fpixel = 1, nbuffer, npixels, ii;
	double datamin, datamax, nullval = 0;

	fits_open_file(&fptr, path, READONLY, &status);
	fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);		/* read the NAXIS1 and NAXIS2 keyword to get image size */
	npixels = naxes[0] * naxes[1];											/* number of pixels in the image, python_arr[naxes[1], naxes[0]] */
	double *buffer = new double[npixels];										/* create a new array */

	fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval, buffer, &anynull, &status);

	for (ii = 0; ii < npixels; ii++) arr[ii] = buffer[ii];
	delete[] buffer;														/* (have to) delete the array */

	fits_close_file(fptr, &status);
}

void write_img(double *img, int ysize, int xsize, char *filename)
{
	fitsfile *fptr;															 /* pointer to the FITS file; defined in fitsio.h */
	int status, ii, jj;
	long fpixel = 1, naxis = 2, nelements, exposure;
	long naxes[2] ;								/* x, y */
	naxes[0] = xsize;
	naxes[1] = ysize;
																//float *array = new float[xsize*ysize];					/* the same as numpy array */
																	//char path[50] = "E:\\testfile.fits";
	remove(filename);                                                /* overwrite the old file */
	status = 0;																						/* initialize status before calling fitsio routines */
	fits_create_file(&fptr, filename, &status);												/* create new file */
	fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);				 /* Create the primary array image (16-bit short integer pixels */
	exposure = 1500.;
	fits_update_key(fptr, TDOUBLE, "EXPOSURE", &exposure, "Total Exposure Time", &status);  /* Write a keyword; must pass the ADDRESS of the value */
	nelements = xsize*ysize;         /* number of pixels to write */
	fits_write_img(fptr, TDOUBLE, fpixel, nelements, img, &status);     /* Write the array of integers to the image */
	fits_close_file(fptr, &status);              /* close the file */
	fits_report_error(stderr, status);      /* print out any error messages */
}

void pow_spec(double *in_img, double *out_img, int column, int row)
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

void get_radius(double *in_img, para *paras, double scale, int size, int type, double sig_level)
{	 /* will not change the inputted array */
	/* the image should be larger than 12*12 */
	/* setting scale = infinity ,  one can obtain the area of the signal */

	int x, y, xp = 0, yp = 0, num0 = 0, num = 1, nump, p = 1;
	double max = 0, flux = 0, flux_sq=0. ;	

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
	double *cp_img = new double[size*size]();
	for (x = 0; x < size*size; x++)
	{	
		if (in_img[x] > max / scale && in_img[x] > 1.5 * sig_level)
		{
			cp_img[x] = in_img[x];
		}
	}
	
	int *col = new int[size*size];
	int *row = new int[size*size];
	col[0] = xp;
	row[0] = yp;
	flux = max;
	flux_sq = max*max;
	/* the peak coordinates have been put in the coordinate lists */
	cp_img[xp + yp*size] = 0.;

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
		if (max > 1.5 * sig_level) 
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
			points_r[i] = (1. + g1)*(rot1 * points[i] - rot2*points[i + num_p] ) + g2*(rot2 * points[i] + rot1*points[i + num_p] ) + size/2.;
			points_r[i + num_p] = g2*(rot1 * points[i] - rot2*points[i + num_p]) + (1. - g1)*(rot2 * points[i] + rot1*points[i + num_p]) +size / 2.;
		}
	}
	else
	{
		for (i = 0; i < num_p; i++)
		{
			points_r[i] = (1. + g1)* points[i] + g2*points[i + num_p] + size/2.;
			points_r[i + num_p] = g2*points[i]  + (1. - g1)*points[i + num_p] + size / 2.;
		}
	}

	/*  convolve PSF and draw the image */
	flux_g = flux / (2* Pi *scale*scale);     /* 1 / sqrt(2*Pi*sig_x^2)/sqrt(2*Pi*sig_x^2) */
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

void create_points(double *point, int num_p, double radius)
{	/* point is the container of the coordinates of the points created which should has 2*num_p elements ,[x..., y....]*/
	int i;
	double x = 0., y = 0., xm=0., ym=0., theta, step;

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
		point[i] = point[i] - xm/num_p;
		point[i + num_p] = point[i + num_p] - ym / num_p;
	}
}

void create_epoints(double *point, int num_p, double radius, double beta, double ellip)
{
	int i=0, j;
	double theta, r, b;
	b = sqrt((1 - ellip) / (1 + ellip))*radius;

	while (i < num_p)
	{
		r =radius* gsl_rng_uniform(rng);
		theta = 2 * Pi*gsl_rng_uniform(rng);
		if (1. / sqrt(cos(theta + beta)*cos(theta + beta) / radius / radius + sin(theta + beta)*sin(theta + beta) / b / b) > r)
		{
			point[i] = r*cos(theta);
			point[i + num_p] = r*sin(theta);
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

void shear_est(double *gal_img, double *psf_img, double *noise_img, para *paras, int size)
{	 /* will not change the inputted array */
	 /* all the inputted images are the powerspectrums */
	/* if there's no backgroud noise, a list of '0' should be putted in  */
	double mg1 = 0., mg2 = 0., mn = 0., mu = 0., mv = 0., beta, max = 0, thres, alpha, kx, ky, tk;
	double mp1=0., mp2=0.;
	int i, j, k;

	alpha = 16*Pi *Pi*Pi*Pi/ (size*size*size*size);
	/* beta is the beta_square in the estimators */
	beta = 1./ paras->psf_hlr / paras->psf_hlr;
	
	//find the maximum of psf power spectrum and set the threshold of max/10000 above which the pixel value will be taken account
	for (i = 0; i < size*size; i++)
	{
		if (psf_img[i] > max)
			max = psf_img[i];
	}
	thres = max / 10000.;

	for (i = 0; i < size; i++)//y coordinates
	{
		ky = i - size*0.5;
		for (j = 0; j < size; j++) // x coordinates
		{
			kx = j - size*0.5;
			if (psf_img[i*size + j] > thres)
			{
				tk = exp( - ( kx*kx + ky*ky ) * beta ) / psf_img[i*size + j] * (gal_img[i*size + j] - noise_img[i*size+j]) * alpha;
				mg1 += -0.5 * ( kx*kx - ky*ky ) * tk;
				mg2 += -kx*ky*tk;
				mn += ( kx*kx + ky*ky - 0.5*beta*(kx*kx + ky*ky)*(kx*kx + ky*ky) ) * tk;
				mu += (kx*kx*kx*kx - 6 * kx*kx*ky*ky + ky*ky*ky*ky)*tk * (-0.5*beta);
				mv += (kx*kx*kx*ky - ky*ky*ky*kx)* tk * (-2.* beta);
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
	//paras->dp1 = 0;//mp1;
	//paras->dp2 = 0;// mp2;

}

void initialize(double *in_img, int length)
{/* will set all the elements to zero */
	int i;
	for (i = 0; i < length; i++)
		in_img[i] = 0.;
}

void stack(double *container, double *stamp, int tag, int size, int row, int col)
{
	int i,j,m;
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

void segment(double *chip, double *stamp, int tag, int size, int row, int col)
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

void addnoise(double *image, int pixel_num,   double sigma)
{
	int i;

	for (i = 0; i < pixel_num; i++)
	{
		image[i] = image[i] + gsl_ran_gaussian(rng, sigma); 
	}
}

void f_snr(double *image, para *paras, int size, int edge)
{	/* will not change the inputted array */
	/* estimate the snr in Fourier space */
	double noise=0, signal;
	int i,k=0;
	for (i = 0; i < size; i++) //y coordinates
	{
		for (k = 0; k < size; k++) // x coordinates
		{
			if(i< edge||i> size-edge-1 )
			{
				noise += image[i*size + k];
			}
			else if (k < edge || k>size - edge-1)
			{
				noise += image[i*size + k];
			}
		}
	}
	if (image[size / 2 * size + size / 2] > image[(size / 2 - 1)* size + size / 2] && image[(size / 2 * size) + size / 2] > image[(size / 2 + 1) * size + size / 2]
		&& image[size / 2 * size + size / 2] > image[size / 2 * size + size / 2 - 1] && image[size / 2 * size + size / 2] > image[size / 2 * size + size / 2 + 1])
	{
		signal = image[size / 2 * size + size / 2];
	}
	else
	{
		signal = (image[(size / 2 - 1)* size + size / 2] + image[(size / 2 + 1) * size + size / 2] + image[size / 2 * size + size / 2 + 1] + image[size / 2 * size + size / 2 - 1]) *0.25;
	}
	noise = noise*0.25 / ((size - edge)*edge);
	paras->gal_fsnr = sqrt(signal/noise);
}
