#include "FQlib.h"

const gsl_rng_type *T0, *T1, *T2, *T3;
gsl_rng *rng0, *rng1, *rng2, *rng3;
int GSL_SETUP_LABEL=0;


/********************************************************************************************************************************************/
/* operations on the image */
/********************************************************************************************************************************************/

void stack(double *container, const double *stamp, const int tag, const int size, const int row, const int col)
{
	int i, j, m;
	/* 'container' is the container of 'row'*'col' stamps which is array with 'size'*'size' pixels.
	'tag' is the tag of the stamp in all stamps.
	the (i,j)'th pixel of stamp will locate at '(tag - tag%col)/col*size*col*size + i*n*size + tag%col8size + j' */
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



void source_detector(const double *source_img, int *source_x, int*source_y, double*source_paras, fq_paras*paras, bool cross, int &detection, std::string &info)
{	/* remember to add the peak detection! to locate the source */
	/* it will not change the inputted array */
	int i, j, k, m, c, ix, iy, tx, ty, x, y, y_size = paras->img_y, x_size = paras->img_x;
	int s = y_size*x_size;
	int peak = 0, yp, xp, area = 0, half_light_area = 0;
	double flux = 0, flux_sq = 0, half_light_flux = 0;
	int  len0 = 0, len=0, s_num = 0, num0, num, num_new;
	double detect_thresh = paras->detect_thresh;
	double *cp_img = new double[s] {};
	int *temp_x = new int[s] {};
	int *temp_y = new int[s] {};
	int cross_x[4]{ -1,1,0,0 };
	int cross_y[4]{ 0,0,-1,1 };
	
	/* copy and mask the candidates pixels */
	for (i = 0; i < s; i++)
	{
			if (source_img[i] >= detect_thresh)
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
				peak = 0;
				half_light_area = 0;
				half_light_flux = 0;
				flux = cp_img[i*x_size + j];
				flux_sq = flux*flux;

				// if(i==29 and j ==30)
				// {
				// 	std::cout<<"Begin peak: "<<peak<<std::endl;
				// }

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
						else /* search the nearest eight pixels */
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

											// if(i==29 and j ==30)
											// {
											// 	std::cout<<"Finding: "<<peak<<" "<<yp<<" "<<xp<<" "<<ty<<" "<<tx<<std::endl;
											// }
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

				if (len >= paras->area_thresh)
				{	
					//std::cout<<len<<" "<<yp<<" "<<xp<<std::endl;
					if (s_num >= paras->max_source)
					{
						std::cout << "Too many source!" << std::endl;
						info = "Too many source!";
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
	if (info != "Too many source!")
	{
		info = "Normal.";
	}
	delete[] cp_img;
	delete[] temp_x;
	delete[] temp_y;
	detection = s_num;
}

void source_detector(const float *source_img, int *source_x, int*source_y, float*source_paras, fq_paras_float*paras, bool cross, int &detection, std::string &info)
{	/* remember to add the peak detection! to locate the source */
	/* it will not change the inputted array */
	int i, j, k, m, c, ix, iy, tx, ty, x, y, y_size = paras->img_y, x_size = paras->img_x;
	int s = y_size * x_size;
	int peak = 0, yp, xp, area = 0, half_light_area = 0;
	float flux = 0, flux_sq = 0, half_light_flux = 0;
	int  len0 = 0, len = 0, s_num = 0, num0, num, num_new;
	float detect_thresh = paras->detect_thresh;
	float *cp_img = new float[s] {};
	int *temp_x = new int[s] {};
	int *temp_y = new int[s] {};
	int cross_x[4]{ -1,1,0,0 };
	int cross_y[4]{ 0,0,-1,1 };

	/* copy and mask the candidates pixels */
	for (i = 0; i < s; i++)
	{
		if (source_img[i] >= detect_thresh)
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
				peak = 0;
				half_light_area = 0;
				half_light_flux = 0;
				flux = cp_img[i*x_size + j];
				flux_sq = flux * flux;

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

				if (len >= paras->area_thresh)
				{
					if (s_num >= paras->max_source)
					{
						std::cout << "Too many source!" << std::endl;
						info = "Too many source!";
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
	if (info != "Too many source!")
	{
		info = "Normal.";
	}
	delete[] cp_img;
	delete[] temp_x;
	delete[] temp_y;
	detection = s_num;
}



void galaxy_finder(const double *stamp_arr, int *check_mask, fq_paras *paras, bool cross, int &detect_label, std::string &info)
{	
	int elem_unit = 8; // the number of parameters for each source detected 
	int source_num, area = 0, hlr_area, yp, xp;
	int size = paras->stamp_size, pix_num = size * size;
	int xc = size / 2, yc = size / 2;
	int detect = -1;

	int tag_s = 0, tag_e, i, j;

	double hlr, flux, radius, max_distance = paras->max_distance*paras->max_distance;
	double dy, dx;
	int *source_x = new int[pix_num] {};
	int *source_y = new int[pix_num] {};

	double *source_para = new double[elem_unit*paras->max_source]{}; // determined by 'max_sources' in paras.
	initialize_arr(check_mask, pix_num, 0);
	source_detector(stamp_arr, source_x, source_y, source_para, paras, cross, source_num, info);
	//std::cout << source_num << std::endl;

	int max_tag=-1, max_area=0;
	for ( i = 0; i < source_num; i++)
	{		
		dy = source_para[i * elem_unit + 1] - yc;
		dx = source_para[i * elem_unit + 2] - xc;
		radius = dy*dy+ dx*dx;
		area = source_para[i * elem_unit];
		if (radius <= max_distance and area >= max_area)
		{
			max_tag = i;
			max_area = area;
		}
		// std::cout<<dy<<" "<<dx<<" "<<area<<" "<<max_tag<<" "<<max_area<<std::endl;
	}
	detect = max_tag;
	// std::cout<<"Max tag "<<detect<<std::endl;
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
			// detection mask
			check_mask[source_y[j] * size + source_x[j]] = 1;
		}
	
	}

	if (detect > -1)
	{
		// mask the target galaxy
		tag_s = 0;
		int tag_e;
		for (i = 0; i < detect; i++)
		{
			tag_s += source_para[i*elem_unit];
		}
		tag_e = tag_s + source_para[detect*elem_unit];
		for (i = tag_s; i < tag_e; i++)
		{
			check_mask[source_y[i] * size + source_x[i]] = 2;
		}

		paras->gal_size = max_area;
		paras->gal_py = source_para[detect * elem_unit + 1];
		paras->gal_px = source_para[detect * elem_unit + 2];
		paras->gal_peak = source_para[detect * elem_unit + 3];
		paras->gal_hsize = source_para[detect * elem_unit + 4];
		paras->gal_flux = source_para[detect * elem_unit + 5];
		paras->gal_hflux = source_para[detect * elem_unit + 6];

		paras->gal_effective_radius = sqrt(max_area / Pi);
		paras->gal_hlr = sqrt(paras->gal_hsize / Pi);
		paras->gal_snr = sqrt(source_para[detect * elem_unit + 7]) / paras->gal_noise_sig;
		paras->gal_osnr = source_para[detect * elem_unit + 5] / sqrt(max_area) / paras->gal_noise_sig;
	}
	else
	{
		initialize_para(paras); // set the relative parameters to zero
	}
	delete[] source_x;
	delete[] source_y;
	delete[] source_para;
	detect_label = detect;
}

void galaxy_finder(const float *stamp_arr, int *check_mask, fq_paras_float *paras, bool cross, int &detect_label, std::string &info)
{
	int elem_unit = 8; // the number of parameters for each source detected 
	int source_num, area = 0, hlr_area, yp, xp;
	int size = paras->stamp_size, pix_num = size * size;
	int xc = size / 2, yc = size / 2;
	int detect = -1;

	int tag_s = 0, tag_e, i, j;

	double hlr, flux, radius, max_distance = paras->max_distance*paras->max_distance;
	double dy, dx;
	int *source_x = new int[pix_num] {};
	int *source_y = new int[pix_num] {};

	float *source_para = new float[elem_unit*paras->max_source]{}; // determined by 'max_sources' in paras.
	initialize_arr(check_mask, pix_num, 0);
	source_detector(stamp_arr, source_x, source_y, source_para, paras, cross, source_num, info);
	//std::cout << source_num << std::endl;

	int max_tag=-1, max_area=0;
	for ( i = 0; i < source_num; i++)
	{		
		dy = source_para[i * elem_unit + 1] - yc;
		dx = source_para[i * elem_unit + 2] - xc;
		radius = dy*dy+ dx*dx;
		area = source_para[i * elem_unit];
		if (radius <= max_distance and area >= max_area)
		{
			max_tag = i;
			max_area = area;
		}
		// std::cout<<dy<<" "<<dx<<" "<<area<<" "<<max_tag<<" "<<max_area<<std::endl;
	}
	detect = max_tag;
	// std::cout<<"Max tag "<<detect<<std::endl;
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
			// detection mask
			check_mask[source_y[j] * size + source_x[j]] = 1;
		}
	
	}

	if (detect > -1)
	{
		// mask the target galaxy
		tag_s = 0;
		int tag_e;
		for (i = 0; i < detect; i++)
		{
			tag_s += source_para[i*elem_unit];
		}
		tag_e = tag_s + source_para[detect*elem_unit];
		for (i = tag_s; i < tag_e; i++)
		{
			check_mask[source_y[i] * size + source_x[i]] = 2;
		}

		paras->gal_size = max_area;
		paras->gal_py = source_para[detect * elem_unit + 1];
		paras->gal_px = source_para[detect * elem_unit + 2];
		paras->gal_peak = source_para[detect * elem_unit + 3];
		paras->gal_hsize = source_para[detect * elem_unit + 4];
		paras->gal_flux = source_para[detect * elem_unit + 5];
		paras->gal_hflux = source_para[detect * elem_unit + 6];

		paras->gal_effective_radius = sqrt(max_area / Pi);
		paras->gal_hlr = sqrt(paras->gal_hsize / Pi);
		paras->gal_snr = sqrt(source_para[detect * elem_unit + 7]) / paras->gal_noise_sig;
		paras->gal_osnr = source_para[detect * elem_unit + 5] / sqrt(max_area) / paras->gal_noise_sig;
	}
	else
	{
		initialize_para(paras); // set the relative parameters to zero
	}
	delete[] source_x;
	delete[] source_y;
	delete[] source_para;
	detect_label = detect;
}

int edge_extend(int *mask, const int *source_y, const int* source_x, const int source_id, const int source_len, fq_paras *paras, const int iters)
{
	int size = paras->stamp_size, pix_len=0, pix_len_0, pix_new,ix, iy, i, j, m,n,sub=2;
	int *cp_y = new int[size*size]{};
	int *cp_x = new int[size*size]{};
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



void create_points(double *point, const int num_p, const double radius, const double step, const gsl_rng *gsl_rand_rng)
{	/* point is the container of the coordinates of the points created which should has 2*num_p elements ,[x..., y....]*/
	int i;
	double x = 0., y = 0., xm = 0., ym = 0., theta;

	for (i = 0; i < num_p; i++)
	{
		theta = 2.*Pi*gsl_rng_uniform(gsl_rand_rng);
		//step = 1*gsl_rng_uniform(rng);
		x += cos(theta)*step;
		y += sin(theta)*step;
		if (x*x + y*y > radius*radius)
		{
			x = cos(theta)*step;
			y = sin(theta)*step;
		}
		point[i] = x;
		point[i + num_p] = y;
		xm += x;
		ym += y;
	}
	
	//theta = 2.*Pi*gsl_rng_uniform(gsl_rand_rng);
	//x = cos(theta);
	//y = sin(theta);

	/* move the mass center of the points cluster to (0,0) */
	xm /= num_p;
	ym /= num_p;
	for (i = 0; i < num_p; i++)
	{
		point[i] = point[i] - xm;
		point[i + num_p] = point[i + num_p] - ym;
	}
}


void create_points(float *point, const int num_p, const float radius, const float step, const gsl_rng *gsl_rand_rng)
{	/* point is the container of the coordinates of the points created which should has 2*num_p elements ,[x..., y....]*/
	int i;
	float x = 0., y = 0., xm = 0., ym = 0., theta;

	for (i = 0; i < num_p; i++)
	{
		theta = 2.*Pi*gsl_rng_uniform(gsl_rand_rng);
		//step = 1*gsl_rng_uniform(rng);
		x += cos(theta)*step;
		y += sin(theta)*step;
		if (x*x + y*y > radius*radius)
		{
			x = cos(theta)*step;
			y = sin(theta)*step;
		}
		point[i] = x;
		point[i + num_p] = y;
		xm += x;
		ym += y;
	}
	
	//theta = 2.*Pi*gsl_rng_uniform(gsl_rand_rng);
	//x = cos(theta);
	//y = sin(theta);

	/* move the mass center of the points cluster to (0,0) */
	xm /= num_p;
	ym /= num_p;
	for (i = 0; i < num_p; i++)
	{
		point[i] = point[i] - xm;
		point[i + num_p] = point[i + num_p] - ym;
	}
}

void create_epoints(double *point, const int num_p, const double ellip, const gsl_rng *gsl_rand_rng)
{
	int i = 0;
	double theta, beta, radius, r1, r2, b, y, x;
	radius = 5 + 4 * gsl_rng_uniform(gsl_rand_rng);
	b = sqrt((1 - ellip) / (1 + ellip))*radius;
	beta = 2 * Pi*gsl_rng_uniform(gsl_rand_rng);
	r1 = cos(beta);
	r2 = sin(beta);
	while (i < num_p)
	{
		x = gsl_ran_gaussian(gsl_rand_rng, radius);
		y = gsl_ran_gaussian(gsl_rand_rng, b);
		if (x*x / (radius*radius) + y*y / (b*b) <= 1)
		{
			point[i] = x*r1 - y*r2;
			point[i + num_p] = x*r2 + r1*y;
			i++;
		}
	}
}

void create_epoints(float *point, const int num_p, const float ellip, const gsl_rng *gsl_rand_rng)
{
	int i = 0;
	float theta, beta, radius, r1, r2, b, y, x;
	radius = 5 + 4 * gsl_rng_uniform(gsl_rand_rng);
	b = sqrt((1 - ellip) / (1 + ellip))*radius;
	beta = 2 * Pi*gsl_rng_uniform(gsl_rand_rng);
	r1 = cos(beta);
	r2 = sin(beta);
	while (i < num_p)
	{
		x = gsl_ran_gaussian(gsl_rand_rng, radius);
		y = gsl_ran_gaussian(gsl_rand_rng, b);
		if (x*x / (radius*radius) + y*y / (b*b) <= 1)
		{
			point[i] = x*r1 - y*r2;
			point[i + num_p] = x*r2 + r1*y;
			i++;
		}
	}
}


void coord_rotation(const double*xy, const int pts_num, const double theta, double *xy_r)
{
    double sin_theta, cos_theta;
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    for (int i = 0; i < pts_num; i++)
    {
        xy_r[i] = cos_theta * xy[i] - sin_theta*xy[i + pts_num];
        xy_r[i + pts_num] = sin_theta * xy[i] + cos_theta*xy[i + pts_num];
        // std::cout<<xy[i]<<" "<<xy[i+pts_num]<<" "<<xy_r[i]<<" "<<xy_r[i+pts_num]<<std::endl;
    }
}

void coord_rotation(const float*xy, const int pts_num, const float theta, float *xy_r)
{
    float sin_theta, cos_theta;
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    for (int i = 0; i < pts_num; i++)
    {
        xy_r[i] = cos_theta * xy[i] - sin_theta*xy[i + pts_num];
        xy_r[i + pts_num] = sin_theta * xy[i] + cos_theta*xy[i + pts_num];
        // std::cout<<xy[i]<<" "<<xy[i+pts_num]<<" "<<xy_r[i]<<" "<<xy_r[i+pts_num]<<std::endl;
    }
}


void create_psf(double*in_img, const double scale, const int size, const double img_cent, const int psf)
{
	int i, j;
	double rs, r1, val, flux_g, flux_m, rd;
	
	flux_g = 1. / (2 * Pi *scale*scale);     /* 1 / sqrt(2*Pi*sig_x^2)/sqrt(2*Pi*sig_x^2) */
	flux_m = 1. / (Pi*scale*scale*(1. - pow(10, -2.5))*0.4); /* 1 / ( Pi*scale^2*( (1 + alpha^2)^(1-beta) - 1) /(1-beta)), where alpha = 3, beta = 3.5 */

	initialize_arr(in_img, size*size, 0);
	rd = 1. / scale / scale;
	for (i = 0; i < size; i++)
	{
		r1 = (i - img_cent)*(i - img_cent)*rd;
		for (j = 0; j < size; j++)
		{
			rs = r1 + (j - img_cent)*(j - img_cent)*rd;
			if (psf == 1)  // Gaussian PSF
			{
				if (rs <= 9) in_img[i*size + j] = flux_g*exp(-rs*0.5);
			}
			else              // Moffat PSF
			{
				if (rs <= 9.) in_img[i*size + j] = flux_m*pow(1. + rs, -3.5);
			}
		}
	}
}

void create_psf(float*in_img, const float scale, const int size, const float img_cent, const int psf)
{
	int i, j;
	float rs, r1, val, flux_g, flux_m, rd;
	
	flux_g = 1. / (2 * Pi *scale*scale);     /* 1 / sqrt(2*Pi*sig_x^2)/sqrt(2*Pi*sig_x^2) */
	flux_m = 1. / (Pi*scale*scale*(1. - pow(10, -2.5))*0.4); /* 1 / ( Pi*scale^2*( (1 + alpha^2)^(1-beta) - 1) /(1-beta)), where alpha = 3, beta = 3.5 */

	initialize_arr(in_img, size*size, 0);
	rd = 1. / scale / scale;
	for (i = 0; i < size; i++)
	{
		r1 = (i - img_cent)*(i - img_cent)*rd;
		for (j = 0; j < size; j++)
		{
			rs = r1 + (j - img_cent)*(j - img_cent)*rd;
			if (psf == 1)  // Gaussian PSF
			{
				if (rs <= 9) in_img[i*size + j] = flux_g*exp(-rs*0.5);
			}
			else              // Moffat PSF
			{
				if (rs <= 9.) in_img[i*size + j] = flux_m*pow(1. + rs, -3.5);
			}
		}
	}
}



void create_psf_e(double*in_img, const double scale, const int size, const double img_cent, const double ellip, const double theta, const int psf)
{
	int i, j;
	double rs, val, flux_norm, rd;
	double a2_inv, b2_inv;
	double q2, rot_1, rot_2;
	double r1, r2, ry1, ry2;
	double nscale;

	initialize_arr(in_img, size*size, 0);

	rot_1 = cos(theta);
	rot_2 = -sin(theta);
	// [ cos \theta, -sin \theta]
	// [sin \theta, cos \theta]
	q2 =  (1 - ellip)/(1 + ellip);

	// x^2/a^2 + y^2/b^2 = 1, (a > b)
	// q = b/a, q^2 = (1-e)/(1+e) 
	// => q^2 x^2/b^2 + y^2/b^2 = 1

	b2_inv = 1. / scale / scale;
	a2_inv = b2_inv*q2;

	// just for Moffat PSF now!!!
	// plot the PSF to 3*rd, default
	nscale = 9;
	// integrate to sqrt(nscale)*rd
	flux_norm = sqrt(q2)/ (Pi*scale*scale*(1. - pow(1+nscale, -2.5))*0.4);

	// xr = cos\theta * x - sin\theta *y
	// yr = sin\theta * x + cos\theta *y
	for (i = 0; i < size; i++) // y
	{
		// -sin\theta *y, the sign is correct!
		ry1 =  -rot_2 * (i - img_cent);
		// cos\theta *y
		ry2 =  rot_1 * (i - img_cent);

		for (j = 0; j < size; j++) // x
		{
			// xr = cos\theta * x - sin\theta *y, the sign is correct!
			r1 = rot_1 * (j - img_cent) + ry1;
			// yr = sin\theta * x + cos\theta *y, the sign is correct!
			r2 = rot_2 * (j - img_cent) + ry2;
			// xr^2/a^2 + yr^2/b^2
			rs = r1 * r1*a2_inv + r2 * r2*b2_inv;

			if (psf == 1)  // Gaussian PSF
			{
				// the flux normalization may be wrong for Gaussian PSF!!!
				if (rs <= nscale) in_img[i*size + j] = flux_norm * exp(-rs * 0.5);
			}
			else // Moffat PSF
			{
				if (rs <= nscale) in_img[i*size + j] = flux_norm * pow(1. + rs, -3.5);
			}
		}
	}
}

void create_psf_e(float*in_img, const float scale, const int size, const float img_cent, const float ellip, const float theta, const int psf)
{
	int i, j;
	float rs, val, flux_norm, rd;
	float a2_inv, b2_inv;
	float q2, rot_1, rot_2;
	float r1, r2, ry1, ry2;
	float nscale;

	initialize_arr(in_img, size*size, 0);

	rot_1 = cos(theta);
	rot_2 = -sin(theta);
	// [ cos \theta, -sin \theta]
	// [sin \theta, cos \theta]
	q2 =  (1 - ellip)/(1 + ellip);

	// x^2/a^2 + y^2/b^2 = 1, (a > b)
	// q = b/a, q^2 = (1-e)/(1+e) 
	// => q^2 x^2/b^2 + y^2/b^2 = 1

	b2_inv = 1. / scale / scale;
	a2_inv = b2_inv*q2;

	// just for Moffat PSF now!!!
	// plot the PSF to 3*rd, default
	nscale = 9;
	// integrate to sqrt(nscale)*rd
	flux_norm = sqrt(q2)/ (Pi*scale*scale*(1. - pow(1+nscale, -2.5))*0.4);

	// xr = cos\theta * x - sin\theta *y
	// yr = sin\theta * x + cos\theta *y
	for (i = 0; i < size; i++) // y
	{
		// -sin\theta *y
		ry1 =  -rot_2 * (i - img_cent);
		// cos\theta *y
		ry2 =  rot_1 * (i - img_cent);

		for (j = 0; j < size; j++) // x
		{
			// xr = cos\theta * x - sin\theta *y
			r1 = rot_1 * (j - img_cent) + ry1;
			// yr = sin\theta * x + cos\theta *y
			r2 = rot_2 * (j - img_cent) + ry2;
			// xr^2/a^2 + yr^2/b^2
			rs = r1 * r1*a2_inv + r2 * r2*b2_inv;

			if (psf == 1)  // Gaussian PSF
			{
				if (rs <= nscale) in_img[i*size + j] = flux_norm * exp(-rs * 0.5);
			}
			else              // Moffat PSF
			{
				if (rs <= nscale) in_img[i*size + j] = flux_norm * pow(1. + rs, -3.5);
			}
		}
	}
}



void convolve(const double * points, const int num_p, const double flux_per_pts, const double g1, const double g2,
 				double *in_img, const int size, const double img_cent, const double psf_scale, const int psf_type)
{	 /* will not change the inputted array */
	 /* in_img is the container of the final image,
	 points is the array of points' coordinates,
	 rotate is the radian in units of pi/4,
	 psf_scale is the scale length of PSF,
	 psf_type=1 means the Gaussian PSF, psf_type=2 means the Moffat PSF	*/
	int i, j, k, m;
	double r1, r2, n, flux_g, flux_m;
	// |rot1, - rot2 |
	// |rot2,  rot1  |
	double rot1, rot2, val, rs, rd;
	rd = 1. / psf_scale / psf_scale;  // scale of PSF	

	double *points_r = new double[num_p * 2];

	initialize_arr(in_img, size*size, 0);

	/* shear the profile and move the center to image center */
	for (i = 0; i < num_p; i++)
	{
		points_r[i] = (1. + g1)* points[i] + g2*points[i + num_p] + img_cent;
		points_r[i + num_p] = g2*points[i] + (1. - g1)*points[i + num_p] + img_cent;
	}
	

	/*  convolve PSF and draw the image */
	flux_g = flux_per_pts / (2 * Pi *psf_scale*psf_scale);     /* 1 / sqrt(2*Pi*sig_x^2)/sqrt(2*Pi*sig_x^2) */
	/* \int_{0}^{alpha*rd}  PSF(r)_{moffat} = Pi*scale^2*( (1 + alpha^2)^(1-beta) - 1) /(1-beta)), beta = 3.5 for moffat*/
	flux_m = flux_per_pts / (Pi*psf_scale*psf_scale*(1. - pow(10, -2.5))*0.4);   /* alpha = 3, beta = 3.5 */
	//flux_m = flux / (Pi*psf_scale*psf_scale*(1. - pow(17, -2.5))*0.4);   /*  alpha = 4, beta = 3.5 */
	
	for (k = 0; k < num_p; k++)
	{
		for (i = 0; i < size; i++)  /* y coordinate */
		{
			r1 = (i - points_r[k + num_p])*(i - points_r[k + num_p])*rd;
			for (j = 0; j < size; j++) /* x coordinate */
			{
				rs = r1 + (j - points_r[k])*(j - points_r[k])*rd;
				if (psf_type == 1)  // Gaussian PSF
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

void convolve(const float * points, const int num_p, const float flux_per_pts, const float g1, const float g2, 
				float *in_img, const int size, const float img_cent, const float psf_scale, const int psf_type)
{	 /* will not change the inputted array */
	 /* in_img is the container of the final image,
	 points is the array of points' coordinates,
	 rotate is the radian in units of pi/4,
	 psf_scale is the scale length of PSF,
	 psf_type=1 means the Gaussian PSF, psf_type=2 means the Moffat PSF	*/
	int i, j, k, m;
	float r1, r2, n, flux_g, flux_m;
	// |rot1, - rot2 |
	// |rot2,  rot1  |
	float rot1, rot2, val, rs, rd;
	rd = 1. / psf_scale / psf_scale;  // scale of PSF	

	float *points_r = new float[num_p * 2];

	initialize_arr(in_img, size*size, 0);


	/* shear the profile and move the center to image center */
	for (i = 0; i < num_p; i++)
	{
		points_r[i] = (1. + g1)* points[i] + g2*points[i + num_p] + img_cent;
		points_r[i + num_p] = g2*points[i] + (1. - g1)*points[i + num_p] + img_cent;
	}


	/*  convolve PSF and draw the image */
	flux_g = flux_per_pts / (2 * Pi *psf_scale*psf_scale);     /* 1 / sqrt(2*Pi*sig_x^2)/sqrt(2*Pi*sig_x^2) */
	/* \int_{0}^{alpha*rd}  PSF(r)_{moffat} = Pi*scale^2*( (1 + alpha^2)^(1-beta) - 1) /(1-beta)), beta = 3.5 for moffat*/
	flux_m = flux_per_pts / (Pi*psf_scale*psf_scale*(1. - pow(10, -2.5))*0.4);   /* alpha = 3, beta = 3.5 */
	//flux_m = flux / (Pi*psf_scale*psf_scale*(1. - pow(17, -2.5))*0.4);   /*  alpha = 4, beta = 3.5 */
	
	for (k = 0; k < num_p; k++)
	{
		for (i = 0; i < size; i++)  /* y coordinate */
		{
			r1 = (i - points_r[k + num_p])*(i - points_r[k + num_p])*rd;
			for (j = 0; j < size; j++) /* x coordinate */
			{
				rs = r1 + (j - points_r[k])*(j - points_r[k])*rd;
				if (psf_type == 1)  // Gaussian PSF
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



void convolve_e(const double * points, const int num_p, const double flux_per_pts, const double g1, const double g2,
				double *in_img, const int size, const double img_cent, const double psf_scale, const int psf_type, const double ellip, const double theta)
{	 /* will not change the inputted array */
	 /* in_img is the container of the final image,
	 points is the array of points' coordinates,
	 rotate is the radian in units of pi/4,
	 psf_scale is the scale length of PSF,
	 psf_type=1 means the Gaussian PSF, psf_type=2 means the Moffat PSF	*/
	int i, j, k, m;
	double r1, r2, n, flux_norm;
	double ry1, ry2, rx1, rx2;
	double a2_inv, b2_inv;
	double psf_rot_1, psf_rot_2, q2;
	double nscale;

	double val, rs;

	double *points_r = new double[num_p * 2];

	initialize_arr(in_img, size*size, 0);

	// rotation of psf
	psf_rot_1 = cos(theta);
	psf_rot_2 = -sin(theta);

	q2 =  (1 - ellip)/(1 + ellip);
	b2_inv = 1. / psf_scale / psf_scale;
	a2_inv = b2_inv*q2;
	
	/* shear the profile and move the center to image center */
	for (i = 0; i < num_p; i++)
	{
		points_r[i] = (1. + g1)* points[i] + g2 * points[i + num_p] + img_cent;
		points_r[i + num_p] = g2 * points[i] + (1. - g1)*points[i + num_p] + img_cent;
	}


	/*  convolve PSF and draw the image */

	// now, just for Moffat PSF
	// plot the PSF to 3*rd, default
	nscale = 9;
	// integrate to sqrt(nscale)*rd
	flux_norm = flux_per_pts*sqrt(q2)/ (Pi*psf_scale*psf_scale*(1. - pow(1+nscale, -2.5))*0.4);

	for (k = 0; k < num_p; k++)
	{
		for (i = 0; i < size; i++)  /* y coordinate */
		{
			ry1 = -psf_rot_2 * (i - points_r[k + num_p]);
			ry2 = psf_rot_1 * (i - points_r[k + num_p]);
			for (j = 0; j < size; j++) /* x coordinate */
			{
				r1 = psf_rot_1 * (j - points_r[k]) + ry1;
				r2 = psf_rot_2 * (j - points_r[k]) + ry2;
				rs = r1*r1*a2_inv + r2*r2*b2_inv;
				if (psf_type == 1)  // Gaussian PSF
				{
					if (rs <= nscale) in_img[i*size + j] += flux_norm * exp(-rs * 0.5);
				}
				else				// Moffat PSF
				{
					if (rs <= nscale) in_img[i*size + j] += flux_norm * pow(1. + rs, -3.5);
				}
			}
		}
	}
	delete[] points_r;
}

void convolve_e(const float * points, const int num_p, const float flux_per_pts, const float g1, const float g2, 
				float *in_img, const int size, const float img_cent, const float psf_scale, const int psf_type, const float ellip, const float theta)
{	 /* will not change the inputted array */
	 /* in_img is the container of the final image,
	 points is the array of points' coordinates,
	 rotate is the radian in units of pi/4,
	 psf_scale is the scale length of PSF,
	 psf_type=1 means the Gaussian PSF, psf_type=2 means the Moffat PSF	*/
	int i, j, k, m;
	float r1, r2, n, flux_norm;
	float ry1, ry2, rx1, rx2;
	float a2_inv, b2_inv;
	float psf_rot_1, psf_rot_2, q2;
	float nscale;

	float val, rs;

	float *points_r = new float[num_p * 2];

	initialize_arr(in_img, size*size, 0);

	// rotation of psf
	psf_rot_1 = cos(theta);
	psf_rot_2 = -sin(theta);

	q2 =  (1 - ellip)/(1 + ellip);
	b2_inv = 1. / psf_scale / psf_scale;
	a2_inv = b2_inv*q2;
	

	/* shear the profile and move the center to image center */
	for (i = 0; i < num_p; i++)
	{
		points_r[i] = (1. + g1)* points[i] + g2 * points[i + num_p] + img_cent;
		points_r[i + num_p] = g2 * points[i] + (1. - g1)*points[i + num_p] + img_cent;
	}


	/*  convolve PSF and draw the image */

	// now, just for Moffat PSF
	// plot the PSF to 3*rd, default
	nscale = 9;
	// integrate to sqrt(nscale)*rd
	flux_norm = flux_per_pts*sqrt(q2)/ (Pi*psf_scale*psf_scale*(1. - pow(1+nscale, -2.5))*0.4);

	for (k = 0; k < num_p; k++)
	{
		for (i = 0; i < size; i++)  /* y coordinate */
		{
			ry1 = -psf_rot_2 * (i - points_r[k + num_p]);
			ry2 = psf_rot_1 * (i - points_r[k + num_p]);
			for (j = 0; j < size; j++) /* x coordinate */
			{
				r1 = psf_rot_1 * (j - points_r[k]) + ry1;
				r2 = psf_rot_2 * (j - points_r[k]) + ry2;
				rs = r1*r1*a2_inv + r2*r2*b2_inv;
				if (psf_type == 1)  // Gaussian PSF
				{
					if (rs <= nscale) in_img[i*size + j] += flux_norm * exp(-rs * 0.5);
				}
				else				// Moffat PSF
				{
					if (rs <= nscale) in_img[i*size + j] += flux_norm * pow(1. + rs, -3.5);
				}
			}
		}
	}
	delete[] points_r;
}



void deconvolution(const double *gal_pow, const double *psf_pow, double *out_img, const int column, const int row)
{
	int i, j;
	j=column*row;
	for(i=0;i<j;i++)
	{
		out_img[i]=gal_pow[i]/psf_pow[i];
	}
}

void deconvolution(const double *gal_pow, const double*gal_pow_real, const double*gal_pow_imag, const double *psf_pow, const double*psf_pow_real, const double*psf_pow_imag, 
double *out_img, const int column, const int row)
{
	int i,j;
	j = column*row;
	for(i=0;i<j;i++)
	{
		out_img[i] = gal_pow[i]/psf_pow[i] + 4*gal_pow_real[i]*gal_pow_imag[i]*psf_pow_real[i]*psf_pow_imag[i]/psf_pow[i]/psf_pow[i];
	}
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
		in[i][0] = in_img[i];
		in[i][1] = 0;
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
		in[i][0] = in_img[i];
		in[i][1] = 0;
	}
	// out[i][0] real, out[i][1] imagine
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

void pow_spec(const double *in_img, double *out_img, double *out_img_real, double *out_img_imag, double *phase_arg, const int column, const int row)
{   /* will not change the inputted array */
	/* in_img is the inputted array and the out_img is the container of the outputted image */
	fftwl_complex *in, *out;
	fftwl_plan p;
	long int i, j, m, n;
	// out[i][0] real, out[i][1] imagine

	in = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex) *(row*column));
	for (i = 0; i < (row*column); i++)
	{
		in[i][0] = in_img[i];
		in[i][1] = 0;
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
			phase_arg[n] = atan2(out[m][1], out[m][0]);
			out_img_real[n] = out[m][0];
			out_img_imag[n] = out[m][1];
		}
	}
	/* shift the signal to the center */

	fftwl_destroy_plan(p);
	fftwl_free(in);
	fftwl_free(out);
}

void pow_spec(const float *in_img, float *out_img, float *out_img_real, float *out_img_imag,float *phase_arg, const int column, const int row)
{   /* will not change the inputted array */
	/* in_img is the inputted array and the out_img is the container of the outputted image */
	fftwf_complex *in, *out;
	fftwf_plan p;
	long int i, j, m, n;

	in = (fftwf_complex*)fftwl_malloc(sizeof(fftwf_complex) *(row*column));
	for (i = 0; i < (row*column); i++)
	{
		in[i][0] = in_img[i];
		in[i][1] = 0;
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
			phase_arg[n] = atan2(out[m][1], out[m][0]);
			out_img_real[n] = out[m][0];
			out_img_imag[n] = out[m][1];
		}
	}
	/* shift the signal to the center */

	fftwf_destroy_plan(p);
	fftwf_free(in);
	fftwf_free(out);
}


void pow_spec(const double *in_img, double *out_img_pow, double *out_img_real, double *out_img_imag, const int column, const int row)
{
    /* will not change the inputted array */
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

			out_img_pow[n] = out[m][0] * out[m][0] + out[m][1] * out[m][1];
			out_img_real[n] = out[m][0];
			out_img_imag[n] = out[m][1];
		}
	}
	/* shift the signal to the center */

	fftwl_destroy_plan(p);
	fftwl_free(in);
	fftwl_free(out);
}

void pow_spec(const float *in_img, float *out_img, float *out_img_real, float *out_img_imag,const int column, const int row)
{   /* will not change the inputted array */
	/* in_img is the inputted array and the out_img is the container of the outputted image */
	fftwf_complex *in, *out;
	fftwf_plan p;
	long int i, j, m, n;

	in = (fftwf_complex*)fftwl_malloc(sizeof(fftwf_complex) *(row*column));
	for (i = 0; i < (row*column); i++)
	{
		in[i][0] = in_img[i];
		in[i][1] = 0;
	}
	// out[i][0] real, out[i][1] imaginary
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
			out_img_real[n] = out[m][0];
			out_img_imag[n] = out[m][1];
		}
	}
	/* shift the signal to the center */

	fftwf_destroy_plan(p);
	fftwf_free(in);
	fftwf_free(out);
}


void get_radius(double *in_img, fq_paras *paras, double scale, int type, double sig_level)
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


void get_psf_radius(const double *psf_pow, fq_paras*paras, const double scale)
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
	paras->psf_pow_thresh = max / 100000.;
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

void get_psf_radius(const float *psf_pow, fq_paras_float *paras, const float scale)
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
	paras->psf_pow_thresh = max / 100000.;
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



void get_quad(const double *img, const int img_size, const double img_cent, const double weight_sigma_sq, double &quad_size)
{
    // calculate the gaussian-weighted quadrupole of a galaxy or PSF image 
    double temp_quad, temp_norm;
    int i,j,m;
    double ry_sq, r_sq;
    double wei_coeff,wei_img;
    

    wei_coeff = 0.5/weight_sigma_sq;

    temp_quad = 0;
    temp_norm = 0;
    for(i=0; i<img_size; i++)
    {   
        ry_sq = (i - img_cent)*(i-img_cent);
        m = i*img_size;
        for(j=0; j<img_size; j++)
        {
            r_sq = (j - img_cent)*(j-img_cent) + ry_sq;

            wei_img = exp(-r_sq*wei_coeff)*img[m+j];
   
            temp_quad += wei_img*r_sq;
            temp_norm += wei_img;
        }
    }
    if(temp_norm == 0){quad_size = 0;}
	else{quad_size = temp_quad/temp_norm;}    
}

void get_quad(const float *img, const int img_size, const float img_cent, const float weight_sigma_sq, float &quad_size)
{
    // calculate the gaussian-weighted quadrupole of a galaxy or PSF image 
    float temp_quad, temp_norm;
    int i,j,m;
    float ry_sq, r_sq;
    float wei_coeff,wei_img;
    

    wei_coeff = 0.5/weight_sigma_sq;

    temp_quad = 0;
    temp_norm = 0;
    for(i=0; i<img_size; i++)
    {   
        ry_sq = (i - img_cent)*(i-img_cent);
        m = i*img_size;
        for(j=0; j<img_size; j++)
        {
            r_sq = (j - img_cent)*(j-img_cent) + ry_sq;

            wei_img = exp(-r_sq*wei_coeff)*img[m+j];
   
            temp_quad += wei_img*r_sq;
            temp_norm += wei_img;
        }
    }
    if(temp_norm == 0){quad_size = 0;}
	else{quad_size = temp_quad/temp_norm;}    
}





void initialize_arr(long *arr, const int length, const long x)
{/* will set all the elements to x */
	for (int i = 0; i < length; i++)
	{
		arr[i] = x;
	}
}

void initialize_arr(double *arr, const int length, const double x)
{/* will set all the elements to x */
	for (int i = 0; i < length; i++)
	{
		arr[i] = x;
	}
}

void initialize_arr(float *arr, const int length, const float x)
{/* will set all the elements to x */
	for (int i = 0; i < length; i++)
	{
		arr[i] = x;
	}
}

void initialize_arr(int *arr, const int length, const int x)
{/* will set all the elements to x */
	for (int i = 0; i < length; i++)
	{
		arr[i] = x;
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

void snr_est(const double *image, fq_paras *paras, int fit)
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
	paras->gal_flux2_ext[0] = sqrt(image[xc*size + xc]);

	if (fit == 2)
	{
		for (i = 0; i < 20; i++)
		{
			fz[i] = image[(xc + y[i])*size + xc + x[i]];
		}
		hyperfit_5(fz, fit_paras, paras);
		paras->gal_flux_alt = sqrt(pow(10, fit_paras[0]))/ noise;
		paras->gal_flux2_ext[1] = sqrt(pow(10, fit_paras[0]));
		paras->gal_flux2_ext[2] = std::max(paras->gal_flux2_ext[1], paras->gal_flux2_ext[0]);
		paras->gal_flux2_ext[3] = std::max(paras->gal_flux2, paras->gal_flux_alt);
	}
}

void snr_est(const float *image, fq_paras_float *paras, int fit)
{	/* will not change the inputted array */
	/* estimate the snr in Fourier space */
	float n = 0, noise;
	int size = paras->stamp_size;
	int i, k, edge = 1, xc = size / 2;
	int x[20]{ -1,  0,  1, -2, -1,  0,  1,  2, -2, -1,  1,  2, -2, -1,  0,  1,  2, -1,  0,  1 };
	int y[20]{ -2, -2, -2, -1, -1, -1, -1, -1,  0,  0,  0,  0,  1,  1,  1,  1,  1,  2, 2, 2 };
	float fz[20], fit_paras[6];

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
	paras->gal_flux2_ext[0] = sqrt(image[xc*size + xc]);

	if (fit == 2)
	{
		for (i = 0; i < 20; i++)
		{
			fz[i] = image[(xc + y[i])*size + xc + x[i]];
		}
		hyperfit_5(fz, fit_paras, paras);
		paras->gal_flux_alt = sqrt(pow(10, fit_paras[0]))/ noise;
		paras->gal_flux2_ext[1] = sqrt(pow(10, fit_paras[0]));
		paras->gal_flux2_ext[2] = std::max(paras->gal_flux2_ext[1], paras->gal_flux2_ext[0]);
		paras->gal_flux2_ext[3] = std::max(paras->gal_flux2, paras->gal_flux_alt);
	}
}



void possion_subtraction(double *arr, fq_paras *paras, int edge)
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

void possion_subtraction(float *arr, fq_paras_float *paras, int edge)
{
	int i, j,size = paras->stamp_size;
	float noise = 0;
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


void noise_subtraction(double *image_pow, double *noise_pow, fq_paras *paras, const int edge, const int possion)
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

void noise_subtraction(float *image_pow, float *noise_pow, fq_paras_float *paras, const int edge, const int possion)
{
	int i, j, size = paras->stamp_size;
	float inoise = 0, pnoise=0;

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



void shear_est(double *gal_img, double *psf_img, fq_paras *paras)
{	 /* will not change the inputted array */
	 /* all the inputted images are the powerspectrums */
	/* if there's no background noise, a array of '0' should be inputted */
	double mg1 = 0., mg2 = 0., mn = 0., mu = 0., mv = 0., beta, thresh, alpha, kx, kx2, ky2, ky, tk, k2, k4;
	double mp1=0., mp2=0.;
	int i, j, k, size;
	size = paras->stamp_size;

	alpha = 16*Pi*Pi*Pi*Pi/ (size*size*size*size);
	/* beta is the beta_square in the estimators */
	// use the hlr/1.414 of the PSF as the sigma of the target PSF
	// then the hlr of the target PSF will be smaller than the PSF
	// HLR = hlr/1.414*1.177 < hlr
	beta = 1./ paras->psf_hlr / paras->psf_hlr;
	
	//find the maximum of psf power spectrum and set the threshold of max/10000 above which the pixel value will be taken into account
	thresh = paras->psf_pow_thresh;

	for (i = 0; i < size; i++)//y coordinates
	{
		ky = i - size*0.5;
		for (j = 0; j < size; j++) // x coordinates
		{
			kx = j - size*0.5;
			if (psf_img[i*size + j] >= thresh)
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

void shear_est(float *gal_img, float *psf_img, fq_paras_float *paras)
{	 /* will not change the inputted array */
	 /* all the inputted images are the powerspectrums */
	/* if there's no background noise, a array of '0' should be inputted */
	float mg1 = 0., mg2 = 0., mn = 0., mu = 0., mv = 0., beta, thresh, alpha, kx, kx2, ky2, ky, tk, k2, k4;
	float mp1=0., mp2=0.;
	int i, j, k, size;
	size = paras->stamp_size;

	alpha = 16*Pi*Pi*Pi*Pi/ (size*size*size*size);
	/* beta is the beta_square in the estimators */
	// use the hlr/1.414 of the PSF as the sigma of the target PSF
	// then the hlr of the target PSF will be smaller than the PSF
	// HLR = hlr/1.414*1.177 < hlr
	beta = 1./ paras->psf_hlr / paras->psf_hlr;
	
	//find the maximum of psf power spectrum and set the threshold of max/10000 above which the pixel value will be taken into account
	thresh = paras->psf_pow_thresh;

	for (i = 0; i < size; i++)//y coordinates
	{
		ky = i - size*0.5;
		for (j = 0; j < size; j++) // x coordinates
		{
			kx = j - size*0.5;
			if (psf_img[i*size + j] >= thresh)
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

void shear_est(const double *gal_pow, const double *gal_pow_real, const double *gal_pow_imag,
               double *psf_pow,const double *psf_pow_real,const double *psf_pow_imag, fq_paras *paras)
{	 /* will not change the inputted array */
	 /* all the inputted images are the powerspectrums */
	/* if there's no background noise, a array of '0' should be inputted */
	double mg1 = 0., mg2 = 0., mn = 0., mu = 0., mv = 0., beta, thresh, alpha, kx, kx2, ky2, ky, tk, k2, k4;
	double mp1=0., mp2=0.;
	int i, j, k, size, tag;
	size = paras->stamp_size;

	alpha = 16*Pi*Pi*Pi*Pi/ (size*size*size*size);
	/* beta is the beta_square in the estimators */
	// use the hlr/1.414 of the PSF as the sigma of the target PSF
	// then the hlr of the target PSF will be smaller than the PSF
	// HLR = hlr/1.414*1.177 < hlr
	beta = 1./ paras->psf_hlr / paras->psf_hlr;
	
	//find the maximum of psf power spectrum and set the threshold of max/10000 above which the pixel value will be taken into account
	thresh = paras->psf_pow_thresh;

	for (i = 0; i < size; i++)//y coordinates
	{
		ky = i - size*0.5;
		for (j = 0; j < size; j++) // x coordinates
		{
			kx = j - size*0.5;
			tag = i*size+j;
			if (psf_pow[tag] >= thresh)
			{	
				k2 = kx*kx + ky*ky;
				k4 = k2*k2;
				kx2 = kx*kx;
				ky2 = ky*ky;			

				tk = exp(-k2*beta)*(gal_pow[tag]+4*gal_pow_real[tag]*gal_pow_imag[tag]*psf_pow_real[tag]*psf_pow_imag[tag]/psf_pow[tag])/psf_pow[tag]* alpha;
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



void ellip_est(const double *gal_img, const int size, fq_paras*paras)
{
	int i, j, sizeh = size*0.5;
	double x, y, y2, xg, q11 =0, q12=0, q22=0;
	for (i = 0; i < size; i++)
	{	
		y = i - sizeh;
		y2 = y * y;
		for (j = 0; j < size; j++)
		{
			x = j - sizeh;
			xg = x * gal_img[i*size + j];
			q11 += x * xg;
			q12 += y * xg;
			q22 += y2 * gal_img[i*size + j];
		}
	}
	paras->gal_e1 = (q11 - q22) / (q11 + q22);
	paras->gal_e2 = 2 * q12 / (q11 + q22);
}

void ellip_est(const float *gal_img, const int size, fq_paras_float*paras)
{
	int i, j, sizeh = size*0.5;
	float x, y, y2, xg, q11 =0, q12=0, q22=0;
	for (i = 0; i < size; i++)
	{	
		y = i - sizeh;
		y2 = y * y;
		for (j = 0; j < size; j++)
		{
			x = j - sizeh;
			xg = x * gal_img[i*size + j];
			q11 += x * xg;
			q12 += y * xg;
			q22 += y2 * gal_img[i*size + j];
		}
	}
	paras->gal_e1 = (q11 - q22) / (q11 + q22);
	paras->gal_e2 = 2 * q12 / (q11 + q22);
}


void block_bound(const double scale, const int ny, const int nx, double *bound_y, double *bound_x)
{
	int i, j, k;
	for (i = 0; i < ny; i++)
	{
		for (j = 0; j < nx; j++)
		{
			// the sequence of the vertexes of the blocks：
			// | (y1,x1), (y1,x2) |
			// | (y2,x1), (y2,x2) |
			bound_y[4 * (i * nx + j)] = i * scale;
			bound_y[4 * (i * nx + j) + 1] = i * scale;
			bound_y[4 * (i * nx + j) + 2] = (i + 1) * scale;
			bound_y[4 * (i * nx + j) + 3] = (i + 1) * scale;

			bound_x[4 * (i * nx + j)] = j * scale;
			bound_x[4 * (i * nx + j) + 1] = (j + 1) * scale;
			bound_x[4 * (i * nx + j) + 2] = j * scale;
			bound_x[4 * (i * nx + j) + 3] = (j + 1) * scale;
		}
	}
}

void find_block(const pts_info *infos, const double radius_s, const double radius_e, const double *bound_y, const double *bound_x, int *block_mask)
{
	int i, j, k, lb, lby, lb_d, lb_d_, seq = 0;
	int idy = infos->idy;// block id of the point
	int idx = infos->idx;
	double y = infos->y;// the coordinates of the point
	double x = infos->x;
	double scale = infos->scale;// the length of the side of the square blocks
	int ny = infos->ny; // the number of blocks along each axis 
	int nx = infos->nx;
	int num = infos->blocks_num;// the numbers of total blocks

	int nx_left, nx_right, ny_up; // the max blocks number in each direction
												// away from the point within radius_e
	int nx_s, nx_e, ny_e;

	double rs_sq = radius_s * radius_s;
	double re_sq = radius_e * radius_e;
	double dx, dy;
	double distance[4];
	// "distance" stores the distance of the four vertexes of each block,
	// the sequence of the vertexes：
	// | (y1,x1), (y1,x2) |
	// | (y2,x1), (y2,x2) |

	// find the minimum square that contains the target blocks
	nx_left = (int)((radius_e - x + bound_x[0]) / scale) + idx + 1;
	nx_right = (int)((radius_e + x - bound_x[0]) / scale) - idx;
	ny_up = (int)((radius_e + y - bound_y[0]) / scale) - idy;

	nx_s = std::max(idx - nx_left, 0);
	nx_e = std::min(idx + nx_right + 1, nx);
	ny_e = std::min(idy + ny_up + 1, ny);

	// initialiize the mask
	for (i = 0; i < num; i++)
	{
		block_mask[i] = -1;
	}

	for (i = idy; i < ny_e; i++)
	{
		lby = i * nx; // for speed
		for (j = nx_s; j < nx_e; j++)
		{
			lb = lby + j; // label of the block
			lb_d = lb * 4;
			for (k = 0; k < 4; k++)
			{
				lb_d_ = lb_d + k;
				dy = bound_y[lb_d_] - y;
				dx = bound_x[lb_d_] - x;
				distance[k] = dy * dy + dx * dx;
			}
			sort_arr(distance, 4, 1); // ascending order
			if (distance[3] < rs_sq or (distance[0] > re_sq and i not_eq idy and j not_eq idx))
			{
				//"distance[3] < rs_sq": 
				// if the max distance between the vertex and the point is smaller than radius_s,
				//	it is not the target block.
				//	"distance[0] > re_sq and i not_eq idy and j not_eq idx":  
				//  if the minimum distance between the vertexes of a block ,
				//	which is not on the cross centers on  (idy, idx) , is larger 
				//	than radius_e, it is not the target one. 
				//	The blocks on the cross centers on  (idy, idx) are the 
				//	targets.
				continue;
			}
			else
			{
				if (i > idy or (i == idy and j >= idx))
				{
					block_mask[seq] = lb;
					seq++;
				}
			}
		}
	}
}

void find_block(const pts_info *infos, const double radius, const double *bound_y, const double *bound_x, int *block_mask)
{
	int i, j, k, lb, lby, lb_d, lb_d_, seq = 0;
	int idy = infos->idy;// block id of the point
	int idx = infos->idx;
	double y = infos->y;// the coordinates of the point
	double x = infos->x;
	double scale = infos->scale;// the length of the side of the square blocks
	int ny = infos->ny; // the number of blocks along each axis 
	int nx = infos->nx;
	int num = infos->blocks_num;// the numbers of total blocks

	int nx_left, nx_right, ny_up, ny_down; // the max blocks number in each direction
												// away from the point within radius_e
	int nx_s, nx_e, ny_s, ny_e;

	double re_sq = radius * radius;
	double dx, dy;
	double distance[4];
	// "distance" stores the distance of the four vertexes of each block,
	// the sequence of the vertexes：
	// | (y1,x1), (y1,x2) |
	// | (y2,x1), (y2,x2) |

	// find the minimum square that contains the target blocks
	nx_left = (int)((radius - x + bound_x[0]) / scale) + idx + 1;
	nx_right = (int)((radius + x - bound_x[0]) / scale) - idx;
	ny_up = (int)((radius + y - bound_y[0]) / scale) - idy;
	ny_down = (int)((radius - y + bound_y[0]) / scale) + idy + 1;

	nx_s = std::max(idx - nx_left, 0);
	nx_e = std::min(idx + nx_right + 1, nx);
	ny_s = std::max(idy - ny_down, 0);
	ny_e = std::min(idy + ny_up + 1, ny);

	// initialiize the mask
	for (i = 0; i < num; i++)
	{
		block_mask[i] = -1;
	}
	//std::cout<<nx_s<<" "<<nx_e<<" "<<ny_s<<" "<<ny_e<<std::endl;
	for (i = ny_s; i < ny_e; i++)
	{
		lby = i * nx; // for speed
		for (j = nx_s; j < nx_e; j++)
		{
			lb = lby + j; // label of the block
			lb_d = lb * 4;
			for (k = 0; k < 4; k++)
			{
				lb_d_ = lb_d + k;
				dy = bound_y[lb_d_] - y;
				dx = bound_x[lb_d_] - x;
				distance[k] = dy * dy + dx * dx;
			}
			sort_arr(distance, 4, 1); // ascending order
			if (distance[0] > re_sq and i not_eq idy and j not_eq idx)
			{
				//	"distance[0] > re_sq and i not_eq idy and j not_eq idx":  
				//  if the minimum distance between the vertexes of a block ,
				//	which is not on the cross centers on  (idy, idx) , is larger 
				//	than radius_e, it is not the target one. 
				//	The blocks on the cross centers on  (idy, idx) are the 
				//	targets.
				continue;
			}
			else
			{
				block_mask[seq] = lb;
				seq++;
			}
		}
	}
}



void chisq_Gbin_1d(const double *mg, const double *mnu, const int data_num, const double *bins, const int bin_num, const double gh, double &result)
{
	int i, j, k;
	double *temp = new double[data_num];
	int *num_in_bin = new int[bin_num];

	for (i = 0; i < data_num; i++)
	{
		temp[i] = mg[i] - gh * mnu[i];
	}
	histogram(temp, bins, num_in_bin, data_num, bin_num);
	try
	{
		cal_chisq_1d(num_in_bin, bin_num, result);
	}
	catch(const char *msg)
	{	
		std::cout << "g_guess: " << gh << std::endl;
		std::cout << "Num: " << std::endl;
		show_arr(num_in_bin, 1, bin_num);
		std::cout << "Bin: " << std::endl;
		show_arr(bins, 1, bin_num + 1);
		char err_inform[200];
		sprintf(err_inform,"%s, Chi square divided by zero (chisq_Gbin_1d -> cal_chisq_1d)!!!",msg);
		throw err_inform;
	}

	delete[] num_in_bin;
	delete[] temp;
}

void chisq_Gbin_1d(const float *mg, const float *mnu, const int data_num, const float *bins, const int bin_num, const float gh, float &result)
{
	int i, j, k;
	float *temp = new float[data_num];
	int *num_in_bin = new int[bin_num];

	for (i = 0; i < data_num; i++)
	{
		temp[i] = mg[i] - gh * mnu[i];
	}
	histogram(temp, bins, num_in_bin, data_num, bin_num);
	try
	{
		cal_chisq_1d(num_in_bin, bin_num, result);
	}
	catch(const char *msg)
	{	
		std::cout << "g_guess: " << gh << std::endl;
		std::cout << "Num: " << std::endl;
		show_arr(num_in_bin, 1, bin_num);
		std::cout << "Bin: " << std::endl;
		show_arr(bins, 1, bin_num + 1);
		char err_inform[200];
		sprintf(err_inform,"%s, Chi square divided by zero (chisq_Gbin_1d -> cal_chisq_1d)!!!",msg);
		throw err_inform;
	}

	delete[] num_in_bin;
	delete[] temp;
}


void chisq_Gbin_1d(const double *mg, const int data_num, const double *bins, const int bin_num, double &result)
{
	int i, j, k;

	int *num_in_bin = new int[bin_num];

	histogram(mg, bins, num_in_bin, data_num, bin_num);
	try
	{
		cal_chisq_1d(num_in_bin, bin_num, result);
	}
	catch(const char *msg)
	{	
		std::cout << "Num: " << std::endl;
		show_arr(num_in_bin, 1, bin_num);
		std::cout << "Bin: " << std::endl;
		show_arr(bins, 1, bin_num + 1);
		char err_inform[200];
		sprintf(err_inform,"%s, Chi square divided by zero (chisq_Gbin_1d -> cal_chisq_1d)!!!",msg);
		throw err_inform;
	}
	delete[] num_in_bin;
}


void cal_chisq_2d(const double *hist_arr, const int size, double &result)
{
	int h = size / 2, i, j, s1, s2 = size * size;
	double chi = 0, n, m;
	for (i = 0; i < h; i++)
	{
		s1 = i * size;
		for (j = 0; j < h; j++)
		{
			m = hist_arr[s1 + j] + hist_arr[s2 - s1 - j - 1] - (hist_arr[s1 + size - j - 1] + hist_arr[s2 - s1 - size + j]);
			n = hist_arr[s1 + j] + hist_arr[s2 - s1 - j - 1] + hist_arr[s1 + size - j - 1] + hist_arr[s2 - s1 - size + j];
			if (hist_arr[s1 + j] + hist_arr[s2 - s1 - j - 1] + hist_arr[s1 + size - j - 1] + hist_arr[s2 - s1 - size + j] == 0)
			{
				std::cout << "Chi square divided by zero!!!" << std::endl;
				exit(0);
			}
			chi += m * m / n;
		}
	}
	result = chi * 0.5;

}

void cal_chisq_2d(const long *hist_arr, const int size, double &result)
{
	int h = size / 2, i, j, s1, s2 = size * size;
	double chi = 0, n, m;

	for (i = 0; i < h; i++)
	{
		s1 = i * size;
		for (j = 0; j < h; j++)
		{
			m = hist_arr[s1 + j] + hist_arr[s2 - s1 - j - 1] - (hist_arr[s1 + size - j - 1] + hist_arr[s2 - s1 - size + j]);
			n = hist_arr[s1 + j] + hist_arr[s2 - s1 - j - 1] + hist_arr[s1 + size - j - 1] + hist_arr[s2 - s1 - size + j];
			if (hist_arr[s1 + j] + hist_arr[s2 - s1 - j - 1] + hist_arr[s1 + size - j - 1] + hist_arr[s2 - s1 - size + j] == 0)
			{
				std::cout << "Chi square divided by zero!!!" << std::endl;
				exit(0);
			}
			chi += m * m / n;
		}
	}
	result = chi * 0.5;
}

void cal_chisq_2d(const int *hist_arr, const int size, double &result)
{
	int h = size / 2, i, j, s1, s2 = size * size;
	double chi = 0, n, m;

	for (i = 0; i < h; i++)
	{
		s1 = i * size;
		for (j = 0; j < h; j++)
		{
			m = hist_arr[s1 + j] + hist_arr[s2 - s1 - j - 1] - (hist_arr[s1 + size - j - 1] + hist_arr[s2 - s1 - size + j]);
			n = hist_arr[s1 + j] + hist_arr[s2 - s1 - j - 1] + hist_arr[s1 + size - j - 1] + hist_arr[s2 - s1 - size + j];
			if (hist_arr[s1 + j] + hist_arr[s2 - s1 - j - 1] + hist_arr[s1 + size - j - 1] + hist_arr[s2 - s1 - size + j] == 0)
			{
				std::cout << "Chi square divided by zero!!!" << std::endl;
				exit(0);
			}
			chi += m * m / n;
		}
	}
	result = chi * 0.5;
}


void cal_chisq_1d(const double *hist_num, const int bin_num, double &result)
{
	// the size must be an even number
	int i, j;
	int mid = bin_num / 2;
	double chi_count = 0;
	double dn, sn;
	for (i = mid; i < bin_num; i++)
	{
		dn = hist_num[i] - hist_num[bin_num - i - 1];
		sn = hist_num[i] + hist_num[bin_num - i - 1];
		if (hist_num[i] + hist_num[bin_num - i - 1] == 0)
		{
			std::cout << "ERROR in cal_chisq_1d !"<<std::endl;
			std::cout << "ERROR Data : ";
			show_arr(hist_num, 1, bin_num);
			throw "ERROR chi square divided by zero !!!";
		}
		chi_count += dn * dn / sn;
	}
	result = chi_count * 0.5;
}

void cal_chisq_1d(const long *hist_num, const int bin_num, double &result)
{
	// the size must be an even number
	int i, j;
	int mid = bin_num / 2;
	double chi_count = 0;
	double dn, sn;
	for (i = mid; i < bin_num; i++)
	{
		dn = hist_num[i] - hist_num[bin_num - i - 1];
		sn = hist_num[i] + hist_num[bin_num - i - 1];
		if (hist_num[i] + hist_num[bin_num - i - 1] == 0)
		{
			std::cout << "ERROR in cal_chisq_1d !" << std::endl;
			std::cout << "ERROR Data : ";
			show_arr(hist_num, 1, bin_num);
			throw "ERROR chi square divided by zero !!!";
		}
		chi_count += dn * dn / sn;
	}
	result = chi_count * 0.5;
}

void cal_chisq_1d(const int *hist_num, const int bin_num, double &result)
{
	// the size must be an even number
	int i, j;
	int mid = bin_num / 2;
	double chi_count = 0;
	double dn, sn;
	for (i = mid; i < bin_num; i++)
	{
		dn = hist_num[i] - hist_num[bin_num - i - 1];
		sn = hist_num[i] + hist_num[bin_num - i - 1];
		if (hist_num[i] + hist_num[bin_num - i - 1] == 0)
		{
			std::cout << "ERROR in cal_chisq_1d !" << std::endl;
			std::cout<<"ERROR Data : ";
			show_arr(hist_num, 1, bin_num);
			throw "ERROR chi square divided by zero !!!";
		}
		chi_count += dn * dn / sn;
	}
	result = chi_count * 0.5;
}

void cal_chisq_1d(const int *hist_num, const int bin_num, float &result)
{
	// the size must be an even number
	int i, j;
	int mid = bin_num / 2;
	float chi_count = 0;
	float dn, sn;
	for (i = mid; i < bin_num; i++)
	{
		dn = hist_num[i] - hist_num[bin_num - i - 1];
		sn = hist_num[i] + hist_num[bin_num - i - 1];
		if (hist_num[i] + hist_num[bin_num - i - 1] == 0)
		{
			std::cout << "ERROR in cal_chisq_1d !" << std::endl;
			std::cout<<"ERROR Data : ";
			show_arr(hist_num, 1, bin_num);
			throw "ERROR chi square divided by zero !!!";
		}
		chi_count += dn * dn / sn;
	}
	result = chi_count * 0.5;
}

void cal_chisq_1d(const int *hist_num, const int bin_num, const int num, double &result)
{
	// the size must be an even number
	int i, j;
	int mid = bin_num / 2;
	double chi_count = 0;
	double dn, sn;
	for (i = mid; i < mid+num+1; i++)
	{
		dn = hist_num[i] - hist_num[bin_num - i - 1];
		sn = hist_num[i] + hist_num[bin_num - i - 1];
		if (hist_num[i] + hist_num[bin_num - i - 1] == 0)
		{
			std::cout << "ERROR in cal_chisq_1d !";
			std::cout << "ERROR Data : ";
			show_arr(hist_num, 1, bin_num);
			throw "ERROR chi square divided by zero !!!";
		}
		chi_count += dn * dn / sn;
	}
	result = chi_count * 0.5;
}



void find_shear_mean(const double *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col,
						 double &gh, double &gh_sig, const int sub_block_num, const double scale)
{
	// if the data array is very large, say > 10^7\, then summing it directly may cause numerical problem
	// it's better to sum the sub-block and then add the these quantities together.

	int i,j, m,n;

	double mg_mean, mg_sq_mean, mn_mean, mg_sum, mn_sum;
	double sub_sum1, sub_sum2, sub_sum3;;

	double s1= 1./scale;
	
	int *num_in_block = new int[sub_block_num];
	int *block_st = new int[sub_block_num];

	m = data_row/sub_block_num;
	n = data_row%sub_block_num;

	for(i=0;i<sub_block_num;i++)
	{
		block_st[i] = 0;

		num_in_block[i] = m;
		if(i<n)num_in_block[i] += 1;

		for(j=0;j<i;j++)
		{
			block_st[i] += num_in_block[j];
		}
	}

	mg_mean = 0;
	mg_sq_mean = 0;
	mn_mean = 0;
	for(i=0;i<sub_block_num;i++)
	{
		sub_sum1 = 0;
		sub_sum2 = 0;
		sub_sum3 = 0;

		for(j=block_st[i]; j<block_st[i]+num_in_block[i]; j++)
		{	
			m = j*data_col;
			sub_sum1 += data_arr[m + mg_col]*s1;
			sub_sum2 += data_arr[m + mn_col]*s1;
			sub_sum3 += (data_arr[m + mg_col]*s1)*(data_arr[m + mg_col]*s1);
		}
		mg_mean += sub_sum1*(scale/data_row);
		mn_mean += sub_sum2*(scale/data_row);
		mg_sq_mean += sub_sum3*(scale*scale/data_row);
	}

	gh = mg_mean/mn_mean;
	gh_sig = sqrt(mg_sq_mean/mn_mean/mn_mean/data_row);

	delete[] num_in_block;
	delete[] block_st;
}

void find_shear_mean(const float *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col, 
						float &gh, float &gh_sig, const int sub_block_num, const float scale)
{
	// if the data array is very large, say > 10^7\, then summing it directly may cause numerical problem
	// it's better to sum the sub-block and then add the these quantities together.

	int i,j, m,n;

	float mg_mean, mg_sq_mean, mn_mean, mg_sum, mn_sum;
	float sub_sum1, sub_sum2, sub_sum3;;

	float s1= 1./scale;
	
	int *num_in_block = new int[sub_block_num];
	int *block_st = new int[sub_block_num];

	m = data_row/sub_block_num;
	n = data_row%sub_block_num;

	for(i=0;i<sub_block_num;i++)
	{
		block_st[i] = 0;

		num_in_block[i] = m;
		if(i<n)num_in_block[i] += 1;

		for(j=0;j<i;j++)
		{
			block_st[i] += num_in_block[j];
		}
	}

	mg_mean = 0;
	mg_sq_mean = 0;
	mn_mean = 0;
	for(i=0;i<sub_block_num;i++)
	{
		sub_sum1 = 0;
		sub_sum2 = 0;
		sub_sum3 = 0;

		for(j=block_st[i];j<block_st[i]+num_in_block[i];j++)
		{	
			m = j*data_col;
			sub_sum1 += data_arr[m+mg_col]*s1;
			sub_sum2 += data_arr[m+mn_col]*s1;
			sub_sum3 += (data_arr[m+mg_col]*s1)*(data_arr[m+mg_col]*s1);
		}
		mg_mean += sub_sum1*(scale/data_row);
		mn_mean += sub_sum2*(scale/data_row);
		mg_sq_mean += sub_sum3*(scale*scale/data_row);
	}

	gh = mg_mean/mn_mean;
	gh_sig = sqrt(mg_sq_mean/mn_mean/mn_mean/data_row);

	delete[] num_in_block;
	delete[] block_st;
}


void find_shear_mean(const double *mg, const double *mn, const int data_row, double &gh, double &gh_sig, const int sub_block_num, const double scale)
{
	// if the data array is very large, say > 10^7\, then summing it directly may cause numerical problem
	// it's better to sum the sub-block and then add the these quantities together.

	int i,j, m,n;

	double mg_mean, mg_sq_mean, mn_mean, mg_sum, mn_sum;
	double sub_sum1, sub_sum2, sub_sum3;;

	double s1= 1./scale;
	
	int *num_in_block = new int[sub_block_num];
	int *block_st = new int[sub_block_num];

	m = data_row/sub_block_num;
	n = data_row%sub_block_num;

	for(i=0;i<sub_block_num;i++)
	{
		block_st[i] = 0;

		num_in_block[i] = m;
		if(i<n)num_in_block[i] += 1;

		for(j=0;j<i;j++)
		{
			block_st[i] += num_in_block[j];
		}
	}

	mg_mean = 0;
	mg_sq_mean = 0;
	mn_mean = 0;
	for(i=0;i<sub_block_num;i++)
	{
		sub_sum1 = 0;
		sub_sum2 = 0;
		sub_sum3 = 0;

		for(j=block_st[i]; j<block_st[i]+num_in_block[i]; j++)
		{	
			sub_sum1 += mg[j]*s1;
			sub_sum2 += mn[j]*s1;
			sub_sum3 += (mg[j]*s1)*(mg[j]*s1);
		}
		mg_mean += sub_sum1*(scale/data_row);
		mn_mean += sub_sum2*(scale/data_row);
		mg_sq_mean += sub_sum3*(scale*scale/data_row);
	}

	gh = mg_mean/mn_mean;
	gh_sig = sqrt(mg_sq_mean/mn_mean/mn_mean/data_row);

	delete[] num_in_block;
	delete[] block_st;
}

void find_shear_mean(const float *mg, const float *mn, const int data_row, float &gh, float &gh_sig, const int sub_block_num, const float scale)
{
	// if the data array is very large, say > 10^7\, then summing it directly may cause numerical problem
	// it's better to sum the sub-block and then add the these quantities together.

	int i,j, m,n;

	float mg_mean, mg_sq_mean, mn_mean, mg_sum, mn_sum;
	float sub_sum1, sub_sum2, sub_sum3;;

	float s1= 1./scale;
	
	int *num_in_block = new int[sub_block_num];
	int *block_st = new int[sub_block_num];

	m = data_row/sub_block_num;
	n = data_row%sub_block_num;

	for(i=0;i<sub_block_num;i++)
	{
		block_st[i] = 0;

		num_in_block[i] = m;
		if(i<n)num_in_block[i] += 1;

		for(j=0;j<i;j++)
		{
			block_st[i] += num_in_block[j];
		}
	}

	mg_mean = 0;
	mg_sq_mean = 0;
	mn_mean = 0;
	for(i=0;i<sub_block_num;i++)
	{
		sub_sum1 = 0;
		sub_sum2 = 0;
		sub_sum3 = 0;

		for(j=block_st[i]; j<block_st[i]+num_in_block[i]; j++)
		{	
			sub_sum1 += mg[j]*s1;
			sub_sum2 += mn[j]*s1;
			sub_sum3 += (mg[j]*s1)*(mg[j]*s1);
		}
		mg_mean += sub_sum1*(scale/data_row);
		mn_mean += sub_sum2*(scale/data_row);
		mg_sq_mean += sub_sum3*(scale*scale/data_row);
	}

	gh = mg_mean/mn_mean;
	gh_sig = sqrt(mg_sq_mean/mn_mean/mn_mean/data_row);

	delete[] num_in_block;
	delete[] block_st;
}


void find_shear(const double *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col, const int mu_col, const int bin_num, 
				int g_label, double &gh, double &gh_sig, double &chisq_min_fit, double *chi_check, 
				const int chi_fit_num, const int choice, const double max_scale, const double ini_left, const double ini_right,const double chi_gap)
{
	int i, j;
	int max_iters = 12;
	int temp_num;

	int *num_in_bin = new int[bin_num];
	double *bins = new double[bin_num + 1];
	double *gh_fit = new double[chi_fit_num];
	double *chisq_fit = new double[chi_fit_num];
	double chisq_fit_coeff[3];

	// record the each g_left, chisq_left, g_right, chisq_right
	int record_col = 4;
	double fit_max_chisq, new_end;
	// double *search_vals = new double[(max_iters+1)*record_col]{};

	int same = 0, iters = 0, change = 1;
	double left = ini_left, right = ini_right, step;
	double chi_left, chi_right, chi_mid;
	double gh_left, gh_right, gh_mid;
	
	temp_num = data_row;
	if(choice > 0)
	{
		temp_num = choice;
	}

	double *mg = new double[temp_num];
	int ch_step = data_row / temp_num;
	for(i=0;i<temp_num;i++){mg[i]=data_arr[i*ch_step*data_col + mg_col];}
	
	set_bin(mg, temp_num, bins, bin_num, max_scale);
	delete[] mg;
	//double st1, st2, st3, st4, st5, st6;
	//st1 = clock();
	// set the bins for G1(2)
	
	// show_arr(bins, 1, bin_num + 1);
	//st2 = clock();
	while (change == 1)
	{		
		change = 0;
		gh_mid = (left + right) *0.5;
		gh_left = left;
		gh_right = right;

		try
		{
			fourier_hist(data_arr, data_row, data_col, mg_col, mn_col, mu_col, gh_left, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_left);

			fourier_hist(data_arr, data_row, data_col, mg_col, mn_col, mu_col, gh_mid, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_mid);

			fourier_hist(data_arr, data_row, data_col, mg_col, mn_col, mu_col, gh_right, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		// search_vals[iters*record_col] = gh_left;
		// search_vals[iters*record_col+1] = chi_left;
		// search_vals[iters*record_col+2] = gh_right;
		// search_vals[iters*record_col+3] = chi_right;

		//std::cout << left << " "<< gh_left<<" "<< gh_mid<<" "<< gh_right <<" "<< right << std::endl;

		if (chi_left > chi_mid + chi_gap)
		{
			left = (gh_mid + gh_left) *0.5;
			change = 1;
		}
		if (chi_right > chi_mid + chi_gap)
		{
			right = (gh_mid + gh_right)*0.5;
			change = 1;
		}

		iters += 1;
		if (iters > max_iters)
		{
			break;
		}
	}
	
	// fit_max_chisq = std::max(chi_left, chi_right)+chi_gap;
	
	temp_num = 7;
	step = (right - left) / (temp_num - 1);
	for (i = 0; i < temp_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < temp_num; i++)
	{	
		try
		{
			fourier_hist(data_arr, data_row, data_col, mg_col, mn_col, mu_col, gh_fit[i], g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
	}
	//st4 = clock();
	fit_shear(gh_fit, chisq_fit, temp_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff, -1);

	// to get a more symmetrical interval for fitting
	new_end = std::max(gh-left, right-gh);
	
	left = gh - new_end;
	step = 2*new_end / (chi_fit_num - 1);

	for (i = 0; i < chi_fit_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < chi_fit_num; i++)
	{	
		try
		{
			fourier_hist(data_arr, data_row, data_col, mg_col, mn_col, mu_col, gh_fit[i], g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
		if (chi_check)
		{	// for checking
			chi_check[i] = gh_fit[i];
			chi_check[chi_fit_num + i] = chi_right;
		}
	}
	
	//st4 = clock();
	fit_shear(gh_fit, chisq_fit, chi_fit_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff, -1);

	//st5 = clock();
	//std::cout << gh << " " << gh_sig << std::endl;
	//std::cout <<"Time: "<< (st2 - st1) / CLOCKS_PER_SEC << " " << (st3 - st2) / CLOCKS_PER_SEC << " " << (st4 - st3) / CLOCKS_PER_SEC << " " << (st5 - st4) / CLOCKS_PER_SEC << std::endl;
	delete[] gh_fit;
	delete[] chisq_fit;
	delete[] bins;
	// delete[] search_vals;
	delete[] num_in_bin;
}

void find_shear(const float *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col, const int mu_col, const int bin_num, 
				const int g_label, float &gh, float &gh_sig, float &chisq_min_fit, float *chi_check,	const int chi_fit_num, 
				const int choice, const float max_scale, const float ini_left, const float ini_right, const float chi_gap)
{
	int i, j;
	int max_iters = 12;
	int temp_num;

	int *num_in_bin = new int[bin_num];

	float *bins = new float[bin_num + 1];

	float *gh_fit = new float[chi_fit_num];
	float *chisq_fit = new float[chi_fit_num];
	float chisq_fit_coeff[3];

	// record the each g_left, chisq_left, g_right, chisq_right
	int record_col = 4;
	int left_tag=-1, right_tag=-1;
	float fit_max_chisq, new_end;
	// float *search_vals = new float[(max_iters+1)*record_col]{};

	int same = 0, iters = 0, change = 1;
	float left = ini_left, right = ini_right, step;
	float chi_left, chi_right, chi_mid;
	float gh_left, gh_right, gh_mid;

	//double st1, st2, st3, st4, st5, st6;
	//st1 = clock();
	// set the bins for G1(2)
	temp_num = data_row;
	if(choice > 0)
	{
		temp_num = choice;
	}

	float *mg = new float[temp_num];
	int ch_step = data_row / temp_num;
	for(i=0;i<temp_num;i++){mg[i]=data_arr[i*ch_step*data_col + mg_col];}
	
	set_bin(mg, temp_num, bins, bin_num, max_scale);
	delete[] mg;
	// show_arr(bins, 1, bin_num + 1);
	//st2 = clock();
	while (change == 1)
	{		
		change = 0;
		gh_mid = (left + right) *0.5;
		gh_left = left;
		gh_right = right;

		try
		{	
			fourier_hist(data_arr, data_row, data_col, mg_col, mn_col, mu_col, gh_left, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_left);

			fourier_hist(data_arr, data_row, data_col, mg_col, mn_col, mu_col, gh_mid, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_mid);

			fourier_hist(data_arr, data_row, data_col, mg_col, mn_col, mu_col, gh_right, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		// search_vals[iters*record_col] = gh_left;
		// search_vals[iters*record_col+1] = chi_left;
		// search_vals[iters*record_col+2] = gh_right;
		// search_vals[iters*record_col+3] = chi_right;

		//std::cout << left << " "<< gh_left<<" "<< gh_mid<<" "<< gh_right <<" "<< right << std::endl;

		if (chi_left > chi_mid + chi_gap)
		{
			left = (gh_mid + gh_left) *0.5;
			change = 1;
		}
		if (chi_right > chi_mid + chi_gap)
		{
			right = (gh_mid + gh_right)*0.5;
			change = 1;
		}

		iters += 1;
		if (iters > max_iters)
		{
			break;
		}
	}
	
	// fit_max_chisq = std::max(chi_left, chi_right)+chi_gap;
	// std::cout<<right<<" "<<left<<std::endl;
	// std::cout<<"chi max "<<fit_max_chisq<<std::endl;

	temp_num = 7;
	step = (right - left) / (temp_num - 1);
	for (i = 0; i < temp_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < temp_num; i++)
	{	
		try
		{
			fourier_hist(data_arr, data_row, data_col, mg_col, mn_col, mu_col, gh_fit[i], g_label, bins, num_in_bin, bin_num);

			// std::cout<<i<<" "<<gh_fit[i]<<std::endl;
			// show_arr(num_in_bin, 1, bin_num);

			cal_chisq_1d(num_in_bin, bin_num, chi_right);
			// std::cout<<chi_right<<std::endl;
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
	}
	//st4 = clock();

	fit_shear(gh_fit, chisq_fit, temp_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff, -1);
	// std::cout<<gh<<" "<<gh_sig<<std::endl;

	// to get a more symmetrical interval for fitting
	new_end = std::max(gh-left, right-gh);
	
	left = gh - new_end;
	step = 2*new_end / (chi_fit_num - 1);

	for (i = 0; i < chi_fit_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < chi_fit_num; i++)
	{	
		try
		{
			fourier_hist(data_arr, data_row, data_col, mg_col, mn_col, mu_col, gh_fit[i], g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
		if (chi_check)
		{	// for checking
			chi_check[i] = chi_right;
			chi_check[chi_fit_num + i] = gh_fit[i];
		}
	}
	
	//st4 = clock();
	fit_shear(gh_fit, chisq_fit, chi_fit_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff,-1);

	//st5 = clock();
	//std::cout << gh << " " << gh_sig << std::endl;
	//std::cout <<"Time: "<< (st2 - st1) / CLOCKS_PER_SEC << " " << (st3 - st2) / CLOCKS_PER_SEC << " " << (st4 - st3) / CLOCKS_PER_SEC << " " << (st5 - st4) / CLOCKS_PER_SEC << std::endl;
	delete[] gh_fit;
	delete[] chisq_fit;
	delete[] bins;
	// delete[] search_vals;
	delete[] num_in_bin;
}



void find_shear(const double *mg, const double *mn, const double *mu, const int data_num, const int bin_num, const double *bins, int g_label, double &gh, 
				double &gh_sig, double &chisq_min_fit, double *chi_check, 
				const int chi_fit_num, const double ini_left, const double ini_right, const double chi_gap)
{
	int i, j, k;
	int max_iters = 12;
	int temp_num;

	int *num_in_bin = new int[bin_num];
	double *gh_fit = new double[chi_fit_num];
	double *chisq_fit = new double[chi_fit_num];
	double chisq_fit_coeff[3];

	// record the each g_left, chisq_left, g_right, chisq_right
	int record_col = 4;
	int left_tag=-1, right_tag=-1;
	double fit_max_chisq, new_end;
	

	int same = 0, iters = 0, change = 1;
	double left = ini_left, right = ini_right, step;
	double chi_left, chi_right, chi_mid;
	double gh_left, gh_right, gh_mid;

	//double st1, st2, st3, st4, st5, st6;
	//st1 = clock();
	//st2 = clock();
	while (change == 1)
	{		
		change = 0;
		gh_mid = (left + right) *0.5;
		gh_left = left;
		gh_right = right;

		try
		{
			fourier_hist(mg, mn, mu, data_num, gh_left, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_left);

			fourier_hist(mg, mn, mu, data_num, gh_mid, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_mid);

			fourier_hist(mg, mn, mu, data_num, gh_right, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		//std::cout << left << " "<< gh_left<<" "<< gh_mid<<" "<< gh_right <<" "<< right << std::endl;

		if (chi_left > chi_mid + chi_gap)
		{
			left = (gh_mid + gh_left) *0.5;
			change = 1;
		}
		if (chi_right > chi_mid + chi_gap)
		{
			right = (gh_mid + gh_right)*0.5;
			change = 1;
		}

		iters += 1;
		if (iters > max_iters)
		{
			break;
		}
	}
	
	
	temp_num = 7;
	step = (right - left) / (temp_num - 1);
	for (i = 0; i < temp_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < temp_num; i++)
	{	
		try
		{
			fourier_hist(mg, mn, mu, data_num, gh_fit[i], g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
	}
	//st4 = clock();
	fit_shear(gh_fit, chisq_fit, temp_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff,-1);

	// to get a more symmetrical interval for fitting
	new_end = std::max(gh-left, right-gh);
	
	left = gh - new_end;
	step = 2*new_end / (chi_fit_num - 1);

	for (i = 0; i < chi_fit_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < chi_fit_num; i++)
	{	
		try
		{
			fourier_hist(mg, mn, mu, data_num, gh_fit[i], g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
		if (chi_check)
		{	// for checking
			chi_check[i] = chi_right;
			chi_check[chi_fit_num + i] = gh_fit[i];
		}
	}
	
	//st4 = clock();
	fit_shear(gh_fit, chisq_fit, chi_fit_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff, -1);

	//st5 = clock();
	//std::cout << gh << " " << gh_sig << std::endl;
	//std::cout <<"Time: "<< (st2 - st1) / CLOCKS_PER_SEC << " " << (st3 - st2) / CLOCKS_PER_SEC << " " << (st4 - st3) / CLOCKS_PER_SEC << " " << (st5 - st4) / CLOCKS_PER_SEC << std::endl;
	delete[] gh_fit;
	delete[] chisq_fit;
	delete[] num_in_bin;
}

void find_shear(const float *mg, const float *mn,  const float *mu, const int data_num, const int bin_num, const float *bins, int g_label, 
				float &gh, float &gh_sig, float &chisq_min_fit, float *chi_check,
				const int chi_fit_num, const float ini_left, const float ini_right, const float chi_gap)
{
	int i, j, k;
	int max_iters = 12;
	int temp_num;

	int *num_in_bin = new int[bin_num];
	float *gh_fit = new float[chi_fit_num];
	float *chisq_fit = new float[chi_fit_num];
	float chisq_fit_coeff[3];

	// record the each g_left, chisq_left, g_right, chisq_right
	int record_col = 4;
	int left_tag=-1, right_tag=-1;
	float fit_max_chisq, new_end;


	int same = 0, iters = 0, change = 1;
	float left = ini_left, right = ini_right, step;
	float chi_left, chi_right, chi_mid;
	float gh_left, gh_right, gh_mid;

	//double st1, st2, st3, st4, st5, st6;
	//st1 = clock();
	// set the bins for G1(2)

	//show_arr(bins, 1, bin_num + 1);
	//st2 = clock();
	while (change == 1)
	{		
		change = 0;
		gh_mid = (left + right) *0.5;
		gh_left = left;
		gh_right = right;

		try
		{	
			fourier_hist(mg, mn, mu, data_num, gh_left, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_left);

			fourier_hist(mg, mn, mu, data_num, gh_mid, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_mid);

			fourier_hist(mg, mn, mu, data_num, gh_right, g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);

		}
		catch(const char *msg)
		{
			throw msg;
		}

		//std::cout << left << " "<< gh_left<<" "<< gh_mid<<" "<< gh_right <<" "<< right << std::endl;

		if (chi_left > chi_mid + chi_gap)
		{
			left = (gh_mid + gh_left) *0.5;
			change = 1;
		}
		if (chi_right > chi_mid + chi_gap)
		{
			right = (gh_mid + gh_right)*0.5;
			change = 1;
		}

		iters += 1;
		if (iters > max_iters)
		{
			break;
		}
	}
	
	temp_num = 7;
	step = (right - left) / (temp_num - 1);
	for (i = 0; i < temp_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < temp_num; i++)
	{	
		try
		{
			fourier_hist(mg, mn, mu, data_num, gh_fit[i], g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
	}
	//st4 = clock();

	fit_shear(gh_fit, chisq_fit, temp_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff, -1);

	// to get a more symmetrical interval for fitting
	new_end = std::max(gh-left, right-gh);
	
	left = gh - new_end;
	step = 2*new_end / (chi_fit_num - 1);

	for (i = 0; i < chi_fit_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < chi_fit_num; i++)
	{	
		try
		{
			fourier_hist(mg, mn, mu, data_num, gh_fit[i], g_label, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
		if (chi_check)
		{	// for checking
			chi_check[i] = chi_right;
			chi_check[chi_fit_num + i] = gh_fit[i];
		}
	}
	
	//st4 = clock();
	fit_shear(gh_fit, chisq_fit, chi_fit_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff, -1);

	//st5 = clock();
	//std::cout << gh << " " << gh_sig << std::endl;
	//std::cout <<"Time: "<< (st2 - st1) / CLOCKS_PER_SEC << " " << (st3 - st2) / CLOCKS_PER_SEC << " " << (st4 - st3) / CLOCKS_PER_SEC << " " << (st5 - st4) / CLOCKS_PER_SEC << std::endl;
	delete[] gh_fit;
	delete[] chisq_fit;
	delete[] num_in_bin;
}



void find_shear_NU(const double *mg, const double *mn, const double *mu, const int data_num, const int bin_num, const double *bins, int NU_label, 
					double fix_g, double &gh, double &gh_sig, double &chisq_min_fit, double *chi_check, 
					const int chi_fit_num, const double ini_left, const double ini_right, const double chi_gap)
{
	int i, j, k;
	int max_iters = 12;
	int temp_num;

	int *num_in_bin = new int[bin_num];
	double *gh_fit = new double[chi_fit_num];
	double *chisq_fit = new double[chi_fit_num];
	double chisq_fit_coeff[3];

	// record the each g_left, chisq_left, g_right, chisq_right
	int record_col = 4;
	int left_tag=-1, right_tag=-1;
	double fit_max_chisq, new_end;
	

	int same = 0, iters = 0, change = 1;
	double left = ini_left, right = ini_right, step;
	double chi_left, chi_right, chi_mid;
	double gh_left, gh_right, gh_mid;
	double a1, a2,a3, b1, b2,b3;
	//double st1, st2, st3, st4, st5, st6;
	//st1 = clock();
	//st2 = clock();
	while (change == 1)
	{		
		change = 0;
		gh_mid = (left + right) *0.5;
		gh_left = left;
		gh_right = right;
		if(NU_label == 1)
		{
			a1 = fix_g;	a2 = fix_g;	a3 = fix_g;
			b1 = gh_left;  b2 = gh_mid;  b3 = gh_right;
		}
		else
		{
			b1 = fix_g;	b2 = fix_g;	b3 = fix_g;
			a1 = gh_left; a2 = gh_mid; a3 = gh_right;
		}
		try
		{	
			fourier_hist_NU(mg, mn, mu, data_num, a1, b1, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_left);
			// std::cout<<a1<<" "<<b1<<" "<<chi_left<<std::endl;
			// show_arr(num_in_bin, 1, bin_num);

			fourier_hist_NU(mg, mn, mu, data_num, a2, b2, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_mid);
			// std::cout<<a2<<" "<<b2<<" "<<chi_left<<std::endl;
			// show_arr(num_in_bin, 1, bin_num);

			fourier_hist_NU(mg, mn, mu, data_num, a3, b3, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
			// std::cout<<a3<<" "<<b3<<" "<<chi_left<<std::endl;
			// show_arr(num_in_bin, 1, bin_num);

		}
		catch(const char *msg)
		{
			throw msg;
		}
		//std::cout << left << " "<< gh_left<<" "<< gh_mid<<" "<< gh_right <<" "<< right << std::endl;

		if (chi_left > chi_mid + chi_gap)
		{
			left = (gh_mid + gh_left) *0.5;
			change = 1;
		}
		if (chi_right > chi_mid + chi_gap)
		{
			right = (gh_mid + gh_right)*0.5;
			change = 1;
		}

		iters += 1;
		if (iters > max_iters)
		{
			break;
		}
	}
	
	temp_num = 7;
	step = (right - left) / (temp_num - 1);
	for (i = 0; i < temp_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < temp_num; i++)
	{	
		if(NU_label == 1)
		{
			a1 = fix_g;
			b1 = gh_fit[i];
		}
		else
		{
			b1 = fix_g;
			a1 = gh_fit[i];
		}

		try
		{
			fourier_hist_NU(mg, mn, mu, data_num, a1, b1, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
	}
	//st4 = clock();
	fit_shear(gh_fit, chisq_fit, temp_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff, -1);

	// to get a more symmetrical interval for fitting
	new_end = std::max(gh-left, right-gh);
	
	left = gh - new_end;
	step = 2*new_end / (chi_fit_num - 1);

	for (i = 0; i < chi_fit_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < chi_fit_num; i++)
	{	
		if(NU_label == 1)
		{
			a1 = fix_g;
			b1 = gh_fit[i];
		}
		else
		{
			b1 = fix_g;
			a1 = gh_fit[i];
		}

		try
		{
			fourier_hist_NU(mg, mn, mu, data_num, a1, b1, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
		if (chi_check)
		{	// for checking
			chi_check[i] = chi_right;
			chi_check[chi_fit_num + i] = gh_fit[i];
		}
	}
	
	//st4 = clock();
	fit_shear(gh_fit, chisq_fit, chi_fit_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff, -1);

	//st5 = clock();
	//std::cout << gh << " " << gh_sig << std::endl;
	//std::cout <<"Time: "<< (st2 - st1) / CLOCKS_PER_SEC << " " << (st3 - st2) / CLOCKS_PER_SEC << " " << (st4 - st3) / CLOCKS_PER_SEC << " " << (st5 - st4) / CLOCKS_PER_SEC << std::endl;
	delete[] gh_fit;
	delete[] chisq_fit;
	delete[] num_in_bin;
}

void find_shear_NU(const float *mg, const float *mn, const float *mu, const int data_num, const int bin_num, const float *bins, int NU_label, 
					float fix_g, float &gh, float &gh_sig, float &chisq_min_fit, float *chi_check, 
					const int chi_fit_num, const float ini_left, const float ini_right, const float chi_gap)
{
	int i, j, k;
	int max_iters = 12;
	int temp_num;

	int *num_in_bin = new int[bin_num];
	float *gh_fit = new float[chi_fit_num];
	float *chisq_fit = new float[chi_fit_num];
	float chisq_fit_coeff[3];
	
	// record the each g_left, chisq_left, g_right, chisq_right
	int record_col = 4;
	int left_tag=-1, right_tag=-1;
	float fit_max_chisq, new_end;
	

	int same = 0, iters = 0, change = 1;
	float left = ini_left, right = ini_right, step;
	float chi_left, chi_right, chi_mid;
	float gh_left, gh_right, gh_mid;
	float a1, a2,a3, b1, b2,b3;
	//double st1, st2, st3, st4, st5, st6;
	//st1 = clock();
	//st2 = clock();
	while (change == 1)
	{		
		change = 0;
		gh_mid = (left + right) *0.5;
		gh_left = left;
		gh_right = right;
		if(NU_label == 1)
		{
			a1 = fix_g;	a2 = fix_g;	a3 = fix_g;
			b1 = gh_left;  b2 = gh_mid;  b3 = gh_right;
		}
		else
		{
			b1 = fix_g;	b2 = fix_g;	b3 = fix_g;
			a1 = gh_left; a2 = gh_mid; a3 = gh_right;
		}
		try
		{	
			fourier_hist_NU(mg, mn, mu, data_num, a1, b1, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_left);

			fourier_hist_NU(mg, mn, mu, data_num, a2, b2, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_mid);

			fourier_hist_NU(mg, mn, mu, data_num, a3, b3, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		//std::cout << left << " "<< gh_left<<" "<< gh_mid<<" "<< gh_right <<" "<< right << std::endl;

		if (chi_left > chi_mid + chi_gap)
		{
			left = (gh_mid + gh_left) *0.5;
			change = 1;
		}
		if (chi_right > chi_mid + chi_gap)
		{
			right = (gh_mid + gh_right)*0.5;
			change = 1;
		}

		iters += 1;
		if (iters > max_iters)
		{
			break;
		}
	}
	
	
	temp_num = 7;
	step = (right - left) / (temp_num - 1);
	for (i = 0; i < temp_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < temp_num; i++)
	{	
		if(NU_label == 1)
		{
			a1 = fix_g;
			b1 = gh_fit[i];
		}
		else
		{
			b1 = fix_g;
			a1 = gh_fit[i];
		}

		try
		{
			fourier_hist_NU(mg, mn, mu, data_num, a1, b1, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
	}
	//st4 = clock();
	fit_shear(gh_fit, chisq_fit, temp_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff, -1);

	// to get a more symmetrical interval for fitting
	new_end = std::max(gh-left, right-gh);
	
	left = gh - new_end;
	step = 2*new_end / (chi_fit_num - 1);

	for (i = 0; i < chi_fit_num; i++)
	{
		gh_fit[i] = left + step * i;
	}
	for (i = 0; i < chi_fit_num; i++)
	{	
		if(NU_label == 1)
		{
			a1 = fix_g;
			b1 = gh_fit[i];
		}
		else
		{
			b1 = fix_g;
			a1 = gh_fit[i];
		}

		try
		{
			fourier_hist_NU(mg, mn, mu, data_num, a1, b1, bins, num_in_bin, bin_num);
			cal_chisq_1d(num_in_bin, bin_num, chi_right);
		}
		catch(const char *msg)
		{
			throw msg;
		}
		chisq_fit[i] = chi_right;
		if (chi_check)
		{	// for checking
			chi_check[i] = chi_right;
			chi_check[chi_fit_num + i] = gh_fit[i];
		}
	}
	
	//st4 = clock();
	fit_shear(gh_fit, chisq_fit, chi_fit_num, gh, gh_sig, chisq_min_fit, chisq_fit_coeff, -1);

	//st5 = clock();
	//std::cout << gh << " " << gh_sig << std::endl;
	//std::cout <<"Time: "<< (st2 - st1) / CLOCKS_PER_SEC << " " << (st3 - st2) / CLOCKS_PER_SEC << " " << (st4 - st3) / CLOCKS_PER_SEC << " " << (st5 - st4) / CLOCKS_PER_SEC << std::endl;
	delete[] gh_fit;
	delete[] chisq_fit;
	delete[] num_in_bin;
}

void find_shear_iter(const double *mg, const double *mn, const double *mu, const int data_num, const int bin_num, const double *bins, int g_label, const int iters, double *result, 
						const int chi_fit_num, const double ini_left, const double ini_right, const double chi_gap)
{
	int result_col=6, i,j,k;
	double chimin;
	double gN, gN_sig, gU, gU_sig;
	double *chi_check = new double[chi_fit_num*2];
	find_shear(mg, mn, mu, data_num, bin_num, bins, g_label, gN, gN_sig, chimin, chi_check, chi_fit_num, ini_left, ini_right, chi_gap);
	
	// gN, gN_sig, gU, gU_sig, chimin_gN, chimin_gU
	result[0] = gN;
	result[1] = gN_sig;
	result[2] = gN;
	result[3] = gN_sig;
	result[4] = chimin;
	result[5] = chimin;

	for(i=1; i<iters; i++)
	{	
		find_shear_NU(mg, mn, mu, data_num, bin_num, bins, 1, gN, gU, gU_sig, chimin, chi_check, chi_fit_num, ini_left, ini_right, chi_gap);
		result[i*result_col + 2] = gU;
		result[i*result_col + 3] = gU_sig;
		result[i*result_col + 5] = chimin;

		find_shear_NU(mg, mn, mu, data_num, bin_num, bins, 2, gU, gN, gN_sig, chimin, chi_check, chi_fit_num, ini_left, ini_right, chi_gap);
		result[i*result_col] = gN;
		result[i*result_col + 1] = gN_sig;
		result[i*result_col + 4] = chimin;
	}

	delete[] chi_check;
}

void find_shear_iter(const float *mg, const float *mn, const float *mu, const int data_num, const int bin_num, const float *bins, int g_label, const int iters, float *result, 
						const int chi_fit_num, const float ini_left, const float ini_right, const float chi_gap)
{
	int result_col=6, i,j,k;
	float chimin;
	float gN, gN_sig, gU, gU_sig;
	float *chi_check = new float[chi_fit_num*2];
	find_shear(mg, mn, mu, data_num, bin_num, bins, g_label, gN, gN_sig, chimin, chi_check, chi_fit_num, ini_left, ini_right, chi_gap);
	// gN, gN_sig, gU, gU_sig, chimin_gN, chimin_gU
	result[0] = gN;
	result[1] = gN_sig;
	result[2] = gN;
	result[3] = gN_sig;
	result[4] = chimin;
	result[5] = chimin;

	for(i=1; i<iters; i++)
	{
		find_shear_NU(mg, mn, mu, data_num, bin_num, bins, 1, gN, gU, gU_sig, chimin, chi_check, chi_fit_num, ini_left, ini_right, chi_gap);
		result[i*result_col + 2] = gU;
		result[i*result_col + 3] = gU_sig;
		result[i*result_col + 5] = chimin;

		find_shear_NU(mg, mn, mu, data_num, bin_num, bins, 2, gU, gN, gN_sig, chimin, chi_check, chi_fit_num, ini_left, ini_right, chi_gap);
		result[i*result_col] = gN;
		result[i*result_col + 1] = gN_sig;
		result[i*result_col + 4] = chimin;
	}

	delete[] chi_check;
}


void fit_shear(const double *shear, const double *chisq, const int num, double &gh, double &gh_sig, double &chisq_min_fit, double *chisq_fit_coeff, const int extend_chi, const double d_chi)
{
	// fit a 2nd order 1-D curve for estimate the shear.
	// y = ax^2 + bx + c

	int i, count = 0;
	double min_chi = 10000;
	double coeff[3];
	if (d_chi > 0)
	{
		int *mask = new int[num] {};
		// find the minimum
		for (i = 0; i < num; i++)
		{
			if (chisq[i] < min_chi and chisq[i] >=0)
			{
				min_chi = chisq[i];
			}
		}
		// find the width for fitting
		for (i = 0; i < num; i++)
		{
			if (chisq[i] >= 0 and chisq[i] < min_chi + d_chi)
			{
				count++;
				mask[i] = 1;
			}
		}
		//std::cout << min_chi << std::endl;
		//show_arr(chisq, 1, num);
		// for fitting
		if (count < 5)
		{	
			if(extend_chi > 0)
			{
				while(count < 6)
				{
					for(i=0; i<num-1; i++)
					{
						if(mask[i] < 1 and mask[i+1] == 1)
						{
							mask[i] = 1;
							count ++;
						}
					}
					for(i=num-1; i>0; i--)
					{
						if(mask[i] == 0 and mask[i-1] == 1)
						{
							mask[i] = 1;
							count ++;
						}
					}
				}
			}
			else
			{
				char err_log[100];
				sprintf(err_log, "Too less points ( %d (%d) ) for fitting!!!", count, num);
				std::cout<<err_log<<std::endl;
				show_arr(chisq, 1, num);
				show_arr(mask, 1, num);
				throw err_log;
			}
		}

		double *new_chisq = new double[count] {};
		double *new_shear = new double[count] {};
		count = 0;
		for (i = 0; i < num; i++)
		{
			if (mask[i] == 1)
			{
				new_chisq[count] = chisq[i];
				new_shear[count] = shear[i];
				count++;
			}
		}

		// g`= a1 + a2*g + a3*g^2
		poly_fit_1d(new_shear, new_chisq, count, 2, coeff);

		if (coeff[2] < 0)
		{
			char err_log[35];
			sprintf(err_log, "Bad shear fitting !!!");
			show_arr(new_chisq, 1, count);
			show_arr(coeff, 1, 3);
			std::cout<<err_log<<std::endl;
			throw err_log;
		}


		delete[] mask;
		delete[] new_chisq;
		delete[] new_shear;
	}
	else
	{
		poly_fit_1d(shear, chisq, num, 2, coeff);

		if (coeff[2] < 0)
		{
			char err_log[35];
			sprintf(err_log, "Bad shear fitting !!!");
			show_arr(chisq, 1, num);
			std::cout<<err_log<<std::endl;
			throw err_log;
		}
	}

	gh = -coeff[1] / coeff[2]*0.5;
	gh_sig = sqrt(0.5 / coeff[2]);

	chisq_min_fit = coeff[0] + coeff[1]*gh + coeff[2]*gh*gh;
	chisq_fit_coeff[0] = coeff[0];
	chisq_fit_coeff[1] = coeff[1];
	chisq_fit_coeff[2] = coeff[2];

}

void fit_shear(const float *shear, const float *chisq, const int num, float &gh, float &gh_sig, float &chisq_min_fit, float* chisq_fit_coeff,const int extend_chi, const float d_chi)
{
	// fit a 2nd order 1-D curve for estimate the shear.
	// y = ax^2+bx + c

	int i, count = 0;
	float min_chi = 10000;
	float coeff[3];
	if (d_chi > 0)
	{
		int *mask = new int[num] {};
		// find the minimum
		for (i = 0; i < num; i++)
		{
			if (chisq[i] < min_chi and chisq[i] >=0)
			{
				min_chi = chisq[i];
			}
		}
		// find the width for fitting
		for (i = 0; i < num; i++)
		{
			if (chisq[i] >= 0 and chisq[i] < min_chi + d_chi)
			{
				count++;
				mask[i] = 1;
			}
		}
		//std::cout << min_chi << std::endl;
		//show_arr(chisq, 1, num);
		// for fitting
		if (count < 5)
		{
			if(extend_chi > 0)
			{
				while(count < 6)
				{
					for(i=0; i<num; i++)
					{
						if(mask[i] < 1 and mask[i+1] == 1)
						{
							mask[i] = 1;
							count ++;
						}
					}
					for(i=num-1; i>0; i--)
					{
						if(mask[i] == 0 and mask[i-1] == 1)
						{
							mask[i] = 1;
							count ++;
						}
					}
				}
			}
			else
			{
				char err_log[100];
				sprintf(err_log, "Too less points ( %d (%d) ) for fitting!!!", count, num);
				std::cout<<err_log<<std::endl;
				show_arr(chisq, 1, num);
				show_arr(mask, 1, num);
				throw err_log;
			}
		}
		float *new_chisq = new float[count] {};
		float *new_shear = new float[count] {};
		count = 0;
		for (i = 0; i < num; i++)
		{
			if (mask[i] == 1)
			{
				new_chisq[count] = chisq[i];
				new_shear[count] = shear[i];
				count++;
			}
		}

		// g`= a1 + a2*g + a3*g^2
		poly_fit_1d(new_shear, new_chisq, count, 2, coeff);

		delete[] mask;
		delete[] new_chisq;
		delete[] new_shear;
	}
	else
	{
		poly_fit_1d(shear, chisq, num, 2, coeff);
		// show_arr(shear,1, num);
		// show_arr(chisq,1, num);
	}
	if (coeff[2] < 0)
	{
		char err_log[35];
		sprintf(err_log, "Bad shear fitting !!!");
		throw err_log;
	}
	gh = -coeff[1] / coeff[2]*0.5;
	gh_sig = sqrt(0.5 / coeff[2]);

	chisq_min_fit = coeff[0] + coeff[1]*gh + coeff[2]*gh*gh;
	chisq_fit_coeff[0] = coeff[0];
	chisq_fit_coeff[1] = coeff[1];
	chisq_fit_coeff[2] = coeff[2];
}



void fourier_hist_NU(const double *mg, const double *mn,const double *mu, const int data_row, const double gN, const double gU, const double *bins, int *num_in_bin, const int bin_num)
{
	int i,j, k;
	double temp;

	initialize_arr(num_in_bin, bin_num, 0);


	for(i=0;i<data_row;i++)
	{	
		temp = mg[i] - gN*mn[i] - gU*mu[i];
		
		for(j=0;j<bin_num;j++)
		{
			if(temp >= bins[j] and temp < bins[j+1])
			{
				num_in_bin[j] += 1;
				continue;
			}
		}
	}
	
}

void fourier_hist_NU(const float *mg, const float *mn,const float *mu, const int data_row, const float gN, const float gU, const float *bins, int *num_in_bin, const int bin_num)
{
	int i,j, k;
	float temp;

	initialize_arr(num_in_bin, bin_num, 0);


	for(i=0;i<data_row;i++)
	{	
		temp = mg[i] - gN*mn[i] - gU*mu[i];;
		
		for(j=0;j<bin_num;j++)
		{
			if(temp >= bins[j] and temp < bins[j+1])
			{
				num_in_bin[j] += 1;
				continue;
			}
		}
	}	
}



void fourier_hist(const double *mg, const double *mn,const double *mu, const int data_row, const double gh, const int g_label, const double *bins, int *num_in_bin, const int bin_num)
{
	int i,j, k;
	double temp;

	initialize_arr(num_in_bin, bin_num, 0);

	if(g_label==1)
	{
		for(i=0;i<data_row;i++)
		{	
			
			temp = mg[i] - gh*(mn[i] + mu[i]);
			
			for(j=0;j<bin_num;j++)
			{
				if(temp >= bins[j] and temp < bins[j+1])
				{
					num_in_bin[j] += 1;
					continue;
				}
			}
		}
	}
	else
	{
		for(i=0;i<data_row;i++)
		{	
			temp = mg[i] - gh*(mn[i] - mu[i]);
			
			for(j=0;j<bin_num;j++)
			{
				if(temp >= bins[j] and temp < bins[j+1])
				{
					num_in_bin[j] += 1;
					continue;
				}
			}
		}
	}
}

void fourier_hist(const float *mg, const float *mn,const float *mu, const int data_row, const float gh, const int g_label, const float *bins, int *num_in_bin, const int bin_num)
{
	int i,j, k;
	float temp;

	initialize_arr(num_in_bin, bin_num, 0);

	if(g_label==1)
	{
		for(i=0;i<data_row;i++)
		{	
			
			temp = mg[i] - gh*(mn[i] + mu[i]);
			
			for(j=0;j<bin_num;j++)
			{
				if(temp >= bins[j] and temp < bins[j+1])
				{	
					// if(i <10)
					// {std::cout<<j<<" "<<bins[j]<<" "<<temp<<" "<<mg[i]<<" "<<gh<<" "<<mn[i]<<" "<<mu[i]<<" "<<bins[j+1]<<std::endl;}
					num_in_bin[j] += 1;
					continue;
				}
			}
		}
	}
	else
	{
		for(i=0;i<data_row;i++)
		{	
			temp = mg[i] - gh*(mn[i] - mu[i]);
			
			for(j=0;j<bin_num;j++)
			{
				if(temp >= bins[j] and temp < bins[j+1])
				{
					num_in_bin[j] += 1;
					continue;
				}
			}
		}
	}
}



void fourier_hist(const double *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col, const int mu_col, 
					const double gh, const int g_label, const double *bins, int *num_in_bin, const int bin_num)
{
	int i,j, k;
	double temp;

	initialize_arr(num_in_bin, bin_num, 0);

	if(g_label==1)
	{
		for(i=0;i<data_row;i++)
		{	
			k = i*data_col;
			temp = data_arr[k+mg_col] - gh*(data_arr[k+mn_col] + data_arr[k+mu_col]);
			
			for(j=0;j<bin_num;j++)
			{
				if(temp >= bins[j] and temp < bins[j+1])
				{
					num_in_bin[j] += 1;
					continue;
				}
			}
		}
	}
	else
	{
		for(i=0;i<data_row;i++)
		{	
			k = i*data_col;
			temp = data_arr[k+mg_col] - gh*(data_arr[k+mn_col] - data_arr[k+mu_col]);
			
			for(j=0;j<bin_num;j++)
			{
				if(temp >= bins[j] and temp < bins[j+1])
				{
					num_in_bin[j] += 1;
					continue;
				}
			}
		}
	}
}

void fourier_hist(const float *data_arr, const int data_row, const int data_col, const int mg_col, const int mn_col, const int mu_col,
					const float gh, const int g_label, const float *bins, int *num_in_bin, const int bin_num)
{
	int i,j, k;
	float temp;

	initialize_arr(num_in_bin, bin_num, 0);

	if(g_label==1)
	{
		for(i=0;i<data_row;i++)
		{	
			k = i*data_col;
			temp = data_arr[k+mg_col] - gh*(data_arr[k+mn_col] + data_arr[k+mu_col]);
			
			for(j=0;j<bin_num;j++)
			{
				if(temp >= bins[j] and temp < bins[j+1])
				{
					num_in_bin[j] += 1;
					continue;
				}
			}
		}
	}
	else
	{
		for(i=0;i<data_row;i++)
		{	
			k = i*data_col;
			temp = data_arr[k+mg_col] - gh*(data_arr[k+mn_col] - data_arr[k+mu_col]);
			
			for(j=0;j<bin_num;j++)
			{
				if(temp >= bins[j] and temp < bins[j+1])
				{
					num_in_bin[j] += 1;
					continue;
				}
			}
		}
	}

}



void estimator_rotation(const double theta,const double mg1, const double mg2, const double mn, const double mu, const double mv, double *output)
{
	double cos2theta, sin2theta, cos4theta, sin4theta;
	cos2theta = cos(2*theta);
	sin2theta = sin(2*theta);
	cos4theta = cos(4*theta);
	sin4theta = sin(4*theta);
	output[0] = mg1*cos2theta - mg2*sin2theta;
	output[1] = mg1*sin2theta + mg2*cos2theta;
	output[2] = mn;
	output[3] = mu*cos4theta - mv*sin4theta;
	output[4] = mu*sin4theta + mv*cos4theta;

}

void estimator_rotation(const float theta,const float mg1, const float mg2, const float mn, const float mu, const float mv, float *output)
{
	float cos2theta, sin2theta, cos4theta, sin4theta;
	cos2theta = cos(2*theta);
	sin2theta = sin(2*theta);
	cos4theta = cos(4*theta);
	sin4theta = sin(4*theta);
	output[0] = mg1*cos2theta - mg2*sin2theta;
	output[1] = mg1*sin2theta + mg2*cos2theta;
	output[2] = mn;
	output[3] = mu*cos4theta - mv*sin4theta;
	output[4] = mu*sin4theta + mv*cos4theta;

}

void image_rotation(const double *image_in, double *image_out, const int size)
{
	int i, j, m1,n1, m2, n2;
	int size2;
	size2 = size/2;
	
	for(i=-size2; i<size2+1; i++) // y
	{
		if(i == 0) continue;
		for(j=-size2; j<size2+1; j++) // x
		{	
			if(j == 0) continue;
			m1 = i + size2;
			if(i > 0){m1 -= 1;}
			n1 = j + size2;
			if(j > 0){n1 -= 1;}

			m2 = j + size2;
			n2 = size2 - i;
			if(j > 0){m2 -= 1;}
			if(i < 0){n2 -= 1;}

			image_out[m2*size + n2] = image_in[m1*size+n1];
		}
	}
}

void image_rotation(const float *image_in, float *image_out, const int size)
{
	int i, j, m1,n1, m2, n2;
	int size2;
	size2 = size/2;
	
	for(i=-size2; i<size2+1; i++) // y
	{	
		if(i == 0) continue;
		for(j=-size2; j<size2+1; j++) // x
		{	
			if(j == 0) continue;
			m1 = i + size2;
			if(i > 0){m1 -= 1;}
			n1 = j + size2;
			if(j > 0){n1 -= 1;}

			m2 = j + size2;
			n2 = size2 - i;
			if(j > 0){m2 -= 1;}
			if(i < 0){n2 -= 1;}

			image_out[m2*size + n2] = image_in[m1*size+n1];
		}
	}
}

/********************************************************************************************************************************************/
/* cosmology */
/********************************************************************************************************************************************/

void log_bin(const double start, const double end, const int num, double * bins)
{
	double st, et, dp;
	int i;
	if (start <= 0 or end <= 0 or start >= end)
	{
		std::cout << "Wrong start- or end-point !!!" << std::endl;
		exit(0);
	}
	st = log10(start);
	et = log10(end);
	dp = (et - st) / (num - 1);

	for (i = 0; i < num; i++)
	{
		bins[i] = pow(10, st + i * dp);
	}
}

void linspace(const double start, const double end, const int num, double *bins)
{
	int i, j, k;
	double diff;
	if (end <= start)
	{
		std::cout << "The end point must be larger than the start" << std::endl;
		exit(0);
	}
	diff = (end - start) / (num - 1);
	for (i = 0; i < num; i++)
	{
		bins[i] = start + i * diff;
	}
}


/********************************************************************************************************************************************/
/* fitting */
/********************************************************************************************************************************************/
void image_convole(double *image_in, double *image_out, const int img_size, const double *kernel, const int kernel_size)
{
	int i,j,m,n,p,q;
	
	int kernel_width;
	int alpha;
	double temp;

	alpha = img_size - kernel_size/2;

	for(i=0;i<img_size;i++) // y
	{
		for(j=0;j<img_size;j++) // x
		{
			temp = 0;
			for(m=0; m<kernel_size; m++) //y kernel
			{
				p = (i + m + alpha)%img_size*img_size;

				for(n=0; n<kernel_size; n++) // x kernel
				{
					q = (j + n + alpha)%img_size;
					temp += image_in[p + q]*kernel[m*kernel_size + n];
				}
			}
			image_out[i*img_size + j] = temp;
		}
	}
}

void image_convole(float *image_in, float *image_out, const int img_size, const float *kernel, const int kernel_size)
{
	int i,j,m,n,p,q;
	
	int alpha;
	float temp;
	alpha = img_size-kernel_size/2;

	for(i=0;i<img_size;i++) // y
	{
		for(j=0;j<img_size;j++) // x
		{
			temp = 0;
			for(m=0; m<kernel_size;m++) //y kernel
			{
				p = (i + m + alpha)%img_size*img_size;

				for(n=0; n<kernel_size;n++) // x kernel
				{
					q = (j + n + alpha)%img_size;
					temp += image_in[p + q]*kernel[m*kernel_size+n];
				}
			}
			image_out[i*img_size + j] = temp;
		}
	}
}

void smooth(double *image,  const double *coeffs, fq_paras*paras)//be careful of the memset()
{
	/*  to fit the curve: a1 + a2*x +a3*y + a4*x^2 +a5*x*y + a6*y^2  */
	int i, j, m, n, q, p, pk = 0, tag, cen, coe, jx, iy, size = paras->stamp_size;
	double fz[6]{}, z[25]{}, fit_para_6, max = 0., thresh;
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

void smooth(double *image, const double* psf_pow, const double *coeffs, fq_paras*paras)//be careful of the memset()
{
	/*  to fit the curve: a1 + a2*x +a3*y + a4*x^2 +a5*x*y + a6*y^2  */
	int i, j, m, n, q, p, pk = 0, tag, cen, coe, jx, iy, size = paras->stamp_size;
	double fz[6]{}, z[25]{}, fit_para_6, max = 0., thresh = paras->psf_pow_thresh/1.2; // the lower thres (divided by 1.2) for safety!!
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
			if (psf_pow[i*size + j] >= thresh)
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

void smooth_real(double *image, const double *coeffs, fq_paras*paras)//be careful of the memset()
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

void hyperfit_5(const double *data, double *fit_paras, fq_paras *paras)
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

void hyperfit_5(const float *data, float *fit_paras, fq_paras_float *paras)
{
	float temp = 0;

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


double fval_at_xy(const double x, const double y, const int order, const double *coeffs)
{		
	int terms = (order + 1)*(order + 2) / 2;
	int i, j, k;
	int *pow_y = new int[terms] {};
	int *pow_x = new int[terms] {};
	k = 0;
	for (i = 0; i < order + 1; i++)
	{
		for (j = 0; j < i + 1; j++)
		{
			pow_x[k] = i - j;
			pow_y[k] = j;
			k++;
		}
	}

	double fvals = 0;
	for (i = 0; i < terms; i++)
	{
		fvals += coeffs[i] * pow(x, pow_x[i])*pow(y, pow_y[i]);
	}
	
	delete[] pow_y;
	delete[] pow_x;
	return fvals;
}



void poly_fit_1d(const double *x, const double *fx, const int data_num, const int order, double *coeffs)
{
	if (order < 1)
	{
		std::cout << "Order < 1 !!!" << std::endl;
		exit(0);
	}
	int i, j, s;
	int terms = order + 1;
	double *cov_matrix = new double[terms*terms]{};
	double *f_vector = new double[terms] {};
	
	cov_matrix_1d(x, fx, data_num, order, cov_matrix, f_vector);

	gsl_matrix *cov_mat = gsl_matrix_alloc(terms, terms);
	gsl_vector * vect_b = gsl_vector_alloc(terms);
	gsl_matrix *mat_inv = gsl_matrix_alloc(terms, terms);
	gsl_permutation *permu = gsl_permutation_alloc(terms);
	gsl_vector *pamer = gsl_vector_alloc(terms);

	for (j = 0; j < terms; j++)
	{
		gsl_vector_set(vect_b, j, f_vector[j]);
	}

	for (i = 0; i < terms; i++)
	{
		for (j = 0; j < terms; j++)
		{
			gsl_matrix_set(cov_mat, i, j, cov_matrix[i*terms + j]);
		}
	}

	gsl_linalg_LU_decomp(cov_mat, permu, &s);
	//gsl_linalg_LU_invert(cov_mat, permu, mat_inv);
	gsl_linalg_LU_solve(cov_mat, permu, vect_b, pamer);

	for (i = 0; i < terms; i++)
	{
		coeffs[i] = gsl_vector_get(pamer, i);
	}

	gsl_matrix_free(cov_mat);
	gsl_vector_free(vect_b);
	gsl_matrix_free(mat_inv);
	gsl_vector_free(pamer);
	gsl_permutation_free(permu);

	delete[] cov_matrix;
	delete[] f_vector;
}

void poly_fit_1d(const float *x, const float *fx, const int data_num, const int order, float *coeffs)
{
	if (order < 1)
	{
		std::cout << "Order < 1 !!!" << std::endl;
		exit(0);
	}
	int i, j, s;
	int terms = order + 1;
	double *doub_coeffs = new double[terms];	
	double *doub_x = new double[data_num];
	double *doub_fx = new double[data_num];
	for(i=0; i<data_num; i++)
	{
		doub_x[i] = x[i];
		doub_fx[i] = fx[i];
	}
	poly_fit_1d(doub_x, doub_fx, data_num, order, doub_coeffs);
	for(i=0;i<terms;i++){coeffs[i] = doub_coeffs[i];}
	delete[] doub_coeffs;
	delete[] doub_x;
	delete[] doub_fx;

}


void poly_fit_1d(const double *x, const double *fx, const double *fx_err, const int data_num, double *coeffs, int weight)
{
	double chi, c0, c1, cov00, cov01, cov11;
	if (1 == weight)
	{
		double *wi = new double[data_num];
		for (int i = 0; i < data_num; i++)
		{
			if (0 == fx_err[i])
			{
				wi[i] = 1;
			}
			else
			{
				wi[i] = 1. / fx_err[i] / fx_err[i];
			}
		}
		gsl_fit_wlinear(x, 1, wi,1,fx,1,data_num,&c0,&c1,&cov00,&cov01,&cov11,&chi);
		delete[] wi;
	}
	else
	{
		gsl_fit_linear(x, 1, fx, 1, data_num, &c0, &c1, &cov00, &cov01, &cov11, &chi);
	}
	coeffs[0] = c0;
	coeffs[1] = sqrt(cov00);
	coeffs[2] = c1;
	coeffs[3] = sqrt(cov11);
}

void poly_fit_1d(const float *x, const float *fx, const float *fx_err, const int data_num, float *coeffs, int weight)
{	
	double *doub_x = new double[data_num];
	double *doub_fx = new double[data_num];
	double *doub_fx_err = new double[data_num];
	double *doub_coeffs = new double[4];
	
	int i;
	for(i=0; i<data_num; i++)
	{
		doub_x[i] = x[i];
		doub_fx[i] = fx[i];
		doub_fx_err[i] = fx_err[i];
	}
	poly_fit_1d(doub_x, doub_fx, doub_fx_err, data_num, doub_coeffs, weight);
	for(i=0;i<4;i++){coeffs[i] = doub_coeffs[i];}

	delete[] doub_x;
	delete[] doub_fx;
	delete[] doub_fx_err;
	delete[] doub_coeffs;
}


void poly_fit_2d(const double *x, const double *y, const double *fxy, const int data_num, const int order, double *coeffs)
{
	// order >= 1
	if (order < 1)
	{
		std::cout << "Order < 1 !!!" << std::endl;
		exit(0);
	}

	int terms = (order + 1)*(order + 2) / 2;
	int i, s, j, k;
	double temp_scale = 1;

	double *cov_matrix = new double[terms*terms]{};
	double *f_vector = new double[terms] {};
	double *check = new double[terms*terms]{};
	
	cov_matrix_2d(x, y, fxy, data_num, order, cov_matrix, f_vector);

	gsl_matrix *cov_mat = gsl_matrix_alloc( terms, terms);
	gsl_vector * vect_b = gsl_vector_alloc(terms);
	gsl_matrix *mat_inv = gsl_matrix_alloc( terms, terms);
	gsl_permutation *permu = gsl_permutation_alloc(terms);
	gsl_vector *pamer = gsl_vector_alloc(terms);	

	for (j = 0; j < terms; j++)
	{
		gsl_vector_set(vect_b, j, f_vector[j]/temp_scale);
	}

	for (i = 0; i < terms; i++)
	{
		for (j = 0; j < terms; j++)
		{
			gsl_matrix_set(cov_mat, i, j, cov_matrix[i*terms + j]/temp_scale);
		}
	}

	gsl_linalg_LU_decomp (cov_mat, permu, &s);
	gsl_linalg_LU_invert (cov_mat, permu, mat_inv);	
	gsl_linalg_LU_solve (cov_mat, permu, vect_b, pamer);

	for (i = 0; i < terms; i++)
	{
		coeffs[i] = gsl_vector_get(pamer, i);
	}

	gsl_matrix_free(cov_mat);
	gsl_vector_free(vect_b);
	gsl_matrix_free(mat_inv);
	gsl_vector_free(pamer);
	gsl_permutation_free(permu);	

	delete[] cov_matrix;
	delete[] f_vector;
	delete[] check;
}



void cov_matrix_1d(const double *x, const double *fx, const int data_num, const int order, double *cov_mat, double *f_vector)
{
	if (order < 1)
	{
		std::cout << "Order < 1 !!!" << std::endl;
		exit(0);
	}
	int i, j, k, m;
	int terms = order + 1;
	double *xn = new double[(2*order+1)*data_num];
	double elem_sum = 0;

	for (j = 0; j < data_num; j++)
	{
		xn[j] = 1;
		xn[data_num + j] = x[j];
	}
	for (i = 2; i < 2 * order + 1; i++)
	{
		for (j = 0; j < data_num; j++)
		{
			xn[i*data_num + j] = pow(x[j], i);
		}	
	}

	for (i = 0; i < terms; i++)
	{
		elem_sum = 0;
		for (m = 0; m < data_num; m++)
		{
			elem_sum += fx[m] * xn[i*data_num + m];
		}
		f_vector[i] = elem_sum;

		for (j = 0; j < terms; j++)
		{
			elem_sum = 0;
			for (k = 0; k < data_num; k++)
			{
				elem_sum += xn[(i + j)*data_num + k];
			}
			cov_mat[i*terms+ j] = elem_sum;
		}
	}

	delete[] xn;
}

void cov_matrix_1d(const float *x, const float *fx, const int data_num, const int order, float *cov_mat, float *f_vector)
{
	if (order < 1)
	{
		std::cout << "Order < 1 !!!" << std::endl;
		exit(0);
	}
	int i, j, k, m;
	int terms = order + 1;
	float *xn = new float[(2*order+1)*data_num];
	float elem_sum = 0;

	for (j = 0; j < data_num; j++)
	{
		xn[j] = 1;
		xn[data_num + j] = x[j];
	}
	for (i = 2; i < 2 * order + 1; i++)
	{
		for (j = 0; j < data_num; j++)
		{
			xn[i*data_num + j] = pow(x[j], i);
		}	
	}

	for (i = 0; i < terms; i++)
	{
		elem_sum = 0;
		for (m = 0; m < data_num; m++)
		{
			elem_sum += fx[m] * xn[i*data_num + m];
		}
		f_vector[i] = elem_sum;

		for (j = 0; j < terms; j++)
		{
			elem_sum = 0;
			for (k = 0; k < data_num; k++)
			{
				elem_sum += xn[(i + j)*data_num + k];
			}
			cov_mat[i*terms+ j] = elem_sum;
		}
	}

	delete[] xn;
}

void cov_matrix_2d(const double *x, const double *y, const double *fxy, const int data_num, const int order, double *cov_mat, double *f_vector)
{
	// order >= 1;
	if (order < 1)
	{
		std::cout << "Order < 1 !!!" << std::endl;
		exit(0);
	}

	int terms = (order + 1)*(order + 2) / 2;
	// rows * data number, "row+1" represent the power of x(y)
	int size = 2 * order * data_num;
	// the length of the square matrix of powers of x and y
	int mask_size = (2 * order + 1)*(2 * order + 1);
	// the number of possible "x^n*y^m" term
	int xys_size = (2 * order - 1)*order * data_num, xy_seq;
	int i, j, k, m, n;
	int y_odr, x_odr, xy_count = 0;
	double dn_sum = 0;

	// the powers of x and y of each turn in the polynomial
	int *pow_y = new int[terms] {};
	int *pow_x = new int[terms] {};
	double *ys = new double[size] {};
	double *xs = new double[size] {};

	// xy_pow_mask stores the location of "x^n*y^m" in the array "xys"
	// xy_pow_mask is used as 2 dimessional array , (y order, x order)
	int *xy_pow_mask = new int[mask_size] {};
	double *xys = new double[xys_size] {};

	for (i = 0; i < mask_size; i++)
	{
		xy_pow_mask[i] = -1;
	}

	k = 0;
	for (i = 0; i < order + 1; i++)
	{
		for (j = 0; j < i + 1; j++)
		{
			pow_x[k] = i - j;
			pow_y[k] = j;
			k++;
		}
	}

	// calculate each order of x and y for the covariance matrix
	// to x^(order*order), y^(order*order)
	for (i = 0; i < 2*order; i++)
	{
		
		if (0 == i) // the first order
		{
			for (k = 0; k < data_num; k++)
			{
				ys[k] = y[k];
				xs[k] = x[k];
			}
		}
		else // the higher order  
		{
			m = i*data_num;
			n = (i - 1)*data_num;
			for (k = 0; k < data_num; k++)
			{
				ys[m + k] = ys[n + k] * y[k];
				xs[m + k] = xs[n + k] * x[k];
			}
		}
	}

	// calculate the covariance matrix		
	//row
	for (i = 0; i < terms; i++)
	{
		//column
		for (j = 0; j < terms; j++)
		{
			y_odr = pow_y[j] + pow_y[i];
			x_odr = pow_x[j] + pow_x[i];

			//std::cout << x_odr << " " << y_odr << ", ";
			//if (j == terms - 1)
			//{
			//	std::cout << std::endl;
			//}

			// the terms without x^n
			if (0 == x_odr)
			{
				//the cov[0,0] = number of data points
				if (0 == y_odr)
				{
					dn_sum = data_num;
				}
				// the y^n terms
				else
				{
					sum_arr(ys, size, (y_odr - 1)*data_num, y_odr*data_num, dn_sum);
				}
			}
			// the terms with x^n
			else
			{
				// the x^n terms
				if (0 == y_odr)
				{
					sum_arr(xs, size, (x_odr - 1)*data_num, x_odr*data_num, dn_sum);
				}
				// the x^n*y^m terms
				else
				{
					// if this term has been gotten
					if (xy_pow_mask[y_odr*(2*order + 1) + x_odr] > -1)
					{
						xy_seq = xy_pow_mask[y_odr*(2*order + 1) + x_odr];
						sum_arr(xys, xys_size, xy_seq*data_num, (xy_seq + 1)*data_num, dn_sum);
					}
					// if not
					else
					{
						xy_pow_mask[y_odr*(2*order + 1) + x_odr] = xy_count;
						for (k = 0; k < data_num; k++)
						{
							xys[xy_count*data_num + k] = xs[(x_odr - 1)*data_num + k] * ys[(y_odr - 1)*data_num + k];
						}
						sum_arr(xys, xys_size, xy_count*data_num, (xy_count + 1)*data_num, dn_sum);
						xy_count++;
					}
				}
			}
			cov_mat[i*terms + j] = dn_sum;
		}
	}

	// the vector of the right side
	for (i = 0; i < terms; i++)
	{
		dn_sum = 0;
		y_odr = pow_y[i];
		x_odr = pow_x[i];

		//std::cout << x_odr << " " << y_odr << ", ";
		// the terms without x^n
		if (0 == x_odr)
		{
			if (0 == y_odr)
			{
				sum_arr(fxy, data_num, 0, data_num, dn_sum);
			}
			// the f(x,y)*y^n terms
			else
			{
				for (j = 0; j < data_num; j++)
				{
					dn_sum += fxy[j] * ys[(y_odr - 1)*data_num + j];
				}
			}
		}
		// the terms with x^n
		else
		{
			// the f(x,y)*x^n terms
			if (0 == y_odr)
			{
				for (j = 0; j < data_num; j++)
				{
					dn_sum += fxy[j] * xs[(x_odr - 1)*data_num + j];
				}
			}
			// the f(x,y)* x^n*y^m terms
			else
			{
				xy_seq = xy_pow_mask[y_odr*(2*order + 1) + x_odr];
				for (j = 0; j < data_num; j++)
				{
					dn_sum += fxy[j] * xys[xy_seq*data_num + j];
				}
			}
		}
		f_vector[i] = dn_sum;
	}

	delete[] xys;
	delete[] pow_x;
	delete[] pow_y;
	delete[] ys;
	delete[] xs;
	delete[] xy_pow_mask;
}



void sum_arr(const double *arr, const int size, const int start_t, const int stop_t, double &total)
{
	double temp = 0;
	if (stop_t > size)
	{	
		std::cout << "Cross the boundary of array!!! Stop_t: " <<stop_t<<",   Size: "<<size<< std::endl;
		exit(0);
	}
	for (int i = start_t; i < stop_t; i++)
	{
		temp += arr[i];
	}
	total = temp;
}

void sum_arr(const long *arr, const int size, const int start_t, const int stop_t, long &total)
{
	long temp = 0;
	if (stop_t > size)
	{
		std::cout << "Cross the boundary of array!!! Stop_t: " << stop_t << ",   Size: " << size << std::endl;
		exit(0);
	}
	for (int i = start_t; i < stop_t; i++)
	{
		temp += arr[i];
	}
	total = temp;
}

void sum_arr(const int *arr, const int size, const int start_t, const int stop_t, int &total)
{
	int temp = 0;
	if (stop_t > size)
	{
		std::cout << "Cross the boundary of array!!! Stop_t: " << stop_t << ",   Size: " << size << std::endl;
		exit(0);
	}
	for (int i = start_t; i < stop_t; i++)
	{
		temp += arr[i];
	}
	total = temp;
}

void arr_pow(const double *arr, double *arr_out, const int size, const int alpha, const int beta, const double power)
{	
	int i;
	if (1 == alpha)
	{
		if (0 == beta)
		{
			for ( i = 0; i < size; i++)
			{
				arr_out[i] = pow(arr[i], power);
			}
		}
		else
		{
			for ( i = 0; i < size; i++)
			{
				arr_out[i] = pow(arr[i]+beta, power);
			}
		}
	}
	else
	{
		if (0 == beta)
		{
			for (i = 0; i < size; i++)
			{
				arr_out[i] = pow(arr[i]*alpha, power);
			}
		}
		else
		{
			for (i = 0; i < size; i++)
			{
				arr_out[i] = pow(arr[i]*alpha + beta, power);
			}
		}
	}
}

void arr_rescale(double *x, const double dx, const double scale, const int num)
{
	for (int i = 0; i < num; i++)
	{
		x[i] = scale * (x[i] + dx);
	}
}

void matrix_product(const double*arr_left, const int size_1, const int size_2, const int size_3, const double *arr_right, double *result)
{
	gsl_matrix *a = gsl_matrix_alloc(size_1, size_2);
	gsl_matrix *b = gsl_matrix_alloc(size_2, size_3);
	gsl_matrix *c = gsl_matrix_alloc(size_1, size_3);

	int i, j;

	for (i = 0; i < size_1; i++)
	{
		for (j = 0; j < size_2; j++)
		{
			gsl_matrix_set(a, i, j, arr_left[i*size_2 + j]);
		}
	}
	for (i = 0; i < size_2; i++)
	{
		for (j = 0; j < size_3; j++)
		{
			gsl_matrix_set(b, i, j, arr_right[i*size_3 + j]);
		}
	}
	for (i = 0; i < size_1; i++)
	{
		for (j = 0; j < size_3; j++)
		{
			gsl_matrix_set(c, i, j, 0);
		}
	}

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0, a, b, 0.0, c);

	for (i = 0; i < size_1; i++)
	{
		for (j = 0; j < size_3; j++)
		{
			result[i*size_3+j] = gsl_matrix_get(c, i, j);

		}
	}

	gsl_matrix_free(a);
	gsl_matrix_free(b);
	gsl_matrix_free(c);
}

void matrix_inv(const double *arr, const int size, double *arr_inv)
{
	gsl_matrix *gsl_arr = gsl_matrix_alloc(size, size);
	gsl_matrix *gsl_arr_inv = gsl_matrix_alloc(size, size);
	gsl_permutation *permu = gsl_permutation_alloc(size);
	int i, j, s;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			gsl_matrix_set(gsl_arr, i, j, arr[i*size + j]);
		}
	}
	gsl_linalg_LU_decomp(gsl_arr, permu, &s);
	gsl_linalg_LU_invert(gsl_arr, permu, gsl_arr_inv);
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			arr_inv[i*size + j] = gsl_matrix_get(gsl_arr_inv, i, j);
		}
	}

	gsl_matrix_free(gsl_arr);
	gsl_matrix_free(gsl_arr_inv);
	gsl_permutation_free(permu);
}



/********************************************************************************************************************************************/
/* general methods */
/********************************************************************************************************************************************/



void find_near(const double *arr, const double tar_val, const int arr_len, int &label)
{
	double near;
	int sl, sm, sr, ds, tag;
	sl = 0;
	sr = arr_len - 1;
	sm = int ((sr - sl)*0.5);
	ds = sr - sl;
	if (ds >= 4)
	{
		while (true)
		{
			if (arr[sm] <= tar_val)
			{
				sl = sm;
				sm = int((sr + sm)*0.5);
				//std::cout << sl << " " << sm << " " << sr << std::endl;
			}
			else
			{
				sr = sm;
				sm = int((sl +sr)*0.5);
				//std::cout << sl << " " << sm << " " << sr << std::endl;
			}
			ds = sr - sl;
			if (ds <= 4)
			{
				break;
			}
		}
	}

	near = fabs(arr[sr] - tar_val);
	tag = sr;
	for (int i = sl; i < sr + 1; i++)
	{
		if (fabs(arr[i] - tar_val) < near)
		{
			near = fabs(arr[i] - tar_val);
			tag = i;
		}
	}
	label = tag;
}

void check_buffer(double *target_arr, double *buffer, const int start_t, const int buffer_size, int & count, int count_line)
{
	if (count > count_line)
	{
		int i;
		for (i = 0; i < buffer_size; i++)
		{
			target_arr[start_t + i] += buffer[i];
			buffer[i] = 0;
		}
		count = 0; //re-count
	}
}

void task_alloc(const int *label_list, const int total_task_num, const int portion, const int portion_label, int *allocated_list)
{
	int i, j, k, m, n;
	m = total_task_num / portion;
	n = total_task_num % portion;
	for (i = 0; i < total_task_num; i++)
	{
		allocated_list[i] = -1;
	}
	
	for (i = 0; i < m; i++)
	{
		allocated_list[i] = label_list[m*portion_label + i];
	}
	
	if (portion_label< n)
	{
		allocated_list[m] = label_list[m*portion + portion_label];
	}
}

void task_alloc(const int total_task_num, const int portion, const int my_part_id, int &my_start, int &my_end, int *task_count)
{	
	int i,j;
	int sub_num;
	int st, ed;

	sub_num = total_task_num / portion;
	j = total_task_num%portion;
	for(i=0;i<portion;i++)
	{
		task_count[i] = sub_num;
		if(i<j){task_count[i] +=1;}
	}
	st = 0;
	for(i=0;i<my_part_id;i++)
	{
		st += task_count[i];
	}
	my_start = st;
	my_end = st + task_count[my_part_id];
}

void task_alloc(const int total_task_num, const int division_num, const int my_part_id, int &my_st_id, int &my_ed_id, int *task_count, int *entry_for_gather)
{
    int i,j,m,n;
    m = total_task_num/division_num;
    n = total_task_num%division_num;

    for(i=0;i<division_num;i++)
    {
        task_count[i] = m;
        if(i<n){task_count[i] +=1;}
    }
    m=0;
    n=0;
    for(i=0;i<division_num;i++)
    {   
        m=0;
        for(j=0;j<i;j++)
        {
            m+=task_count[j];
        }
        entry_for_gather[i]=m;
    }
    n = m+task_count[my_part_id];
    my_st_id = entry_for_gather[my_part_id];
    my_ed_id = my_st_id + task_count[my_part_id];
}

void show_arr(const double*arr, const int rows, const int cols)
{
	int i, j;
	if (rows >= 1)
	{
		for (i = 0; i < rows; i++)
		{
			for (j = 0; j < cols; j++)
			{
				std::cout << arr[i*cols + j] << "  ";
			}
			std::cout << std::endl;
		}
	}
	else
	{
		std::cout << "rows must >= 1 !!!" << std::endl;
		exit(0);
	}
}
void show_arr(const long*arr, const int rows, const int cols)
{
	int i, j;
	if (rows >= 1)
	{
		for (i = 0; i < rows; i++)
		{
			for (j = 0; j < cols; j++)
			{
				std::cout << arr[i*cols + j] << "  ";
			}
			std::cout << std::endl;
		}
	}
	else
	{
		std::cout << "rows must >= 1 !!!" << std::endl;
		exit(0);
	}
}

void show_arr(const float*arr, const int rows, const int cols)
{
	int i, j;
	if (rows >= 1)
	{
		for (i = 0; i < rows; i++)
		{
			for (j = 0; j < cols; j++)
			{
				std::cout << arr[i*cols + j] << "  ";
			}
			std::cout << std::endl;
		}
	}
	else
	{
		std::cout << "rows must >= 1 !!!" << std::endl;
		exit(0);
	}
}

void show_arr(const int*arr, const int rows, const int cols)
{
	int i, j;
	if (rows >= 1)
	{
		for (i = 0; i < rows; i++)
		{
			for (j = 0; j < cols; j++)
			{
				std::cout << arr[i*cols + j] << "  ";
			}
			std::cout << std::endl;
		}
	}
	else
	{
		std::cout << "rows must >= 1 !!!" << std::endl;
		exit(0);
	}
}

void initialize_para(fq_paras *paras)
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

void initialize_para(fq_paras_float *paras)
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

void set_bin(double *data, const int data_num, double * bins, const int bin_num, const double max_scale)
{
	int i;
	int mid = bin_num / 2, step, num;
	
	for (i = 0; i < data_num; i++)
	{
		data[i] = fabs(data[i]);
	}

	step = data_num / bin_num * 2;
	sort_arr(data, data_num, 1);
	
	// make the boundary big enough to enclose all the data
	bins[0] = -data[data_num-1] * max_scale;
	bins[bin_num] = data[data_num-1] * max_scale;
	bins[mid] = 0;
	for (i = 1; i < bin_num / 2; i++)
	{
		bins[mid + i] = data[step*i];
		bins[mid - i] = -data[step*i];
	}

}
void set_bin(float *data, const int data_num, float * bins, const int bin_num, const float max_scale)
{
	int i;
	int mid = bin_num / 2, step, num;
	
	for (i = 0; i < data_num; i++)
	{
		data[i] = fabs(data[i]);
	}

	step = data_num / bin_num * 2;
	sort_arr(data, data_num, 1);
	
	// make the boundary big enough to enclose all the data
	bins[0] = -data[data_num-1] * max_scale;
	bins[bin_num] = data[data_num-1] * max_scale;
	bins[mid] = 0;
	for (i = 1; i < bin_num / 2; i++)
	{
		bins[mid + i] = data[step*i];
		bins[mid - i] = -data[step*i];
	}
}

void set_bin(const double *data, const int data_num, const int bin_num, double * bins, const double max_scale, int choice)
{
	int i;
	int mid = bin_num / 2, step, num;
	double *data_cp, data_max = data[0], data_min = data[0], bound;
	for (i = 0; i < data_num; i++)
	{
		if (data[i] > data_max)
		{
			data_max = data[i];
		}
		if (data[i] < data_min)
		{
			data_min = data[i];
		}
	}
	data_min = fabs(data_min);
	bound = std::max(data_min, data_max);
	// std::cout<<data_min<<" "<<data_max<<std::endl;

	if (0 >= choice)
	{
		num = data_num;
		data_cp = new double[num];		
		for (i = 0; i < num; i++)
		{
			data_cp[i] = fabs(data[i]);
		}
	}
	else
	{	
		// choice the data "randomly" to save time
		num = choice;
		int ch_step = data_num / choice;		
		data_cp = new double[num];
		for (i = 0; i < choice; i++)
		{
			data_cp[i] = fabs(data[i*ch_step]);
		}
	}

	step = num / bin_num * 2;
	sort_arr(data_cp, num, 1);
	
	// make the boundary big enough to enclose all the data
	bins[0] = -bound * max_scale;
	bins[bin_num] = bound * max_scale;
	bins[mid] = 0;
	for (i = 1; i < bin_num / 2; i++)
	{
		bins[mid + i] = data_cp[step*i];
		bins[mid - i] = -data_cp[step*i];
	}
	delete[] data_cp;
}

void set_bin(const float *data, const int data_num, const int bin_num, float * bins, const float max_scale, int choice)
{
	int i;
	int mid = bin_num / 2, step, num;
	float *data_cp, data_max = data[0], data_min = data[0], bound;
	for (i = 0; i < data_num; i++)
	{
		if (data[i] > data_max)
		{
			data_max = data[i];
		}
		if (data[i] < data_min)
		{
			data_min = data[i];
		}
	}
	data_min = fabs(data_min);
	bound = std::max(data_min, data_max);

	if (choice < 0)
	{
		std::cout << "choice must be non-negative!!!" << std::endl;
		exit(0);
	}
	else if (0 == choice)
	{
		num = data_num;
		data_cp = new float[num];
		for (i = 0; i < num; i++)
		{
			data_cp[i] = fabs(data[i]);
		}
	}
	else
	{
		// choice the data "randomly" to save time
		num = choice;
		int ch_step = data_num / choice;
		data_cp = new float[num];
		for (i = 0; i < choice; i++)
		{
			data_cp[i] = fabs(data[i*ch_step]);
		}
	}

	step = num / bin_num * 2;
	sort_arr(data_cp, num, 1);

	// make the boundary big enough to enclose all the data
	bins[0] = -bound * max_scale;
	bins[bin_num] = bound * max_scale;
	bins[mid] = 0;
	for (i = 1; i < bin_num / 2; i++)
	{
		bins[mid + i] = data_cp[step*i];
		bins[mid - i] = -data_cp[step*i];
	}
	delete[] data_cp;
}

void set_bin(const int *data, const int data_num, int * bins, const int bin_num, const int max_scale, int choice)
{
	int i;
	int mid = bin_num / 2, step, num;
	int *data_cp, data_max = data[0], data_min = data[0], bound;
	for (i = 0; i < data_num; i++)
	{
		if (data[i] > data_max)
		{
			data_max = data[i];
		}
		if (data[i] < data_min)
		{
			data_min = data[i];
		}
	}
	data_min = fabs(data_min);
	bound = std::max(data_min, data_max);

	if (choice < 0)
	{
		std::cout << "choice must be non-negative!!!" << std::endl;
		exit(0);
	}
	else if (0 == choice)
	{
		num = data_num;
		data_cp = new int[num];
		for (i = 0; i < num; i++)
		{
			data_cp[i] = fabs(data[i]);
		}
	}
	else
	{
		// choice the data "randomly" to save time
		num = choice;
		int ch_step = data_num / choice;
		data_cp = new int[num];
		for (i = 0; i < choice; i++)
		{
			data_cp[i] = fabs(data[i*ch_step]);
		}
	}

	step = num / bin_num * 2;
	sort_arr(data_cp, num, 1);

	// make the boundary big enough to enclose all the data
	bins[0] = -bound * max_scale;
	bins[bin_num] = bound * max_scale;
	bins[mid] = 0;
	for (i = 1; i < bin_num / 2; i++)
	{
		bins[mid + i] = data_cp[step*i];
		bins[mid - i] = -data_cp[step*i];
	}
	delete[] data_cp;
}



void histogram(const double *data, const double *bins, int *num_in_bin, const int data_num, const int bin_num)
{
	// initialize
	int i, j;
	initialize_arr(num_in_bin, bin_num, 0);
	for ( i = 0; i < data_num; i++)
	{
		for ( j = 0; j < bin_num; j++)
		{
			if (data[i] < bins[j + 1] && data[i] >= bins[j])
			{
				num_in_bin[j] += 1;
				break;
			}
		}
	}
}

void histogram(const double *data, const double *bins, long *num_in_bin, const int data_num, const int bin_num)
{
	// initialize
	int i, j;
	initialize_arr(num_in_bin, bin_num, 0);
	for (i = 0; i < data_num; i++)
	{
		for (j = 0; j < bin_num; j++)
		{
			if (data[i] < bins[j + 1] && data[i] >= bins[j])
			{
				num_in_bin[j] += 1;
				break;
			}
		}
	}
}

void histogram(const float *data, const float *bins, int *num_in_bin, const int data_num, const int bin_num)
{	
	// initialize
	int i, j;
	initialize_arr(num_in_bin, bin_num, 0);
	for ( i = 0; i < data_num; i++)
	{
		for ( j = 0; j < bin_num; j++)
		{
			if (data[i] < bins[j + 1] && data[i] >= bins[j])
			{
				num_in_bin[j] += 1;
				break;
			}
		}
	}
}

void histogram(const float *data, const float *bins, long *num_in_bin, const int data_num, const int bin_num)
{
	// initialize
	int i, j;
	initialize_arr(num_in_bin, bin_num, 0);
	for (i = 0; i < data_num; i++)
	{
		for (j = 0; j < bin_num; j++)
		{
			if (data[i] < bins[j + 1] && data[i] >= bins[j])
			{
				num_in_bin[j] += 1;
				break;
			}
		}
	}
}

void histogram(const int *data, const  int *bins, int *num_in_bin, const  int data_num, const  int bin_num)
{
	// initialize
	int i, j;
	initialize_arr(num_in_bin, bin_num,0);
	for ( i = 0; i < data_num; i++)
	{
		for ( j = 0; j < bin_num; j++)
		{
			if (data[i] < bins[j + 1] && data[i] >= bins[j])
			{
				num_in_bin[j] += 1;
				break;
			}
		}
	}
}

void histogram(const int *data, const  int *bins, long *num_in_bin, const  int data_num, const  int bin_num)
{
	// initialize
	int i, j;
	initialize_arr(num_in_bin, bin_num, 0);
	for (i = 0; i < data_num; i++)
	{
		for (j = 0; j < bin_num; j++)
		{
			if (data[i] < bins[j + 1] && data[i] >= bins[j])
			{
				num_in_bin[j] += 1;
				break;
			}
		}
	}
}



void histogram2d(const double *data_y, const double*data_x, const double *bin_y, const double *bin_x, int *num, const int data_num, const int ybin_num, const int xbin_num)
{
	int i, j, k;
	// initialize
	initialize_arr(num, ybin_num*xbin_num, 0);
	// bin_num = len(bins) - 1
	for ( k = 0; k < data_num; k++)
	{
		// loop y-bins
		for ( i = 0; i < ybin_num; i++)
		{
			if (data_y[k] < bin_y[i + 1] && data_y[k] >= bin_y[i])
			{
				// loop x-bins
				for ( j = 0; j < xbin_num; j++)
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
	int i, j, k;
	// initialize
	initialize_arr(num, ybin_num*xbin_num, 0);
	// bin_num = len(bins) - 1
	for ( k = 0; k < data_num; k++)
	{
		// loop y-bins
		for ( i = 0; i < ybin_num; i++)
		{
			if (data_y[k] < bin_y[i + 1] && data_y[k] >= bin_y[i])
			{
				// loop x-bins
				for ( j = 0; j < xbin_num; j++)
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
	// initialize
	int i, j, k;
	initialize_arr(num, ybin_num*xbin_num, 0);
	// bin_num = len(bins) - 1
	for ( k = 0; k < data_num; k++)
	{
		// loop y-bins
		for ( i = 0; i < ybin_num; i++)
		{
			if (data_y[k] < bin_y[i + 1] && data_y[k] >= bin_y[i])
			{
				// loop x-bins
				for ( j = 0; j < xbin_num; j++)
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

void histogram_s(const double data, const double *data_bin, const int bin_num, int &bin_label)
{
	int i;
	// loop y-bins

	for (i = 0; i < bin_num; i++)
	{
		if (data >= data_bin[i] and data < data_bin[i + 1])
		{
			bin_label = i;
			break;
		}
	}

}

void histogram2d_s(const double data_y, const double data_x, const double *bin_y, const double *bin_x, const int ybin_num, const  int xbin_num, int &bin_label)
{
	int i, j;
	// loop y-bins

	for ( i = 0; i < ybin_num; i++)
	{
		if (data_y >= bin_y[i] and data_y < bin_y[i + 1])
		{
			// loop x-bins
			for ( j = 0; j < xbin_num; j++)
			{
				if (data_x >= bin_x[j] and data_x < bin_x[j + 1])
				{
					bin_label = i * xbin_num + j;
					break;
				}
			}
			break;
		}
	}

}



void sort_arr(double* arr, int size, int order=1)
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

void sort_arr(float *arr, int size, int order=1)
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

void sort_arr(int *arr, int size, int order=1)
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


void get_time(char *str_time, int length)
{	
	/* the length of str_time should be larger than 40 */
	time_t timer;
	time(&timer);
	strftime(str_time, length, "%Y-%m-%d %H:%M:%S", localtime(&timer));
}


/******************************************* random **********************************************************/

void rand_shuffle(double *seq, int length, const gsl_rng *gsl_rand_rng)
{
	int rand_i;
	double rand_n;
	for (int i = 0; i < length-1; i++)
	{	
		rand_uniform(0, length, rand_n, gsl_rand_rng);
		rand_i = int(rand_n);
		std::swap(seq[i], seq[rand_i]);
	}
}

void rand_shuffle(float *seq, int length, const gsl_rng *gsl_rand_rng)
{
	int rand_i;
	double rand_n;
	for (int i = 0; i < length-1; i++)
	{
		rand_uniform(0, length, rand_n, gsl_rand_rng);
		rand_i = int(rand_n);
		std::swap(seq[i], seq[rand_i]);
	}
}

void rand_shuffle(int *seq, int length, const gsl_rng *gsl_rand_rng)
{
	int rand_i;
	double rand_n;
	for (int i = 0; i < length - 1; i++)
	{
		rand_uniform(0, length, rand_n, gsl_rand_rng);
		rand_i = int(rand_n);
		std::swap(seq[i], seq[rand_i]);
	}
}

/************************* GSL ***********************************/
// if use GSL, you must call gsl_setup() firtstly
// then initialize the gsl_rng, gsl_initialize(gsl_rng)
// finally, you must free the gsl_rng by calling gsl_free(gsl_rng)
void gsl_setup()
{	
	if(GSL_SETUP_LABEL == 0)
	{
		gsl_rng_env_setup();
		GSL_SETUP_LABEL = 1;
	}
}

void gsl_initialize(int seed, const int rng_label)
{
	gsl_setup();

	if(rng_label==0)
	{
		T0 = gsl_rng_ranlxs0;
		rng0 = gsl_rng_alloc(T0);
		gsl_rng_set(rng0, seed);
	}
	if(rng_label==1)
	{
		T1 = gsl_rng_ranlxs1;
		rng1 = gsl_rng_alloc(T1);
		gsl_rng_set(rng1, seed);
	}
	if(rng_label==2)
	{
		T2 = gsl_rng_ranlxs2;
		rng2 = gsl_rng_alloc(T2);
		gsl_rng_set(rng2, seed);
	}
	if(rng_label==3)
	{
		T3 = gsl_rng_mt19937;
		rng3 = gsl_rng_alloc(T3);
		gsl_rng_set(rng3, seed);
	}
}

void gsl_free(const int rng_label)
{	
	if(rng_label==0)
	{
		gsl_rng_free(rng0);
	}
	if(rng_label==1)
	{
		gsl_rng_free(rng1);
	}
	if(rng_label==2)
	{
		gsl_rng_free(rng2);
	}
	if(rng_label==3)
	{
		gsl_rng_free(rng3);
	}	
}

void noise2pow(double *pow_img, const int size, const double sigma, gsl_rng *gsl_rand_rng)
{
    int i,j, kx1, ky1, kx2, ky2;
    int cent;
    double noise;
    cent = size*0.5;

    for(i=1;i<cent;i++) // y
    {
        ky1 = (cent + i)*size;
        ky2 = (cent - i)*size;

        // ky = size/2
        rand_gauss(sigma, 0, noise, gsl_rand_rng);
        pow_img[cent*size + cent+i] += noise;
        pow_img[cent*size + cent- i] += noise;

        // kx = size/2
        rand_gauss(sigma, 0, noise, gsl_rand_rng);
        pow_img[ky1 + cent] += noise;
        pow_img[ky2 + cent] += noise;

        // ky = 0, kx = [1, size]
        rand_gauss(sigma, 0, noise, gsl_rand_rng);
        pow_img[cent + i] += noise;
        pow_img[cent - i] += noise;

        // kx=0, ky = [1, size]
        rand_gauss(sigma, 0, noise, gsl_rand_rng);
        pow_img[ky1] += noise;
        pow_img[ky2] += noise;

        for(j=1;j<cent;j++) // x
        {
            kx1 = cent + j;
            kx2 = cent - j;

            rand_gauss(sigma, 0, noise, gsl_rand_rng);
            pow_img[ky1 + kx1] += noise;
            pow_img[ky2 + kx2] += noise;
            
            rand_gauss(sigma, 0, noise, gsl_rand_rng);
            pow_img[ky1 + kx2] += noise;
            pow_img[ky2 + kx1] += noise;
        }
    }

    rand_gauss(sigma, 0, noise, gsl_rand_rng);
    pow_img[0] += noise; // [0,0]
    rand_gauss(sigma, 0, noise, gsl_rand_rng);
    pow_img[cent] += noise; // [0, size/2]
    rand_gauss(sigma, 0, noise, gsl_rand_rng);
    pow_img[cent*size] += noise; // [size/2, 0]
    rand_gauss(sigma, 0, noise, gsl_rand_rng);
    pow_img[cent*size + cent] += noise; // [size/2, size/2]
}

void addnoise(double *image, const int pixel_num, const double sigma, const gsl_rng *gsl_rand_rng)
{
	for (int i = 0; i < pixel_num; i++)
	{
		image[i] = image[i] + gsl_ran_gaussian(gsl_rand_rng, sigma);
	}
}

void addnoise(float *image, const int pixel_num, const float sigma, const gsl_rng *gsl_rand_rng)
{
	for (int i = 0; i < pixel_num; i++)
	{
		image[i] = image[i] + gsl_ran_gaussian(gsl_rand_rng, sigma);
	}
}


// void rand_multi_gauss(const double*cov, const double *mu, const int num, double *result, const gsl_rng *gsl_rand_rng)
// {	
// 	int i, j;
// 	gsl_matrix *covs = gsl_matrix_alloc(num, num);
// 	gsl_vector *mus = gsl_vector_alloc(num);
// 	gsl_vector *res = gsl_vector_alloc(num);

// 	for (i = 0; i < num; i++)
// 	{
// 		gsl_vector_set(mus, i, mu[i]);
// 		for (j = 0; j < num; j++)
// 		{
// 			gsl_matrix_set(covs, i, j, cov[i*num + j]);
// 			//std::cout << gsl_matrix_get(covs, i, j) << ", ";
// 		}
// 		//std::cout << std::endl;
// 	}

// 	gsl_linalg_cholesky_decomp(covs);
// 	gsl_ran_multivariate_gaussian(gsl_rand_rng, mus, covs, res);

// 	for (i = 0; i < num; i++)
// 	{
// 		result[i] = gsl_vector_get(res, i);
// 	}

// 	gsl_matrix_free(covs);
// 	gsl_vector_free(mus);
// 	gsl_vector_free(res);
// }

void rand_gauss(double sigma, double mean, double &rand_n, const gsl_rng *gsl_rand_rng)
{
	rand_n = gsl_ran_gaussian(gsl_rand_rng, sigma) + mean;
}

void rand_uniform(double start, double end, double &rand_n, const gsl_rng *gsl_rand_rng)
{
	rand_n = gsl_ran_flat(gsl_rand_rng, start, end);
}





