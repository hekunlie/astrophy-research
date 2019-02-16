#include<FQlib.h>
//void mgrid(const int size_y, const int size_x, int *arr)
//{
//	int i, j;
//	for (i = 0; i < size_y; i++)
//	{
//		for (j = 0; j < size_x; j++)
//		{
//			arr[i*size_x+j] 
//		}
//	}
//}

void trim_img(const double *img, const int size_y, const int size_x, const int edge, double *img_t)
{
	// cfht: the widthes of the useless pixel around the image is (2, 32, 32, 34) 
	//			for (bottom(y->0), left(x->0), right(x->size_x), top(y->size_y))
	// the same witd of border for segmentation

	int i, j, m, n, k=0;
	for (i = edge; i < size_y-edge; i++)
	{	
		m = i * size_x;
		for (j = edge; j < size_x-edge; j++)
		{
			n = m + j;
			img_t[k] = img[n];
			k++;
		}
	}
}

void background_remove(const double *img, const int size_y, const int size_x, const int pts_num, double *img_sub)
{	
	int i, j, k, m, n, num = size_x*size_y;
	double thres_1, thres_2, pix;
	int order = 1;
	int terms = (order + 1)*(order + 2) / 2;
	double *coeffs = new double[terms];
	double *sample = new double[pts_num*3] {};
	int *ch_x = new int[pts_num*3] {};
	int *ch_y = new int[pts_num*3] {};
	int *ch = new int[num]{};

	double *fit_vals = new double[pts_num] {};
	double *fit_x = new double[pts_num] {};
	double *fit_y = new double[pts_num] {};

	for (i = 0; i < num; i++)
	{
		ch[i] = i;
	}
	rand_shuffle(ch, num);
	std::cout << "SHUFFLE" << std::endl;

	for (i = 0; i < pts_num*3; i++)
	{
		m = ch[i] / size_x;
		n = ch[i] % size_x;
		sample[i] = img[m*size_x + n];
		ch_y[i] = m;
		ch_x[i] = n;
	}
	sort_arr(sample, 3 * pts_num, 1);
	thres_1 = sample[pts_num];
	thres_2 = sample[pts_num*2];
	
	std::cout << thres_1 << ", " << thres_2 << std::endl;
	char path[100], set[20];
	sprintf(path, "/home/hkli/work/cpp/test/ch.hdf5");
	sprintf(set, "/data");
	write_h5(path, set, sample, 3*pts_num, 1);
	sprintf(path, "/home/hkli/work/cpp/test/ch_xy.hdf5");
	sprintf(set, "/data");
	write_h5(path, set, ch, num, 1);

	delete[] sample;// saving memory
	k = 0;
	for (i = 0; i < pts_num * 3; i++)
	{	
		pix = img[ch_y[i]*size_x + n];
		if ( k<pts_num)
		{
			if (pix >= thres_1 and pix < thres_2)
			{
				fit_vals[k] = pix;
				fit_y[k] = ch_y[i];
				fit_x[k] = ch_x[i];
				k++;
			}
		}
		else
		{
			k--;
			break;
		}
	}
	// check
	if (k !=pts_num-1)
	{
		char buf[80];
		sprintf(buf, "Background Removing: pixel number is %d (must be %d!!!) ", k+1, pts_num);
		std::cout << buf << std::endl;
		exit(0);
	}

	poly_fit_2d(fit_x, fit_y, fit_vals, pts_num, order, coeffs);
	for (i = 0; i < size_y; i++)
	{
		m = i * size_x;
		for (j = 0; j < size_x; j++)
		{
			n = m + j;
			img_sub[n] = img[n] - fval_at_xy(j, i, order, coeffs);
		}
	}
	for (i = 0; i < terms; i++)
	{
		std::cout << coeffs[i] << ", ";
	}

	delete[] coeffs;
	delete[] fit_x;
	delete[] fit_y;
	delete[] ch_x;
	delete[] ch_y;
	delete[] ch;
}


int main(int argc, char *argv[])
{
	gsl_initialize(1230);
	int i, j, k, m, n;
	int edge = 34;
	int imgx = 2112;
	int imgy = 4644;
	double *img = new double[imgx*imgy]{};
	double *img_t = new double[(imgx - edge * 2)*(imgy - edge * 2)]{};
	double *img_sub = new double[(imgx - edge * 2)*(imgy - edge * 2)]{};
	char img_path[100];
	sprintf(img_path, "/home/hkli/work/cpp/test/data/img.fits");
	read_fits(img_path, img);

	trim_img(img, imgy, imgx, edge, img_t);
	sprintf(img_path, "!/home/hkli/work/cpp/test/data/img_t.fits");
	write_fits(img_path, img_t, imgy-2*edge, imgx - 2 * edge);

	std::cout << "READ" << std::endl;
	background_remove(img_t, imgy - 2 * edge, imgx - 2 * edge, 50000, img_sub);

	sprintf(img_path, "!/home/hkli/work/cpp/test/data/img_sub.fits");
	write_fits(img_path, img_sub, imgy, imgx);
	delete[] img;
	delete[] img_t;
	delete[] img_sub;
	gsl_free();
	return 0;
}