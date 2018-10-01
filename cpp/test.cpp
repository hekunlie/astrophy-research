#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include "FQlib.h"

using namespace std;

int main()
{
	int size = 40, s_num, area_s=0, area_e;
	para paras;
	paras.img_x = size;
	paras.img_y = size;
	paras.detect_thres = 4.5;
	paras.stamp_size = size;

	double t1, t2;
	double *img = new double[size*size]{};
	double* img_t = new double[size*size]{};
	int *s_x = new int[size*size]{};
	int *s_y = new int[size*size]{};
	double *s_c = new double[8*paras.max_source]{};
	int detect;
	char buffer[60];
	cout << "starting..." << endl;
	sprintf(buffer, "/home/hkli/temp/test.fits");
	read_img(img, buffer);
	cout << "read img..." << endl;
	s_num = source_detector(img, s_x, s_y, s_c, &paras, false);
	detect = galaxy_finder(img, &paras, false);
	cout << paras.gal_py << " " << paras.gal_px << endl;
	cout << detect<<"finished, find: "<< s_num << endl;
	for (int m = 0; m < s_num; m++)
	{
		cout << "AREA: "<<s_c[8 * m] << endl;
		area_e = area_s + s_c[m * 8];

		for (int i = area_s; i < area_e; i++)
		{
			img_t[s_x[i] + s_y[i] * size] = 1;
		}
		img_t[(int)s_c[8 * m + 1] * size + (int)s_c[8 * m + 2]] += 2;
		area_s = area_e;
	}
	img_t[paras.gal_py*size + paras.gal_px] += 1;
	sprintf(buffer, "!/home/hkli/temp/t_img.fits");
	write_img(img_t, size, size, buffer);
	
	/*double t1, t2;
	t1 = clock();
	f_snr(img, &paras);
	cout << paras.gal_flux2 << " " << paras.gal_flux_alt << endl;
	int x[20]{ -1,  0,  1, -2, -1,  0,  1,  2, -2, -1,  1,  2, -2, -1,  0,  1,  2, -1,  0,  1 };
	int y[20]{ -2, -2, -2, -1, -1, -1, -1, -1,  0,  0,  0,  0,  1,  1,  1,  1,  1,  2, 2, 2 };
	double fz[20], fit_paras[6];
	int xc = size / 2, i;
	cout << xc << endl;
	for (i = 0; i < 20; i++)
	{
		fz[i] = img[(xc + y[i])*size + xc + x[i]];
	}

	hyperfit_5(fz, fit_paras,&paras);
	for (i = 0; i < 6; i++)
	{
		cout << fit_paras[i] << " ";
	}
	cout << endl;
	cout << pow(10, fit_paras[0]) <<" "<<img[xc*size+xc]<< endl;
	/*
	}*/
	t2 = clock();
	cout << (t2 - t1) / CLOCKS_PER_SEC;
	return 0;
}
