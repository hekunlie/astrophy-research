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
	int size = 50;
	para paras;
	paras.img_x = size;
	paras.img_y = size;
	paras.detect_thres = 0.3;
	paras.area_thres = 5;
	paras.img_size = size;
	
	double *img = new double[size*size]{};
	double* img_t = new double[size*size]{};
	int *s_x = new int[size*size]{};
	int *s_y = new int[size*size]{};
	double *s_c = new double[7*paras.max_source]{};

	char buffer[60];
	cout << "starting..." << endl;
	sprintf(buffer, "/home/hkli/temp/img.fits");
	read_img(img, buffer);
	cout << "read img..." << endl;
	//for (int k = 0; k < 2112 * 4644; k++) 
	//{
	//	img[k] = img[k] - 4000;
	//}
	//detector(img, s_x, s_y, s_c, &paras,FALSE);
	//cout << "finished" << endl;
	//int count = 0, s=0;
	double t1, t2;
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
	/*for (int m = 0; m < 1; m++)
	{
		for (int i = 0; i <paras.max_source; i++)
		{
			
			if (s_c[7*i] != 0)
			{	
				for (int k = 0; k < 7; k++)
				{
					cout << s_c[7 * i+k] << endl;
				}
				for (int j = 0; j < s_c[7*i]; j++)
				{
					img_t[s_x[s + j] + s_y[s + j] * size] = 1;
				}
				s += s_c[7*i];
			}
			else
				break;

		}

		sprintf(buffer, "!/home/hkli/temp/t_img.fits");
		write_img(img_t, size, size, buffer);
	}*/
	t2 = clock();
	cout << (t2 - t1) / CLOCKS_PER_SEC;
	return 0;
}
