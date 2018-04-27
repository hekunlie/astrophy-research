#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include<FQlib.h>

using namespace std;

int main()
{
	double *img = new double[4606 * 2043]{};
	double* img_t = new double[4606 * 2043]{};
	int *s_x = new int[4606 * 2043]{};
	int *s_y = new int[4606 * 2043]{};
	int *s_c = new int[4606 * 2043]{};
	char buffer[60];
	cout << "starting..." << endl;
	sprintf(buffer, "/home/hklee/temp/831549p_1OFCBC.fits");
	read_img(img, buffer);
	cout << "read img..." << endl;
	//for (int k = 0; k < 2112 * 4644; k++) 
	//{
	//	img[k] = img[k] - 4000;
	//}
	detector(img, s_x, s_y, s_c, 45., 4606, 2043);
	cout << "finished" << endl;
	int count = 0, s=0;
	double t1, t2;
	t1 = clock();
	for (int m = 0; m < 10; m++)
	{
		for (int i = 0; i < 4606 * 2043; i++)
		{
			if (s_c[i] != 0)
			{
				for (int j = 0; j < s_c[i]; j++)
				{
					img_t[s_x[s + j] + s_y[s + j] * 2043] = 1;
				}
				s += s_c[i];
			}
		}

		sprintf(buffer, "!/home/hklee/temp/t_img.fits");
		write_img(img_t, 4606, 2043, buffer);
	}
	t2 = clock();
	cout << (t2 - t1) / CLOCKS_PER_SEC;
	return 0;
}
