#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
//#include<Eigen/Core>
//#include<Eigen/Dense>
#include<FQlib.h>

using namespace std;
//para all_paras;
//using namespace Eigen;
struct MyStruct
{
	double t1, t2, t3, t4;
};
int main()
{
	//MyStruct myst;
	//MatrixXd mat(5, 5),mat1(5,5),mat2;

	para all_paras;
	all_paras.img_size = 90;
	//mat.setZero();
	//mat1.setZero();
	//mat2.setZero();
	//void poly_fit(int *x, int*y, double *f_vals,int order, int num);
	all_paras.t1 = 0.;
	all_paras.t2 = 0.;
	all_paras.t3 = 0.;
	double coeffs[150 * 25]{};
	double st, ed;
	int i, j;
	double img[90 * 90], pimg[90 * 90], fimg[90 * 90], psf[90 * 90], ppsf[90 * 90];
	char buffer[100],setname[30];
	double one[5]{ 1,1,1,1,1 };
	double two[5]{ 2,2,2,2,2 };
	double y;
	//for (j = 0; j < 3; j++)
	//{
	//	y = cblas_ddot(5, one, 1, two, 1);
	//	cout << y << endl;
	//	for (i = 0; i < 5; i++)
	//	{
	//		cout << one[i] << " " << two[i] << endl;
	//	}
	//}
	sprintf(buffer, "coeffs.hdf5");
	sprintf(setname, "/data");
	read_h5(buffer, setname, coeffs, NULL, NULL, NULL, NULL);

	sprintf(buffer, "/home/hkli/gal.fits");
	read_img(img,buffer);
	pow_spec(img, pimg, 90, 90);

	sprintf(buffer, "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/psf.fits");
	read_img(psf, buffer);
	pow_spec(psf, ppsf, 90, 90);
	get_psf_thres(ppsf, &all_paras);
	int max = 0;
	for (i = 0; i < 90 * 90; i++)
	{
		if(ppsf[i] > all_paras.psf_pow_thres)
		{
			max++;
		}
	}
	cout << max << endl;
	cout <<"PSF THRES: "<< all_paras.psf_pow_thres << endl;

	sprintf(buffer, "!/home/hkli/galp.fits");
	write_img(pimg, 90, 90, buffer);

	sprintf(buffer, "!/home/hkli/psfp.fits");
	write_img(ppsf, 90, 90, buffer);
	st = clock();
	for (i = 0; i < 20000; i++)
	{
		smooth(pimg, fimg, ppsf, coeffs, &all_paras);
		if (i == 0)
		{
			sprintf(buffer, "!/home/hkli/galpf.fits");
			write_img(fimg, 90, 90, buffer);
		}
	}


	ed = clock();
	cout << (ed - st) / CLOCKS_PER_SEC << endl;
	//cin >> i;
	cout << all_paras.t1 / CLOCKS_PER_SEC << " " << all_paras.t2 / CLOCKS_PER_SEC << " " << all_paras.t3 / CLOCKS_PER_SEC << endl;
	return 0;
}
//void poly_fit(int *x, int *y, double*f_vals, int order, int num)
//{
//	/* the expression to m-th-order polynomials: SUM_{n: 0~m} SUM_{i: 0~n} x^{n-i}*y^i */
//	/* then it's easy to write down the matrix of the equations come from the least square method*/
//	//double a1, a2, a3, a4, a5, a6;
//	int i, j, k, turns, px, py, z;
//	turns = (order + 1)*(order + 2)*0.5;
//	Array<double,Dynamic, Dynamic> xy_pow(2,turns),fx(1,num),fy(1,num),fxy(1,num), fx2(1, num), fy2(1, num);
//	Array<double,Dynamic, Dynamic> coeffs(turns, turns),chi_xy(turns,num);
//	Array<double, Dynamic, Dynamic> fvals(1, num), fz(1,turns);
//
//	for (i = 0; i < num; i++)
//	{
//		fx(0, i) = x[i];
//		fy(0, i) = y[i];
//		fvals(0, i) = f_vals[i];
//	}
//	
//	/* the power of x and y of each row of the matrix */
//	k = 0;
//	for (i = 0; i <= order; i++)
//	{
//		for (j = 0; j <= i; j++)
//		{
//			xy_pow(0, k) = i - j; //power of x
//			xy_pow(1, k) = j;    // power of y
//			if (i == j)
//			{
//				if (j == 0)
//				{
//					chi_xy.block(k, 0, 1, num).setOnes();
//				}
//				else if (j == 1)
//				{
//					chi_xy.block(k, 0, 1, num) = fy;
//				}
//				else
//				{
//					chi_xy.block(k, 0, 1, num) = fy.pow(j);
//				}
//			}
//			else if (i ==j+1)
//			{
//				if (j == 0)
//				{
//					chi_xy.block(k, 0, 1, num) = fx;
//				}
//				else if (j == 1)
//				{
//					chi_xy.block(k, 0, 1, num) = fx*fy;
//				}
//				else
//				{
//					chi_xy.block(k, 0, 1, num) = fx*fy.pow(j);
//				}
//			}
//
//			else
//			{
//				chi_xy.block(k, 0, 1, num) = fx.pow(i - j)*fy.pow(j);
//			}
//
//			k++;
//		}
//	}
//	//cout << "X: "<<fx1 << endl;
//	//cout << "Chi_x: " << endl<< chi_x << endl;
//	//cout << "Y: "<<fy1 << endl;
//	//cout << "Chi_y: "<<endl<<chi_y << endl;
//	//cout << "Chi_xy: " << chi_ << endl;
//	//cout << xy_pow << endl;
//	for (i = 0; i < turns; i++)
//	{	
//		px = xy_pow(0, i);
//		py = xy_pow(1, i);
//
//		if (px == 0)
//		{
//			fx2.setOnes();
//		}
//		else if (px == 1)
//		{
//			fx2 = fx;
//		}
//		else
//		{
//			fx2 = fx.pow(px);
//		}
//
//		if (py == 0)
//		{
//			fy2.setOnes();
//		}
//		else if (py == 1)
//		{
//			fy2 = fy;
//		}
//		else
//		{
//			fy2 = fy.pow(py);
//		}
//		fxy = fx2*fy2;
//		fz(0, i) = (fvals*fxy).sum();
//
//		for (j = i; j < turns; j++)
//		{	
//			if (i == 0 && j == 0)
//			{
//				coeffs(i, j) = num;
//			}
//			else
//			{
//				coeffs(i, j) = (chi_xy.block(j, 0, 1, num)*fxy).sum();
//				coeffs(j, i) = coeffs(i, j);
//			}
//		}
//	}
//	//cout << coeffs << endl;
//}

/*sprintf(chip_path, "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/0/gal_chip_000%d.fits", myid);
read_img(gals,chip_path);
sprintf(buffer, "Rank: %d, open gal_chip_000%d.fits", myid, myid);
cout << buffer << endl;
for (i = 0; i < stamp_num; i++)
{
segment(gals, gal, i, size, stamp_nx, stamp_nx);

pow_spec(gal, ppow, size, size);
stack(n_fit_pow, ppow, i, size, stamp_nx, stamp_nx);

paraboloid_fit(ppow, noise, &all_paras, size);
stack(fit_pow, noise, i, size, stamp_nx, stamp_nx);

initialize(gal, size*size);
initialize(ppow, size*size);
initialize(noise, size*size);
}
sprintf(chip_path, "!/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/gal_chip_000%d_pow.fits", myid);
for (i = 0; i < stamp_num*size*size; i++)
{
n_fit_pow[i] = log10(n_fit_pow[i]);
}
write_img(n_fit_pow, 100 * size, 100 * size, chip_path);
sprintf(buffer, "Rank: %d, write gal_chip_000%d_pow.fits", myid, myid);
cout << buffer << endl;
sprintf(chip_path, "!/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/gal_chip_000%d_fit_pow.fits", myid);
write_img(fit_pow, 100 * size, 100 * size, chip_path);
sprintf(buffer, "Rank: %d, write gal_chip_000%d_fit_pow.fits", myid, myid);
cout << buffer << endl;
ed = clock();*/