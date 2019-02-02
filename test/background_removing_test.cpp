#include"FQlib.h"

main(int argc, char* argv[])
{
	int i, j, k,m,n;
	int order = std::atoi(argv[1]);
	int img_num = std::atoi(argv[2]);
	std::cout << order << ", " << img_num << std::endl;
	int num = 40000, terms;
	terms = (order + 1)*(order + 2) / 2;
	double *x = new double[num] {};
	double *y = new double[num] {};
	double *fxy = new double[num] {};
	double *coeff = new double[terms] {};
	char path[150], set1[20], set2[20], set3[20];
	sprintf(set1, "/x");
	sprintf(set2, "/y");
	sprintf(set3, "/fxy");

	double temp = 0;
	double t1, t2;
	/*double *a1 = new double[20]{ 0.22822754, 0.42790551, 0.06420378, 0.75472387, 0.9865212,  0.21199389,
 0.79058317, 0.94770809, 0.58544705, 0.60150541, 0.78670175, 0.73167472,
 0.39444194, 0.60943937, 0.14885859, 0.28801111, 0.59735845, 0.10567613,
 0.2499499,  0.30942367 };
	double *a2 = new double[25]{ 0.32453876, 0.14855579, 0.60610805, 0.61943799, 0.89532018, 0.2736469,
 0.82622208, 0.10691952, 0.10825967, 0.4469363,  0.33545435, 0.24134618,
 0.98546919, 0.27013476, 0.50424156, 0.69110948, 0.91134749, 0.09294967,
 0.49239249, 0.00615354, 0.72758652, 0.4069141,  0.44794632, 0.47482358,
 0.50322338 };
	double *a3 = new double[20]{};
	matrix_product(a1, 4, 5, 5, a2, a3);
	for (i = 0; i < 20; i++)
	{
		std::cout << a3[i] << ", ";
	}



	delete[] a1;
	delete[] a2;
	delete[] a3;
*/
	for (i = 0; i < img_num; i++)
	{
		initialize_arr(x, num);
		initialize_arr(y, num);
		initialize_arr(fxy, num);
		initialize_arr(coeff, terms);
		
	
		sprintf(path, "/home/hkli/work/cpp/test/data/%d.hdf5", i);
		read_h5(path, set1, x);
		read_h5(path, set2, y);
		read_h5(path, set3, fxy);
		temp = 0;
		//for (m = 0; m < 20; m++)
		//{
		//	for (n = 0; n < 20; n++)
		//	{
		//		std::cout << x[m * 200 + n] << ", ";
		//	}
		//	std::cout << std::endl;
		//}

		//for (j = 0; j < num; j++)
		//{
		//	temp += fxy[j];

		//}
		//std::cout << "SUM: " << temp << std::endl;
		t1 = clock();
		poly_fit_2d(x, y, fxy, num, order, coeff);
		t2 = clock();
		temp += t2 - t1;
		for (j = 0; j < terms; j++)
		{
			std::cout << coeff[j] << ",  ";
		}
		std::cout << std::endl;
	}
	std::cout << temp / CLOCKS_PER_SEC;

	delete[] x;
	delete[] y;
	delete[] fxy;
	delete[] coeff;
	return 0;
}