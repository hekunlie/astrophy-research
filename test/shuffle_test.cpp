#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<stdio.h>
#include<ctime>
#include <stdlib.h>
#include <FQlib.h>

int main()
{	
	int seed = 1230, i;
	gsl_rng_initialize(seed);
	int num = 100;
	double *arr1 = new double[num];
	float *arr2 = new float[num];
	int *arr3 = new int[num];
	for (i = 0; i < num; i++)
	{
		arr1[i] = i*1.1;
		arr2[i] = i * 1.2;
		arr3[i] = i;
	}

	rand_shuffle(arr1, num);
	rand_shuffle(arr2, num);
	rand_shuffle(arr3, num);
	for (i = 0; i < 10; i++)
	{
		std::cout << arr1[i] << ",  " << arr2[i] << ",  " << arr3[i] << std::endl;
	}
	sort_double(arr1, num, 1);
	sort_float(arr2, num, 1);
	sort_int(arr3, num, 1);
	std::cout << "SORT:" << std::endl;
	for (i = 0; i < num; i++)
	{
		std::cout << arr1[i] << ",  " << arr2[i] << ",  " << arr3[i] << std::endl;
	}
	gsl_rng_free();
	delete[] arr1;
	delete[] arr2;
	delete[] arr3;
	return 0;
}