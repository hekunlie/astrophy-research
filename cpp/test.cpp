#include <cmath>
#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include "mpi.h"
#include "FQlib.h"
#include<hdf5.h>
#include<stdio.h>
#include<string>
#include<algorithm>
#include<cstdio>
using namespace std;

int com(const void *a, const void *b)
{
	return (*(double *)a - *(double *)b);
}
int main(int argc, char*argv[])
{	
	char buffer[100], set[20];
	sprintf(buffer, "test.hdf5");
	sprintf(set, "/data");
	double t1, t2, t = 0;
	
	for (int i = 0; i < 1; i++)
	{
		double dou[10]= { 0,1,8,7,6,5,4,3,2,1 };
		float fl[10] = { 0,1,8,7,6,5,4,3,2,1 };
		int in[10] = { 0,1,8,7,6,5,4,3,2,1 };
		qsort_double(dou, 10, 2);
		qsort_float(fl, 10, 2);
		qsort_int(in, 10, 2);
		for (int j = 0; j < 10; j++)
		{
			cout << "DOUBLE: "<<dou[j]<<" FLOAT: "<<fl[j]<<" INT: "<<in[j] <<endl;
		}

	}
	cout << "total: "<<t / CLOCKS_PER_SEC << endl;
	return 0;
}
