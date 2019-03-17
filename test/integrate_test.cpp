#include<stdlib.h>
#include<iostream>
#include<ctime>
#include<cmath>
#include<FQlib.h>

void distance(double low_z, double high_z, double precision_thresh, double omg_m, double omg_lam, double &result)
{
	int i, j, k, num;
	int max_run = 30;
	double result_1, result_2, d_res;
	double dz, z;
	double *fbin1, *fbin2;
	int tag1, iters=1;

	num = 20;
	
	fbin1 = new double[int(pow(2,max_run)*(num-1)+1)];
	tag1 = 1;

	result_1 = 0;
	for (i = 0; i < max_run; i++)
	{		
		dz = (high_z - low_z) / num;
		
		if (iters == 0)
		{
			for (j = 0; j < num; j++)
			{
				z = low_z + j * dz;
				fbin1[j] = 1. / sqrt(pow(1 + z, 3)*omg_m + omg_lam);
			}
		}
		else
		{
			for (j = 0; j < num - 1; j++)
			{
				k = 2 * j + 1;
				z = low_z + k * dz;
				fbin1[k] = 1. / sqrt(pow(1 + z, 3)*omg_m + omg_lam);				
			}
		}

		result_2 = 0;
		for (j = 0; j < num-1; j++)
		{
			result_2 += fbin1[j] + fbin1[j+1];
		}
		result_2 = result_2 * dz*0.5;
	
		if (fabs(result_2 - result_1) <= precision_thresh)
		{
			if (tag1 == 1)
			{
				delete[] fbin1;
			}
			result = (result_2 + result_1) *0.5;
			break;
		}
		else
		{
			iters += 1;
			result_1 = result_2;

			num = num + num - 1;
			fbin2 = new double[num];
			for (j = 0; j = num; j++)
			{
				fbin2[j*2] = fbin1[j];
			}
			delete[] fbin1;
			tag1 = 0;
			fbin1 = new double[num];
			tag1 = 1;
			for (j= 0; j= num; j++)
			{
				fbin1[j] = fbin2[j];
			}
			delete[] fbin2;			
		}
	}
	if (i == max_run)
	{
		if (tag1 == 1)
		{
			delete[] fbin1;
		}
		std::cout << "Max iteration, doesn't converge, " << result_2 - result_1 << "." << std::endl;
		exit(0);
	}
}
int main(int argc, char *argv[])
{
	double a1, a2, z1,z2;
	double t1, t2;
	int num;
	char data_path[150];
	char set_name[50];
	int i, j, k;
	int tag;
	double redshift[1001]{};
	double distance[1001]{};
	char log_inform[200];
	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/redshift.hdf5");
	sprintf(set_name, "/redshift");
	read_h5(data_path, set_name, redshift);
	sprintf(set_name, "/distance");
	read_h5(data_path, set_name, distance);

	
	z1 = atof(argv[1]);
	z2 = atof(argv[2]);

	a1 = 1. / (z1 + 1);
	a2 = 1. / (z2 + 1);
	
	num = atoi(argv[3]);
	std::cout << z1 << " " << z2 << std::endl;
	double *bins = new double[num];
	log_bin(z1+0.1, z2+1, num, bins);
	//show_arr(bins, 1, num);
	for ( i = 0; i < num; i++)
	{
		std::cout << log10(bins[i]) << ", ";
	}
	std::cout << std::endl;
	double r,r1;
	t1 = clock();
	for ( i = 0; i < 1; i++)
	{
		com_distance(z1, z2, 0.000001, 0.31, 0.69, r);
		//distance(z1, z2, 0.00001, 0.31, 0.69, r1);
		find_near(redshift, z2, 1001, tag);
	}
	t2 = clock();
	sprintf(log_inform, "z1: %.3f, z2: %.3f. distance(z1, z2): %.7f. found z: %d %.3f (%.2f), distance(z1, z): %.7f. time: %4f sec", z1, z2, r, tag, redshift[tag], z2, distance[tag], (t2 - t1) / CLOCKS_PER_SEC);
	std::cout <<log_inform << std::endl;



	delete[] bins;
	
	return 0;
}