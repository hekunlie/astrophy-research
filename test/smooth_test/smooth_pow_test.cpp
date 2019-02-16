#include <cstring>
#include<sstream>
#include <ctime>
#include <stdlib.h>
#include "mpi.h"
#include<hdf5.h>
#include<stdio.h>
#include<string>
#include <FQlib.h>

int main()
{
	para all_paras;
	int size = 64, i, j, k;
	double *gal = new double[size*size]{};
	double *pgal = new double[size*size]{};
	double *img = new double[size*size * 10000]{};
	double *data = new double[60000]{};
	all_paras.stamp_size = size;

	double temp_flux = 0;
	char img_path[150], set_name[20];
	
	sprintf(img_path, "/mnt/perc/hklee/simu_test/0/gal_chip_0000.fits");
	read_fits(img_path, img);
	
	for (i = 0; i < 10000; i++)
	{
		segment(img, gal, i, size, 100, 100);
		pow_spec(gal, pgal, size, size);
		snr_est(pgal, &all_paras, 2);
		temp_flux = 0;
		
		for (k = 0; k < size*size; k++)
		{
			temp_flux += gal[k];
		}

		data[i * 6 + 0] = pgal[size / 2 + size * size / 2];
		data[i * 6 + 1] = all_paras.gal_flux2;
		data[i * 6 + 2] = all_paras.gal_flux_alt;
		data[i * 6 + 3] = all_paras.gal_flux2_ext[0];
		data[i * 6 + 4] = all_paras.gal_flux2_ext[1];
		data[i * 6 + 5] = fabs(temp_flux);
	}
	sprintf(img_path, "/home/hklee/work/data.hdf5");
	sprintf(set_name, "/data");
	write_h5(img_path, set_name, data, 10000, 6);

	delete[] data;
	delete[] gal;
	delete[] pgal;
	delete[] img;
	return 0;
}