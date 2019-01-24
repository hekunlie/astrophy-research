#include "FQlib.h"

main()
{
	char path[100], set_name[20];
	int num = 100000, i, j, k, bin_num = 10;

	double *ddata = new double[num]{};
	float *fdata = new float[num]{};
	int *idata = new int[num]{};

	double *dbin = new double[bin_num + 1]{};
	float *fbin = new float[bin_num + 1]{};
	int *ibin = new int[bin_num + 1]{};

	sprintf(path, "/home/hkli/work/test/test.hdf5");
	sprintf(set_name, "/ddata");
	read_h5(path, set_name, ddata);

	sprintf(set_name, "/fdata");
	read_h5(path, set_name, fdata);

	sprintf(set_name, "/idata");
	read_h5(path, set_name, idata);
	

	std::cout << ddata[0] << ",  " << fdata[0] << std::endl;
	ddata[0] = fdata[0];
	std::cout << ddata[0] << ",  " << fdata[0] << std::endl;

	std::cout << ddata[1] << ",  " << fdata[1] << std::endl;
	fdata[1] = ddata[1];
	std::cout << ddata[1] << ",  " << fdata[1] << std::endl;

	// sort() test
	sort_arr(ddata, num, 1);
	sort_arr(fdata, num, 1);
	sort_arr(idata, num, 1);

	for (i = 0; i < num-1; i++)
	{
		if (ddata[i] > ddata[i + 1])
		{
			std::cout << "DOUBLE SORT() FAILD: " << ddata[i] << " " << ddata[i + 1] << std::endl;
		}
		if (fdata[i] > fdata[i + 1])
		{
			std::cout << "FLOAT SORT() FAILD: " << fdata[i] << " " << fdata[i + 1] << std::endl;
		}
		if (idata[i] > idata[i + 1])
		{
			std::cout << "INT SORT() FAILD: " << idata[i] << " " << idata[i + 1] << std::endl;
		}
	}
	for (i = 0; i < 20; i++)
	{
		//std::cout << ddata[i] << ", " << fdata[i] << ", " << idata[i] << std::endl;
	}
	std::cout << "DESCEND" << std::endl;
	sort_arr(ddata, num, 2);
	sort_arr(fdata, num, 2);
	sort_arr(idata, num, 2);
	for (i = 0; i < num - 1; i++)
	{
		if (ddata[i] < ddata[i + 1])
		{
			std::cout << "DOUBLE SORT() FAILD: " << ddata[i] << " " << ddata[i + 1] << std::endl;
		}
		if (fdata[i] < fdata[i + 1])
		{
			std::cout << "FLOAT SORT() FAILD: " << fdata[i] << " " << fdata[i + 1] << std::endl;
		}
		if (idata[i] < idata[i + 1])
		{
			std::cout << "INT SORT() FAILD: " << idata[i] << " " << idata[i + 1] << std::endl;
		}
	}
	for (i = 0; i < 20; i++)
	{
		std::cout << ddata[i] << ", " << fdata[i] << ", " << idata[i] << std::endl;
	}
	// set_bin() test
	set_bin(ddata, num, dbin, bin_num, 1.2);
	set_bin(fdata, num, fbin, bin_num, 1.2);
	set_bin(idata, num, ibin, bin_num, 2);

	for (i = 0; i < bin_num + 1; i++)
	{
		std::cout << dbin[i] << ",  ";
	}
	std::cout << std::endl;
	for (i = 0; i < bin_num + 1; i++)
	{
		std::cout << fbin[i] << ",  ";
	}
	std::cout << std::endl;
	for (i = 0; i < bin_num + 1; i++)
	{
		std::cout << ibin[i] << ",  ";
	}
	std::cout << std::endl;



	// histogram() test


	delete[] dbin;
	delete[] fbin;
	delete[] ibin;
	delete[] ddata;
	delete[] fdata;
	delete[] idata;
	return 0;
}