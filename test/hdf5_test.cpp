#include<FQlib.h>

int main(int argc, char *argv[])
{
	char name[100], set_name[100];
	char logs[100];
	sprintf(name, "test.hdf5");
	sprintf(set_name, "/a");
	creat_h5_group(name, set_name, TRUE);
	sprintf(set_name, "/a/b");
	creat_h5_group(name, set_name, FALSE);

	int i;
	double data[100000]{}, a1=0;
	float fdata[100000]{}, a2 = 0;
	int idata[100000]{}, a3 = 0;
	long ldata[100000]{}, a4 = 0;

	for (i = 0; i < 100000; i++)
	{
		data[i] = i;
		fdata[i] = i;
		idata[i] = i;
		ldata[i] = i;
	}
	for (i = 0; i < 100000; i++)
	{
		a1 += data[i];
		a2 += fdata[i];
		a3 += idata[i];
		a4 += ldata[i];
	}
	sprintf(logs, "%.8f, %.8f, %d, %ld", a1, a2, a3, a4);
	std::cout << logs << std::endl;
	sprintf(set_name, "/a/b/c/c1");
	write_h5(name, set_name, data, 100000, 1, FALSE);
	sprintf(set_name, "/a/b/c/c2");
	write_h5( name, set_name, fdata, 100000, 1, FALSE);
	sprintf(set_name, "/a/b/c/c3");
	write_h5(name, set_name, idata, 100000, 1,FALSE);
	sprintf(set_name, "/a/b/c4");
	write_h5(name, set_name, ldata, 100000, 1, FALSE);

	for (i = 0; i < 100000; i++)
	{
		data[i] = 0;
		fdata[i] = 0;
		idata[i] = 0;
		ldata[i] = 0;
	}
	a1 = 0;
	a2 = 0;
	a3 = 0;
	a4 = 0;
	sprintf(set_name, "/a/b/c/c1");
	read_h5(name, set_name, data);
	sprintf(set_name, "/a/b/c/c2");
	read_h5(name, set_name, fdata);
	sprintf(set_name, "/a/b/c/c3");
	read_h5(name, set_name, idata);
	sprintf(set_name, "/a/b/c4");
	read_h5(name, set_name, ldata);
	for (i = 0; i < 100000; i++)
	{
		a1 += data[i];
		a2 += fdata[i];
		a3 += idata[i];
		a4 += ldata[i];
	}
	sprintf(logs, "%.8f, %.8f, %d, %ld", a1, a2, a3, a4);
	std::cout << logs << std::endl;

	return 0;
}
