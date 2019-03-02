#include<FQlib.h>

int main(int argc, char *argv[])
{
	char file_name[100], set_name[100], attrs_name[50];;
	char logs[200];
	sprintf(file_name, "test.hdf5");
	sprintf(set_name, "/a/b/c");
	creat_h5_group(file_name, set_name, TRUE);
	//sprintf(set_name, "/a/b");
	//creat_h5_group(name, set_name, FALSE);
	//sprintf(set_name, "/a/b/c/c1");
	//creat_h5_group(name, set_name, FALSE);
	//sprintf(set_name, "/a/b/cc");
	//creat_h5_group(name, set_name, FALSE);
	int i, num = 3,a_num=5;
	double *data= new double[num]{}, a1=0;
	float *fdata = new float[num]{}, a2 = 0;
	int *idata = new int[num]{}, a3 = 0;
	long *ldata= new long[num]{}, a4 = 0;

	double *attr_d = new double[a_num]{};
	int *attr_i = new int[a_num] {};

	for (i = 0; i < a_num; i++)
	{
		attr_d[i] = i;
		attr_i[i] = i;

	}

	for (i = 0; i < num; i++)
	{
		data[i] = i;
		fdata[i] = i;
		idata[i] = i;
		ldata[i] = i;
	}
	for (i = 0; i < num; i++)
	{
		a1 += data[i];
		a2 += fdata[i];
		a3 += idata[i];
		a4 += ldata[i];
	}
	//sprintf(logs, "%.8f, %.8f, %d, %ld", a1, a2, a3, a4);
	//std::cout << logs << std::endl;

	sprintf(set_name, "/a/b/c");
	sprintf(attrs_name, "group attrs");
	write_h5_attrs(file_name, set_name, attrs_name, attr_i, a_num, "g");

	sprintf(set_name, "/a/b/c/c1");
	write_h5(file_name, set_name, data, num, 1, FALSE);
	sprintf(attrs_name, "double attrs");
	write_h5_attrs(file_name, set_name, attrs_name, attr_d, a_num, "d");

	sprintf(set_name, "/a/b/c/c2");
	write_h5(file_name, set_name, fdata, num, 1, FALSE);

	sprintf(set_name, "/a/b/c/c3");
	write_h5(file_name, set_name, idata, num, 1, FALSE);
	sprintf(attrs_name, "int attrs");
	write_h5_attrs(file_name, set_name, attrs_name, attr_i, a_num, "d");

	sprintf(set_name, "/a/b/cc");
	write_h5(file_name, set_name, ldata, num, 1, FALSE);

	for (i = 0; i < num; i++)
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
	read_h5(file_name, set_name, data);
	sprintf(set_name, "/a/b/c/c2");
	read_h5(file_name, set_name, fdata);
	sprintf(set_name, "/a/b/c/c3");
	read_h5(file_name, set_name, idata);
	sprintf(set_name, "/a/b/cc");
	read_h5(file_name, set_name, ldata);
	for (i = 0; i < num; i++)
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
