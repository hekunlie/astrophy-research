#include<FQlib.h>
void creat_h5_group_(const char *filename, const char *set_name, const bool trunc, int row, int col)
{
	int i, j, m, count = 0, s_count = 0;

	int *slash;
	char *names;

	hid_t file_id, group_id;
	herr_t status;
	unsigned rank = 2;
	hsize_t dims[2];
	hid_t dataspace_id;
	hid_t dataset_id;

	if (file_exist(filename))
	{
		if (trunc == TRUE)
		{
			remove(filename);
		}
		else
		{
			file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
		}
		std::cout << "FOUND" << std::endl;
	}
	else
	{
		file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		std::cout << "NOT FOUND" << std::endl;
	}

	// the length of the "set_name"
	for (count = 0; ; count++)
	{
		if (set_name[count] == '\0')
		{
			break;
		}
	}

	slash = new int[count];
	names = new char[count + 1];


	for (i = 0; i < count; i++)
	{
		if (set_name[i] == '/')
		{
			slash[s_count] = i;
			s_count++;
		}
	}
	slash[s_count] = count;

	for (i = 0; i < s_count; i++)
	{

		for (j = 0; j < slash[i + 1]; j++)
		{
			names[j] = set_name[j];
		}
		// label the end
		names[j] = '\0';

		status = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
		status = H5Lget_info(file_id, names, NULL, H5P_DEFAULT);
		std::cout << names << " " << status << std::endl;
		// if the group doesn't exist, create it
		if (status != 0)
		{
			if (i == 0)
			{
				group_id = H5Gcreate(file_id, names, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			}
			else
			{
				group_id = H5Gcreate(group_id, names, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			}
			std::cout << "creat_group " << names << std::endl;
		}
		// if it exist, open it
		else
		{
			if (i == 0)
			{
				group_id = H5Gopen(file_id, names, H5P_DEFAULT);
			}
			else
			{
				group_id = H5Gopen(group_id, names, H5P_DEFAULT);
			}
		}

	}

	dims[0] = row;
	dims[1] = col;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	
	dataset_id = H5Dcreate(group_id, set_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);

	status = H5Gclose(group_id);

	status = H5Eset_auto(H5E_DEFAULT, (H5E_auto2_t)H5Eprint2, stderr);
	status = H5Fclose(file_id);
	delete[] slash;
	delete[] names;
}
void write_h5_(const char *filename, const char *set_name,const double*arr, const int row, const int column)
{
	hid_t file_id;
	herr_t status;

	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

	unsigned rank = 2;
	hsize_t dims[2];
	dims[0] = row;
	dims[1] = column;
	hid_t dataspace_id;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	hid_t dataset_id;
	dataset_id = H5Dcreate(file_id, set_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);

	status = H5Fclose(file_id);

}

int main(int argc, char *argv[])
{
	char name[100], set_name[100];
	sprintf(name, "test.hdf5");
	sprintf(set_name, "/a");
	//creat_h5_group_(name, set_name, TRUE,5,1);
	//sprintf(set_name, "/a/b/c/d/e");
	//creat_h5_group(name, set_name, FALSE);
	//sprintf(set_name, "/a/b/c/f/h");
	//creat_h5_group(name, set_name, FALSE);
	double data[5]{};
	write_h5(name, set_name, data, 5, 1);

	return 0;
}