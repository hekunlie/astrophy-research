#include<FQlib.h>
void creat_h5_group_(const char *filename, const char *set_name, const bool trunc)
{
	int i, j, m, count = 0, s_count = 0;

	int *slash;
	char *name;

	hid_t file_id, group_id;
	herr_t status;
	unsigned rank = 2;
	hsize_t dims[2];
	hid_t dataspace_id;
	hid_t dataset_id;

	// the length of the "set_name"
	for (count = 0; ; count++)
	{
		if (set_name[count] == '\0')
		{
			break;
		}
	}
	slash = new int[count];
	name = new char[count + 1];

	// find the position of "/"
	for (i = 0; i < count; i++)
	{
		if (set_name[i] == '/')
		{
			slash[s_count] = i;
			s_count++;
		}
	}
	// the position of '\0', the end
	slash[s_count] = count;

	if (file_exist(filename))
	{
		// trucate it
		if (trunc == TRUE)
		{
			remove(filename);
			std::cout << "REMOVE" << std::endl;
			file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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

	for (i = 0; i < s_count; i++)
	{

		for (j = 0; j < slash[i + 1]; j++)
		{
			name[j] = set_name[j];
		}
		name[j] = '\0';// label the end
		status = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
		status = H5Lget_info(file_id, name, NULL, H5P_DEFAULT);
		if (status != 0)
		{
			group_id = H5Gcreate(file_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Gclose(group_id);
			std::cout << "creat_group " << name << std::endl;
		}
		else
		{
			std::cout << "group exit " << name << std::endl;
			continue;
		}
	}

	status = H5Eset_auto(H5E_DEFAULT, (H5E_auto2_t)H5Eprint2, stderr);
	status = H5Fclose(file_id);
	delete[] slash;
	delete[] name;

}
void write_h5_(const char *filename, const char *set_name,const double*arr, const int row, const int column)
{
	int i, count, s_count;

	for (count = 0; ; count++)
	{
		if (set_name[count] == '\0')
		{
			break;
		}
	}

	int *slash = new int[count] {};
	char *names = new char[count + 1];
	char *new_names = new char[count + 1];
	// find the slashes in the set_name
	s_count = 0;
	for (i = 0; i < count; i++)
	{
		if (set_name[i] == '/')
		{
			slash[s_count] = i;
			s_count++;
		}
	}
	// label the end of the set_name
	slash[s_count] = count;
	// the set_name is like /a/b/c
	// the dataset will be under c,while /a/b must be created before
	// 
	// the "c" is in the "new_names",  "/a/b" in "names"
	for (i = 0; i < slash[s_count - 1]; i++)
	{
		names[i] = set_name[i];
	}
	names[slash[s_count - 1] + 1] = '\0';

	for (i = slash[s_count - 1]+1; i < slash[s_count]; i++)
	{
		new_names[i - slash[s_count - 1]-1] = set_name[i];
	}
	new_names[slash[s_count] - slash[s_count - 1]-1] = '\0';	
	// try to create /a/b if it doesn't exist
	creat_h5_group(filename, names, FALSE);	

	hid_t file_id, group_id;
	herr_t status;

	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

    group_id = H5Gopen1(file_id, names);

	unsigned rank = 2;
	hsize_t dims[2];
	dims[0] = row;
	dims[1] = column;
	hid_t dataspace_id;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	hid_t dataset_id;

	dataset_id = H5Dcreate(group_id, new_names, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);

    status = H5Gclose(group_id);
	status = H5Fclose(file_id);

	delete[] slash;
	delete[] names;
	delete[] new_names;

}

int main(int argc, char *argv[])
{
	hid_t file_id, group_id_1, group_id_2;
	herr_t status;
	char name[100], set_name[100];
	sprintf(name, "test.hdf5");
	sprintf(set_name, "/a/b/c");
	creat_h5_group_(name, set_name, FALSE);
	sprintf(set_name, "/a/b/c/d/e/f");
	creat_h5_group_(name, set_name, FALSE);
	//file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	//group_id_1 = H5Gcreate(file_id, set_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//status = H5Gclose(group_id_1);
	//std::cout << "CREATE /a" << std::endl;
	//sprintf(set_name, "/a/b");
	//group_id_1 = H5Gcreate(file_id, set_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//status = H5Gclose(group_id_1);
	//status = H5Fclose(file_id);
	//std::cout << "CREATE /a/b" << std::endl;
	//sprintf(set_name, "/a11/b2/c3");
	//creat_h5_group(name, set_name, TRUE,5,1);
	//sprintf(set_name, "/a/b/c/d/e");
	//creat_h5_group(name, set_name, FALSE);
	//sprintf(set_name, "/a11/b2/c3/dddd");
	//creat_h5_group(name, set_name, FALSE);
	double data[5]{};
	//write_h5_( name, set_name, data, 5, 1);
	//sprintf(set_name, "/a11/b2/c3/ffff");
	//creat_h5_group(name, set_name, FALSE);
	//write_h5_(name, set_name, data, 5, 1);

	return 0;
}
