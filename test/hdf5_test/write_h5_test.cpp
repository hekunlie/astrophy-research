#include <FQlib.h>

void write_h5_test(const char *filename, const char *set_name, const double*data, const int row, const int column, const bool trunc)
{
	int i, count, s_count;
	int *slash;
	char *name, *new_name;
	// the length of the set_name
	for (count = 0; ; count++)
	{
		if (set_name[count] == '\0')
		{
			break;
		}
	}
	slash = new int[count] {};
	name = new char[count + 1];
	new_name = new char[count + 1];
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
	// the set_name is something like /a/b/c
	// the dataset will be under c,while /a/b must be created before
	// the "c" is in the "new_names",  "/a/b" in "names"
	if (s_count > 1)
	{
		// the case like /a/b, /a/b/c...
		for (i = 0; i < slash[s_count - 1]; i++)
		{
			name[i] = set_name[i];
		}
		name[slash[s_count - 1]] = '\0';
	}
	else
	{
		// the case like /a
		name[0] = '/';
		name[1] = '\0';
	}

	for (i = slash[s_count - 1] + 1; i < slash[s_count]; i++)
	{
		new_name[i - slash[s_count - 1] - 1] = set_name[i];
	}
	new_name[slash[s_count] - slash[s_count - 1] - 1] = '\0';
	//std::cout << set_name << std::endl;
	//std::cout << name << std::endl;
	//std::cout << new_name << std::endl;
	//show_arr(slash, 1, count);
	//std::cout << s_count << std::endl;
	// try to create /a/b
	create_h5_group(filename, name, trunc);

	hid_t file_id, group_id, dataset_id, dataspace_id;
	herr_t status;
    unsigned rank;

    file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	group_id = H5Gopen1(file_id, name);

    if(row == 1 or column == 1)
    {
        hsize_t dims[1];
	    rank = 1;
        dims[0] = column*row;
        dataspace_id = H5Screate_simple(rank, dims, NULL);
    }
    else
	{
        hsize_t dims[2];
	    rank = 2;
	    dims[0] = row;
	    dims[1] = column;
        dataspace_id = H5Screate_simple(rank, dims, NULL);
    }

	dataset_id = H5Dcreate(group_id, new_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Gclose(group_id);
	status = H5Fclose(file_id);

	delete[] slash;
	delete[] name;
	delete[] new_name;
}

int main(int argc, char**argv)
{
    int num = 100000, i, row, col,data_size;
    double *data_1 = new double[num];

    for(i=0;i<num;i++)
    {
        data_1[i] = i;
    }
    row = num;
    col = 1;

    char data_path[100], set_name[50];

    sprintf(data_path,"test.hdf5");
    sprintf(set_name,"/data");
    write_h5(data_path, set_name, data_1, row, col, true);


    read_h5_datasize(data_path, set_name, data_size);
    double *data_2 = new double[data_size];
    read_h5(data_path, set_name, data_2);

    std::cout<<"Read "<<data_size<<" elements"<<std::endl;

    double diff=0;
    for(i=0;i<data_size;i++)
    {
        diff += data_1[i] - data_2[i];
    }
    std::cout<<diff<<std::endl;

    return 0;
}