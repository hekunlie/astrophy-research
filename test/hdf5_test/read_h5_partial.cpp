#include<FQlib.h>
#include<mpi.h>

int main(int argc, char**argv)
{
	char filename[50];
	char setname[20];
	sprintf(filename, "test.hdf5");
	sprintf(setname, "/data");
	double *data, *data_partial;
	int num, num_partial, num_s, num_e;

	read_h5_datasize(filename, setname, num);
	
	data = new double[num];
	read_h5(filename, setname, data);

	num_s = atoi(argv[1]);
	num_e = atoi(argv[2]);
	num_partial = num_e - num_s;
	data_partial = new double[num_partial];
	
	std::cout << "The whole data: " << std::endl;
	show_arr(data, 1, num);
	std::cout << "Select: " << num_s << " to " << num_e << std::endl;

	/* get the file dataspace. */
	hsize_t offset[1], count[1], dimsm[2], count_out[2], offset_out[2];
	hid_t file_id;
	herr_t status;
	hid_t dataset_id, dataspace_id, memspace_id;

	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset_id = H5Dopen(file_id, setname, H5P_DEFAULT);
	dataspace_id = H5Dget_space(dataset_id); /* dataspace  identifier */
	std::cout << "Get dataspace_id" << std::endl;

	/* Define hyperslab in the dataset.	*/
	offset[0] = num_s;
	count[0] = num_partial;
	status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET,offset, NULL, count, NULL);
	std::cout << "Define hyperslab in the dataset" << std::endl;

	/* Define memory dataspace.	*/
	dimsm[0] = 1;
	dimsm[1] = num_s;
	memspace_id = H5Screate_simple(2, dimsm, NULL);

	/*  Define memory hyperslab. 	*/
	offset_out[0] = 0;
	offset_out[1] = 0;
	count_out[0] = num_s;
	count_out[1] = num_s;
	status = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET,	offset_out, NULL, count_out, NULL);
	std::cout << "Define memory hyperslab" << std::endl;

	H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, data_partial);

	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Fclose(file_id);

	char inform[100];
	sprintf(inform, "Select [%d, %d)", num_s, num_e);
	std::cout << inform << std::endl;
	show_arr(data_partial, 1, num_partial);

	return 0;
}