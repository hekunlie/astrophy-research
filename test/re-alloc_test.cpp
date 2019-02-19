#include<FQlib.h>
#include"mpi.h"

int main(int argc, char *argv[])
{
	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	int i, j, k;
	double *mask;
	double temp = 0, temp_=0;
	double *buffer_ptr;
	int *num_ptr;
	int *err_ptr;
	int arr_len=500000;
	int grid_num = 100;
	int disp_unit, disp_unit_num, disp_unit_err;
	MPI_Win win, win_num, win_err;
	MPI_Aint buffer_size, size_num, size_err;

	char data_path[100], set_name[30];
	char informs[200];

	sprintf(data_path, "/home/hkli/work/cpp/test/test.hdf5");
	sprintf(set_name, "/data");
	double *writedat = new double[arr_len] {};
	initialize_arr(writedat, arr_len, 1);
	write_h5(data_path, set_name, writedat, arr_len, 1);

	long *long_test_w = new long[arr_len] {};
	long *long_test_r = new long[arr_len] {};
	long long_temp = 0;
	for (i = 0; i < arr_len; i++)
	{
		long_test_w[i] = 1;
	}
	sprintf(data_path, "/home/hkli/work/cpp/test/long_test.hdf5");
	sprintf(set_name, "/data");
	write_h5(data_path, set_name, long_test_w, arr_len, 1);

	for (i = 0; i < 10; i++)
	{
		mask = new double[arr_len]{};
		temp = 0;
		temp_ = 0;
		long_temp = 0;
		if (rank == 0)
		{
			buffer_size = arr_len * sizeof(double);
			MPI_Win_allocate_shared(buffer_size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &buffer_ptr, &win);
			size_num = grid_num * sizeof(int);
			MPI_Win_allocate_shared(size_num, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &num_ptr, &win_num);
			size_err = sizeof(int);
			MPI_Win_allocate_shared(size_err, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &err_ptr, &win_err);
		}
		else
		{
			MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &buffer_ptr, &win);
			MPI_Win_shared_query(win, 0, &buffer_size, &disp_unit, &buffer_ptr);

			MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &num_ptr, &win_num);
			MPI_Win_shared_query(win_num, 0, &size_num, &disp_unit_num, &num_ptr);

			MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &err_ptr, &win_err);
			MPI_Win_shared_query(win_err, 0, &size_err, &disp_unit_err, &err_ptr);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		// check: read the hdf5 by the buffer itself
		sprintf(data_path, "/home/hkli/work/cpp/test/test.hdf5");
		read_h5(data_path, set_name, buffer_ptr);
		sprintf(data_path, "/home/hkli/work/cpp/test/long_test.hdf5");
		read_h5(data_path, set_name, long_test_r);
		for (j = 0; j < arr_len; j++)
		{
			temp += buffer_ptr[j];
			long_temp += long_test_r[i];
			if (long_test_r[i] != 1)
			{
				std::cout << "Wrong " << long_test_r << std::endl;
			}
		}
		for (j = 0; j < numprocs; j++)
		{
			if (j == rank)
			{
				sprintf(informs, "%d reads hdf5 & accumulates: %g (%d), Long: %d", rank, temp, arr_len, long_temp);
				std::cout << informs<< std::endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		temp = 0;
		
		
		if (rank == 0)
		{
			for (j = 0; j < arr_len; j++)
			{
				buffer_ptr[j] = 1;				
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (j = 0; j < arr_len; j++)
		{
			temp +=buffer_ptr[j];
			mask[j] = 2;
		}
		for (j = 0; j < numprocs; j++)
		{
			if (j == rank)
			{
				std::cout << rank << ", " << temp << std::endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

		temp = 0;
		for (j = 0; j < arr_len; j++)
		{
			temp += mask[j];
		}
		for (j = 0; j < numprocs; j++)
		{
			if (j == rank)
			{
				std::cout << rank << ", " << temp << std::endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		if (rank == 0)
		{
			std::cout << std::endl;
		}
		MPI_Win_free(&win);
		MPI_Win_free(&win_num);
		MPI_Win_free(&win_err);
		delete[] mask;
		
	}

	delete[] writedat;
	delete[] long_test_w;
	delete[] long_test_r;
	return 0;

}