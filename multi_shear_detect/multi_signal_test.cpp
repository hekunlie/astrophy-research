#include<FQlib.h>
#include<mpi.h>

#define signal_num 2

int main(int argc, char **argv)
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);
	
	int i, j, k;
	int seed = 1223;
	char data_path[200], set_name[20];

	// the data number of each data set
	int bin_num = atoi(argv[1]);
	int gh_num = atoi(argv[2]);
	int data_num[signal_num]{ 6*10000000, 6*10000000 };	
	double sigma[signal_num]{ 0.05, 0.05 };
	double signal[signal_num]{ -0.05, 0.05 };
	int total_data_num = 0;

	for (i = 0; i < signal_num; i++)
	{
		total_data_num += data_num[i];
	}
	double *data = new double[total_data_num]{};
	double *temp = new double[total_data_num] {};
	double *gh = new double[gh_num];
	for (i = 0; i < gh_num; i++)
	{
		gh[i] = -0.1 + 0.2 / gh_num * i;
	}

	// chi square
	double *chisq, chisq_i;
	int *num_in_bin = new int[bin_num] {};
	double *data_bin = new double[bin_num+1];

	MPI_Win win_chisq;
	MPI_Aint chisq_size;

	if (0 == rank)
	{
		MPI_Win_allocate_shared(gh_num * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chisq, &win_chisq);
	}
	else
	{
		int disp_unit;
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &chisq, &win_chisq);
		MPI_Win_shared_query(win_chisq, 0, &chisq_size, &disp_unit, &chisq);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// the Gaussian data
	gsl_initialize(seed);
	if (0 == rank)
	{		
		for (i = 0; i < signal_num; i++)
		{
			int count = 0;
			for (j = 0; j < i; j++)
			{
				count += data_num[j];
			}
			std::cout << count << std::endl;
			for (j = 0; j< data_num[i]; j++)
			{
				rand_gauss(sigma[i], signal[i], data[j+count]);
			}

			sprintf(set_name, "/data");
			sprintf(data_path, "/mnt/perc/hklee/CFHT/multi_shear/data.hdf5");
			write_h5(data_path, set_name, data, 1, total_data_num, TRUE);
		}		
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(data, total_data_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// set bin
	set_bin(data, total_data_num, data_bin, bin_num, 1000);

	int range_s, range_e, step;
	step = gh_num / numprocs;
	range_s = step * rank;
	if (rank < numprocs - 1)
	{
		range_e = step * (rank + 1);
	}
	else
	{
		range_e = step * (rank + 1) + gh_num % numprocs;
	}
	// calculate the chi square
	for (i = range_s; i < range_e; i++)
	{
		for (j = 0; j < total_data_num; j++)
		{
			temp[j] = data[j] - gh[i];
		}
		initialize_arr(num_in_bin, bin_num, 0);
		histogram(temp, data_bin, num_in_bin, total_data_num, bin_num);
		cal_chisq_1d(num_in_bin, bin_num, chisq[i]);
	}

	//for (i = 0; i < numprocs; i++)
	//{
	//	if (i == rank)
	//	{
	//		show_arr(data_bin, 1, bin_num + 1);
	//		std::cout << std::endl;
	//		std::cout << std::endl;

	//		std::cout << rank << ", " << std::endl;
	//		for (j = 0; j < 20; j++)
	//		{
	//			std::cout << data[j * 5000] << " ";
	//		}
	//		std::cout << std::endl;

	//	}
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}

	if (0 == rank)
	{
		double est_g, est_g_sig;
		sprintf(set_name, "/chisq");
		sprintf(data_path, "/mnt/perc/hklee/CFHT/multi_shear/chisq.hdf5");
		write_h5(data_path, set_name, chisq, 1, gh_num, TRUE);
		sprintf(set_name, "/shear");
		sprintf(data_path, "/mnt/perc/hklee/CFHT/multi_shear/chisq.hdf5");
		write_h5(data_path, set_name, gh, 1, gh_num, FALSE);
		fit_shear(gh, chisq, gh_num, est_g, est_g_sig, 50);
		std::cout << est_g << " " << est_g_sig << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	

	delete[] data;
	delete[] data_bin;
	MPI_Win_free(&win_chisq);
	gsl_free();
	MPI_Finalize();
	return 0;
}