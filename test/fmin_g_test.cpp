#include<FQlib.h>
#include<mpi.h>

int main(int argc, char**argv)
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	int i, j, k;
	char data_path[200], set_name[50], inform[250];
	std::string data_path_s;
	data_path_s = "/mnt/ddnfs/data_users/hkli/selection_bias_64/parameters/shear.dat";
	int data_row = 10000000, data_col = 7;
	int shear_num = 10;
	
	double st1, st2, st3;

	double gh1, gh1_sig, gh2, gh2_sig;
	
	double *data = new double[data_row*data_col]{};
	double *shear = new double[shear_num * 2]{};
	double *mg1 = new double[data_row] {};
	double *mg2 = new double[data_row] {};
	double *mnu1 = new double[data_row] {};
	double *mnu2 = new double[data_row] {};

	read_text(data_path_s, shear, shear_num * 2);

	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/selection_bias_64_bright/result/data/data_%d.hdf5", rank);
	sprintf(set_name, "/data");
	read_h5(data_path, set_name, data);

	for (i = 0; i < data_row; i++)
	{
		mg1[i] = data[i*data_col + 2];
		mg2[i] = data[i*data_col + 3];
		mnu1[i] = data[i*data_col + 4] + data[i*data_col + 5];
		mnu2[i] = data[i*data_col + 4] - data[i*data_col + 5];
	}
	sprintf(inform, "%d My g1: %7.5f, g2: %7.5f", rank,  shear[rank], shear[shear_num + rank]);
	for (i = 0; i < numprocs; i++)
	{
		if (i == rank)
		{
			std::cout << inform << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	st1 = clock();
	find_shear(mg1, mnu1, data_row, 12, gh1, gh1_sig);
	st2 = clock();
	find_shear(mg2, mnu2, data_row, 12, gh2, gh2_sig);
	st3 = clock();
	sprintf(inform, "%d True g1: %7.5f, True g2: %7.5f. \nEst g1: %7.5f (%7.5f). Est g2: %7.5f (%7.5f). \nTime: %.2f. %.2f", rank, shear[rank], shear[rank + shear_num], gh1, gh1_sig,  gh2, gh2_sig, (st2-st1)/CLOCKS_PER_SEC, (st3-st2)/CLOCKS_PER_SEC);
	for (i = 0; i < numprocs; i++)
	{
		if(i==rank)
		{
			std::cout << inform << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	

	delete[] mg2;
	delete[] mg1;
	delete[] mnu1;
	delete[] mnu2;
	delete[] shear;
	delete[] data;
	MPI_Finalize();
	return 0;
}