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
	int cut_num = 1;
	double st1, st2, st3, st4, st5, st6;

	int shear_num = 10;
	if (numprocs > shear_num)
	{
		std::cout << "Cpu "<<numprocs << " > shear_num" << shear_num << std::endl;
		exit(0);
	}

	int data_num = 10000000, data_col = 7;
	int size = data_num * data_col;
	char data_path[200], set_name[50], log_inform[250];

	double gh1, gh2, gh1_sig, gh2_sig;
	double m1, m1_sig, m2, m2_sig, c1, c1_sig, c2, c2_sig;

	double *shear = new double[shear_num * 2];
	double *g1 = new double[shear_num];
	double *g2 = new double[shear_num];
	double *result;

	double *data = new double[size];
	double *mg1 = new double[data_num];
	double *mg2 = new double[data_num];
	double *mnu1 = new double[data_num];
	double *mnu2 = new double[data_num];
	//double *pk0 = new double[data_num];


	std::string data_path_s;
	std::stringstream str_1, str_2;

	MPI_Win win_shear;
	MPI_Aint shear_size;
	shear_size = shear_num * 4 * sizeof(double);

	data_path_s = "/mnt/ddnfs/data_users/hkli/simu_test1/parameters/shear.dat";
	read_text(data_path_s, shear, 2 * shear_num);

	if (rank == 0)
	{

		for (i = 0; i < shear_num; i++)
		{
			g1[i] = shear[i];
			g2[i] = shear[i + shear_num];
		}
		
		MPI_Win_allocate_shared(shear_size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &result, &win_shear);
	}
	else
	{
		int disp_unit;
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &result, &win_shear);
		MPI_Win_shared_query(win_shear, 0, &shear_size, &disp_unit, &result);
	}

	sprintf(set_name, "/data");
	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/simu_test1/result/data/data_%d.hdf5", rank);
	read_h5(data_path, set_name, data);

	for (i = 0; i < data_num; i++)
	{
		mg1[i] = data[i*data_col + 2];
		mg2[i] = data[i*data_col + 3];
		mnu1[i] = data[i*data_col + 4] + data[i*data_col + 5];
		mnu2[i] = data[i*data_col + 4] - data[i*data_col + 5];
	}
	st1 = clock();
	/* calculate the g1 & g2 */
	find_shear(mg1, mnu1, data_num, 12, gh1, gh1_sig, 10000,-0.1, 0.1);
	find_shear(mg2, mnu2, data_num, 12, gh2, gh2_sig, 10000, -0.1, 0.1);
	st2 = clock();
	result[rank] = gh1;
	result[rank +shear_num] = gh1_sig;
	result[rank + shear_num*2] = gh2;
	result[rank + shear_num*3] = gh2_sig;
	sprintf(log_inform, "Rank %d. True g1: %8.5f Est: %8.5f (%8.5f). True g2: %8.5f Est: %8.5f (%8.5f). %.2f sec", rank, shear[rank], gh1, gh1_sig, shear[rank + shear_num], gh2, gh2_sig, (st2-st1)/CLOCKS_PER_SEC);
	std::cout << log_inform << std::endl;

	MPI_Barrier(MPI_COMM_WORLD);
	
	// least square to fit the m & c
	if (rank == 0)
	{	
		double coeff[4];
		double *measured_g1 = new double[shear_num];
		double *measured_g2 = new double[shear_num];
		double *measured_g1_sig = new double[shear_num];
		double *measured_g2_sig = new double[shear_num];
			 
		for (i = 0; i < shear_num; i++)
		{
			measured_g1[i] = result[i];
			measured_g1_sig[i] = result[i + shear_num];
			measured_g2[i] = result[i + shear_num * 2];
			measured_g2_sig[i] = result[i + shear_num * 3];
		}
		sprintf(data_path, "/home/hkli/work/shear_result.hdf5");
		sprintf(set_name, "/data");
		write_h5(data_path, set_name, result, 4, shear_num, TRUE);
		//show_arr(measured_g1, 1, shear_num);
		//show_arr(measured_g1_sig, 1, shear_num);
		//show_arr(measured_g2, 1, shear_num);
		//show_arr(measured_g2_sig, 1, shear_num);

		st3 = clock();
		poly_fit_1d(g1, measured_g1, measured_g1_sig, shear_num, coeff, 1);
		st4 = clock();
		show_arr(coeff, 1, 4);
		sprintf(log_inform, "m1: %8.6f (%8.6f), c1: %9.6f (%9.6f). %.2f", coeff[2]-1, coeff[3], coeff[0], coeff[1], (st4-st3)/CLOCKS_PER_SEC);
		std::cout << log_inform << std::endl;

		st5 = clock();
		poly_fit_1d(g2, measured_g2, measured_g2_sig, shear_num, coeff, 1);
		st6 = clock();
		show_arr(coeff, 1, 4);
		sprintf(log_inform, "m2: %6.4f (%6.4f), c2: %9.6f (%9.6f)", coeff[2] - 1, coeff[3], coeff[0], coeff[1], (st6 - st5) / CLOCKS_PER_SEC);
		std::cout << log_inform << std::endl;

		delete[] measured_g1;
		delete[] measured_g2;
		delete[] measured_g1_sig;
		delete[] measured_g2_sig; 
	}
	
	delete[] mnu2;
	delete[] mnu1;
	delete[] mg2;
	delete[] mg1;
	delete[] data;
	delete[] g1;
	delete[] g2;
	delete[] shear;
	MPI_Finalize();
	return 0;
}