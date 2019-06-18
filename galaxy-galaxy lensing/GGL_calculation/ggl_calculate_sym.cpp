#include<FQlib.h>
#include<mpi.h>

#define DATA_COL 5
#define MAX_AREA 20
#define CRIT_COEFF 388.2833518
#define MG_BIN_NUM 8

int main(int argc, char ** argv)
{
	/* crit_coeff is the coefficient of the critical density, see the note for detials */

	/* argv[1]: 1for calculate	gamma_t, 2 for calculate the differential surface mass density	*/
	/* argv[2]: the file name of foreground source		*/
	/* argv[3] -- argv[n]: area labels, like 1,3, 4			*/

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	char total_path[300], data_path[300], set_name[30], result_path[300];
	char fore_source[50];
	char log[300];

	int i, j, k;
	int choice, shear_cmd, chi_fit_num=20;

	int area_id[MAX_AREA]{}, area_num, radius_id;
	int total_data_num, total_pair_num;
	int data_num[MAX_AREA], pair_num[MAX_AREA];

	double *data[MAX_AREA];
	double *mgt, *mgx, *mnu1, *mnu2, *crit;
	double *total_result, *total_chi_check;
	double *chi_check[2];
	double gt, gx, gt_sig, gx_sig, coeff;
	double st1, st2;

	shear_cmd= atoi(argv[1]);
	strcpy(fore_source, argv[2]);
	area_num = argc - 3;
	for (i = 3; i < argc; i++)
	{
		area_id[i - 3] = atoi(argv[i]);
	}
	if (0 == rank)
	{
		std::cout << "Foreground: " << fore_source << std::endl;
		std::cout << "Area: ";
		for (i = 0; i < area_num; i++)
		{
			std::cout << area_id[i] << " ";
		}
		std::cout << std::endl;
	}

	sprintf(total_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/");

	total_data_num = 0;
	total_pair_num = 0;
	radius_id = rank;


	// read data from each files
	for (i = 0; i < area_num; i++)
	{	
		sprintf(data_path, "%sresult/%s/w_%d/radius_%d.hdf5", total_path, fore_source, area_id[i], radius_id);
		sprintf(set_name, "/pair_data");
		read_h5_datasize(data_path, set_name, data_num[i]);
		pair_num[i] = data_num[i] / DATA_COL;

		total_data_num += data_num[i];
		total_pair_num += pair_num[i];

		data[i] = new double[data_num[i]]{};
		read_h5(data_path, set_name, data[i]);
	}
	for (i = 0; i < numprocs; i++)
	{
		if (i == rank)
		{
			for (j = 0; j < area_num; j++)
			{
				std::cout << "In area " << area_id[j] << ", " << pair_num[j] << " galaxy pairs are found in " << radius_id << " radius bin." << std::endl;
			}
			std::cout << "Total " << total_pair_num << " galaxy pairs are found" << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	mgt = new double[total_pair_num];
	mgx = new double[total_pair_num];
	mnu1 = new double[total_pair_num];
	mnu2 = new double[total_pair_num];
	//crit = new double[total_pair_num];
	chi_check[0] = new double[2*chi_fit_num];
	chi_check[1] = new double[2*chi_fit_num];

	// stack the data into one array
	for (i = 0; i < area_num; i++)
	{
		k = 0;
		for (j = 0; j < i; j++)
		{
			k += pair_num[j];
		}
		for (j = 0; j < pair_num[i]; j++)
		{
			mgt[k + j] = data[i][j*DATA_COL] * data[i][j*DATA_COL + 4]* CRIT_COEFF;
			mgx[k + j] = data[i][j*DATA_COL + 1] * data[i][j*DATA_COL + 4]* CRIT_COEFF;
			if (shear_cmd == 1)
			{
				mnu1[k + j] = (data[i][j*DATA_COL + 2] + data[i][j*DATA_COL + 3])* data[i][j*DATA_COL + 4]* CRIT_COEFF;
				mnu2[k + j] = (data[i][j*DATA_COL + 2] -  data[i][j*DATA_COL + 3])* data[i][j*DATA_COL + 4]* CRIT_COEFF;
			}
			else
			{
				mnu1[k + j] = data[i][j*DATA_COL + 2] + data[i][j*DATA_COL + 3];
				mnu2[k + j] = data[i][j*DATA_COL + 2] - data[i][j*DATA_COL + 3];
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	
	choice = 0;
	if (total_pair_num > 200000)
	{
		choice = 100000;
	}
	//double bins[MG_BIN_NUM];
	//set_bin(mgt, total_pair_num, bins, MG_BIN_NUM, 1, choice);
	//for (i = 0; i < numprocs; i++)
	//{
	//	if (i == rank)
	//	{
	//		//for (k = 0; k < 100; k++)
	//		//{
	//		//	std::cout << mgt[k] << std::endl;
	//		//}
	//		show_arr(bins, 1, MG_BIN_NUM);
	//	}
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}
	//exit(0);

	MPI_Win win_result, win_chi;
	MPI_Aint size_result, size_chi;

	size_result = 4 * numprocs;
	size_chi = 2 * chi_fit_num*numprocs;

	if (0 == rank)
	{
		MPI_Win_allocate_shared(size_result*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &total_result, &win_result);
		MPI_Win_allocate_shared(size_chi * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &total_chi_check, &win_chi);
	}
	else
	{
		int disp_unit;
		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &total_result, &win_result);
		MPI_Win_shared_query(win_result, 0, &size_result, &disp_unit, &total_result);

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &total_chi_check, &win_chi);
		MPI_Win_shared_query(win_chi, 0, &size_chi, &disp_unit, &total_chi_check);
	}

	double gh_left, gh_right;
	gh_left = -0.5*CRIT_COEFF;
	gh_right = 0.5*CRIT_COEFF;
	if (shear_cmd == 1)
	{
		gh_left = -0.1;
		gh_right = 0.1;
	}
	st1 = clock();
	try
	{
		find_shear(mgt, mnu1, total_pair_num, MG_BIN_NUM, gt, gt_sig, chi_check[0], chi_fit_num, choice, 1000, gh_left, gh_right, 50);
	}
	catch (const char* msg)
	{
		std::cerr << "Rank: " << rank << " g1  " << msg << std::endl;
		gt = 0;
		gt_sig = 0;
	}
	try
	{
		find_shear(mgx, mnu2, total_pair_num, MG_BIN_NUM, gx, gx_sig, chi_check[1], chi_fit_num, choice, 1000, gh_left, gh_right, 50);
	}
	catch (const char* msg)
	{
		std::cerr << "Rank: " << rank <<" g2 "<< msg << std::endl;
		gx = 0;
		gx_sig = 0;
	}
	st2 = clock();

	coeff = 1;
	//if (shear_cmd > 1)
	//{
	//	coeff = CRIT_COEFF;
	//}
	total_result[4 * rank] = gt* coeff;
	total_result[4 * rank + 1] = gt_sig* coeff;
	total_result[4 * rank + 2] = gx* coeff;
	total_result[4 * rank + 3] = gx_sig* coeff;
	for (i = 0; i < chi_fit_num; i++)
	{
		total_chi_check[2*chi_fit_num*rank + i] = chi_check[0][i];
		total_chi_check[2*chi_fit_num*rank + chi_fit_num + i] = chi_check[1][i];
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for (i = 0; i < numprocs; i++)
	{
		if (i == rank)
		{
			std::cout << "Radius: " << radius_id << ".  Time: " << (st2 - st1) / CLOCKS_PER_SEC << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}


	if (0 == rank)
	{
		char result_name[40];
		if (shear_cmd == 1)
		{
			sprintf(result_name, "result_gamma.hdf5");
		}
		else
		{
			sprintf(result_name, "result_crit.hdf5");
		}
		if (area_num == 1)
		{
			sprintf(result_path, "%sresult/%s/w_%d/%s", total_path, fore_source, area_id[0],result_name);
		}
		else
		{
			sprintf(result_path, "%sresult/%s/total/%s", total_path, fore_source, result_name);
		}
		sprintf(set_name, "/result");
		write_h5(result_path, set_name, total_result, numprocs, 4, TRUE);
		sprintf(set_name, "/chisq");
		write_h5(result_path, set_name, total_result, numprocs, 2*chi_fit_num, FALSE);
	}


	for (i = 0; i < area_num; i++)
	{
		delete[] data[i];
	}
	MPI_Win_free(&win_chi);
	MPI_Win_free(&win_result);

	delete[] mgt;
	delete[] mgx;
	delete[] mnu1;
	delete[] mnu2;

	MPI_Finalize();
	return 0;
}