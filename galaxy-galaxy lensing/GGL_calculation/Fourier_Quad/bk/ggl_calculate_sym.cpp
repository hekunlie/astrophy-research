#include<FQlib.h>
#include<mpi.h>

#define DATA_COL 6
#define RESULT_ROW 11
#define MAX_AREA 20
#define CRIT_COEFF 174648.0703379
#define MG_BIN_NUM 8

int main(int argc, char ** argv)
{
	/* crit_coeff is the coefficient of the critical density, see the note for detials */

	/* argv[1]: the file name of foreground source		*/
	/* argv[2] -- argv[n]: area labels, like 1,3, 4			*/

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
	int choice, chi_fit_num=20;

	int area_id[MAX_AREA]{}, area_num, radius_id;
	int total_data_num, total_pair_num;
	int data_num[MAX_AREA], pair_num[MAX_AREA];

	int file_mask[1];

	double *data[MAX_AREA];
	double *mgt, *mgx, *mnu1, *mnu2, *mnu1_crit, *mnu2_crit, *crit, *dist_radius;
	double *total_result, *total_chi_check;
	double *chi_check[2];
	double gt, gx, gt_sig, gx_sig, coeff;
	double sigma_t, sigma_x, sigma_t_sig, sigma_x_sig;

	double radius_mean, radius_mean_temp;

	double gh_left, gh_right;

	double st1, st2;


	strcpy(fore_source, argv[1]);
	area_num = argc - 2;
	for (i = 2; i < argc; i++)
	{
		area_id[i - 2] = atoi(argv[i]);
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

	//sprintf(total_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/");
	sprintf(total_path, "/mnt/perc/hklee/CFHT/gg_lensing/");

	total_data_num = 0;
	total_pair_num = 0;
	radius_id = rank;


	// read data from each files
	initialize_arr(pair_num, MAX_AREA, 0);
	initialize_arr(data_num, MAX_AREA, 0);
	for (i = 0; i < area_num; i++)
	{	
		sprintf(data_path, "%sresult/%s/fourier/w_%d/radius_%d.hdf5", total_path, fore_source, area_id[i], radius_id);
		sprintf(set_name, "/mask");
		read_h5(data_path, set_name, file_mask);
		if (rank == 0)
		{
			std::cout << file_mask[0] << std::endl;
		}
		if (file_mask[0] == 1)
		{
			sprintf(set_name, "/pair_data");
			read_h5_datasize(data_path, set_name, data_num[i]);
			pair_num[i] = data_num[i] / DATA_COL;

			total_data_num += data_num[i];
			total_pair_num += pair_num[i];

			data[i] = new double[data_num[i]]{};
			read_h5(data_path, set_name, data[i]);
		}
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

	MPI_Win win_result, win_chi;
	MPI_Aint size_result, size_chi;

	size_result = RESULT_ROW * numprocs;
	size_chi = 2 * chi_fit_num*numprocs;

	if (0 == rank)
	{
		MPI_Win_allocate_shared(size_result * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &total_result, &win_result);
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

	mgt = new double[total_pair_num];
	mgx = new double[total_pair_num];
	mnu1 = new double[total_pair_num];
	mnu2 = new double[total_pair_num];
	mnu1_crit = new double[total_pair_num];
	mnu2_crit = new double[total_pair_num];
	dist_radius = new double[total_pair_num];
	//crit = new double[total_pair_num];
	chi_check[0] = new double[2 * chi_fit_num];
	chi_check[1] = new double[2 * chi_fit_num];

	if (total_pair_num > 10)
	{   
		
		// stack the data into one array
		for (i = 0; i < area_num; i++)
		{
			// if the file contains data
			k = 0;
			for (j = 0; j < i; j++)
			{
				k += pair_num[j];
			}
			for (j = 0; j < pair_num[i]; j++)
			{
				mgt[k + j] = data[i][j*DATA_COL] * data[i][j*DATA_COL + 4] * CRIT_COEFF;
				mgx[k + j] = data[i][j*DATA_COL + 1] * data[i][j*DATA_COL + 4] * CRIT_COEFF;

				mnu1_crit[k + j] = (data[i][j*DATA_COL + 2] + data[i][j*DATA_COL + 3])* data[i][j*DATA_COL + 4] * CRIT_COEFF;
				mnu2_crit[k + j] = (data[i][j*DATA_COL + 2] - data[i][j*DATA_COL + 3])* data[i][j*DATA_COL + 4] * CRIT_COEFF;

				mnu1[k + j] = data[i][j*DATA_COL + 2] + data[i][j*DATA_COL + 3];
				mnu2[k + j] = data[i][j*DATA_COL + 2] - data[i][j*DATA_COL + 3];

				dist_radius[k + j] = data[i][j*DATA_COL + 5];
			}

		}

		// the mean radius from the center of each raidus bin
		radius_mean = 0;
		radius_mean_temp = 0;
		for (i = 0; i < total_pair_num; i++)
		{
			radius_mean_temp += dist_radius[i];

			if (radius_mean_temp > 100000)
			{
				radius_mean += radius_mean_temp;
				radius_mean_temp = 0;
			}
		}
		radius_mean += radius_mean_temp;
		radius_mean = radius_mean / total_pair_num;
		std::cout << rank<<" "<<radius_mean <<"Mpc"<< std::endl;

		choice = 0;
		if (total_pair_num > 200000)
		{
			choice = 100000;
		}
		gh_left = -0.5*CRIT_COEFF;
		gh_right = 0.5*CRIT_COEFF;

		st1 = clock();
		try
		{
			find_shear(mgt, mnu1, total_pair_num, MG_BIN_NUM, sigma_t, sigma_t_sig, chi_check[0], chi_fit_num, choice, 1000, gh_left, gh_right, 50);
		}
		catch (const char* msg)
		{
			std::cerr << "Rank: " << rank << " Sigma_t  " << msg << std::endl;
			sigma_t = 0;
			sigma_t_sig = 0;
		}
		try
		{
			find_shear(mgx, mnu2, total_pair_num, MG_BIN_NUM, sigma_x, sigma_x_sig, chi_check[1], chi_fit_num, choice, 1000, gh_left, gh_right, 50);
		}
		catch (const char* msg)
		{
			std::cerr << "Rank: " << rank << " Sigma_x " << msg << std::endl;
			sigma_x = 0;
			sigma_x_sig = 0;
		}

		gh_left = -0.1;
		gh_right = 0.1;

		/*try
		{
			find_shear(mgt, mnu1_crit, total_pair_num, MG_BIN_NUM, gt, gt_sig, chi_check[0], chi_fit_num, choice, 1000, gh_left, gh_right, 50);
		}
		catch (const char* msg)
		{
			std::cerr << "Rank: " << rank << " g_t  " << msg << std::endl;
			gt = 0;
			gt_sig = 0;
		}
		try
		{
			find_shear(mgx, mnu2_crit, total_pair_num, MG_BIN_NUM, gx, gx_sig, chi_check[1], chi_fit_num, choice, 1000, gh_left, gh_right, 50);
		}
		catch (const char* msg)
		{
			std::cerr << "Rank: " << rank << " g_x " << msg << std::endl;
			gx = 0;
			gx_sig = 0;
		}
		*/
		st2 = clock();

		total_result[rank] = gt;
		total_result[1 * numprocs + rank] = gt_sig;
		total_result[2 * numprocs + rank] = gx;
		total_result[3 * numprocs + rank] = gx_sig;

		total_result[4 * numprocs + rank] = sigma_t;
		total_result[5 * numprocs + rank] = sigma_t_sig;
		total_result[6 * numprocs + rank] = sigma_x;
		total_result[7 * numprocs + rank] = sigma_x_sig;

		total_result[8 * numprocs + rank] = sigma_t * radius_mean;
		total_result[9 * numprocs + rank] = sigma_t_sig * radius_mean;

		total_result[10 * numprocs + rank] = radius_mean;

		for (i = 0; i < chi_fit_num; i++)
		{
			total_chi_check[2 * chi_fit_num*rank + i] = chi_check[0][i];
			total_chi_check[2 * chi_fit_num*rank + chi_fit_num + i] = chi_check[1][i];
		}

		std::cout << "Radius: " << radius_id << ".  Time: " << (st2 - st1) / CLOCKS_PER_SEC << std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (numprocs - 1 == rank)
	{
		if (area_num == 1)
		{
			sprintf(result_path, "%sresult/%s/fourier/result/fourier_%s_result_w_%d.hdf5", total_path, fore_source, fore_source, area_id[0]);
		}
		else
		{
			sprintf(result_path, "%sresult/%s/fourier/result/fourier_%s_result_total.hdf5", total_path, fore_source, fore_source);
		}
		sprintf(set_name, "/data");
		write_h5(result_path, set_name, total_result, RESULT_ROW, numprocs, TRUE);

		char times[50];
		get_time(times, 50);
		std::cout << times << std::endl;

	}

	for (i = 0; i < area_num; i++)
	{
		if (pair_num[i] > 0)
		{
			delete[] data[i];
		}
	}
	MPI_Win_free(&win_chi);
	MPI_Win_free(&win_result);

	
	delete[] mgt;
	delete[] mgx;
	delete[] mnu1;
	delete[] mnu2;
	delete[] mnu1_crit;
	delete[] mnu2_crit;
	delete[] dist_radius;

	MPI_Finalize();
	return 0;
}