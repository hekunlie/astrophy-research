#include<FQlib.h>
#include<mpi.h>

#define data_col 6
#define max_area 20
#define crit_coeff 388.2833518
#define mg_bin_num 12

/* argv[1]: the file name of foreground source		*/
/* argv[2]: radius bin label										*/
/* argv[3] -- argv[n]: area labels, like 1,3, 4			*/

int main(int argc, char ** argv)
{
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
	double st1, st2;
	int area_id[max_area]{}, area_num, radius_id;
	int total_data_num, total_pair_num;
	int data_num[max_area], pair_num[max_area];


	strcpy(fore_source, argv[1]);
	radius_id = atoi(argv[2]);
	area_num = argc - 3;
	for (i = 3; i < argc; i++)
	{
		area_id[i - 3] = atoi(argv[i]);
	}
	if (0 == rank)
	{
		std::cout << "Foreground: " << fore_source << std::endl;
		std::cout << "Radius bin: " << radius_id << std::endl;
		std::cout << "Area: ";
		for (i = 0; i < area_num; i++)
		{
			std::cout << area_id[i]<<" ";
		}
		std::cout << std::endl;
	}

	sprintf(total_path, "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/");
	total_data_num = 0;
	total_pair_num = 0;
	double *data[max_area];
	double *mgt, *mgx, *mnu1, *mnu2, *crit, *mg1;
	int choice;
	double gh_left, gh_right, gstep, *gh;
	int gh_num, gh_s, gh_e;
	double *mg_bin;
	double*total_chisq, chisq_temp;
	double gt, gx, gt_sig, gx_sig;

	gh_left = -130 + 10 * rank;
	gh_right = fabs(gh_left);
	
	gstep = 0.5;
	gh_num = int(gh_right * 2/ gstep);

	gh_s = gh_num / numprocs * rank;
	gh_e = gh_num / numprocs * (rank + 1);
	if (rank == numprocs - 1)
	{
		gh_e = gh_num / numprocs * (rank + 1) + gh_num % numprocs;
	}

	gh = new double[gh_num] {};
	for (i = 0; i < gh_num; i++)
	{
		gh[i] = gh_left + i* gstep;
	}

	for (i = 0; i < area_num; i++)
	{
		sprintf(data_path, "%sresult/%s/w_%d/%d.hdf5", total_path, fore_source, area_id[i], radius_id);
		sprintf(set_name, "/data_0");//////
		read_h5_datasize(data_path, set_name, data_num[i]);
		pair_num[i] = data_num[i] / data_col;

		total_data_num += data_num[i];
		total_pair_num += pair_num[i];

		data[i] = new double[data_num[i]] {};
		read_h5(data_path, set_name, data[i]);
	}

	mgt = new double[total_pair_num];
	mgx = new double[total_pair_num];
	mnu1 = new double[total_pair_num];
	mnu2 = new double[total_pair_num];
	mg1 = new double[total_pair_num];
	
	for (i = 0; i < area_num; i++)
	{
		k = 0;
		for (j = 0; j < i; j++)
		{
			k += pair_num[j];
		}
		for (j = 0; j < pair_num[i]; j++)
		{
			mgt[k+ j] = data[i][j*data_col] * data[i][j*data_col + 4]* crit_coeff;
			mgx[k + j] = data[i][j*data_col + 1] * data[i][j*data_col + 4]* crit_coeff;
			mnu1[k + j] = data[i][j*data_col + 2] + data[i][j*data_col + 3];
			mnu2[k+ j] = data[i][j*data_col + 2] - data[i][j*data_col + 3];
			mg1[k + j] = data[i][j*data_col + 5];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for (i = 0; i < numprocs; i++)
	{
		if (i == rank)
		{	
			std::cout << "gh_num: " << gh_num << std::endl;
			std::cout << "gh: " << gh_s << " ~ " << gh_e << std::endl;
			std::cout << "Radius: " << radius_id << ".  Pair num: ";
			for (j = 0; j < area_num; j++)
			{
				std::cout << pair_num[j] << " ";
			}
			std::cout << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}


	MPI_Win win_bin, win_chisq;
	MPI_Aint bin_size, chisq_size;
	bin_size = mg_bin_num+1;
	chisq_size = 2*gh_num;

	if (rank == 0)
	{
		MPI_Win_allocate_shared(bin_size * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD,
			&mg_bin, &win_bin);
		MPI_Win_allocate_shared(chisq_size * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, 
			&total_chisq, &win_chisq);
	}
	else
	{
		int disp_unit, disp_unit_c;
		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &mg_bin, &win_bin);
		MPI_Win_shared_query(win_bin, 0, &bin_size, &disp_unit, &mg_bin);

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &total_chisq, &win_chisq);
		MPI_Win_shared_query(win_chisq, 0, &chisq_size, &disp_unit_c, &total_chisq);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// set up bins for PDF-SYM
	if (0 == rank)
	{
		if (total_pair_num < 200000)
		{
			choice = 0;
		}
		else
		{
			choice = 100000;
		}
		set_bin(mgt, total_pair_num, mg_bin, mg_bin_num, 10000, choice);
	}
	MPI_Barrier(MPI_COMM_WORLD);


	st1 = clock();
	for(i = gh_s; i < gh_e; i++)
	{
		try 
		{
			chisq_Gbin_1d(mgt, mnu1, total_pair_num, mg_bin, mg_bin_num, gh[i], chisq_temp);
			total_chisq[i] = chisq_temp;
		}
		catch (const char* msg)
		{
			std::cerr << "Rank: " << rank <<". g_gusee: "<<i<< ". Tangential. " << msg << std::endl;
			exit(0);
		}
		try
		{
			chisq_Gbin_1d(mgx, mnu2, total_pair_num, mg_bin, mg_bin_num, gh[i], chisq_temp);
			total_chisq[gh_num + i] = chisq_temp;
		}
		catch (const char* msg)
		{
			std::cerr << "Rank: " << rank << ". g_gusee: " << i<<". Cross. " << msg << std::endl;
			exit(0);
		}		
	}
	st2 = clock();

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
		double *chisq_fit = new double[gh_num];
		double signal[4];

		sprintf(result_path, "%sresult/%s/%d.hdf5", total_path, fore_source, radius_id);
		sprintf(set_name, "/chisq");
		write_h5(result_path, set_name, total_chisq, 2, gh_num, TRUE);

		for (j = 0; j < 2; j++)
		{
			for (i = 0; i < gh_num; i++)
			{
				chisq_fit[i] = total_chisq[i + j * gh_num];
			}
			fit_shear(gh, chisq_fit, gh_num, signal[j * 2], signal[j * 2 + 1], 60);
		}
		sprintf(set_name, "/signal");
		write_h5(result_path, set_name, signal, 4, 1, FALSE);

		delete[] chisq_fit;
	}
	MPI_Win_free(&win_bin);
	MPI_Win_free(&win_chisq);
	for (i = 0; i < area_num; i++)
	{
		delete[] data[i];
	}

	delete[] gh;
	delete[] mgt;
	delete[] mgx;
	delete[] mnu1;
	delete[] mnu2;
	delete[] mg1;

	MPI_Finalize();
	return 0;
}