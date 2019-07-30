#include<FQlib.h>
#include<mpi.h>

#define DATA_COL 5
#define MAX_AREA 20
#define MG_BIN_NUM 12
#define ALLOWED_SIZE 300000000

//#define SMALL_CATA

#ifdef SMALL_CATA
#define MY_INT int
#else
#define MY_INT long
#endif 

int main(int argc, char ** argv)
{
	/* crit_coeff is the coefficient of the critical density, see the note for detials */
	
	/* argv[1]: the file name of foreground source												*/
	/* argv[2]: radius bin label																				*/
	/* argv[3] -- argv[n]: area labels, like 1,3, 4														*/

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	char total_path[300], data_path[300], set_name[30], result_path[300];
	char fore_source[50];
	char logs[300];

	int i, j, k;
	double st1, st2;

	int area_id[MAX_AREA]{}, area_num, radius_id;

	MY_INT m, n, my_gal_s, my_gal_e;
	MY_INT total_data_num, total_pair_num, my_num, temp_num;
	MY_INT data_num[MAX_AREA]{}, pair_num[MAX_AREA]{};
	MY_INT *pair_count;
	MY_INT *total_chisq, *my_chisq;
	MY_INT data_count[MG_BIN_NUM];

	double *data[MAX_AREA];
	double *All_mgt, *All_mgx, *All_mnut, *All_mnux;
	double *my_mgt, *my_mgx, *my_mnut, *my_mnux, *mg;

	int choice;
	double mg1_bin[MG_BIN_NUM+1]{}, mg2_bin[MG_BIN_NUM+1]{};
	double *mg_bin;

	int gh_num;
	double gh_s, gh_e, dg, *gh;
	double gt, gx, gt_sig, gx_sig;
	double mg_max, mg_min;



	// the foreground name
	strcpy(fore_source, argv[1]);
	// the radius label
	radius_id = atoi(argv[2]);
	area_num = argc - 3;
	for (i = 3; i < argc; i++)
	{
		area_id[i - 3] = atoi(argv[i]);
	}

	sprintf(total_path, "/mnt/perc/hklee/CFHT/gg_lensing/result/%s/fourier_cata_new/", fore_source);

	if (0 == rank)
	{
		int radius_num;
		double *radius_bin;
		sprintf(data_path, "%sradius_bin.hdf5", total_path);
		sprintf(set_name, "/radius_bin");
		read_h5_datasize(data_path, set_name, radius_num);
		radius_bin = new double[radius_num];
		read_h5(data_path, set_name, radius_bin);

		sprintf(logs, "Foreground: %s.\nRaidus_bin: [%.5f, %.5f].\nArea: ", fore_source, radius_bin[radius_id], radius_bin[radius_id + 1]);
		std::cout << logs;
		for (i = 0; i < area_num; i++)
		{
			std::cout << area_id[i]<<" ";
		}
		std::cout << std::endl;
	}
	
	gh_s = -150 + 10 * radius_id;
	gh_e = fabs(gh_s);

	dg = 0.2;

	gh_num = int((gh_e - gh_s)/ dg) + 1;
	gh = new double[gh_num] {};
	for (i = 0; i < gh_num; i++)
	{
		gh[i] = gh_s + i * dg;
	}
	if (rank == 0)
	{
		//std::cout << "The gh: ";
		//show_arr(gh, 1, gh_num);
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// rank 0 collects the data in each area
	// then scatters them to each rank
	total_data_num = 0;
	total_pair_num = 0;

	// collect the data from each file
	for (i = 0; i < area_num; i++)
	{
		sprintf(data_path, "%sw_%d/radius_%d.hdf5", total_path, area_id[i], radius_id);
		sprintf(set_name, "/pair_data");
		read_h5_datasize(data_path, set_name, data_num[i]);
		// each row contains G_t, G_x, N, U, Crit_ij
		pair_num[i] = data_num[i] / DATA_COL;

		total_data_num += data_num[i];
		total_pair_num += pair_num[i];
		if (rank == 0)
		{	
			data[i] = new double[data_num[i]]{};
			read_h5(data_path, set_name, data[i]);
			std::cout << "Read " << pair_num[i] << " pairs in " << i << "'th area" << std::endl;
			if (i == area_num - 1)
			{
				std::cout <<"Total pair: "<< total_pair_num << std::endl;
			}
		}
	}

	if (rank == 0)
	{
		All_mgt = new double[total_pair_num];
		All_mgx = new double[total_pair_num];
		All_mnut = new double[total_pair_num];
		All_mnux = new double[total_pair_num];

		// stack each column from each file
		for (i = 0; i < area_num; i++)
		{
			k = 0;
			for (j = 0; j < i; j++)
			{
				k += pair_num[j];
			}
			for (j = 0; j < pair_num[i]; j++)
			{
				All_mgt[k + j] = data[i][j*DATA_COL] * data[i][j*DATA_COL + 4];
				All_mgx[k + j] = data[i][j*DATA_COL + 1] * data[i][j*DATA_COL + 4];

				All_mnut[k + j] = data[i][j*DATA_COL + 2] + data[i][j*DATA_COL + 3];
				All_mnux[k + j] = data[i][j*DATA_COL + 2] - data[i][j*DATA_COL + 3];
			}
		}
		// set up the G bins
		choice = total_pair_num;
		if (total_pair_num > 500000)
		{
			choice = 200000;
		}
		set_bin(All_mgt, total_pair_num, mg1_bin, MG_BIN_NUM, 1000, choice);
		set_bin(All_mgx, total_pair_num, mg2_bin, MG_BIN_NUM, 1000, choice);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// distribute the source
	// it's safe for LONG/INT
	m = total_pair_num / numprocs;
	n = total_pair_num % numprocs;

	my_gal_s = m * rank;
	my_gal_e = m * (rank + 1);
	if (rank == numprocs - 1)
	{
		my_gal_e += n;
	}

	my_num = my_gal_e - my_gal_s;

	my_chisq = new MY_INT[2 * gh_num*MG_BIN_NUM]{};

	MPI_Win win_chisq, win_mg_bin, win_pair_count;
	MPI_Aint  size_chisq, size_mg_bin, size_pair_count;

	size_chisq = 2*gh_num*MG_BIN_NUM;
	size_mg_bin = 2*(MG_BIN_NUM + 1);
	size_pair_count = numprocs;

	if (0 == rank)
	{
		MPI_Win_allocate_shared(size_chisq * sizeof(MY_INT), sizeof(MY_INT), MPI_INFO_NULL, MPI_COMM_WORLD, &total_chisq, &win_chisq);

		MPI_Win_allocate_shared(numprocs * sizeof(MY_INT), sizeof(MY_INT), MPI_INFO_NULL, MPI_COMM_WORLD, &pair_count, &win_pair_count);

		MPI_Win_allocate_shared(size_mg_bin * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &mg_bin, &win_mg_bin);
	}
	else
	{
		int dispu_total;
		MPI_Win_allocate_shared(0, sizeof(MY_INT), MPI_INFO_NULL, MPI_COMM_WORLD, &total_chisq, &win_chisq);
		MPI_Win_shared_query(win_chisq, 0, &size_chisq, &dispu_total, &total_chisq);

		MPI_Win_allocate_shared(0, sizeof(MY_INT), MPI_INFO_NULL, MPI_COMM_WORLD, &pair_count, &win_pair_count);
		MPI_Win_shared_query(win_pair_count, 0, &size_pair_count, &dispu_total, &pair_count);

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &mg_bin, &win_mg_bin);
		MPI_Win_shared_query(win_mg_bin, 0, &size_mg_bin, &dispu_total, &mg_bin);

	}
	MPI_Barrier(MPI_COMM_WORLD);

	pair_count[rank] = my_num;
	if (rank == 0)
	{
		initialize_arr(total_chisq, 2*gh_num*MG_BIN_NUM, 0);
		for (i = 0; i < MG_BIN_NUM+1; i++)
		{
			mg_bin[i] = mg1_bin[i];
			mg_bin[i+MG_BIN_NUM+1] = mg2_bin[i];
		}
	}
	initialize_arr(my_chisq, 2*gh_num*MG_BIN_NUM, 0);
	MPI_Barrier(MPI_COMM_WORLD);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	mg = new double[pair_count[rank]];
	my_mgt = new double[pair_count[rank]];
	my_mgx = new double[pair_count[rank]];
	my_mnut = new double[pair_count[rank]];
	my_mnux = new double[pair_count[rank]];

	for (i = 0; i < MG_BIN_NUM+1; i++)
	{
		mg1_bin[i] = mg_bin[i];
		mg2_bin[i] = mg_bin[i + MG_BIN_NUM+1];
	}
	if (rank == 0)
	{
		for (i = 0; i < pair_count[0];i++)
		{
			my_mgt[i] = All_mgt[i];
			my_mgx[i] = All_mgx[i];
			my_mnut[i] = All_mnut[i];
			my_mnux[i] = All_mnux[i];
		}
		mg_min = *std::min_element(All_mgt, All_mgt + total_pair_num);
		mg_max = *std::max_element(All_mgt, All_mgt + total_pair_num);
		std::cout << "G1[" << mg_min << ", " << mg_max << "]" << " bin: " << std::endl;;
		show_arr(mg1_bin, 1, MG_BIN_NUM+1);
		show_arr(mg_bin, 1, MG_BIN_NUM + 1);

		mg_min = *std::min_element(All_mgx, All_mgx + total_pair_num);
		mg_max = *std::max_element(All_mgx, All_mgx + total_pair_num);
		std::cout << "G1[" << mg_min << ", " << mg_max << "]" << " bin: " << std::endl;;
		show_arr(mg2_bin, 1, MG_BIN_NUM+1);
		show_arr(mg_bin+MG_BIN_NUM+1, 1, MG_BIN_NUM + 1);

		std::cout << "Num of each thread: ";
		show_arr(pair_count, 1, numprocs);
		sum_arr(pair_count, numprocs, 0, numprocs, temp_num);
		std::cout << "Total num: " << temp_num << " [" << total_pair_num << "]." << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// distribute the data to each rank
	// (developing) if the catalog is too big, the pairs will be to many to store in the memory.
	// split the whole cata into many sub-sets, then loop the sub-sets,
	// in each loop, distribute the set to each rank.

	MPI_Status status;
#define mysyn( a, b ) \
	{ \
		if (rank == 0) \
		{\
			for (i = 1; i < numprocs; i++)\
			{\
				MPI_Send((a) + pair_count[i - 1], pair_count[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);\
			}\
		}\
		else\
		{\
			MPI_Recv((b),  pair_count[rank], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);\
		}\
	}

	mysyn(All_mgt, my_mgt);
	mysyn(All_mgx, my_mgx);
	mysyn(All_mnut, my_mnut);
	mysyn(All_mnux, my_mnux);
#undef mysyn

	if (rank == 0)
	{
		std::cout << "Begin chi calculation." << std::endl;
	}
	// calculate the chi square
	MPI_Barrier(MPI_COMM_WORLD);

	st1 = clock();
	for(i = 0; i < gh_num; i++)
	{
		for (j = 0; j < pair_count[rank]; j++)
		{
			mg[j] = my_mgt[j] - gh[i] * my_mnut[j];			
		}
		histogram(mg, mg1_bin, data_count, pair_count[rank], MG_BIN_NUM);
		for (k = 0; k < MG_BIN_NUM; k++)
		{
			my_chisq[i*MG_BIN_NUM + k] += data_count[k];
		}

		for (j = 0; j < pair_count[rank]; j++)
		{
			mg[j] = my_mgx[j] - gh[i] * my_mnux[j];	
		}
		histogram(mg, mg2_bin, data_count, pair_count[rank], MG_BIN_NUM);
		for (k = 0; k < MG_BIN_NUM; k++)
		{
			my_chisq[(i + gh_num)*MG_BIN_NUM + k] += data_count[k];
		}		
	}
	if (rank == 0)
	{
		std::cout << "Finish chi calculation." << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for (i = 0; i < numprocs; i++)
	{
		if (i == rank)
		{
			for (j = 0; j < 2 * gh_num*MG_BIN_NUM; j++)
			{
				total_chisq[j] += my_chisq[j];
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if (rank == 0)
	{
		std::cout << "Begin chi square calculation." << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	st2 = clock();

	// fit the chis square
	if (rank == 0)
	{
		double *final_chisq_t = new double[gh_num] {};
		double *final_chisq_x = new double[gh_num] {};
		MY_INT temp_chi_t[MG_BIN_NUM]{}, temp_chi_x[MG_BIN_NUM]{};
		double signal[4]{ -1, -1,-1,-1 };

		for (i = 0; i < gh_num; i++)
		{
			for (j = 0; j < MG_BIN_NUM; j++)
			{
				temp_chi_t[j] = total_chisq[i*MG_BIN_NUM + j];
				temp_chi_x[j] = total_chisq[(i + gh_num)*MG_BIN_NUM + j];
			}
			try
			{
				cal_chisq_1d(temp_chi_t, MG_BIN_NUM, final_chisq_t[i]);
			}
			catch (const char* msg)
			{
				std::cout << "Faliure in tangential shear chi sqaure calculation. " << std::endl;
				std::cout << gh[i]<<" "<< msg << std::endl;
				final_chisq_t[i] = -99;
			}
			try
			{
				cal_chisq_1d(temp_chi_x, MG_BIN_NUM, final_chisq_x[i]);
			}
			catch (const char* msg)
			{
				std::cout << "Faliure in tangential shear chi sqaure calculation. "<< std::endl;
				std::cout << gh[i] << " " << msg << std::endl;
				final_chisq_x[i] = -99;
			}

		}
		//std::cout << "The tan chi square:";
		//show_arr(final_chisq_t, 1, gh_num);
		//std::cout << std::endl;
		//std::cout << "The tan chi square:";
		//show_arr(final_chisq_x, 1, gh_num);
		try
		{
			fit_shear(gh, final_chisq_t, gh_num, signal[0], signal[1], 50);
		}
		catch (const char* msg)
		{
			std::cerr << "Faliure in tangential shear fitting. " << msg << std::endl;
			exit(0);
		}
		try
		{
			fit_shear(gh, final_chisq_x, gh_num, signal[2], signal[3], 50);
		}
		catch (const char* msg)
		{
			std::cerr << "Faliure in cross shear fitting. " << msg << std::endl;
			exit(0);
		}
		if (area_num > 1)
		{
			sprintf(result_path, "%sresult/total_%d.hdf5", total_path, radius_id);
		}
		else
		{
			sprintf(result_path, "%sresult/w_%d_%d.hdf5", total_path, area_id[0], radius_id);
		}

		sprintf(set_name, "/chisq_t");
		write_h5(result_path, set_name, final_chisq_t, 1, gh_num, TRUE);
		sprintf(set_name, "/chisq_x");
		write_h5(result_path, set_name, final_chisq_x, 1, gh_num, FALSE);

		sprintf(set_name, "/shear");
		write_h5(result_path, set_name, gh, 1, gh_num, FALSE);

		sprintf(set_name, "/signal");
		write_h5(result_path, set_name, signal, 4, 1, FALSE);
		std::cout << "Signal: ";
		show_arr(signal, 1, 4);
	}


	MPI_Win_free(&win_mg_bin);
	MPI_Win_free(&win_chisq);
	MPI_Win_free(&win_pair_count);

	if (rank == 0)
	{
		for (i = 0; i < area_num; i++)
		{
			delete[] data[i];
		}
		delete All_mgt;
		delete All_mgx;
		delete All_mnut;
		delete All_mnux;
	}

	delete[] gh;
	delete[] my_mgt;
	delete[] my_mgx;
	delete[] my_mnut;
	delete[] my_mnux;
	delete[] mg;

	MPI_Finalize();
	return 0;
}