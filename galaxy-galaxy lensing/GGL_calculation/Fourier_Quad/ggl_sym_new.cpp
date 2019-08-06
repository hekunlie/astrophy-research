#include<FQlib.h>
#include<mpi.h>
#include<vector>

#define MAX_DATA_COL 20
#define MAX_AREA 20
#define MG_BIN_NUM 8
#define ALLOT_THRESH 500000 /* it's smaller than the number threshold, 20000000, for saving data separately */
#define CRIT_DENSITY 554.682135528


int main(int argc, char ** argv)
{
	/* crit_coeff is the coefficient of the critical density, see the note for detials */
	
	/* argv[1]: the file name of foreground source												*/											
	/* argv[2] -- argv[n]: area labels, like 1,3, 4														*/

	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);
	MPI_Status mpi_status;

	char total_path[300], data_path[300], set_name[30], result_path[300];
	char fore_source[50];
	char logs[300];

	int i, j, k;
	double st1, st2, st3, st4, st5, st6, dt1, dt2, dt3;

	int ia, areas[MAX_AREA]{}, area_num, area_id, tag;

	int radius_num, radius_id;
	double *radius_bin;

	int *pair_num[500], *subset_num, *pair_each_radius,  dataset_num, set_id, pair_count;
	int data_num, data_col, data_row;
	int small_num;
	int temp[1];
	int mgt_id, mgx_id, mnut_id, mnux_id, dist_id, z_id, crit_id;
	double *data[MAX_DATA_COL], *total_data[MAX_AREA];

	int *total_chisq;
	double mg_bin[MG_BIN_NUM+1];
	int num_in_mg_bin[MG_BIN_NUM];
	double *mgt, *mgx, *mnut, *mnux, *crit, *dist_r, *redshift, mg_temp; 
	double gh_s, gh_e, dg, *gh;
	int gh_num;
	int my_gh_s, my_gh_e, ig;
	double gt, gt_sig, gx, gx_sig;
	double *final_result, *final_chisq_t, *final_chisq_x, *mean_dist;

	// the foreground name
	strcpy(fore_source, argv[1]);
	// the radius label
	area_num = argc - 2;
	if (rank == 0) std::cout << fore_source << ": ";

	for (i = 2; i < argc; i++)
	{
		areas[i - 2] = atoi(argv[i]);

		if (rank == 0) std::cout << areas[i - 2] << " ";
	}

	sprintf(total_path, "/mnt/perc/hklee/CFHT/gg_lensing/result/%s/fourier_cata_old/", fore_source);	

	mgt_id = 0;
	mgx_id = 1;
	mnut_id = 2;
	mnux_id = 3;
	crit_id = 4;
	dist_id = 5;
	z_id = 6;


	// read the radius bin
	sprintf(data_path, "%sw_1/radius_0.hdf5", total_path);
	sprintf(set_name, "/radius_bin");
	read_h5_datasize(data_path, set_name, radius_num);
	radius_bin = new double[radius_num];
	read_h5(data_path, set_name, radius_bin);
	radius_num = radius_num - 1;

	if (rank == 0)
	{
		if (area_num > 1)
		{
			sprintf(result_path, "%sresult/%s_result_total.hdf5", total_path, fore_source);
		}
		else
		{
			sprintf(result_path, "%sresult/%s_result_w_%d.hdf5", total_path, fore_source, areas[0]);
		}
		write_h5(result_path, set_name, radius_bin, 1, radius_num + 1, TRUE);
	}

	// the total pair count in each radius bin of all areas
	pair_each_radius = new int[radius_num];
	initialize_arr(pair_each_radius, radius_num, 0);

	// how many subset of each radius bin in each area
	subset_num = new int[area_num*radius_num]{};
	
	// count the total pair num in each radius bin
	for(radius_id =0; radius_id <radius_num; radius_id++)
	{
		for (ia = 0; ia < area_num; ia++)
		{
			tag = radius_id * area_num + ia;

			sprintf(data_path, "%sw_%d/radius_%d.hdf5", total_path, areas[ia], radius_id);
			sprintf(set_name, "/pair_count");
			read_h5_datasize(data_path, set_name, data_num);

			pair_num[tag] = new int[data_num];

			read_h5(data_path, set_name, pair_num[tag]);

			sum_arr(pair_num[tag], data_num, 0, data_num, pair_count);
			pair_each_radius[radius_id] += pair_count;

			sprintf(set_name, "/subset_num");
			read_h5(data_path, set_name, temp);
			subset_num[tag] = temp[0];

			sprintf(set_name, "/data_col");
			read_h5(data_path, set_name, temp);
			data_col = temp[0];

		}
	}
	if (rank == 0)
	{

		std::cout <<std::endl<< "Radius bin:" << std::endl;
		show_arr(radius_bin, 1, radius_num + 1);
		std::cout << "Pair count:" << std::endl;
		show_arr(pair_each_radius, 1, radius_num);
		std::cout << std::endl;
		//show_arr(subset_num, radius_num, area_num);
		//std::cout << data_col << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);


	//////////////////////////////////////////////// calculate ////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// first, calculate the signal of the radius bin in which pair is too little
	// it will allot those to each threads
	small_num = 0;
	for (radius_id = 0; radius_id < radius_num; radius_id++)
	{
		if (pair_each_radius[radius_id] <= ALLOT_THRESH)
		{
			small_num++;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0)
	{
		// gt, gt_sig, gx, gx_sig, 4 x n
		final_result = new double[4 * radius_num]{};
		mean_dist = new double[radius_num] {};
		initialize_arr(final_result, 4 * radius_num, 0);
		initialize_arr(mean_dist,  radius_num, 0);
	}

	double gh_guess[20]{ 130, 60, 50, 40, 30, 20, 20, 15,15,15,10,10,10,10,10,10,10,10,10,10};
	
	dg = 0.2;
	gh_num = int((gh_guess[0] *2) / dg) + 1;

	MPI_Win win_chisq;
	MPI_Aint  size_chisq;
	size_chisq = 2 * gh_num*MG_BIN_NUM;
	if (0 == rank)
	{
		MPI_Win_allocate_shared(size_chisq * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &total_chisq, &win_chisq);
	}
	else
	{
		int dispu_total;
		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &total_chisq, &win_chisq);
		MPI_Win_shared_query(win_chisq, 0, &size_chisq, &dispu_total, &total_chisq);
	}

	for (radius_id=0; radius_id<radius_num; radius_id++)
	{
		st1 = clock();

		if (rank == 0)
		{
			initialize_arr(total_chisq, 2 * gh_num*MG_BIN_NUM, 0);

			final_chisq_t = new double[gh_num] {};
			final_chisq_x = new double[gh_num] {};
		}
		MPI_Barrier(MPI_COMM_WORLD);


		gh_s = -gh_guess[radius_id];
		gh_e = fabs(gh_s);
		gh_num = int((gh_e - gh_s) / dg) + 1;
		gh = new double[gh_num] {};
		for (i = 0; i < gh_num; i++)
		{
			gh[i] = gh_s + i * dg;
		}

		i = gh_num / numprocs;
		j = gh_num % numprocs;
		my_gh_s = rank * i;
		my_gh_e = (rank + 1)*i;
		if (rank == numprocs - 1)
		{
			my_gh_e += j;
		}
//		std::cout << rank << " " << my_gh_s << " " << my_gh_e << " "<<gh_num<<std::endl;

		if (radius_id < small_num)
		{
			st3 = clock();
			// the data column and row of each area
			int *temp_row = new int[area_num];
			data_num = 0;
			for (ia = 0; ia < area_num; ia++)
			{
				sprintf(data_path, "%sw_%d/radius_%d.hdf5", total_path, areas[ia], radius_id);
				sprintf(set_name, "/pair_data_0");
				read_h5_datasize(data_path, set_name, dataset_num);

				data_row = dataset_num / data_col;
				temp_row[ia] = data_row;

				data_num += data_row;
				total_data[ia] = new double[dataset_num];

				read_h5(data_path, set_name, total_data[ia]);
			}
			st4 = clock();
			// stack data
			mgt = new double[data_num];
			mgx = new double[data_num];
			mnut = new double[data_num];
			mnux = new double[data_num];
			dist_r = new double[data_num];
			redshift = new double[data_num];	
			for (ia = 0; ia < area_num; ia++)
			{
				data_row = 0;
				for (i = 0; i < ia; i++)
				{
					data_row += temp_row[i];
				}
				for (i = 0; i < temp_row[ia]; i++)
				{
					mgt[data_row+i] = total_data[ia][i*data_col + mgt_id]* total_data[ia][i*data_col + crit_id]* CRIT_DENSITY;
					mgx[data_row + i] = total_data[ia][i*data_col + mgx_id] * total_data[ia][i*data_col + crit_id] * CRIT_DENSITY;
					mnut[data_row + i] = total_data[ia][i*data_col + mnut_id];
					mnux[data_row + i] = total_data[ia][i*data_col + mnux_id];
					dist_r[data_row + i] = total_data[ia][i*data_col + dist_id];
					redshift[data_row + i] = total_data[ia][i*data_col + z_id];

				}
			}
			st5 = clock();
			if (rank == 0)
			{
				for (i = 0; i < data_num; i++)
				{
					// times 1000 for precision
					mean_dist[radius_id] += dist_r[i] * 1000;
				}
				mean_dist[radius_id] = mean_dist[radius_id] / pair_each_radius[radius_id] / 1000.;
			}
			// set up bins for PDF_SYM
			initialize_arr(mg_bin, 0, MG_BIN_NUM + 1);
			set_bin(mgt, pair_each_radius[radius_id], mg_bin, MG_BIN_NUM, 10000);

			//std::cout << "Calculation" << std::endl;
			//show_arr(mg_bin, 1, MG_BIN_NUM + 1);


			// calculate the signal
			for (ig = my_gh_s; ig < my_gh_e; ig++)
			{
				// tangential
				for (i = 0; i < data_num; i++)
				{
					mg_temp = mgt[i] - gh[ig] * mnut[i];
					histogram_s(mg_temp, mg_bin, MG_BIN_NUM, k);
					total_chisq[ig*MG_BIN_NUM + k] += 1;
				}
				// cross
				for (i = 0; i < data_num; i++)
				{
					mg_temp = mgx[i] - gh[ig] * mnux[i];
					histogram_s(mg_temp, mg_bin, MG_BIN_NUM, k);
					total_chisq[(gh_num+ig)*MG_BIN_NUM + k] += 1;
				}
			}
			st6 = clock();
			delete[] temp_row;
			for (ia = 0; ia < area_num; ia++)
			{
				delete[] total_data[ia];
			}
			delete[] mgt;
			delete[] mgx;
			delete[] mnut;
			delete[] mnux;
			delete[] dist_r;
			delete[] redshift;
			if (rank == 0)
			{
				std::cout << (st4 - st3) / CLOCKS_PER_SEC << " " << (st5 - st4) / CLOCKS_PER_SEC << " " << (st6 - st5) / CLOCKS_PER_SEC << std::endl;
			}
		}
		else // if the dataset is too big
		{
			// calculate in each subset
			dt1 = 0;
			dt2 = 0;
			for (ia = 0; ia < area_num; ia++)
			{
				sprintf(data_path, "%sw_%d/radius_%d.hdf5", total_path, areas[ia], radius_id);
				tag = radius_id * area_num + ia;

				// loop the sub-dataset
				for (set_id = 0; set_id < subset_num[tag]; set_id++)
				{
					st3 = clock();

					sprintf(set_name, "/pair_data_%d", set_id);
					read_h5_datasize(data_path, set_name, dataset_num);
					data_row = dataset_num / data_col;
					total_data[ia] = new double[dataset_num];

					if (rank == 0)
					{
						read_h5(data_path, set_name, total_data[ia]);
						for (j = 0; j < data_row; j++)
						{
							mean_dist[radius_id] += total_data[ia][j*data_col + dist_id] * 100;
						}
						for (i = 1; i < numprocs; i++)
						{
							MPI_Send(total_data[ia], dataset_num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
						}
					}
					else
					{
						MPI_Recv(total_data[ia], dataset_num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &mpi_status);
					}
					st4 = clock();
					// set up bins for PDF_SYM
					if (set_id == 0)
					{
						mgt = new double[data_row];
						for (j = 0; j < data_row; j++)
						{
							mgt[j] = total_data[ia][j*data_col+mgt_id] * total_data[ia][j*data_col + crit_id] * CRIT_DENSITY;
						}

						initialize_arr(mg_bin, MG_BIN_NUM + 1,0);
						set_bin(mgt, data_row, mg_bin, MG_BIN_NUM, 10000);

						//show_arr(mgt, 1, 100);
						//show_arr(mg_bin, 1, MG_BIN_NUM + 1);
						delete[] mgt;
					}
					st5 = clock();
					//calculate the signal
					for (ig = my_gh_s; ig < my_gh_e; ig++)
					{
						// tangential
						for (j = 0; j < data_row; j++)
						{
							mg_temp = total_data[ia][j*data_col + mgt_id] * total_data[ia][j*data_col + crit_id] * CRIT_DENSITY - gh[ig]* total_data[ia][j*data_col + mnut_id];
							histogram_s(mg_temp, mg_bin, MG_BIN_NUM, k);
							total_chisq[ig*MG_BIN_NUM + k] += 1;
						}
						// cross
						for (j = 0; j < data_row; j++)
						{
							mg_temp = total_data[ia][j*data_col+mgx_id] * total_data[ia][j*data_col + crit_id] * CRIT_DENSITY - gh[ig] * total_data[ia][j*data_col + mnux_id];
							histogram_s(mg_temp, mg_bin, MG_BIN_NUM, k);
							total_chisq[(gh_num + ig)*MG_BIN_NUM + k] += 1;
						}
					}
					st6 = clock();
					dt1 = (st4 - st3) / CLOCKS_PER_SEC;
					dt2 = (st5 - st4) / CLOCKS_PER_SEC;
					dt3 = (st6 - st5) / CLOCKS_PER_SEC;
					delete[] total_data[ia];
				}

				if (rank == 0)
				{
					mean_dist[radius_id] = mean_dist[radius_id] / pair_each_radius[radius_id] / 100;
				}
			}
			//for (i = 0; i < numprocs; i++)
			//{
			//	if (rank == i)
			//	{
			//		std::cout << dt1 << " " << dt2 <<" "<<dt3<< std::endl;
			//	}
			//	MPI_Barrier(MPI_COMM_WORLD);
			//}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		// write down the results
		if (rank == 0)
		{
			// calculate the final chi squared and fitting
			for (ig = 0; ig < gh_num; ig++)
			{
				for (i = 0; i < MG_BIN_NUM; i++)
				{
					num_in_mg_bin[i] = total_chisq[ig*MG_BIN_NUM + i];
				}
				try
				{
					cal_chisq_1d(num_in_mg_bin, MG_BIN_NUM, final_chisq_t[ig]);
				}
				catch (const char* msg)
				{
					std::cout << "Faliure in tangential shear chi sqaure calculation. " << std::endl;
					std::cout << gh[ig] << " " << msg << std::endl;
					final_chisq_t[ig] = -99;
				}
				for (i = 0; i < MG_BIN_NUM; i++)
				{
					num_in_mg_bin[i] = total_chisq[(ig + gh_num)*MG_BIN_NUM + i];
				}
				try
				{
					cal_chisq_1d(num_in_mg_bin, MG_BIN_NUM, final_chisq_x[ig]);
				}
				catch (const char* msg)
				{
					std::cout << "Faliure in cross shear chi sqaure calculation. " << std::endl;
					std::cout << gh[ig] << " " << msg << std::endl;
					final_chisq_x[ig] = -99;
				}
			}
			try
			{
				fit_shear(gh, final_chisq_t, gh_num, final_result[radius_id], final_result[radius_num + radius_id], 30);
			}
			catch (const char* msg)
			{
				std::cout << "Faliure in tangential shear chi sqaure fitting. " << std::endl;
			}
			try
			{
				fit_shear(gh, final_chisq_x, gh_num, final_result[radius_num*2 + radius_id], final_result[radius_num*3 + radius_id], 30);
			}
			catch (const char* msg)
			{
				std::cout << "Faliure in cross shear chi sqaure fitting. " << std::endl;
			}
			sprintf(set_name, "/num_%d", radius_id);
			write_h5(result_path, set_name, total_chisq, 2, gh_num*MG_BIN_NUM, FALSE);

			sprintf(set_name, "/chisq_t_%d", radius_id);
			write_h5(result_path, set_name, final_chisq_t, 1, gh_num, FALSE);
			sprintf(set_name, "/chisq_x_%d", radius_id);
			write_h5(result_path, set_name, final_chisq_x, 1, gh_num, FALSE);

			delete[] final_chisq_t;
			delete[] final_chisq_x;
		}
			
		delete[] gh;
		st2 = clock();

		if (rank == 0)
		{
			gt = final_result[radius_id];
			gt_sig = final_result[radius_num + radius_id];
			gx = final_result[radius_num*2 + radius_id];
			gx_sig = final_result[radius_num*3 + radius_id];

			sprintf(logs, "Radius: [%.4f, %.4f]. Time: %.2f sec\ngt: %.3f (%.4f), gx: %.3f (%.4f). mean R:%.5f Mpc/h\n", radius_bin[radius_id], radius_bin[radius_id+1] ,(st2-st1)/CLOCKS_PER_SEC,gt, gt_sig, gx, gx_sig, mean_dist[radius_id]);
			std::cout << logs << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if (rank == 0)
	{
		sprintf(set_name, "/result");
		write_h5(result_path, set_name, final_result, 4, radius_num, FALSE);

		sprintf(set_name, "/mean_dist");
		write_h5(result_path, set_name, mean_dist, 1, radius_num, FALSE);

		std::cout << "Radius:" << std::endl;
		show_arr(mean_dist, 1, radius_num);
		std::cout << "Result:" << std::endl;
		show_arr(final_result, 4, radius_num);
	}
	MPI_Win_free(&win_chisq);

	MPI_Finalize();
	return 0;
}