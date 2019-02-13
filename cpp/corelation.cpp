#include<FQlib.h>
#include<mpi.h>

int main(int argc, char *argv[])
{
	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	gsl_initialize(123);

	std::string str_shear_num = "shear_num", str_total_num = "total_num";
	char data_path[100], log_inform[250], set_1[50], npy_data_path[200];
	char data_path_1[200], data_path_2[200];

	sprintf(data_path, "/mnt/ddnfs/data_users/hkli/correlation/simu/");
	std::string str_data_path = "/mnt/ddnfs/data_users/hkli/correlation/simu/";
	std::string str_paraf_path = str_data_path + "parameters/para.ini";
	int total_chips, shear_pairs, shear_num;
	read_para(str_paraf_path, str_total_num, total_chips);
	read_para(str_paraf_path, str_shear_num, shear_num);
	shear_pairs = shear_num / 2;

	int i, j, k;
	sprintf(data_path_1, "%sresult/data/data_%d.hdf5", data_path, rank);
	sprintf(data_path_2, "%sresult/data/data_%d.hdf5", data_path, rank+shear_pairs);
	sprintf(npy_data_path, "%smgauss_%d.hdf5", data_path, rank);
	sprintf(set_1, "/data");
	if (rank < shear_pairs)
	{	
		// for multi_gaussian()
		int fit_num = 40, data_num = total_chips*10000, bin_num = 6;
		double cov11;
		double *covs = new double[4]{};
		double *cor_gs = new double[2]{};
		double *mus = new double[2]{};
		double st1, st2, st3;

		double *data_1 = new double[data_num * 7]{};
		double *data_2 = new double[data_num * 7]{};

		double *cor_gs_npy = new double[fit_num*data_num * 2];

		double *gch11 = new double[data_num]{};
		double *gch12 = new double[data_num]{};

		double *gch21 = new double[data_num]{};
		double *gch22 = new double[data_num]{};

		double *chisqs = new double[fit_num*2]{};

		double *bins = new double[bin_num + 1]{};
		int *nums_1 = new int[bin_num*bin_num]{};
		int *nums_2 = new int[bin_num*bin_num]{};

		read_h5(data_path_1, set_1, data_1);
		read_h5(data_path_2, set_1, data_2);
		read_h5(npy_data_path, set_1, cor_gs_npy);

		if (rank == 0)
		{
			std::cout << numprocs << ", " << data_num << std::endl;
		}

		// set bins for histogram
		for (i = 0; i < data_num; i++)
		{
			gch11[i] = data_1[i * 7 + 2];
		}
		set_bin(gch11, data_num, bins, bin_num, 100.);

		for (i =0; i < fit_num; i++)
		{	
			// covariance 
			cov11 = i * 0.00005 - 0.001;
			covs[0] = fabs(4 * cov11);
			covs[1] = cov11;
			covs[2] = cov11;
			covs[3] = fabs(4 * cov11);
			if(rank == 0)
			{	
				std::cout << covs[0] * covs[3] - cov11 * cov11 << std::endl;
				show_arr(covs, 2, 2);
			}
			st1 = clock();

			for (j = 0; j < data_num; j++)
			{
				rand_multi_gauss(covs, mus, 2, cor_gs); // check agnist numpy 
				//cor_gs[0] = cor_gs_npy[i*data_num * 2 + j * 2];
				//cor_gs[1] = cor_gs_npy[i*data_num * 2 + j * 2 + 1];
				// correlation for g1
				gch11[j] = data_1[j * 7 + 2] - cor_gs[0] * (data_1[j * 7 + 4] + data_1[j * 7 + 5]);
				gch12[j] = data_2[j * 7 + 2] - cor_gs[1] * (data_2[j * 7 + 4] + data_2[j * 7 + 5]);
				// correlation for g2
				gch21[j] = data_1[j * 7 + 3] - cor_gs[0] * (data_1[j * 7 + 4] - data_1[j * 7 + 5]);
				gch22[j] = data_2[j * 7 + 3] - cor_gs[1] * (data_2[j * 7 + 4] - data_2[j * 7 + 5]);
			}
			st2 = clock();

			histogram2d(gch11, gch12, bins, bins, nums_1, data_num, bin_num, bin_num);
			histogram2d(gch21, gch22, bins, bins, nums_2, data_num, bin_num, bin_num);

			chisqs[i] = chisq_2d(nums_1, bin_num);
			chisqs[i+fit_num] = chisq_2d(nums_2, bin_num);

			//for (j = 0; j < numprocs; j++)
			//{
			//	if (j == rank)
			//	{
			//		std::cout << rank << std::endl;

			//		std::cout << chisqs[i] << ", " << chisqs[i+fit_num] << std::endl;
			//		show_arr(bins, 1, bin_num + 1);
			//		show_arr(nums_1, bin_num, bin_num);
			//		std::cout << std::endl;
			//		show_arr(nums_2, bin_num, bin_num);
			//	}
			//	MPI_Barrier(MPI_COMM_WORLD);
			//}


			st3 = clock();
			if (rank == 0)
			{
				std::cout <<i<<", "<< (st2 - st1) / CLOCKS_PER_SEC<<", "<<(st3 - st2) / CLOCKS_PER_SEC << std::endl;
			}

		}
		sprintf(data_path_1, "%schisq_%d.hdf5", data_path, rank);
		write_h5(data_path_1, set_1, chisqs, 2, fit_num);

		delete[] data_1;
		delete[] data_2;

		delete[] gch11;
		delete[] gch12;
		delete[] gch21;
		delete[] gch22;
		delete[] chisqs;
		delete[] bins;
		delete[] nums_1;
		delete[] nums_2;
		delete[] covs;
		delete[] mus;
		delete[] cor_gs;
		delete[] cor_gs_npy;
	}

	gsl_free();
	MPI_Finalize();
	return 0;
}