#include<FQlib.h>
#include<mpi.h>

/* argv[1]: the the name of total directory of the data				 */
/* argv[2]: filter name, like sex2_2 ...										     */
/* argv[3]: sigma, like 1.5, 2, ..													     */
/* argv[4]: selection name, like mag_auto, flux2_ex1, ..			     */


// an array of struct will store the index of each cutoff criterion
struct cutoff_info
{
	std::string cut_name; // flux2_ex1,..., mag_auto,snr,...
	std::string source;// sex or fourier
	int idx_in_cata, data_col;
	int cut_idx;
};

int main(int argc, char**argv)
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	double st1, st2, st3, st4;
	st1 = clock();

	int i, j, k;
	int my_shear, my_cut, shear_change;
	int shear_num, cut_num, total_cutoff_cells, sub_cutoff_cells;

	int source_count = 0, label, cut_step, cata_label;

	int total_data_num = 10000000;// the total number of source in one shear point
	int data_col = 7; // column of the shear mesurement results
	int data_s_col = 4, data_f_col = 10;// column of the sex mesurement results and our own "snr" measurement results
	int mg1_idx = 2, mg2_idx = 3, mn_idx = 4, mu_idx = 5;

	double chi_check[30];
	int chi_fit_num  = 30;

	shear_num = 10;
	cut_num = 10;
	total_cutoff_cells = shear_num * cut_num;

	if (numprocs > total_cutoff_cells)
	{
		std::cout << "Only " << total_cutoff_cells << " cpus are needed!!!" << std::endl;
		exit(0);
	}

	char total_path[200], log_inform[200], set_name[30];
	char data_path[200], sex_path[200], cut_result_path[200];
	char sigma[5], filter_name[20], select_name[20];
	

	std::string total_path_s, shear_path_s;
	std::string section_name, para_name;
	std::string filter_name_s, sigma_s, select_name_s;

	char_to_str(argv[1], para_name);
	char_to_str(argv[2], filter_name_s);
	char_to_str(argv[3], sigma_s);
	char_to_str(argv[4], select_name_s);

	sprintf(filter_name, argv[2]);
	sprintf(sigma, argv[3]);
	sprintf(select_name, argv[4]);

	read_config("/home/hkli/work/envs/envs.dat", "selection_bias", para_name+"_path", total_path_s);
	shear_path_s = total_path_s + "parameters/shear.dat";
	std::strcpy(total_path, total_path_s.c_str());


	double *shear = new double[shear_num * 2];
	double *g1_true, *g2_true;
	g1_true = new double[shear_num];
	g2_true = new double[shear_num];
	read_text(shear_path_s, shear, 2 * shear_num);
	for (i = 0; i < shear_num; i++)
	{
		g1_true[i] = shear[i];
		g2_true[i] = shear[i + shear_num];
	}	

	   	  
	int criterion_label;
	cutoff_info cutoffs[8];
	for (i = 0; i < 8; i++)
	{
		if (i < 5)
		{
			cutoffs[i].cut_name = "flux2_ex" + std::to_string( i + 1);
			cutoffs[i].source = "fourier";	
			cutoffs[i].idx_in_cata = 4 + i;
			cutoffs[i].data_col = data_f_col;
		}
		else
		{
			cutoffs[i].source = "sex";
			cutoffs[i].data_col = data_s_col;
		}
		cutoffs[i].cut_idx = i;
	}
	cutoffs[5].cut_name = "mag_auto";
	cutoffs[5].idx_in_cata = 0;
	cutoffs[6].cut_name = "sex_snr";
	cutoffs[6].idx_in_cata = 1;
	cutoffs[7].cut_name = "snr_auto";
	cutoffs[7].idx_in_cata = 2;

	// find the selection criterion
	for (i = 0; i < 8; i++)
	{
		if (select_name_s == cutoffs[i].cut_name)
		{
			criterion_label = i;
		}
	}


	
	// task distribution
	//					cut_1		cut_2		cut_3 ...
	//  shear_1
	//  shear_2
	//		..

	int *tasks = new int[total_cutoff_cells];
	int *my_task = new int[total_cutoff_cells];
	for (i = 0; i < total_cutoff_cells; i++)
	{
		tasks[i] = i;
	}
	task_alloc(tasks, total_cutoff_cells, numprocs, rank, my_task);

	double selects[8];
	double *data = new double[total_data_num*data_col];
	double *data_select[2];
	data_select[0] = new double[total_data_num*data_f_col];// the fourier catalog
	data_select[1] = new double[total_data_num*data_s_col];// the sextractor catalog
	int *mask = new int[total_data_num];
	int *mask_s = new int[total_data_num];

	double *mg1, *mg2, *mnu1, *mnu2;
	int exist_tag=0;
	double gh1, gh1_sig, gh2, gh2_sig;

	double *cut_scale, *shear_result;
	int *source_num;

	MPI_Win win_cut_scale, win_shear, win_num;
	MPI_Aint scale_size, shear_size, num_size;

	shear_size = shear_num * cut_num * 4 * sizeof(double);
	scale_size = cut_num * sizeof(double);
	num_size = shear_num * cut_num * sizeof(int);

	if (0 == rank)
	{	
		// win_shear stores the [ [g1], [g1_sig], [g2], [g2_sig]], each sub "[ ]" is a block of size shear_num*cut_num
		MPI_Win_allocate_shared(shear_size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shear_result, &win_shear);
		// win_cut_scale stores the cutoff thresholds 10 elements.
		MPI_Win_allocate_shared(scale_size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &cut_scale, &win_cut_scale);

		MPI_Win_allocate_shared(num_size, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &source_num, &win_num);
	}
	else
	{
		int disp_unit_s, disp_unit_c, disp_unit_num;

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shear_result, &win_shear);
		MPI_Win_shared_query(win_shear, 0, &shear_size, &disp_unit_s, &shear_result);

		MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &cut_scale, &win_cut_scale);
		MPI_Win_shared_query(win_cut_scale, 0, &scale_size, &disp_unit_c, &cut_scale);

		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &source_num, &win_num);
		MPI_Win_shared_query(win_num, 0, &num_size, &disp_unit_num, &source_num);
	}

	// read all the data and calculate the cutoff thresholds 
	sprintf(set_name, "/data");
	if (0 == rank)
	{
		double *data_cut = new double[total_data_num*shear_num];		

		//read all the data
		for (i = 0; i < shear_num; i++)
		{
			// read the detection mask from sextractor
			sprintf(data_path, "%sresult/data/%s/mask_%d.hdf5", total_path, filter_name, i);
			read_h5(data_path, set_name, mask);
			
			if ("fourier" == cutoffs[criterion_label].source)
			{
				// read the fourier quad results
				sprintf(data_path, "%sresult/data/data_%.1fsig/data_%d.hdf5", total_path, atof(sigma), i);
				cata_label = 0;
			}
			else
			{
				// read the sextractor results
				sprintf(data_path, "%sresult/data/%s/sex_%d.hdf5", total_path, filter_name, i);
				cata_label = 1;
			}

			read_h5(data_path, set_name, data_select[cata_label]);
			std::cout << data_path << std::endl;
			for (j = 0; j < total_data_num; j++)
			{
				if (mask[j] > 0)
				{
					// the column index & the data column
					label = cutoffs[criterion_label].idx_in_cata + j * cutoffs[criterion_label].data_col;
					data_cut[source_count] = data_select[cata_label][label];
					source_count++;
				}
			}
		}

		// sort the selection criterion and calculate the cutoff thresholds
		sort_arr(data_cut, source_count, 1);
		cut_step = source_count / cut_num;
		for (i = 0; i < cut_num; i++)
		{
			cut_scale[i] = data_cut[i*cut_step];
		}

		std::cout << "The selection criterion: " << select_name_s << " " << criterion_label << std::endl;
		//show_arr(cut_scale, 1, cut_num);
		delete[] data_cut;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	st2 = clock();
	
	shear_change = -1;
	for (i = 0; i < total_cutoff_cells; i++)
	{
		if (my_task[i] > -1)
		{
			source_count = 0;
			my_shear = my_task[i] / cut_num;
			my_cut = my_task[i] % cut_num;
			//std::cout << my_task[i] << " " << my_shear << " " << my_cut << std::endl;
			initialize_arr(mask_s, total_data_num, 0);

			// if in different shear point, reload the file in that point
			if (my_shear != shear_change)
			{
				shear_change = my_shear;

				// read detection mask
				sprintf(data_path, "%sresult/data/%s/mask_%d.hdf5", total_path, filter_name, my_shear);
				read_h5(data_path, set_name, mask);
				
				// read the shear measurement result G1, G2 ....
				sprintf(data_path, "%sresult/data/data_%d.hdf5", total_path, my_shear);
				read_h5(data_path, set_name, data);

				if ("fourier" == cutoffs[criterion_label].source)
				{
					// read the fourier quad results, selction criteria
					sprintf(data_path, "%sresult/data/data_%.1fsig/data_%d.hdf5", total_path, atof(sigma), my_shear);
					cata_label = 0;
				}
				else
				{
					// read the sextractor results, selction criteria
					sprintf(data_path, "%sresult/data/%s/sex_%d.hdf5", total_path, filter_name, my_shear);
					cata_label = 1;
				}
				read_h5(data_path, set_name, data_select[cata_label]);

				for (j = 0; j < total_data_num; j++)
				{
					label = cutoffs[criterion_label].idx_in_cata + j * cutoffs[criterion_label].data_col;
					if (mask[j] > 0 and data_select[cata_label][label] >= cut_scale[my_cut])
					{
						mask_s[j] = 1;
						source_count++;
					}
				}
			}
			else
			{
				for (j = 0; j < total_data_num; j++)
				{
					label = cutoffs[criterion_label].idx_in_cata + j * cutoffs[criterion_label].data_col;
					if (mask[j] > 0 and data_select[cata_label][label] >= cut_scale[my_cut])
					{
						mask_s[j] = 1;
						source_count++;
					}
				}
			}
			
			mg1 = new double[source_count];
			mg2 = new double[source_count];
			mnu1 = new double[source_count];
			mnu2 = new double[source_count];

			source_count = 0;
			for (j = 0; j < total_data_num; j++)
			{
				if (mask_s[j] > 0)
				{
					mg1[source_count] = data[j*data_col + mg1_idx];
					mg2[source_count] = data[j*data_col + mg2_idx];
					mnu1[source_count] = data[j*data_col + mn_idx] + data[j*data_col + mu_idx];
					mnu2[source_count] = data[j*data_col + mn_idx] - data[j*data_col + mu_idx];
					source_count++;
				}
			}

			find_shear(mg1, mnu1, source_count, 8, gh1, gh1_sig, chi_check, chi_fit_num);
			find_shear(mg2, mnu2, source_count, 8, gh2, gh2_sig, chi_check, chi_fit_num);
			//std::cout << my_task[i] << " " << my_shear << " " << my_cut << g1_true[my_shear] << " " << gh1 << g2_true[my_shear]<<" "<<gh2<<" "<<source_count << std::endl;


			shear_result[my_shear * 4 * cut_num + my_cut] = gh1;
			shear_result[my_shear * 4 * cut_num + cut_num + my_cut] = gh1_sig;
			shear_result[my_shear * 4 * cut_num + cut_num * 2 + my_cut] = gh2;
			shear_result[my_shear * 4 * cut_num + cut_num * 3 + my_cut] = gh2_sig;
			
			source_num[my_shear * cut_num + my_cut] = source_count;

			delete[] mg1;
			delete[] mg2;
			delete[] mnu1;
			delete[] mnu2;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	st3 = clock();

	if (0 == rank)
	{
		double *fit_gh1 = new double[shear_num];
		double *fit_gh1_sig = new double[shear_num];

		double *fit_gh2 = new double[shear_num];
		double *fit_gh2_sig = new double[shear_num];

		double *mc1_array = new double[4 * cut_num];
		double *mc2_array = new double[4 * cut_num];
		double coeff[4];

		for (i = 0; i < cut_num; i++)
		{
			for (j = 0; j < shear_num; j++)
			{
				fit_gh1[j] = shear_result[j * 4 * cut_num + i];
				fit_gh1_sig[j] = shear_result[j * 4 * cut_num + cut_num+i];

				fit_gh2[j] = shear_result[j * 4 * cut_num + cut_num * 2 + i];
				fit_gh2_sig[j] = shear_result[j * 4 * cut_num + cut_num * 3 + i];
			}

			poly_fit_1d(g1_true, fit_gh1, fit_gh1_sig, shear_num, coeff, 1);
			mc1_array[i] = coeff[2] - 1;// m
			mc1_array[i + cut_num] = coeff[3];//m_sig
			mc1_array[i + cut_num * 2] = coeff[0];//c
			mc1_array[i + cut_num * 3] = coeff[1];//c_sig

			poly_fit_1d(g2_true, fit_gh2, fit_gh2_sig, shear_num, coeff, 1);
			mc2_array[i] = coeff[2] - 1;// m
			mc2_array[i + cut_num] = coeff[3];//m_sig
			mc2_array[i + cut_num * 2] = coeff[0];//c
			mc2_array[i + cut_num * 3] = coeff[1];//c_sig
		}
		
		// save the cutoff scales
		sprintf(data_path, "%sresult/cuts/sym/%s/%s/total.hdf5", total_path, filter_name, select_name);
		sprintf(set_name, "/cut_scale");
		write_h5(data_path, set_name, cut_scale, 1, cut_num, TRUE);
		sprintf(set_name, "/mc1");
		write_h5(data_path, set_name, mc1_array, 4, cut_num, FALSE);
		sprintf(set_name, "/mc2");
		write_h5(data_path, set_name, mc2_array, 4, cut_num, FALSE);
		sprintf(set_name, "/shear");
		write_h5(data_path, set_name, shear_result, shear_num, 4*cut_num, FALSE);
		sprintf(set_name, "/num");
		write_h5(data_path, set_name, source_num, shear_num, cut_num, FALSE);

		delete[] fit_gh1;
		delete[] fit_gh1_sig;
		delete[] fit_gh2;
		delete[] fit_gh2_sig;
		delete[] mc1_array;
		delete[] mc2_array;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	st4 = clock();
	if (0 == rank)
	{
		std::cout << (st2 - st1) / CLOCKS_PER_SEC << " " << (st3 - st2) / CLOCKS_PER_SEC << " " << (st4 - st3) / CLOCKS_PER_SEC << std::endl;
	}
	delete[] mask;
	delete[] data_select[0];
	delete[] data_select[1];
	delete[] data;
	delete[] tasks;
	delete[] my_task;
	MPI_Finalize();
	return 0;
}