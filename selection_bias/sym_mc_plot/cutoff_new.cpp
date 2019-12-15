#include<FQlib.h>
#include<hk_iolib.h>
#include<hk_mpi.h>
#define PDF_SYM

int main(int argc, char**argv)
{
    /* Resolution factor cutoff */
    /* argv[1]: the the name of total directory of the data				 */
    /* argv[2]: filter name, like sex2_2 ...							*/
	/* argv[3]: name of selection criterion 							*/
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	double st1, st2, st3, st4;
	st1 = clock();

	int i, j, k;
	int my_shear, my_cut, cell_st, cell_ed, shear_change;
	int shear_num, cut_num, total_cutoff_cells, sub_cutoff_cells;

	int source_count, label, cut_step, cata_label;

	int total_data_num;
	int data_col; 
	int mg1_idx = 2, mg2_idx = 3, mn_idx = 4, mu_idx = 5;

	double chi_check[40];// be the 2Xchi_fit_num
	int chi_fit_num;

    total_data_num = 10000000;// the total number of source in one shear point
    data_col = 7;// column of the shear mesurement results
    chi_fit_num  = 20;

	shear_num = 10;
	cut_num = 10;
	total_cutoff_cells = shear_num * cut_num;

	if (numprocs > total_cutoff_cells)
	{
		std::cout << "Only " << total_cutoff_cells << " cpus are needed!!!" << std::endl;
		exit(0);
	}

	char total_path[300], mask_path[300], selection_path[300], data_path[300], result_path[300], shear_path[300],log_inform[300], set_name[30];
    char source_name[50], filter_name[50], select_name[50];
	std::string select_name_s;

	std::strcpy(total_path, argv[1]);
	std::strcpy(filter_name, argv[2]);
	std::strcpy(select_name, argv[3]);
	char_to_str(select_name, select_name_s);

    //sprintf(total_path,"/mnt/perc/hklee/selection_bias/%s", source_name);
    sprintf(shear_path,"%s/parameters/shear.dat",total_path);
	if(rank==0)
    std::cout<<total_path<<" "<<shear_path<<" "<<std::endl;

	double *shear = new double[shear_num * 2];
	double *g1_true, *g2_true;
	g1_true = new double[shear_num];
	g2_true = new double[shear_num];
    std::string shear_path_s;
    char_to_str(shear_path,shear_path_s);


	read_text(shear_path_s, shear, 2 * shear_num);
    
	for (i = 0; i < shear_num; i++)
	{
		g1_true[i] = shear[i];
		g2_true[i] = shear[i + shear_num];
	}	

	////////////////////////////////////////////////////////////////////////////////
	///////////////////////////// task distribution ///////////////////////////////
	//					cut_1		cut_2		cut_3 ...
	//  shear_1
	//  shear_2
	//		..
    i = total_cutoff_cells/numprocs;
    j = total_cutoff_cells%numprocs;
    cell_st = i*rank;
    cell_ed = i*(rank+1);
    if(rank == numprocs-1)
    {
        cell_ed += j;
    }
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////


	double *data;
	int *mask; 
    double *selection; 

	double *mg1, *mg2, *mnu1, *mnu2;
	double gh1, gh1_sig, gh2, gh2_sig;
	double left, right;

	double *cut_scale, *shear_result;
	int *source_num;

	MPI_Win win_cut_scale, win_shear, win_num;
	MPI_Aint scale_size, shear_size, num_size;

	shear_size = shear_num * cut_num * 4 * sizeof(double);
	scale_size = cut_num * sizeof(double);
	num_size = shear_num * cut_num * sizeof(int);

	if (0 == rank)
	{	
		// win_shear stores the block [ [g1_cut1, g1_cut2, g1_cut3, ...], 
        //                              [g1_sig,......], 
        //                              [g2,..........], 
        //                              [g2_sig,......]], 4*cut_num elements, there are "shear_num" blocks
		MPI_Win_allocate_shared(shear_size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shear_result, &win_shear);
		// win_cut_scale stores the cutoff thresholds 10 elements.
		MPI_Win_allocate_shared(scale_size, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &cut_scale, &win_cut_scale);
        // stores the source number of each cutoff step of each shear point
        // [[num_cut1, num_cut2, ...], of shear point 0
        //  [num_cut1, num_cut2, ...], of shear point 1
        //  ....]
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

	////////////////// read all the data and calculate the cutoff thresholds /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	sprintf(set_name, "/data");
	
	data = new double[total_data_num*data_col];
	mask = new int[total_data_num];
    selection = new double[total_data_num];

	if (0 == rank)
	{
		double *data_cut = new double[total_data_num*shear_num];		

		//read all the data
	    source_count = 0;
    	for (i = 0; i < shear_num; i++)
		{   
			// read the detection mask from sextractor
			sprintf(mask_path, "%s/result/data/%s/mask_%d.hdf5", total_path, filter_name, i);
			read_h5(mask_path, set_name, mask);

			sprintf(selection_path, "%s/result/data/%s/%s_%d.hdf5", total_path, filter_name, select_name, i);
			read_h5(selection_path, set_name, selection);

			std::cout << selection_path << std::endl;
			if(select_name_s == "mag_auto" or select_name_s == "mag_true")
			{
				for (j = 0; j < total_data_num; j++)
				{	
					if (mask[j] > 0)
					{
						data_cut[source_count] = selection[j];
						source_count++;
					}				
				}
			}
			else
			{
				for (j = 0; j < total_data_num; j++)
				{	
					if (mask[j] > 0 and selection[j] > 0)
					{
						data_cut[source_count] = selection[j];
						source_count++;
					}				
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

		show_arr(cut_scale, 1, cut_num);
        //std::cout<<source_count<<std::endl;
		delete[] data_cut;

	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0)
	{
		//show_arr(cut_scale,1, cut_num);
		//std::cout<<cell_st<<" "<<cell_ed<<std::endl;
		std::cout<<"Cutting off"<<std::endl;

	}

	st2 = clock();

	shear_change = -1;
	for (i = cell_st; i < cell_ed; i++)
	{	
        my_shear = i/cut_num;
        my_cut = i % cut_num;

        //if(rank == 0)
        //{std::cout << i << " " << my_shear << " " << my_cut << std::endl;}

        // if in different shear point, reload the file of that point
        if (my_shear != shear_change)
        {
			initialize_arr(mask, total_data_num, 0);

            shear_change = my_shear;

            // read detection mask
            sprintf(mask_path, "%s/result/data/%s/mask_%d.hdf5", total_path, filter_name, my_shear);
            read_h5(mask_path, set_name, mask);
            //if(rank==0){std::cout<<mask_path<<std::endl;}

            // read the shear measurement result G1, G2 ....
            sprintf(data_path, "%s/result/data/data_%d.hdf5", total_path, my_shear);
            read_h5(data_path, set_name, data);
            //if(rank==0){std::cout<<data_path<<std::endl;}

            // read selection criterion
			sprintf(selection_path, "%s/result/data/%s/%s_%d.hdf5", total_path, filter_name,select_name, my_shear);
			read_h5(selection_path, set_name, selection);
			//if(rank==0){std::cout<<selection_path<<std::endl;}

        }

        source_count = 0;
		for (j = 0; j < total_data_num; j++)
		{
			if (mask[j] > 0 and selection[j] >= cut_scale[my_cut])
			{
				source_count++;
			}
		}
		//if(rank ==0){std::cout<<cut_scale[my_cut]<<" "<<source_count<<std::endl;}
		        
        mg1 = new double[source_count];
        mg2 = new double[source_count];
        mnu1 = new double[source_count];
        mnu2 = new double[source_count];

        source_count = 0;
        for (j = 0; j < total_data_num; j++)
        {
            if (mask[j] > 0 and selection[j] >= cut_scale[my_cut])
            {
                mg1[source_count] = data[j*data_col + mg1_idx];
                mg2[source_count] = data[j*data_col + mg2_idx];
                mnu1[source_count] = data[j*data_col + mn_idx] + data[j*data_col + mu_idx];
                mnu2[source_count] = data[j*data_col + mn_idx] - data[j*data_col + mu_idx];
                source_count++;
            }
        }

		try
        {
#ifdef PDF_SYM_RANGE
			if(i == 0){std::cout<<"PDF_SYM Range fit"<<std::endl;}
			find_shear_mean(mg1, mnu1, source_count, gh1, gh1_sig, 100);
			left = gh1 - 0.015;
			right = gh1 + 0.015;
			find_shear_fit(mg1, mnu1, source_count, 40, chi_fit_num, chi_check, left, right, gh1, gh1_sig);
#endif
#ifdef PDF_SYM
			if(i == 0){std::cout<<"PDF_SYM"<<std::endl;}	
			find_shear(mg1, mnu1, source_count, 40, gh1, gh1_sig, chi_check, chi_fit_num);
#endif
#ifdef MEAN	
			if(i == 0){std::cout<<"MEAN"<<std::endl;}		
			find_shear_mean(mg1, mnu1, source_count, gh1, gh1_sig, 100);
#endif
		}
		catch(const char *msg)
		{
			char err_inform[300];
			sprintf(err_inform,"g1, %d, %d, %d, %s",my_shear, my_cut, source_count, msg);
			std::cerr<<err_inform<<std::endl;
		}
		try
        {
#ifdef PDF_SYM_RANGE
			find_shear_mean(mg2, mnu2, source_count, gh2, gh2_sig, 100);
			left = gh2 - 0.015;
			right = gh2 + 0.015;
			find_shear_fit(mg2, mnu2, source_count, 40, chi_fit_num, chi_check, left, right, gh2, gh2_sig);
#endif
#ifdef PDF_SYM
			find_shear(mg2, mnu2, source_count, 40, gh2, gh2_sig, chi_check, chi_fit_num);
#endif
#ifdef MEAN
			find_shear_mean(mg2, mnu2, source_count, gh2, gh2_sig, 100);
#endif
		}
		catch(const char *msg)
		{
			char err_inform[300];
			sprintf(err_inform,"g2, %d, %d, %d, %s",my_shear, my_cut, source_count, msg);
			std::cerr<<err_inform<<std::endl;
		}
        // //std::cout << my_task[i] << " " << my_shear << " " << my_cut << g1_true[my_shear] << " " << gh1 << g2_true[my_shear]<<" "<<gh2<<" "<<source_count << std::endl;

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
	delete[] mask;
	delete[] selection;
	delete[] data;
	MPI_Barrier(MPI_COMM_WORLD);
	st3 = clock();

	if (0 == rank)
	{			
		show_arr(source_num, shear_num, cut_num);
		std::cout<<(st2-st1)/CLOCKS_PER_SEC<<" "<<(st3-st2)/CLOCKS_PER_SEC<<std::endl;
		std::cout<<"Fitting"<<std::endl;

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
		
		std::cout<<"Result: "<<std::endl;
		show_arr(mc1_array, 4, cut_num);
		show_arr(mc2_array, 4, cut_num);

		sprintf(result_path, "%s/result/cuts/sym/%s/%s/total.hdf5", total_path, filter_name, select_name);

		// save the cutoff scales
		sprintf(set_name, "/cut_scale");
		write_h5(result_path, set_name, cut_scale, 1, cut_num, TRUE);
		sprintf(set_name, "/mc1");
		write_h5(result_path, set_name, mc1_array, 4, cut_num, FALSE);
		sprintf(set_name, "/mc2");
		write_h5(result_path, set_name, mc2_array, 4, cut_num, FALSE);
		sprintf(set_name, "/shear");
		write_h5(result_path, set_name, shear_result, shear_num, 4*cut_num, FALSE);
		sprintf(set_name, "/num");
		write_h5(result_path, set_name, source_num, shear_num, cut_num, FALSE);
		sprintf(set_name, "/g1");
		write_h5(result_path, set_name, g1_true, 1, shear_num, FALSE);
		sprintf(set_name, "/g2");
		write_h5(result_path, set_name, g2_true, 1, shear_num, FALSE);
		std::cout<<"Write to "<<result_path<<std::endl;

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
    
    MPI_Win_free(&win_shear);
    MPI_Win_free(&win_cut_scale);
    MPI_Win_free(&win_num);

	MPI_Finalize();
	return 0;
}