#include<FQlib.h>
#include<hk_iolib.h>
#include<hk_mpi.h>

void task_allot(const int total_task_num, const int division_num, const int my_part_id, int &my_st_id, int &my_ed_id, int *task_count)
{
    int i,j,m,n;
    m = total_task_num/division_num;
    n = total_task_num%division_num;

    for(i=0;i<division_num;i++)
    {
        task_count[i] = m;
        if(i<n){task_count[i] +=1;}
    }
    m=0;
    n=0;
    for(i=0;i<my_part_id;i++)
    {
        m+=task_count[i];
    }
    n = m+task_count[my_part_id];
    my_st_id = m;
    my_ed_id = n;
}

int main(int argc, char**argv)
{   
    int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

    char parent_path[400], shear_path[400], result_path[400];
    char set_name[30], inform[300], time_now[40];
    char result_nm[60];
    char data_path_1[400], data_path_2[400],data_path_3[400], data_type_1[40], data_type_2[40], data_type_3[40];

    int i,j,k;

    int shear_num, shear_st, shear_ed;
    int*shear_point;
    double *shears, *g1_t, *g2_t;
    double gh1, gh1_sig, gh2, gh2_sig;

    int mg_bin_num;
    double *mg_bins;
    int chi_check_num;
    double *chi_check;
    double chisq_min_fit;
    double left_guess, right_guess;
  
    double *total_chi_check_g1, *sub_chi_check_g1;
    double *total_chi_check_g2, *sub_chi_check_g2;
    int *chi_send_count;

    int data_col, data_row;
    int result_col;
    double *data, *mg1, *mg2,*mn,*mnu1, *mnu2;
    double weight;

    // for results
    int sub_num;
    double *result_all_mean, *result_sub_mean;
    double *result_all_pdf, *result_sub_pdf;
    int *send_count;
    double *rotation;

    // data shape
    strcpy(parent_path, argv[1]);
    strcpy(result_nm, argv[2]);
    strcpy(data_type_1, argv[3]);
    strcpy(data_type_2, argv[4]);
    int shear_rotation_cmd = atoi(argv[5]);
    shear_num = atoi(argv[6]);
    data_row = atoi(argv[7])*10000;

    data_col = 5;// G1, G2, N, U
    result_col = 4;// g1, g1_sig, g2, g2_sig
    mg_bin_num = 10;//atoi(argv[4]);
    chi_check_num =20;
    left_guess = -0.1;
    right_guess = 0.1;
    
    rotation = new double[5];

    //sprintf(result_path, "%s/shear_result_%s_fit_range_%.4f.hdf5", parent_path, data_type, fit_range[fit_range_label]);
    // sprintf(result_path, "%s/shear_result_%s_rotate_%s_minus_%s.hdf5", parent_path, data_type_1,data_type_2, data_type_3);
    sprintf(result_path, "%s/shear_result_%s.hdf5", parent_path, result_nm);

    data = new double[data_row*data_col]{};
    mg1 = new double[data_row]{};
    mg2 = new double[data_row]{};
    mn = new double[data_row]{};
    mnu1 = new double[data_row]{};
    mnu2 = new double[data_row]{};

    shear_point = new int[numprocs]{};
    send_count = new int[numprocs]{};
    chi_send_count = new int[numprocs]{};

    mg_bins = new double[mg_bin_num+1]{};
    chi_check = new double[2*chi_check_num]{};
 
    // shear point distribution
    //exit(0);
    i = shear_num/numprocs;
    j = shear_num%numprocs;
    for(k=0;k<numprocs;k++)
    {
        shear_point[k] = i;
    }
    for(k=0;k<j;k++)
    {
        shear_point[k] += 1;
    }
    shear_st=0;
    shear_ed=0;
    for(k=0;k<rank;k++)
    {
        shear_st += shear_point[k];
    }
    shear_ed = shear_st + shear_point[rank];

    for(i=0;i<numprocs;i++)
    {
        send_count[i] = shear_point[i]*4;
        chi_send_count[i] = shear_point[i]*40;
    }

    // the measured g1, sig1, g2, sig2 of each thread
    // measured from all source
    // will be sent to rank 0 stack into result_all
    result_sub_mean = new double[(shear_ed - shear_st)*result_col]{};
    result_sub_pdf = new double[(shear_ed - shear_st)*result_col]{};

    if(rank == 0)
    {
        show_arr(shear_point, 1, numprocs);
        std::cout<<"Bin_num "<<mg_bin_num<<" data:"<<parent_path<<" "<<std::endl;
    }

    sub_chi_check_g1 = new double[(shear_ed - shear_st)*chi_check_num*2];
    sub_chi_check_g2 = new double[(shear_ed - shear_st)*chi_check_num*2];


    MPI_Barrier(MPI_COMM_WORLD);
    for(i=0;i<numprocs;i++)
    {
        if(i == rank)
        {
            std::cout<<rank<<" "<<shear_st<<" "<<shear_ed<<" "<<data_row<<" "<<data_col<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    
    sprintf(shear_path, "%s/shear.hdf5",parent_path);
    g1_t = new double[shear_num]{};
    g2_t = new double[shear_num]{};
    sprintf(set_name,"/g1");
    read_h5(shear_path, set_name,g1_t);
    sprintf(set_name,"/g2");
    read_h5(shear_path, set_name,g2_t);

    // read and calculate
    sprintf(set_name,"/data");
    for(i=shear_st;i<shear_ed;i++)
    {   
        sprintf(data_path_1,"%s/data_%s_%d.hdf5",parent_path, data_type_1, i);
        read_h5(data_path_1, set_name, data);
        // read data
        if(rank == 0)
        {
            std::cout<<data_path_1<<std::endl;
        }
        for(k=0;k<data_row;k++)
        {
            mg1[k] = data[k*data_col];
            mg2[k] = data[k*data_col + 1];
            mn[k] = data[k*data_col + 2];
            mnu1[k] = data[k*data_col + 2] + data[k*data_col + 3];
            mnu2[k] = data[k*data_col + 2] - data[k*data_col + 3];
        }

        sprintf(data_path_2,"%s/data_%s_%d.hdf5",parent_path, data_type_2, i);
        read_h5(data_path_2, set_name, data);
        // read data
        if(rank == 0)
        {
            std::cout<<data_path_2<<std::endl;
        }
        for(k=0;k<data_row;k++)
        {
            if(shear_rotation_cmd == 0)
            {
                mg1[k] = mg1[k] + data[k*data_col];
                mg2[k] = mg2[k] + data[k*data_col + 1];
                mn[k] = mn[k] + data[k*data_col + 2];
                mnu1[k] = mnu1[k] + (data[k*data_col + 2] + data[k*data_col + 3]);
                mnu2[k] = mnu2[k] + (data[k*data_col + 2] - data[k*data_col + 3]);
            }
            else
            {
                estimator_rotation(Pi/2, data[k*data_col],data[k*data_col+1],data[k*data_col+2],data[k*data_col+3],data[k*data_col+4],rotation);
                mg1[k] = mg1[k] + rotation[0];
                mg2[k] = mg2[k] + rotation[1];
                mn[k] = mn[k] + rotation[2];
                mnu1[k] = mnu1[k] + rotation[2] + rotation[3];
                mnu2[k] = mnu2[k] + rotation[2] - rotation[3];
            }
        }
        
        // sprintf(data_path_3,"%s/data_%d_%s.hdf5",parent_path, i, data_type_3);
        // read_h5(data_path_3, set_name, data);
        // // read data
        // if(rank == 0)
        // {
        //     std::cout<<data_path_3<<std::endl;
        // }
        // for(k=0;k<data_row;k++)
        // {   
        //     estimator_rotation(Pi/2, data[k*data_col],data[k*data_col+1],data[k*data_col+2],data[k*data_col+3],data[k*data_col+4],rotation);
        //     mg1[k] = mg1[k] + rotation[0];
        //     mg2[k] = mg2[k] + rotation[1];
        //     mn[k] = mn[k] + rotation[2];
        //     mnu1[k] = mnu1[k] + (rotation[2] + rotation[3]);
        //     mnu2[k] = mnu2[k] + (rotation[2] - rotation[3]);
        // }


        // MEAN
        find_shear_mean(mg1, mn, data_row, gh1, gh1_sig, 1000,100);
        find_shear_mean(mg2, mn, data_row, gh2, gh2_sig, 1000,100);

        result_sub_mean[(i-shear_st)*result_col] = gh1;
        result_sub_mean[(i-shear_st)*result_col + 1] = gh1_sig;
        result_sub_mean[(i-shear_st)*result_col + 2] = gh2;
        result_sub_mean[(i-shear_st)*result_col + 3] = gh2_sig;
        sprintf(inform, "%03d, Ave. True g1: %9.6f, Est.: %9.6f (%8.6f), True g2: %9.6f, Est.: %9.6f (%8.6f).", rank, g1_t[i], gh1, gh1_sig,g2_t[i],gh2,gh2_sig);
        std::cout<<inform<<std::endl;

        // PDF_SYM
        try
        {
            find_shear(mg1, mnu1, data_row, mg_bin_num, gh1, gh1_sig, chisq_min_fit,chi_check, chi_check_num);
            //find_shear_fit(mg1, mnu1, data_row, mg_bin_num, chi_check_num,chi_check ,left_guess, right_guess, gh1, gh1_sig);
            for(k=0;k<2*chi_check_num;k++)
            {
                sub_chi_check_g1[(i-shear_st)*2*chi_check_num + k] = chi_check[k];
            }
        }
        catch(const char*img)
        {
            std::cout<<rank<<" "<<i<<" PDF of g1 is going wrong"<<std::endl;
        }

        try
        {
            find_shear(mg2, mnu2, data_row, mg_bin_num, gh2, gh2_sig, chisq_min_fit,chi_check, chi_check_num);
            //find_shear_fit(mg2, mnu2, data_row, mg_bin_num, chi_check_num, chi_check, left_guess, right_guess,  gh2, gh2_sig);
            for(k=0;k<2*chi_check_num;k++)
            {
                sub_chi_check_g2[(i-shear_st)*2*chi_check_num + k] = chi_check[k];
            }
        }
        catch(const char*img)
        {
            std::cout<<rank<<" "<<i<<" PDF of g2 is going wrong"<<std::endl;
        } 

        result_sub_pdf[(i-shear_st)*result_col] = gh1;
        result_sub_pdf[(i-shear_st)*result_col + 1] = gh1_sig;
        result_sub_pdf[(i-shear_st)*result_col + 2] = gh2;
        result_sub_pdf[(i-shear_st)*result_col + 3] = gh2_sig;
        
        sprintf(inform, "%03d, PDF. True g1: %9.6f, Est.: %9.6f (%8.6f), True g2: %9.6f, Est.: %9.6f (%8.6f).",rank, g1_t[i], gh1,gh1_sig,g2_t[i],gh2,gh2_sig);
        std::cout<<inform<<std::endl;

    }

    // for(i=0;i<numprocs;i++)
    // {
    //     if(rank == i)
    //     {
    //         std::cout<<rank<<std::endl;
    //         show_arr(result_sub_pdf, shear_ed - shear_st, result_col);
    //         show_arr(result_sub_mean, shear_ed - shear_st, result_col);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
    {
        result_all_mean = new double[shear_num*result_col]{};
        result_all_pdf = new double[shear_num*result_col]{};

        total_chi_check_g1 = new double[shear_num*2*chi_check_num]{};
        total_chi_check_g2 = new double[shear_num*2*chi_check_num]{};
    }

    my_Gatherv(result_sub_mean, send_count, result_all_mean, numprocs, rank, 0);
    my_Gatherv(result_sub_pdf, send_count, result_all_pdf, numprocs, rank, 0);

    my_Gatherv(sub_chi_check_g1, chi_send_count, total_chi_check_g1, numprocs, rank, 0);
    my_Gatherv(sub_chi_check_g2, chi_send_count, total_chi_check_g2, numprocs, rank, 0);



    if(rank == 0)
    {   
        double *mc = new double[4];
        double *pdf_mc = new double[8];
        double *mean_mc = new double[8];
        double *fit_val = new double[shear_num];
        double *fit_err = new double[shear_num];
        double *result_arr = new double[6*shear_num]{};

        // mean
        for(k=0;k<2;k++)
        {
            for(j=0;j<shear_num;j++)
            {
                fit_val[j] = result_all_mean[j*result_col + k*2];
                fit_err[j] = result_all_mean[j*result_col + k*2 +1];
            }
            if(k ==0){poly_fit_1d(g1_t, fit_val, fit_err, shear_num, mc, 1);}
            else{poly_fit_1d(g2_t, fit_val, fit_err, shear_num, mc, 1);}

            mean_mc[k*4] = mc[2]-1;// m
            mean_mc[k*4 + 1] = mc[3];// m_sig
            mean_mc[k*4 + 2] = mc[0];// c
            mean_mc[k*4 + 3] = mc[1];// c_sig
        }
        for (i=0;i<shear_num;i++)
        {
            result_arr[i] = g1_t[i];          
            result_arr[shear_num + i] = result_all_mean[i*result_col];
            result_arr[2*shear_num + i] = result_all_mean[i*result_col+1];
            result_arr[3*shear_num + i] = g2_t[i];
            result_arr[4*shear_num + i] = result_all_mean[i*result_col+2];
            result_arr[5*shear_num + i] = result_all_mean[i*result_col+3];
        }  

        std::cout<<"AVERAGE: m & c"<<std::endl;
        show_arr(mean_mc, 2, 4);
        sprintf(set_name,"/mean_result");
        write_h5(result_path, set_name, result_arr, 6, shear_num, true);
        sprintf(set_name,"/mean_mc");
        write_h5(result_path,set_name, mean_mc, 2, 4, false);

        // PDF_SYM
        for(k=0;k<2;k++)
        {
            for(j=0;j<shear_num;j++)
            {
                fit_val[j] = result_all_pdf[j*result_col + k*2];
                fit_err[j] = result_all_pdf[j*result_col + k*2 +1];
            }
            if(k ==0){poly_fit_1d(g1_t, fit_val, fit_err, shear_num, mc, 1);}
            else{poly_fit_1d(g2_t, fit_val, fit_err, shear_num, mc, 1);}

            pdf_mc[k*4] = mc[2]-1;// m
            pdf_mc[k*4 + 1] = mc[3];// m_sig
            pdf_mc[k*4 + 2] = mc[0];// c
            pdf_mc[k*4 + 3] = mc[1];// c_sig
        }

        std::cout<<"PDF_SYM: m & c"<<std::endl;
        show_arr(pdf_mc, 2, 4);

        
        //sprintf(result_path, "%s/result/result_%s_bin_num_%d.hdf5", parent_path, data_type, mg_bin_num);
        for (i=0;i<shear_num;i++)
        {
            result_arr[i] = g1_t[i];          
            result_arr[shear_num + i] = result_all_pdf[i*result_col];
            result_arr[2*shear_num + i] = result_all_pdf[i*result_col+1];
            result_arr[3*shear_num + i] = g2_t[i];
            result_arr[4*shear_num + i] = result_all_pdf[i*result_col+2];
            result_arr[5*shear_num + i] = result_all_pdf[i*result_col+3];
        }
        sprintf(set_name,"/sym_result");
        write_h5(result_path, set_name, result_arr, 6, shear_num, false);
        sprintf(set_name,"/sym_mc");
        write_h5(result_path,set_name, pdf_mc, 2, 4, false);
        sprintf(set_name,"/mc1");
        write_h5(result_path,set_name, pdf_mc, 1, 4, false);
        sprintf(set_name,"/mc2");
        write_h5(result_path,set_name, &pdf_mc[4], 1, 4, false);
        sprintf(set_name,"/chisq_g1");
        write_h5(result_path,set_name, total_chi_check_g1, shear_num, 2*chi_check_num, false);
        sprintf(set_name,"/chisq_g2");
        write_h5(result_path,set_name, total_chi_check_g2, shear_num, 2*chi_check_num, false);

        delete[] result_arr;
        delete[] mean_mc;
        delete[] pdf_mc;
        delete[] fit_val;
        delete[] fit_err;
        delete[] result_all_mean;
        delete[] result_all_pdf;
        delete[] total_chi_check_g1;
        delete[] total_chi_check_g2;

        get_time(time_now, 40);
        std::cout<<time_now<<std::endl;
        std::cout<<"Bin_num "<<mg_bin_num<<" "<<result_path<<std::endl;

    }
    
    delete[] send_count;
    delete[] result_sub_mean;
    delete[] result_sub_pdf;

    delete[] data;
    delete[] mnu2;
    delete[] mnu1;
    delete[] mn;
    delete[] mg2;
    delete[] mg1;
    delete[] shear_point;
    delete[] g1_t;
    delete[] g2_t;

    MPI_Finalize();
    return 0;
}
