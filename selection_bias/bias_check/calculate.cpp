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

    char parent_path[200], shear_path[200], data_path[200], result_path[200];
    char set_name[30], data_type[40], inform[300];

    int i,j,k;

    int shear_num, shear_st, shear_ed;
    int*shear_point;
    double *shears, *g1_t, *g2_t;
    double gh1, gh1_sig, gh2, gh2_sig;

    int mg_bin_num;
    int chi_check_num;
    double *chi_check;
    double left_guess, right_guess;

    double *total_chi_check_g1, *sub_chi_check_g1;
    double *total_chi_check_g2, *sub_chi_check_g2;
    int *chi_send_count;

    int data_col, data_row;
    int result_col;
    double *data, *mg1, *mg2,*mn,*mnu1, *mnu2;

    // for results
    int sub_num;
    double *result_all_mean, *result_sub_mean;
    double *result_all_pdf, *result_sub_pdf;
    int *send_count;

    // data shape
    strcpy(parent_path, argv[1]);
    shear_num = atoi(argv[2]);
    strcpy(data_type, argv[3]); // noise_free or noisy
    data_row = atoi(argv[4])*10000;
    data_col = 4;// G1, G2, N, U
    result_col = 4;// g1, g1_sig, g2, g2_sig
    mg_bin_num = 8;
    chi_check_num =20;
    left_guess = -0.04;
    right_guess = 0.04;

    // MPI_Barrier(MPI_COMM_WORLD);
    // for(i=0;i<numprocs;i++)
    // {
    //     if(i == rank)
    //     {
    //         std::cout<<rank<<" "<<data_row<<" "<<data_col<<std::endl;
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    data = new double[data_row*data_col];
    mg1 = new double[data_row];
    mg2 = new double[data_row];
    mn = new double[data_row];
    mnu1 = new double[data_row];
    mnu2 = new double[data_row];

    shear_point = new int[numprocs]{};
    send_count = new int[numprocs]{};
    chi_send_count = new int[numprocs]{};

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
        sprintf(data_path,"%s/data_%d_%s.hdf5",parent_path, i, data_type);
        read_h5(data_path, set_name, data);

        for(k=0;k<numprocs;k++)
        {
            if(k == rank)
            {
                std::cout<<rank<<" shear: "<<k<<" data:"<<data_path<<" g1: "<<g1_t[k]<<" g2: "<<g2_t[k]<<std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        // read data
        for(k=0;k<data_row;k++)
        {
            mg1[k] = data[k*data_col];
            mg2[k] = data[k*data_col + 1];
            mn[k] = data[k*data_col + 2];
            mnu1[k] = data[k*data_col + 2] + data[k*data_col + 3];
            mnu2[k] = data[k*data_col + 2] - data[k*data_col + 3];
        }

        // MEAN
        find_shear_mean(mg1, mn, data_row, gh1, gh1_sig, 1000,100);
        find_shear_mean(mg2, mn, data_row, gh2, gh2_sig, 1000,100);

        result_sub_mean[(i-shear_st)*result_col] = gh1;
        result_sub_mean[(i-shear_st)*result_col + 1] = gh1_sig;
        result_sub_mean[(i-shear_st)*result_col + 2] = gh2;
        result_sub_mean[(i-shear_st)*result_col + 3] = gh2_sig;
        sprintf(inform, "%d, Ave. g1: %.6f, %.6f (%.6f), g2: %.6f, %.6f (%.6f).", rank,g1_t[i], gh1,gh1_sig,g2_t[i],gh2,gh2_sig);
        std::cout<<inform<<std::endl;

        // PDF_SYM
        try
        {
            //find_shear(mg=mg1, mnu=mnu1, data_num=data_row, bin_num=8, gh=gh1, gh_sig=gh1_sig, chi_check=chi_check, chi_fit_num=20, ini_left=-0.06, ini_right=0.06);
            find_shear(mg1, mnu1, data_row, mg_bin_num, gh1, gh1_sig, chi_check, chi_check_num, 0, 100, left_guess, right_guess, 30);
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
            find_shear(mg2, mnu2, data_row, mg_bin_num, gh2, gh2_sig, chi_check, chi_check_num, 0, 100, left_guess, right_guess,30);
            //find_shear(mg=mg2, mnu=mnu2, data_num=data_row, bin_num=8, gh=gh2, gh_sig=gh2_sig, chi_check=chi_check, chi_fit_num=20, ini_left=-0.06, ini_right=0.06);
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

        sprintf(inform, "%d, PDF. g1: %.6f, %.6f (%.6f), g2: %.6f, %.6f (%.6f).",rank, g1_t[i], gh1,gh1_sig,g2_t[i],gh2,gh2_sig);
        std::cout<<inform<<std::endl;

    }

    for(i=0;i<numprocs;i++)
    {
        if(rank == i)
        {
            std::cout<<rank<<std::endl;
            show_arr(result_sub_pdf, shear_ed - shear_st, result_col);
            show_arr(result_sub_mean, shear_ed - shear_st, result_col);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

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
        std::cout<<"AVERAGE: m & c"<<std::endl;
        show_arr(mean_mc, 2, 4);
        std::cout<<"PDF_SYM: m & c"<<std::endl;
        show_arr(pdf_mc, 2, 4);

        sprintf(result_path, "%s/result_%s.hdf5", parent_path, data_type);
        sprintf(set_name,"/mean_result");
        write_h5(result_path, set_name, result_all_mean, shear_num, result_col, true);
        sprintf(set_name,"/sym_result");
        write_h5(result_path, set_name, result_all_pdf, shear_num, result_col, false);
        sprintf(set_name,"/mean_mc");
        write_h5(result_path,set_name, mean_mc, 2, 4, false);
        sprintf(set_name,"/sym_mc");
        write_h5(result_path,set_name, pdf_mc, 2, 4, false);
        sprintf(set_name,"/chisq_g1");
        write_h5(result_path,set_name, total_chi_check_g1, shear_num, chi_check_num, false);
        sprintf(set_name,"/chisq_g2");
        write_h5(result_path,set_name, total_chi_check_g2, shear_num, chi_check_num, false);

        delete[] mean_mc;
        delete[] pdf_mc;
        delete[] fit_val;
        delete[] fit_err;
        delete[] result_all_mean;
        delete[] result_all_pdf;
        delete[] total_chi_check_g1;
        delete[] total_chi_check_g2;

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
    MPI_Finalize();
    return 0;
}