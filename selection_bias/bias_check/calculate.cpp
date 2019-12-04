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
    char set_name[30];

    int i,j,k;
    int shear_num, shear_st, shear_ed;
    int*shear_point;
    double *shears, *g1_t, *g2_t;
    double gh1, gh1_sig, gh2, gh2_sig;

    double chi_check[40];
    double *total_chi_check, *sub_chi_check;
    int *chi_send_count;

    int data_col, data_row;
    double *data, *mg1, *mg2,*mnu1, *mnu2;

    // for results
    int sub_num;
    double *result_all, *result_sub;
    int *send_count;

    // data shape
    strcpy(parent_path, argv[1]);
    data_row = atoi(argv[2])*10000;
    data_col = atoi(argv[3]);
    shear_num = 20;

    data = new double[data_row*data_col];
    mg1 = new double[data_row];
    mg2 = new double[data_row];
    mnu1 = new double[data_row];
    mnu2 = new double[data_row];

    shear_point = new int[numprocs]{};
    send_count = new int[numprocs]{};
    chi_send_count = new int[numprocs]{};

    // shear point distribution
    
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
        send_count[i] = shear_point[i]*8;
        chi_send_count[i] = shear_point[i]*160;
    }

    // the measured g1, sig1, g2, sig2 of each thread
    // measured from all source
    // will be sent to rank 0 stack into result_all
    result_sub = new double[(shear_ed - shear_st)*8]{};

    if(rank == 0)
    {
        show_arr(shear_point, 1, numprocs);
    }

    sub_chi_check = new double[(shear_ed - shear_st)*160];

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
        sprintf(data_path,"%s/data_%d_noise_free.hdf5",parent_path, i);
        read_h5(data_path, set_name, data);

        for(k=0;k<numprocs;k++)
        {
            if(k == rank)
            {
                std::cout<<rank<<" shear: "<<k<<" data:"<<data_path<<" g1: "<<g1_t[k]<<" g2: "<<g2_t[k]<<std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        // noise free
        for(k=0;k<data_row;k++)
        {
            mg1[k] = data[k*data_col];
            mg2[k] = data[k*data_col + 1];
            mnu1[k] = data[k*data_col + 2] + data[k*data_col + 3];
            mnu2[k] = data[k*data_col + 2] - data[k*data_col + 3];
        }
        try
        {
            //find_shear(mg=mg1, mnu=mnu1, data_num=data_row, bin_num=8, gh=gh1, gh_sig=gh1_sig, chi_check=chi_check, chi_fit_num=20, ini_left=-0.06, ini_right=0.06);
            find_shear(mg1, mnu1, data_row, 8, gh1, gh1_sig, chi_check, 20, 0, 100, -0.06,0.06,30);
            for(k=0;k<40;k++)
            {
                sub_chi_check[(i-shear_st)*160+k] = chi_check[k];
            }
        }
        catch(const char*img)
        {
            std::cout<<rank<<" "<<i<<" wrong, noise free g1"<<std::endl;
        }

        try
        {
            find_shear(mg2, mnu2, data_row, 8, gh2, gh2_sig, chi_check, 20, 0, 100, -0.06,0.06,30);
            //find_shear(mg=mg2, mnu=mnu2, data_num=data_row, bin_num=8, gh=gh2, gh_sig=gh2_sig, chi_check=chi_check, chi_fit_num=20, ini_left=-0.06, ini_right=0.06);
            for(k=0;k<40;k++)
            {
                sub_chi_check[(i-shear_st)*160+ 40 + k] = chi_check[k];
            }
        }
        catch(const char*img)
        {
            std::cout<<rank<<" "<<i<<" wrong, noise free g2"<<std::endl;
        } 

        result_sub[(i-shear_st)*8] = gh1;
        result_sub[(i-shear_st)*8 + 1] = gh1_sig;
        result_sub[(i-shear_st)*8 + 2] = gh2;
        result_sub[(i-shear_st)*8 + 3] = gh2_sig;
        
    
        // noisy
        for(k=0;k<data_row;k++)
        {
            mg1[k] = data[k*data_col + 4];
            mg2[k] = data[k*data_col + 5];
            mnu1[k] = data[k*data_col + 6] + data[k*data_col + 7];
            mnu2[k] = data[k*data_col + 6] - data[k*data_col + 7];
        }

        try
        {
            find_shear(mg1, mnu1, data_row, 8, gh1, gh1_sig, chi_check, 20, 0, 100, -0.06,0.06,30);
            //find_shear(mg=mg1, mnu=mnu1, data_num=data_row, bin_num=8, gh=gh1, gh_sig=gh1_sig, chi_check=chi_check, chi_fit_num=20, ini_left=-0.06, ini_right=0.06);
            for(k=0;k<40;k++)
            {
                sub_chi_check[(i-shear_st)*160+ 80 + k] = chi_check[k];
            }
        }
        catch(const char*img)
        {
            std::cout<<rank<<" "<<i<<" wrong, g1"<<std::endl;
        }
        try
        {
            find_shear(mg2, mnu2, data_row, 8, gh2, gh2_sig, chi_check, 20, 0, 100, -0.06,0.06,30);
            //find_shear(mg=mg2, mnu=mnu2, data_num=data_row, bin_num=8, gh=gh2, gh_sig=gh2_sig, chi_check=chi_check, chi_fit_num=20, ini_left=-0.06, ini_right=0.06);
            for(k=0;k<40;k++)
            {
                sub_chi_check[(i-shear_st)*160+ 120 + k] = chi_check[k];
            }
        }
        catch(const char*img)
        {
            std::cout<<rank<<" "<<i<<" wrong, g1"<<std::endl;
        } 

        //find_shear(mg1, mnu1, data_row, 8, gh1, gh1_sig, chi_check, 20);
        //find_shear(mg2, mnu2, data_row, 8, gh2, gh2_sig, chi_check, 20);

        result_sub[(i-shear_st)*8 + 4] = gh1;
        result_sub[(i-shear_st)*8 + 5] = gh1_sig;
        result_sub[(i-shear_st)*8 + 6] = gh2;
        result_sub[(i-shear_st)*8 + 7] = gh2_sig;

    }
    for(i=0;i<numprocs;i++)
    {
        if(rank == i)
        {
            std::cout<<rank<<std::endl;
            show_arr(result_sub, shear_ed - shear_st, 8);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if(rank == 0)
    {
        result_all = new double[shear_num*8]{};
        total_chi_check = new double[shear_num*160]{};
    }

    my_Gatherv(result_sub, send_count, result_all, numprocs, rank, 0);
    my_Gatherv(sub_chi_check, chi_send_count, total_chi_check, numprocs, rank, 0);


    if(rank == 0)
    {   
        double *mc = new double[4];
        double *mc_all = new double[16];
        double *fit_val = new double[shear_num];
        double *fit_err = new double[shear_num];
        
        // noise free g1
        for(k=0;k<shear_num;k++)
        {
            fit_val[k] = result_all[k*8];
            fit_err[k] = result_all[k*8 + 1];
        }
        poly_fit_1d(g1_t, fit_val, fit_err, shear_num, mc, 1);
        for(i=0;i<4;i++)
        {
            mc_all[i] = mc[i];
        }

        // noise free g2
        for(k=0;k<shear_num;k++)
        {
            fit_val[k] = result_all[k*8 + 2];
            fit_err[k] = result_all[k*8 + 3];
        }
        poly_fit_1d(g2_t, fit_val, fit_err, shear_num, mc, 1);
        for(i=0;i<4;i++)
        {
            mc_all[i+4] = mc[i];
        }
        // noisy g1
        for(k=0;k<shear_num;k++)
        {
            fit_val[k] = result_all[k*8 + 4];
            fit_err[k] = result_all[k*8 + 4 + 1];
        }
        poly_fit_1d(g1_t, fit_val, fit_err, shear_num, mc, 1);
        for(i=0;i<4;i++)
        {
            mc_all[i+8] = mc[i];
        }
        // noisy g2
        for(k=0;k<shear_num;k++)
        {
            fit_val[k] = result_all[k*8 + 4 + 2];
            fit_err[k] = result_all[k*8 + 4 + 3];
        }
        poly_fit_1d(g2_t, fit_val, fit_err, shear_num, mc, 1);
        for(i=0;i<4;i++)
        {
            mc_all[i+12] = mc[i];
        }

        sprintf(result_path, "%s/result.hdf5", parent_path);
        sprintf(set_name,"/data");
        write_h5(result_path, set_name, result_all, shear_num, 8, true);
        sprintf(set_name,"/mc");
        write_h5(result_path,set_name, mc_all, 2, 8, false);
        sprintf(set_name,"/chisq");
        write_h5(result_path,set_name, total_chi_check, shear_num, 160, false);
        delete[] mc;
        delete[] mc_all;
        delete[] fit_val;
        delete[] fit_err;
        delete[] result_all;
    }
    
    delete[] send_count;
    delete[] result_sub;
    delete[] data;
    delete[] mnu2;
    delete[] mnu1;
    delete[] mg2;
    delete[] mg1;
    delete[] shear_point;
    MPI_Finalize();
    return 0;
}