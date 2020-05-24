#include<FQlib.h>
#include<hk_iolib.h>
#include<hk_mpi.h>

#define MY_FLOAT double

void task_alloc(const int total_task_num, const int division_num, const int my_part_id, int &my_st_id, int &my_ed_id, int *task_count)
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

    char parent_path[300], data_path[300], set_name[30], data_name[40], result_path[300];
    char inform[300], *temp_inform[50];
    char data_nm[30];
    char time_now[50],time_now_1[50];
    double st1, st2;

    st1 = clock();

    MY_FLOAT *mg[5];
    char *mg_name[5];


    MY_FLOAT *mg1_bin, *mg2_bin;
    int *num_in_bin, bin_num;

    MY_FLOAT *g1t, *g2t;
    MY_FLOAT *shear_m, *shear_result, *shear_nu, *shear_result_nu;
    int result_col_nu;

    int grid_row, grid_col, total_grid, my_grid_st, my_grid_ed;

    int *shear_count, *send_count, *send_count_nu;   

    MY_FLOAT *chi_check, chi_min;
    int chi_fit_num;

    int i,j,k;
    int shear_num, shear_id;
    int data_row, data_col;

    MY_FLOAT g1, g1_sig, g2, g2_sig;
    int mg1_idx, mg2_idx, mn_idx, mu_idx, mv_idx;

    MY_FLOAT *mg_temp, *mnu_temp;
    MY_FLOAT g1_temp, g1_temp_sig, g2_temp, g2_temp_sig;
    MY_FLOAT gN, gU, gN_sig, gU_sig,delta_g, chisq_N, chisq_U;
    int my_shear_st, my_shear_ed;
    int iters;

    int ratio;
    int add_time, add_scale;
    
    add_time = 100;//atoi(argv[5]);
    add_scale = 100;//atoi(argv[6]);
    bin_num = 20;
    data_col = 5;
    result_col_nu = 6;
    iters = 4;
    chi_fit_num = 20;

    mg1_idx = 0;
    mg2_idx = 1;
    mn_idx = 2;
    mu_idx = 3;
    mv_idx = 4;

    strcpy(parent_path, argv[1]);
    strcpy(data_nm, argv[2]);
    data_row = atoi(argv[3])*10000;
    shear_num = atoi(argv[4]);
    

    chi_check = new MY_FLOAT[2*chi_fit_num];

    for(i=0; i<iters*2; i++)
    {
        temp_inform[i] = new char[200];
    }
    
    // std::cout<<"Get task"<<std::endl;
    // show_arr(grid_count,1,numprocs);
    

    num_in_bin = new int[bin_num]{};
    mg1_bin = new MY_FLOAT[bin_num+1]{};
 

    for(i=0;i<5;i++)
    {
        mg[i] = new MY_FLOAT[data_row];
        mg_name[i] = new char[50];
    }
    sprintf(mg_name[0], "/mg1");
    sprintf(mg_name[1], "/mg2");
    sprintf(mg_name[2], "/mn");
    sprintf(mg_name[3], "/mu");
    sprintf(mg_name[4], "/mv");



    sprintf(result_path, "%s/shear_result_pdf_iter_%s.hdf5",parent_path, data_nm);


    sprintf(data_path, "%s/shear.hdf5", parent_path);
    
    g1t = new MY_FLOAT[shear_num];
    g2t = new MY_FLOAT[shear_num];

    sprintf(set_name,"/g1");
    read_h5(data_path, set_name, g1t);

    sprintf(set_name,"/g2");
    read_h5(data_path, set_name, g2t);
    


    MY_FLOAT *mc = new MY_FLOAT[4];
    MY_FLOAT *result_mc = new MY_FLOAT[8];
    MY_FLOAT *result_mc_iter = new MY_FLOAT[8*iters];

    MY_FLOAT *fit_val = new MY_FLOAT[shear_num];
    MY_FLOAT *fit_err = new MY_FLOAT[shear_num];

    shear_m = new MY_FLOAT[4];
    shear_nu = new MY_FLOAT[iters*result_col_nu];
    if(rank == 0)
    {
        shear_result = new MY_FLOAT[shear_num*4];
        shear_result_nu = new MY_FLOAT[iters*result_col_nu*shear_num];
    }
    shear_count = new int[shear_num];
    send_count = new int[shear_num];
    send_count_nu = new int[shear_num];

    
    task_alloc(shear_num, numprocs, rank, my_shear_st, my_shear_ed, shear_count);

    for(i=0;i<shear_num;i++)
    {
        send_count[i] = 4;
        send_count_nu[i] = 6*iters;
    }
    if(rank == 0)
    {   
        show_arr(g1t, 1, shear_num);
        show_arr(g2t, 1, shear_num);
        std::cout<<"Bin_num "<<bin_num<<" data:"<<parent_path<<" "<<std::endl;
    }
    
    // for(i=0;i<numprocs;i++)
    // {
    //     if(i==rank){std::cout<<my_grid_st<<" "<<my_grid_ed<<std::endl;}
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    for(shear_id=my_shear_st; shear_id<my_shear_ed; shear_id++)
    {   
        // sprintf(data_name, data_nm, shear_id);
        sprintf(data_path, "%s/data_%s_%d.hdf5",parent_path, data_nm, shear_id);
        for(j=0;j<5;j++)
        {
            read_h5(data_path, mg_name[j], mg[j]);
        }
        std::cout<<"Read data "<<data_path<<std::endl;
        
        //////// average  ////////
        find_shear_mean(mg[0], mg[2], data_row, g1, g1_sig, 1000, 1);
        find_shear_mean(mg[1], mg[2], data_row, g2, g2_sig, 1000, 1);

        sprintf(temp_inform[0], "Ave g1(true %9.6f): %9.6f (%9.6f)", g1t[shear_id], g1, g1_sig);
        sprintf(temp_inform[1], "Ave g2(true %9.6f): %9.6f (%9.6f)", g2t[shear_id], g2, g2_sig);
        
        shear_m[0] = g1;
        shear_m[1] = g1_sig;
        shear_m[2] = g2;
        shear_m[3] = g2_sig;
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        for(i=0;i<numprocs;i++)
        {
            if(i==rank)
            {
                std::cout<<temp_inform[0]<<std::endl<<temp_inform[1]<<std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        my_Gatherv(shear_m, send_count, shear_result, numprocs, rank);
        if(rank == 0)
        {
            sprintf(set_name, "/average");
            write_h5(result_path, set_name, shear_result, shear_num, 4, true);
            // fitting
            for(i=0;i<shear_num;i++)
            {
                fit_val[i] = shear_result[i*4];
                fit_err[i] = shear_result[i*4+1];
            }
            poly_fit_1d(g1t, fit_val, fit_err, shear_num, mc, 1);
            result_mc[0] = mc[2]-1;// m
            result_mc[1] = mc[3];// m_sig
            result_mc[2] = mc[0];// c
            result_mc[3] = mc[1];// c_sig
            for(i=0;i<shear_num;i++)
            {
                fit_val[i] = shear_result[i*4+2];
                fit_err[i] = shear_result[i*4+3];
            }
            poly_fit_1d(g2t, fit_val, fit_err, shear_num, mc, 1);
            result_mc[4] = mc[2]-1;// m
            result_mc[5] = mc[3];// m_sig
            result_mc[6] = mc[0];// c
            result_mc[7] = mc[1];// c_sig
            sprintf(set_name, "/mean_mc");
            write_h5(result_path, set_name, result_mc, 2, 4, false);
            std::cout<<"Mean: m & c"<<std::endl;
            show_arr(result_mc,2,4);
        }
        MPI_Barrier(MPI_COMM_WORLD);



        set_bin(mg[0], data_row, bin_num, mg1_bin, 100, 0);
        // show_arr(mg1_bin,1,bin_num+1);

        find_shear_iter(mg[0], mg[2], mg[3], data_row, bin_num, mg1_bin, 1, iters, shear_nu, chi_fit_num, -0.1, 0.1, 40);

        MPI_Barrier(MPI_COMM_WORLD);   
        my_Gatherv(shear_nu, send_count_nu, shear_result_nu, numprocs, rank);
        if(rank == 0)
        {
            sprintf(set_name, "/new_PDF/g1");
            write_h5(result_path, set_name, shear_result_nu, shear_num, result_col_nu*iters, false);
                        // fitting
            for(i=0;i<iters;i++)
            {
                for(j=0;j<shear_num;j++)
                {
                    fit_val[j] = shear_result_nu[j*result_col_nu*iters + i*result_col_nu];
                    fit_err[j] = shear_result_nu[j*result_col_nu*iters + 1 + i*result_col_nu];
                }
                poly_fit_1d(g1t, fit_val, fit_err, shear_num, mc, 1);
                result_mc_iter[i*8] = mc[2]-1;// m
                result_mc_iter[i*8+1] = mc[3];// m_sig
                result_mc_iter[i*8+2] = mc[0];// c
                result_mc_iter[i*8+3] = mc[1];// c_sig
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        find_shear_iter(mg[1], mg[2], mg[3], data_row, bin_num, mg1_bin, 2, iters, shear_nu, chi_fit_num, -0.1, 0.1, 40);

        MPI_Barrier(MPI_COMM_WORLD);   
        my_Gatherv(shear_nu, send_count_nu, shear_result_nu, numprocs, rank);
        if(rank == 0)
        {
            sprintf(set_name, "/new_PDF/g2");
            write_h5(result_path, set_name, shear_result_nu, shear_num, result_col_nu*iters, false);

            // fitting
            for(i=0;i<iters;i++)
            {
                for(j=0;j<shear_num;j++)
                {
                    fit_val[j] = shear_result_nu[j*result_col_nu*iters + i*result_col_nu];
                    fit_err[j] = shear_result_nu[j*result_col_nu*iters + 1 + i*result_col_nu];
                }
                poly_fit_1d(g2t, fit_val, fit_err, shear_num, mc, 1);
                result_mc_iter[i*8 + 4] = mc[2]-1;// m
                result_mc_iter[i*8 + 5] = mc[3];// m_sig
                result_mc_iter[i*8 + 6] = mc[0];// c
                result_mc_iter[i*8 + 7] = mc[1];// c_sig
                
                sprintf(set_name, "/new_sym_mc_iter_%d", i);
                write_h5(result_path, set_name, &result_mc_iter[i*8], 2, 4, false);

            }

        }
        MPI_Barrier(MPI_COMM_WORLD);

    }
    st2 = clock();
    if(rank == 0)
    {           
        for(i=0;i<iters;i++)
        {   std::cout<<"PDF_Iter "<<i<<" : m & c"<<std::endl;
            show_arr(&result_mc_iter[i*8], 2, 4);
        }              
        std::cout<<result_path<<std::endl;       
        std::cout<<(st2-st1)/CLOCKS_PER_SEC<<" SEC"<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    

    MPI_Finalize();
    return 0;
}