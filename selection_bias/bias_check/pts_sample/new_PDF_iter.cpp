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

    MY_FLOAT *data, *mg1, *mg2, *mn, *mu, *mnu1, *mnu2, *mg;
    MY_FLOAT *mg1_w, *mg2_w, *mn_w, *weight;
    
    MY_FLOAT *mg1_bin, *mg2_bin;
    int *num_in_bin, bin_num;

    MY_FLOAT *g1t, *g2t;
    MY_FLOAT *shear_m, *shear_result, *shear_nu, *shear_result_nu;
    int result_col_nu;

    int grid_row, grid_col, total_grid, my_grid_st, my_grid_ed;

    int *shear_count, *send_count, *send_count_nu;   

    MY_FLOAT *check1;
    MY_FLOAT *check2;    

    int i,j,k;
    int shear_num, shear_id;
    int data_row, sub_row, data_col;
    MY_FLOAT g1, g1_sig, g2, g2_sig;

    MY_FLOAT *check_temp;
    MY_FLOAT *mg_temp, *mnu_temp;
    MY_FLOAT g1_temp, g1_temp_sig, g2_temp, g2_temp_sig;
    MY_FLOAT gN, gU, gN_sig, gU_sig,delta_g, chisq_N, chisq_U;
    int my_shear_st, my_shear_ed;
    int iters;

    int ratio;
    int add_time, add_scale;
    
    add_time = 100;//atoi(argv[5]);
    add_scale = 100;//atoi(argv[6]);
    bin_num = 8;
    data_col = 5;
    result_col_nu = 6;
    iters = 3;

    strcpy(parent_path, argv[1]);
    strcpy(data_nm, argv[2]);
    data_row = atoi(argv[3])*10000;
    
    ratio = 1;
    
    sub_row = data_row/ratio;

    grid_row = 15;

    check1 = new MY_FLOAT[2*grid_row];
    check2 = new MY_FLOAT[2*grid_row];
    check_temp = new MY_FLOAT[2*grid_row];

    for(i=0;i<iters*2;i++)
    {
        temp_inform[i] = new char[200];
    }
    
    // std::cout<<"Get task"<<std::endl;
    // show_arr(grid_count,1,numprocs);
    

    num_in_bin = new int[bin_num]{};
    mg1_bin = new MY_FLOAT[bin_num+1]{};
    mg2_bin = new MY_FLOAT[bin_num+1]{};    


    data = new MY_FLOAT[data_row*data_col];
    mg = new MY_FLOAT[sub_row];
    mg1 = new MY_FLOAT[sub_row];
    mg2 = new MY_FLOAT[sub_row];
    mn = new MY_FLOAT[sub_row];

    mg_temp = new MY_FLOAT[sub_row];
    mnu_temp = new MY_FLOAT[sub_row];

    mu = new MY_FLOAT[sub_row];
    mnu1 = new MY_FLOAT[sub_row];
    mnu2 = new MY_FLOAT[sub_row];


    sprintf(result_path, "%s/shear_result_pdf_iter_%s.hdf5",parent_path, data_nm);

    sprintf(data_path, "%s/shear.hdf5", parent_path);
    
    sprintf(set_name,"/g1");
    read_h5_datasize(data_path, set_name,shear_num);
    g1t = new MY_FLOAT[shear_num];
    g2t = new MY_FLOAT[shear_num];
    read_h5(data_path, set_name, g1t);

    sprintf(set_name,"/g2");
    read_h5_datasize(data_path, set_name,shear_num);
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
        std::cout<<rank<<" "<<sub_row<<" "<< my_shear_st<<" "<<my_shear_ed<<std::endl;
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
        sprintf(set_name,"/data");
        read_h5(data_path, set_name, data);

        for(i=0;i<sub_row;i++)
        {   
            //weight[i] = data[i*data_col + 7]*data[i*data_col + 7];

            mg1[i] = data[i*data_col + 0];
            mg2[i] = data[i*data_col + 1];
            mn[i] = data[i*data_col + 2];
            mu[i] = data[i*data_col + 3];
            mnu1[i] = mn[i] + mu[i];
            mnu2[i] = mn[i] - mu[i];

            // mg1_w[i] = mg1[i]/data[i*data_col]/data[i*data_col];
            // mg2_w[i] = mg2[i]/data[i*data_col]/data[i*data_col];
            // mn_w[i] = mn[i]/data[i*data_col]/data[i*data_col];

        }
        std::cout<<"Read data "<<data_path<<std::endl;
        
        //////// average  ////////
        find_shear_mean(mg1, mn, sub_row, g1, g1_sig, add_time, add_scale);
        find_shear_mean(mg2, mn, sub_row, g2, g2_sig, add_time, add_scale);

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

        }
        MPI_Barrier(MPI_COMM_WORLD);


        //////// PDF ////////
        set_bin(mg1, sub_row, mg1_bin, bin_num, 100);
        set_bin(mg2, sub_row, mg2_bin, bin_num, 100);

        find_shear(mg1, mnu1, sub_row, 10, g1, g1_sig, chisq_N, check1, grid_row, 0,100,-0.1, 0.1, 40);
        find_shear(mg2, mnu2, sub_row, 10, g2, g2_sig, chisq_U, check2, grid_row, 0,100,-0.1, 0.1, 40);
        
        shear_m[0] = g1;
        shear_m[1] = g1_sig;
        shear_m[2] = g2;
        shear_m[3] = g2_sig;

        MPI_Barrier(MPI_COMM_WORLD);
        // for(i=0;i<numprocs;i++)
        // {
        //     if(i==rank)
        //     {
        //         sprintf(inform, "PDF g1(true %9.6f): %9.6f (%9.6f), chisq: %.4f", g1t[shear_id], g1, g1_sig, chisq_N);
        //         std::cout<<inform<<std::endl;            
        //         sprintf(inform, "PDF g2(true %9.6f): %9.6f (%9.6f), chisq: %.4f", g2t[shear_id], g2, g2_sig, chisq_U);
        //         std::cout<<inform<<std::endl;
        //     }
        //     MPI_Barrier(MPI_COMM_WORLD);
        // }
        // MPI_Barrier(MPI_COMM_WORLD);   
        // my_Gatherv(shear_m, send_count, shear_result, numprocs, rank);
        // if(rank == 0)
        // {
        //     sprintf(set_name, "/PDF");
        //     write_h5(result_path, set_name, shear_result, shear_num, 4, false);
        //     // fitting
        //     for(i=0;i<shear_num;i++)
        //     {
        //         fit_val[i] = shear_result[i*4];
        //         fit_err[i] = shear_result[i*4+1];
        //     }
        //     poly_fit_1d(g1t, fit_val, fit_err, shear_num, mc, 1);
        //     result_mc[0] = mc[2]-1;// m
        //     result_mc[1] = mc[3];// m_sig
        //     result_mc[2] = mc[0];// c
        //     result_mc[3] = mc[1];// c_sig
        //     for(i=0;i<shear_num;i++)
        //     {
        //         fit_val[i] = shear_result[i*4+2];
        //         fit_err[i] = shear_result[i*4+3];
        //     }
        //     poly_fit_1d(g2t, fit_val, fit_err, shear_num, mc, 1);
        //     result_mc[4] = mc[2]-1;// m
        //     result_mc[5] = mc[3];// m_sig
        //     result_mc[6] = mc[0];// c
        //     result_mc[7] = mc[1];// c_sig
        //     sprintf(set_name, "/sym_mc");
        //     write_h5(result_path, set_name, result_mc, 2, 4, false);

        // }
        // MPI_Barrier(MPI_COMM_WORLD);
        


        //////// new PDF ////////
        delta_g = 1;
        gN = g1;
        gU = g1;

        for(j=0;j<iters;j++)
        {   
            // fix g_N, find the g_U
            for(i=0;i<sub_row;i++){mg_temp[i] = mg1[i] - gN*mn[i];}
            find_shear(mg_temp, mu, sub_row, 10, gU, gU_sig, chisq_U, check_temp, grid_row, 0, 100,-0.07, 0.07, 60);
            // fix g_U, find the g_N
            for(i=0;i<sub_row;i++){mg_temp[i] = mg1[i] - gU*mu[i];}
            find_shear(mg_temp, mn, sub_row, 10, gN, gN_sig, chisq_N, check_temp, grid_row, 0, 100,-0.07, 0.07, 60);
            
            sprintf(temp_inform[j], "Iter %d. (true %9.6f) g1_N: %9.6f (%9.6f). g1_U: %9.6f (%9.6f). chi^2: %.4f, %.4f", j, g1t[shear_id], gN, gN_sig, gU, gU_sig, chisq_U, chisq_N);

            shear_nu[j*result_col_nu] = gN;
            shear_nu[j*result_col_nu + 1] = gN_sig;
            shear_nu[j*result_col_nu + 2] = gU;
            shear_nu[j*result_col_nu + 3] = gU_sig;
            shear_nu[j*result_col_nu + 4] = chisq_N;
            shear_nu[j*result_col_nu + 5] = chisq_U;
        }
        
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
        
        gN = g2;
        gU = g2;
        for(j=0;j<iters;j++)
        {   
            // fix g_N, find the g_U
            for(i=0;i<sub_row;i++){mg_temp[i] = mg2[i] - gN*mn[i];}
            find_shear(mg_temp, mu, sub_row, 10, gU, gU_sig, chisq_U, check_temp, grid_row, 0, 100,-0.07, 0.07, 60);
            // fix g_U, find the g_N
            for(i=0;i<sub_row;i++){mg_temp[i] = mg2[i] - gU*mu[i];}
            find_shear(mg_temp, mn, sub_row, 10, gN, gN_sig, chisq_N, check_temp, grid_row, 0, 100,-0.07, 0.07, 60);
            
            sprintf(temp_inform[j+iters], "Iter %d. (true %9.6f) g2_N: %9.6f (%9.6f). g2_U: %9.6f (%9.6f). chi^2: %.4f, %.4f", j, g2t[shear_id], gN, gN_sig, gU, gU_sig, chisq_U, chisq_N);
            
            shear_nu[j*result_col_nu] = gN;
            shear_nu[j*result_col_nu + 1] = gN_sig;
            shear_nu[j*result_col_nu + 2] = gU;
            shear_nu[j*result_col_nu + 3] = gU_sig;
            shear_nu[j*result_col_nu + 4] = chisq_N;
            shear_nu[j*result_col_nu + 5] = chisq_U;

        }
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
                show_arr(&result_mc_iter[i*8], 2, 4);
            }

        }
        MPI_Barrier(MPI_COMM_WORLD);

        for(i=0;i<numprocs;i++)
        {
            if(i==rank)
            {   
                for(j=0;j<2*iters;j++)
                {
                    std::cout<<temp_inform[j]<<std::endl;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        
      

    }

    MPI_Barrier(MPI_COMM_WORLD);


    MPI_Finalize();
    return 0;
}