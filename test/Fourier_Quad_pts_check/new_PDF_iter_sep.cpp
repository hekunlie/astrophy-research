#include<FQlib.h>
#include<hk_iolib.h>
#include<hk_mpi.h>

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

    char parent_path[200], data_path[200], set_name[30], data_name[40], result_path[200];
    char inform[300], *temp_inform[50];
    char data_nm[30];
    char time_now[50],time_now_1[50];

    double *data, *mg1, *mg2, *mn, *mu, *mnu1, *mnu2, *mg;
    double *mg1_w, *mg2_w, *mn_w, *weight;
    
    double *mg1_bin, *mg2_bin;
    int *num_in_bin, bin_num;

    double *g1t, *g2t;
    double *shear_m, *shear_result, *shear_nu, *shear_result_nu;
    int result_col_nu;

    int grid_row, grid_col, total_grid, my_grid_st, my_grid_ed;

    int *shear_count, *send_count, *send_count_nu;   

    double *check1;
    double *check2;    

    int i, j, k, tag;
    int shear_num, shear_id;
    int data_row = 20000000, sub_row, data_col;
    double g1, g1_sig, g2, g2_sig;

    double *check_temp;
    double *mg_temp, *mg_temp_sub[50], *mnu_temp, *mu_temp_sub[50];
    double g1_temp, g1_temp_sig, g2_temp, g2_temp_sig;
    double gN, gU, gN_sig, gU_sig,delta_g, chisq_N, chisq_U;
    int my_shear_st, my_shear_ed;
    int iters;
    
    double *scale_vals;
    double *scale_bins;
    int scale_bin_num, *scale_num_counts, *scale_bin_label, temp;

    int ratio;
    int add_time, add_scale;
    add_time = 20;//atoi(argv[5]);
    add_scale = 20;//atoi(argv[6]);
    bin_num = 8;
    data_col = 7;
    
    iters = 4;

    scale_bin_num = 5;

    result_col_nu = 3 + scale_bin_num*3;

    strcpy(parent_path, argv[1]);
    strcpy(data_nm, argv[2]);
    
    ratio = atoi(argv[3]);
    
    sub_row = data_row/ratio;

    grid_row = 15;

    check1 = new double[2*grid_row];
    check2 = new double[2*grid_row];
    check_temp = new double[2*grid_row];

    for(i=0;i<iters*2;i++)
    {
        temp_inform[i] = new char[200];
    }
    
    // std::cout<<"Get task"<<std::endl;
    // show_arr(grid_count,1,numprocs);
    

    num_in_bin = new int[bin_num]{};
    mg1_bin = new double[bin_num+1]{};
    mg2_bin = new double[bin_num+1]{};    

    scale_bins = new double[scale_bin_num+1]{0, 3, 5, 7, 10, 100000};
    scale_vals = new double[sub_row];
    scale_bin_label = new int[sub_row];
    scale_num_counts = new int[scale_bin_num];

    data = new double[data_row*data_col];
    mg = new double[sub_row];
    mg1 = new double[sub_row];
    mg2 = new double[sub_row];
    mn = new double[sub_row];

    mg_temp = new double[sub_row];
    mnu_temp = new double[sub_row];

    mu = new double[sub_row];
    mnu1 = new double[sub_row];
    mnu2 = new double[sub_row];

    // weight = new double[sub_row];
    
    sprintf(result_path, "%s/new_pdf_%d.hdf5",parent_path, scale_bin_num);

    sprintf(data_path, "%s/shear.hdf5", parent_path);
    
    sprintf(set_name,"/g1");
    read_h5_datasize(data_path, set_name,shear_num);
    g1t = new double[shear_num];
    g2t = new double[shear_num];
    read_h5(data_path, set_name, g1t);

    sprintf(set_name,"/g2");
    read_h5_datasize(data_path, set_name,shear_num);
    read_h5(data_path, set_name, g2t);
    

    shear_m = new double[6];
    shear_nu = new double[iters*result_col_nu];
    if(rank == 0)
    {
        shear_result = new double[shear_num*6];
        shear_result_nu = new double[iters*result_col_nu*shear_num];
    }
    shear_count = new int[shear_num];
    send_count = new int[shear_num];
    send_count_nu = new int[shear_num];

    
    task_alloc(shear_num, numprocs, rank, my_shear_st, my_shear_ed, shear_count);

    for(i=0;i<shear_num;i++)
    {
        send_count[i] = 6;
        send_count_nu[i] = result_col_nu*iters;
    }
    if(rank == 0)
    {   
        show_arr(g1t, 1, shear_num);
        show_arr(g2t, 1, shear_num);
        std::cout<<rank<<" "<<sub_row<<" "<< my_shear_st<<" "<<my_shear_ed<<std::endl;
        show_arr(scale_bins, 1, scale_bin_num+1);
    }
    
    // for(i=0;i<numprocs;i++)
    // {
    //     if(i==rank){std::cout<<my_grid_st<<" "<<my_grid_ed<<std::endl;}
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    for(shear_id=my_shear_st; shear_id<my_shear_ed; shear_id++)
    {   
        // sprintf(data_name, data_nm, shear_id);
        sprintf(data_path, "%s/data_%d.hdf5",parent_path, shear_id);
        sprintf(set_name,"/data");
        read_h5(data_path, set_name, data);

        sprintf(data_path, "%s/sex2_1.5/flux2_ex1_%d.hdf5",parent_path, shear_id);
        sprintf(set_name,"/data");
        read_h5(data_path, set_name, scale_vals);

        for(i=0;i<sub_row;i++)
        {   
            //weight[i] = data[i*data_col + 7]*data[i*data_col + 7];

            mg1[i] = data[i*data_col + 2];
            mg2[i] = data[i*data_col + 3];
            mn[i] = data[i*data_col + 4];
            mu[i] = data[i*data_col + 5];
            mnu1[i] = mn[i] + mu[i];
            mnu2[i] = mn[i] - mu[i];

            scale_vals[i] = scale_vals[i]/48./60.;
            histogram_s(scale_vals[i], scale_bins, scale_bin_num, k);
            scale_bin_label[i] = k;
            if(rank == 0 and i <= 20){std::cout<<scale_bin_label[i]<<std::endl;}
            // mg1_w[i] = mg1[i]/data[i*data_col]/data[i*data_col];
            // mg2_w[i] = mg2[i]/data[i*data_col]/data[i*data_col];
            // mn_w[i] = mn[i]/data[i*data_col]/data[i*data_col];

        }
        histogram(scale_vals, scale_bins, scale_num_counts, sub_row, scale_bin_num);

        std::cout<<"Read data "<<data_path<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

        for(i=0;i<numprocs;i++)
        {
            if(i==rank)
            {
                show_arr(scale_num_counts, 1, scale_bin_num);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);  


        //////// PDF ////////
        // set_bin(mg1, sub_row, mg1_bin, bin_num, 100);
        // set_bin(mg2, sub_row, mg2_bin, bin_num, 100);

        find_shear(mg1, mnu1, sub_row, 10, g1, g1_sig, chisq_N, check1, grid_row, 0,100,-0.1, 0.1, 40);
        find_shear(mg2, mnu2, sub_row, 10, g2, g2_sig, chisq_U, check2, grid_row, 0,100,-0.1, 0.1, 40);
        
        shear_m[0] = g1;
        shear_m[1] = g1_sig;
        shear_m[2] = chisq_N;
        shear_m[3] = g2;
        shear_m[4] = g2_sig;
        shear_m[5] = chisq_U;

        MPI_Barrier(MPI_COMM_WORLD);
        for(i=0;i<numprocs;i++)
        {
            if(i==rank)
            {
                sprintf(inform, "PDF g1(true %9.6f): %9.6f (%9.6f), chisq: %.4f", g1t[shear_id], g1, g1_sig, chisq_N);
                std::cout<<inform<<std::endl;            
                sprintf(inform, "PDF g2(true %9.6f): %9.6f (%9.6f), chisq: %.4f", g2t[shear_id], g2, g2_sig, chisq_U);
                std::cout<<inform<<std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);   
        my_Gatherv(shear_m, send_count, shear_result, numprocs, rank);
        if(rank == 0)
        {
            sprintf(set_name, "/PDF");
            write_h5(result_path, set_name, shear_result, shear_num, 6, true);
        }
        MPI_Barrier(MPI_COMM_WORLD);


        //////// new PDF ////////
        for(i=0; i<scale_bin_num; i++)
        {
            mg_temp_sub[i] = new double[scale_num_counts[i]];
            mu_temp_sub[i] = new double[scale_num_counts[i]];
        }


        gN = g1;
        gU = g1;


        for(i=0; i<iters; i++)
        {   
            // fix g_N
            for(k=0; k<sub_row; k++){mg_temp[k] = mg1[k] - gN*mn[k];}
            
            // update g_U_i
            for(j=0; j<scale_bin_num; j++)
            {   
                tag = 0;
                for(k=0; k<sub_row; k++)
                {
                    if(scale_bin_label[k] == j)
                    {
                        mg_temp_sub[j][tag] = mg_temp[k];
                        mu_temp_sub[j][tag] = mu[k];
                        tag ++;
                    }
                }
                find_shear(mg_temp_sub[j], mu_temp_sub[j], scale_num_counts[j], 10, gU, gU_sig, chisq_U, check_temp, grid_row, 0, 100,-0.07, 0.07, 60);
                
                shear_nu[i*result_col_nu + 3 + j*3] = gU;
                shear_nu[i*result_col_nu + 3 + j*3 + 1] = gU_sig;
                shear_nu[i*result_col_nu + 3 + j*3 + 2] = chisq_U;

                // update g_U_i
                for(k=0; k<sub_row; k++)
                {
                    if(scale_bin_label[k] == j)
                    {
                        mg_temp[k] = mg1[k] - gU*mu[k];
                    }
                }
            }
            // update g_N
            find_shear(mg_temp, mn, sub_row, 10, gN, gN_sig, chisq_N, check_temp, grid_row, 0, 100,-0.07, 0.07, 60);
            
            sprintf(temp_inform[i], "Iter %d. (true %9.6f) g1_N: %9.6f (%9.6f). chi^2: %.4f", i, g1t[shear_id], gN, gN_sig, chisq_N);

            shear_nu[i*result_col_nu] = gN;
            shear_nu[i*result_col_nu + 1] = gN_sig;
            shear_nu[i*result_col_nu + 2] = chisq_N;
        }
        
        MPI_Barrier(MPI_COMM_WORLD);   
        my_Gatherv(shear_nu, send_count_nu, shear_result_nu, numprocs, rank);
        if(rank == 0)
        {
            sprintf(set_name, "/new_PDF/g1");
            write_h5(result_path, set_name, shear_result_nu, shear_num, result_col_nu*iters, false);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        gN = g2;
        gU = g2;

        for(i=0; i<iters; i++)
        {   
            // fix g_N
            for(k=0; k<sub_row; k++){mg_temp[k] = mg2[k] - gN*mn[k];}
            
            // update g_U_i
            for(j=0; j<scale_bin_num; j++)
            {   
                tag = 0;
                for(k=0; k<sub_row; k++)
                {
                    if(scale_bin_label[k] == j)
                    {
                        mg_temp_sub[j][tag] = mg_temp[k];
                        mu_temp_sub[j][tag] = mu[k];
                        tag ++;
                    }
                }
                find_shear(mg_temp_sub[j], mu_temp_sub[j], scale_num_counts[j], 10, gU, gU_sig, chisq_U, check_temp, grid_row, 0, 100,-0.07, 0.07, 60);
                
                shear_nu[i*result_col_nu + 3 + j*3] = gU;
                shear_nu[i*result_col_nu + 3 + j*3 + 1] = gU_sig;
                shear_nu[i*result_col_nu + 3 + j*3 + 2] = chisq_U;

                // update g_U_i
                for(k=0; k<sub_row; k++)
                {
                    if(scale_bin_label[k] == j)
                    {
                        mg_temp[k] = mg2[k] - gU*mu[k];
                    }
                }
            }
            // update g_N
            find_shear(mg_temp, mn, sub_row, 10, gN, gN_sig, chisq_N, check_temp, grid_row, 0, 100,-0.07, 0.07, 60);
            
            sprintf(temp_inform[iters+i], "Iter %d. (true %9.6f) g2_N: %9.6f (%9.6f). chi^2: %.4f", i, g2t[shear_id], gN, gN_sig, chisq_N);

            shear_nu[i*result_col_nu] = gN;
            shear_nu[i*result_col_nu + 1] = gN_sig;
            shear_nu[i*result_col_nu + 2] = chisq_N;
        }
        
        MPI_Barrier(MPI_COMM_WORLD);   
        my_Gatherv(shear_nu, send_count_nu, shear_result_nu, numprocs, rank);
        if(rank == 0)
        {
            sprintf(set_name, "/new_PDF/g2");
            write_h5(result_path, set_name, shear_result_nu, shear_num, result_col_nu*iters, false);
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

        for(i=0; i<scale_bin_num; i++)
        {
            delete[] mg_temp_sub[i];
            delete[] mu_temp_sub[i];
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);


    MPI_Finalize();
    return 0;
}