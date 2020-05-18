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
    char inform[300];
    char data_nm[30];
    char time_now[50],time_now_1[50];

    double *data, *mg1, *mg2, *mn, *mu, *mnu1, *mnu2, *mg, *mg1_w, *mg2_w, *mn_w;
    
    double *mg1_bin, *mg2_bin;
    int *num_in_bin, bin_num;

    double *g1t, *g2t, *shear_result;
    double *shear_m;

    int grid_row, grid_col, total_grid, my_grid_st, my_grid_ed;
    double *my_gh1_grid, *my_gh2_grid, *my_chisq_grid;

    int *grid_count;
    
    double chisq_min_fit;
    double *check1;
    double *check2;
    double *gh1_grid;
    double *gh2_grid;
    double *chisq_grid;
    double chisq_ij;

    double dgn, dgu;
    double *weight;

    int i,j,k;
    int shear_num, shear_id;
    int data_row = 7000000, sub_row, data_col;
    double g1, g1_sig, g2, g2_sig;

    double *check_temp;
    double *mg_temp, *mnu_temp;
    double g1_temp, g1_temp_sig, g2_temp, g2_temp_sig;
    double gN, gU, gN_sig, gU_sig,delta_g;

    int ratio;
    int add_time, add_scale;
    add_time = 100;//atoi(argv[5]);
    add_scale = 100;//atoi(argv[6]);
    bin_num = 8;
    data_col = 8;

    strcpy(parent_path, argv[1]);
    strcpy(data_nm, argv[2]);
    
    ratio = atoi(argv[3]);
    grid_row = atoi(argv[4]);  
    grid_col = atoi(argv[5]);    

    dgn = atof(argv[6]);
    dgu = atof(argv[7]);

    
    sub_row = data_row/ratio;

    total_grid = grid_row*grid_col;

    check1 = new double[2*grid_row];
    check2 = new double[2*grid_row];
    check_temp = new double[2*grid_row];

    gh1_grid = new double[total_grid];
    gh2_grid = new double[total_grid];
    if(rank == 0)
    {
        chisq_grid = new double[total_grid];
    }
    grid_count = new int[numprocs];
    task_alloc(total_grid, numprocs, rank, my_grid_st, my_grid_ed,grid_count);
    
    // std::cout<<"Get task"<<std::endl;
    // show_arr(grid_count,1,numprocs);

    my_gh1_grid = new double[grid_count[rank]];
    my_gh2_grid = new double[grid_count[rank]];
    my_chisq_grid = new double[grid_count[rank]];

    

    num_in_bin = new int[bin_num]{};
    mg1_bin = new double[bin_num+1]{};
    mg2_bin = new double[bin_num+1]{};    


    data = new double[data_row*data_col];
    mg = new double[sub_row];
    mg1 = new double[sub_row];
    mg2 = new double[sub_row];
    mn = new double[sub_row];

    mg1_w = new double[sub_row];
    mg2_w = new double[sub_row];
    mn_w = new double[sub_row];

    mg_temp = new double[sub_row];
    mnu_temp = new double[sub_row];

    mu = new double[sub_row];
    mnu1 = new double[sub_row];
    mnu2 = new double[sub_row];

    weight = new double[sub_row];

    sprintf(result_path, "%s/new_pdf_%.4f_%.4f_%s_wU_trueflux.hdf5",parent_path, dgn, dgu, data_nm);

    sprintf(data_path, "%s/shear.hdf5", parent_path);
    
    sprintf(set_name,"/g1");
    read_h5_datasize(data_path, set_name,shear_num);
    g1t = new double[shear_num];
    g2t = new double[shear_num];
    shear_m = new double[8*shear_num];

    read_h5(data_path, set_name, g1t);

    sprintf(set_name,"/g2");
    read_h5_datasize(data_path, set_name,shear_num);
    read_h5(data_path, set_name, g2t);
    if(rank == 0)
    {   
        show_arr(g1t, 1, shear_num);
        show_arr(g2t, 1, shear_num);
        std::cout<<rank<<" "<<sub_row<<" "<<total_grid<<std::endl;
    }
    
    // for(i=0;i<numprocs;i++)
    // {
    //     if(i==rank){std::cout<<my_grid_st<<" "<<my_grid_ed<<std::endl;}
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    for(shear_id=0; shear_id<shear_num; shear_id++)
    {   
        // sprintf(data_name, data_nm, shear_id);
        sprintf(data_path, "%s/data_%d_%s.hdf5",parent_path, shear_id, data_nm);
        sprintf(set_name,"/data");
        read_h5(data_path, set_name, data);

        for(i=0;i<sub_row;i++)
        {   
            weight[i] = data[i*data_col + 7]*data[i*data_col + 7];

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
        if(rank==0){std::cout<<"Read data "<<data_path<<std::endl;}
        
        
        set_bin(mg1, sub_row, mg1_bin, bin_num, 100);
        set_bin(mg2, sub_row, mg2_bin, bin_num, 100);

        find_shear(mg1, mnu1, sub_row, 10, g1, g1_sig, chisq_min_fit, check1, grid_row, 0,100,-0.1, 0.1, 40);
        find_shear(mg2, mnu2, sub_row, 10, g2, g2_sig, chisq_min_fit, check2, grid_row, 0,100,-0.1, 0.1, 40);
        if(rank == 0)
        {            
            sprintf(inform, "PDF g1: (true %9.6f) %9.6f (%9.6f),", g1t[shear_id],g1, g1_sig);
            std::cout<<inform<<std::endl;            
            sprintf(inform, "PDF g2: (true %9.6f) %9.6f (%9.6f)", g2t[shear_id], g2, g2_sig);
            std::cout<<inform<<std::endl;
            
            shear_m[shear_id] = g1;
            shear_m[shear_id + shear_num] = g1_sig;
            shear_m[shear_id + 2*shear_num] = g2;
            shear_m[shear_id + 3*shear_num] = g2_sig;
        }     

        MPI_Barrier(MPI_COMM_WORLD);
        
        delta_g = 1;
        gN = g1;
        gU = g1;
        // while(delta_g>g1_sig*2)
        // {
        //     for(i=0;i<sub_row;i++){mg1_temp[i] = mg1[i] - gN*mn[i];}
        //     find_shear(mg1_temp, mu, sub_row, 10, g1_temp, g1_temp_sig, check_temp, grid_row, 0, 100,-0.1, 0.1, 40);
        // }

        for(j=0;j<5;j++)
        {
            for(i=0;i<sub_row;i++){mg_temp[i] = mg1[i] - gN*mn[i];}
            find_shear(mg_temp, mu, sub_row, 10, gU, gU_sig, chisq_min_fit, check_temp, grid_row, 0, 100,-0.07, 0.07, 60);

            for(i=0;i<sub_row;i++){mg_temp[i] = mg1[i] - gU*mu[i];}
            find_shear(mg_temp, mn, sub_row, 10, gN, gN_sig, chisq_min_fit, check_temp, grid_row, 0, 100,-0.07, 0.07, 60);

            sprintf(inform, "Iter %d. g1_N: %9.6f (%9.6f). g1_U: %9.6f (%9.6f).", j, gN, gN_sig, gU, gU_sig);
            if(rank == 0){std::cout<<inform<<std::endl;}
        }
        MPI_Barrier(MPI_COMM_WORLD);




        find_shear_mean(mg1, mn, sub_row, g1, g1_sig, add_time, add_scale);
        find_shear_mean(mg2, mn, sub_row, g2, g2_sig, add_time, add_scale);

        sprintf(inform, "Ave g1: (true %9.6f) %9.6f (%9.6f)", g1t[shear_id],g1, g1_sig);
        if(rank==0){std::cout<<inform<<std::endl;}
        sprintf(inform, "Ave g2: (true %9.6f) %9.6f (%9.6f)", g2t[shear_id], g2, g2_sig);
        if(rank==0)
        {
            std::cout<<inform<<std::endl;
            shear_m[shear_id + 4*shear_num] = g1;
            shear_m[shear_id + 5*shear_num] = g1_sig;
            shear_m[shear_id + 6*shear_num] = g2;
            shear_m[shear_id + 7*shear_num] = g2_sig;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // g1
        // grid
        for(i=0;i<grid_row;i++)
        {   
            for(j=0;j<grid_col;j++)
            {   
                k = i*grid_col + j;
                gh1_grid[k] = g1 - dgn + dgn*2/(grid_row-1)*i;
                gh2_grid[k] = g1 - dgu + dgu*2/(grid_col-1)*j;
                
                // gh1_grid[k] = check1[grid_row + i];
                // gh2_grid[k] = check1[grid_row + i] - dgu + dgu*2/(grid_col-1)*j;

                // gh1_grid[k] = g1 - dgn + dgn*2/(gh_grid_num-1)*j;
                // gh2_grid[k] = g1 - dgu + dgu*2/(gh_grid_num-1)*j;
            }
        }


        for(i=my_grid_st; i<my_grid_ed; i++)
        {
            for(k=0;k<sub_row;k++){mg[k] = mg1[k] - gh1_grid[i]*mn[k] - gh2_grid[i]*mu[k]/weight[k];}            
            // histogram(mg, mg1_bin, num_in_bin, sub_row, bin_num);
            // cal_chisq_1d(num_in_bin, bin_num, chisq_ij);   
            chisq_Gbin_1d(mg, sub_row, mg1_bin, bin_num, chisq_ij);
            my_chisq_grid[i-my_grid_st] = chisq_ij;
        }
        // if(rank == 0)
        // {
        //     show_arr(gh1_grid,gh_grid_num, gh_grid_num);
        //     show_arr(gh2_grid,gh_grid_num, gh_grid_num);
        //     show_arr(my_chisq_grid,gh_grid_num, gh_grid_num);
        // }
        MPI_Barrier(MPI_COMM_WORLD);
        my_Gatherv(my_chisq_grid, grid_count, chisq_grid, numprocs, rank);

        if(rank == 0)
        {
            sprintf(set_name,"/PDF_new/%d/g1/grid1", shear_id);
            if(shear_id == 0){write_h5(result_path, set_name, gh1_grid, grid_row, grid_col,true);}
            else{write_h5(result_path, set_name, gh1_grid, grid_row, grid_col,false);}

            sprintf(set_name,"/PDF_new/%d/g1/grid2", shear_id);
            write_h5(result_path, set_name, gh2_grid, grid_row, grid_col,false);

            sprintf(set_name,"/PDF_new/%d/g1/chisq_grid", shear_id);
            write_h5(result_path, set_name, chisq_grid, grid_row, grid_col,false);

            sprintf(set_name,"/PDF_new/%d/g1/check", shear_id);
            write_h5(result_path, set_name, check1, 2, grid_row,false);

        }
        MPI_Barrier(MPI_COMM_WORLD);


        //////////////////////////////////////////////////////////////////////////////// 
        //g2
        // grid
        for(i=0;i<grid_row;i++)
        {   
            for(j=0;j<grid_col;j++)
            {   
                k = i*grid_col + j;
                gh1_grid[k] = g2 - dgn + dgn*2/(grid_row-1)*i;
                gh2_grid[k] = g2 - dgu + dgu*2/(grid_col-1)*j;

                // gh1_grid[k] = check2[grid_row + i];
                // gh2_grid[k] = check2[grid_row + i] - dgu + dgu*2/(grid_col-1)*j;

                // gh1_grid[k] = g2 - dgn + dgn*2/(gh_grid_num-1)*j;
                // gh2_grid[k] = g2 - dgu + dgu*2/(gh_grid_num-1)*j;
            }
        }

        for(i=my_grid_st; i<my_grid_ed; i++)
        {
            for(k=0;k<sub_row;k++){mg[k] = mg2[k]- gh1_grid[i]*mn[k] + gh2_grid[i]*mu[k]/weight[k];}            
            // histogram(mg, mg2_bin, num_in_bin, sub_row, bin_num);
            // cal_chisq_1d(num_in_bin, bin_num, chisq_ij);
            chisq_Gbin_1d(mg, sub_row, mg2_bin, bin_num, chisq_ij);   
            my_chisq_grid[i-my_grid_st] = chisq_ij;
        }
        // if(rank==0)
        // {
        //     std::cout<<std::endl;
        //     show_arr(gh1_grid,gh_grid_num, gh_grid_num);
        //     show_arr(gh2_grid,gh_grid_num, gh_grid_num);
        //     show_arr(my_chisq_grid,gh_grid_num, gh_grid_num);
        // }
        MPI_Barrier(MPI_COMM_WORLD);
        
        my_Gatherv(my_chisq_grid, grid_count, chisq_grid, numprocs, rank);

        if(rank == 0)
        {
            sprintf(set_name,"/PDF_new/%d/g2/grid1", shear_id);
            write_h5(result_path, set_name, gh1_grid, grid_row, grid_col,false);

            sprintf(set_name,"/PDF_new/%d/g2/grid2", shear_id);
            write_h5(result_path, set_name, gh2_grid, grid_row, grid_col,false);

            sprintf(set_name,"/PDF_new/%d/g2/chisq_grid", shear_id);
            write_h5(result_path, set_name, chisq_grid, grid_row, grid_col,false);

            sprintf(set_name,"/PDF_new/%d/g2/check", shear_id);
            write_h5(result_path, set_name, check2, 2, grid_row,false);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        

    }

    if(rank == 0)
    {
        sprintf(set_name,"/PDF_new/shear");
        write_h5(result_path, set_name, shear_m, 8, shear_num, false);
        std::cout<<"Finish"<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // delete[] data;
    // delete[] mg1;
    // delete[] mg2;
    // delete[] mn;
    // delete[] mu;
    // delete[] mnu1;
    // delete[] mnu2;
    // delete[] g2t;
    // delete[] g1t;
    // delete[] check1;
    // delete[] check2;
    // delete[] gh1_grid;
    // delete[] gh2_grid;
    // delete[] chisq_grid;
    // delete[] grid_count;
    // delete[] my_gh1_grid;
    // delete[] my_gh2_grid;
    // delete[] my_chisq_grid;
    // delete[] shear_m;
    // delete[] num_in_bin;
    // delete[] mg1_bin;
    // delete[] mg2_bin;

    MPI_Finalize();
    return 0;
}