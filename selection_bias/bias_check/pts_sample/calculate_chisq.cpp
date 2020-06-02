#include<FQlib.h>
#include<hk_iolib.h>
#include<hk_mpi.h>

#define SAVE_MEM

#ifdef SAVE_MEM
#define MY_FLOAT float
#else
#define MY_FLOAT double
#endif

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

    char parent_path[200], shear_path[200], result_path[200];
    char set_name[30], inform[300], time_now[40];

    char data_path[300], *data_type[10], data_comb[100], *mg_name[5];
    int data_type_num;

    int i,j,k,m;

    int shear_num, shear_st, shear_ed;
    int*shear_point;

    
    int mg_bin_num, *num_in_bin;
    int chisq_num;
    MY_FLOAT *mg_bins;
    MY_FLOAT *chisq1, *chisq2, *shear_for_chi, chisq_min_fit;
    MY_FLOAT left_guess, right_guess;
    MY_FLOAT gh1, gh1_sig, gh2, gh2_sig;
    MY_FLOAT *result;

    int data_col, data_row;
    MY_FLOAT *mg_data[5], *temp_read[5];
    MY_FLOAT *rotation;

    // data shape
    strcpy(parent_path, argv[1]);
    shear_num = atoi(argv[2]);
    data_row = atoi(argv[3])*10000;
    mg_bin_num = atoi(argv[4]);
    for(i=5; i<argc; i++)
    {
        data_type[i-5] = new char[40];
        strcpy(data_type[i-5], argv[i]);
        if(rank == 0){std::cout<<data_type[i-5]<<std::endl;}
    }
    data_type_num = argc-5;
    data_col = 5;// G1, G2, N, U, V
    chisq_num = 101;
    left_guess = -0.1;
    right_guess = 0.1;
    rotation = new MY_FLOAT[5];

    for(i=0; i<data_col; i++)
    {
        mg_name[i] = new char[50];
        mg_data[i] = new MY_FLOAT[data_row]{};
        temp_read[i] = new MY_FLOAT[data_row]{};
        
    }

    sprintf(mg_name[0], "/mg1");
    sprintf(mg_name[1], "/mg2");
    sprintf(mg_name[2], "/mn");
    sprintf(mg_name[3], "/mu");
    sprintf(mg_name[4], "/mv");

    char_stack(data_type, data_type_num, data_comb);

    result = new MY_FLOAT[4];

    shear_point = new int[numprocs]{};

    mg_bins = new MY_FLOAT[mg_bin_num+1]{};
    num_in_bin = new int[mg_bin_num]{};
    chisq1 = new MY_FLOAT[chisq_num]{};
    chisq2 = new MY_FLOAT[chisq_num]{};
    shear_for_chi = new MY_FLOAT[chisq_num]{};
    for(i=0;i<chisq_num;i++)
    {
        shear_for_chi[i] = 0.2/(chisq_num-1)*i - 0.1;
    }
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

    MPI_Barrier(MPI_COMM_WORLD);
    for(i=0;i<numprocs;i++)
    {
        if(i == rank)
        {
            std::cout<<rank<<" "<<shear_st<<" "<<shear_ed<<" "<<data_row<<" "<<data_col<<" "<<mg_bin_num<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // read and calculate
    sprintf(set_name,"/data");
    for(i=shear_st;i<shear_ed;i++)
    {   
        // read data
        for(j=0;j<data_type_num;j++)
        {
            sprintf(data_path,"%s/data_%s_%d.hdf5",parent_path, data_type[j], i);
            for(k=0;k<data_col;k++)
            {
                read_h5(data_path, mg_name[k], temp_read[k]);
            }

            if(rank== 0){std::cout<<"Read data "<<data_path<<" "<<data_type[j]<<std::endl;}
            for(k=0;k<data_col;k++)
            {
                for(m=0;m<data_row;m++)
                {
                    // estimator_rotation(Pi/2, data[k*data_col],data[k*data_col+1],data[k*data_col+2],data[k*data_col+3],data[k*data_col+4],rotation);
                    // mg1[k] += rotation[0];
                    // mg2[k] += rotation[1];
                    // mn[k] += rotation[2];
                    // mnu1[k] += rotation[2] + rotation[3];
                    // mnu2[k] += rotation[2] - rotation[3];

                    mg_data[k][m] += temp_read[k][m];
                }
            }           
            
        }

        // calculate chisq
        set_bin(mg_data[0], data_row, mg_bins, mg_bin_num, 100);

        for (k = 0; k < chisq_num; k++)
        {	
            try
            {
                fourier_hist(mg_data[0], mg_data[2], mg_data[3], data_row, shear_for_chi[k], 1, mg_bins, num_in_bin, mg_bin_num);
			    cal_chisq_1d(num_in_bin, mg_bin_num, left_guess);
            }
            catch(const char *msg)
            {
                throw msg;
            }
            chisq1[k] = left_guess;
            
            try
            {
                fourier_hist(mg_data[1], mg_data[2], mg_data[3], data_row, shear_for_chi[k], 2, mg_bins, num_in_bin, mg_bin_num);
			    cal_chisq_1d(num_in_bin, mg_bin_num, left_guess);
            }
            catch(const char *msg)
            {
                throw msg;
            }
            chisq2[k] = left_guess;
        }
        
        sprintf(data_path, "%s/chisq_%d/chisq_%d_%s.hdf5", parent_path, mg_bin_num,rank, data_comb);
        sprintf(set_name,"/chisq1");
        write_h5(data_path, set_name, chisq1, chisq_num, 1, true);
        sprintf(set_name,"/chisq2");
        write_h5(data_path, set_name, chisq2, chisq_num, 1, false);
        sprintf(set_name,"/shear");
        write_h5(data_path, set_name, shear_for_chi, chisq_num, 1, false);    


    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){std::cout<<"Finish"<<std::endl;}

    MPI_Finalize();
    return 0;
}
