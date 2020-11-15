#include<hk_mpi.h>
#include<functions_expo_wise.h>
#include<ctime>
#include<vector>


int main(int argc, char **argv)
{
    int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);


    char inform[400];
    char time_now[40];
    
    int i;
    double st1, st2, st3, st4, st5;
    double tt1, tt2;

    corr_cal all_paras;


    strcpy(all_paras.parent_path, argv[1]);
    all_paras.resample_num = atoi(argv[2]);
    all_paras.corr_cal_result_file_num = atoi(argv[3]);
    all_paras.corr_cal_thread_num = numprocs;
    all_paras.corr_cal_rank = rank;

    st1 = clock();

    sprintf(inform,"Prepare data");
    if(rank == 0){std::cout<<inform<<std::endl;}
    
    read_para(&all_paras);

    // create index file
    if(rank == 0){prepare_data(&all_paras, 0);}
    MPI_Barrier(MPI_COMM_WORLD);

    // read index file
    if(rank > 0){prepare_data(&all_paras, 1);}
    MPI_Barrier(MPI_COMM_WORLD);


    sprintf(inform,"Calculate %d ~ %d",all_paras.jackknife_resample_st[rank],all_paras.jackknife_resample_ed[rank]);
    if(all_paras.corr_cal_rank == 0){std::cout<<inform<<std::endl;}

    for(i=all_paras.jackknife_resample_st[rank]; i<all_paras.jackknife_resample_ed[rank]; i++)
    {   
        st2 = clock();
        get_time(time_now, 40);
        if(rank == 0)
        {     
            sprintf(inform, "%s. Jackknife %d. start.", time_now, i);
            std::cout<<inform<<std::endl;
        }


        resample_jackknife(&all_paras,i);
        corr_calculate(&all_paras, i);

        st3 = clock();
        get_time(time_now, 40);

        if(rank == 0)
        {   
            tt1 = (st3-st2)/CLOCKS_PER_SEC;    
            sprintf(inform, "%s. Jackknife %d. %.2f sec", time_now, i,tt1);
            std::cout<<inform<<std::endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    st4 = clock();

    // write down the result
    for(i=0; i<numprocs; i++)
    {
        if(i==rank)
        {
            save_result(&all_paras);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    st5 = clock();

    if(rank == 0)
    {   
        tt1 = (st4-st1)/CLOCKS_PER_SEC;
        tt2 = (st5-st4)/CLOCKS_PER_SEC;

        sprintf(inform,"Finish. %.2f sec, %.2f sec.", tt1, tt2);
        std::cout<<inform<<std::endl;
    }


    MPI_Finalize();
    return 0;
}