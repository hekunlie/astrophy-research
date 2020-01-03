#include<FQlib.h>
#include<hk_mpi.h>


int main(int argc, char**argv)
{
    int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

    double t1,t2;
    int i,j,k;

    char filter_name[50];
    char cut_name[50];
    char time_now[50];

    int buff_len= 50;
    double *buff;
    double *recv_buff = new double[buff_len]{};

    if(rank == 0)
    {
        buff = new double[buff_len]{};
        for(i=0;i<buff_len;i++){buff[i]=1;}
    }

    my_Send_all(buff, recv_buff, buff_len, numprocs, 0, rank);
    
    // strcpy(cut_name, argv[1]);
    // strcpy(filter_name, argv[2]);

    char inform[200];
    //sprintf(inform, "%d, %s, %s",rank,cut_name, filter_name);
    t1 = clock();

    for(int i=0;i<numprocs;i++)
    {
        if(rank == i)
        {   
            std::cout<<"Rank "<<rank<<" ";
            show_arr(recv_buff,1,buff_len);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // for(int j=0;j<rank*1000000+1;j++)
    // {
    //     sin(10)*exp(10);
    // }
    // t2 = clock();

    // std::cout<<inform<<" "<<(t2-t1)/CLOCKS_PER_SEC<<std::endl;

    delete[] recv_buff;
    if(rank==0)
    delete[] buff;
    MPI_Finalize();
    return 0;
}