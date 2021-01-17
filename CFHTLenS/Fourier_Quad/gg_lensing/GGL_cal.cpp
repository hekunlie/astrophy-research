#include<GG_lensing_functions_expo_wise.h>
#include<iostream>


int main(int argc, char *argv[])
{

    int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    ggl_data_info data_info;
    strcpy(data_info.ggl_total_path, argv[1]);
    data_info.jack_num = atoi(argv[2]);
    data_info.rank = rank;
    data_info.numprocs = numprocs;


    ggl_initialize(&data_info);

    
    for( long i=0;i<10000000;i++)
    {atan2(tan(i*i*i*i*i*i),cos(i*i*i*i))* cos(i*i*i*i)*tan(i*i*i*i*i*i);}
    return 0;
}