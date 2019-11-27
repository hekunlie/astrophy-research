#include<mpi.h>
#include<ctime>
#include<fstream>
#include<iostream>

void write_log(char*filename, char *inform)
{   
    std::ofstream default_loggers;
	char time_now[40];
	time_t timer;
	time(&timer);
	strftime(time_now, 40, "%Y-%m-%d %H:%M:%S", localtime(&timer));

	default_loggers.open(filename, std::ios::out | std::ios::app);
	default_loggers << time_now << " ---- " << inform << std::endl;
	default_loggers.close();
}

int main(int argc, char **argv)
{
    int myid, numprocs, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Get_processor_name(processor_name, &namelen);

    char inform[200];
    char log_path[200];
    sprintf(log_path,"logs.dat");

    MPI_Finalize();
    return 0;
}