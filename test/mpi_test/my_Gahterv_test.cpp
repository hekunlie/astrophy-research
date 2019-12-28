#include<hk_mpi.h>
#include<FQlib.h>

int main(int argc, char **argv)
{
	int rank, numprocs, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	MPI_Status status;

	double t1, t2, t3, s;
	int i, j, num, sub_num;
    int my_st, my_ed;
	num = atoi(argv[1]);

	double *send_buf, *re_gather_buf;
	double *recv_buf = new double[num] {};
    double *sub_buf[6];

	int *send_num = new int[numprocs];
    int *entry_for_gather = new int[numprocs];

    task_alloc(num, numprocs, rank, my_st, my_ed, send_num, entry_for_gather);
    for(i=0;i<6;i++)
    {
        sub_buf[i] = new double[send_num[rank]];
        for(j=0;j<send_num[rank];j++)
        {
            sub_buf[i][j]=rank+i;
        }
    }
    for(i=0;i<6;i++)
    {
        my_Gatherv(sub_buf[i], send_num, recv_buf, numprocs, rank);
        if(rank==0)
        {
            show_arr(recv_buf,1,num);
        }
    }
    MPI_Finalize();
	return 0;
}
