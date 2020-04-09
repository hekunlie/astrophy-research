#include<FQlib.h>
#include<hk_iolib.h>
#include<hk_mpi.h>
#define SKYAREA_NUM_MAX 30

int main(int argc, char *argv[])
{
    /*
        argv[1]: the parent path that includes "./data", "./result"
        argv[2] ~ argv[i]: the sky area name in the "catalog.hdf5", like "w1", "w2"...
    */

	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    char parent_path[300], data_path[300], result_path[300];
    char set_name[30];
    // for the set_name of each area
    char *area_name[SKYAREA_NUM_MAX];
    
    strcpy(parent_path, argv[1]);
    sprintf(data_path,"%s/data/catalog.hdf5", parent_path);
    sprintf(result_path,"%s/result/result.hdf5", parent_path);

    int i,j,k,m,n;
    int sky_area_num;
    sky_area_num = argc - 2;
    for(i=0;i<sky_area_num;i++)
    {
        area_name[i] = new char[30];
        strcpy(area_name[i], argv[i+2]);
    }
    // the 0'th is the program name, the 1'th is the parent_path
    


    return 0;
}
