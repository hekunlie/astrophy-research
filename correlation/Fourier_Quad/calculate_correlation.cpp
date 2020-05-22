// #include<FQlib.h>
#include<hk_mpi.h>
#include<functions.h>



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

    char source_list[300];
    char data_path[300], result_path[300];
    char set_name[50], temp_char[50];

    data_info field_info;

    int i,j,k,m,n;
    int total_expos_num, total_field_num;

    strcpy(source_list, argv[1]);
    strcpy(data_path,argv[2]);
    total_expos_num = atoi(argv[3]);
    total_field_num = atoi(argv[4]);



    // read the information of each exposure file
    initialize(source_list, &field_info, total_expos_num, numprocs, rank);
    std::cout<<"Initialization"<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);


    find_pairs(&field_info);



    for(j=0; j<numprocs;j++)
    {
        if(j==rank)
        {
            std::cout<<total_expos_num<<std::endl;

            for(i=0;i<total_expos_num;i++)
            {
                std::cout<<field_info.field_name[i]<<" "<<field_info.exposure_name[i]
                <<" "<<field_info.field_label[i]
                <<" "<<field_info.exposure_label[i]<<" "<<field_info.exposure_num_of_field[i]
                <<" "<<field_info.field_cen_ra[i]<<" "<<field_info.field_cen_dec[i]
                <<" "<<field_info.delta_ra[i]<<" "<<field_info.delta_dec[i]
                <<" "<<field_info.my_exposure_num<<" "<<field_info.my_exposure_st
                <<" "<<field_info.my_exposure_ed<<std::endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    
    return 0;
}
