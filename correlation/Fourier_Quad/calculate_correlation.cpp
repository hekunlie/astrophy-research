// #include<FQlib.h>
#include<hk_mpi.h>
#include<functions.h>
#include<ctime>
#include<FQlib.h>

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
    int total_expos_num;
    int red_shift_bin_num, z1,z2;

    int ir, radius_bin_num;
    int my_fnm, target_fnm, total_field_num;

    strcpy(source_list, argv[1]);
    strcpy(data_path,argv[2]);
    z1 = atoi(argv[3]);
    z2 = atoi(argv[4]);

    total_field_num = atoi(argv[5]);



    // read the information of each exposure file
    initialize(source_list, &field_info, total_field_num, numprocs, rank);
    if(rank == 0){std::cout<<"Initialization"<<std::endl;}
    MPI_Barrier(MPI_COMM_WORLD);


    if(rank==0)
    {
        for(i=0;i<total_field_num;i++)
        {
            std::cout<<field_info.field_name[i]<<"  "<<field_info.exposure_num_of_field[i]
            <<"  "<<field_info.field_cen_ra[i]<<"  "<<field_info.field_cen_dec[i]
            <<"  "<<field_info.delta_ra[i]<<"  "<<field_info.delta_dec[i]
            <<"  "<<field_info.delta_len[i]<<"  "<<field_info.field_cen_cos_dec[i]<<std::endl;
        }
        for(i=0; i<numprocs; i++){std::cout<<field_info.field_num_each_rank[i]<<" ";}
        std::cout<<std::endl;
    }

    for(j=0; j<numprocs;j++)
    {
        if(j==rank)
        {
            std::cout<<rank<<" "<<field_info.my_field_st<<" "<<field_info.my_field_ed<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0){std::cout<<"Redshifit: "<<z1<<" "<<z2<<std::endl;}


    read_field_data(&field_info, z1, z2);
    if(rank == 0)
    {
        for(i=0; i<total_field_num; i++)
        {
            std::cout<<field_info.field_data_col<<"  "<<field_info.total_gal_num_z1[i]<<"  "<<field_info.total_gal_num_z2[i]<<std::endl;
        }
        for(i=0; i<field_info.zbin_num+1; i++)
        {std::cout<<field_info.zbin[i]<<"  ";}
        std::cout<<std::endl;
        
        for(i=0; i<field_info.theta_bin_num+1; i++)
        {std::cout<<field_info.theta_bin[i]<<"  ";}
        std::cout<<std::endl;
    }


    // loop the fields
    for(my_fnm=0; my_fnm < total_field_num; my_fnm++)
    {
        // loop the radius bin
        for(ir=0;ir<radius_bin_num;ir++)
        {
            // search the target fields
            for(target_fnm=my_fnm; target_fnm<total_field_num; target_fnm++)
            {
                // find the needed fields
                // dx^2 + dy^2 - sqrt(dra^2+ddec^2)*2 <= max radius
                ;
            }
        }

    }

    
    MPI_Finalize();
    
    return 0;
}
