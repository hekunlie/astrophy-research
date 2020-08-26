#include<hk_mpi.h>
#include<functions.h>
#include<ctime>


#define NOT_PRINT_INFO

int main(int argc, char *argv[])
{
    /*

    */

	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    char source_list[300];
    char data_path[300], result_path[300];
    char set_name[50], temp_char[50];
    char log_inform[500], log_path[600];

    double st1, st2, st3, st4, st5, tt;

    data_info field_info;

    int i,j;
    int fnm_1, fnm_2, total_field_num, label, num_label;


    strcpy(source_list, argv[1]);
    strcpy(field_info.parent_path,argv[2]);

    field_info.zbin_label_0 = atoi(argv[3]);
    field_info.zbin_label_1 = atoi(argv[4]);

    total_field_num = atoi(argv[5]);

    sprintf(log_path, "%s/log/%d_log.dat",field_info.parent_path, rank);
    

    // read the information of each exposure file
    initialize(source_list, &field_info, total_field_num);

    // read the catalog of redshift bin z1 & z2
    read_field_data(&field_info);

    // find all the potential field pair for calculation (i, j),
    // does not include the field itself i!= j
    task_prepare(numprocs, rank, &field_info);

    

    ////////////////////////////////  PRINT INFO  ////////////////////////////////////////////
    if(rank == 0){std::cout<<"Initialization"<<std::endl;}
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank==0)
    {
#ifdef PRINT_INFO
        for(i=0;i<total_field_num;i++)
        {
            std::cout<<field_info.field_name_path[i]<<"  "<<field_info.field_name[i]<<"  "<<field_info.exposure_num_of_field[i]
            <<"  "<<field_info.field_cen_ra[i]<<"  "<<field_info.field_cen_dec[i]
            <<"  "<<field_info.field_delta_ra[i]<<"  "<<field_info.field_delta_dec[i]
            <<"  "<<field_info.field_delta_len[i]<<"  "<<field_info.field_cen_cos_dec[i]<<std::endl;
            
            std::cout<<"block num: "<<field_info.block_num[i]
            <<"  total gal nun z1: "<<field_info.total_gal_num_z1[i]
            <<"  total gal nun z2: "<<field_info.total_gal_num_z2[i]<<std::endl;
        }
#endif
        std::cout<<"Field pair for each cpu:"<<std::endl;
        show_arr(field_info.field_pair_num_each_rank,1,numprocs);
        std::cout<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    sprintf(log_inform, "Rank %d. Total field pairs: %d. My field pairs: %d ~ %d.", rank, field_info.field_pair_num,field_info.my_field_pair_st,field_info.my_field_pair_ed);
    write_log(log_path, log_inform);
    for(j=0; j<numprocs;j++)
    {
        if(j==rank)
        {
            std::cout<<log_inform<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);


    if(rank == 0)
    {
        std::cout<<std::endl<<"Calculate Redshift bin: "<<field_info.zbin_label_0<<" "<<field_info.zbin_label_1<<std::endl;

        std::cout<<"Redshift bin:"<<std::endl;
        show_arr(field_info.zbin, 1, field_info.zbin_num+1);
        sprintf(log_inform,"Z1: %.2f ~ %.2f, Z2: %.2f ~ %.2f",field_info.zbin[field_info.zbin_label_0],
        field_info.zbin[field_info.zbin_label_0+1],field_info.zbin[field_info.zbin_label_1],field_info.zbin[field_info.zbin_label_1+1]);
        std::cout<<log_inform<<std::endl;
        
        std::cout<<"Radius bin:"<<std::endl;
        show_arr(field_info.theta_bin, 1, field_info.theta_bin_num+1);
        std::cout<<std::endl;


        std::cout<<"Chi guess: "<<field_info.chi_guess_num<<" points"<<std::endl;
        show_arr(field_info.chi_guess, 1, field_info.chi_guess_num);
        std::cout<<std::endl;
        std::cout<<"Bin num for PDF_SYM: "<<field_info.mg_bin_num<<" "<<field_info.mg_bin_num1<<" "<<field_info.mg_bin_num2<<" "<<field_info.mg_bin_num3<<std::endl;
        show_arr(field_info.mg_bin, 1, field_info.mg_bin_num+1);

        std::cout<<std::endl<<field_info.mg_bin_num<<" "<<field_info.chi_block_len<<" "<<field_info.ir_chi_block_len<<" "<<
        field_info.iexpo_chi_block_len<<std::endl<<std::endl;

        std::cout<<std::endl<<field_info.parent_path<<std::endl;
        std::cout<<std::endl<<field_info.gg_len<<std::endl<<std::endl;     
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ////////////////////////////////  PRINT INFO-end  ////////////////////////////////////////////



    initialize_total_chi_block(&field_info);
    // initialize_field_chi_block(&field_info,my_fnm);

    // loop the field pairs
    for(i=field_info.my_field_pair_st; i < field_info.my_field_pair_ed; i++)
    {   
        // field pair label
        fnm_1 = field_info.field_pair_label_1[i];
        fnm_2 = field_info.field_pair_label_2[i];

        
        ////////////////  search pairs in the current field  ////////////////////
        st1 = clock();
        sprintf(log_inform,"Start: %d-%s(%d) <-> %d-%s(%d)",fnm_1, field_info.field_name[fnm_1], field_info.total_gal_num_z1[fnm_1],
        fnm_2, field_info.field_name[fnm_2],field_info.total_gal_num_z2[fnm_2]);
        if(rank == 0){std::cout<<log_inform<<std::endl;}
        write_log(log_path, log_inform);
        
        // if (fnm_1 == fnm_2)
        // {
        //     find_pairs_same_field(&field_info, fnm_1);
        // }
        // else
        // {   
        //     ;//find_pairs_diff_field(&field_info, fnm_1, fnm_2);
        // }

        find_pairs(&field_info, fnm_1, fnm_2);

        st2 = clock();
        tt =  (st2 - st1)/CLOCKS_PER_SEC;

        sprintf(log_inform,"Finish: %d-%s(%d) <-> %d-%s(%d) in %.2f sec.",fnm_1, field_info.field_name[fnm_1], field_info.total_gal_num_z1[fnm_1],
        fnm_2, field_info.field_name[fnm_2],field_info.total_gal_num_z2[fnm_2], tt);
        if(rank == 0)
        {
            std::cout<<log_inform<<std::endl;
            std::cout<<"================================================================="<<std::endl;
        }
        write_log(log_path, log_inform);
        if(num_label == 1)
        {
            ;// collect_chi_block(&field_info, my_fnm);
            // save_field_chi_block(&field_info, my_fnm);
        }        
        
    }

    
    MPI_Finalize();
    
    return 0;
}
