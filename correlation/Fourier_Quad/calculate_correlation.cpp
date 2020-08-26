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
    int my_fnm, target_fnm,fnm_st, total_field_num, label, num_label;


    strcpy(source_list, argv[1]);
    strcpy(field_info.parent_path,argv[2]);

    field_info.zbin_label_0 = atoi(argv[3]);
    field_info.zbin_label_1 = atoi(argv[4]);

    total_field_num = atoi(argv[5]);

    sprintf(log_path, "%s/log/%d_log.dat",field_info.parent_path, rank);
    

    // read the information of each exposure file
    initialize(source_list, &field_info, total_field_num, numprocs, rank);

    // read the catalog of redshift bin z1 & z2
    read_field_data(&field_info);


    

    if(rank == 0){std::cout<<"Initialization"<<std::endl;}
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef PRINT_INFO
    if(rank==0)
    {
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

        std::cout<<"Field each cpu:"<<std::endl;
        show_arr(field_info.field_num_each_rank,1,numprocs);
        std::cout<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for(j=0; j<numprocs;j++)
    {
        if(j==rank)
        {
            std::cout<<rank<<" My field: "<<field_info.my_field_st<<"~"<<field_info.my_field_ed<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

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

    initialize_total_chi_block(&field_info);

    // loop the fields
    for(my_fnm=field_info.my_field_st; my_fnm < field_info.my_field_ed; my_fnm++)
    {   
        initialize_field_chi_block(&field_info,my_fnm);
        
        ////////////////  search pairs in the current field  ////////////////////
        st1 = clock();
        sprintf(log_inform,"Start: %d-%s num: %d",my_fnm, field_info.field_name[my_fnm], field_info.total_gal_num_z1[my_fnm]);
        if(rank == 0){std::cout<<log_inform<<std::endl;}
        write_log(log_path, log_inform);
        
        sprintf(log_inform,"Search own field: %d-%s",my_fnm,field_info.field_name[my_fnm]);
        if(rank == 0){std::cout<<log_inform<<std::endl;}
        write_log(log_path, log_inform);

        if (field_info.total_gal_num_z1[my_fnm] > 0)
        {
            find_pairs_same_field(&field_info, my_fnm);
            num_label = 1;
        }
        else
        {   
            num_label= 0;
            continue;
        }
        st2 = clock();
        tt =  (st2 - st1)/CLOCKS_PER_SEC;

        sprintf(log_inform,"Finish own field: %d-%s. %.2f sec", my_fnm, field_info.field_name[my_fnm], tt);
        if(rank == 0){std::cout<<log_inform<<std::endl;}
        write_log(log_path, log_inform);

        // collect_chi_block(&field_info, my_fnm);
        // fnm_st = field_info.theta_bin_num*field_info.chi_guess_num*field_info.mg_bin_num;
        // show_arr(field_info.total_num_count_chit,fnm_st,10);
        // MPI_Barrier(MPI_COMM_WORLD);
        
        // if(num_label == 1)
        // {   
        //     if(field_info.zbin_label_0 ==  field_info.zbin_label_1){fnm_st = my_fnm + 1;}
        //     else{fnm_st = 0;}
        //     ////////////////// search pairs in the target fields  //////////////////////
        //     for(target_fnm=fnm_st; target_fnm<total_field_num; target_fnm++)
        //     {   
        //         if(target_fnm != my_fnm)
        //         {
        //             st3 = clock();

        //             field_distance(&field_info, my_fnm, target_fnm, label);

        //             if(label == 1 and field_info.total_gal_num_z2[target_fnm]>0)
        //             {

        //                 sprintf(log_inform,"Search: (%d-%s) %d-%s num: %d",my_fnm,field_info.field_name[my_fnm],
        //                 target_fnm, field_info.field_name[target_fnm],field_info.total_gal_num_z2[target_fnm]);
        //                 if(rank == 0){std::cout<<log_inform<<std::endl;}
        //                 write_log(log_path, log_inform);

        //                 find_pairs_diff_field(&field_info, my_fnm, target_fnm);

        //                 st4 = clock();
        //                 tt = (st4-st3)/CLOCKS_PER_SEC;
        //                 sprintf(log_inform,"Finish search: (%d-%s) %d-%s. %.2f sec",my_fnm,field_info.field_name[my_fnm],
        //                 target_fnm,field_info.field_name[target_fnm], tt);
        //                 if(rank == 0){std::cout<<log_inform<<std::endl;}
        //                 write_log(log_path, log_inform);
        //             }
        //             else
        //             {
        //                 continue;
        //             }
        //         }
        //     }
        // }

        st5 = clock();
        tt = (st5-st1)/CLOCKS_PER_SEC;
        sprintf(log_inform,"Finish: %d-%s. %.2f sec",my_fnm, field_info.field_name[my_fnm], tt);
        if(rank == 0)
        {
            std::cout<<log_inform<<std::endl;
            std::cout<<"===================================================="<<std::endl;
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
