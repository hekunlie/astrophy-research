#include<hk_mpi.h>
#include<functions_expo_wise.h>
#include<ctime>
#include<vector>


void read_para(corr_cal *all_paras)
{   
    int i, j, m, n;
    char set_name[60], data_path[600];

    sprintf(data_path,"%s/cata/gg_cor.hdf5", all_paras->parent_path);
    // read radius bin
    sprintf(set_name,"/theta_bin");
    read_h5_datasize(data_path, set_name,j);
    all_paras->theta_bin = new MY_FLOAT[j]{};
    all_paras->theta_bin_num = j -1;
    read_h5(data_path, set_name, all_paras->theta_bin);


    // read redshift bin
    sprintf(set_name,"/redshift_bin");
    read_h5_datasize(data_path, set_name,j);
    all_paras->zbin = new MY_FLOAT[j]{};
    all_paras->zbin_num = j -1;
    read_h5(data_path, set_name, all_paras->zbin);


    sprintf(set_name, "/chi_guess");
    read_h5_datasize(data_path, set_name, all_paras->chi_guess_num);
    all_paras->chi_guess = new MY_FLOAT[all_paras->chi_guess_num];
    all_paras->corr_cal_chi_guess = new double[all_paras->chi_guess_num];
    read_h5(data_path, set_name, all_paras->chi_guess);
    for(j=0; j<all_paras->chi_guess_num; j++){all_paras->corr_cal_chi_guess[j] = all_paras->chi_guess[j];}


    sprintf(set_name, "/mg_bin");
    read_h5_datasize(data_path, set_name, j);
    all_paras->mg_bin = new MY_FLOAT[j];
    all_paras->mg_bin_num = j -1;
    read_h5(data_path, set_name, all_paras->mg_bin);


    std::cout<<"G bin: "<<all_paras->mg_bin_num<<std::endl;
    show_arr(all_paras->mg_bin,1,all_paras->mg_bin_num+1);
    std::cout<<std::endl;

    std::cout<<"Chi guess bin: "<<all_paras->chi_guess_num<<std::endl;
    show_arr(all_paras->corr_cal_chi_guess,1,all_paras->chi_guess_num);
    std::cout<<std::endl;
    
    std::cout<<"theta bin: "<<all_paras->theta_bin_num<<std::endl;
    show_arr(all_paras->theta_bin,1,all_paras->theta_bin_num+1);
    std::cout<<std::endl;

    std::cout<<"Z bin: "<<all_paras->zbin_num<<std::endl;
    show_arr(all_paras->zbin,1,all_paras->zbin_num+1);
    std::cout<<std::endl;


    // the fundmental block size for number counting in the PDF_SYM
    // the smallest block, mg_bin_num x mg_bin_num, for each chi guess point
    all_paras->chi_block_len = all_paras->mg_bin_num*all_paras->mg_bin_num;
    // the upper level of the above block, chi_guess_num x mg_bin_num x mg_bin_num,
    // for each theta_bin, there're ''theta_bin_num'' of such blocks in each exposure
    all_paras->ir_chi_block_len = all_paras->chi_guess_num*all_paras->chi_block_len;
    // the biggest block, for each exposure, there're ''zbin_num*zbin_num'' of such blocks in each expo
    all_paras->iz_chi_block_len = all_paras->ir_chi_block_len*all_paras->theta_bin_num;
    // tomography, z[i,j] = z[j,i], save zbin_num*zbin_num blocks for saving the time in calculation
    // [i,j] = [j,i], they will be sum after the calculation
    all_paras->expo_chi_block_len = all_paras->iz_chi_block_len*all_paras->zbin_num*all_paras->zbin_num;
    // because [j,i] will be added to [i,j], zbin pair [i,j] = [j, i]
    all_paras->expo_chi_block_len_true = all_paras->iz_chi_block_len*((all_paras->zbin_num*all_paras->zbin_num + all_paras->zbin_num)/2);

    all_paras->theta_accum_len_true = all_paras->theta_bin_num*((all_paras->zbin_num*all_paras->zbin_num + all_paras->zbin_num)/2);


    // tangential and cross components
    all_paras->expo_num_count_chit = new double[all_paras->expo_chi_block_len_true]{};
    all_paras->expo_num_count_chix = new double[all_paras->expo_chi_block_len_true]{};

    all_paras->theta_accum = new double[all_paras->theta_accum_len_true]{};
    all_paras->theta_num_accum = new double[all_paras->theta_accum_len_true]{};


    // stack all the exposure data
    all_paras->corr_cal_stack_num_count_chit = new double[all_paras->expo_chi_block_len_true]{};
    all_paras->corr_cal_stack_num_count_chix = new double[all_paras->expo_chi_block_len_true]{};

    all_paras->corr_cal_stack_expo_theta_accum = new double[all_paras->theta_accum_len_true]{};
    all_paras->corr_cal_stack_expo_theta_num_accum = new double[all_paras->theta_accum_len_true]{};

    // for the PDF calculation
    all_paras->corr_cal_chi_num = all_paras->chi_guess_num*all_paras->theta_bin_num*((all_paras->zbin_num*all_paras->zbin_num + all_paras->zbin_num)/2);
    all_paras->corr_cal_final_data_num = all_paras->theta_bin_num*((all_paras->zbin_num*all_paras->zbin_num + all_paras->zbin_num)/2);
}

void prepare_data(corr_cal *all_paras) 
{   
    int i, j, k, m ,n;
    char data_path[550];
    char set_name[50];
    int *temp[2];

    all_paras->corr_cal_total_pair_num = 0;
    // no 0'th file, because the CPU 0 is the master for the task distribution
    for(i=1; i<all_paras->corr_cal_result_file_num; i++)
    {
        sprintf(data_path, "%s/result/core_%d_num_count.hdf5", all_paras->parent_path,i);
        sprintf(set_name, "/pair_1");
        read_h5_datasize(data_path, set_name, j);
        all_paras->corr_cal_total_pair_num += j;
    }
    all_paras->corr_cal_expo_pair_label[0] = new int[all_paras->corr_cal_total_pair_num]{};
    all_paras->corr_cal_expo_pair_label[1] = new int[all_paras->corr_cal_total_pair_num]{};
    all_paras->corr_cal_expo_pair_file_label = new int[all_paras->corr_cal_total_pair_num]{};

    n = 0;
    for(i=1; i<all_paras->corr_cal_result_file_num; i++)
    {
        sprintf(data_path, "%s/result/core_%d_num_count.hdf5", all_paras->parent_path,i);
        sprintf(set_name, "/pair_1");
        read_h5_datasize(data_path, set_name, j);
        
        temp[0] = new int[j];
        temp[1] = new int[j];

        read_h5(data_path, set_name, temp[0]);
        sprintf(set_name, "/pair_2");
        read_h5(data_path, set_name, temp[1]);

        for(m=0;m<j;m++)
        {
            all_paras->corr_cal_expo_pair_label[0][n] = temp[0][m];
            all_paras->corr_cal_expo_pair_label[1][n] = temp[1][m];
            all_paras->corr_cal_expo_pair_file_label[n] = i;
            n++;
        }

        delete[] temp[0];
        delete[] temp[1];
    }

    // prepare for jackknife 
    pre_jackknife(all_paras);
    
}

void pre_jackknife(corr_cal *all_paras)
{      
    int i, j, k, m, n;
    int *temp;

    // decide the start & end label of each sub-sample
    // each time, the expo-pairs of which the label > start and <=end
    // will be abandoned from the calculation
    temp = new int[all_paras->resample_num];
    all_paras->jackknife_expo_pair_st = new int[all_paras->resample_num+1];
    all_paras->jackknife_expo_pair_ed = new int[all_paras->resample_num+1];
    all_paras->jackknife_expo_pair_st[0] = -1;
    all_paras->jackknife_expo_pair_ed[0] = -1;

    m = all_paras->corr_cal_total_pair_num/all_paras->resample_num;
    n = all_paras->corr_cal_total_pair_num%all_paras->resample_num;
    for(i=0; i<all_paras->resample_num; i++)
    {
        temp[i] = m;
        if(i<n){temp[i] += 1;}
    }
    m = 0;
    for(i=0; i<all_paras->resample_num; i++)
    {
        for(j=0;j<i;j++){m += temp[i];}
        all_paras->jackknife_expo_pair_st[i+1] = m;
        all_paras->jackknife_expo_pair_ed[i+1] = m + temp[i];
    }

    delete[] temp;


    // distribute the resample task to each thread
    // there're "resample_num+1" tasks for "corr_cal_thread_num" CPUs
    temp = new int[all_paras->corr_cal_thread_num];

    m = (all_paras->resample_num+1)/all_paras->corr_cal_thread_num;
    n = (all_paras->resample_num+1)%all_paras->corr_cal_thread_num;

    for(i=0; i<all_paras->corr_cal_thread_num; i++)
    {
        temp[i] = m;
        if(i<n){temp[i] += 1;}
    }
    m = 0;
    for(i=0; i<all_paras->corr_cal_rank; i++)
    {
        for(j=0;j<i;j++){m += temp[i];}
    }
    all_paras->my_resample_label_st = m;
    all_paras->my_resample_label_ed = m + temp[all_paras->corr_cal_rank]; 
    
    for(j=0; j<temp[all_paras->corr_cal_rank]; j++)
    {
        
        all_paras->corr_cal_chi_tt[j] = new double[all_paras->corr_cal_chi_num];
        all_paras->corr_cal_chi_xx[j] = new double[all_paras->corr_cal_chi_num];

        all_paras->corr_cal_gtt[j]  = new double[all_paras->corr_cal_final_data_num];
        all_paras->corr_cal_gxx[j]  = new double[all_paras->corr_cal_final_data_num];
        all_paras->corr_cal_gtt_sig[j]  = new double[all_paras->corr_cal_final_data_num];
        all_paras->corr_cal_gxx_sig[j]  = new double[all_paras->corr_cal_final_data_num];

        all_paras->corr_cal_mean_theta[j] = new double[all_paras->expo_chi_block_len_true]{};
    }

    delete[] temp;
}

void resample_jackknife(corr_cal *all_paras,int resample_label)
{
    int i,j,k;
    int file_tag;
    int abort_st, abort_ed;
    char set_name[60], data_path[600];
    initialize_arr(all_paras->corr_cal_stack_num_count_chit, all_paras->expo_chi_block_len_true, 0);
    initialize_arr(all_paras->corr_cal_stack_num_count_chix, all_paras->expo_chi_block_len_true, 0);
    initialize_arr(all_paras->corr_cal_stack_expo_theta_accum, all_paras->theta_accum_len_true, 0);
    initialize_arr(all_paras->corr_cal_stack_expo_theta_num_accum, all_paras->theta_accum_len_true, 0);

    abort_st = all_paras->jackknife_expo_pair_st[resample_label];
    abort_ed = all_paras->jackknife_expo_pair_ed[resample_label];

    for(i=0; i<all_paras->corr_cal_total_pair_num; i++)
    {   
        // read the data
        if(i>= abort_st and i < abort_ed){continue;}
        file_tag = all_paras->corr_cal_expo_pair_file_label[i];

        sprintf(data_path, "%s/result/core_%d_num_count.hdf5", all_paras->parent_path,file_tag);

        sprintf(set_name, "/t");
        read_h5(data_path, set_name, all_paras->expo_num_count_chit);
        sprintf(set_name, "/x");
        read_h5(data_path, set_name, all_paras->expo_num_count_chix);
        sprintf(set_name, "/theta");
        read_h5(data_path, set_name, all_paras->theta_accum);
        sprintf(set_name, "/theta_num");
        read_h5(data_path, set_name, all_paras->theta_num_accum);

        // stack 
        for(j=0; j<all_paras->expo_chi_block_len_true; j++)
        {
            all_paras->corr_cal_stack_num_count_chit[j] += all_paras->expo_num_count_chit[j];
            all_paras->corr_cal_stack_num_count_chix[j] += all_paras->expo_num_count_chix[j];
        }
        for(j=0; j<all_paras->theta_accum_len_true; j++)
        {
            all_paras->corr_cal_stack_expo_theta_accum[j] += all_paras->theta_accum[j];
            all_paras->corr_cal_stack_expo_theta_num_accum[j] += all_paras->theta_num_accum[j];
        }
    }    
}

void chisq_2d(double *num_count, int mg_bin_num, double &chisq)
{
    int i, j, k;
    int tag[4];
    int bin_num2 = mg_bin_num/2;
    double n1, n2, n;
    n1 = 0;
    n2 = 0;
    n = 0;

    for(i=bin_num2; i<mg_bin_num;i++)
    {   // row

        for(j=bin_num2; j<mg_bin_num; j++)
        {   // col

            tag[0] = i*mg_bin_num + j;// [i,j]
            tag[1] = i*mg_bin_num + mg_bin_num - j - 1;//[i,-j]
            tag[2] = (mg_bin_num-i-1)*mg_bin_num + mg_bin_num - j - 1;//[-i,-j]
            tag[3] = (mg_bin_num-i-1)*mg_bin_num + j;//[-i,j]

            n1 = num_count[tag[0]] + num_count[tag[1]] + num_count[tag[2]] + num_count[tag[3]];
            n2 = num_count[tag[0]] + num_count[tag[2]] - num_count[tag[1]] - num_count[tag[3]];
            n += n2*n2/n1;
        }
    }
    chisq = n/2;
}

void corr_calculate(corr_cal *all_paras, int resample_label)
{   
    int i, j, k, tag;
    int mg_bin_num = all_paras->mg_bin_num;
    int nz = (all_paras->zbin_num*all_paras->zbin_num + all_paras->zbin_num)/2;

    double chisq_tt, chisq_xx;
    double *temp_tt = new double[mg_bin_num*mg_bin_num]{};
    double *temp_xx = new double[mg_bin_num*mg_bin_num]{};

    double *chi_gtt_fit = new double[all_paras->chi_guess_num];
    double *chi_gxx_fit = new double[all_paras->chi_guess_num];
    double gh, gh_sig, chi_min_fit;

    // calculate chi squared
    for(i=0;i<all_paras->corr_cal_chi_num;i++)
    {
        for(j=0;j<mg_bin_num*mg_bin_num;j++)
        {
            tag = i*mg_bin_num*mg_bin_num + j;
            temp_tt[j] = all_paras->corr_cal_stack_num_count_chit[tag];
            temp_xx[j] = all_paras->corr_cal_stack_num_count_chix[tag];
        }
        
        chisq_2d(temp_tt,mg_bin_num, chisq_tt);
        chisq_2d(temp_xx,mg_bin_num, chisq_xx);

        all_paras->corr_cal_chi_tt[resample_label][i] = chisq_tt;
        all_paras->corr_cal_chi_xx[resample_label][i] = chisq_xx;
    }
    // fitting
    for(i=0; i<all_paras->corr_cal_final_data_num;i++)
    {   
        tag = i*all_paras->chi_guess_num;
        // std::cout<<i<<std::endl;
        for(j=0;j<all_paras->chi_guess_num;j++)
        {
            chi_gtt_fit[j] = all_paras->corr_cal_chi_tt[resample_label][tag + j];
            chi_gxx_fit[j] = all_paras->corr_cal_chi_xx[resample_label][tag + j];
        }
        if(resample_label == 1)
        {
            show_arr(chi_gtt_fit,1,all_paras->chi_guess_num);
            show_arr(chi_gxx_fit,1,all_paras->chi_guess_num);
        }

        fit_shear(all_paras->corr_cal_chi_guess, chi_gtt_fit, all_paras->chi_guess_num, gh, gh_sig, chi_min_fit, 150);
        all_paras->corr_cal_gtt[resample_label][i] = gh;
        all_paras->corr_cal_gtt_sig[resample_label][i] = gh_sig;

        fit_shear(all_paras->corr_cal_chi_guess, chi_gxx_fit, all_paras->chi_guess_num, gh, gh_sig, chi_min_fit, 150);
        all_paras->corr_cal_gxx[resample_label][i] = gh;
        all_paras->corr_cal_gxx_sig[resample_label][i] = gh_sig;
    }


    // calculate the mean theta
    for(i=0;i<all_paras->theta_accum_len_true;i++)
    {
        all_paras->corr_cal_mean_theta[resample_label][i] = all_paras->corr_cal_stack_expo_theta_accum[i]/all_paras->corr_cal_stack_expo_theta_num_accum[i];
        std::cout<<all_paras->corr_cal_stack_expo_theta_accum[i]<<" "<<all_paras->corr_cal_stack_expo_theta_num_accum[i]<<std::endl;
    }
    
    std::cout<<"chi block num: "<<all_paras->corr_cal_chi_num<<" Final data point num: "<<all_paras->corr_cal_final_data_num<<" Theta point num: "<<all_paras->theta_accum_len_true<<std::endl;
    delete[] temp_tt;
    delete[] temp_xx;
    delete[] chi_gtt_fit;
    delete[] chi_gxx_fit;

}



int main(int argc, char **argv)
{
    int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    char source_list[300];
    char parent_path[400], cat_path[450], result_path[450], inform[400];
    char *result_file_path[4000];

    int i, j,k;
    corr_cal all_paras;


    strcpy(all_paras.parent_path, argv[1]);
    all_paras.resample_num = atoi(argv[2]);
    all_paras.corr_cal_result_file_num = atoi(argv[3]);
    all_paras.corr_cal_thread_num = numprocs;
    all_paras.corr_cal_rank = rank;


    sprintf(inform,"Reading parameters");
    if(rank == 0){std::cout<<inform<<std::endl;}
    read_para(&all_paras);

    prepare_data(&all_paras);


    sprintf(inform,"Calculate");
    if(rank == 0){std::cout<<inform<<std::endl;}

    for(i=all_paras.my_resample_label_st; i<all_paras.my_resample_label_ed; i++)
    {
        std::cout<<"Sampling "<<i<<std::endl;
        resample_jackknife(&all_paras,i);
        corr_calculate(&all_paras, i);
    }

    // double chisq_test;
    // int nx;
    // nx = atoi(argv[1]);

    // sprintf(parent_path,"test.hdf5");
    // sprintf(cat_path,"/data");
    // double *num = new double[nx*nx];
    // read_h5(parent_path, cat_path,num);
    // // show_arr(num, nx,nx);

    // chisq_2d(num, nx, chisq_test);
    // std::cout<<chisq_test<<std::endl;

    MPI_Finalize();
    return 0;
}