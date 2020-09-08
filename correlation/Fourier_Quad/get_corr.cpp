#include<hk_mpi.h>
#include<functions_expo_wise.h>
#include<ctime>
#include<vector>


void read_result_list(data_info *all_paras, int read_col_idx)
{
    std::ifstream infile;
	std::string str, str_;
	std::stringstream strs;
    char temp[50];
    char temp_path[600];
    char cat_path[600];
    char result_path[600];

    sprintf(cat_path, "%s/cata/source_list.dat", all_paras->parent_path);
    sprintf(result_path, "%s/result", all_paras->parent_path);


	int i, j;

    int line_count;

    infile.open(cat_path);
            
    strs << str;

	line_count = 0;

    all_paras->total_expo_num = 0;
    while (!infile.eof())
    {
        str.clear();
        strs.clear();
        getline(infile, str);
                
        strs << str;
        i = 0;

        while(strs >> str_)
        {   
            if(i == read_col_idx)
            {      
                strcpy(temp, str_.c_str());

                sprintf(temp_path,"%s/%s_num_count.hdf5", result_path, temp);
                if(file_exist(temp_path))
                {   
                    all_paras->expo_name_path[all_paras->total_expo_num] = new char[600];

                    sprintf(all_paras->expo_name_path[all_paras->total_expo_num], "%s", temp_path);

                    all_paras->total_expo_num++;
                }
                else
                {
                    std::cout<<"Cannot find "<<temp_path<<std::endl;
                }                
            }
            str_.clear();

            i++;
        }
        line_count ++;
    }
    std::cout<<"All files: "<<line_count-1<<" Read "<<all_paras->total_expo_num<<" files"<<std::endl;
    std::cout<<std::endl;
}

void read_para(data_info *all_paras)
{   
    int j;
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
    read_h5(data_path, set_name, all_paras->chi_guess);

    sprintf(set_name, "/mg_bin");
    read_h5_datasize(data_path, set_name, j);
    all_paras->mg_bin = new MY_FLOAT[j];
    all_paras->mg_bin_num = j -1;
    read_h5(data_path, set_name, all_paras->mg_bin);

    std::cout<<"G bin: "<<all_paras->mg_bin_num<<std::endl;
    show_arr(all_paras->mg_bin,1,all_paras->mg_bin_num+1);
    std::cout<<std::endl;

    std::cout<<"Chi guess bin: "<<all_paras->chi_guess_num<<std::endl;
    show_arr(all_paras->chi_guess,1,all_paras->chi_guess_num);
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


    // tangential and cross components
    all_paras->expo_num_count_chit = new double[all_paras->expo_chi_block_len]{};
    all_paras->expo_num_count_chix = new double[all_paras->expo_chi_block_len]{};

    all_paras->theta_accum_len = all_paras->theta_bin_num*all_paras->zbin_num*all_paras->zbin_num;
    all_paras->theta_accum_len_true = all_paras->theta_bin_num*((all_paras->zbin_num*all_paras->zbin_num + all_paras->zbin_num)/2);
    all_paras->theta_accum = new double[all_paras->theta_accum_len]{};
    all_paras->theta_num_accum = new double[all_paras->theta_accum_len]{};


    // stack all the exposure data
    all_paras->corr_cal_stack_num_count_chit = new double [all_paras->expo_chi_block_len_true]{};
    all_paras->corr_cal_stack_num_count_chix = new double [all_paras->expo_chi_block_len_true]{};
    all_paras->corr_cal_mean_theta = new double[all_paras->expo_chi_block_len_true]{};

    all_paras->corr_cal_stack_expo_theta_accum = new double[all_paras->theta_accum_len_true]{};
    all_paras->corr_cal_stack_expo_theta_num_accum = new double[all_paras->theta_accum_len_true]{};


    // for the PDF calculation
    all_paras->corr_cal_chi_len = all_paras->chi_guess_num*all_paras->theta_bin_num*((all_paras->zbin_num*all_paras->zbin_num + all_paras->zbin_num)/2);
    all_paras->corr_cal_final_data_num = all_paras->theta_bin_num*((all_paras->zbin_num*all_paras->zbin_num + all_paras->zbin_num)/2);
    for(j=0; j<all_paras->resample_num+1; j++)
    {
        
        all_paras->corr_cal_chi_tt[j] = new double[all_paras->corr_cal_chi_len];
        all_paras->corr_cal_chi_xx[j] = new double[all_paras->corr_cal_chi_len];

        all_paras->corr_cal_gtt[j]  = new double[all_paras->corr_cal_final_data_num];
        all_paras->corr_cal_gxx[j]  = new double[all_paras->corr_cal_final_data_num];
        all_paras->corr_cal_gtt_sig[j]  = new double[all_paras->corr_cal_final_data_num];
        all_paras->corr_cal_gxx_sig[j]  = new double[all_paras->corr_cal_final_data_num];
    }

}

void read_result_data(data_info*all_paras)
{    
    char set_name[60], data_path[600];
    int i, j, k, tag;
    int st, st_ij, st_ji, m;
    
    for(i=0;i<all_paras->total_expo_num; i++)
    {   

        all_paras->corr_cal_expo_num_count_chit[i] = new double [all_paras->expo_chi_block_len_true];
        all_paras->corr_cal_expo_num_count_chix[i] = new double [all_paras->expo_chi_block_len_true];
        all_paras->corr_cal_expo_theta_accum[i] = new double[all_paras->theta_accum_len_true];
        all_paras->corr_cal_expo_theta_num_accum[i] = new double[all_paras->theta_accum_len_true];

        sprintf(set_name, "/t");
        read_h5(all_paras->expo_name_path[i], set_name, all_paras->expo_num_count_chit);
        sprintf(set_name, "/x");
        read_h5(all_paras->expo_name_path[i], set_name, all_paras->expo_num_count_chix);
        sprintf(set_name, "/theta");
        read_h5(all_paras->expo_name_path[i], set_name, all_paras->theta_accum);
        sprintf(set_name, "/theta_num");
        read_h5(all_paras->expo_name_path[i], set_name, all_paras->theta_num_accum);
        // std::cout<<all_paras->expo_name_path[i]<<std::endl;


        for(j=0; j<all_paras->zbin_num; j++)
        {   
            //////////////////////////////////  theta  //////////////////////////////////////
            // z[i,i] part
            st = (j*all_paras->zbin_num + j)*all_paras->theta_bin_num;
            tag = (j*all_paras->zbin_num - j*(j-1)/2)*all_paras->theta_bin_num;
            // if(i == 0 or i == 10)
            // {std::cout<<j<<" "<<tag<<" "<<st<<std::endl;}

            for(m=0; m<all_paras->theta_bin_num; m++)
            {
                all_paras->corr_cal_expo_theta_accum[i][tag+m] = all_paras->theta_accum[st+m];
                all_paras->corr_cal_expo_theta_num_accum[i][tag+m] = all_paras->theta_num_accum[st+m];
            }

            // z[i,j] part, i != j, z[j,i] will be added to z[i,j], for j>i
            for(k=j+1; k<all_paras->zbin_num; k++)
            {   
                tag = (j*all_paras->zbin_num + k - (j*j+j)/2)*all_paras->theta_bin_num;
                if(i == 0 or i == 10)
                {std::cout<<j<<" "<<k<<" "<<tag<<" "<<std::endl;}

                st_ij = (j*all_paras->zbin_num + k)*all_paras->theta_bin_num;
                st_ji = (k*all_paras->zbin_num + j)*all_paras->theta_bin_num;

                for(m=0; m<all_paras->theta_bin_num; m++)
                {
                    all_paras->corr_cal_expo_theta_accum[i][tag+m] = all_paras->theta_accum[st_ij+m]+all_paras->theta_accum[st_ji+m];
                    all_paras->corr_cal_expo_theta_num_accum[i][tag+m] = all_paras->theta_num_accum[st_ij+m]+all_paras->theta_num_accum[st_ji+m];
                }
            }
            //////////////////////////////////////////////////////////////////////////////////


            ///////////////////////////////////  chi  ////////////////////////////////////////
            // z[i,i] part
            st = (j*all_paras->zbin_num + j)*all_paras->iz_chi_block_len;
            tag = (j*all_paras->zbin_num - j*(j-1)/2)*all_paras->iz_chi_block_len;
            // if(i == 0 or i == 10)
            // {std::cout<<j<<" "<<tag<<" "<<st<<std::endl;}

            for(m=0; m<all_paras->iz_chi_block_len; m++)
            {
                all_paras->corr_cal_expo_num_count_chit[i][tag+m] = all_paras->expo_num_count_chit[st+m];
                all_paras->corr_cal_expo_num_count_chix[i][tag+m] = all_paras->expo_num_count_chix[st+m];
            }

            // z[i,j] part, i != j, z[j,i] will be added to z[i,j], for j>i
            for(k=j+1; k<all_paras->zbin_num; k++)
            {   
                tag = (j*all_paras->zbin_num + k - (j*j+j)/2)*all_paras->iz_chi_block_len;
                if(i == 0 or i == 10)
                {std::cout<<j<<" "<<k<<" "<<tag<<" "<<std::endl;}

                st_ij = (j*all_paras->zbin_num + k)*all_paras->iz_chi_block_len;
                st_ji = (k*all_paras->zbin_num + j)*all_paras->iz_chi_block_len;

                for(m=0; m<all_paras->iz_chi_block_len; m++)
                {
                    all_paras->corr_cal_expo_num_count_chit[i][tag+m] = all_paras->expo_num_count_chit[st_ij+m]+all_paras->expo_num_count_chit[st_ji+m];
                    all_paras->corr_cal_expo_num_count_chix[i][tag+m] = all_paras->expo_num_count_chix[st_ij+m]+all_paras->expo_num_count_chix[st_ji+m];
                }
            }
            //////////////////////////////////////////////////////////////////////////////////
        }

        // stack all the exposure data
        for(j=0;j<all_paras->expo_chi_block_len_true;j++)
        {
            all_paras->corr_cal_stack_num_count_chit[j] += all_paras->corr_cal_expo_num_count_chit[i][j];
            all_paras->corr_cal_stack_num_count_chix[j] += all_paras->corr_cal_expo_num_count_chix[i][j];

        }
        for(j=0;j<all_paras->theta_accum_len_true;j++)
        {
            all_paras->corr_cal_stack_expo_theta_accum[j] += all_paras->corr_cal_expo_theta_accum[i][j];
            all_paras->corr_cal_stack_expo_theta_num_accum[j] += all_paras->corr_cal_expo_theta_num_accum[i][j];
        }
    }
}

void resample(data_info *all_paras,int resample_label)
{
    ;
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

void corr_calculate(data_info *all_paras, int resample_label)
{   
    int i, j, k, tag;
    int mg_bin_num = all_paras->mg_bin_num;
    int nz = (all_paras->zbin_num*all_paras->zbin_num + all_paras->zbin_num)/2;

    double chisq_tt, chisq_xx;
    double *temp_tt = new double[mg_bin_num*mg_bin_num]{};
    double *temp_xx = new double[mg_bin_num*mg_bin_num]{};
    double *chi_gtt_fit = new double[all_paras->chi_guess_num];
    double *chi_gxx_fit = new double[all_paras->chi_guess_num];

    // calculate chi squared
    for(i=0;i<all_paras->corr_cal_chi_len;i++)
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
        for(j=0;j<all_paras->chi_guess_num;j++)
        {
            chi_gtt_fit[j] = all_paras->corr_cal_chi_tt[resample_label][tag + j];
            chi_gxx_fit[j] = all_paras->corr_cal_chi_xx[resample_label][tag + j];
        }
    }


    // calculate the mean theta
    for(i=0;i<all_paras->theta_accum_len_true;i++)
    {
        all_paras->corr_cal_mean_theta[i] = all_paras->corr_cal_stack_expo_theta_accum[i]/all_paras->corr_cal_stack_expo_theta_num_accum[i];
    }

    std::cout<<all_paras->corr_cal_chi_len<<" "<<all_paras->corr_cal_final_data_num<<" "<<all_paras->theta_accum_len_true<<std::endl;
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
    char parent_path[400], cat_path[450], result_path[450];
    char *result_file_path[4000];

    int result_file_num;
    int i, j,k;
    data_info all_paras;

    int resample_num;
    
    strcpy(all_paras.parent_path, argv[1]);
    all_paras.resample_num = atoi(argv[2]);
    
    read_result_list(&all_paras, 1);

    read_para(&all_paras);
    read_result_data(&all_paras);


    for(i=0; i<all_paras.resample_num+1; i++)
    {
        // 0 means the result will be stored in the first row,
        // the signal from the whole sample
        // else the i'th row is the result from i'th jackknife or bootstrap
        if(i == 0){corr_calculate(&all_paras, i);}
        // else
        // {
        //     resample(&all_paras,i-1);
        //     corr_calculate(&all_paras, i);
        // }
        
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