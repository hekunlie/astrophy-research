#include"hk_correlation_functions_expo_wise.h"

void initialize(data_info *expo_info)
{   
    int i, j, k;
    char set_name[100];
    char data_path[600];
    char log_informs[200];

    expo_info->mg1_idx = 0;
    expo_info->mg2_idx = 1;
    expo_info->mn_idx = 2;
    expo_info->mu_idx = 3;
    expo_info->mv_idx = 4;
    expo_info->ra_idx = 5;
    expo_info->dec_idx = 6;
    expo_info->cos_dec_idx = 7;
    expo_info->redshift_idx = 8;

    sprintf(data_path,"%s/source_list.dat", expo_info->cata_path);

    line_count(data_path, expo_info);
    if(expo_info->my_rank == 0)
    {
        sprintf(log_informs,"Find %d exposures", expo_info->total_expo_num);
        std::cout<<log_informs<<std::endl;
    }

    for(i=0; i<expo_info->total_expo_num; i++)
    {   
        // the expo file directories
        expo_info->expo_name_path[i] = new char[400];
        expo_info->expo_name[i] = new char[100];
    }

    expo_info->expo_cen_ra = new MY_FLOAT[expo_info->total_expo_num]{};  
    expo_info->expo_cen_dec = new MY_FLOAT[expo_info->total_expo_num]{}; 
    expo_info->expo_delta_ra = new MY_FLOAT[expo_info->total_expo_num]{};  
    expo_info->expo_delta_dec = new MY_FLOAT[expo_info->total_expo_num]{};
    expo_info->expo_delta_len = new MY_FLOAT[expo_info->total_expo_num]{};
    expo_info->expo_cen_cos_dec = new MY_FLOAT[expo_info->total_expo_num]{};
    expo_info->expo_gal_num = new int[expo_info->total_expo_num]{};

    if(expo_info->my_rank == 0)
    {
        sprintf(log_informs,"Allocate memory for %d exposures", expo_info->total_expo_num);
        std::cout<<log_informs<<std::endl;
    }

    // read the infomation of each expo
    if(expo_info->my_rank == 0)
    {
        sprintf(log_informs,"Read informs of %d exposures", expo_info->total_expo_num);
        std::cout<<log_informs<<std::endl;
    }

    read_list(data_path, expo_info, i);

    // array length, all elements
    sprintf(set_name, "/data");
    read_h5_datasize(expo_info->expo_name_path[0], set_name, j);
    // data col = data length/gal_num, it's the same for all expos
    expo_info->expo_data_col = j / expo_info->expo_gal_num[0];


    ///////////////// read the inform of the PDF_SYM  ////////////////
    if(expo_info->my_rank == 0)
    {
        sprintf(log_informs,"Read the parameters for PDF_SYM");
        std::cout<<log_informs<<std::endl;
    }

    sprintf(data_path,"%s/gg_cor.hdf5", expo_info->cata_path);
    // read radius bin
    sprintf(set_name,"/theta_bin");
    read_h5_datasize(data_path, set_name,j);
    expo_info->theta_bin = new MY_FLOAT[j]{};
    expo_info->theta_bin_num = j -1;
    read_h5(data_path, set_name, expo_info->theta_bin);

    // read redshift bin
    sprintf(set_name,"/redshift_bin");
    read_h5_datasize(data_path, set_name,j);
    expo_info->zbin = new MY_FLOAT[j]{};
    expo_info->zbin_num = j -1;
    read_h5(data_path, set_name, expo_info->zbin);

    sprintf(set_name, "/chi_guess");
    read_h5_datasize(data_path, set_name, expo_info->chi_guess_num);
    expo_info->chi_guess = new MY_FLOAT[expo_info->chi_guess_num];
    read_h5(data_path, set_name, expo_info->chi_guess);

    sprintf(set_name, "/mg_bin");
    read_h5_datasize(data_path, set_name, j);
    expo_info->mg_bin = new MY_FLOAT[j];
    expo_info->mg_bin_num = j -1;
    expo_info->mg_bin_num2 = (j-1)/2;
    expo_info->mg_bin_num1 = expo_info->mg_bin_num2/2;
    expo_info->mg_bin_num3 = expo_info->mg_bin_num1 + expo_info->mg_bin_num2;

    read_h5(data_path, set_name, expo_info->mg_bin);
    

    if(expo_info->my_rank == 0)
    {
        sprintf(log_informs,"Allocate memory for PDF_SYM");
        std::cout<<log_informs<<std::endl;
    }
    // the fundmental block size for number counting in the PDF_SYM
    // the smallest block, mg_bin_num x mg_bin_num, for each chi guess point
    expo_info->chi_block_len = expo_info->mg_bin_num*expo_info->mg_bin_num;
    // the upper level of the above block, chi_guess_num x mg_bin_num x mg_bin_num,
    // for each theta_bin, there're ''theta_bin_num'' of such blocks in each exposure
    expo_info->ir_chi_block_len = expo_info->chi_guess_num*expo_info->chi_block_len;
    // the biggest block, for each exposure, there're ''zbin_num*zbin_num'' of such blocks in each expo
    expo_info->iz_chi_block_len = expo_info->ir_chi_block_len*expo_info->theta_bin_num;
    // tomography, z[i,j] = z[j,i], save zbin_num*zbin_num blocks for saving the time in calculation
    // [i,j] = [j,i], they will be sum after the calculation
    expo_info->expo_chi_block_len = expo_info->iz_chi_block_len*expo_info->zbin_num*expo_info->zbin_num;
    expo_info->theta_accum_len = expo_info->theta_bin_num*expo_info->zbin_num*expo_info->zbin_num;

    expo_info->expo_chi_block_len_true = expo_info->iz_chi_block_len*((expo_info->zbin_num*expo_info->zbin_num + expo_info->zbin_num)/2);
    expo_info->theta_accum_len_true = expo_info->theta_bin_num*((expo_info->zbin_num*expo_info->zbin_num + expo_info->zbin_num)/2);


    // tangential and cross components
    expo_info->expo_num_count_chit = new double[expo_info->expo_chi_block_len]{};
    expo_info->expo_num_count_chix = new double[expo_info->expo_chi_block_len]{};

    expo_info->theta_accum = new double[expo_info->theta_accum_len]{};
    expo_info->theta_num_accum = new double[expo_info->theta_accum_len]{};

    // for the data stack in the data saving step, z[i,j] = z[i, j] =  z[j, i]
    expo_info->corr_cal_stack_num_count_chit = new double [expo_info->expo_chi_block_len_true]{};
    expo_info->corr_cal_stack_num_count_chix = new double [expo_info->expo_chi_block_len_true]{};
    expo_info->corr_cal_stack_expo_theta_accum = new double[expo_info->theta_accum_len_true]{};
    expo_info->corr_cal_stack_expo_theta_num_accum = new double[expo_info->theta_accum_len_true]{};

    // read the correlated shear pairs generated before for time-saving
    sprintf(set_name,"/g11");
    read_h5_datasize(data_path, set_name, expo_info->gg_len);

    expo_info->gg_1 = new MY_FLOAT[expo_info->gg_len]{};
    expo_info->gg_2 = new MY_FLOAT[expo_info->gg_len]{};
    
    sprintf(set_name,"/g11");
    read_h5(data_path, set_name, expo_info->gg_1);
    sprintf(set_name,"/g22");
    read_h5(data_path, set_name, expo_info->gg_2);

    expo_info->loop_label = 0;
    expo_info->data_read_label_1 = 0;
    expo_info->data_read_label_2 = 0;
    expo_info->result_file_tag = 0;
    expo_info->task_complete = 0;


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // the buffer, up to 1 GB, for the data of N exposure pairs, 
    // it will be written into the result file, to save the IO in Jackknife process
    // each block (one line): theta, theta_num, chi_tt, chi_xx 
    expo_info->max_buffer_size = 2*1024*1024*(1024/8);// max element num
    expo_info->block_size_in_buffer = 2*(expo_info->expo_chi_block_len_true + expo_info->theta_accum_len_true);
    expo_info->max_block_in_buffer = expo_info->max_buffer_size/expo_info->block_size_in_buffer;
    expo_info->actual_buffer_size = expo_info->max_block_in_buffer*expo_info->block_size_in_buffer;
    expo_info->men_buffer = new double[expo_info->actual_buffer_size]{};
    // layout of jack label: [i,j,m,n.....], pair [i,j], [m,n]....
    expo_info->task_expo_pair_jack_label = new int[expo_info->max_block_in_buffer*2]{};
    expo_info->total_buffer_num = 0;
    expo_info->block_count = 0;
    /////////////////////////////////////////////////////////////////////////////////////////////////////

}

void line_count(char *file_path, data_info *expo_info)
{
    std::ifstream infile;
	std::string str;
	std::stringstream strs;

	int i, j;
    int line_count;

    infile.open(file_path);
	line_count = 0;

    while (!infile.eof())
    {
        str.clear();
        getline(infile, str);
        line_count += 1;
        // std::cout<<str<<std::endl;
    }
    infile.close();
    expo_info->total_expo_num = line_count-1;
    // std::cout<<"Total "<<line_count-1<<" lines in sources list"<<std::endl;
}

void read_list(char *file_path, data_info *expo_info, int &read_file_num)
{
    std::ifstream infile;
	std::string str;
	std::stringstream strs;

	int i, j;
    int line_count;

    infile.open(file_path);
	line_count = 0;
    while (!infile.eof())
    {
        str.clear();
        strs.clear();
        getline(infile, str);
        
        strs << str;

        strs >> expo_info->expo_name_path[line_count];
        strs >> expo_info->expo_name[line_count];

        strs >> expo_info->expo_gal_num[line_count];

        strs >> expo_info->expo_cen_ra[line_count];
        strs >> expo_info->expo_cen_dec[line_count];

        strs >> expo_info->expo_delta_ra[line_count];
        strs >> expo_info->expo_delta_dec[line_count];
        strs >> expo_info->expo_delta_len[line_count];
        strs >> expo_info->expo_cen_cos_dec[line_count];

        // std::cout << str << std::endl;
        line_count += 1;
    }
    infile.close();
}


void read_expo_data_1(data_info *expo_info,int expo_label)
{
    char set_name[100];
    if(expo_info->data_read_label_1 == 1)
    {
        delete[] expo_info->expo_data_1;
        delete[] expo_info->expo_zbin_label_1;
        delete[] expo_info->obs_expo_label_1;
        expo_info->data_read_label_1 = 0;
    }

    expo_info->expo_zbin_label_1 = new int[expo_info->expo_gal_num[expo_label]]{};
    sprintf(set_name, "/redshift_label");
    read_h5(expo_info->expo_name_path[expo_label], set_name, expo_info->expo_zbin_label_1);

    expo_info->expo_data_1 = new MY_FLOAT[expo_info->expo_gal_num[expo_label]*expo_info->expo_data_col]{};
    sprintf(set_name, "/data");
    read_h5(expo_info->expo_name_path[expo_label], set_name, expo_info->expo_data_1);

    sprintf(set_name, "/group_label");
    read_h5(expo_info->expo_name_path[expo_label], set_name, &expo_info->jack_label_1);

    expo_info->obs_expo_label_1 = new int[expo_info->expo_gal_num[expo_label]]{};
    sprintf(set_name, "/expos_label");
    read_h5(expo_info->expo_name_path[expo_label], set_name, expo_info->obs_expo_label_1);

    expo_info->data_read_label_1 = 1;
    
}

void read_expo_data_2(data_info *expo_info,int expo_label)
{
    char set_name[100];
    if(expo_info->data_read_label_2 == 1)
    {
        delete[] expo_info->expo_data_2;
        delete[] expo_info->expo_zbin_label_2;
        delete[] expo_info->obs_expo_label_2;
        expo_info->data_read_label_2 = 0;
    }   

    expo_info->expo_zbin_label_2 = new int[expo_info->expo_gal_num[expo_label]]{};
    sprintf(set_name, "/redshift_label");
    read_h5(expo_info->expo_name_path[expo_label], set_name, expo_info->expo_zbin_label_2);

    expo_info->expo_data_2 = new MY_FLOAT[expo_info->expo_gal_num[expo_label]*expo_info->expo_data_col]{};
    sprintf(set_name, "/data");
    read_h5(expo_info->expo_name_path[expo_label], set_name, expo_info->expo_data_2);

    sprintf(set_name, "/group_label");
    read_h5(expo_info->expo_name_path[expo_label], set_name, &expo_info->jack_label_2);

    expo_info->obs_expo_label_2 = new int[expo_info->expo_gal_num[expo_label]]{};
    sprintf(set_name, "/expos_label");
    read_h5(expo_info->expo_name_path[expo_label], set_name, expo_info->obs_expo_label_2);

    expo_info->data_read_label_2 = 1;
    
}


void task_prepare(int numprocs, int rank, data_info *expo_info)
{
    int i, j, k, m, n;
    // int expo_pair_count = 0;
    // k = expo_info->total_expo_num*expo_info->total_expo_num;

    
    // expo_info->expo_pair_num_each_rank = new int[numprocs]{};
    // m = expo_info->total_expo_num/numprocs;
    // n = expo_info->total_expo_num%numprocs;

    // for(i=0; i<numprocs; i++)
    // {
    //     expo_info->expo_pair_num_each_rank[i] = m;
    //     if(i<n){expo_info->expo_pair_num_each_rank[i]+=1;}
    // }
    
    // expo_info->task_expo_num = expo_info->expo_pair_num_each_rank[rank];
    // expo_info->task_expo_label = new int[expo_info->task_expo_num];
    

    // for(i=0;i<m;i++)
    // {
    //     expo_info->task_expo_label[i] = i*numprocs + rank;
    // }
    // if(rank < n)
    // {
    //     expo_info->task_expo_label[m] = m*numprocs + rank;
    // }
    // // task distribution

    // // task_distribution(numprocs, rank, expo_info);

    for(i=0; i<expo_info->total_expo_num; i++)
    {
        for(j=i; j<expo_info->total_expo_num; j++)
        {
            expo_distance(expo_info, i,j, k);
            if(k == 1)
            {
                expo_info->task_expo_pair_labels_1.push_back(i);
                expo_info->task_expo_pair_labels_2.push_back(j);
            }
        }
    }
    expo_info->task_expo_num = expo_info->task_expo_pair_labels_1.size();
    // std::cout<<expo_info->task_expo_num<<std::endl;
}

void expo_distance(data_info *expo_info, int expo_label_0, int expo_label_1, int &label)
{
    MY_FLOAT ra_1, dec_1, cos_dec_1, delta_len_1;
    MY_FLOAT ra_2, dec_2, cos_dec_2, delta_len_2;
    MY_FLOAT delta_radius, delta_ra, delta_dec;

    ra_1 = expo_info->expo_cen_ra[expo_label_0];
    dec_1 = expo_info->expo_cen_dec[expo_label_0];
    cos_dec_1 = expo_info->expo_cen_cos_dec[expo_label_0];
    delta_len_1 = expo_info->expo_delta_len[expo_label_0];

    ra_2 = expo_info->expo_cen_ra[expo_label_1];
    dec_2 = expo_info->expo_cen_dec[expo_label_1];
    cos_dec_2 = expo_info->expo_cen_cos_dec[expo_label_1];
    delta_len_2 = expo_info->expo_delta_len[expo_label_1];

    // the seperation angle (arc minute)
    delta_ra = (ra_2 - ra_1)*cos_dec_1;
    delta_dec = dec_2 - dec_1;
    delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec) - delta_len_1 - delta_len_2;   
    // std::cout<<ra_1<<" "<<dec_1<<std::endl;
    label = 0;
    if(delta_radius <= 1.2*expo_info->theta_bin[expo_info->theta_bin_num])
    {   
        // std::cout<<ra_1<<" "<<dec_1<<std::endl;
        // std::cout<<ra_2<<" "<<dec_2<<" "<<delta_radius<<" "<< expo_info->theta_bin[expo_info->theta_bin_num]<<std::endl;
        label = 1;
    }
    // std::cout<<expo_info->expo_name[expo_label_0]<<" "<<label<<" "<<expo_info->expo_name[expo_label_1]<<std::endl;
    // std::cout<<ra_1<<" "<<ra_2<<std::endl;
    // std::cout<<dec_1<<" "<<dec_2<<std::endl;
    // std::cout<<delta_radius<<" "<<sqrt(delta_ra*delta_ra + delta_dec*delta_dec)<<std::endl;
}

void initialize_thread_pool(data_info*expo_info,int numprocs)
{
    for(int i=1; i <numprocs; i++)
    {
        expo_info->thread_pool.push_back(i);
    }
}

void thread_pool_resize(data_info *expo_info)
{   
    int i, j;
    
    for(i=0; i<expo_info->thread_del.size(); i++)
    {
        for(j=0; j<expo_info->thread_pool.size(); j++)
        {
            if(expo_info->thread_del[i]==expo_info->thread_pool[i])
            {
                expo_info->thread_pool.erase(expo_info->thread_pool.begin()+j);
            }
        }
    }
    for(i=0; i<expo_info->thread_del.size(); i++){expo_info->thread_del.pop_back();}
}


void find_pairs_diff_expo_dev(data_info *expo_info, int expo_label_0, int expo_label_1)
{
    int i, j, m, n, k;
    int ig1, ig2, ig_st;

    int ir, theta_tag, theta_accum_tag, ic;
    int iz1, iz2, iz1_;

    int gal_num_1, gal_num_2;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;
    MY_FLOAT diff;
    MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
    MY_FLOAT mg1_z1_, mg2_z1_, mn_z1_, mu_z1_, mv_z1_;;

    MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;
    MY_FLOAT mg1_z2_, mg2_z2_, mn_z2_, mu_z2_,mv_z2_;
    MY_FLOAT mgt_corr_1, mgt_corr_2, mgx_corr_1, mgx_corr_2;

    int ix_tt, iy_tt, ix_xx, iy_xx;
    int pre_ix_tt, pre_iy_tt, pre_ix_xx, pre_iy_xx;

    int theta_bin_num;
    int chi_guess_num;
    int expo_chi_block_len, ir_chi_block_len,chi_block_len;
    int ir_len, ic_len;
    MY_FLOAT gg_1, gg_2, gg_len;

    MY_FLOAT delta_ra, delta_dec, delta_radius, delta_radius_check;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    double st1, st2;
    double pairs = 0;

    int loop_label;
    loop_label = expo_info->loop_label;
    loop_label = 0;
    
    int mg_bin_num, mg_bin_num2;
    mg_bin_num = expo_info->mg_bin_num;
    mg_bin_num2 = expo_info->mg_bin_num2;

    chi_guess_num = expo_info->chi_guess_num;
    expo_chi_block_len = expo_info->expo_chi_block_len;
    ir_chi_block_len = expo_info->ir_chi_block_len;
    chi_block_len = expo_info->chi_block_len;
    gg_len = expo_info->gg_len;
    theta_bin_num = expo_info->theta_bin_num;

    gal_num_1 = expo_info->expo_gal_num[expo_label_0];
    gal_num_2 = expo_info->expo_gal_num[expo_label_1];
    
    // int *mask = new int[gal_num_1]{};

    st1 = clock();
    for(ig1=0; ig1<gal_num_1; ig1++)
    {   
        
        // loop the grid in the first zbin, zbin_label_0
        m = ig1*expo_info->expo_data_col;

        ra_z1 = expo_info->expo_data_1[m+expo_info->ra_idx];
        dec_z1 = expo_info->expo_data_1[m+expo_info->dec_idx];
        cos_dec_z1 = expo_info->expo_data_1[m+expo_info->cos_dec_idx];

        mg1_z1_ = expo_info->expo_data_1[m+expo_info->mg1_idx];
        mg2_z1_ = expo_info->expo_data_1[m+expo_info->mg2_idx];

        mn_z1_ = expo_info->expo_data_1[m+expo_info->mn_idx];
        mu_z1_ = expo_info->expo_data_1[m+expo_info->mu_idx];
        mv_z1_ = expo_info->expo_data_1[m+expo_info->mv_idx];

        iz1 = expo_info->expo_zbin_label_1[ig1];
        iz1_ = iz1*expo_info->zbin_num;

        if(expo_label_0 == expo_label_1){ig_st = ig1+1;}
        else{ig_st = 0;}

        for(ig2=ig_st; ig2<gal_num_2; ig2++)
        {   
            // use the pairs from different exposures
            if(expo_info->obs_expo_label_1[ig1] == expo_info->obs_expo_label_2[ig2]){continue;}

            n = ig2*expo_info->expo_data_col;

            ra_z2 = expo_info->expo_data_2[n+expo_info->ra_idx];
            dec_z2 = expo_info->expo_data_2[n+expo_info->dec_idx];
            cos_dec_z2 = expo_info->expo_data_2[n+expo_info->cos_dec_idx];

            // the seperation angle (arc minute)
            delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
            delta_dec = dec_z2 - dec_z1;
            delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);

            // separation(ra_z1/60, dec_z1/60, ra_z2/60, dec_z2/60, delta_radius_check);
            // diff = fabs(delta_radius - delta_radius_check/Pi*180*60)/delta_radius;
            // if( diff > 0.05)
            // {   
            //     std::cout<<ra_z1<<" "<<ra_z2<<" "<<dec_z1<<" "<<dec_z2<<std::endl;
            //     std::cout<<delta_radius<<" "<<delta_radius_check/Pi*180*60<<" "<<diff<<std::endl;
            // }

            theta_tag = -1;
            for(ir=0; ir<theta_bin_num; ir++)
            {
                if(delta_radius >= expo_info->theta_bin[ir] and delta_radius < expo_info->theta_bin[ir+1])
                {theta_tag=ir;break;}
            }
            // std::cout<<delta_radius<<" "<<expo_info->theta_bin[theta_tag]<<" "<<expo_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
            if(theta_tag > -1)
            {   
                pairs+= 1;

                mg1_z2_ = expo_info->expo_data_2[n+expo_info->mg1_idx];
                mg2_z2_ = expo_info->expo_data_2[n+expo_info->mg2_idx];

                mn_z2_ = expo_info->expo_data_2[n+expo_info->mn_idx];
                mu_z2_ = expo_info->expo_data_2[n+expo_info->mu_idx];
                mv_z2_ = expo_info->expo_data_2[n+expo_info->mv_idx];

                // shear estimators rotation (position angle defined as East of North)
                sin_theta = delta_ra/delta_radius;
                cos_theta = delta_dec/delta_radius;

                sin_2theta = 2*sin_theta*cos_theta;
                cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

                sin_4theta = 2*sin_2theta*cos_2theta;
                cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;

                // rotate gal z1
                mg1_z1 = mg1_z1_*cos_2theta - mg2_z1_*sin_2theta;
                mg2_z1 = mg1_z1_*sin_2theta + mg2_z1_*cos_2theta;

                mnu1_z1 = mu_z1_*cos_4theta - mv_z1_*sin_4theta;
                mnu2_z1 = mn_z1_ - mnu1_z1;
                mnu1_z1 = mn_z1_ + mnu1_z1;
                // rotate gal z2
                mg1_z2 = mg1_z2_*cos_2theta - mg2_z2_*sin_2theta;
                mg2_z2 = mg1_z2_*sin_2theta + mg2_z2_*cos_2theta;

                mnu1_z2 = mu_z2_*cos_4theta - mv_z2_*sin_4theta;
                mnu2_z2 = mn_z2_ - mnu1_z2;
                mnu1_z2 = mn_z2_ + mnu1_z2;
                
                // there're zbin_num *zbin_num blocks, iz1 is row, iz2 is the col, each block
                // has a length of mg_bin_num*mg_bin_num*chi_guess_num*theta_bin_num.
                iz2 = expo_info->expo_zbin_label_2[ig2];

                // record the pair separation for the mean separation in that bin
                // which will be the x postion in the last figure
                theta_accum_tag = (iz1_ + iz2)*theta_bin_num + theta_tag;
        
                expo_info->theta_accum[theta_accum_tag] += delta_radius;
                expo_info->theta_num_accum[theta_accum_tag] += 1;


                ////////////////////// the key part of PDF_SYM //////////////////////////////
                ic_len = theta_tag*ir_chi_block_len + (iz1_ + iz2)*expo_info->iz_chi_block_len;

                pre_ix_tt = mg_bin_num2;
                pre_ix_xx = mg_bin_num2;
                pre_iy_tt = mg_bin_num2;
                pre_iy_xx = mg_bin_num2;

                for(ic=0; ic<chi_guess_num; ic++)
                {   
                    gg_1 = expo_info->gg_1[loop_label];
                    gg_2 = expo_info->gg_2[loop_label];

                    mgt_corr_1 = mg1_z1 - gg_1*mnu1_z1; 
                    mgt_corr_2 = mg1_z2 - gg_2*mnu1_z2;   
  
                    mgx_corr_1 = mg2_z1 - gg_1*mnu2_z1; 
                    mgx_corr_2 = mg2_z2 - gg_2*mnu2_z2;  

                    hist2d_fast_dev(mgt_corr_1, mgt_corr_2, expo_info->mg_bin, mg_bin_num, pre_ix_tt, ix_tt, pre_iy_tt, iy_tt);

                    expo_info->expo_num_count_chit[ic_len + iy_tt*mg_bin_num+ix_tt] += 1;
                    
                    pre_ix_tt = ix_tt;
                    pre_iy_tt = iy_tt;

                    hist2d_fast_dev(mgx_corr_1, mgx_corr_2, expo_info->mg_bin, mg_bin_num,pre_ix_xx, ix_xx, pre_iy_xx, iy_xx);

                    expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx] += 1;

                    pre_ix_xx = ix_xx;
                    pre_iy_xx = iy_xx;

                    loop_label += 1;

                    ic_len += chi_block_len;
                    // std::cout<<ic<<" "<<temp_tt[2]<<" "<<temp_tt[3]<<" "<<ix_tt<<" "<<iy_tt<<" "<<gg_1<<std::endl;
                    // std::cout<<ic<<" "<<temp_xx[2]<<" "<<temp_xx[3]<<" "<<ix_xx<<" "<<iy_xx<<" "<<gg_2<<std::endl;
                    // if(ic_len == 0)
                    // {std::cout<<theta_tag<<" "<<iz1<<" "<<iz2<<" "<<iy_xx<<" "<<ix_xx<<" "<<iy_xx*mg_bin_num+ix_xx<<" "<<
                    // mg_bin_num<<" "<<expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx]<<std::endl;}
                    
                }
                if(loop_label >= gg_len){loop_label = 0;}
                ////////////////////// the key part of PDF_SYM -end  //////////////////////////////
                
            }

        }
    }
    // st2 = clock();
    // std::cout<<pairs<<" pairs "<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
    expo_info->gg_pairs = pairs;
    expo_info->loop_label = loop_label;
}


void find_pairs_same_expo_dev(data_info *expo_info, int expo_label_0, int expo_label_1)
{
    int i, j, m, n, k;
    int ig1, ig2, ig_st;

    int ir, theta_tag, theta_accum_tag, ic;
    int iz1, iz2, iz1_;

    int gal_num_1, gal_num_2;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;
    MY_FLOAT diff;
    MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
    MY_FLOAT mg1_z1_, mg2_z1_, mn_z1_, mu_z1_, mv_z1_;;

    MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;
    MY_FLOAT mg1_z2_, mg2_z2_, mn_z2_, mu_z2_,mv_z2_;
    MY_FLOAT mgt_corr_1, mgt_corr_2, mgx_corr_1, mgx_corr_2;

    int ix_tt, iy_tt, ix_xx, iy_xx;
    int pre_ix_tt, pre_iy_tt, pre_ix_xx, pre_iy_xx;

    int theta_bin_num;
    int chi_guess_num;
    int expo_chi_block_len, ir_chi_block_len,chi_block_len;
    int ir_len, ic_len;
    MY_FLOAT gg_1, gg_2, gg_len;

    MY_FLOAT delta_ra, delta_dec, delta_radius, delta_radius_check;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    double st1, st2;
    double pairs = 0;

    int loop_label;
    loop_label = expo_info->loop_label;

    int mg_bin_num, mg_bin_num2;
    mg_bin_num = expo_info->mg_bin_num;
    mg_bin_num2 = expo_info->mg_bin_num2;

    chi_guess_num = expo_info->chi_guess_num;
    expo_chi_block_len = expo_info->expo_chi_block_len;
    ir_chi_block_len = expo_info->ir_chi_block_len;
    chi_block_len = expo_info->chi_block_len;
    gg_len = expo_info->gg_len;
    theta_bin_num = expo_info->theta_bin_num;

    gal_num_1 = expo_info->expo_gal_num[expo_label_0];
    gal_num_2 = expo_info->expo_gal_num[expo_label_1];
    
    // int *mask = new int[gal_num_1]{};

    st1 = clock();
    for(ig1=0; ig1<gal_num_1; ig1++)
    {   
        
        // loop the grid in the first zbin, zbin_label_0
        m = ig1*expo_info->expo_data_col;

        ra_z1 = expo_info->expo_data_1[m+expo_info->ra_idx];
        dec_z1 = expo_info->expo_data_1[m+expo_info->dec_idx];
        cos_dec_z1 = expo_info->expo_data_1[m+expo_info->cos_dec_idx];

        mg1_z1_ = expo_info->expo_data_1[m+expo_info->mg1_idx];
        mg2_z1_ = expo_info->expo_data_1[m+expo_info->mg2_idx];

        mn_z1_ = expo_info->expo_data_1[m+expo_info->mn_idx];
        mu_z1_ = expo_info->expo_data_1[m+expo_info->mu_idx];
        mv_z1_ = expo_info->expo_data_1[m+expo_info->mv_idx];

        iz1 = expo_info->expo_zbin_label_1[ig1];
        iz1_ = iz1*expo_info->zbin_num;

        if(expo_label_0 == expo_label_1){ig_st = ig1+1;}
        else{ig_st = 0;}

        for(ig2=ig_st; ig2<gal_num_2; ig2++)
        {   
            // use the pairs from same exposure
            if(expo_info->obs_expo_label_1[ig1] != expo_info->obs_expo_label_2[ig2]){continue;}

            n = ig2*expo_info->expo_data_col;

            ra_z2 = expo_info->expo_data_2[n+expo_info->ra_idx];
            dec_z2 = expo_info->expo_data_2[n+expo_info->dec_idx];
            cos_dec_z2 = expo_info->expo_data_2[n+expo_info->cos_dec_idx];

            // the seperation angle (arc minute)
            delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
            delta_dec = dec_z2 - dec_z1;
            delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);

            // separation(ra_z1/60, dec_z1/60, ra_z2/60, dec_z2/60, delta_radius_check);
            // diff = fabs(delta_radius - delta_radius_check/Pi*180*60)/delta_radius;
            // if( diff > 0.05)
            // {   
            //     std::cout<<ra_z1<<" "<<ra_z2<<" "<<dec_z1<<" "<<dec_z2<<std::endl;
            //     std::cout<<delta_radius<<" "<<delta_radius_check/Pi*180*60<<" "<<diff<<std::endl;
            // }

            theta_tag = -1;
            for(ir=0; ir<theta_bin_num; ir++)
            {
                if(delta_radius >= expo_info->theta_bin[ir] and delta_radius < expo_info->theta_bin[ir+1])
                {theta_tag=ir;break;}
            }
            // std::cout<<delta_radius<<" "<<expo_info->theta_bin[theta_tag]<<" "<<expo_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
            if(theta_tag > -1)
            {   
                pairs+= 1;

                mg1_z2_ = expo_info->expo_data_2[n+expo_info->mg1_idx];
                mg2_z2_ = expo_info->expo_data_2[n+expo_info->mg2_idx];

                mn_z2_ = expo_info->expo_data_2[n+expo_info->mn_idx];
                mu_z2_ = expo_info->expo_data_2[n+expo_info->mu_idx];
                mv_z2_ = expo_info->expo_data_2[n+expo_info->mv_idx];

                // shear estimators rotation (position angle defined as East of North)
                sin_theta = delta_ra/delta_radius;
                cos_theta = delta_dec/delta_radius;

                sin_2theta = 2*sin_theta*cos_theta;
                cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

                sin_4theta = 2*sin_2theta*cos_2theta;
                cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;

                // rotate gal z1
                mg1_z1 = mg1_z1_*cos_2theta - mg2_z1_*sin_2theta;
                mg2_z1 = mg1_z1_*sin_2theta + mg2_z1_*cos_2theta;

                mnu1_z1 = mu_z1_*cos_4theta - mv_z1_*sin_4theta;
                mnu2_z1 = mn_z1_ - mnu1_z1;
                mnu1_z1 = mn_z1_ + mnu1_z1;
                // rotate gal z2
                mg1_z2 = mg1_z2_*cos_2theta - mg2_z2_*sin_2theta;
                mg2_z2 = mg1_z2_*sin_2theta + mg2_z2_*cos_2theta;

                mnu1_z2 = mu_z2_*cos_4theta - mv_z2_*sin_4theta;
                mnu2_z2 = mn_z2_ - mnu1_z2;
                mnu1_z2 = mn_z2_ + mnu1_z2;
                
                // there're zbin_num *zbin_num blocks, iz1 is row, iz2 is the col, each block
                // has a length of mg_bin_num*mg_bin_num*chi_guess_num*theta_bin_num.
                iz2 = expo_info->expo_zbin_label_2[ig2];

                // record the pair separation for the mean separation in that bin
                // which will be the x postion in the last figure
                theta_accum_tag = (iz1_ + iz2)*theta_bin_num + theta_tag;
        
                expo_info->theta_accum[theta_accum_tag] += delta_radius;
                expo_info->theta_num_accum[theta_accum_tag] += 1;


                ////////////////////// the key part of PDF_SYM //////////////////////////////
                ic_len = theta_tag*ir_chi_block_len + (iz1_ + iz2)*expo_info->iz_chi_block_len;

                pre_ix_tt = mg_bin_num2;
                pre_ix_xx = mg_bin_num2;
                pre_iy_tt = mg_bin_num2;
                pre_iy_xx = mg_bin_num2;

                for(ic=0; ic<chi_guess_num; ic++)
                {   
                    
                    gg_1 = expo_info->gg_1[loop_label];
                    gg_2 = expo_info->gg_2[loop_label];

                    mgt_corr_1 = mg1_z1 - gg_1*mnu1_z1; 
                    mgt_corr_2 = mg1_z2 - gg_2*mnu1_z2;   
  
                    mgx_corr_1 = mg2_z1 - gg_1*mnu2_z1; 
                    mgx_corr_2 = mg2_z2 - gg_2*mnu2_z2;  

                    hist2d_fast_dev(mgt_corr_1, mgt_corr_2,expo_info->mg_bin, mg_bin_num, pre_ix_tt, ix_tt, pre_iy_tt, iy_tt);

                    expo_info->expo_num_count_chit[ic_len + iy_tt*mg_bin_num+ix_tt] += 1;
                    
                    pre_ix_tt = ix_tt;
                    pre_iy_tt = iy_tt;

                    hist2d_fast_dev(mgx_corr_1, mgx_corr_2, expo_info->mg_bin, mg_bin_num,pre_ix_xx, ix_xx, pre_iy_xx, iy_xx);

                    expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx] += 1;

                    pre_ix_xx = ix_xx;
                    pre_iy_xx = iy_xx;

                    loop_label += 1;
                    ic_len += chi_block_len;

                    // std::cout<<ic<<" "<<temp_tt[2]<<" "<<temp_tt[3]<<" "<<ix_tt<<" "<<iy_tt<<" "<<gg_1<<std::endl;
                    // std::cout<<ic<<" "<<temp_xx[2]<<" "<<temp_xx[3]<<" "<<ix_xx<<" "<<iy_xx<<" "<<gg_2<<std::endl;
                    // if(ic_len == 0)
                    // {std::cout<<theta_tag<<" "<<iz1<<" "<<iz2<<" "<<iy_xx<<" "<<ix_xx<<" "<<iy_xx*mg_bin_num+ix_xx<<" "<<
                    // mg_bin_num<<" "<<expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx]<<std::endl;}
                    
                }
                if(loop_label >= gg_len){loop_label = 0;}
                ////////////////////// the key part of PDF_SYM -end  //////////////////////////////
                
            }

        }
    }
    // st2 = clock();
    // std::cout<<pairs<<" pairs "<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
    expo_info->gg_pairs = pairs;
    expo_info->loop_label = loop_label;
}



void find_pairs_stack_expo_dev(data_info *expo_info, int expo_label_0, int expo_label_1)
{
    int i, j, m, n, k;
    int ig1, ig2, ig_st;

    int ir, theta_tag, theta_accum_tag, ic;
    int iz1, iz2, iz1_;

    int gal_num_1, gal_num_2;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;
    MY_FLOAT diff;
    MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
    MY_FLOAT mg1_z1_, mg2_z1_, mn_z1_, mu_z1_, mv_z1_;;

    MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;
    MY_FLOAT mg1_z2_, mg2_z2_, mn_z2_, mu_z2_,mv_z2_;
    MY_FLOAT mgt_corr_1, mgt_corr_2, mgx_corr_1, mgx_corr_2;

    int ix_tt, iy_tt, ix_xx, iy_xx;
    int pre_ix_tt, pre_iy_tt, pre_ix_xx, pre_iy_xx;

    int theta_bin_num;
    int chi_guess_num;
    int expo_chi_block_len, ir_chi_block_len,chi_block_len;
    int ir_len, ic_len;
    MY_FLOAT gg_1, gg_2, gg_len;

    MY_FLOAT delta_ra, delta_dec, delta_radius, delta_radius_check;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    double st1, st2;
    double pairs = 0;

    int loop_label;
    loop_label = expo_info->loop_label;
    loop_label = 0;
    
    int mg_bin_num, mg_bin_num2;
    mg_bin_num = expo_info->mg_bin_num;
    mg_bin_num2 = expo_info->mg_bin_num2;

    chi_guess_num = expo_info->chi_guess_num;
    expo_chi_block_len = expo_info->expo_chi_block_len;
    ir_chi_block_len = expo_info->ir_chi_block_len;
    chi_block_len = expo_info->chi_block_len;
    gg_len = expo_info->gg_len;
    theta_bin_num = expo_info->theta_bin_num;

    gal_num_1 = expo_info->expo_gal_num[expo_label_0];
    gal_num_2 = expo_info->expo_gal_num[expo_label_1];
    
    // int *mask = new int[gal_num_1]{};

    st1 = clock();
    for(ig1=0; ig1<gal_num_1; ig1++)
    {   
        
        // loop the grid in the first zbin, zbin_label_0
        m = ig1*expo_info->expo_data_col;

        ra_z1 = expo_info->expo_data_1[m+expo_info->ra_idx];
        dec_z1 = expo_info->expo_data_1[m+expo_info->dec_idx];
        cos_dec_z1 = expo_info->expo_data_1[m+expo_info->cos_dec_idx];

        mg1_z1_ = expo_info->expo_data_1[m+expo_info->mg1_idx];
        mg2_z1_ = expo_info->expo_data_1[m+expo_info->mg2_idx];

        mn_z1_ = expo_info->expo_data_1[m+expo_info->mn_idx];
        mu_z1_ = expo_info->expo_data_1[m+expo_info->mu_idx];
        mv_z1_ = expo_info->expo_data_1[m+expo_info->mv_idx];

        iz1 = expo_info->expo_zbin_label_1[ig1];
        iz1_ = iz1*expo_info->zbin_num;

        if(expo_label_0 == expo_label_1){ig_st = ig1+1;}
        else{ig_st = 0;}

        for(ig2=ig_st; ig2<gal_num_2; ig2++)
        {   
            n = ig2*expo_info->expo_data_col;

            ra_z2 = expo_info->expo_data_2[n+expo_info->ra_idx];
            dec_z2 = expo_info->expo_data_2[n+expo_info->dec_idx];
            cos_dec_z2 = expo_info->expo_data_2[n+expo_info->cos_dec_idx];

            // the seperation angle (arc minute)
            delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
            delta_dec = dec_z2 - dec_z1;
            delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);

            // separation(ra_z1/60, dec_z1/60, ra_z2/60, dec_z2/60, delta_radius_check);
            // diff = fabs(delta_radius - delta_radius_check/Pi*180*60)/delta_radius;
            // if( diff > 0.05)
            // {   
            //     std::cout<<ra_z1<<" "<<ra_z2<<" "<<dec_z1<<" "<<dec_z2<<std::endl;
            //     std::cout<<delta_radius<<" "<<delta_radius_check/Pi*180*60<<" "<<diff<<std::endl;
            // }

            theta_tag = -1;
            for(ir=0; ir<theta_bin_num; ir++)
            {
                if(delta_radius >= expo_info->theta_bin[ir] and delta_radius < expo_info->theta_bin[ir+1])
                {theta_tag=ir;break;}
            }
            // std::cout<<delta_radius<<" "<<expo_info->theta_bin[theta_tag]<<" "<<expo_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
            if(theta_tag > -1)
            {   
                pairs+= 1;

                mg1_z2_ = expo_info->expo_data_2[n+expo_info->mg1_idx];
                mg2_z2_ = expo_info->expo_data_2[n+expo_info->mg2_idx];

                mn_z2_ = expo_info->expo_data_2[n+expo_info->mn_idx];
                mu_z2_ = expo_info->expo_data_2[n+expo_info->mu_idx];
                mv_z2_ = expo_info->expo_data_2[n+expo_info->mv_idx];

                // shear estimators rotation (position angle defined as East of North)
                sin_theta = delta_ra/delta_radius;
                cos_theta = delta_dec/delta_radius;

                sin_2theta = 2*sin_theta*cos_theta;
                cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

                sin_4theta = 2*sin_2theta*cos_2theta;
                cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;

                // rotate gal z1
                mg1_z1 = mg1_z1_*cos_2theta - mg2_z1_*sin_2theta;
                mg2_z1 = mg1_z1_*sin_2theta + mg2_z1_*cos_2theta;

                mnu1_z1 = mu_z1_*cos_4theta - mv_z1_*sin_4theta;
                mnu2_z1 = mn_z1_ - mnu1_z1;
                mnu1_z1 = mn_z1_ + mnu1_z1;
                // rotate gal z2
                mg1_z2 = mg1_z2_*cos_2theta - mg2_z2_*sin_2theta;
                mg2_z2 = mg1_z2_*sin_2theta + mg2_z2_*cos_2theta;

                mnu1_z2 = mu_z2_*cos_4theta - mv_z2_*sin_4theta;
                mnu2_z2 = mn_z2_ - mnu1_z2;
                mnu1_z2 = mn_z2_ + mnu1_z2;
                
                // there're zbin_num *zbin_num blocks, iz1 is row, iz2 is the col, each block
                // has a length of mg_bin_num*mg_bin_num*chi_guess_num*theta_bin_num.
                iz2 = expo_info->expo_zbin_label_2[ig2];

                // record the pair separation for the mean separation in that bin
                // which will be the x postion in the last figure
                theta_accum_tag = (iz1_ + iz2)*theta_bin_num + theta_tag;
        
                expo_info->theta_accum[theta_accum_tag] += delta_radius;
                expo_info->theta_num_accum[theta_accum_tag] += 1;


                ////////////////////// the key part of PDF_SYM //////////////////////////////
                ic_len = theta_tag*ir_chi_block_len + (iz1_ + iz2)*expo_info->iz_chi_block_len;

                pre_ix_tt = mg_bin_num2;
                pre_ix_xx = mg_bin_num2;
                pre_iy_tt = mg_bin_num2;
                pre_iy_xx = mg_bin_num2;

                for(ic=0; ic<chi_guess_num; ic++)
                {   
                    gg_1 = expo_info->gg_1[loop_label];
                    gg_2 = expo_info->gg_2[loop_label];

                    mgt_corr_1 = mg1_z1 - gg_1*mnu1_z1; 
                    mgt_corr_2 = mg1_z2 - gg_2*mnu1_z2;   
  
                    mgx_corr_1 = mg2_z1 - gg_1*mnu2_z1; 
                    mgx_corr_2 = mg2_z2 - gg_2*mnu2_z2;  

                    hist2d_fast_dev(mgt_corr_1, mgt_corr_2, expo_info->mg_bin, mg_bin_num, pre_ix_tt, ix_tt, pre_iy_tt, iy_tt);

                    expo_info->expo_num_count_chit[ic_len + iy_tt*mg_bin_num+ix_tt] += 1;
                    
                    pre_ix_tt = ix_tt;
                    pre_iy_tt = iy_tt;

                    hist2d_fast_dev(mgx_corr_1, mgx_corr_2, expo_info->mg_bin, mg_bin_num,pre_ix_xx, ix_xx, pre_iy_xx, iy_xx);

                    expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx] += 1;

                    pre_ix_xx = ix_xx;
                    pre_iy_xx = iy_xx;

                    loop_label += 1;

                    ic_len += chi_block_len;
                    // std::cout<<ic<<" "<<temp_tt[2]<<" "<<temp_tt[3]<<" "<<ix_tt<<" "<<iy_tt<<" "<<gg_1<<std::endl;
                    // std::cout<<ic<<" "<<temp_xx[2]<<" "<<temp_xx[3]<<" "<<ix_xx<<" "<<iy_xx<<" "<<gg_2<<std::endl;
                    // if(ic_len == 0)
                    // {std::cout<<theta_tag<<" "<<iz1<<" "<<iz2<<" "<<iy_xx<<" "<<ix_xx<<" "<<iy_xx*mg_bin_num+ix_xx<<" "<<
                    // mg_bin_num<<" "<<expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx]<<std::endl;}
                    
                }
                if(loop_label >= gg_len){loop_label = 0;}
                ////////////////////// the key part of PDF_SYM -end  //////////////////////////////
                
            }

        }
    }
    // st2 = clock();
    // std::cout<<pairs<<" pairs "<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
    expo_info->gg_pairs = pairs;
    expo_info->loop_label = loop_label;
}



void find_pairs_diff_expo(data_info *expo_info, int expo_label_0, int expo_label_1)
{
    int i, j, m, n, k;
    int ig1, ig2, ig_st;

    int ir, theta_tag, theta_accum_tag, ic;
    int iz1, iz2, iz1_;

    int gal_num_1, gal_num_2;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;
    MY_FLOAT diff;
    MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
    MY_FLOAT mg1_z1_, mg2_z1_, mn_z1_, mu_z1_, mv_z1_;;

    MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;
    MY_FLOAT mg1_z2_, mg2_z2_, mn_z2_, mu_z2_,mv_z2_;

    int ix_tt, iy_tt, ix_xx, iy_xx;
    MY_FLOAT temp_tt[4], temp_xx[4];
    int bin_para_tt[2], bin_para_xx[2];

    int theta_bin_num;
    int chi_guess_num;
    int expo_chi_block_len, ir_chi_block_len,chi_block_len;
    int ir_len, ic_len;
    MY_FLOAT gg_1, gg_2, gg_len;

    MY_FLOAT delta_ra, delta_dec, delta_radius, delta_radius_check;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    double st1, st2;
    double pairs = 0;

    int loop_label;
    loop_label = expo_info->loop_label;
    loop_label = 0; 

    int mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3;
    mg_bin_num = expo_info->mg_bin_num;
    mg_bin_num1 = expo_info->mg_bin_num1;
    mg_bin_num2 = expo_info->mg_bin_num2;
    mg_bin_num3 = expo_info->mg_bin_num3;

    chi_guess_num = expo_info->chi_guess_num;
    expo_chi_block_len = expo_info->expo_chi_block_len;
    ir_chi_block_len = expo_info->ir_chi_block_len;
    chi_block_len = expo_info->chi_block_len;
    gg_len = expo_info->gg_len;
    theta_bin_num = expo_info->theta_bin_num;

    gal_num_1 = expo_info->expo_gal_num[expo_label_0];
    gal_num_2 = expo_info->expo_gal_num[expo_label_1];
    
    // int *mask = new int[gal_num_1]{};

    st1 = clock();
    for(ig1=0; ig1<gal_num_1; ig1++)
    {   
        
        // loop the grid in the first zbin, zbin_label_0
        m = ig1*expo_info->expo_data_col;

        ra_z1 = expo_info->expo_data_1[m+expo_info->ra_idx];
        dec_z1 = expo_info->expo_data_1[m+expo_info->dec_idx];
        cos_dec_z1 = expo_info->expo_data_1[m+expo_info->cos_dec_idx];

        mg1_z1_ = expo_info->expo_data_1[m+expo_info->mg1_idx];
        mg2_z1_ = expo_info->expo_data_1[m+expo_info->mg2_idx];

        mn_z1_ = expo_info->expo_data_1[m+expo_info->mn_idx];
        mu_z1_ = expo_info->expo_data_1[m+expo_info->mu_idx];
        mv_z1_ = expo_info->expo_data_1[m+expo_info->mv_idx];

        iz1 = expo_info->expo_zbin_label_1[ig1];
        iz1_ = iz1*expo_info->zbin_num;

        if(expo_label_0 == expo_label_1){ig_st = ig1+1;}
        else{ig_st = 0;}

        for(ig2=ig_st; ig2<gal_num_2; ig2++)
        {   
            // use the pairs from different exposures
            if(expo_info->obs_expo_label_1[ig1] == expo_info->obs_expo_label_2[ig2]){continue;}

            n = ig2*expo_info->expo_data_col;

            ra_z2 = expo_info->expo_data_2[n+expo_info->ra_idx];
            dec_z2 = expo_info->expo_data_2[n+expo_info->dec_idx];
            cos_dec_z2 = expo_info->expo_data_2[n+expo_info->cos_dec_idx];

            // the seperation angle (arc minute)
            delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
            delta_dec = dec_z2 - dec_z1;
            delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);

            // separation(ra_z1/60, dec_z1/60, ra_z2/60, dec_z2/60, delta_radius_check);
            // diff = fabs(delta_radius - delta_radius_check/Pi*180*60)/delta_radius;
            // if( diff > 0.05)
            // {   
            //     std::cout<<ra_z1<<" "<<ra_z2<<" "<<dec_z1<<" "<<dec_z2<<std::endl;
            //     std::cout<<delta_radius<<" "<<delta_radius_check/Pi*180*60<<" "<<diff<<std::endl;
            // }

            theta_tag = -1;
            for(ir=0; ir<theta_bin_num; ir++)
            {
                if(delta_radius >= expo_info->theta_bin[ir] and delta_radius < expo_info->theta_bin[ir+1])
                {theta_tag=ir;break;}
            }
            // std::cout<<delta_radius<<" "<<expo_info->theta_bin[theta_tag]<<" "<<expo_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
            if(theta_tag > -1)
            {   
                pairs+= 1;

                mg1_z2_ = expo_info->expo_data_2[n+expo_info->mg1_idx];
                mg2_z2_ = expo_info->expo_data_2[n+expo_info->mg2_idx];

                mn_z2_ = expo_info->expo_data_2[n+expo_info->mn_idx];
                mu_z2_ = expo_info->expo_data_2[n+expo_info->mu_idx];
                mv_z2_ = expo_info->expo_data_2[n+expo_info->mv_idx];

                // shear estimators rotation (position angle defined as East of North)
                sin_theta = delta_ra/delta_radius;
                cos_theta = delta_dec/delta_radius;

                sin_2theta = 2*sin_theta*cos_theta;
                cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

                sin_4theta = 2*sin_2theta*cos_2theta;
                cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;

                // rotate gal z1
                mg1_z1 = mg1_z1_*cos_2theta - mg2_z1_*sin_2theta;
                mg2_z1 = mg1_z1_*sin_2theta + mg2_z1_*cos_2theta;

                mnu1_z1 = mu_z1_*cos_4theta - mv_z1_*sin_4theta;
                mnu2_z1 = mn_z1_ - mnu1_z1;
                mnu1_z1 = mn_z1_ + mnu1_z1;
                // rotate gal z2
                mg1_z2 = mg1_z2_*cos_2theta - mg2_z2_*sin_2theta;
                mg2_z2 = mg1_z2_*sin_2theta + mg2_z2_*cos_2theta;

                mnu1_z2 = mu_z2_*cos_4theta - mv_z2_*sin_4theta;
                mnu2_z2 = mn_z2_ - mnu1_z2;
                mnu1_z2 = mn_z2_ + mnu1_z2;
                
                // there're zbin_num *zbin_num blocks, iz1 is row, iz2 is the col, each block
                // has a length of mg_bin_num*mg_bin_num*chi_guess_num*theta_bin_num.
                iz2 = expo_info->expo_zbin_label_2[ig2];

                // record the pair separation for the mean separation in that bin
                // which will be the x postion in the last figure
                theta_accum_tag = (iz1_ + iz2)*theta_bin_num + theta_tag;
        
                expo_info->theta_accum[theta_accum_tag] += delta_radius;
                expo_info->theta_num_accum[theta_accum_tag] += 1;


                ////////////////////// the key part of PDF_SYM //////////////////////////////
                ic_len = theta_tag*ir_chi_block_len + (iz1_ + iz2)*expo_info->iz_chi_block_len;

                gg_1 = expo_info->gg_1[loop_label];
                gg_2 = expo_info->gg_2[loop_label];

                temp_tt[2] = mg1_z1 - gg_1*mnu1_z1;
                temp_tt[3] = mg1_z2 - gg_2*mnu1_z2;
                hist_2d_new(temp_tt[2], temp_tt[3], expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_tt, iy_tt);
                
                
                expo_info->expo_num_count_chit[ic_len + iy_tt*mg_bin_num+ix_tt] += 1;
                
                temp_xx[2] = mg2_z1 - gg_1*mnu2_z1;
                temp_xx[3] = mg2_z2 - gg_2*mnu2_z2;

                hist_2d_new(temp_xx[2], temp_xx[3],  expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_xx, iy_xx);
                expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx] += 1;
                loop_label += 1;

                // if(ic_len == 0)
                // {std::cout<<theta_tag<<" "<<iz1<<" "<<iz2<<" "<<iy_xx<<" "<<ix_xx<<" "<<iy_xx*mg_bin_num+ix_xx<<" "<<
                // mg_bin_num<<" "<<expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx]<<std::endl;}
                // std::cout<<0<<" "<<temp_tt[2]<<" "<<temp_tt[3]<<" "<<ix_tt<<" "<<iy_tt<<" "<<gg_1<<std::endl;
                // std::cout<<0<<" "<<temp_xx[2]<<" "<<temp_xx[3]<<" "<<ix_xx<<" "<<iy_xx<<" "<<gg_2<<std::endl;
                for(ic=1; ic<chi_guess_num; ic++)
                {   
                    ic_len += chi_block_len;

                    gg_1 = expo_info->gg_1[loop_label];
                    gg_2 = expo_info->gg_2[loop_label];
                                    
                    bin_para_tt[0] = ix_tt;
                    bin_para_tt[1] = iy_tt;

                    temp_tt[0] = temp_tt[2];
                    temp_tt[1] = temp_tt[3];

                    temp_tt[2] = mg1_z1 - gg_1*mnu1_z1;
                    temp_tt[3] = mg1_z2 - gg_2*mnu1_z2;
                    
                    hist_2d_new(expo_info->mg_bin, mg_bin_num, temp_tt, bin_para_tt, ix_tt, iy_tt);

                    expo_info->expo_num_count_chit[ic_len + iy_tt*mg_bin_num+ix_tt] += 1;
                    
                    bin_para_xx[0] = ix_xx;
                    bin_para_xx[1] = iy_xx;
                    
                    temp_xx[0] = temp_xx[2];
                    temp_xx[1] = temp_xx[3];

                    temp_xx[2] = mg2_z1 - gg_1*mnu2_z1;
                    temp_xx[3] = mg2_z2 - gg_2*mnu2_z2;

                
                    hist_2d_new(expo_info->mg_bin, mg_bin_num, temp_xx, bin_para_xx, ix_xx, iy_xx);
                    expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx] += 1;
                    loop_label += 1;


                    // std::cout<<ic<<" "<<temp_tt[2]<<" "<<temp_tt[3]<<" "<<ix_tt<<" "<<iy_tt<<" "<<gg_1<<std::endl;
                    // std::cout<<ic<<" "<<temp_xx[2]<<" "<<temp_xx[3]<<" "<<ix_xx<<" "<<iy_xx<<" "<<gg_2<<std::endl;
                    // if(ic_len == 0)
                    // {std::cout<<theta_tag<<" "<<iz1<<" "<<iz2<<" "<<iy_xx<<" "<<ix_xx<<" "<<iy_xx*mg_bin_num+ix_xx<<" "<<
                    // mg_bin_num<<" "<<expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx]<<std::endl;}
                    
                }
                if(loop_label >= gg_len){loop_label = 0;}
                //////////////////// the key part of PDF_SYM -end  //////////////////////////////
                
            }

        }
    }
    // st2 = clock();
    // std::cout<<pairs<<" pairs "<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
    expo_info->gg_pairs = pairs;
    expo_info->loop_label = loop_label;
}


void find_pairs_same_expo(data_info *expo_info, int expo_label_0, int expo_label_1)
{
    int i, j, m, n, k;
    int ig1, ig2, ig_st;

    int ir, theta_tag, theta_accum_tag, ic;
    int iz1, iz2, iz1_;

    int gal_num_1, gal_num_2;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;
    MY_FLOAT diff;
    MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
    MY_FLOAT mg1_z1_, mg2_z1_, mn_z1_, mu_z1_, mv_z1_;;

    MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;
    MY_FLOAT mg1_z2_, mg2_z2_, mn_z2_, mu_z2_,mv_z2_;

    int ix_tt, iy_tt, ix_xx, iy_xx;
    MY_FLOAT temp_tt[4], temp_xx[4];
    int bin_para_tt[2], bin_para_xx[2];

    int theta_bin_num;
    int chi_guess_num;
    int expo_chi_block_len, ir_chi_block_len,chi_block_len;
    int ir_len, ic_len;
    MY_FLOAT gg_1, gg_2, gg_len;

    MY_FLOAT delta_ra, delta_dec, delta_radius, delta_radius_check;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    double st1, st2;
    double pairs = 0;

    int loop_label;
    loop_label = expo_info->loop_label;

    int mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3;
    mg_bin_num = expo_info->mg_bin_num;
    mg_bin_num1 = expo_info->mg_bin_num1;
    mg_bin_num2 = expo_info->mg_bin_num2;
    mg_bin_num3 = expo_info->mg_bin_num3;

    chi_guess_num = expo_info->chi_guess_num;
    expo_chi_block_len = expo_info->expo_chi_block_len;
    ir_chi_block_len = expo_info->ir_chi_block_len;
    chi_block_len = expo_info->chi_block_len;
    gg_len = expo_info->gg_len;
    theta_bin_num = expo_info->theta_bin_num;

    gal_num_1 = expo_info->expo_gal_num[expo_label_0];
    gal_num_2 = expo_info->expo_gal_num[expo_label_1];
    
    // int *mask = new int[gal_num_1]{};

    st1 = clock();
    for(ig1=0; ig1<gal_num_1; ig1++)
    {   
        
        // loop the grid in the first zbin, zbin_label_0
        m = ig1*expo_info->expo_data_col;

        ra_z1 = expo_info->expo_data_1[m+expo_info->ra_idx];
        dec_z1 = expo_info->expo_data_1[m+expo_info->dec_idx];
        cos_dec_z1 = expo_info->expo_data_1[m+expo_info->cos_dec_idx];

        mg1_z1_ = expo_info->expo_data_1[m+expo_info->mg1_idx];
        mg2_z1_ = expo_info->expo_data_1[m+expo_info->mg2_idx];

        mn_z1_ = expo_info->expo_data_1[m+expo_info->mn_idx];
        mu_z1_ = expo_info->expo_data_1[m+expo_info->mu_idx];
        mv_z1_ = expo_info->expo_data_1[m+expo_info->mv_idx];

        iz1 = expo_info->expo_zbin_label_1[ig1];
        iz1_ = iz1*expo_info->zbin_num;

        if(expo_label_0 == expo_label_1){ig_st = ig1+1;}
        else{ig_st = 0;}

        for(ig2=ig_st; ig2<gal_num_2; ig2++)
        {   
            
            // use the pairs from the same exposure
            if(expo_info->obs_expo_label_1[ig1] != expo_info->obs_expo_label_2[ig2]){continue;}

            n = ig2*expo_info->expo_data_col;

            ra_z2 = expo_info->expo_data_2[n+expo_info->ra_idx];
            dec_z2 = expo_info->expo_data_2[n+expo_info->dec_idx];
            cos_dec_z2 = expo_info->expo_data_2[n+expo_info->cos_dec_idx];

            // the seperation angle (arc minute)
            delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
            delta_dec = dec_z2 - dec_z1;
            delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);

            // separation(ra_z1/60, dec_z1/60, ra_z2/60, dec_z2/60, delta_radius_check);
            // diff = fabs(delta_radius - delta_radius_check/Pi*180*60)/delta_radius;
            // if( diff > 0.05)
            // {   
            //     std::cout<<ra_z1<<" "<<ra_z2<<" "<<dec_z1<<" "<<dec_z2<<std::endl;
            //     std::cout<<delta_radius<<" "<<delta_radius_check/Pi*180*60<<" "<<diff<<std::endl;
            // }

            theta_tag = -1;
            for(ir=0; ir<theta_bin_num; ir++)
            {
                if(delta_radius >= expo_info->theta_bin[ir] and delta_radius < expo_info->theta_bin[ir+1])
                {theta_tag=ir;break;}
            }
            // std::cout<<delta_radius<<" "<<expo_info->theta_bin[theta_tag]<<" "<<expo_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
            if(theta_tag > -1)
            {   
                pairs+= 1;

                mg1_z2_ = expo_info->expo_data_2[n+expo_info->mg1_idx];
                mg2_z2_ = expo_info->expo_data_2[n+expo_info->mg2_idx];

                mn_z2_ = expo_info->expo_data_2[n+expo_info->mn_idx];
                mu_z2_ = expo_info->expo_data_2[n+expo_info->mu_idx];
                mv_z2_ = expo_info->expo_data_2[n+expo_info->mv_idx];

                // shear estimators rotation (position angle defined as East of North)
                sin_theta = delta_ra/delta_radius;
                cos_theta = delta_dec/delta_radius;

                sin_2theta = 2*sin_theta*cos_theta;
                cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

                sin_4theta = 2*sin_2theta*cos_2theta;
                cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;

                // rotate gal z1
                mg1_z1 = mg1_z1_*cos_2theta - mg2_z1_*sin_2theta;
                mg2_z1 = mg1_z1_*sin_2theta + mg2_z1_*cos_2theta;

                mnu1_z1 = mu_z1_*cos_4theta - mv_z1_*sin_4theta;
                mnu2_z1 = mn_z1_ - mnu1_z1;
                mnu1_z1 = mn_z1_ + mnu1_z1;
                // rotate gal z2
                mg1_z2 = mg1_z2_*cos_2theta - mg2_z2_*sin_2theta;
                mg2_z2 = mg1_z2_*sin_2theta + mg2_z2_*cos_2theta;

                mnu1_z2 = mu_z2_*cos_4theta - mv_z2_*sin_4theta;
                mnu2_z2 = mn_z2_ - mnu1_z2;
                mnu1_z2 = mn_z2_ + mnu1_z2;
                
                // there're zbin_num *zbin_num blocks, iz1 is row, iz2 is the col, each block
                // has a length of mg_bin_num*mg_bin_num*chi_guess_num*theta_bin_num.
                iz2 = expo_info->expo_zbin_label_2[ig2];

                // record the pair separation for the mean separation in that bin
                // which will be the x postion in the last figure
                theta_accum_tag = (iz1_ + iz2)*theta_bin_num + theta_tag;
        
                expo_info->theta_accum[theta_accum_tag] += delta_radius;
                expo_info->theta_num_accum[theta_accum_tag] += 1;


                ////////////////////// the key part of PDF_SYM //////////////////////////////
                ic_len = theta_tag*ir_chi_block_len + (iz1_ + iz2)*expo_info->iz_chi_block_len;

                gg_1 = expo_info->gg_1[loop_label];
                gg_2 = expo_info->gg_2[loop_label];

                temp_tt[2] = mg1_z1 - gg_1*mnu1_z1;
                temp_tt[3] = mg1_z2 - gg_2*mnu1_z2;
                hist_2d_new(temp_tt[2], temp_tt[3], expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_tt, iy_tt);
                
                
                expo_info->expo_num_count_chit[ic_len + iy_tt*mg_bin_num+ix_tt] += 1;
                
                temp_xx[2] = mg2_z1 - gg_1*mnu2_z1;
                temp_xx[3] = mg2_z2 - gg_2*mnu2_z2;

                hist_2d_new(temp_xx[2], temp_xx[3],  expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_xx, iy_xx);
                expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx] += 1;
                loop_label += 1;

                // if(ic_len == 0)
                // {std::cout<<theta_tag<<" "<<iz1<<" "<<iz2<<" "<<iy_xx<<" "<<ix_xx<<" "<<iy_xx*mg_bin_num+ix_xx<<" "<<
                // mg_bin_num<<" "<<expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx]<<std::endl;}
                // std::cout<<0<<" "<<temp_tt[2]<<" "<<temp_tt[3]<<" "<<ix_tt<<" "<<iy_tt<<" "<<gg_1<<std::endl;
                // std::cout<<0<<" "<<temp_xx[2]<<" "<<temp_xx[3]<<" "<<ix_xx<<" "<<iy_xx<<" "<<gg_2<<std::endl;
                for(ic=1; ic<chi_guess_num; ic++)
                {   
                    ic_len += chi_block_len;

                    gg_1 = expo_info->gg_1[loop_label];
                    gg_2 = expo_info->gg_2[loop_label];
                                    
                    bin_para_tt[0] = ix_tt;
                    bin_para_tt[1] = iy_tt;

                    temp_tt[0] = temp_tt[2];
                    temp_tt[1] = temp_tt[3];

                    temp_tt[2] = mg1_z1 - gg_1*mnu1_z1;
                    temp_tt[3] = mg1_z2 - gg_2*mnu1_z2;
                    // hist_2d_new(temp_tt[2], temp_tt[3], field_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_tt, iy_tt);
                    hist_2d_new(expo_info->mg_bin, mg_bin_num, temp_tt, bin_para_tt, ix_tt, iy_tt);

                    expo_info->expo_num_count_chit[ic_len + iy_tt*mg_bin_num+ix_tt] += 1;
                    
                    bin_para_xx[0] = ix_xx;
                    bin_para_xx[1] = iy_xx;
                    
                    temp_xx[0] = temp_xx[2];
                    temp_xx[1] = temp_xx[3];

                    temp_xx[2] = mg2_z1 - gg_1*mnu2_z1;
                    temp_xx[3] = mg2_z2 - gg_2*mnu2_z2;

                    // hist_2d_new(temp_xx[2], temp_xx[3],  field_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_xx, iy_xx);
                    hist_2d_new(expo_info->mg_bin, mg_bin_num, temp_xx, bin_para_xx, ix_xx, iy_xx);
                    expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx] += 1;
                    loop_label += 1;


                    // std::cout<<ic<<" "<<temp_tt[2]<<" "<<temp_tt[3]<<" "<<ix_tt<<" "<<iy_tt<<" "<<gg_1<<std::endl;
                    // std::cout<<ic<<" "<<temp_xx[2]<<" "<<temp_xx[3]<<" "<<ix_xx<<" "<<iy_xx<<" "<<gg_2<<std::endl;
                    // if(ic_len == 0)
                    // {std::cout<<theta_tag<<" "<<iz1<<" "<<iz2<<" "<<iy_xx<<" "<<ix_xx<<" "<<iy_xx*mg_bin_num+ix_xx<<" "<<
                    // mg_bin_num<<" "<<expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx]<<std::endl;}
                    
                }
                if(loop_label >= gg_len){loop_label = 0;}
                ////////////////////// the key part of PDF_SYM -end  //////////////////////////////
                
            }

        }
    }
    // st2 = clock();
    // std::cout<<pairs<<" pairs "<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
    expo_info->gg_pairs = pairs;
    expo_info->loop_label = loop_label;
}

void find_pairs_stack_expo(data_info *expo_info, int expo_label_0, int expo_label_1)
{
    int i, j, m, n, k;
    int ig1, ig2, ig_st;

    int ir, theta_tag, theta_accum_tag, ic;
    int iz1, iz2, iz1_;

    int gal_num_1, gal_num_2;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;
    MY_FLOAT diff;
    MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
    MY_FLOAT mg1_z1_, mg2_z1_, mn_z1_, mu_z1_, mv_z1_;;

    MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;
    MY_FLOAT mg1_z2_, mg2_z2_, mn_z2_, mu_z2_,mv_z2_;

    int ix_tt, iy_tt, ix_xx, iy_xx;
    MY_FLOAT temp_tt[4], temp_xx[4];
    int bin_para_tt[2], bin_para_xx[2];

    int theta_bin_num;
    int chi_guess_num;
    int expo_chi_block_len, ir_chi_block_len,chi_block_len;
    int ir_len, ic_len;
    MY_FLOAT gg_1, gg_2, gg_len;

    MY_FLOAT delta_ra, delta_dec, delta_radius, delta_radius_check;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    double st1, st2;
    double pairs = 0;

    int loop_label;
    loop_label = expo_info->loop_label;
    // loop_label = 0;

    int mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3;
    mg_bin_num = expo_info->mg_bin_num;
    mg_bin_num1 = expo_info->mg_bin_num1;
    mg_bin_num2 = expo_info->mg_bin_num2;
    mg_bin_num3 = expo_info->mg_bin_num3;

    chi_guess_num = expo_info->chi_guess_num;
    expo_chi_block_len = expo_info->expo_chi_block_len;
    ir_chi_block_len = expo_info->ir_chi_block_len;
    chi_block_len = expo_info->chi_block_len;
    gg_len = expo_info->gg_len;
    theta_bin_num = expo_info->theta_bin_num;

    gal_num_1 = expo_info->expo_gal_num[expo_label_0];
    gal_num_2 = expo_info->expo_gal_num[expo_label_1];
    
    // int *mask = new int[gal_num_1]{};

    st1 = clock();
    for(ig1=0; ig1<gal_num_1; ig1++)
    {   
        
        // loop the grid in the first zbin, zbin_label_0
        m = ig1*expo_info->expo_data_col;

        ra_z1 = expo_info->expo_data_1[m+expo_info->ra_idx];
        dec_z1 = expo_info->expo_data_1[m+expo_info->dec_idx];
        cos_dec_z1 = expo_info->expo_data_1[m+expo_info->cos_dec_idx];

        mg1_z1_ = expo_info->expo_data_1[m+expo_info->mg1_idx];
        mg2_z1_ = expo_info->expo_data_1[m+expo_info->mg2_idx];

        mn_z1_ = expo_info->expo_data_1[m+expo_info->mn_idx];
        mu_z1_ = expo_info->expo_data_1[m+expo_info->mu_idx];
        mv_z1_ = expo_info->expo_data_1[m+expo_info->mv_idx];

        iz1 = expo_info->expo_zbin_label_1[ig1];
        iz1_ = iz1*expo_info->zbin_num;

        if(expo_label_0 == expo_label_1){ig_st = ig1+1;}
        else{ig_st = 0;}

        for(ig2=ig_st; ig2<gal_num_2; ig2++)
        {   
            // use the pairs from different exposures

            n = ig2*expo_info->expo_data_col;

            ra_z2 = expo_info->expo_data_2[n+expo_info->ra_idx];
            dec_z2 = expo_info->expo_data_2[n+expo_info->dec_idx];
            cos_dec_z2 = expo_info->expo_data_2[n+expo_info->cos_dec_idx];

            // the seperation angle (arc minute)
            delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
            delta_dec = dec_z2 - dec_z1;
            delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);

            // separation(ra_z1/60, dec_z1/60, ra_z2/60, dec_z2/60, delta_radius_check);
            // diff = fabs(delta_radius - delta_radius_check/Pi*180*60)/delta_radius;
            // if( diff > 0.05)
            // {   
            //     std::cout<<ra_z1<<" "<<ra_z2<<" "<<dec_z1<<" "<<dec_z2<<std::endl;
            //     std::cout<<delta_radius<<" "<<delta_radius_check/Pi*180*60<<" "<<diff<<std::endl;
            // }

            theta_tag = -1;
            for(ir=0; ir<theta_bin_num; ir++)
            {
                if(delta_radius > expo_info->theta_bin[ir] and delta_radius <= expo_info->theta_bin[ir+1])
                {theta_tag=ir;break;}
            }
            // std::cout<<delta_radius<<" "<<expo_info->theta_bin[theta_tag]<<" "<<expo_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
            if(theta_tag > -1)
            {   
                pairs+= 1;

                mg1_z2_ = expo_info->expo_data_2[n+expo_info->mg1_idx];
                mg2_z2_ = expo_info->expo_data_2[n+expo_info->mg2_idx];

                mn_z2_ = expo_info->expo_data_2[n+expo_info->mn_idx];
                mu_z2_ = expo_info->expo_data_2[n+expo_info->mu_idx];
                mv_z2_ = expo_info->expo_data_2[n+expo_info->mv_idx];

                // shear estimators rotation (position angle defined as East of North)
                sin_theta = delta_ra/delta_radius;
                cos_theta = delta_dec/delta_radius;

                sin_2theta = 2*sin_theta*cos_theta;
                cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

                sin_4theta = 2*sin_2theta*cos_2theta;
                cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;

                // rotate gal z1
                mg1_z1 = mg1_z1_*cos_2theta - mg2_z1_*sin_2theta;
                mg2_z1 = mg1_z1_*sin_2theta + mg2_z1_*cos_2theta;

                mnu1_z1 = mu_z1_*cos_4theta - mv_z1_*sin_4theta;
                mnu2_z1 = mn_z1_ - mnu1_z1;
                mnu1_z1 = mn_z1_ + mnu1_z1;
                // rotate gal z2
                mg1_z2 = mg1_z2_*cos_2theta - mg2_z2_*sin_2theta;
                mg2_z2 = mg1_z2_*sin_2theta + mg2_z2_*cos_2theta;

                mnu1_z2 = mu_z2_*cos_4theta - mv_z2_*sin_4theta;
                mnu2_z2 = mn_z2_ - mnu1_z2;
                mnu1_z2 = mn_z2_ + mnu1_z2;
                
                // there're zbin_num *zbin_num blocks, iz1 is row, iz2 is the col, each block
                // has a length of mg_bin_num*mg_bin_num*chi_guess_num*theta_bin_num.
                iz2 = expo_info->expo_zbin_label_2[ig2];

                // record the pair separation for the mean separation in that bin
                // which will be the x postion in the last figure
                theta_accum_tag = (iz1_ + iz2)*theta_bin_num + theta_tag;
        
                expo_info->theta_accum[theta_accum_tag] += delta_radius;
                expo_info->theta_num_accum[theta_accum_tag] += 1;


                ////////////////////// the key part of PDF_SYM //////////////////////////////
                ic_len = theta_tag*ir_chi_block_len + (iz1_ + iz2)*expo_info->iz_chi_block_len;

                gg_1 = expo_info->gg_1[loop_label];
                gg_2 = expo_info->gg_2[loop_label];

                temp_tt[2] = mg1_z1 - gg_1*mnu1_z1;
                temp_tt[3] = mg1_z2 - gg_2*mnu1_z2;
                hist_2d_new(temp_tt[2], temp_tt[3], expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_tt, iy_tt);
                
                
                expo_info->expo_num_count_chit[ic_len + iy_tt*mg_bin_num+ix_tt] += 1;
                
                temp_xx[2] = mg2_z1 - gg_1*mnu2_z1;
                temp_xx[3] = mg2_z2 - gg_2*mnu2_z2;

                hist_2d_new(temp_xx[2], temp_xx[3],  expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_xx, iy_xx);
                expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx] += 1;
                loop_label += 1;

                // if(ic_len == 0)
                // {std::cout<<theta_tag<<" "<<iz1<<" "<<iz2<<" "<<iy_xx<<" "<<ix_xx<<" "<<iy_xx*mg_bin_num+ix_xx<<" "<<
                // mg_bin_num<<" "<<expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx]<<std::endl;}
                // std::cout<<0<<" "<<temp_tt[2]<<" "<<temp_tt[3]<<" "<<ix_tt<<" "<<iy_tt<<" "<<gg_1<<std::endl;
                // std::cout<<0<<" "<<temp_xx[2]<<" "<<temp_xx[3]<<" "<<ix_xx<<" "<<iy_xx<<" "<<gg_2<<std::endl;
                for(ic=1; ic<chi_guess_num; ic++)
                {   
                    ic_len += chi_block_len;

                    gg_1 = expo_info->gg_1[loop_label];
                    gg_2 = expo_info->gg_2[loop_label];
                                    
                    bin_para_tt[0] = ix_tt;
                    bin_para_tt[1] = iy_tt;

                    temp_tt[0] = temp_tt[2];
                    temp_tt[1] = temp_tt[3];

                    temp_tt[2] = mg1_z1 - gg_1*mnu1_z1;
                    temp_tt[3] = mg1_z2 - gg_2*mnu1_z2;
                    // hist_2d_new(temp_tt[2], temp_tt[3], field_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_tt, iy_tt);
                    hist_2d_new(expo_info->mg_bin, mg_bin_num, temp_tt, bin_para_tt, ix_tt, iy_tt);

                    expo_info->expo_num_count_chit[ic_len + iy_tt*mg_bin_num+ix_tt] += 1;
                    
                    bin_para_xx[0] = ix_xx;
                    bin_para_xx[1] = iy_xx;
                    
                    temp_xx[0] = temp_xx[2];
                    temp_xx[1] = temp_xx[3];

                    temp_xx[2] = mg2_z1 - gg_1*mnu2_z1;
                    temp_xx[3] = mg2_z2 - gg_2*mnu2_z2;

                    // hist_2d_new(temp_xx[2], temp_xx[3],  field_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_xx, iy_xx);
                    hist_2d_new(expo_info->mg_bin, mg_bin_num, temp_xx, bin_para_xx, ix_xx, iy_xx);
                    expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx] += 1;
                    loop_label += 1;


                    // std::cout<<ic<<" "<<temp_tt[2]<<" "<<temp_tt[3]<<" "<<ix_tt<<" "<<iy_tt<<" "<<gg_1<<std::endl;
                    // std::cout<<ic<<" "<<temp_xx[2]<<" "<<temp_xx[3]<<" "<<ix_xx<<" "<<iy_xx<<" "<<gg_2<<std::endl;
                    // if(ic_len == 0)
                    // {std::cout<<theta_tag<<" "<<iz1<<" "<<iz2<<" "<<iy_xx<<" "<<ix_xx<<" "<<iy_xx*mg_bin_num+ix_xx<<" "<<
                    // mg_bin_num<<" "<<expo_info->expo_num_count_chix[ic_len + iy_xx*mg_bin_num+ix_xx]<<std::endl;}
                    
                }
                if(loop_label >= gg_len){loop_label = 0;}
                ////////////////////// the key part of PDF_SYM -end  //////////////////////////////
                
            }

        }
    }
    // st2 = clock();
    // std::cout<<pairs<<" pairs "<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
    expo_info->gg_pairs = pairs;
    expo_info->loop_label = loop_label;
}

void initialize_expo_chi_block(data_info *expo_info)
{   
    int i;
    for(i=0; i<expo_info->expo_chi_block_len; i++)
    {   
        expo_info->expo_num_count_chit[i] = 0;
        expo_info->expo_num_count_chix[i] = 0;
    }
    for(i=0; i<expo_info->theta_accum_len; i++)
    {   
        expo_info->theta_accum[i] = 0;
        expo_info->theta_num_accum[i] = 0;
    }
}

void collect_chi_block(data_info *expo_info, int expo_label)
{
    ;
}

void save_expo_data_new(data_info *expo_info, int rank, int task_end_tag)
{
    int row, col;
    int st, i,j;
    char result_file_path[600], set_name[50];

    sprintf(result_file_path, "%s/core_%d_num_count.hdf5", expo_info->result_path, rank);
    
    if(task_end_tag == 1)
    {   
        // if the progress comes to the end, force it to save the buffer (may not be full)
        if(expo_info->block_count > 0)
        {
            // for that case in which the tasks have been completed,
            // but the buffer has not come to the end
            row = expo_info->block_count;
            col = expo_info->block_size_in_buffer;
            sprintf(set_name, "/%d/data", expo_info->total_buffer_num);
            if(expo_info->total_buffer_num == 0)
            {
                write_h5(result_file_path, set_name, expo_info->men_buffer, row, col, true);
            }
            else
            {
                write_h5(result_file_path, set_name, expo_info->men_buffer, row, col, false);
            }
            
            col = 2*expo_info->block_count;
            sprintf(set_name, "/%d/jack_label", expo_info->total_buffer_num);
            write_h5(result_file_path, set_name, expo_info->task_expo_pair_jack_label, 1, col, false);
            
            expo_info->total_buffer_num += 1;
            sprintf(set_name, "/buffer_num");
            write_h5(result_file_path, set_name, &expo_info->total_buffer_num, 1, 1, false);

        }
    }
    else
    {   
        //if it is full, write into the disk, save it first
        if(expo_info->block_count == expo_info->max_block_in_buffer)
        {   
            row = expo_info->max_block_in_buffer;
            col = expo_info->block_size_in_buffer;

            // save the main data, theta, theta_num, chi_tt, chi_xx
            sprintf(set_name, "/%d/data", expo_info->total_buffer_num);
            if(expo_info->total_buffer_num == 0)
            {
                write_h5(result_file_path, set_name, expo_info->men_buffer, row, col, true);
            }
            else
            {
                write_h5(result_file_path, set_name, expo_info->men_buffer, row, col, false);
            }
            // save the jack id of each block
            col = 2*expo_info->max_block_in_buffer;
            sprintf(set_name, "/%d/jack_label", expo_info->total_buffer_num);
            write_h5(result_file_path, set_name, expo_info->task_expo_pair_jack_label, 1, col, false);

            expo_info->total_buffer_num ++;
            expo_info->block_count = 0;
        }

        // assign the buffer of each expo pair to the buffer for the stacked data in the emory
        // the buffer contains data from z_{ij}, however, z_{ij} and z_{ji} should contribute to z_{ij}
        // that is the meaning of merge_data().
        merge_data(expo_info);

        // theta, theta_num
        st = expo_info->block_count*expo_info->block_size_in_buffer;
        for(i=0; i<expo_info->theta_accum_len_true; i++)
        {
            expo_info->men_buffer[st + i] = expo_info->corr_cal_stack_expo_theta_accum[i];
            expo_info->men_buffer[st + expo_info->theta_accum_len_true + i] = expo_info->corr_cal_stack_expo_theta_num_accum[i];
        }

        // chi_tt, chi_xx
        st = expo_info->block_count*expo_info->block_size_in_buffer + 2*expo_info->theta_accum_len_true;
        for(i=0; i<expo_info->expo_chi_block_len_true; i++)
        {
            expo_info->men_buffer[st + i] = expo_info->corr_cal_stack_num_count_chit[i];
            expo_info->men_buffer[st + expo_info->expo_chi_block_len_true + i] = expo_info->corr_cal_stack_num_count_chix[i];
        }

        expo_info->task_expo_pair_jack_label[expo_info->block_count*2] = expo_info->jack_label_1;
        expo_info->task_expo_pair_jack_label[expo_info->block_count*2 + 1] = expo_info->jack_label_2;

        expo_info->block_count ++;
        
    }    
}


void save_expo_data(data_info *expo_info, int expo_label_1, int expo_label_2, int rank)
{   
    int row, col, i;
    char result_file_path[600], set_name[50];

    expo_info->expo_pair_label_1.push_back(expo_label_1);
    expo_info->expo_pair_label_2.push_back(expo_label_2);

    col = expo_info->mg_bin_num;
    row = expo_info->expo_chi_block_len_true/col;

    sprintf(result_file_path, "%s/core_%d_num_count.hdf5", expo_info->result_path, rank);
    sprintf(set_name, "/%d-%d/tt",expo_label_1, expo_label_2);
    
    merge_data(expo_info);

    if(expo_info->result_file_tag == 0)
    {
        write_h5(result_file_path, set_name, expo_info->corr_cal_stack_num_count_chit, row, col, true);
        expo_info->result_file_tag=1;
    }
    else{write_h5(result_file_path, set_name, expo_info->corr_cal_stack_num_count_chit, row, col, false);}

    sprintf(set_name, "/%d-%d/xx",expo_label_1, expo_label_2);
    write_h5(result_file_path, set_name, expo_info->corr_cal_stack_num_count_chix, row, col, false);


    col = expo_info->theta_bin_num;
    row = expo_info->theta_accum_len_true/col;

    sprintf(set_name, "/%d-%d/theta",expo_label_1, expo_label_2);
    write_h5(result_file_path, set_name, expo_info->corr_cal_stack_expo_theta_accum, row, col, false);
    sprintf(set_name, "/%d-%d/theta_num",expo_label_1, expo_label_2);
    write_h5(result_file_path, set_name, expo_info->corr_cal_stack_expo_theta_num_accum, row, col, false);
}

void save_expo_pair_label(data_info *expo_info, int rank)
{
    int i, col, row;
    char result_file_path[600], set_name[50];

    sprintf(result_file_path, "%s/core_%d_num_count.hdf5", expo_info->result_path, rank);

    col = expo_info->expo_pair_label_1.size();
    
    int *labels = new int[col];

    for(i=0; i<col; i++){labels[i] = expo_info->expo_pair_label_1[i];}
    sprintf(set_name, "/pair_1");
    write_h5(result_file_path, set_name, labels, 1, col, false);
    
    for(i=0; i<col; i++){labels[i] = expo_info->expo_pair_label_2[i];}
    sprintf(set_name, "/pair_2");
    write_h5(result_file_path, set_name, labels, 1, col, false);
    
    delete[] labels;
    
}

void save_expo_data(data_info *expo_info, int expo_label, char *file_name)
{   
    int row, col;
    char set_name[50];

    col = expo_info->mg_bin_num;
    row = expo_info->expo_chi_block_len/col;

    sprintf(set_name, "/t");
    write_h5(file_name, set_name, expo_info->expo_num_count_chit, row, col, true);
    sprintf(set_name, "/x");
    write_h5(file_name, set_name, expo_info->expo_num_count_chix, row, col, false);
}


void merge_data(data_info *expo_info)
{   
    int i, j, k, m,n;
    int st, st_ij, st_ji, tag;

    for(j=0; j<expo_info->zbin_num; j++)
    {   
        //////////////////////////////////  theta  //////////////////////////////////////
        // z[i,i] part
        st = (j*expo_info->zbin_num + j)*expo_info->theta_bin_num;
        tag = (j*expo_info->zbin_num - j*(j-1)/2)*expo_info->theta_bin_num;
        // if(i == 0 or i == 10)
        // {std::cout<<j<<" "<<tag<<" "<<st<<std::endl;}

        for(m=0; m<expo_info->theta_bin_num; m++)
        {
            expo_info->corr_cal_stack_expo_theta_accum[tag+m] = expo_info->theta_accum[st+m];
            expo_info->corr_cal_stack_expo_theta_num_accum[tag+m] = expo_info->theta_num_accum[st+m];
        }

        // z[i,j] part, i != j, z[j,i] will be added to z[i,j], for j>i
        for(k=j+1; k<expo_info->zbin_num; k++)
        {   
            tag = (j*expo_info->zbin_num + k - (j*j+j)/2)*expo_info->theta_bin_num;
            // if(i == 0 or i == 10)
            // {std::cout<<j<<" "<<k<<" "<<tag<<" "<<std::endl;}

            st_ij = (j*expo_info->zbin_num + k)*expo_info->theta_bin_num;
            st_ji = (k*expo_info->zbin_num + j)*expo_info->theta_bin_num;

            for(m=0; m<expo_info->theta_bin_num; m++)
            {
                expo_info->corr_cal_stack_expo_theta_accum[tag+m] = expo_info->theta_accum[st_ij+m]+expo_info->theta_accum[st_ji+m];
                expo_info->corr_cal_stack_expo_theta_num_accum[tag+m] = expo_info->theta_num_accum[st_ij+m]+expo_info->theta_num_accum[st_ji+m];
            }
        }
        //////////////////////////////////////////////////////////////////////////////////


        ///////////////////////////////////  PDF number count  ////////////////////////////////////////
        // z[i,i] part
        st = (j*expo_info->zbin_num + j)*expo_info->iz_chi_block_len;
        tag = (j*expo_info->zbin_num - j*(j-1)/2)*expo_info->iz_chi_block_len;
        // if(i == 0 or i == 10)
        // {std::cout<<j<<" "<<tag<<" "<<st<<std::endl;}

        for(m=0; m<expo_info->iz_chi_block_len; m++)
        {
            expo_info->corr_cal_stack_num_count_chit[tag+m] = expo_info->expo_num_count_chit[st+m];
            expo_info->corr_cal_stack_num_count_chix[tag+m] = expo_info->expo_num_count_chix[st+m];
        }

        // z[i,j] part, i != j, z[j,i] will be added to z[i,j], for j>i
        for(k=j+1; k<expo_info->zbin_num; k++)
        {   
            tag = (j*expo_info->zbin_num + k - (j*j+j)/2)*expo_info->iz_chi_block_len;
            // if(i == 0 or i == 10)
            // {std::cout<<j<<" "<<k<<" "<<tag<<" "<<std::endl;}

            st_ij = (j*expo_info->zbin_num + k)*expo_info->iz_chi_block_len;
            st_ji = (k*expo_info->zbin_num + j)*expo_info->iz_chi_block_len;

            for(m=0; m<expo_info->iz_chi_block_len; m++)
            {
                expo_info->corr_cal_stack_num_count_chit[tag+m] = expo_info->expo_num_count_chit[st_ij+m]+expo_info->expo_num_count_chit[st_ji+m];
                expo_info->corr_cal_stack_num_count_chix[tag+m] = expo_info->expo_num_count_chix[st_ij+m]+expo_info->expo_num_count_chix[st_ji+m];
            }
        }
        //////////////////////////////////////////////////////////////////////////////////
    }
}


void hist_2d_fast(MY_FLOAT x, MY_FLOAT y, MY_FLOAT*bins, int bin_num, int bin_num2, int &ix, int &iy)
{
    int i;
    if(x < 0)
    {
        for(i=0; i<bin_num2; i++)
        {
            if(x >= bins[i] and x < bins[i+1]){ix=i;break;}
        }
    }
    else
    {
        for(i=bin_num2; i<bin_num; i++)
        {
            if(x >= bins[i] and x < bins[i+1]){ix=i;break;}
        }
    }
    if(y < 0)
    {
        for(i=0; i<bin_num2; i++)
        {
            if(y >= bins[i] and y < bins[i+1]){iy=i;break;}
        }
    }
    else
    {
        for(i=bin_num2; i<bin_num; i++)
        {
            if(y >= bins[i] and y < bins[i+1]){iy=i;break;}
        }
    }
}

void hist_2d_new(MY_FLOAT x, MY_FLOAT y, MY_FLOAT*bins, int bin_num, int bin_num1,int bin_num2, int bin_num3,int &ix, int &iy)
{   
    int im, im1, im2;
    if(x < 0)
    {
        if(x >= bins[bin_num1])
        { im1 = bin_num1; im2=bin_num2;}
        else
        { im1 = 0; im2=bin_num1;}
    }
    else
    {
        if(x < bins[bin_num3])
        { im1 = bin_num2; im2=bin_num3;}
        else
        { im1 = bin_num3; im2=bin_num;}
    }
    for(im=im1; im<im2; im++)
    {if(x >= bins[im] and x < bins[im+1]){ix=im;break;}}
    
    if(y < 0)
    {
        if(y >= bins[bin_num1])
        { im1 = bin_num1; im2=bin_num2;}
        else
        { im1 = 0; im2=bin_num1;}
    }
    else
    {
        if(y < bins[bin_num3])
        { im1 = bin_num2; im2=bin_num3;}
        else
        { im1 = bin_num3; im2=bin_num;}
    }
    for(im=im1; im<im2; im++)
    {if(y >= bins[im] and y < bins[im+1]){iy=im;break;}}
}

void hist_2d_new(MY_FLOAT x, MY_FLOAT y, MY_FLOAT*bins, int *bin_num,int &ix, int &iy)
{   
    int im, im1, im2;
    if(x < 0)
    {
        if(x >= bins[bin_num[1]])
        { im1 = bin_num[1]; im2=bin_num[2];}
        else
        { im1 = 0; im2=bin_num[1];}
    }
    else
    {
        if(x < bins[bin_num[3]])
        { im1 = bin_num[2]; im2=bin_num[3];}
        else
        { im1 = bin_num[3]; im2=bin_num[0];}
    }
    for(im=im1; im<im2; im++)
    {if(x >= bins[im] and x < bins[im+1]){ix=im;break;}}
    
    if(y < 0)
    {
        if(y >= bins[bin_num[1]])
        { im1 = bin_num[1]; im2=bin_num[2];}
        else
        { im1 = 0; im2=bin_num[1];}
    }
    else
    {
        if(y < bins[bin_num[3]])
        { im1 = bin_num[2]; im2=bin_num[3];}
        else
        { im1 = bin_num[3]; im2=bin_num[0];}
    }
    for(im=im1; im<im2; im++)
    {if(y >= bins[im] and y < bins[im+1]){iy=im;break;}}
}

void hist_2d_new(MY_FLOAT*bins, int bin_num, MY_FLOAT *xy, int *bin_para, int &ix, int &iy)
{
    int im, im1, im2;
    if(xy[2]>xy[0])
    {
        for(im=bin_para[0]; im<bin_num; im++)
        {/*std::cout<<im<<" "<<std::endl;*/if(xy[2] >= bins[im] and xy[2] < bins[im+1]){ix=im;break;}}
    }
    else
    {
        for(im=bin_para[0]; im > -1; im--)
        {if(xy[2] >= bins[im] and xy[2] < bins[im+1]){ix=im;break;}}
    }
    if(xy[3]>xy[1])
    {
        for(im=bin_para[1]; im<bin_num; im++)
        {if(xy[3] >= bins[im] and xy[3] < bins[im+1]){iy=im;break;}}
    }
    else
    {
        for(im=bin_para[1]; im > -1; im--)
        {if(xy[3] >= bins[im] and xy[3] < bins[im+1]){iy=im;break;}}
    }
    // std::cout<<"found "<<ix<<" "<<iy<<std::endl;
}

void hist2d_fast_dev(MY_FLOAT x,  MY_FLOAT y, MY_FLOAT *bins, int bin_num, int pre_xbin_tag, int &xbin_tag, int pre_ybin_tag, int &ybin_tag)
{   
    // because the correction term is very small (g*NU or \xi*NU), 
    // therefore, the corrected (by g_i or \xi_i) G doesn't move too far away from the previous one (g_{i-1} or \xi_{i-1})
    int i, j;
    // int tag;
    if(x >= bins[pre_xbin_tag])
    {
        for(i=pre_xbin_tag; i<bin_num; i++)
        {
            if(x >= bins[i] and x < bins[i+1]){xbin_tag = i;break;}
            // std::cout<<x<<" "<<bins[i]<<" "<<bins[i+1]<<" "<<pre_xbin_tag<<" "<<i<<std::endl;
        }
    }
    else
    {
        for(i=pre_xbin_tag; i>0; i--)
        {   
            j = i-1;
            if(x >= bins[j] and x < bins[i]){xbin_tag = j;break;}
        }
    }
    // xbin_tag = tag;

    if(y >= bins[pre_ybin_tag])
    {
        for(i=pre_ybin_tag; i<bin_num; i++)
        {
            if(y >= bins[i] and y < bins[i+1]){ybin_tag = i;break;}
        }
    }
    else
    {
        for(i=pre_ybin_tag; i>0; i--)
        {
            j = i - 1;
            if(y >= bins[j] and y < bins[i]){ybin_tag = j;break;}
        }
    }
    // ybin_tag = tag;
}


///////////////////////////////////// for the last step, \Xi^2 calculation and estimation of correlation function //////////////////////////////// 
void read_para(corr_cal *all_paras)
{   
    int i, j, m, n;
    char set_name[60], data_path[600];

    sprintf(all_paras->log_path, "%s/log/cal_j%d_%d.dat", all_paras->parent_path, all_paras->resample_num, all_paras->corr_cal_rank);

    sprintf(data_path,"%s/source_list.dat", all_paras->cata_path);

    line_count(data_path, all_paras->corr_cal_expo_num);
    // sprintf(data_path,"Rank %d read %d lines\n", all_paras->corr_cal_rank,all_paras->corr_cal_expo_num);
    // std::cout<<data_path;

    sprintf(data_path,"%s/gg_cor.hdf5", all_paras->cata_path);
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


    all_paras->expo_block_len_in_buffer = 2*(all_paras->expo_chi_block_len_true+all_paras->theta_accum_len_true);

    // for the PDF calculation
    all_paras->corr_cal_chi_num = all_paras->chi_guess_num*all_paras->theta_bin_num*((all_paras->zbin_num*all_paras->zbin_num + all_paras->zbin_num)/2);
    all_paras->corr_cal_final_data_num = all_paras->theta_bin_num*((all_paras->zbin_num*all_paras->zbin_num + all_paras->zbin_num)/2);

    pre_jackknife(all_paras);
    
    if(all_paras->corr_cal_rank == 0)
    {   
        std::cout<<all_paras->corr_cal_expo_num<<" exposures"<<std::endl;
        
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

        std::cout<<"chi^2 num: "<<all_paras->corr_cal_chi_num<<std::endl;
        std::cout<<"expo data len: "<<all_paras->expo_chi_block_len_true<<std::endl;
        std::cout<<"theta num: "<<all_paras->corr_cal_final_data_num<<std::endl;
        std::cout<<"block len: "<<all_paras->expo_block_len_in_buffer<<std::endl;
        std::cout<<"final data num: "<<all_paras->corr_cal_final_data_num<<std::endl;
        
    }

}

void prepare_data(corr_cal *all_paras, int tag)
{   
    int i, j, k, m ,n;
    int buffer_num, block_num;
    char data_path[600];
    char set_name[60];
    int *temp[2];

    all_paras->buffer_num_in_file = new int[all_paras->corr_cal_result_file_num]{};

    if(tag == 0)
    {     
        //std::cout<<all_paras->corr_cal_result_file_num<<" files"<<std::endl;
        all_paras->corr_cal_total_pair_num = 0;
        // no 0'th file, because the CPU 0 is the master for the task distribution
        for(i=1; i<all_paras->corr_cal_result_file_num; i++)
        {
            sprintf(data_path, "%s/core_%d_num_count.hdf5", all_paras->result_path,i);
            sprintf(set_name, "/buffer_num");
            read_h5(data_path, set_name, &buffer_num);
            all_paras->buffer_num_in_file[i] = buffer_num;
            for(j=0;j<buffer_num;j++)
            {
                sprintf(set_name, "/%d/jack_label", j);
                read_h5_datasize(data_path, set_name, block_num);
                all_paras->corr_cal_total_pair_num += block_num/2;
            }
        }
        
        std::cout<<all_paras->corr_cal_total_pair_num<<" Pairs\n";

        sprintf(data_path, "%s/pair_index.hdf5", all_paras->result_path);
        sprintf(set_name, "/buffer_num");
        write_h5(data_path, set_name,  all_paras->buffer_num_in_file, 1, all_paras->corr_cal_result_file_num,true);

        sprintf(all_paras->inform, "Write jack_id list");
        write_log(all_paras->log_path, all_paras->inform);

    }
    else
    {
        sprintf(data_path, "%s/pair_index.hdf5", all_paras->result_path);

        sprintf(set_name, "/buffer_num");
        read_h5(data_path, set_name, all_paras->buffer_num_in_file);

        sprintf(all_paras->inform, "Read jack_id list");
        write_log(all_paras->log_path, all_paras->inform);
    }
    
}

void pre_jackknife(corr_cal *all_paras)
{      
    int i, m, n;
    
    // distribute the resample task to each thread
    // there're "resample_num+1" tasks for "corr_cal_thread_num" CPUs
    all_paras->jackknife_resample_st = new int[all_paras->corr_cal_thread_num];
    all_paras->jackknife_resample_ed = new int[all_paras->corr_cal_thread_num];
    corr_task_alloc(all_paras->resample_num+1, all_paras->corr_cal_thread_num, all_paras->jackknife_resample_st, all_paras->jackknife_resample_ed);
    
    // if(all_paras->corr_cal_rank == 0)
    // {
    //     show_arr(all_paras->jackknife_subsample_pair_st, 1, all_paras->resample_num+1);
    //     show_arr(all_paras->jackknife_subsample_pair_ed, 1, all_paras->resample_num+1);
    //     show_arr(all_paras->jackknife_resample_st, 1, all_paras->corr_cal_thread_num);
    //     show_arr(all_paras->jackknife_resample_ed, 1, all_paras->corr_cal_thread_num);
    // }    

    all_paras->my_jack_st = all_paras->jackknife_resample_st[all_paras->corr_cal_rank];
    all_paras->my_jack_ed = all_paras->jackknife_resample_ed[all_paras->corr_cal_rank];
    if(all_paras->corr_cal_rank == 0)
    {   
        std::cout<<"Resample label st&ed: "<<std::endl;
        show_arr(all_paras->jackknife_resample_st,1,all_paras->corr_cal_thread_num);
        show_arr(all_paras->jackknife_resample_ed,1,all_paras->corr_cal_thread_num);
    }

    for(i=all_paras->my_jack_st; i<all_paras->my_jack_ed; i++)
    {
        all_paras->corr_cal_chi_tt[i] = new double[all_paras->corr_cal_chi_num]{};
        all_paras->corr_cal_chi_xx[i] = new double[all_paras->corr_cal_chi_num]{};

        all_paras->corr_cal_gtt[i]  = new double[all_paras->corr_cal_final_data_num]{};
        all_paras->corr_cal_gxx[i]  = new double[all_paras->corr_cal_final_data_num]{};
        all_paras->corr_cal_gtt_sig[i]  = new double[all_paras->corr_cal_final_data_num]{};
        all_paras->corr_cal_gxx_sig[i]  = new double[all_paras->corr_cal_final_data_num]{};

        all_paras->corr_cal_mean_theta[i] = new double[all_paras->theta_accum_len_true]{};

        // stack the exposure data for jackknife
        all_paras->corr_cal_stack_num_count_chit[i] = new double[all_paras->expo_chi_block_len_true]{};
        all_paras->corr_cal_stack_num_count_chix[i] = new double[all_paras->expo_chi_block_len_true]{};

        all_paras->corr_cal_stack_expo_theta_accum[i] = new double[all_paras->theta_accum_len_true]{};
        all_paras->corr_cal_stack_expo_theta_num_accum[i] = new double[all_paras->theta_accum_len_true]{};
    }

}

void resample_jackknife(corr_cal *all_paras)
{
    // read the data
    // if the expo pair label is the sub-sample label interval,
    // it will be skipped. Throw a sub-sample away each time.
    // the rest will be stack for the calculation.
    // sprintf(all_paras->inform, "Jackknife %d resample-start",resample_label);
    // write_log(all_paras->log_path, all_paras->inform);

    int i,j, k, kk, q;
    int m, n;
    int st, ed;
    int buffer_tag, file_tag, buffer_num;
    int block_num, block_need, block_throw;
    int abort_st, abort_ed;
    char set_name[60], data_path[600];
    double *temp_read;
    int *jack_labels;

    char time_now[40];
    
            
    get_time(time_now, 40);
    sprintf(all_paras->inform, "%s. Jackkinfe resample %d ~ %d starts.",time_now, all_paras->my_jack_st, all_paras->my_jack_ed);
    write_log(all_paras->log_path, all_paras->inform);

    for(i=all_paras->my_jack_st; i<all_paras->my_jack_ed; i++)
    {
        initialize_arr(all_paras->corr_cal_stack_num_count_chit[i], all_paras->expo_chi_block_len_true, 0);
        initialize_arr(all_paras->corr_cal_stack_num_count_chix[i], all_paras->expo_chi_block_len_true, 0);
        initialize_arr(all_paras->corr_cal_stack_expo_theta_accum[i], all_paras->theta_accum_len_true, 0);
        initialize_arr(all_paras->corr_cal_stack_expo_theta_num_accum[i], all_paras->theta_accum_len_true, 0);
    }
    

    for(file_tag=1; file_tag<all_paras->corr_cal_result_file_num; file_tag++)
    {      
        buffer_num = all_paras->buffer_num_in_file[file_tag];
        sprintf(data_path, "%s/core_%d_num_count.hdf5", all_paras->result_path, file_tag);
        // std::cout<<buffer_tag<<std::endl;
        
        get_time(time_now, 40);
        sprintf(all_paras->inform, "%s. read file %d",time_now, file_tag);
        write_log(all_paras->log_path, all_paras->inform);

        for(buffer_tag=0; buffer_tag<buffer_num; buffer_tag++)
        {
            sprintf(set_name, "/%d/jack_label", buffer_tag);
            read_h5_datasize(data_path, set_name, block_num);
            
            // scan the jack_label first, 
            // if there're block needed, read the file
            jack_labels = new int[block_num];
            read_h5(data_path, set_name, jack_labels);
            block_num = block_num/2;
            block_need = 0;
            for(j=0; j<block_num; j++)
            {
                m = jack_labels[2*j];
                n = jack_labels[2*j+1];

                for(k=all_paras->my_jack_st; k<all_paras->my_jack_ed; k++)
                {
                    if(m != k and n != k){block_need++;}
                }
            }

            get_time(time_now, 40);
            sprintf(all_paras->inform, "%s. file %d, group %d, %d blocks needed",time_now, file_tag, buffer_tag, block_need);
            write_log(all_paras->log_path, all_paras->inform);

            if(block_need > 0)
            {
                sprintf(set_name, "/%d/data", buffer_tag);
                read_h5_datasize(data_path,set_name, k);
                temp_read = new double[k];
                read_h5(data_path, set_name, temp_read);
                block_num = k/all_paras->expo_block_len_in_buffer;

                // stack
                for(j=0; j<block_num; j++)
                {
                    m = jack_labels[2*j];
                    n = jack_labels[2*j+1];

                    for(k=all_paras->my_jack_st; k<all_paras->my_jack_ed; k++)
                    {
                        if(k != m and k != n)
                        {
                            // if the exposure pair labels does not contain the jack_id, add the data to the buffer
                            st = j*all_paras->expo_block_len_in_buffer;
                            for(kk=0; kk<all_paras->theta_accum_len_true; kk++)
                            {
                                all_paras->corr_cal_stack_expo_theta_accum[k][kk] += temp_read[st + kk];
                                all_paras->corr_cal_stack_expo_theta_num_accum[k][kk] += temp_read[st + all_paras->theta_accum_len_true + kk];
                            }

                            st = j*all_paras->expo_block_len_in_buffer + 2*all_paras->theta_accum_len_true;
                            for(kk=0; kk<all_paras->expo_chi_block_len_true; kk++)
                            {
                                all_paras->corr_cal_stack_num_count_chit[k][kk] += temp_read[st + kk];
                                all_paras->corr_cal_stack_num_count_chix[k][kk] += temp_read[st + all_paras->expo_chi_block_len_true + kk];
                            }
                        }
                    }
                }
                delete[] temp_read;

                get_time(time_now, 40);
                sprintf(all_paras->inform, "%s. file %d, group %d -- end ",time_now, file_tag, buffer_tag);
                write_log(all_paras->log_path, all_paras->inform);

                if(all_paras->corr_cal_rank == 0){std::cout<<all_paras->inform<<std::endl;}
            }
            delete[] jack_labels;

        }

        // if(all_paras->corr_cal_rank == all_paras->corr_cal_thread_num - 1)
        // {
        //     if(file_tag == 5 or file_tag == 10 or file_tag == 20)
        //     {   
        //         std::cout<<"All stacked: "<<std::endl;
        //         for(i=0; i<all_paras->theta_accum_len_true; i++)
        //         {
        //             std::cout<<all_paras->corr_cal_stack_expo_theta_accum[200][i]/all_paras->corr_cal_stack_expo_theta_num_accum[200][i]<<" ";
        //         }
        //         std::cout<<std::endl;
        //     }
        // }

        get_time(time_now, 40);
        sprintf(all_paras->inform, "%s. file %d -- end",time_now, file_tag);
        write_log(all_paras->log_path, all_paras->inform);
    }
    get_time(time_now, 40);
    sprintf(all_paras->inform, "%s. Jackkinfe resample %d ~ %d end.",time_now, all_paras->my_jack_st, all_paras->my_jack_ed);
    write_log(all_paras->log_path, all_paras->inform);
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


void corr_calculate(corr_cal *all_paras)
{   
    int i, j, k, tag;
    int jacK_tag;
    int mg_bin_num = all_paras->mg_bin_num;
 
    double chisq_tt, chisq_xx;
 
    double *temp_tt = new double[mg_bin_num*mg_bin_num]{};
    double *temp_xx = new double[mg_bin_num*mg_bin_num]{};

    double *chi_gtt_fit = new double[all_paras->chi_guess_num];
    double *chi_gxx_fit = new double[all_paras->chi_guess_num];
    double gh, gh_sig, chi_min_fit;

    double *chisq_fit_coeff = new double[3];

    // calculate chi squared
    for(jacK_tag=all_paras->my_jack_st; jacK_tag<all_paras->my_jack_ed; jacK_tag++)
    {
        for(i=0;i<all_paras->corr_cal_chi_num;i++)
        {
            for(j=0;j<mg_bin_num*mg_bin_num;j++)
            {
                tag = i*mg_bin_num*mg_bin_num + j;
                temp_tt[j] = all_paras->corr_cal_stack_num_count_chit[jacK_tag][tag];
                temp_xx[j] = all_paras->corr_cal_stack_num_count_chix[jacK_tag][tag];
            }
            
            chisq_2d(temp_tt,mg_bin_num, chisq_tt);
            chisq_2d(temp_xx,mg_bin_num, chisq_xx);

            all_paras->corr_cal_chi_tt[jacK_tag][i] = chisq_tt;
            all_paras->corr_cal_chi_xx[jacK_tag][i] = chisq_xx;
        }

        // fitting
        for(i=0; i<all_paras->corr_cal_final_data_num;i++)
        {   
            for(j=0;j<all_paras->chi_guess_num;j++)
            {   
                tag = i*all_paras->chi_guess_num + j;
                chi_gtt_fit[j] = all_paras->corr_cal_chi_tt[jacK_tag][tag];
                chi_gxx_fit[j] = all_paras->corr_cal_chi_xx[jacK_tag][tag];
            }
            // if(resample_label == 1)
            // {
            //     show_arr(chi_gtt_fit,1,all_paras->chi_guess_num);
            //     show_arr(chi_gxx_fit,1,all_paras->chi_guess_num);
            // }

            fit_shear(all_paras->corr_cal_chi_guess, chi_gtt_fit, all_paras->chi_guess_num, gh, gh_sig, chi_min_fit, chisq_fit_coeff, 150);
            all_paras->corr_cal_gtt[jacK_tag][i] = gh;
            all_paras->corr_cal_gtt_sig[jacK_tag][i] = gh_sig;

            fit_shear(all_paras->corr_cal_chi_guess, chi_gxx_fit, all_paras->chi_guess_num, gh, gh_sig, chi_min_fit, chisq_fit_coeff, 150);
            all_paras->corr_cal_gxx[jacK_tag][i] = gh;
            all_paras->corr_cal_gxx_sig[jacK_tag][i] = gh_sig;
        }


        // calculate the mean theta
        for(i=0;i<all_paras->theta_accum_len_true;i++)
        {
            all_paras->corr_cal_mean_theta[jacK_tag][i] = all_paras->corr_cal_stack_expo_theta_accum[jacK_tag][i]/all_paras->corr_cal_stack_expo_theta_num_accum[jacK_tag][i];
            // std::cout<<all_paras->corr_cal_stack_expo_theta_accum[i]<<" "<<all_paras->corr_cal_stack_expo_theta_num_accum[i]<<std::endl;
        }
        
        // std::cout<<"chi block num: "<<all_paras->corr_cal_chi_num<<" Final data point num: "<<all_paras->corr_cal_final_data_num<<" Theta point num: "<<all_paras->theta_accum_len_true<<std::endl;


        // sprintf(all_paras->inform, "Jackknife %d end",resample_label);
        // write_log(all_paras->log_path, all_paras->inform);
    }
    delete[] temp_tt;
    delete[] temp_xx;
    delete[] chi_gtt_fit;
    delete[] chi_gxx_fit;
}

void corr_task_alloc(int total_task_num, int portion, int *label_st, int *label_ed)
{   
    int i, j, k, m, n;
    int *temp;

    temp = new int[portion];

    m = total_task_num/portion;
    n = total_task_num%portion;

    for(i=0; i<portion; i++)
    {
        temp[i] = m;
        if(i<n){temp[i] += 1;}
    }

    for(i=0; i<portion; i++)
    {   
        m = 0;
        for(j=0;j<i;j++){m += temp[j];}
        label_st[i] = m;
        label_ed[i] = m + temp[i];
    }

    delete[] temp;
}

void save_result(corr_cal *all_paras)
{
    char data_path[600];
    char set_name[100];

    int i, j;
    int m, n;
    int row, col;
    bool overwrite;

    sprintf(data_path, "%s/result_%d.hdf5", all_paras->result_path, all_paras->resample_num);

    for(i=all_paras->my_jack_st; i<all_paras->my_jack_ed; i++)
    {   
        // rank 0 will start write first
        if(all_paras->corr_cal_rank == 0 and i == 0){overwrite = true;}
        else{overwrite=false;}

        // the \chi squared  
        col = all_paras->chi_guess_num;
        row = all_paras->corr_cal_chi_num/all_paras->chi_guess_num;

        sprintf(set_name, "/%d/chi_tt",i);
        write_h5(data_path, set_name, all_paras->corr_cal_chi_tt[i], row, col, overwrite);

        sprintf(set_name, "/%d/chi_xx",i);
        write_h5(data_path, set_name, all_paras->corr_cal_chi_xx[i], row, col, false);

        // the signal
        col = all_paras->theta_bin_num;
        row = all_paras->corr_cal_final_data_num/all_paras->theta_bin_num;

        sprintf(set_name, "/%d/tt",i);
        write_h5(data_path, set_name, all_paras->corr_cal_gtt[i], row, col, false);
        sprintf(set_name, "/%d/tt_sig",i);
        write_h5(data_path, set_name, all_paras->corr_cal_gtt_sig[i], row, col, false);

        sprintf(set_name, "/%d/xx",i);
        write_h5(data_path, set_name, all_paras->corr_cal_gxx[i], row, col, false);
        sprintf(set_name, "/%d/xx_sig",i);
        write_h5(data_path, set_name, all_paras->corr_cal_gxx_sig[i], row, col, false);

        sprintf(set_name, "/%d/theta",i);
        write_h5(data_path, set_name, all_paras->corr_cal_mean_theta[i], row, col, false);

        sprintf(set_name, "/%d/total_gal_count",i);
        write_h5(data_path, set_name, all_paras->corr_cal_stack_expo_theta_num_accum[i], row, col, false);

        // row = 1;
        // col = all_paras->expo_chi_block_len_true;

        // sprintf(set_name, "/%d/chi_tt_count",i);
        // write_h5(data_path, set_name, all_paras->corr_cal_stack_num_count_chit[i], row, col, false);

        // sprintf(set_name, "/%d/chi_xx_count",i);
        // write_h5(data_path, set_name, all_paras->corr_cal_stack_num_count_chix[i], row, col, false);
    }
    
}
