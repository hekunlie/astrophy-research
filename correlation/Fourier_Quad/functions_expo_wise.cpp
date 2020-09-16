#include"functions_expo_wise.h"

void initialize(data_info *expo_info, int total_expo_num)
{
    int i, j, k;
    char set_name[100];
    char data_path[600];

    expo_info->total_expo_num = total_expo_num;

    for(i=0;i<total_expo_num;i++)
    {   
        // the expo file directories
        expo_info->expo_name_path[i] = new char[400];
        expo_info->expo_name[i] = new char[100];
    }
    
    expo_info->expo_cen_ra = new MY_FLOAT[total_expo_num]{};  
    expo_info->expo_cen_dec = new MY_FLOAT[total_expo_num]{}; 
    expo_info->expo_delta_ra = new MY_FLOAT[total_expo_num]{};  
    expo_info->expo_delta_dec = new MY_FLOAT[total_expo_num]{};
    expo_info->expo_delta_len = new MY_FLOAT[total_expo_num]{};
    expo_info->expo_cen_cos_dec = new MY_FLOAT[total_expo_num]{};
    expo_info->expo_gal_num = new int[total_expo_num]{};

    sprintf(data_path,"%s/cata/source_list.dat", expo_info->parent_path);
    // read the infomation of each expo
    read_list(data_path, expo_info, i);

    // array length, all elements
    sprintf(set_name, "/data");
    read_h5_datasize(expo_info->expo_name_path[0], set_name, j);
    // data col = data length/gal_num, it's the same for all expos
    expo_info->expo_data_col = j / expo_info->expo_gal_num[0];


    ///////////////// read the inform of the PDF_SYM  ////////////////
    sprintf(data_path,"%s/cata/gg_cor.hdf5", expo_info->parent_path);
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


    expo_info->theta_accum = new double[expo_info->theta_accum_len];
    expo_info->theta_num_accum = new double[expo_info->theta_accum_len];

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
    read_file_num = line_count;

}

void read_data(data_info *expo_info)
{
    int i, j;
    char set_name[100];
    for(i=0; i<expo_info->total_expo_num; i++)
    {   
        ///////////////// read all the expos  //////////////////////////
        expo_info->expo_zbin_label[i] = new int[expo_info->expo_gal_num[i]]{};
        sprintf(set_name, "/redshift_label");
        read_h5(expo_info->expo_name_path[i], set_name, expo_info->expo_zbin_label[i]);

        expo_info->expo_data[i] = new MY_FLOAT[expo_info->expo_gal_num[i]*expo_info->expo_data_col]{};
        sprintf(set_name, "/data");
        read_h5(expo_info->expo_name_path[i], set_name, expo_info->expo_data[i]);
       
        // // read the gal num in each zbin stored in each expo file
        // expo_info->expo_num_in_zbin[i] = new int[expo_info->zbin_num];
        // sprintf(set_name, "/num_in_zbin");
        // read_h5(expo_info->expo_name_path[i], set_name, expo_info->expo_num_in_zbin[i]);


        // expo_info->expo_zbin_st[i] = new int[expo_info->zbin_num];
        // sprintf(set_name, "/zbin_st");
        // read_h5(expo_info->expo_name_path[i], set_name, expo_info->expo_zbin_st[i]);

        // expo_info->expo_zbin_ed[i] = new int[expo_info->zbin_num];
        // sprintf(set_name, "/zbin_ed");
        // read_h5(expo_info->expo_name_path[i], set_name, expo_info->expo_zbin_ed[i]);
    }
}

void read_expo_data_1(data_info *expo_info,int expo_label)
{
    char set_name[100];
    if(expo_info->data_read_label_1 == 1)
    {
        delete[] expo_info->expo_data_1;
        delete[] expo_info->expo_zbin_label_1;
        expo_info->data_read_label_1 = 0;
    }

    expo_info->expo_zbin_label_1 = new int[expo_info->expo_gal_num[expo_label]]{};
    sprintf(set_name, "/redshift_label");
    read_h5(expo_info->expo_name_path[expo_label], set_name, expo_info->expo_zbin_label_1);

    expo_info->expo_data_1 = new MY_FLOAT[expo_info->expo_gal_num[expo_label]*expo_info->expo_data_col]{};
    sprintf(set_name, "/data");
    read_h5(expo_info->expo_name_path[expo_label], set_name, expo_info->expo_data_1);

    expo_info->data_read_label_1 = 1;
    
}

void read_expo_data_2(data_info *expo_info,int expo_label)
{
    char set_name[100];
    if(expo_info->data_read_label_2 == 1)
    {
        delete[] expo_info->expo_data_2;
        delete[] expo_info->expo_zbin_label_2;
        expo_info->data_read_label_2 = 0;
    }   

    expo_info->expo_zbin_label_2 = new int[expo_info->expo_gal_num[expo_label]]{};
    sprintf(set_name, "/redshift_label");
    read_h5(expo_info->expo_name_path[expo_label], set_name, expo_info->expo_zbin_label_2);

    expo_info->expo_data_2 = new MY_FLOAT[expo_info->expo_gal_num[expo_label]*expo_info->expo_data_col]{};
    sprintf(set_name, "/data");
    read_h5(expo_info->expo_name_path[expo_label], set_name, expo_info->expo_data_2);

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
        for(j=i+1; j<expo_info->total_expo_num; j++)
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

    label = 0;
    if(delta_radius <= expo_info->theta_bin[expo_info->theta_bin_num])
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

void find_pairs(data_info *expo_info, int expo_label_0, int expo_label_1)
{
    int i, j, m, n, k;
    int ig1, ig2;

    int ir, theta_tag, ic;
    int iz1, iz2;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;

    MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
    MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;
    MY_FLOAT temp_x, temp_y;
    int ix_tt, iy_tt, ix_xx, iy_xx;
    MY_FLOAT temp_tt[4], temp_xx[4];
    int bin_para_tt[2], bin_para_xx[2];

    int ix, iy, im, im1,im2;
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

    st1 = clock();
    for(ig1=0; ig1<expo_info->expo_gal_num[expo_label_0]; ig1++)
    {   
        
        // loop the grid in the first zbin, zbin_label_0
        m = ig1*expo_info->expo_data_col;

        ra_z1 = expo_info->expo_data[expo_label_0][m+expo_info->ra_idx];
        dec_z1 = expo_info->expo_data[expo_label_0][m+expo_info->dec_idx];
        cos_dec_z1 = expo_info->expo_data[expo_label_0][m+expo_info->cos_dec_idx];

        mg1_z1 = expo_info->expo_data[expo_label_0][m+expo_info->mg1_idx];
        mg2_z1 = expo_info->expo_data[expo_label_0][m+expo_info->mg2_idx];

        mnu1_z1 = expo_info->expo_data[expo_label_0][m+expo_info->mn_idx] +
                    expo_info->expo_data[expo_label_0][m+expo_info->mu_idx];
        mnu2_z1 = expo_info->expo_data[expo_label_0][m+expo_info->mn_idx] -
                    expo_info->expo_data[expo_label_0][m+expo_info->mu_idx];
        
        iz1 = expo_info->expo_zbin_label[expo_label_0][ig1]*expo_info->zbin_num;

        for(ig2=0; ig2<expo_info->expo_gal_num[expo_label_1]; ig2++)
        {   
            n = ig2*expo_info->expo_data_col;

            ra_z2 = expo_info->expo_data[expo_label_1][n+expo_info->ra_idx];
            dec_z2 = expo_info->expo_data[expo_label_1][n+expo_info->dec_idx];
            cos_dec_z2 = expo_info->expo_data[expo_label_1][n+expo_info->cos_dec_idx];

            // the seperation angle (arc minute)
            delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
            delta_dec = dec_z2 - dec_z1;
            delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);
            // separation(ra_z1/60, dec_z1/60, ra_z2/60, dec_z2/60, delta_radius_check);
            // std::cout<<delta_radius<<" "<<delta_radius_check/Pi*180*60<<std::endl;

            theta_tag = -1;
            for(ir=0; ir<expo_info->theta_bin_num; ir++)
            {
                if(delta_radius > expo_info->theta_bin[ir] and delta_radius <= expo_info->theta_bin[ir+1])
                {theta_tag=ir;break;}
            }
            // std::cout<<delta_radius<<" "<<expo_info->theta_bin[theta_tag]<<" "<<expo_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
            if(theta_tag > -1)
            {   
                pairs+= 1;

                // shear estimators rotation (position angle defined as East of North)
                sin_theta = delta_ra/delta_radius;
                cos_theta = delta_dec/delta_radius;

                sin_2theta = 2*sin_theta*cos_theta;
                cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

                sin_4theta = 2*sin_2theta*cos_2theta;
                cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;


                mg1_z2 = expo_info->expo_data[expo_label_1][n+expo_info->mg1_idx]*cos_2theta - 
                        expo_info->expo_data[expo_label_1][n+expo_info->mg2_idx]*sin_2theta;
                mg2_z2 = expo_info->expo_data[expo_label_1][n+expo_info->mg1_idx]*sin_2theta + 
                        expo_info->expo_data[expo_label_1][n+expo_info->mg2_idx]*cos_2theta;

                mnu1_z2 = expo_info->expo_data[expo_label_1][n+expo_info->mu_idx]*cos_4theta -
                        expo_info->expo_data[expo_label_1][n+expo_info->mv_idx]*sin_4theta;
                mnu2_z2 = mnu1_z2;

                mnu1_z2 = expo_info->expo_data[expo_label_1][n+expo_info->mn_idx] + mnu2_z2;
                mnu2_z2 = expo_info->expo_data[expo_label_1][n+expo_info->mn_idx] - mnu2_z2;
                
                // there're zbin_num *zbin_num blocks, iz1 is row, iz2 is the col, each block
                // has a length of mg_bin_num*mg_bin_num*chi_guess_num*theta_bin_num.
                iz2 = expo_info->expo_zbin_label[expo_label_1][ig2];

                ////////////////////// the key part of PDF_SYM //////////////////////////////
                ic_len = theta_tag*ir_chi_block_len + (iz1 + iz2)*expo_info->iz_chi_block_len;

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

void find_pairs_new(data_info *expo_info, int expo_label_0, int expo_label_1)
{
    int i, j, m, n, k;
    int ig1, ig2;

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

        for(ig2=0; ig2<gal_num_2; ig2++)
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

void save_expo_data(data_info *expo_info, int expo_label_1, int expo_label_2, int rank)
{   
    int row, col, i;
    char result_path[600], set_name[50];

    expo_info->expo_pair_label_1.push_back(expo_label_1);
    expo_info->expo_pair_label_2.push_back(expo_label_2);

    col = expo_info->mg_bin_num;
    row = expo_info->expo_chi_block_len_true/col;

    sprintf(result_path, "%s/result/core_%d_num_count.hdf5", expo_info->parent_path, rank);
    sprintf(set_name, "/%d-%d/tt",expo_label_1, expo_label_2);
    
    merge_data(expo_info);

    if(expo_info->result_file_tag == 0)
    {
        write_h5(result_path, set_name, expo_info->corr_cal_stack_num_count_chit, row, col, true);
        expo_info->result_file_tag=1;
    }
    else{write_h5(result_path, set_name, expo_info->corr_cal_stack_num_count_chit, row, col, false);}

    sprintf(set_name, "/%d-%d/xx",expo_label_1, expo_label_2);
    write_h5(result_path, set_name, expo_info->corr_cal_stack_num_count_chix, row, col, false);


    col = expo_info->theta_bin_num;
    row = expo_info->theta_accum_len_true/col;

    sprintf(set_name, "/%d-%d/theta",expo_label_1, expo_label_2);
    write_h5(result_path, set_name, expo_info->corr_cal_stack_expo_theta_accum, row, col, false);
    sprintf(set_name, "/%d-%d/theta_num",expo_label_1, expo_label_2);
    write_h5(result_path, set_name, expo_info->corr_cal_stack_expo_theta_num_accum, row, col, false);
}

void save_expo_pair_label(data_info *expo_info, int rank)
{
    int i, col, row;
    char result_path[600], set_name[50];

    sprintf(result_path, "%s/result/core_%d_num_count.hdf5", expo_info->parent_path, rank);

    col = expo_info->expo_pair_label_1.size();
    
    int *labels = new int[col];

    for(i=0; i<col; i++){labels[i] = expo_info->expo_pair_label_1[i];}
    sprintf(set_name, "/pair_1");
    write_h5(result_path, set_name, labels, 1, col, false);
    
    for(i=0; i<col; i++){labels[i] = expo_info->expo_pair_label_2[i];}
    sprintf(set_name, "/pair_2");
    write_h5(result_path, set_name, labels, 1, col, false);
    
    delete[] labels;
    
}

void save_expo_data(data_info *expo_info, int expo_label, char *file_name)
{   
    int row, col;
    char result_path[600], set_name[50];

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


        ///////////////////////////////////  number count  ////////////////////////////////////////
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



///////////////////////////////////// for the last step, \Xi^2 calculation and estimation of correlation function //////////////////////////////// 
void read_para(corr_cal *all_paras)
{   
    int i, j, m, n;
    char set_name[60], data_path[600];

    sprintf(all_paras->log_path, "%s/log/cal_%d.dat", all_paras->parent_path, all_paras->corr_cal_rank);

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

    if(all_paras->corr_cal_rank == 0)
    {
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
        std::cout<<"expo len: "<<all_paras->expo_chi_block_len_true<<std::endl;
        std::cout<<"final data num: "<<all_paras->corr_cal_final_data_num<<std::endl;
    }

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
    
    if(all_paras->corr_cal_rank == 0){std::cout<<"All expo pairs: "<<all_paras->corr_cal_total_pair_num<<std::endl;}
    if(all_paras->corr_cal_rank == 0)
    {
        sprintf(data_path, "%s/pairs.hdf5",all_paras->parent_path);
        sprintf(set_name, "/label_1");
        write_h5(data_path, set_name, all_paras->corr_cal_expo_pair_label[0], 1, all_paras->corr_cal_total_pair_num, true);
        sprintf(set_name, "/label_2");
        write_h5(data_path, set_name, all_paras->corr_cal_expo_pair_label[1], 1, all_paras->corr_cal_total_pair_num, false);

        sprintf(set_name, "/file_tag");
        write_h5(data_path, set_name, all_paras->corr_cal_expo_pair_file_label, 1, all_paras->corr_cal_total_pair_num, false);
    }
    
    // prepare for jackknife 
    pre_jackknife(all_paras);
    
    // if(all_paras->corr_cal_rank == 0)
    // {
    //     show_arr(all_paras->jackknife_subsample_pair_st, 1, all_paras->resample_num+1);
    //     show_arr(all_paras->jackknife_subsample_pair_ed, 1, all_paras->resample_num+1);
    //     show_arr(all_paras->jackknife_resample_st, 1, all_paras->corr_cal_thread_num);
    //     show_arr(all_paras->jackknife_resample_ed, 1, all_paras->corr_cal_thread_num);
    // }
}

void pre_jackknife(corr_cal *all_paras)
{      
    int i, m, n;

    all_paras->jackknife_subsample_pair_st = new int[all_paras->resample_num+1];
    all_paras->jackknife_subsample_pair_ed = new int[all_paras->resample_num+1];
    all_paras->jackknife_resample_st = new int[all_paras->corr_cal_thread_num];
    all_paras->jackknife_resample_ed = new int[all_paras->corr_cal_thread_num];
    
    // the first one is the all pairs
    all_paras->jackknife_subsample_pair_st[0] = -1;
    all_paras->jackknife_subsample_pair_ed[0] = -1;
    
    // decide the start & end label of each sub-sample
    // each time, the expo-pairs of which the label > start and <=end
    // will be abandoned from the calculation
    corr_task_alloc(all_paras->corr_cal_total_pair_num, all_paras->resample_num, &all_paras->jackknife_subsample_pair_st[1],&all_paras->jackknife_subsample_pair_ed[1]);
    
    // distribute the resample task to each thread
    // there're "resample_num+1" tasks for "corr_cal_thread_num" CPUs
    corr_task_alloc(all_paras->resample_num+1, all_paras->corr_cal_thread_num, all_paras->jackknife_resample_st, all_paras->jackknife_resample_ed);
    
    // if(all_paras->corr_cal_rank == 0)
    // {
    //     show_arr(all_paras->jackknife_subsample_pair_st, 1, all_paras->resample_num+1);
    //     show_arr(all_paras->jackknife_subsample_pair_ed, 1, all_paras->resample_num+1);
    //     show_arr(all_paras->jackknife_resample_st, 1, all_paras->corr_cal_thread_num);
    //     show_arr(all_paras->jackknife_resample_ed, 1, all_paras->corr_cal_thread_num);
    // }    

    m = all_paras->jackknife_resample_st[all_paras->corr_cal_rank];
    n = all_paras->jackknife_resample_ed[all_paras->corr_cal_rank];
    
    for(i=m; i<n; i++)
    {
        all_paras->corr_cal_chi_tt[i] = new double[all_paras->corr_cal_chi_num];
        all_paras->corr_cal_chi_xx[i] = new double[all_paras->corr_cal_chi_num];

        all_paras->corr_cal_gtt[i]  = new double[all_paras->corr_cal_final_data_num];
        all_paras->corr_cal_gxx[i]  = new double[all_paras->corr_cal_final_data_num];
        all_paras->corr_cal_gtt_sig[i]  = new double[all_paras->corr_cal_final_data_num];
        all_paras->corr_cal_gxx_sig[i]  = new double[all_paras->corr_cal_final_data_num];

        all_paras->corr_cal_mean_theta[i] = new double[all_paras->expo_chi_block_len_true]{};
    }

}

void resample_jackknife(corr_cal *all_paras,int resample_label)
{
    // read the data
    // if the expo pair label is the sub-sample label interval,
    // it will be skipped. Throw a sub-sample away each time.
    // the rest will be stack for the calculation.
    sprintf(all_paras->inform, "Jackknife %d resample-start",resample_label);
    write_log(all_paras->log_path, all_paras->inform);

    int i,j, k, q;
    int m, n;
    int file_tag;
    int abort_st, abort_ed;
    char set_name[60], data_path[600];

    initialize_arr(all_paras->corr_cal_stack_num_count_chit, all_paras->expo_chi_block_len_true, 0);
    initialize_arr(all_paras->corr_cal_stack_num_count_chix, all_paras->expo_chi_block_len_true, 0);
    initialize_arr(all_paras->corr_cal_stack_expo_theta_accum, all_paras->theta_accum_len_true, 0);
    initialize_arr(all_paras->corr_cal_stack_expo_theta_num_accum, all_paras->theta_accum_len_true, 0);

    abort_st = all_paras->jackknife_subsample_pair_st[resample_label];
    abort_ed = all_paras->jackknife_subsample_pair_ed[resample_label];

    for(i=0; i<all_paras->corr_cal_total_pair_num; i++)
    {   
        
        file_tag = all_paras->corr_cal_expo_pair_file_label[i];

        m = all_paras->corr_cal_expo_pair_label[0][i];
        n = all_paras->corr_cal_expo_pair_label[1][i];
        if(m>= abort_st and m < abort_ed){continue;}
        if(n>= abort_st and n < abort_ed){continue;}

        sprintf(data_path, "%s/result/core_%d_num_count.hdf5", all_paras->parent_path,file_tag);

        
        sprintf(set_name, "/%d-%d/tt", m, n);
        read_h5(data_path, set_name, all_paras->expo_num_count_chit);
        sprintf(set_name, "/%d-%d/xx", m, n);
        read_h5(data_path, set_name, all_paras->expo_num_count_chix);
        sprintf(set_name, "/%d-%d/theta", m, n);
        read_h5(data_path, set_name, all_paras->theta_accum);
        sprintf(set_name, "/%d-%d/theta_num", m, n);
        read_h5(data_path, set_name, all_paras->theta_num_accum);

        // stack 
        for(j=0; j<all_paras->expo_chi_block_len_true; j++)
        {
            all_paras->corr_cal_stack_num_count_chit[j] += all_paras->expo_num_count_chit[j];
            all_paras->corr_cal_stack_num_count_chix[j] += all_paras->expo_num_count_chix[j];

            if(j<all_paras->theta_accum_len_true)
            {
                all_paras->corr_cal_stack_expo_theta_accum[j] += all_paras->theta_accum[j];
                all_paras->corr_cal_stack_expo_theta_num_accum[j] += all_paras->theta_num_accum[j];
            }
        }
        // q = int(all_paras->corr_cal_total_pair_num*0.05);
        // k = i/q;
        // if(k%q == 0)
        // {
        //     sprintf(all_paras->inform, "Jackknife %d read file %d%%",resample_label, k);
        //     write_log(all_paras->log_path, all_paras->inform);
        // }
    }


    sprintf(all_paras->inform, "Jackknife %d resample-end",resample_label);
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

void corr_calculate(corr_cal *all_paras, int resample_label)
{   
    int i, j, k, tag;
 
    int mg_bin_num = all_paras->mg_bin_num;
 
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
        // if(resample_label == 1)
        // {
        //     show_arr(chi_gtt_fit,1,all_paras->chi_guess_num);
        //     show_arr(chi_gxx_fit,1,all_paras->chi_guess_num);
        // }

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
        // std::cout<<all_paras->corr_cal_stack_expo_theta_accum[i]<<" "<<all_paras->corr_cal_stack_expo_theta_num_accum[i]<<std::endl;
    }
    
    // std::cout<<"chi block num: "<<all_paras->corr_cal_chi_num<<" Final data point num: "<<all_paras->corr_cal_final_data_num<<" Theta point num: "<<all_paras->theta_accum_len_true<<std::endl;
    delete[] temp_tt;
    delete[] temp_xx;
    delete[] chi_gtt_fit;
    delete[] chi_gxx_fit;

    sprintf(all_paras->inform, "Jackknife %d calculate",resample_label);
    write_log(all_paras->log_path, all_paras->inform);

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

    sprintf(data_path, "%s/result/result.hdf5", all_paras->parent_path);

    m = all_paras->jackknife_resample_st[all_paras->corr_cal_rank];
    n = all_paras->jackknife_resample_ed[all_paras->corr_cal_rank];

    for(i=m; i<n; i++)
    {   
        if(all_paras->corr_cal_rank == 0 and i == m){overwrite = true;}
        else{overwrite=false;}

        // the \chi squared  
        row = all_paras->chi_guess_num;
        col = all_paras->corr_cal_chi_num/all_paras->chi_guess_num;

        sprintf(set_name, "/%d/chi_tt",i);
        write_h5(data_path, set_name, all_paras->corr_cal_chi_tt[i], row, col, overwrite);

        sprintf(set_name, "/%d/chi_xx",i);
        write_h5(data_path, set_name, all_paras->corr_cal_chi_xx[i], row, col, false);

        // the signal
        row = all_paras->theta_bin_num;
        col = all_paras->corr_cal_final_data_num/all_paras->theta_bin_num;

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

    }
    
}
