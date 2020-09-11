#include"functions_expo_wise.h"

void initialize(data_info *expo_info, int total_expo_num)
{
    int i, j;
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
    // tangential and cross components
    expo_info->expo_num_count_chit = new double[expo_info->expo_chi_block_len]{};
    expo_info->expo_num_count_chix = new double[expo_info->expo_chi_block_len]{};

    expo_info->theta_accum_len = expo_info->theta_bin_num*expo_info->zbin_num*expo_info->zbin_num;
    expo_info->theta_accum = new double[expo_info->theta_accum_len];
    expo_info->theta_num_accum = new double[expo_info->theta_accum_len];

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
    int expo_pair_count = 0;
    k = expo_info->total_expo_num*expo_info->total_expo_num;

    
    expo_info->expo_pair_num_each_rank = new int[numprocs]{};
    m = expo_info->total_expo_num/numprocs;
    n = expo_info->total_expo_num%numprocs;

    for(i=0; i<numprocs; i++)
    {
        expo_info->expo_pair_num_each_rank[i] = m;
        if(i<n){expo_info->expo_pair_num_each_rank[i]+=1;}
    }
    
    expo_info->task_expo_num = expo_info->expo_pair_num_each_rank[rank];
    expo_info->task_expo_label = new int[expo_info->task_expo_num];
    

    for(i=0;i<m;i++)
    {
        expo_info->task_expo_label[i] = i*numprocs + rank;
    }
    if(rank < n)
    {
        expo_info->task_expo_label[m] = m*numprocs + rank;
    }
    // task distribution

    // task_distribution(numprocs, rank, expo_info);
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
    int row, col;
    char result_path[600], set_name[50];

    col = expo_info->mg_bin_num;
    row = expo_info->expo_chi_block_len/col;

    sprintf(result_path, "%s/result/core_%s_num_count.hdf5", expo_info->parent_path, rank);
    sprintf(set_name, "/%d-%d/tt",expo_label_1, expo_label_2);
    
    if(expo_info->result_file_tag == 0)
    {
        write_h5(result_path, set_name, expo_info->expo_num_count_chit, row, col, true);
        expo_info->result_file_tag=1;
    }
    else{write_h5(result_path, set_name, expo_info->expo_num_count_chit, row, col, false);}

    sprintf(set_name, "/%d-%d/xx",expo_label_1, expo_label_2);
    write_h5(result_path, set_name, expo_info->expo_num_count_chix, row, col, false);

    col = expo_info->theta_bin_num;
    row = expo_info->zbin_num*expo_info->zbin_num;
    sprintf(set_name, "/%d-%d/theta",expo_label_1, expo_label_2);
    write_h5(result_path, set_name, expo_info->theta_accum, row, col, false);
    sprintf(set_name, "/%d-%d/theta_num",expo_label_1, expo_label_2);
    write_h5(result_path, set_name, expo_info->theta_num_accum, row, col, false);
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




