#include"functions_expo_wise.h"

void initialize(char *file_path, data_info *expo_info, int total_expo_num)
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

    // read the infomation of each expo
    read_list(file_path, expo_info, i);


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
    // the biggest block, for each exposure, there're ''zbin_num'' of such blocks in each expo
    expo_info->iz_chi_block_len = expo_info->ir_chi_block_len*expo_info->theta_bin_num;

    expo_info->expo_chi_block_len = expo_info->iz_chi_block_len*expo_info->zbin_num;
    // tangential and cross components
    for(i=0; i<expo_info->total_expo_num; i++)
    {                      
        expo_info->num_count_chit[i] = new double[expo_info->expo_chi_block_len]{};
        expo_info->num_count_chix[i] = new double[expo_info->expo_chi_block_len]{};
    }
    expo_info->total_num_count_chit = new double[expo_info->expo_chi_block_len]{};
    expo_info->total_num_count_chix = new double[expo_info->expo_chi_block_len]{};

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
        // read the array length to get the data column
        if(i == 0)
        {   
            // array length, all elements, of zbin i
            sprintf(set_name, "/data");
            read_h5_datasize(expo_info->expo_name_path[i], set_name, j);
            // data col = data length/gal_num, it's the same for all expos
            expo_info->expo_data_col = j / expo_info->expo_gal_num[0];
        }


        ///////////////// read all the expos  //////////////////////////
        expo_info->expo_zbin_label[i] = new int[expo_info->expo_gal_num[i]]{};
        sprintf(set_name, "/redshift_label");
        read_h5(expo_info->expo_name_path[i], set_name, expo_info->expo_zbin_label[i]);

        expo_info->expo_data[i] = new MY_FLOAT[expo_info->expo_gal_num[i]*expo_info->expo_data_col]{};
        sprintf(set_name, "/data");
        read_h5(expo_info->expo_name_path[i], set_name, expo_info->expo_data[i]);
       
        // read the gal num in each zbin stored in each expo file
        expo_info->expo_num_in_zbin[i] = new int[expo_info->zbin_num];
        sprintf(set_name, "/num_in_zbin");
        read_h5(expo_info->expo_name_path[i], set_name, expo_info->expo_num_in_zbin[i]);


        expo_info->expo_zbin_st[i] = new int[expo_info->zbin_num];
        sprintf(set_name, "/zbin_st");
        read_h5(expo_info->expo_name_path[i], set_name, expo_info->expo_zbin_st[i]);

        expo_info->expo_zbin_ed[i] = new int[expo_info->zbin_num];
        sprintf(set_name, "/zbin_ed");
        read_h5(expo_info->expo_name_path[i], set_name, expo_info->expo_zbin_ed[i]);
    }
}



void initialize_expo_chi_block(data_info *expo_info, int expo_label)
{
    for(int i =0; i< expo_info->expo_chi_block_len; i++)
    {
        expo_info->num_count_chit[expo_label][i] = 0;
        expo_info->num_count_chix[expo_label][i] = 0;
    }
}

void initialize_total_chi_block(data_info *expo_info)
{   
    int i;
    for(i=0; i<expo_info->expo_chi_block_len; i++)
    {   
        expo_info->total_num_count_chit[i] = 0;
        expo_info->total_num_count_chix[i] = 0;
    }
}

void collect_chi_block(data_info *expo_info, int expo_label)
{
    ;
}

void task_distribution(int portion, int my_id, data_info *expo_info)
{   
    int i, j, m, n;
    m = expo_info->expo_pair_num/portion;
    n = expo_info->expo_pair_num%portion;

    for(i=0; i<portion; i++)
    {
        expo_info->expo_pair_num_each_rank[i] = m;
        if(i<n){expo_info->expo_pair_num_each_rank[i]+=1;}
    }
    

    j = 0;
    for(i=0; i<my_id; i++)
    {
        j += expo_info->expo_pair_num_each_rank[i];
    }
    expo_info->my_expo_pair_st = j;
    expo_info->my_expo_pair_ed = j + expo_info->expo_pair_num_each_rank[my_id];

}

void task_prepare(int numprocs, int rank, data_info *expo_info)
{
    int i, ii, j, k;
    int expo_pair_count = 0;
    k = expo_info->total_expo_num*expo_info->total_expo_num;
    int *expo_pair_1 = new int[k]{};
    int *expo_pair_2 = new int[k]{};
    if(k>= 40000){std::cout<<"Too many exposures! Int overflow!"<<std::endl;}

    for(i=0; i<expo_info->total_expo_num; i++)
    {   
        // to avoid double counting the pairs
        for(j=i+1; j<expo_info->total_expo_num; j++)
        {   
            expo_distance(expo_info, i, j, k);
            if(k == 1)
            {
                expo_pair_1[expo_pair_count] = i;
                expo_pair_2[expo_pair_count] = j;
                expo_pair_count += 1;
            } 
        }
        
    }
    expo_info->expo_pair_label_1 = new int[expo_pair_count];
    expo_info->expo_pair_label_2 = new int[expo_pair_count];    
    expo_info->expo_pair_num = expo_pair_count;

    for(i=0;i<expo_pair_count;i++)
    {
        expo_info->expo_pair_label_1[i] = expo_pair_1[i];
        expo_info->expo_pair_label_2[i] = expo_pair_2[i];
    }

    // task distribution
    expo_info->expo_pair_num_each_rank = new int[numprocs]{};
    task_distribution(numprocs, rank, expo_info);

    delete[] expo_pair_1;
    delete[] expo_pair_2;
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
}

void hist_2d(MY_FLOAT x, MY_FLOAT y, MY_FLOAT*bins, int bin_num, int &ix, int &iy)
{
    int i;
    for(i=0; i<bin_num; i++)
    {
        if(x > bins[i] and x <= bins[i+1]){ix=i;break;}
    }
    for(i=0; i<bin_num; i++)
    {
        if(y > bins[i] and y <= bins[i+1]){iy=i;break;}
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

void hist_2d_new(MY_FLOAT*bins, int bin_num, MY_FLOAT former_x, MY_FLOAT former_y, int former_ix, int former_iy, MY_FLOAT x, MY_FLOAT y, int &ix, int &iy)
{
    int im, im1, im2;
    if(x>former_x)
    {
        for(im=former_ix; im<bin_num; im++)
        {if(x >= bins[im] and x < bins[im+1]){ix=im;break;}}
    }
    else
    {
        for(im=former_ix; im > -1; im--)
        {if(x >= bins[im] and x < bins[im+1]){ix=im;break;}}
    }
    if(y>former_y)
    {
        for(im=former_iy; im<bin_num; im++)
        {if(y >= bins[im] and y < bins[im+1]){iy=im;break;}}
    }
    else
    {
        for(im=former_iy; im > -1; im--)
        {if(y >= bins[im] and y < bins[im+1]){iy=im;break;}}
    }
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



// void find_pairs_diff_expo(data_info *expo_info, int expo_label_0, int expo_label_1)
// {
//     int i, j, m, n, k;

//     int ib1, ib2;
//     int expo_label_0;
//     int ir, theta_tag, ic;
//     MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
//     MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;

//     MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
//     MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;
//     MY_FLOAT temp_x, temp_y;
//     int ix, iy, im, im1,im2;
//     int chi_guess_num;
//     int iexpo_chi_block_len, ir_chi_block_len,chi_block_len;
//     int iexpo_len, ir_len, ic_len;
//     MY_FLOAT gg_1, gg_2, gg_len;

//     MY_FLOAT delta_ra, delta_dec, delta_radius;
//     MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

//     int *target_block_label = new int[expo_info->block_num[expo_label_1]]{};
//     int *block_label_mask = new int[expo_info->block_num[expo_label_1]]{};
//     int target_block_num = 0;
//     double st1, st2;
//     double pairs = 0;

//     int loop_label;
//     loop_label = expo_info->loop_label;

//     int mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3;
//     mg_bin_num = expo_info->mg_bin_num;
//     mg_bin_num1 = expo_info->mg_bin_num1;
//     mg_bin_num2 = expo_info->mg_bin_num2;
//     mg_bin_num3 = expo_info->mg_bin_num3;

//     chi_guess_num = expo_info->chi_guess_num;
//     iexpo_chi_block_len = expo_info->iexpo_chi_block_len;
//     ir_chi_block_len = expo_info->ir_chi_block_len;
//     chi_block_len = expo_info->chi_block_len;
//     gg_len = expo_info->gg_len;

//     // decide which the blocks in the target expo, expo_label_1, to calculation
//     // it will calculate the distance of all the block pairs, if the distance is 
//     // smaller than the biggest separation angle, this block in "expo_label_1"
//     // will be labeled as tagert block, it will be included in the calculation.
//     // it will loop all the galaxy pairs which is labeled by the
//     // block, then this pair will distributed to the right separation angle bin.
//     for(ib2=0; ib2<expo_info->block_num[expo_label_1]; ib2++)
//     {
//         block_label_mask[ib2] = 1;
//     }
//     for(ib1=0; ib1<expo_info->block_num[expo_label_0]; ib1++)
//     {
//         ra_z1 = expo_info->block_cen_ra[expo_label_0][ib1];
//         dec_z1 = expo_info->block_cen_dec[expo_label_0][ib1];
//         cos_dec_z1 = expo_info->block_cen_cos_dec[expo_label_0][ib1];
//         delta_len_z1 = expo_info->block_delta_len[expo_label_0][ib1];

//         for(ib2=0; ib2<expo_info->block_num[expo_label_1]; ib2++)
//         {
//             ra_z2 = expo_info->block_cen_ra[expo_label_1][ib2];
//             dec_z2 = expo_info->block_cen_dec[expo_label_1][ib2];
//             cos_dec_z2 = expo_info->block_cen_cos_dec[expo_label_1][ib2];
//             delta_len_z2 = expo_info->block_delta_len[expo_label_1][ib2];

//             // the seperation angle (arc minute)
//             delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
//             delta_dec = dec_z2 - dec_z1;
//             delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec) - delta_len_z1 - delta_len_z2;

//             if(delta_radius <= expo_info->theta_bin[expo_info->theta_bin_num] and block_label_mask[ib2] == 1)
//             {
//                 target_block_label[target_block_num] = ib2;
//                 block_label_mask[ib2] = 0;
//                 target_block_num += 1;
//             }
//         }
//     }
    
//     for(ib1=0; ib1<expo_info->block_num[expo_label_0]; ib1++)
//     {   
//         // st1 = clock();
//         for(i=expo_info->block_st_z1[expo_label_0][ib1]; i<expo_info->block_ed_z1[expo_label_0][ib1]; i++)
//         {
//             // loop the grid in the first zbin, zbin_label_0
//             m = i*expo_info->expo_data_col;

//             ra_z1 = expo_info->expo_data_z1[expo_label_0][m+expo_info->ra_idx];
//             dec_z1 = expo_info->expo_data_z1[expo_label_0][m+expo_info->dec_idx];
//             cos_dec_z1 = expo_info->expo_data_z1[expo_label_0][m+expo_info->cos_dec_idx];

//             mg1_z1 = expo_info->expo_data_z1[expo_label_0][m+expo_info->mg1_idx];
//             mg2_z1 = expo_info->expo_data_z1[expo_label_0][m+expo_info->mg2_idx];

//             mnu1_z1 = expo_info->expo_data_z1[expo_label_0][m+expo_info->mn_idx] +
//                         expo_info->expo_data_z1[expo_label_0][m+expo_info->mu_idx];
//             mnu2_z1 = expo_info->expo_data_z1[expo_label_0][m+expo_info->mn_idx] -
//                         expo_info->expo_data_z1[expo_label_0][m+expo_info->mu_idx];

//             expo_label_0 = expo_info->expo_expo_label_z1[expo_label_0][i];
            
//             iexpo_len = expo_label_0*iexpo_chi_block_len;

//             for(k=0; k<target_block_num; k++)
//             {
//                 ib2=target_block_label[k];

//                 for(j=expo_info->block_st_z2[expo_label_1][ib2]; j<expo_info->block_ed_z2[expo_label_1][ib2]; j++)
//                 {   
//                     // loop the grid in the second zbin, zbin_label_1
//                     n = j*expo_info->expo_data_col;

//                     ra_z2 = expo_info->expo_data_z2[expo_label_1][n+expo_info->ra_idx];
//                     dec_z2 = expo_info->expo_data_z2[expo_label_1][n+expo_info->dec_idx];
//                     cos_dec_z2 = expo_info->expo_data_z2[expo_label_1][n+expo_info->cos_dec_idx];

//                     // the seperation angle (arc minute)
//                     delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
//                     delta_dec = dec_z2 - dec_z1;
//                     delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);
                    
//                     theta_tag = -1;
//                     for(ir=0; ir<expo_info->theta_bin_num; ir++)
//                     {
//                         if(delta_radius > expo_info->theta_bin[ir] and delta_radius <= expo_info->theta_bin[ir+1])
//                         {theta_tag=ir;break;}
//                     }
//                     // std::cout<<delta_radius<<" "<<expo_info->theta_bin[theta_tag]<<" "<<expo_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
//                     if(theta_tag > -1)
//                     {   
//                         // pairs+= 1;

//                         // shear estimators rotation (position angle defined as East of North)
//                         sin_theta = delta_ra/delta_radius;
//                         cos_theta = delta_dec/delta_radius;

//                         sin_2theta = 2*sin_theta*cos_theta;
//                         cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

//                         sin_4theta = 2*sin_2theta*cos_2theta;
//                         cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;


//                         mg1_z2 = expo_info->expo_data_z2[expo_label_1][n+expo_info->mg1_idx]*cos_2theta - 
//                                 expo_info->expo_data_z2[expo_label_1][n+expo_info->mg2_idx]*sin_2theta;
//                         mg2_z2 = expo_info->expo_data_z2[expo_label_1][n+expo_info->mg1_idx]*sin_2theta + 
//                                 expo_info->expo_data_z2[expo_label_1][n+expo_info->mg2_idx]*cos_2theta;

//                         mnu1_z2 = expo_info->expo_data_z2[expo_label_1][n+expo_info->mu_idx]*cos_4theta -
//                                 expo_info->expo_data_z2[expo_label_1][n+expo_info->mv_idx]*sin_4theta;
//                         mnu2_z2 = mnu1_z2;

//                         mnu1_z2 = expo_info->expo_data_z2[expo_label_1][n+expo_info->mn_idx] + mnu2_z2;
//                         mnu2_z2 = expo_info->expo_data_z2[expo_label_1][n+expo_info->mn_idx] - mnu2_z2;
                        
//                         // the key part of PDF_SYM
//                         ir_len = theta_tag*ir_chi_block_len + iexpo_len;

//                         for(ic=0; ic<chi_guess_num; ic++)
//                         {   
//                             if(loop_label >= gg_len){loop_label = 0;}

//                             ic_len = ic*chi_block_len + ir_len;
//                             gg_1 = expo_info->gg_1[loop_label];
//                             gg_2 = expo_info->gg_2[loop_label];
//                             // std::cout<<theta_tag<<" "<<expo_label_0<<" "<<ic_len<<std::endl;
//                             // std::cout<<loop_label<<std::endl;
//                             temp_x = mg1_z1 - gg_1*mnu1_z1;
//                             temp_y = mg1_z2 - gg_2*mnu1_z2;
                            
//                             // hist_2d_fast(temp_x, temp_y, expo_info->mg_bin, expo_info->mg_bin_num, expo_info->mg_bin_num2,ix, iy);
//                             hist_2d_new(temp_x, temp_y, expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix, iy);
//                             // std::cout<<iy<<" "<<ix<<std::endl;

//                             expo_info->num_count_chit[expo_label_0][ic_len + iy*mg_bin_num+ix] += 1;
                            

//                             temp_x = mg2_z1 - gg_1*mnu2_z1;
//                             temp_y = mg2_z2 - gg_2*mnu2_z2;

//                             // hist_2d_fast(temp_x, temp_y, expo_info->mg_bin, expo_info->mg_bin_num, expo_info->mg_bin_num2,ix, iy);
//                             hist_2d_new(temp_x, temp_y, expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix, iy);

//                             // std::cout<<iy<<" "<<ix<<std::endl;

//                             expo_info->num_count_chix[expo_label_0][ic_len + iy*mg_bin_num+ix] += 1;
//                             // std::cout<<iy<<" "<<ix<<std::endl;
//                             loop_label += 1;
//                         }
//                     }
//                 }
//             }
//         }
//         // st2 = clock();
//         // std::cout<<"Block "<<ib1<<" "<<pairs<<" pairs "<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
//         // pairs = 0;
//     }
//     expo_info->loop_label = loop_label;

//     delete[] block_label_mask;
//     delete[] target_block_label;
// }

// void find_pairs_same_expo(data_info *expo_info, int expo_label)
// {
//     int i, j, m, n, k;

//     int ib1, ib2;
//     int expo_label_0;
//     int ir, theta_tag, ic;
//     MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
//     MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;

//     MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
//     MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;

//     MY_FLOAT temp_x_tt, temp_y_tt, temp_x_xx, temp_y_xx;
//     MY_FLOAT former_temp_x_tt, former_temp_y_tt,former_temp_x_xx, former_temp_y_xx;
//     int ix_tt, iy_tt, ix_xx, iy_xx;
//     int former_ix_tt, former_iy_tt, former_ix_xx, former_iy_xx;
//     MY_FLOAT temp_tt[4], temp_xx[4];
//     int ix, iy, im, im1,im2;


//     int chi_guess_num;
//     int iexpo_chi_block_len, ir_chi_block_len,chi_block_len;
//     int iexpo_len, ir_len, ic_len;
//     MY_FLOAT gg_1, gg_2;
//     int gg_len;
//     MY_FLOAT *mg_bin;
//     MY_FLOAT delta_ra, delta_dec, delta_radius;
//     MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

//     double st1, st2;

//     double pairs = 0;
    
//     int loop_label;
//     loop_label = expo_info->loop_label;

//     int mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3;
//     mg_bin_num = expo_info->mg_bin_num;
//     mg_bin_num1 = expo_info->mg_bin_num1;
//     mg_bin_num2 = expo_info->mg_bin_num2;
//     mg_bin_num3 = expo_info->mg_bin_num3;

//     int bin_para_tt[2], bin_para_xx[2];

//     chi_guess_num = expo_info->chi_guess_num;
//     iexpo_chi_block_len = expo_info->iexpo_chi_block_len;
//     ir_chi_block_len = expo_info->ir_chi_block_len;
//     chi_block_len = expo_info->chi_block_len;
//     gg_len = expo_info->gg_len;

//     for(ib1=0; ib1<expo_info->block_num[expo_label]; ib1++)
//     {   
//         // st1 = clock();
//         for(i=expo_info->block_st_z1[expo_label][ib1]; i<expo_info->block_ed_z1[expo_label][ib1]; i++)
//         {
//             // loop the grid in the first zbin, zbin_label_0
//             m = i*expo_info->expo_data_col;

//             ra_z1 = expo_info->expo_data_z1[expo_label][m+expo_info->ra_idx];
//             dec_z1 = expo_info->expo_data_z1[expo_label][m+expo_info->dec_idx];
//             cos_dec_z1 = expo_info->expo_data_z1[expo_label][m+expo_info->cos_dec_idx];

//             mg1_z1 = expo_info->expo_data_z1[expo_label][m+expo_info->mg1_idx];
//             mg2_z1 = expo_info->expo_data_z1[expo_label][m+expo_info->mg2_idx];

//             mnu1_z1 = expo_info->expo_data_z1[expo_label][m+expo_info->mn_idx] +
//                         expo_info->expo_data_z1[expo_label][m+expo_info->mu_idx];
//             mnu2_z1 = expo_info->expo_data_z1[expo_label][m+expo_info->mn_idx] -
//                         expo_info->expo_data_z1[expo_label][m+expo_info->mu_idx];

//             expo_label_0 = expo_info->expo_expo_label_z1[expo_label][i];
            
//             iexpo_len = expo_label_0*iexpo_chi_block_len;

//             if(expo_info->zbin_label_0 != expo_info->zbin_label_1){k=0;}
//             else{k=ib1;}

//             for(ib2=k; ib2<expo_info->block_num[expo_label]; ib2++)
//             {
//                 for(j=expo_info->block_st_z2[expo_label][ib2]; j<expo_info->block_ed_z2[expo_label][ib2]; j++)
//                 {   
//                     if(expo_info->expo_expo_label_z2[expo_label][j] !=  expo_label_0)
//                     {
//                         // loop the grid in the second zbin, zbin_label_1
//                         n = j*expo_info->expo_data_col;

//                         ra_z2 = expo_info->expo_data_z2[expo_label][n+expo_info->ra_idx];
//                         dec_z2 = expo_info->expo_data_z2[expo_label][n+expo_info->dec_idx];
//                         cos_dec_z2 = expo_info->expo_data_z2[expo_label][n+expo_info->cos_dec_idx];

//                         // the seperation angle (arc minute)
//                         delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
//                         delta_dec = dec_z2 - dec_z1;
//                         delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);
                        
//                         theta_tag = -1;
//                         for(ir=0; ir<expo_info->theta_bin_num; ir++)
//                         {
//                             if(delta_radius > expo_info->theta_bin[ir] and delta_radius <= expo_info->theta_bin[ir+1])
//                             {theta_tag=ir;break;}
//                         }
//                         // std::cout<<delta_radius<<" "<<expo_info->theta_bin[theta_tag]<<" "<<expo_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
//                         if(theta_tag > -1)
//                         {   
//                             pairs += 1;
//                             // shear estimators rotation (position angle defined as East of North)
//                             sin_theta = delta_ra/delta_radius;
//                             cos_theta = delta_dec/delta_radius;

//                             sin_2theta = 2*sin_theta*cos_theta;
//                             cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

//                             sin_4theta = 2*sin_2theta*cos_2theta;
//                             cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;


//                             mg1_z2 = expo_info->expo_data_z2[expo_label][n+expo_info->mg1_idx]*cos_2theta - 
//                                     expo_info->expo_data_z2[expo_label][n+expo_info->mg2_idx]*sin_2theta;
//                             mg2_z2 = expo_info->expo_data_z2[expo_label][n+expo_info->mg1_idx]*sin_2theta + 
//                                     expo_info->expo_data_z2[expo_label][n+expo_info->mg2_idx]*cos_2theta;

//                             mnu1_z2 = expo_info->expo_data_z2[expo_label][n+expo_info->mu_idx]*cos_4theta -
//                                     expo_info->expo_data_z2[expo_label][n+expo_info->mv_idx]*sin_4theta;
//                             mnu2_z2 = mnu1_z2;

//                             mnu1_z2 += expo_info->expo_data_z2[expo_label][n+expo_info->mn_idx];
//                             mnu2_z2 = expo_info->expo_data_z2[expo_label][n+expo_info->mn_idx] - mnu2_z2;
                            
//                             ////////////////////// the key part of PDF_SYM //////////////////////////////
//                             ir_len = theta_tag*ir_chi_block_len + iexpo_len;

//                             gg_1 = expo_info->gg_1[loop_label];// expo_info->gg_1[0][loop_label];
//                             gg_2 = expo_info->gg_2[loop_label];// expo_info->gg_2[0][loop_label];


//                             temp_tt[2] = mg1_z1 - gg_1*mnu1_z1;
//                             temp_tt[3] = mg1_z2 - gg_2*mnu1_z2;
//                             hist_2d_new(temp_tt[2], temp_tt[3], expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_tt, iy_tt);
                            
//                             // std::cout<<0<<" "<<temp_x_tt<<" "<<temp_y_tt<<" "<<ix_tt<<" "<<iy_tt<<std::endl;
//                             expo_info->num_count_chit[expo_label][ir_len + iy_tt*mg_bin_num+ix_tt] += 1;
                            
//                             temp_xx[2] = mg2_z1 - gg_1*mnu2_z1;
//                             temp_xx[3] = mg2_z2 - gg_2*mnu2_z2;

//                             hist_2d_new(temp_xx[2], temp_xx[3],  expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_xx, iy_xx);
//                             expo_info->num_count_chix[expo_label][ir_len + iy_xx*mg_bin_num+ix_xx] += 1;
//                             loop_label += 1;
//                             // std::cout<<ic<<" "<<temp_x_xx<<" "<<temp_y_xx<<" "<<ix_xx<<" "<<iy_xx<<std::endl;
//                             for(ic=1; ic<chi_guess_num; ic++)
//                             {   
//                                 ic_len = ic*chi_block_len + ir_len;

//                                 gg_1 = expo_info->gg_1[loop_label];
//                                 gg_2 = expo_info->gg_2[loop_label];

//                                 bin_para_tt[0] = ix_tt;
//                                 bin_para_tt[1] = iy_tt;

//                                 temp_tt[0] = temp_tt[2];
//                                 temp_tt[1] = temp_tt[3];

//                                 temp_tt[2] = mg1_z1 - gg_1*mnu1_z1;
//                                 temp_tt[3] = mg1_z2 - gg_2*mnu1_z2;
//                                 // hist_2d_new(temp_x_tt, temp_y_tt, expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_tt, iy_tt);
//                                 hist_2d_new(expo_info->mg_bin, mg_bin_num, temp_tt, bin_para_tt, ix_tt, iy_tt);

//                                 expo_info->num_count_chit[expo_label][ic_len + iy_tt*mg_bin_num+ix_tt] += 1;
                                
//                                 bin_para_xx[0] = ix_xx;
//                                 bin_para_xx[1] = iy_xx;
                                
//                                 temp_xx[0] = temp_xx[2];
//                                 temp_xx[1] = temp_xx[3];

//                                 temp_xx[2] = mg2_z1 - gg_1*mnu2_z1;
//                                 temp_xx[3] = mg2_z2 - gg_2*mnu2_z2;

//                                 // hist_2d_new(temp_x_xx, temp_y_xx,  expo_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix_xx, iy_xx);
//                                 hist_2d_new(expo_info->mg_bin, mg_bin_num, temp_xx, bin_para_xx, ix_xx, iy_xx);
//                                 expo_info->num_count_chix[expo_label][ic_len + iy_xx*mg_bin_num+ix_xx] += 1;

//                                 loop_label += 1;
//                             }
//                             if(loop_label >= gg_len){loop_label = 0;}
//                             ////////////////////// the key part of PDF_SYM -end  //////////////////////////////
//                         }
//                     }
//                     else
//                     {
//                         continue;
//                     }
//                 }
//             }
//         }
//         // st2 = clock();
//         // std::cout<<"Block "<<ib1<<", "<<pairs<<" pairs, "<<(st2-st1)/CLOCKS_PER_SEC<<" Sec"<<std::endl;
//         // pairs = 0;
//     }
//     expo_info->loop_label = loop_label;
// }

// void find_pairs_same_expo(data_info *expo_info, int expo_label)
// {

// }
// void find_pairs_diff_expo(data_info *expo_info, int expo_label_0, int expo_label_1)
// {
    
// }

