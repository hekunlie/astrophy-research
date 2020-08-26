#include"functions.h"

void read_inform(char *file_path, data_info *field_info, int &read_file_num)
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

        strs >> field_info->field_name_path[line_count];
        strs >> field_info->field_name[line_count];

        strs >> field_info->exposure_num_of_field[line_count];

        strs >> field_info->field_cen_ra[line_count];
        strs >> field_info->field_cen_dec[line_count];

        strs >> field_info->field_delta_ra[line_count];
        strs >> field_info->field_delta_dec[line_count];
        strs >> field_info->field_delta_len[line_count];
        strs >> field_info->field_cen_cos_dec[line_count];

        // std::cout << str << std::endl;
        line_count += 1;
    }
    read_file_num = line_count;

}

void read_field_data(data_info *field_info)
{
    int i, j;
    char set_name[100];
    for(i=0; i<field_info->total_field_num; i++)
    {   
        // std::cout<<"Read "<<field_info->field_name_path[i]<<std::endl;
        // read the gal num in each zbin stored in each field file
        sprintf(set_name, "/total_num_in_zbin");
        read_h5(field_info->field_name_path[i], set_name, field_info->num_in_zbin);
        field_info->total_gal_num_z1[i] = field_info->num_in_zbin[field_info->zbin_label_0];
        field_info->total_gal_num_z2[i] = field_info->num_in_zbin[field_info->zbin_label_1];

        // read the array length to get the data column
        if(i == 0)
        {   
            // array length, all elements, of zbin i
            sprintf(set_name, "/z%d/field",field_info->zbin_label_0);
            read_h5_datasize(field_info->field_name_path[i], set_name, j);
            // data col = data length/gal_num, it's the same for all fields
            field_info->field_data_col = j / field_info->num_in_zbin[field_info->zbin_label_0];
        }


        ///////////////// read all the fields  //////////////////////////
        // zbin i
        j = field_info->num_in_zbin[field_info->zbin_label_0]*field_info->field_data_col;
        field_info->field_data_z1[i] = new MY_FLOAT[j]{};
        sprintf(set_name, "/z%d/field",field_info->zbin_label_0);
        read_h5(field_info->field_name_path[i], set_name, field_info->field_data_z1[i]);
        // zbin j
        j = field_info->num_in_zbin[field_info->zbin_label_1]*field_info->field_data_col;
        field_info->field_data_z2[i] = new MY_FLOAT[j]{};
        sprintf(set_name, "/z%d/field",field_info->zbin_label_1);
        read_h5(field_info->field_name_path[i], set_name, field_info->field_data_z2[i]);

        field_info->field_expo_label_z1[i] = new int[field_info->total_gal_num_z1[i]]{};
        sprintf(set_name, "/z%d/exposure_label", field_info->zbin_label_0);
        read_h5(field_info->field_name_path[i], set_name, field_info->field_expo_label_z1[i]);
        
        field_info->field_expo_label_z2[i] = new int[field_info->total_gal_num_z2[i]]{};
        sprintf(set_name, "/z%d/exposure_label", field_info->zbin_label_1);
        read_h5(field_info->field_name_path[i], set_name, field_info->field_expo_label_z2[i]);


        //////////// read the block inform in each field  //////////////////////////
        // the block labels 
        field_info->field_block_label_z1[i] = new int[field_info->total_gal_num_z1[i]]{};
        sprintf(set_name, "/z%d/block_label", field_info->zbin_label_0);
        read_h5(field_info->field_name_path[i], set_name, field_info->field_block_label_z1[i]);
        
        field_info->field_block_label_z2[i] = new int[field_info->total_gal_num_z2[i]]{};
        sprintf(set_name, "/z%d/block_label", field_info->zbin_label_1);
        read_h5(field_info->field_name_path[i], set_name, field_info->field_block_label_z2[i]);

        if(i==0){ field_info->block_num = new int[field_info->total_field_num]{}; }
        // number of block in each field
        sprintf(set_name, "/block_cen_ra");
        read_h5_datasize(field_info->field_name_path[i], set_name, j);
        field_info->block_num[i] = j;
        // RA of the center of each block
        field_info->block_cen_ra[i] = new MY_FLOAT[j];
        read_h5(field_info->field_name_path[i], set_name, field_info->block_cen_ra[i]);

        // Dec of the center of each block
        field_info->block_cen_dec[i] = new MY_FLOAT[j];
        sprintf(set_name, "/block_cen_dec");
        read_h5(field_info->field_name_path[i], set_name, field_info->block_cen_dec[i]);
        
        // cos(Dec) of the center of each block
        field_info->block_cen_cos_dec[i] = new MY_FLOAT[j];
        sprintf(set_name, "/block_cen_cos_dec");
        read_h5(field_info->field_name_path[i], set_name, field_info->block_cen_cos_dec[i]);

        // cos(Dec) of the center of each block
        field_info->block_delta_len[i] = new MY_FLOAT[j];
        sprintf(set_name, "/block_delta_len");
        read_h5(field_info->field_name_path[i], set_name, field_info->block_delta_len[i]);


        // block start & end 
        field_info->block_st_z1[i] = new MY_FLOAT[j];
        field_info->block_st_z2[i] = new MY_FLOAT[j];
        field_info->block_ed_z1[i] = new MY_FLOAT[j];
        field_info->block_ed_z2[i] = new MY_FLOAT[j];
        sprintf(set_name, "/z%d/block_st", field_info->zbin_label_0);
        read_h5(field_info->field_name_path[i], set_name, field_info->block_st_z1[i]);
        sprintf(set_name, "/z%d/block_st", field_info->zbin_label_1);
        read_h5(field_info->field_name_path[i], set_name, field_info->block_st_z2[i]);
        sprintf(set_name, "/z%d/block_ed", field_info->zbin_label_0);
        read_h5(field_info->field_name_path[i], set_name, field_info->block_ed_z1[i]);
        sprintf(set_name, "/z%d/block_ed", field_info->zbin_label_1);
        read_h5(field_info->field_name_path[i], set_name, field_info->block_ed_z2[i]);
    }
}


void initialize(char *file_path, data_info *field_info, int total_field_num)
{
    int i, j;
    char set_name[100];
    char data_path[600];

    field_info->total_field_num = total_field_num;

    for(i=0;i<total_field_num;i++)
    {   
        // the field file directories
        field_info->field_name_path[i] = new char[400];
        field_info->field_name[i] = new char[50];  
    }

    field_info->exposure_num_of_field = new int[total_field_num]{};
     
    field_info->field_cen_ra = new MY_FLOAT[total_field_num]{};  
    field_info->field_cen_dec = new MY_FLOAT[total_field_num]{}; 

    field_info->field_delta_ra = new MY_FLOAT[total_field_num]{};  
    field_info->field_delta_dec = new MY_FLOAT[total_field_num]{};
    field_info->field_delta_len = new MY_FLOAT[total_field_num]{};
    field_info->field_cen_cos_dec = new MY_FLOAT[total_field_num]{};

    field_info->total_gal_num_z1 = new int[total_field_num]{};
    field_info->total_gal_num_z2 = new int[total_field_num]{};
    

    // read the infomation of each field
    read_inform(file_path, field_info, i);


    ///////////////// read the inform of the PDF_SYM  ////////////////
    sprintf(data_path,"%s/cata/gg_cor.hdf5", field_info->parent_path);
    // read radius bin
    sprintf(set_name,"/theta_bin");
    read_h5_datasize(data_path, set_name,j);
    field_info->theta_bin = new MY_FLOAT[j]{};
    field_info->theta_bin_num = j -1;
    read_h5(data_path, set_name, field_info->theta_bin);

    // read redshift bin
    sprintf(set_name,"/redshift_bin");
    read_h5_datasize(data_path, set_name,j);
    field_info->zbin = new MY_FLOAT[j]{};
    field_info->zbin_num = j -1;
    field_info->num_in_zbin = new int[j]{};
    read_h5(data_path, set_name, field_info->zbin);


    sprintf(set_name, "/chi_guess");
    read_h5_datasize(data_path, set_name, field_info->chi_guess_num);
    field_info->chi_guess = new MY_FLOAT[field_info->chi_guess_num];
    read_h5(data_path, set_name, field_info->chi_guess);

    sprintf(set_name, "/mg_bin");
    read_h5_datasize(data_path, set_name, j);
    field_info->mg_bin = new MY_FLOAT[j];
    field_info->mg_bin_num = j -1;
    field_info->mg_bin_num2 = (j-1)/2;
    field_info->mg_bin_num1 = field_info->mg_bin_num2/2;
    field_info->mg_bin_num3 = field_info->mg_bin_num1 + field_info->mg_bin_num2;

    read_h5(data_path, set_name, field_info->mg_bin);

    // the fundmental block size for number counting in the PDF_SYM
    // the smallest block, mg_bin_num x mg_bin_num, for each chi guess point
    field_info->chi_block_len = field_info->mg_bin_num*field_info->mg_bin_num;
    // the upper level of the above block, chi_guess_num x mg_bin_num x mg_bin_num,
    // for each theta_bin, there're ''theta_bin_num'' of such blocks in each exposure
    field_info->ir_chi_block_len = field_info->chi_guess_num*field_info->chi_block_len;
    // the biggest block, for each exposure, there're ''exposure_num'' of such blocks in each field
    field_info->iexpo_chi_block_len = field_info->theta_bin_num*field_info->ir_chi_block_len;

    field_info->field_chi_block_len = new int[total_field_num]{};
    // tangential and cross components
    for(i=0; i<field_info->total_field_num; i++)
    {                      
        field_info->field_chi_block_len[i] = field_info->iexpo_chi_block_len*field_info->exposure_num_of_field[i];
        field_info->num_count_chit[i] = new double[field_info->field_chi_block_len[i]]{};
        field_info->num_count_chix[i] = new double[field_info->field_chi_block_len[i]]{};
    }
    field_info->total_num_count_chit = new double[field_info->iexpo_chi_block_len]{};
    field_info->total_num_count_chix = new double[field_info->iexpo_chi_block_len]{};

    // read the correlated shear pairs generated before for time-saving
    for(i=0; i<field_info->chi_guess_num; i++)
    {
        sprintf(set_name,"/%d/g11",i);
        read_h5_datasize(data_path, set_name, field_info->gg_len);

        field_info->gg_1[i] = new MY_FLOAT[field_info->gg_len]{};
        field_info->gg_2[i] = new MY_FLOAT[field_info->gg_len]{};
        
        sprintf(set_name,"/%d/g11",i);
        read_h5(data_path, set_name, field_info->gg_1[i]);
        sprintf(set_name,"/%d/g22",i);
        read_h5(data_path, set_name, field_info->gg_2[i]);
    }
    field_info->loop_label = 0;
}

void initialize_field_chi_block(data_info *field_info, int field_label)
{
    for(int i=0; i<field_info->field_chi_block_len[field_label]; i++)
    {
        field_info->num_count_chit[field_label][i] = 0;
    }
}

void initialize_total_chi_block(data_info *field_info)
{   
    int i;
    
    for(i=0; i<field_info->iexpo_chi_block_len; i++)
    {   
        field_info->total_num_count_chit[i] = 0;
        field_info->total_num_count_chix[i] = 0;
    }
}

void collect_chi_block(data_info *field_info, int field_label)
{
    int i,j, k;
    for(i=0;i<field_info->exposure_num_of_field[field_label]; i++)
    {   
        k = i*field_info->iexpo_chi_block_len;
        for(j=0; j<field_info->iexpo_chi_block_len; j++)
        {
            field_info->total_num_count_chit[j] = field_info->num_count_chit[field_label][j + k];
            field_info->total_num_count_chix[j] = field_info->num_count_chix[field_label][j + k];
        }
    }
}

void task_distribution(int portion, int my_id, data_info *field_info)
{   
    int i, j, m, n;
    m = field_info->field_pair_num/portion;
    n = field_info->field_pair_num%portion;

    for(i=0; i<portion; i++)
    {
        field_info->field_pair_num_each_rank[i] = m;
        if(i<n){field_info->field_pair_num_each_rank[i]+=1;}
    }
    

    j = 0;
    for(i=0; i<my_id; i++)
    {
        j += field_info->field_pair_num_each_rank[i];
    }
    field_info->my_field_pair_st = j;
    field_info->my_field_pair_ed = j + field_info->field_pair_num_each_rank[my_id];

}

void task_prepare(int numprocs, int rank, data_info *field_info)
{
    int i, ii, j, k;
    int field_pair_count = 0;
    k = field_info->total_field_num*field_info->total_field_num;
    int *field_pair_1 = new int[k]{};
    int *field_pair_2 = new int[k]{};
    
    for(i=0; i<field_info->total_field_num; i++)
    {   
        // to avoid double counting the pairs
        ii = i;
        if(field_info->zbin_label_0 != field_info->zbin_label_1){ii = 0;}

        if(field_info->total_gal_num_z1[i] > 0)
        {
            for(j=ii; j<field_info->total_field_num; j++)
            {   
                if(field_info->total_gal_num_z2[j] > 0)
                {
                    field_distance(field_info, i, j, k);
                    if(k == 1)
                    {
                        field_pair_1[field_pair_count] = i;
                        field_pair_2[field_pair_count] = j;
                        field_pair_count += 1;
                    }                  
                }
            }
        }
    }
    field_info->field_pair_label_1 = new int[field_pair_count];
    field_info->field_pair_label_2 = new int[field_pair_count];    
    field_info->field_pair_num = field_pair_count;

    for(i=0;i<field_pair_count;i++)
    {
        field_info->field_pair_label_1[i] = field_pair_1[i];
        field_info->field_pair_label_2[i] = field_pair_2[i];
    }

    // task distribution
    field_info->field_pair_num_each_rank = new int[numprocs]{};
    task_distribution(numprocs, rank, field_info);

    delete[] field_pair_1;
    delete[] field_pair_2;
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

void field_distance(data_info *field_info, int field_label_0, int field_label_1, int &label)
{
    MY_FLOAT ra_1, dec_1, cos_dec_1, delta_len_1;
    MY_FLOAT ra_2, dec_2, cos_dec_2, delta_len_2;
    MY_FLOAT delta_radius, delta_ra, delta_dec;

    ra_1 = field_info->field_cen_ra[field_label_0];
    dec_1 = field_info->field_cen_dec[field_label_0];
    cos_dec_1 = field_info->field_cen_cos_dec[field_label_0];
    delta_len_1 = field_info->field_delta_len[field_label_0];

    ra_2 = field_info->field_cen_ra[field_label_1];
    dec_2 = field_info->field_cen_dec[field_label_1];
    cos_dec_2 = field_info->field_cen_cos_dec[field_label_1];
    delta_len_2 = field_info->field_delta_len[field_label_1];

    // the seperation angle (arc minute)
    delta_ra = (ra_2 - ra_1)*cos_dec_1;
    delta_dec = dec_2 - dec_1;
    delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec) - delta_len_1 - delta_len_2;
    
    label = 0;
    if(delta_radius <= field_info->theta_bin[field_info->theta_bin_num])
    {   
        // std::cout<<ra_1<<" "<<dec_1<<std::endl;
        // std::cout<<ra_2<<" "<<dec_2<<" "<<delta_radius<<" "<< field_info->theta_bin[field_info->theta_bin_num]<<std::endl;
        label = 1;
    }
}

void find_pairs_diff_field(data_info *field_info, int field_label_0, int field_label_1)
{
    int i, j, m, n, k;

    int ib1, ib2;
    int expo_label_0;
    int ir, theta_tag, ic;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;

    MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
    MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;
    MY_FLOAT temp_x, temp_y;
    int ix, iy, im, im1,im2;
    int iexpo_len, ir_len, ic_len;
    MY_FLOAT gg_1, gg_2;

    MY_FLOAT delta_ra, delta_dec, delta_radius;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    int *target_block_label = new int[field_info->block_num[field_label_1]]{};
    int *block_label_mask = new int[field_info->block_num[field_label_1]]{};
    int target_block_num = 0;
    double st1, st2;
    double pairs = 0;

    int loop_label;
    loop_label = field_info->loop_label;
    int mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3;
    mg_bin_num = field_info->mg_bin_num;
    mg_bin_num1 = field_info->mg_bin_num1;
    mg_bin_num2 = field_info->mg_bin_num2;
    mg_bin_num3 = field_info->mg_bin_num3;

    // decide which the blocks in the target field, field_label_1, to calculation
    // it will calculate the distance of all the block pairs, if the distance is 
    // smaller than the biggest separation angle, this block in "field_label_1"
    // will be labeled as tagert block, it will be included in the calculation.
    // it will loop all the galaxy pairs which is labeled by the
    // block, then this pair will distributed to the right separation angle bin.
    for(ib2=0; ib2<field_info->block_num[field_label_1]; ib2++)
    {
        block_label_mask[ib2] = 1;
    }
    for(ib1=0; ib1<field_info->block_num[field_label_0]; ib1++)
    {
        ra_z1 = field_info->block_cen_ra[field_label_0][ib1];
        dec_z1 = field_info->block_cen_dec[field_label_0][ib1];
        cos_dec_z1 = field_info->block_cen_cos_dec[field_label_0][ib1];
        delta_len_z1 = field_info->block_delta_len[field_label_0][ib1];

        for(ib2=0; ib2<field_info->block_num[field_label_1]; ib2++)
        {
            ra_z2 = field_info->block_cen_ra[field_label_1][ib2];
            dec_z2 = field_info->block_cen_dec[field_label_1][ib2];
            cos_dec_z2 = field_info->block_cen_cos_dec[field_label_1][ib2];
            delta_len_z2 = field_info->block_delta_len[field_label_1][ib2];

            // the seperation angle (arc minute)
            delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
            delta_dec = dec_z2 - dec_z1;
            delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec) - delta_len_z1 - delta_len_z2;

            if(delta_radius <= field_info->theta_bin[field_info->theta_bin_num] and block_label_mask[ib2] == 1)
            {
                target_block_label[target_block_num] = ib2;
                block_label_mask[ib2] = 0;
                target_block_num += 1;
            }
        }
    }
    
    for(ib1=0; ib1<field_info->block_num[field_label_0]; ib1++)
    {   
        // st1 = clock();
        for(i=field_info->block_st_z1[field_label_0][ib1]; i<field_info->block_ed_z1[field_label_0][ib1]; i++)
        {
            // loop the grid in the first zbin, zbin_label_0
            m = i*field_info->field_data_col;

            ra_z1 = field_info->field_data_z1[field_label_0][m+field_info->ra_idx];
            dec_z1 = field_info->field_data_z1[field_label_0][m+field_info->dec_idx];
            cos_dec_z1 = field_info->field_data_z1[field_label_0][m+field_info->cos_dec_idx];

            mg1_z1 = field_info->field_data_z1[field_label_0][m+field_info->mg1_idx];
            mg2_z1 = field_info->field_data_z1[field_label_0][m+field_info->mg2_idx];

            mnu1_z1 = field_info->field_data_z1[field_label_0][m+field_info->mn_idx] +
                        field_info->field_data_z1[field_label_0][m+field_info->mu_idx];
            mnu2_z1 = field_info->field_data_z1[field_label_0][m+field_info->mn_idx] -
                        field_info->field_data_z1[field_label_0][m+field_info->mu_idx];

            expo_label_0 = field_info->field_expo_label_z1[field_label_0][i];
            
            iexpo_len = expo_label_0*field_info->iexpo_chi_block_len;

            for(k=0; k<target_block_num; k++)
            {
                ib2=target_block_label[k];

                for(j=field_info->block_st_z2[field_label_1][ib2]; j<field_info->block_ed_z2[field_label_1][ib2]; j++)
                {   
                    // loop the grid in the second zbin, zbin_label_1
                    n = j*field_info->field_data_col;

                    ra_z2 = field_info->field_data_z2[field_label_1][n+field_info->ra_idx];
                    dec_z2 = field_info->field_data_z2[field_label_1][n+field_info->dec_idx];
                    cos_dec_z2 = field_info->field_data_z2[field_label_1][n+field_info->cos_dec_idx];

                    // the seperation angle (arc minute)
                    delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
                    delta_dec = dec_z2 - dec_z1;
                    delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);
                    
                    theta_tag = -1;
                    for(ir=0; ir<field_info->theta_bin_num; ir++)
                    {
                        if(delta_radius > field_info->theta_bin[ir] and delta_radius <= field_info->theta_bin[ir+1])
                        {theta_tag=ir;break;}
                    }
                    // std::cout<<delta_radius<<" "<<field_info->theta_bin[theta_tag]<<" "<<field_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
                    if(theta_tag > -1)
                    {   
                        // pairs+= 1;

                        // shear estimators rotation (position angle defined as East of North)
                        sin_theta = delta_ra/delta_radius;
                        cos_theta = delta_dec/delta_radius;

                        sin_2theta = 2*sin_theta*cos_theta;
                        cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

                        sin_4theta = 2*sin_2theta*cos_2theta;
                        cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;


                        mg1_z2 = field_info->field_data_z2[field_label_1][n+field_info->mg1_idx]*cos_2theta - 
                                field_info->field_data_z2[field_label_1][n+field_info->mg2_idx]*sin_2theta;
                        mg2_z2 = field_info->field_data_z2[field_label_1][n+field_info->mg1_idx]*sin_2theta + 
                                field_info->field_data_z2[field_label_1][n+field_info->mg2_idx]*cos_2theta;

                        mnu1_z2 = field_info->field_data_z2[field_label_1][n+field_info->mu_idx]*cos_4theta -
                                field_info->field_data_z2[field_label_1][n+field_info->mv_idx]*sin_4theta;
                        mnu2_z2 = mnu1_z2;

                        mnu1_z2 = field_info->field_data_z2[field_label_1][n+field_info->mn_idx] + mnu2_z2;
                        mnu2_z2 = field_info->field_data_z2[field_label_1][n+field_info->mn_idx] - mnu2_z2;
                        
                        // the key part of PDF_SYM
                        ir_len = theta_tag*field_info->ir_chi_block_len + iexpo_len;

                        for(ic=0; ic<field_info->chi_guess_num; ic++)
                        {   
                            if(loop_label >= field_info->gg_len){loop_label = 0;}

                            ic_len = ic*field_info->chi_block_len + ir_len;
                            gg_1 = field_info->gg_1[ic][loop_label];
                            gg_2 = field_info->gg_2[ic][loop_label];
                            // std::cout<<theta_tag<<" "<<expo_label_0<<" "<<ic_len<<std::endl;
                            // std::cout<<loop_label<<std::endl;
                            temp_x = mg1_z1 - gg_1*mnu1_z1;
                            temp_y = mg1_z2 - gg_2*mnu1_z2;
                            
                            // hist_2d_fast(temp_x, temp_y, field_info->mg_bin, field_info->mg_bin_num, field_info->mg_bin_num2,ix, iy);
                            hist_2d_new(temp_x, temp_y, field_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix, iy);
                            // std::cout<<iy<<" "<<ix<<std::endl;

                            field_info->num_count_chit[field_label_0][ic_len + iy*mg_bin_num+ix] += 1;
                            

                            temp_x = mg2_z1 - gg_1*mnu2_z1;
                            temp_y = mg2_z2 - gg_2*mnu2_z2;

                            // hist_2d_fast(temp_x, temp_y, field_info->mg_bin, field_info->mg_bin_num, field_info->mg_bin_num2,ix, iy);
                            hist_2d_new(temp_x, temp_y, field_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix, iy);

                            // std::cout<<iy<<" "<<ix<<std::endl;

                            field_info->num_count_chix[field_label_0][ic_len + iy*field_info->mg_bin_num+ix] += 1;
                            // std::cout<<iy<<" "<<ix<<std::endl;
                            loop_label += 1;
                        }
                    }
                }
            }
        }
        // st2 = clock();
        // std::cout<<"Block "<<ib1<<" "<<pairs<<" pairs "<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
        // pairs = 0;
    }
    field_info->loop_label = loop_label;

    delete[] block_label_mask;
    delete[] target_block_label;
}

void find_pairs_same_field(data_info *field_info, int field_label)
{
    int i, j, m, n, k;

    int ib1, ib2;
    int expo_label_0;
    int ir, theta_tag, ic;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;

    MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
    MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;
    MY_FLOAT temp_x, temp_y;
    int ix, iy, im, im1,im2;
    int iexpo_len, ir_len, ic_len;

    MY_FLOAT gg_1, gg_2;
    
    MY_FLOAT delta_ra, delta_dec, delta_radius;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    double st1, st2;

    double pairs = 0;
    
    int loop_label;
    loop_label = field_info->loop_label;
    int mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3;
    mg_bin_num = field_info->mg_bin_num;
    mg_bin_num1 = field_info->mg_bin_num1;
    mg_bin_num2 = field_info->mg_bin_num2;
    mg_bin_num3 = field_info->mg_bin_num3;

    
    for(ib1=0; ib1<field_info->block_num[field_label]; ib1++)
    {   
        // st1 = clock();
        for(i=field_info->block_st_z1[field_label][ib1]; i<field_info->block_ed_z1[field_label][ib1]; i++)
        {
            // loop the grid in the first zbin, zbin_label_0
            m = i*field_info->field_data_col;

            ra_z1 = field_info->field_data_z1[field_label][m+field_info->ra_idx];
            dec_z1 = field_info->field_data_z1[field_label][m+field_info->dec_idx];
            cos_dec_z1 = field_info->field_data_z1[field_label][m+field_info->cos_dec_idx];

            mg1_z1 = field_info->field_data_z1[field_label][m+field_info->mg1_idx];
            mg2_z1 = field_info->field_data_z1[field_label][m+field_info->mg2_idx];

            mnu1_z1 = field_info->field_data_z1[field_label][m+field_info->mn_idx] +
                        field_info->field_data_z1[field_label][m+field_info->mu_idx];
            mnu2_z1 = field_info->field_data_z1[field_label][m+field_info->mn_idx] -
                        field_info->field_data_z1[field_label][m+field_info->mu_idx];

            expo_label_0 = field_info->field_expo_label_z1[field_label][i];
            
            iexpo_len = expo_label_0*field_info->iexpo_chi_block_len;

            if(field_info->zbin_label_0 != field_info->zbin_label_1){k=0;}
            else{k=ib1;}

            for(ib2=k; ib2<field_info->block_num[field_label]; ib2++)
            {
                for(j=field_info->block_st_z2[field_label][ib2]; j<field_info->block_ed_z2[field_label][ib2]; j++)
                {   
                    if(field_info->field_expo_label_z2[field_label][j] !=  expo_label_0)
                    {
                        // loop the grid in the second zbin, zbin_label_1
                        n = j*field_info->field_data_col;

                        ra_z2 = field_info->field_data_z2[field_label][n+field_info->ra_idx];
                        dec_z2 = field_info->field_data_z2[field_label][n+field_info->dec_idx];
                        cos_dec_z2 = field_info->field_data_z2[field_label][n+field_info->cos_dec_idx];

                        // the seperation angle (arc minute)
                        delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
                        delta_dec = dec_z2 - dec_z1;
                        delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);
                        
                        theta_tag = -1;
                        for(ir=0; ir<field_info->theta_bin_num; ir++)
                        {
                            if(delta_radius > field_info->theta_bin[ir] and delta_radius <= field_info->theta_bin[ir+1])
                            {theta_tag=ir;break;}
                        }
                        // std::cout<<delta_radius<<" "<<field_info->theta_bin[theta_tag]<<" "<<field_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
                        if(theta_tag > -1)
                        {   
                            // pairs += 1;
                            // shear estimators rotation (position angle defined as East of North)
                            sin_theta = delta_ra/delta_radius;
                            cos_theta = delta_dec/delta_radius;

                            sin_2theta = 2*sin_theta*cos_theta;
                            cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

                            sin_4theta = 2*sin_2theta*cos_2theta;
                            cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;


                            mg1_z2 = field_info->field_data_z2[field_label][n+field_info->mg1_idx]*cos_2theta - 
                                    field_info->field_data_z2[field_label][n+field_info->mg2_idx]*sin_2theta;
                            mg2_z2 = field_info->field_data_z2[field_label][n+field_info->mg1_idx]*sin_2theta + 
                                    field_info->field_data_z2[field_label][n+field_info->mg2_idx]*cos_2theta;

                            mnu1_z2 = field_info->field_data_z2[field_label][n+field_info->mu_idx]*cos_4theta -
                                    field_info->field_data_z2[field_label][n+field_info->mv_idx]*sin_4theta;
                            mnu2_z2 = mnu1_z2;

                            mnu1_z2 += field_info->field_data_z2[field_label][n+field_info->mn_idx];
                            mnu2_z2 = field_info->field_data_z2[field_label][n+field_info->mn_idx] - mnu2_z2;
                            
                            // the key part of PDF_SYM
                            ir_len = theta_tag*field_info->ir_chi_block_len + iexpo_len;

                            for(ic=0; ic<field_info->chi_guess_num; ic++)
                            {   
                                if(loop_label >= field_info->gg_len){loop_label = 0;}

                                ic_len = ic*field_info->chi_block_len + ir_len;
                                gg_1 = field_info->gg_1[ic][loop_label];
                                gg_2 = field_info->gg_2[ic][loop_label];

                                temp_x = mg1_z1 - gg_1*mnu1_z1;
                                temp_y = mg1_z2 - gg_2*mnu1_z2;
                                
                            //     // hist_2d_fast(temp_x, temp_y, field_info->mg_bin, field_info->mg_bin_num, field_info->mg_bin_num2,ix, iy);
                                hist_2d_new(temp_x, temp_y, field_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix, iy);
                            //     // std::cout<<iy<<" "<<ix<<std::endl;

                                field_info->num_count_chit[field_label][ic_len + iy*mg_bin_num+ix] += 1;
                                

                                temp_x = mg2_z1 - gg_1*mnu2_z1;
                                temp_y = mg2_z2 - gg_2*mnu2_z2;

                            //     // hist_2d_fast(temp_x, temp_y, field_info->mg_bin, field_info->mg_bin_num, field_info->mg_bin_num2,ix, iy);
                                hist_2d_new(temp_x, temp_y, field_info->mg_bin, mg_bin_num,mg_bin_num1, mg_bin_num2, mg_bin_num3, ix, iy);
   
                            //     // std::cout<<iy<<" "<<ix<<std::endl;

                                field_info->num_count_chix[field_label][ic_len + iy*field_info->mg_bin_num+ix] += 1;
                            //     // std::cout<<iy<<" "<<ix<<std::endl;
                                loop_label += 1;
                            }
                        }
                    }
                    else
                    {
                        continue;
                    }
                }
            }
        }
        // st2 = clock();
        // std::cout<<"Block "<<ib1<<", "<<pairs<<" pairs, "<<(st2-st1)/CLOCKS_PER_SEC<<" Sec"<<std::endl;
        // pairs = 0;
    }
    field_info->loop_label = loop_label;
}

void find_pairs(data_info *field_info, int field_label_0, int field_label_1)
{
    int i, j, m, n, k;

    int ib1, ib2;
    int expo_label_0;
    int ir, theta_tag, ic;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1, delta_len_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2, delta_len_z2;

    MY_FLOAT mg1_z1, mg2_z1, mnu1_z1, mnu2_z1;
    MY_FLOAT mg1_z2, mg2_z2, mnu1_z2, mnu2_z2;
    MY_FLOAT temp_x, temp_y;
    int ix, iy, im, im1,im2;
    int iexpo_len, ir_len, ic_len;

    MY_FLOAT delta_ra, delta_dec, delta_radius;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    int *target_block_label = new int[field_info->block_num[field_label_1]]{};
    int *block_label_mask = new int[field_info->block_num[field_label_1]]{};
    int target_block_num = 0;
    double st1, st2;

    int loop_label;
    loop_label = field_info->loop_label;

    // decide which the blocks in the target field, field_label_1, to calculation
    // it will calculate the distance of all the block pairs, if the distance is 
    // smaller than the biggest separation angle, this block in "field_label_1"
    // will be labeled as tagert block, it will be included in the calculation.
    // it will loop all the galaxy pairs which is labeled by the
    // block, then this pair will distributed to the right separation angle bin.
    for(ib2=0; ib2<field_info->block_num[field_label_1]; ib2++)
    {
        block_label_mask[ib2] = 1;
    }
    for(ib1=0; ib1<field_info->block_num[field_label_0]; ib1++)
    {
        ra_z1 = field_info->block_cen_ra[field_label_0][ib1];
        dec_z1 = field_info->block_cen_dec[field_label_0][ib1];
        cos_dec_z1 = field_info->block_cen_cos_dec[field_label_0][ib1];
        delta_len_z1 = field_info->block_delta_len[field_label_0][ib1];

        // if search pairs between different zbin, if should loop all the grids
        // else it should search in the grid that has a label larger than "ib1"
        // to avoid double counting
        if(field_label_0 == field_label_1 and field_info->zbin_label_0 == field_info->zbin_label_1){k=ib1;}
        else{k=0;}

        for(ib2=k; ib2<field_info->block_num[field_label_1]; ib2++)
        {
            ra_z2 = field_info->block_cen_ra[field_label_1][ib2];
            dec_z2 = field_info->block_cen_dec[field_label_1][ib2];
            cos_dec_z2 = field_info->block_cen_cos_dec[field_label_1][ib2];
            delta_len_z2 = field_info->block_delta_len[field_label_1][ib2];

            // the seperation angle (arc minute)
            delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
            delta_dec = dec_z2 - dec_z1;
            delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec) - delta_len_z1 - delta_len_z2;

            if(delta_radius <= field_info->theta_bin[field_info->theta_bin_num] and block_label_mask[ib2] == 1)
            {
                target_block_label[target_block_num] = ib2;
                block_label_mask[ib2] = 0;
                target_block_num += 1;
            }
        }
    }
    
    for(ib1=0; ib1<field_info->block_num[field_label_0]; ib1++)
    {   
        // st1 = clock();
        for(i=field_info->block_st_z1[field_label_0][ib1]; i<field_info->block_ed_z1[field_label_0][ib1]; i++)
        {
            // loop the grid in the first zbin, zbin_label_0
            m = i*field_info->field_data_col;

            ra_z1 = field_info->field_data_z1[field_label_0][m+field_info->ra_idx];
            dec_z1 = field_info->field_data_z1[field_label_0][m+field_info->dec_idx];
            cos_dec_z1 = field_info->field_data_z1[field_label_0][m+field_info->cos_dec_idx];

            mg1_z1 = field_info->field_data_z1[field_label_0][m+field_info->mg1_idx];
            mg2_z1 = field_info->field_data_z1[field_label_0][m+field_info->mg2_idx];

            mnu1_z1 = field_info->field_data_z1[field_label_0][m+field_info->mn_idx] +
                        field_info->field_data_z1[field_label_0][m+field_info->mu_idx];
            mnu2_z1 = field_info->field_data_z1[field_label_0][m+field_info->mu_idx] -
                        field_info->field_data_z1[field_label_0][m+field_info->mu_idx];

            expo_label_0 = field_info->field_expo_label_z1[field_label_0][i];
            
            iexpo_len = expo_label_0*field_info->iexpo_chi_block_len;


            for(k=0; k<target_block_num; k++)
            {
                ib2=target_block_label[k];

                for(j=field_info->block_st_z2[field_label_1][ib2]; j<field_info->block_ed_z2[field_label_1][ib2]; j++)
                {   

                    // loop the grid in the second zbin, zbin_label_1
                    n = j*field_info->field_data_col;

                    ra_z2 = field_info->field_data_z2[field_label_1][n+field_info->ra_idx];
                    dec_z2 = field_info->field_data_z2[field_label_1][n+field_info->dec_idx];
                    cos_dec_z2 = field_info->field_data_z2[field_label_1][n+field_info->cos_dec_idx];

                    // the seperation angle (arc minute)
                    delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
                    delta_dec = dec_z2 - dec_z1;
                    delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);
                    
                    theta_tag = -1;
                    for(ir=0; ir<field_info->theta_bin_num; ir++)
                    {
                        if(delta_radius > field_info->theta_bin[ir] and delta_radius <= field_info->theta_bin[ir+1])
                        {theta_tag=ir;break;}
                    }
                    // std::cout<<delta_radius<<" "<<field_info->theta_bin[theta_tag]<<" "<<field_info->theta_bin[theta_tag+1]<<" "<<theta_tag<<std::endl;
                    if(theta_tag > -1)
                    {
                        // shear estimators rotation (position angle defined as East of North)
                        sin_theta = delta_ra/delta_radius;
                        cos_theta = delta_dec/delta_radius;

                        sin_2theta = 2*sin_theta*cos_theta;
                        cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

                        sin_4theta = 2*sin_2theta*cos_2theta;
                        cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;


                        mg1_z2 = field_info->field_data_z2[field_label_1][n+field_info->mg1_idx]*cos_2theta - 
                                field_info->field_data_z2[field_label_1][n+field_info->mg2_idx]*sin_2theta;
                        mg2_z2 = field_info->field_data_z2[field_label_1][n+field_info->mg1_idx]*sin_2theta + 
                                field_info->field_data_z2[field_label_1][n+field_info->mg2_idx]*cos_2theta;

                        mnu1_z2 = field_info->field_data_z2[field_label_1][n+field_info->mu_idx]*cos_4theta -
                                field_info->field_data_z2[field_label_1][n+field_info->mv_idx]*sin_4theta;
                        mnu2_z2 = mnu1_z2;

                        mnu1_z2 += field_info->field_data_z2[field_label_1][n+field_info->mn_idx];
                        mnu2_z2 = field_info->field_data_z2[field_label_1][n+field_info->mn_idx] - mnu2_z2;
                        
                        // the key part of PDF_SYM
                        ir_len = theta_tag*field_info->ir_chi_block_len + iexpo_len;

                        for(ic=0; ic<field_info->chi_guess_num; ic++)
                        {   
                            if(loop_label >= field_info->gg_len){loop_label = 0;}

                            ic_len = ic*field_info->chi_block_len + ir_len;
                            // std::cout<<theta_tag<<" "<<expo_label_0<<" "<<ic_len<<std::endl;
                            // std::cout<<loop_label<<std::endl;
                            temp_x = mg1_z1 - field_info->gg_1[ic][loop_label]*mnu1_z1;
                            temp_y = mg1_z2 - field_info->gg_2[ic][loop_label]*mnu1_z2;
                            
                            // hist_2d_fast(temp_x, temp_y, field_info->mg_bin, field_info->mg_bin_num, field_info->mg_bin_num2,ix, iy);
                            hist_2d_new(temp_x, temp_y, field_info->mg_bin, field_info->mg_bin_num,field_info->mg_bin_num1, 
                            field_info->mg_bin_num2, field_info->mg_bin_num3, ix, iy);

                            field_info->num_count_chit[field_label_0][ic_len + iy*field_info->mg_bin_num+ix] += 1;
                            

                            temp_x = mg2_z1 - field_info->gg_1[ic][loop_label]*mnu2_z1;
                            temp_y = mg2_z2 - field_info->gg_2[ic][loop_label]*mnu2_z2;

                            // hist_2d_fast(temp_x, temp_y, field_info->mg_bin, field_info->mg_bin_num, field_info->mg_bin_num2,ix, iy);
                            hist_2d_new(temp_x, temp_y, field_info->mg_bin, field_info->mg_bin_num,field_info->mg_bin_num1, 
                            field_info->mg_bin_num2, field_info->mg_bin_num3, ix, iy);
                           
                            // std::cout<<iy<<" "<<ix<<std::endl;

                            field_info->num_count_chix[field_label_0][ic_len + iy*field_info->mg_bin_num+ix] += 1;
                            // std::cout<<iy<<" "<<ix<<std::endl;
                            loop_label += 1;
                        }
                    }
                }
            }
        }
        // st2 = clock();
        // std::cout<<"One block "<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
    }
    field_info->loop_label = loop_label;

    delete[] block_label_mask;
    delete[] target_block_label;
}


void save_field_chi_block(data_info*field_info, int field_label)
{
    int i, j, k, m, n, p;
    char result_path[600], set_name[50];
    double *radius_chi_buf;
    double temp;

    radius_chi_buf = new double[field_info->ir_chi_block_len];

    p = field_info->chi_guess_num*field_info->mg_bin_num;
    
    for(i=0; i<field_info->exposure_num_of_field[field_label]; i++)
    {   
        k = 0;
        temp = 0;
        // std::cout<<i<<" "<<temp<<" "<<field_info->exposure_num_of_field[field_label]<<std::endl;
        for(j=0; j<field_info->iexpo_chi_block_len; j++)
        {
            temp += field_info->num_count_chit[field_label][i*field_info->iexpo_chi_block_len+j];
        }
        // std::cout<<i<<" "<<temp<<std::endl;
        // if there're number counts in chi block of this exposure, it will be saved to file
        if(fabs(temp) > 1)
        {
            sprintf(result_path, "%s/result/%s-%d.hdf5",field_info->parent_path, field_info->field_name[field_label], k);

            for(j=0; j<field_info->theta_bin_num; j++)
            {
                n = i*field_info->iexpo_chi_block_len + j*field_info->ir_chi_block_len;
                for(m=0; m<field_info->ir_chi_block_len; m++)
                {
                    radius_chi_buf[m] = field_info->num_count_chit[field_label][n+m];
                }

                sprintf(set_name,"/theta_%d", j);
                if(file_exist(result_path))
                {
                    write_h5(result_path, set_name, radius_chi_buf, p, field_info->mg_bin_num, false);
                }
                else
                {
                    write_h5(result_path, set_name, radius_chi_buf, p, field_info->mg_bin_num, true);
                }
            }
            k++;
        }

    }
    delete[] radius_chi_buf;
}

void find_pairs_same_field_(data_info *field_info, int field_label)
{
    int i, j, m, n, k;

    int ib1, ib2;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2;

    MY_FLOAT mg1_z1, mg2_z1, mn_z1, mu_z1;
    MY_FLOAT mg1_z2, mg2_z2, mn_z2, mu_z2;

    MY_FLOAT delta_ra, delta_dec, delta_radius;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    if(field_info->zbin_label_0 == field_info->zbin_label_1)
    {   
        // if search pairs between different zbin, if should loop all the grids
        // else it should search in the grid that has a label larger than "ib1"
        // to avoid double count
        for(i=0; i<field_info->total_gal_num_z1[field_label]; i++)
        {
            // loop the grid in the first zbin, zbin_label_0
            m = i*field_info->field_data_col;

            ra_z1 = field_info->field_data_z1[field_label][m+field_info->ra_idx];
            dec_z1 = field_info->field_data_z1[field_label][m+field_info->dec_idx];
            cos_dec_z1 = field_info->field_data_z1[field_label][m+field_info->cos_dec_idx];

            mg1_z1 = field_info->field_data_z1[field_label][m+field_info->mg1_idx];
            mg2_z1 = field_info->field_data_z1[field_label][m+field_info->mg2_idx];

            mn_z1 = field_info->field_data_z1[field_label][m+field_info->mn_idx];
            mu_z1 = field_info->field_data_z1[field_label][m+field_info->mu_idx];

            for(j=0; j<field_info->total_gal_num_z2[field_label]; j++)
            {   
                // search pairs in different exposures
                if(field_info->field_block_label_z1[i] > field_info->field_block_label_z2[j] or 
                    field_info->field_expo_label_z1[i] == field_info->field_expo_label_z2[j])
                {
                    continue;
                }
                // loop the grid in the second zbin, zbin_label_1
                n = j*field_info->field_data_col;

                ra_z2 = field_info->field_data_z2[field_label][n+field_info->ra_idx];
                dec_z2 = field_info->field_data_z2[field_label][n+field_info->dec_idx];
                cos_dec_z2 = field_info->field_data_z2[field_label][n+field_info->cos_dec_idx];

                // the seperation angle (arc minute)
                delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
                delta_dec = dec_z2 - dec_z1;
                delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);

                // shear estimators rotation (position angle defined as East of North)
                sin_theta = delta_ra/delta_radius;
                cos_theta = delta_dec/delta_radius;

                sin_2theta = 2*sin_theta*cos_theta;
                cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

                sin_4theta = 2*sin_2theta*cos_2theta;
                cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;


                mg1_z2 = field_info->field_data_z2[field_label][n+field_info->mg1_idx];
                mg2_z2 = field_info->field_data_z2[field_label][n+field_info->mg2_idx];

                mn_z2 = field_info->field_data_z2[field_label][n+field_info->mn_idx];
                mu_z2 = field_info->field_data_z2[field_label][n+field_info->mu_idx];

            }
        }

    }
}