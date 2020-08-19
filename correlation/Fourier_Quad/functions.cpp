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
        sprintf(set_name, "/block_cen_ra");
        read_h5(field_info->field_name_path[i], set_name, field_info->block_cen_ra[i]);

        // Dec of the center of each block
        field_info->block_cen_dec[i] = new MY_FLOAT[j];
        sprintf(set_name, "/block_cen_dec");
        read_h5(field_info->field_name_path[i], set_name, field_info->block_cen_dec[i]);
        
        // cos(Dec) of the center of each block
        field_info->block_cen_cos_dec[i] = new MY_FLOAT[j];
        sprintf(set_name, "/block_cen_cos_dec");
        read_h5(field_info->field_name_path[i], set_name, field_info->block_cen_cos_dec[i]);

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


void initialize(char *file_path, data_info *field_info, int total_field_num, int numprocs, int rank)
{
    int i;
    char set_name[100];
    char data_path[500];

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


    // read radius bin
    sprintf(set_name,"/theta_bin");
    read_h5_datasize(field_info->field_name_path[0], set_name,i);
    field_info->theta_bin = new MY_FLOAT[i]{};
    field_info->theta_bin_num = i -1;
    read_h5(field_info->field_name_path[0], set_name, field_info->theta_bin);

    // read redshift bin
    sprintf(set_name,"/redshift_bin");
    read_h5_datasize(field_info->field_name_path[0], set_name,i);
    field_info->zbin = new MY_FLOAT[i]{};
    field_info->zbin_num = i -1;
    field_info->num_in_zbin = new int[i]{};
    read_h5(field_info->field_name_path[0], set_name, field_info->zbin);

    // read the inform of the PDF_SYM
    sprintf(set_name, "/chi_guess");
    read_h5_datasize(field_info->field_name_path[0], set_name, i);
    field_info->chi_guess_num = i;
    field_info->chi_guess = new MY_FLOAT[field_info->chi_guess_num];
    read_h5(field_info->field_name_path[0], set_name, field_info->chi_guess);
    // the num count of each bin (for the PDF_SYM)
    // for each theta bin, there're chi_guess_num * chi_bin_num* chi_bin_num elements
    // tangential and cross components
    for(i=0; i<field_info->theta_bin_num; i++)
    {
        field_info->num_count_chit[i] = new double[field_info->chi_guess_num*field_info->chi_bin_num*field_info->chi_bin_num]{};
        field_info->num_count_chix[i] = new double[field_info->chi_guess_num*field_info->chi_bin_num*field_info->chi_bin_num]{};
    }
    for(i=0; i<field_info->chi_guess_num; i++)
    {
        sprintf(set_name,"/gg_t");
        sprintf(data_path,"%s/cata/cor_%d.hdf5", i);
        read_h5_datasize(data_path, set_name, field_info->gg_len);

        field_info->gg_tt[i] = new MY_FLOAT[field_info->gg_len]{};
        field_info->gg_tt[i] = new MY_FLOAT[field_info->gg_len]{};
    }

    // task distribution
    field_info->field_num_each_rank = new int[total_field_num]{};
    task_distribution(numprocs, rank, field_info);

}


void task_distribution(int portion, int my_id, data_info *field_info)
{   
    int i, j, m, n;
    m = field_info->total_field_num/portion;
    n = field_info->total_field_num%portion;

    for(i=0; i<field_info->total_field_num; i++)
    {
        field_info->field_num_each_rank[i] = m;
        if(i<n){field_info->field_num_each_rank[i]+=1;}
    }
    

    j = 0;
    for(i=0; i<my_id; i++)
    {
        j += field_info->field_num_each_rank[i];
    }
    field_info->my_field_st = j;
    field_info->my_field_ed = j + field_info->field_num_each_rank[my_id];

}



void fast_hist(MY_FLOAT data, MY_FLOAT*bins, int *num_in_bin, int bin_num)
{
    int i;
    for(i=0; i<bin_num; i++)
    {
        if(data > bins[i] and data <= bins[i+1]){num_in_bin[i] += 1;break;}
    }

}

void find_pairs_diff_field(data_info *field_info, int field_label_0, int field_label_1)
{
    int i, j, k;

    int ib1, ib2;

    MY_FLOAT mg1_z1, mg2_z1, mn_z1, mu_z1;
    MY_FLOAT mg1_z2, mg2_z2, mn_z2, mu_z2;

    MY_FLOAT delta_ra, delta_dec, delta_radius;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    for(ib1=0; ib1<field_info->block_num[field_label_0]; ib1++)
    {   
        for(ib2=0; ib2<field_info->block_num[field_label_1]; ib2++)
        {
            ;
              // for(i=0; i<field_info.block_num[my_fnm]; i++)
        // {   
        //     m = i*field_info.field_data_col;

        //     ra_z1 = field_info.field_data_z1[my_fnm][m+field_info.ra_idx];
        //     dec_z1 = field_info.field_data_z1[my_fnm][m+field_info.dec_idx];
        //     cos_dec_z1 = field_info.field_data_z1[my_fnm][m+field_info.cos_dec_idx];

        //     mg1_z1 = field_info.field_data_z1[my_fnm][n+field_info.mg1_idx];
        //     mg2_z1 = field_info.field_data_z1[my_fnm][n+field_info.mg1_idx];

        //     mn_z1 = field_info.field_data_z1[my_fnm][n+field_info.mn_idx];
        //     mu_z1 = field_info.field_data_z1[my_fnm][n+field_info.mu_idx];

        //     for(j=i; j<field_info.block_num[my_fnm]; j++)
        //     {   
        //         n = j*field_info.field_data_col;
        //         ra_z2 = field_info.field_data_z2[my_fnm][n+field_info.ra_idx];
        //         dec_z2 = field_info.field_data_z2[my_fnm][n+field_info.dec_idx];
        //         cos_dec_z2 = field_info.field_data_z2[my_fnm][n+field_info.cos_dec_idx];
                
        //         // the seperation angle (arc minute)
        //         delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
        //         delta_dec = dec_z2 - dec_z1;
        //         delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);
 
        //         // shear estimators rotation (position angle defined as East of North)
        //         sin_theta = delta_ra/delta_radius;
        //         cos_theta = delta_dec/delta_radius;

        //         sin_2theta = 2*sin_theta*cos_theta;
        //         cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

        //         sin_4theta = 2*sin_2theta*cos_2theta;
        //         cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;

        //         mg1_z2 = field_info.field_data_z2[my_fnm][n+field_info.mg1_idx]*cos_2theta - 
        //             field_info.field_data_z2[my_fnm][n+field_info.mg2_idx]*sin_2theta;
        //         mg2_z2 = field_info.field_data_z2[my_fnm][n+field_info.mg1_idx]*sin_2theta + 
        //             field_info.field_data_z2[my_fnm][n+field_info.mg2_idx]*cos_2theta;

        //         mn_z2 = field_info.field_data_z2[my_fnm][n+field_info.mn_idx];

        //         mu_z2 = field_info.field_data_z2[my_fnm][n+field_info.mu_idx]*cos_4theta - 
        //             field_info.field_data_z2[my_fnm][n+field_info.mv_idx]*sin_4theta;
        //         for(ic=0; ic<chi_num; ic++)
        //         {
                    
        //         }
        //         // loop the radius bin
        //         for(ir=0;ir<radius_bin_num;ir++)
        //         {

        //         }
        //     }
        // }
        }

    //     for(ib1=0; ib1<field_info->block_num[field_label_0]; ib1++)
    // {   
    //     // if search pairs between different zbin, if should loop all the grids
    //     // else it should search in the grid that has a label larger than "ib1"
    //     // to avoid double count
    //     if(field_info->zbin_label_0 != field_info->zbin_label_1){k=0;}
    //     else{k=ib1;}

    //     for(ib2=k; ib2<field_info->block_num[field_label_1]; ib2++)
    //     {   
    //         for(i=field_info->block_st_z1[field_label_0][ib1]; i<field_info->block_ed_z1[field_label_0][ib1]; i++);
    //         {
    //             // loop the grid in the first zbin, zbin_label_0
    //             m = i*field_info->field_data_col;

    //             ra_z1 = field_info->field_data_z1[field_label_0][m+field_info->ra_idx];
    //             dec_z1 = field_info->field_data_z1[field_label_0][m+field_info->dec_idx];
    //             cos_dec_z1 = field_info->field_data_z1[field_label_0][m+field_info->cos_dec_idx];

    //             mg1_z1 = field_info->field_data_z1[field_label_0][m+field_info->mg1_idx];
    //             mg2_z1 = field_info->field_data_z1[field_label_0][m+field_info->mg2_idx];

    //             mn_z1 = field_info->field_data_z1[field_label_0][m+field_info->mn_idx];
    //             mu_z1 = field_info->field_data_z1[field_label_0][m+field_info->mu_idx];

    //             for(j=field_info->block_st_z2[field_label_1][ib2]; j<field_info->block_ed_z2[field_label_1][ib2]; j++)
    //             {   
    //                 // loop the grid in the second zbin, zbin_label_1
    //                 n = j*field_info->field_data_col;

    //                 ra_z2 = field_info->field_data_z2[field_label_1][n+field_info->ra_idx];
    //                 dec_z2 = field_info->field_data_z2[field_label_1][n+field_info->dec_idx];
    //                 cos_dec_z2 = field_info->field_data_z2[field_label_1][n+field_info->cos_dec_idx];

    //                 // the seperation angle (arc minute)
    //                 delta_ra = (ra_z2 - ra_z1)*cos_dec_z1;
    //                 delta_dec = dec_z2 - dec_z1;
    //                 delta_radius = sqrt(delta_ra*delta_ra + delta_dec*delta_dec);

    //                 // shear estimators rotation (position angle defined as East of North)
    //                 sin_theta = delta_ra/delta_radius;
    //                 cos_theta = delta_dec/delta_radius;

    //                 sin_2theta = 2*sin_theta*cos_theta;
    //                 cos_2theta = cos_theta*cos_theta - sin_theta*sin_theta;

    //                 sin_4theta = 2*sin_2theta*cos_2theta;
    //                 cos_4theta = cos_2theta*cos_2theta - sin_2theta*sin_2theta;


    //                 mg1_z2 = field_info->field_data_z2[field_label][n+field_info->mg1_idx];
    //                 mg2_z2 = field_info->field_data_z2[field_label][n+field_info->mg2_idx];

    //                 mn_z2 = field_info->field_data_z2[field_label][n+field_info->mn_idx];
    //                 mu_z2 = field_info->field_data_z2[field_label][n+field_info->mu_idx];

    //             }
    //         }
    //     }
    // }
    }
}

void find_pairs_same_field(data_info *field_info, int field_label)
{
    int i, j, m, n, k;

    int ib1, ib2;
    int expo_label_0, expo_label_1;
    int ir,ic;
    MY_FLOAT ra_z1, dec_z1, cos_dec_z1;
    MY_FLOAT ra_z2, dec_z2, cos_dec_z2;

    MY_FLOAT mg1_z1, mg2_z1, mn_z1, mu_z1;
    MY_FLOAT mg1_z2, mg2_z2, mn_z2, mu_z2;

    MY_FLOAT delta_ra, delta_dec, delta_radius;
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;

    for(ib1=0; ib1<field_info->block_num[field_label]; ib1++)
    {   
        // if search pairs between different zbin, if should loop all the grids
        // else it should search in the grid that has a label larger than "ib1"
        // to avoid double count
        if(field_info->zbin_label_0 != field_info->zbin_label_1){k=0;}
        else{k=ib1;}

        for(i=field_info->block_st_z1[field_label][ib1]; i<field_info->block_ed_z1[field_label][ib1]; i++)
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

            expo_label_0 = field_info->field_expo_label_z1[field_label][i];

            for(ib2=k; ib2<field_info->block_num[field_label]; ib2++)
            {
                for(j=field_info->block_st_z2[field_label][ib2]; j<field_info->block_ed_z2[field_label][ib2]; j++)
                {   
                    if(field_info->field_expo_label_z2[field_label][j] ==  expo_label_0)
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