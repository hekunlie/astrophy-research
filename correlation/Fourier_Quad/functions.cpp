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

void read_field_data(data_info *field_info, int zbin_label_0, int zbin_label_1)
{
    int i, j;
    char set_name[100];
    for(i=0; i<field_info->total_field_num; i++)
    {   
        // std::cout<<"Read "<<field_info->field_name[i]<<std::endl;
        // read the gal num in each zbin stored in each field file
        sprintf(set_name, "/total_num_in_zbin");
        read_h5(field_info->field_name[i], set_name, field_info->num_in_zbin);
        field_info->total_gal_num_z1[i] = field_info->num_in_zbin[zbin_label_0];
        field_info->total_gal_num_z2[i] = field_info->num_in_zbin[zbin_label_1];

        // read the array length to get the data column
        if(i == 0)
        {   
            // array length, all elements, of zbin i
            sprintf(set_name, "/z%d/field",zbin_label_0);
            read_h5_datasize(field_info->field_name[i], set_name, j);
            // data col = data length/gal_num, it's the same for all fields
            field_info->field_data_col = j / field_info->num_in_zbin[zbin_label_0];
        }


        ///////////////// read all the fields  //////////////////////////
        // zbin i
        j = field_info->num_in_zbin[zbin_label_0]*field_info->field_data_col;
        field_info->field_data_z1[i] = new MY_FLOAT[j]{};
        sprintf(set_name, "/z%d/field",zbin_label_0);
        read_h5(field_info->field_name[i], set_name, field_info->field_data_z1[i]);
        // zbin j
        j = field_info->num_in_zbin[zbin_label_1]*field_info->field_data_col;
        field_info->field_data_z2[i] = new MY_FLOAT[j]{};
        sprintf(set_name, "/z%d/field",zbin_label_1);
        read_h5(field_info->field_name[i], set_name, field_info->field_data_z2[i]);


        //////////// read the block inform in each field  //////////////////////////
        if(i==0){ field_info->block_num = new int[field_info->total_field_num]{}; }
        // number of block in each field
        sprintf(set_name, "/block_cen_ra");
        read_h5_datasize(field_info->field_name[i], set_name, j);
        field_info->block_num[i] = j;

        // RA of the center of each block
        field_info->block_cen_ra[i] = new MY_FLOAT[j];
        sprintf(set_name, "/block_cen_ra");
        read_h5(field_info->field_name[i], set_name, field_info->block_cen_ra[i]);

        // Dec of the center of each block
        field_info->block_cen_dec[i] = new MY_FLOAT[j];
        sprintf(set_name, "/block_cen_dec");
        read_h5(field_info->field_name[i], set_name, field_info->block_cen_dec[i]);
        
        // cos(Dec) of the center of each block
        field_info->block_cen_cos_dec[i] = new MY_FLOAT[j];
        sprintf(set_name, "/block_cen_cos_dec");
        read_h5(field_info->field_name[i], set_name, field_info->block_cen_cos_dec[i]);

        // block start & end 
        field_info->block_st_z1[i] = new MY_FLOAT[j];
        field_info->block_st_z2[i] = new MY_FLOAT[j];
        field_info->block_ed_z1[i] = new MY_FLOAT[j];
        field_info->block_ed_z2[i] = new MY_FLOAT[j];
        sprintf(set_name, "/z%d/block_st", zbin_label_0);
        read_h5(field_info->field_name[i], set_name, field_info->block_st_z1[i]);
        sprintf(set_name, "/z%d/block_st", zbin_label_1);
        read_h5(field_info->field_name[i], set_name, field_info->block_st_z2[i]);
        sprintf(set_name, "/z%d/block_ed", zbin_label_0);
        read_h5(field_info->field_name[i], set_name, field_info->block_ed_z1[i]);
        sprintf(set_name, "/z%d/block_ed", zbin_label_1);
        read_h5(field_info->field_name[i], set_name, field_info->block_ed_z2[i]);
    }
}


void initialize(char *file_path, data_info *field_info, int total_field_num, int numprocs, int rank)
{
    int i;
    char set_name[100];

    field_info->total_field_num = total_field_num;

    for(i=0;i<total_field_num;i++)
    {   
        // the field file directories
        field_info->field_name[i] = new char[400];  
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
    read_h5_datasize(field_info->field_name[0], set_name,i);
    field_info->theta_bin = new MY_FLOAT[i]{};
    field_info->theta_bin_num = i -1;
    read_h5(field_info->field_name[0], set_name, field_info->theta_bin);

    // read redshift bin
    sprintf(set_name,"/redshift_bin");
    read_h5_datasize(field_info->field_name[0], set_name,i);
    field_info->zbin = new MY_FLOAT[i]{};
    field_info->zbin_num = i -1;
    field_info->num_in_zbin = new int[i]{};
    read_h5(field_info->field_name[0], set_name, field_info->zbin);

    // read the inform of the PDF_SYM
    sprintf(set_name, "/chi_guess");
    read_h5_datasize(field_info->field_name[0], set_name, i);
    field_info->chi_guess_num = i;
    field_info->chi_guess = new MY_FLOAT[field_info->chi_guess_num];
    read_h5(field_info->field_name[0], set_name, field_info->chi_guess);
    // the num count of each bin (for the PDF_SYM)
    // for each theta bin, there're chi_guess_num * chi_bin_num elements
    for(i=0; i<field_info->theta_bin_num; i++)
    {
        field_info->num_count[i] = new double[field_info->chi_guess_num*field_info->chi_bin_num]{};
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

void find_pairs(data_info *field_info)
{
    int i, j, k;
    int expos_st, expos_ed, expos_num;

    double radius;
    double expos_ra_1, expos_dec_1, expos_dra_1, expos_ddec_1;
    double expos_ra_2, expos_dec_2, expos_dra_2, expos_ddec_2;


    for(i=expos_st; i<expos_ed; i++)
    {
        ;
    }

}