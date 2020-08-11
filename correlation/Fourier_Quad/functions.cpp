#include"functions.h"

void read_file(char *file_path, data_info *field_info, int &read_file_num)
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

        strs >> field_info->delta_ra[line_count];
        strs >> field_info->delta_dec[line_count];
        strs >> field_info->delta_len[line_count];


        // std::cout << str << std::endl;
        line_count += 1;
    }
    read_file_num = line_count;

}


void initialize(char *file_path, data_info *field_info, int total_field_num, int numprocs, int rank)
{
    int i;

    field_info->total_field_num = total_field_num;

    for(i=0;i<total_field_num;i++)
    {   
        // the field file directories
        field_info->field_name[i] = new char[400];  
    }

    field_info->exposure_num_of_field = new int[total_field_num];
     
    field_info->field_cen_ra = new float[total_field_num];  
    field_info->field_cen_dec = new float[total_field_num]; 

    field_info->delta_ra = new float[total_field_num];  
    field_info->delta_dec = new float[total_field_num];
    field_info->delta_len = new float[total_field_num];

    
    read_file(file_path, field_info, i);

    // tell the cpus from where their tasks start

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