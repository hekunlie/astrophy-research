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
        strs >> field_info->exposure_name[line_count];

        strs >> field_info->field_label[line_count];

        strs >> field_info->exposure_label[line_count];
        strs >> field_info->exposure_num_of_field[line_count];

        strs >> field_info->field_cen_ra[line_count];
        strs >> field_info->field_cen_dec[line_count];

        strs >> field_info->delta_ra[line_count];
        strs >> field_info->delta_dec[line_count];

        // std::cout << str << std::endl;
        line_count += 1;
    }
    read_file_num = line_count;

}


void initialize(char *file_path, data_info *field_info, int total_exposure_num, int numprocs, int rank)
{
    int i;

    for(i=0;i<total_exposure_num;i++)
    {
        field_info->field_name[i] = new char[100];
        field_info->exposure_name[i] = new char[100];  
    }

    field_info->field_label = new int[total_exposure_num];
    field_info->exposure_label = new int[total_exposure_num];
    field_info->exposure_num_of_field = new int[total_exposure_num];
     
    field_info->field_cen_ra = new float[total_exposure_num];  
    field_info->field_cen_dec = new float[total_exposure_num]; 

    field_info->delta_ra = new float[total_exposure_num];  
    field_info->delta_dec = new float[total_exposure_num];

    field_info->total_exposure_num = total_exposure_num;

    
    read_file(file_path, field_info, i);

    // tell the cpus from where their tasks start
    task_distribution(total_exposure_num, numprocs, rank, field_info);

}


void task_distribution(int total_task_num, int portion, int my_id, data_info *field_info)
{   
    int i, j, m, n;
    m = total_task_num/portion;
    n = total_task_num%portion;

    field_info->my_exposure_num = m;
    if(n > 0 and my_id < n){field_info->my_exposure_num += 1;}

    j = 0;
    for(i=0; i<my_id; i++)
    {
        if(i<n){j += m+1;}
        else{j += m;}
    }
    field_info->my_exposure_st = j;
    field_info->my_exposure_ed = j + field_info->my_exposure_num;

}

void find_pairs(data_info *field_info)
{
    int i, j, k;
    int expos_st, expos_ed, expos_num;

    double radius;
    double expos_ra_1, expos_dec_1, expos_dra_1, expos_ddec_1;
    double expos_ra_2, expos_dec_2, expos_dra_2, expos_ddec_2;

    expos_st = field_info->my_exposure_st;
    expos_ed = field_info->my_exposure_ed;
    expos_num = field_info->my_exposure_num;

    for(i=expos_st; i<expos_ed; i++)
    {
        
    }

}