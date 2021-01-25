#include<science_cal_lib.h>

void separation_angle_1(const double RA1, const double DEC1, const double RA2, const double DEC2, double &sep_radian)
{
    double diff_ra_rad = (RA1 - RA2)*DEG2RAD;
    double dec1_rad, dec2_rad;
    dec1_rad = DEC1*DEG2RAD;
    dec2_rad = DEC2*DEG2RAD;

    double m = sin(dec1_rad)*sin(dec2_rad) + cos(dec1_rad)*cos(dec2_rad)*cos(diff_ra_rad);
    sep_radian = acos(m);
}
void separation_angle_1(const float RA1, const float DEC1, const float RA2, const float DEC2, float &sep_radian)
{
    float diff_ra_rad = (RA1 - RA2)*DEG2RAD;
    float dec1_rad, dec2_rad;
    dec1_rad = DEC1*DEG2RAD;
    dec2_rad = DEC2*DEG2RAD;

    float m = sin(dec1_rad)*sin(dec2_rad) + cos(dec1_rad)*cos(dec2_rad)*cos(diff_ra_rad);
    sep_radian = acos(m);
}

void separation_angle_2(const double RA1, const double DEC1, const double RA2, const double DEC2, double &sep_radian)
{
	double dec1_rad, dec2_rad, ra1_rad, ra2_rad;
	double diff_dec_rad, diff_ra_rad;
	double cos_dec1, cos_dec2, sin_dec1, sin_dec2, cos_diff_ra, sin_diff_ra;
	double m, m1,m2, n;

	dec1_rad = DEC1 * DEG2RAD;
	dec2_rad = DEC2 * DEG2RAD;

	diff_ra_rad = (RA1 - RA2)* DEG2RAD;

	cos_dec1 = cos(dec1_rad);
	cos_dec2 = cos(dec2_rad);
	sin_dec1 = sin(dec1_rad);
	sin_dec2 = sin(dec2_rad);

	sin_diff_ra = sin(diff_ra_rad);
	cos_diff_ra = cos(diff_ra_rad);

	m1 = cos_dec2 * sin_diff_ra;
	m2 = cos_dec1 * sin_dec2 - sin_dec1 * cos_dec2*cos_diff_ra;
	m = sqrt(m1*m1 + m2 * m2);

	n = sin_dec1 * sin_dec2 + cos_dec1 * cos_dec2*cos_diff_ra;

	sep_radian = fabs(atan2(m, n));
}

void separation_angle_2(const float RA1, const float DEC1, const float RA2, const float DEC2, float &sep_radian)
{
	float dec1_rad, dec2_rad, ra1_rad, ra2_rad;
	float diff_dec_rad, diff_ra_rad;
	float cos_dec1, cos_dec2, sin_dec1, sin_dec2, cos_diff_ra, sin_diff_ra;
	float m, m1,m2, n;

	dec1_rad = DEC1 * DEG2RAD;
	dec2_rad = DEC2 * DEG2RAD;

	diff_ra_rad = (RA1 - RA2)* DEG2RAD;

	cos_dec1 = cos(dec1_rad);
	cos_dec2 = cos(dec2_rad);
	sin_dec1 = sin(dec1_rad);
	sin_dec2 = sin(dec2_rad);

	sin_diff_ra = sin(diff_ra_rad);
	cos_diff_ra = cos(diff_ra_rad);

	m1 = cos_dec2 * sin_diff_ra;
	m2 = cos_dec1 * sin_dec2 - sin_dec1 * cos_dec2*cos_diff_ra;
	m = sqrt(m1*m1 + m2 * m2);

	n = sin_dec1 * sin_dec2 + cos_dec1 * cos_dec2*cos_diff_ra;

	sep_radian = fabs(atan2(m, n));
}


void com_distance(const double low_z, const double high_z, const double omg_m, const double omg_lam, 
double &result, const double precision_thresh, const bool integ_only)
{
	// Int f(x) = (f_1+f_2)/2*dx + (f_2+f_3)/2*dx + (f_3+f_4)/2*dx...
	//				= (f_1 + 2*f_2 + 2*f_3 + .. 2*f_{n-1} + f_n)/2*dx
	//				= ((f_1+f_n)/2 + f_2 + f_3 + .. + f_{n-1})*dx
	int i, j, bin_num, bin_num_0;
	int max_iters = 50, iters = 0;
	double res_1, res_2, d_res;
	double dz, z;
	double scale = 1000 * C_0_hat;

	// run first time, if it converges, it returns the result
	// else, do the while block
	bin_num = 5;
	dz = (high_z - low_z) / (bin_num - 1);
	res_1 = 0;
	// (f_1 + f_n)/2
	res_2 = (1. / sqrt(pow(1 + low_z, 3)*omg_m + omg_lam) + 1. / sqrt(pow(1 + high_z, 3)*omg_m + omg_lam))*0.5;
	// add the median value
	for (j = 1; j < bin_num - 1; j++)
	{
		z = low_z + j * dz;
		res_2 = res_2 + 1. / sqrt(pow(1 + z, 3)*omg_m + omg_lam);
	}
	while (true)
	{
		d_res = fabs((res_2* dz - res_1)*scale);
		if (d_res < precision_thresh)
		{
			break;
		}
		if (iters > max_iters)
		{
			std::cout << "com_dist() Max iteration, doesn't converge, " << d_res << "( " << precision_thresh << " )." << std::endl;
			exit(0);
		}
		res_1 = res_2 * dz;

		dz = dz * 0.5;

		for (j = 0; j < bin_num - 1; j++)
		{
			z = low_z + (2 * j + 1) * dz;
			res_2 = res_2 + 1. / sqrt(pow(1 + z, 3)*omg_m + omg_lam);
		}
		bin_num = bin_num + bin_num - 1;
		iters++;
	}
	if (integ_only)
	{
		result = (res_2 * dz + res_1)*0.5;
	}
	else
	{
		result = (res_2 * dz + res_1)*0.5*scale;
	}
}

/*void com_distance(const double low_z, const double high_z, const double omg_m,
 const double omg_lam, double &result, const double precision_thresh, int integ_tag)
{
	int i, j, bin_num;
	int max_run = 500;
	double result_1, result_2, d_res;
	double dz, z1, z2;

	bin_num = 20;
	result_1 = 0;
	for (i = 0; i < max_run; i++)
	{
		dz = (high_z - low_z) / bin_num;
		result_2 = 0;
		for (j = 0; j < bin_num; j++)
		{
			z1 = low_z + j * dz;
			z2 = z1 + dz;
			result_2 = result_2 + 1. / sqrt(pow(1 + z1, 3)*omg_m + omg_lam) + 1. / sqrt(pow(1 + z2, 3)*omg_m + omg_lam);
			//result_2 = result_2 + 1. / sqrt(z1*omg_m + pow(z1, 4)*omg_lam) + 1. / sqrt(z2*omg_m + pow(z2, 4)*omg_lam);
		}
		if (integ_tag == 0)
		{
			result_2 = result_2 * dz*0.5 * 1000 * C_0_hat;
		}
		else
		{
			result_2 = result_2 * dz*0.5;
		}
		if (fabs(result_2 - result_1) <= precision_thresh)
		{
			// comoving distance  [ Mpc/h ]
			result = (result_2 + result_1) *0.5;
			break;
		}
		else
		{
			result_1 = result_2;
			bin_num *= 2;
		}
	}
	if (i == max_run)
	{
		std::cout << "Max iteration, doesn't converge, " << result_2 - result_1 << "." << std::endl;
		exit(0);
	}
}*/


/////////////////////////  GGL part  /////////////////////////////////////////
void ggl_initialize(ggl_data_info *data_info)
{   
    int i;
    int cal_tag = 0;

    data_info->len_ra_col = 0;
    data_info->len_dec_col = 1;
    data_info->len_cos_dec_col = 2;
    data_info->len_z_col = 3;
    data_info->len_com_dist_col = 4;
    data_info->len_prop_dist_col = 5;
    data_info->len_jackid_col = 6;
    
    data_info->src_mg1_col = 0;
    data_info->src_mg2_col = 1;
    data_info->src_mn_col = 2;
    data_info->src_mu_col = 3;
    data_info->src_mv_col = 4;

    data_info->src_ra_col = 5;
    data_info->src_dec_col = 6;
    data_info->src_cos_dec_col = 7;
    data_info->src_z_col = 8;
    data_info->src_zerr_col = 9;
    data_info->src_prop_dist_col = 10;
    

    // no len/src data array exists in memory
    data_info->len_expo_read_tag = 0;
    data_info->src_expo_read_tag = 0;

    data_info->back_dz = 0.1;

    data_info->pos_inform_num = 4;

    data_info->crit_coeff = 1662895.2081868195;

    sprintf(data_info->ggl_log_path,"%s/log/log_%d.dat", data_info->ggl_total_path, data_info->rank);
    sprintf(data_info->ggl_result_path,"%s/result/result_cache_%d.hdf5", data_info->ggl_total_path, data_info->rank);

    ggl_read_list(data_info);

    ggl_read_pdf_inform(data_info);

#ifdef GGL_DELTA_SIGMA
    cal_tag ++;
#endif 
#ifdef GGL_GAMMA_T
    cal_tag ++;
#endif
    if(cal_tag == 0)
    {   
        sprintf(data_info->ggl_log_inform,"Nothing to do.\n");
        std::cout<<data_info->ggl_log_inform;
        exit(0);
    }

    if(data_info->rank == 0)
    {  
#ifdef GGL_DELTA_SIGMA
        sprintf(data_info->ggl_log_inform,"\nCalculate DELTA_SIGMA. G bins %d. Guess: %d\n", data_info->mg_sigma_bin_num,data_info->pdf_sigma_num);
        std::cout<<data_info->ggl_log_inform;
        show_arr(data_info->mg_sigma_bin,1,data_info->mg_sigma_bin_num+1);
        // std::cout<<std::endl;
        // show_arr(data_info->delta_sigma_guess,1,data_info->pdf_sigma_num);
        // std::cout<<std::endl;

        sprintf(data_info->ggl_log_inform,"\nDELTA_SIGMA Chi len\nof one pts: %d\nof %d pts: %d\nof %d jackknife: %d\n", 
        data_info->chi_sigma_theta_block_len, 
        data_info->signal_pts_num, 
        data_info->chi_sigma_signal_block_len,
        data_info->jack_num, data_info->total_chi_sigma_len);
        std::cout<<data_info->ggl_log_inform;
#endif 
#ifdef GGL_GAMMA_T
        sprintf(data_info->ggl_log_inform,"\nCalculate Tangential shear. G bins %d. Guess: %d\n", data_info->mg_gt_bin_num, data_info->pdf_gt_num);
        std::cout<<data_info->ggl_log_inform;
        show_arr(data_info->mg_gt_bin,1,data_info->mg_gt_bin_num+1);
        // std::cout<<std::endl;
        // show_arr(data_info->gt_guess,1,data_info->pdf_gt_num);
        // std::cout<<std::endl;

        sprintf(data_info->ggl_log_inform,"\nTangential shear Chi len\nof one pts: %d\nof %d pts: %d\nof %d jackknife: %d\n", 
        data_info->chi_gt_theta_block_len, 
        data_info->signal_pts_num, 
        data_info->chi_gt_signal_block_len,
        data_info->jack_num, data_info->total_chi_gt_len);
        std::cout<<data_info->ggl_log_inform;
#endif  
        sprintf(data_info->ggl_log_inform,"\nSeparation bins: %d\n", data_info->sep_bin_num);
        std::cout<<data_info->ggl_log_inform;
        show_arr(data_info->separation_bin,1,data_info->sep_bin_num+1);
        std::cout<<std::endl;
    }
}

// void ggl_task_prepare(ggl_data_info *data_info)
// {
//     int i, j, k, m, n;
//     MY_FLOAT theta, theta_max;

//     for(i=0; i<data_info->len_expo_num; i++)
//     {   
        
//         for(j=i; j<data_info->src_expo_num; j++)
//         {   
//             theta = data_info->len_nearest_dist[i]/data_info->separation_bin[data_info->sep_bin_num+1]*1.2
//              + data_info->len_width_informs[i];
//             if(k == 1)
//             {
//                 data_info->task_len_expo_labels.push_back(i);
//                 data_info->task_src_expo_labels.push_back(j);
//             }
//         }
//     }
//     data_info->task_expo_num = data_info->task_len_expo_labels.size();
// }

void ggl_read_list(ggl_data_info *data_info)
{   
    int i;
    // read fore/background exposure informs
    sprintf(data_info->ggl_foreground_inform_path, "%s/cata/foreground/foreground_source_list.dat", data_info->ggl_total_path);
    sprintf(data_info->ggl_background_inform_path, "%s/cata/background/background_source_list.dat", data_info->ggl_total_path);
    line_count(data_info->ggl_foreground_inform_path, data_info->len_expo_num);
    line_count(data_info->ggl_background_inform_path, data_info->src_expo_num);
    
    if(data_info->rank == 0)
    {   
        std::cout<<data_info->ggl_total_path<<std::endl;
        sprintf(data_info->ggl_log_inform, "Foreground %d expos.  Background %d expos\n", data_info->len_expo_num, data_info->src_expo_num);
        std::cout<<data_info->ggl_log_inform;
    }

    // data_info->len_nearest_dist = new MY_FLOAT[data_info->len_expo_num];

    for(i=0; i<data_info->len_expo_num; i++)
    {
         data_info->len_expo_path[i] = new char[450];
         data_info->len_expo_name[i] = new char[50];
    }
    data_info->len_data_row = new int[data_info->len_expo_num];
    data_info->len_expo_jackid = new int[data_info->len_expo_num];
    data_info->src_expo_needed_tag = new int[data_info->src_expo_num];
    initialize_arr(data_info->src_expo_needed_tag, data_info->src_expo_num,-1);


    for(i=0; i<data_info->src_expo_num; i++)
    { 
        data_info->src_expo_path[i] = new char[450];
        data_info->src_expo_name[i] = new char[50];
        data_info->src_pos_informs[i] = new MY_FLOAT[data_info->pos_inform_num];
    }
    data_info->src_data_row = new int[data_info->src_expo_num];
    
    ggl_read_len_list(data_info->ggl_foreground_inform_path, data_info);
    ggl_read_src_list(data_info->ggl_background_inform_path, data_info);
}

void ggl_read_len_list(char *file_path, ggl_data_info* data_info)
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

        strs >> data_info->len_expo_path[line_count];
        // file name
        strs >> data_info->len_expo_name[line_count];
        // galaxy count
        strs >> data_info->len_data_row[line_count];
        // strs >> data_info->len_nearest_dist[line_count];
        // jack id
        strs >> data_info->len_expo_jackid[line_count];

        // std::cout << str << std::endl;
        line_count += 1;
    }
    infile.close();
}

void ggl_read_src_list(char *file_path, ggl_data_info* data_info)
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

        strs >> data_info->src_expo_path[line_count];
        // file name
        strs >> data_info->src_expo_name[line_count];
        // galaxy count
        strs >> data_info->src_data_row[line_count];
        // ra of the exposure center
        strs >> data_info->src_pos_informs[line_count][0];
        // dec of the exposure center
        strs >> data_info->src_pos_informs[line_count][1];
        // cos(dec) of the exposure center
        strs >> data_info->src_pos_informs[line_count][2];
        // half of the diagonal 
        strs >> data_info->src_pos_informs[line_count][3];

        // std::cout << str << std::endl;
        line_count += 1;
    }
    infile.close();
}

void ggl_read_len_exp(ggl_data_info *data_info, int len_expo_label)
{
    if(data_info->len_expo_read_tag == 1)
    {delete[] data_info->len_expo_data;}

    char set_name[50];
    int size;
    sprintf(set_name,"/data");
    
    read_h5_datasize(data_info->len_expo_path[len_expo_label], set_name,size);
    data_info->len_expo_data = new MY_FLOAT[size];
    read_h5(data_info->len_expo_path[len_expo_label], set_name, data_info->len_expo_data);

    data_info->len_data_col = size/data_info->len_data_row[len_expo_label];

    data_info->len_expo_read_tag = 1;
}

void ggl_read_src_exp(ggl_data_info *data_info, int src_expo_label)
{
    if(data_info->src_expo_read_tag == 1)
    {delete[] data_info->src_expo_data;}

    char set_name[50];
    int size;
    sprintf(set_name,"/data");
    
    read_h5_datasize(data_info->src_expo_path[src_expo_label], set_name, size);
    data_info->src_expo_data = new MY_FLOAT[size];
    read_h5(data_info->src_expo_path[src_expo_label], set_name, data_info->src_expo_data);
    
    data_info->src_data_col = size/data_info->src_data_row[src_expo_label];

    data_info->src_expo_read_tag = 1;
}


void ggl_read_pdf_inform(ggl_data_info *data_info)
{
    sprintf(data_info->ggl_pdf_inform_path,"%s/cata/pdf_inform.hdf5", data_info->ggl_total_path);

    sprintf(data_info->set_name,"/separation_bin");
    read_h5_datasize(data_info->ggl_pdf_inform_path, data_info->set_name,data_info->sep_bin_num);
    data_info->separation_bin = new MY_FLOAT[data_info->sep_bin_num];
    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->separation_bin);
    data_info->signal_pts_num = data_info->sep_bin_num - 1;
    data_info->sep_bin_num = data_info->sep_bin_num - 1;

    data_info->sub_signal_count_len = 3*data_info->signal_pts_num;
    data_info->total_signal_count_len = data_info->sub_signal_count_len*(data_info->jack_num+1);

    data_info->worker_sub_signal_count = new double[data_info->sub_signal_count_len];
    data_info->worker_total_signal_count = new double[data_info->total_signal_count_len];

    initialize_arr(data_info->worker_total_signal_count, data_info->total_signal_count_len, 0);


#ifdef GGL_DELTA_SIGMA
    // read G bins for delta sigma
    sprintf(data_info->set_name,"/mg_sigma_bin");
    read_h5_datasize(data_info->ggl_pdf_inform_path, data_info->set_name,data_info->mg_sigma_bin_num);
    data_info->mg_sigma_bin = new MY_FLOAT[data_info->mg_sigma_bin_num];
    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->mg_sigma_bin);
    data_info->mg_sigma_bin_num = data_info->mg_sigma_bin_num - 1;

    // read delta sigma guesses
    sprintf(data_info->set_name,"/delta_sigma_guess");
    read_h5_datasize(data_info->ggl_pdf_inform_path, data_info->set_name,data_info->pdf_sigma_num);
    data_info->delta_sigma_guess = new double[data_info->pdf_sigma_num];
    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->delta_sigma_guess);

    data_info->chi_sigma_theta_block_len = data_info->pdf_sigma_num*data_info->mg_sigma_bin_num;
    data_info->chi_sigma_signal_block_len = data_info->signal_pts_num*data_info->chi_sigma_theta_block_len;
    data_info->total_chi_sigma_len = data_info->chi_sigma_signal_block_len*(data_info->jack_num+1);

    data_info->worker_sub_chi_sigma_count = new double[2*data_info->chi_sigma_signal_block_len];
    data_info->worker_total_chi_sigma_count = new double[2*data_info->total_chi_sigma_len];

    initialize_arr(data_info->worker_total_chi_sigma_count, 2*data_info->total_chi_sigma_len, 0);

#endif


#ifdef GGL_GAMMA_T
    sprintf(data_info->set_name,"/mg_gt_bin");
    read_h5_datasize(data_info->ggl_pdf_inform_path, data_info->set_name,data_info->mg_gt_bin_num);
    data_info->mg_gt_bin = new MY_FLOAT[data_info->mg_gt_bin_num];
    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->mg_gt_bin);
    data_info->mg_gt_bin_num = data_info->mg_gt_bin_num - 1;

    // read gt guesses
    sprintf(data_info->set_name,"/gt_guess");
    read_h5_datasize(data_info->ggl_pdf_inform_path, data_info->set_name,data_info->pdf_gt_num);
    data_info->gt_guess = new double[data_info->pdf_gt_num];   
    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->gt_guess);

    data_info->chi_gt_theta_block_len = data_info->pdf_gt_num*data_info->mg_gt_bin_num;
    data_info->chi_gt_signal_block_len = data_info->signal_pts_num*data_info->chi_gt_theta_block_len;
    data_info->total_chi_gt_len = data_info->chi_gt_signal_block_len*(data_info->jack_num+1);

    data_info->worker_sub_chi_gt_count = new double[2*data_info->chi_gt_signal_block_len];
    data_info->worker_total_chi_gt_count = new double[2*data_info->total_chi_gt_len];

    initialize_arr(data_info->worker_total_chi_gt_count, 2*data_info->total_chi_gt_len, 0);

#endif


    if(data_info->rank == 0)
    {
        data_info->total_chi_sigma_count = new double[data_info->total_chi_sigma_len];
        data_info->total_chi_gt_count = new double[data_info->total_chi_gt_len];
        data_info->total_signal_count = new double[data_info->total_signal_count_len];

        initialize_arr(data_info->total_chi_sigma_count, data_info->total_chi_sigma_len, 0);
        initialize_arr(data_info->total_chi_gt_count, data_info->total_chi_gt_len, 0);
        initialize_arr(data_info->total_signal_count, data_info->total_signal_count_len, 0);
    }

}


void ggl_find_src_needed(ggl_data_info *data_info, int len_expo_label)
{
    int ifg, ifg_row, bkg;
    MY_FLOAT max_sep_theta;
    MY_FLOAT dra, ddec, sep_theta;
    MY_FLOAT len_ra, len_dec, len_cos_dec;

    initialize_arr(data_info->src_expo_needed_tag, data_info->src_expo_num, -1);

    //std::cout<<len_expo_label<<" "<<data_info->len_data_row[len_expo_label]<<std::endl;


    for(ifg=0; ifg<data_info->len_data_row[len_expo_label]; ifg++)
    {   

        ifg_row = ifg*data_info->len_data_col;
        len_ra = data_info->len_expo_data[ifg_row + data_info->len_ra_col];
        len_dec = data_info->len_expo_data[ifg_row + data_info->len_dec_col];
        len_cos_dec = data_info->len_expo_data[ifg_row + data_info->len_cos_dec_col];

        // if stack the signal in physical coordinate,
        // 1.5 for safety, enlarges the search radius
#ifdef GGL_PROP_DIST_STACK
        max_sep_theta = 1.5*data_info->separation_bin[data_info->sep_bin_num]/
        data_info->len_expo_data[ifg_row + data_info->len_prop_dist_col]*180/Pi;
#else
        max_sep_theta = 1.5*data_info->separation_bin[data_info->sep_bin_num+1]/
        data_info->len_expo_data[ifg_row + data_info->len_com_dist_col]*180/Pi;
#endif       
        for(bkg=0; bkg<data_info->src_expo_num; bkg++)
        {
            dra = (data_info->src_pos_informs[bkg][0] - len_ra)*len_cos_dec;
            ddec = data_info->src_pos_informs[bkg][1] - len_dec;
            sep_theta = sqrt(dra*dra + ddec*ddec) - data_info->src_pos_informs[bkg][2];
            // std::cout<<data_info->src_pos_informs[bkg][0]<<" "<<data_info->src_pos_informs[bkg][1]<<" "<<len_ra<<" "<<len_dec<<std::endl;
            if(sep_theta <= max_sep_theta){data_info->src_expo_needed_tag[bkg] = 1;}
        }
    }

}

void ggl_rotation_matrix(MY_FLOAT cent_ra, MY_FLOAT cent_dec, MY_FLOAT cent_cos_dec,MY_FLOAT src_ra, MY_FLOAT src_dec, MY_FLOAT*rotation_matrix)
{
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;
    MY_FLOAT dra, ddec, delta_radius;
    
    dra = src_ra - cent_ra;
    ddec = (src_dec - cent_dec)*cent_cos_dec;
    delta_radius = sqrt(dra*dra + ddec*ddec);

    // theta is position angle
    // sin_theta cos_theta
    rotation_matrix[0] = dra/delta_radius;
    rotation_matrix[1] = ddec/delta_radius;
    // sin_2theta cos_2theta
    rotation_matrix[2] = 2*rotation_matrix[0]*rotation_matrix[1];
    rotation_matrix[3] = rotation_matrix[1]*rotation_matrix[1] - rotation_matrix[0]*rotation_matrix[0];
    // sin_4theta cos_4theta
    rotation_matrix[4] = 2*rotation_matrix[2]*rotation_matrix[3];
    rotation_matrix[5] = rotation_matrix[3]*rotation_matrix[3] - rotation_matrix[2]*rotation_matrix[2];
}   

void ggl_rotate_estimator(MY_FLOAT G1, MY_FLOAT G2, MY_FLOAT U, MY_FLOAT V, MY_FLOAT *rotation_matrix, MY_FLOAT &Gt, MY_FLOAT &Gx, MY_FLOAT &Ut)
{
    Gt = G1*rotation_matrix[3] - G2*rotation_matrix[2];
    Gx = G1*rotation_matrix[2] + G2*rotation_matrix[3];
    Ut = U*rotation_matrix[5] - V*rotation_matrix[4];
    // Ux = U*rotation_matrix[4] + V*rotation_matrix[5];
}

void ggl_fast_hist(MY_FLOAT *bins, int bin_num, MY_FLOAT val, int pre_bin_tag, int &bin_tag)
{   
    int i;
    if(val >= bins[pre_bin_tag])
    {
        for(i=pre_bin_tag; i<bin_num; i++)
        {
            if(val>=bins[i] and val<bins[i+1]){bin_tag = i;}
        }
    }
    else
    {
        for(i=pre_bin_tag; i>0; i--)
        {
            if(val>=bins[i-1] and val<bins[i]){bin_tag = i;}
        }
    }
}

void ggl_find_pair(ggl_data_info *data_info, int len_expo_label)
{
    int ibkg, bkg, ibkg_row, ifg, ifg_row;
    int ir, sep_bin_tag;
    int chi_sigma_pos, chi_gt_pos, chi_sigma_pos_i, chi_gt_pos_i;
    int i,j,k;
    double st, ed;

    MY_FLOAT len_ra, len_dec, len_cos_dec, src_ra, src_dec;
    MY_FLOAT len_z, len_dist, src_z, src_dist;
    MY_FLOAT src_z_err, len_z_dz;
    MY_FLOAT dra, ddec, delta_radius;
    MY_FLOAT sep_dist, sep_theta;
    MY_FLOAT sigma_crit, coeff;
    MY_FLOAT src_mg1, src_mg2, src_mn, src_mv, src_mu;
    MY_FLOAT src_mg1_rot, src_mg2_rot, src_mu_rot, src_mv_rot;
    MY_FLOAT temp_mgt, temp_mgx, temp_mnut, temp_mnux,temp_mnut_g, temp_mnux_g;

    int pre_pdf_bin_tag1, pdf_bin_tag1;
    int pre_pdf_bin_tag2, pdf_bin_tag2;

    MY_FLOAT rotation_mat[6];
    int pair_count = 0;

    st = clock();

    initialize_arr(data_info->worker_sub_chi_sigma_count, 2*data_info->chi_sigma_signal_block_len, 0);
    initialize_arr(data_info->worker_sub_chi_gt_count, 2*data_info->chi_gt_signal_block_len, 0);
    initialize_arr(data_info->worker_sub_signal_count, data_info->sub_signal_count_len, 0);

    ggl_read_len_exp(data_info, len_expo_label);

    ggl_find_src_needed(data_info, len_expo_label);

    for(bkg=0; bkg<data_info->src_expo_num; bkg++)
    {   
        if(data_info->src_expo_needed_tag[bkg]< 1){continue;}

        ggl_read_src_exp(data_info, bkg);

        for(ifg=0; ifg<data_info->len_data_row[len_expo_label]; ifg++)
        {    
            ifg_row = ifg*data_info->len_data_col;

            len_ra = data_info->len_expo_data[ifg_row + data_info->len_ra_col];
            len_dec = data_info->len_expo_data[ifg_row + data_info->len_dec_col];
            len_cos_dec = data_info->len_expo_data[ifg_row + data_info->len_cos_dec_col];

            len_z = data_info->len_expo_data[ifg_row + data_info->len_z_col];
            len_z_dz = len_z + data_info->back_dz;

            len_dist = data_info->len_expo_data[ifg_row + data_info->len_prop_dist_col];

            // stacking in physical or comoving coordinate
#ifdef GGL_PROP_DIST_STACK
            coeff = data_info->crit_coeff/len_dist;
#else
            coeff = data_info->crit_coeff/len_dist/(1+len_z)/(1+len_z);
#endif            
            for(ibkg=0; ibkg<data_info->src_data_row[bkg]; ibkg++)
            {   
                ibkg_row = ibkg*data_info->src_data_col;
                src_z = data_info->src_expo_data[ibkg_row + data_info->src_z_col];
                src_z_err = src_z - data_info->src_expo_data[ibkg_row + data_info->src_zerr_col];
                if(len_z >= src_z_err or len_z_dz >= src_z){continue;}

                src_ra = data_info->src_expo_data[ibkg_row + data_info->src_ra_col];
                src_dec = data_info->src_expo_data[ibkg_row + data_info->src_dec_col];
                src_dist = data_info->src_expo_data[ibkg_row + data_info->src_prop_dist_col];

                sigma_crit = coeff*src_dist/(src_dist - len_dist);

                separation_angle_1(len_ra, len_dec, src_ra, src_dec, sep_theta);

#ifdef GGL_PROP_DIST_STACK
                sep_dist = sep_theta*data_info->len_expo_data[ifg_row + data_info->len_prop_dist_col];
#else
                sep_dist = sep_theta*data_info->len_expo_data[ifg_row + data_info->len_com_dist_col];
#endif
                sep_bin_tag = -1;
                for(ir=0; ir<data_info->sep_bin_num; ir++)
                {
                    if(sep_dist >= data_info->separation_bin[ir] and sep_dist< data_info->separation_bin[ir+1])
                    { sep_bin_tag = ir; break; }
                }
                // if(sep_bin_tag<2 and sep_bin_tag> -1)
                // {sprintf(data_info->ggl_log_inform,"%f %f %f %f %f %f %f %f %d\n",
                // len_ra, len_dec, src_ra, src_dec, sep_theta,data_info->len_expo_data[ifg_row + data_info->len_com_dist_col],len_z, sep_dist,sep_bin_tag);
                // std::cout<<data_info->ggl_log_inform;}

                if(sep_bin_tag > -1)
                {   
                    pair_count ++;

                    data_info->worker_sub_signal_count[sep_bin_tag] += sep_theta;
                    data_info->worker_sub_signal_count[data_info->signal_pts_num + sep_bin_tag] += sep_dist;
                    data_info->worker_sub_signal_count[data_info->signal_pts_num*2 + sep_bin_tag] += 1;

                    // if(sep_bin_tag == 0)std::cout<<"Find 0 pairs "<<sep_bin_tag<<" "<<len_expo_label<<std::endl;
                    // rotation, sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta
                    ggl_rotation_matrix(len_ra,len_dec, len_cos_dec, src_ra, src_dec, rotation_mat);
                    src_mg1 = data_info->src_expo_data[ibkg_row + data_info->src_mg1_col];
                    src_mg2 = data_info->src_expo_data[ibkg_row + data_info->src_mg2_col];
                    src_mn = data_info->src_expo_data[ibkg_row + data_info->src_mn_col];
                    src_mu = data_info->src_expo_data[ibkg_row + data_info->src_mu_col];
                    src_mv = data_info->src_expo_data[ibkg_row + data_info->src_mv_col];
                    
                    ggl_rotate_estimator(src_mg1, src_mg2, src_mu, src_mv, rotation_mat, src_mg1_rot, src_mg2_rot, src_mu_rot);

                    src_mg1_rot *= sigma_crit;
                    src_mg2_rot *= sigma_crit;
                    temp_mnut = src_mn + src_mu_rot;
                    temp_mnux = src_mn - src_mu_rot;
                    // pdf
                    
 
#ifdef GGL_DELTA_SIGMA
                    pre_pdf_bin_tag1 = 0;
                    pre_pdf_bin_tag2 = 0;

                    chi_sigma_pos = sep_bin_tag*data_info->chi_sigma_theta_block_len;
                    for(i=0; i<data_info->pdf_sigma_num; i++)
                    { 
                        // \Delta\Sigma(R) & \Delta\Sigma(R)_x
                        chi_sigma_pos_i = chi_sigma_pos + i*data_info->mg_sigma_bin_num;
                        temp_mgt = src_mg1_rot - data_info->delta_sigma_guess[i]*temp_mnut;
                        temp_mgx = src_mg2_rot - data_info->delta_sigma_guess[i]*temp_mnux;
                        ggl_fast_hist(data_info->mg_sigma_bin, data_info->mg_sigma_bin_num, temp_mgt, pre_pdf_bin_tag1, pdf_bin_tag1);
                        ggl_fast_hist(data_info->mg_sigma_bin, data_info->mg_sigma_bin_num, temp_mgx, pre_pdf_bin_tag2, pdf_bin_tag2);
                        
                        data_info->worker_sub_chi_sigma_count[chi_sigma_pos_i + pdf_bin_tag1] +=1;
                        data_info->worker_sub_chi_sigma_count[data_info->chi_sigma_signal_block_len + chi_sigma_pos_i + pdf_bin_tag2] +=1;
                        
                        pre_pdf_bin_tag1 = pdf_bin_tag1;
                        pre_pdf_bin_tag2 = pdf_bin_tag2;
                    }
#endif

#ifdef GGL_GAMMA_T
                    pre_pdf_bin_tag1 = 0;
                    pre_pdf_bin_tag2 = 0;

                    temp_mnut_g = temp_mnut*sigma_crit;
                    temp_mnux_g = temp_mnux*sigma_crit;

                    chi_gt_pos = sep_bin_tag*data_info->chi_gt_theta_block_len;
                    for(i=0; i<data_info->pdf_gt_num; i++)
                    { 
                        // \gamma_t & \gamma_x
                        chi_gt_pos_i = chi_gt_pos + i*data_info->mg_gt_bin_num;
                        temp_mgt = src_mg1_rot - data_info->gt_guess[i]*temp_mnut_g;
                        temp_mgx = src_mg2_rot - data_info->gt_guess[i]*temp_mnux_g;
                        ggl_fast_hist(data_info->mg_gt_bin, data_info->mg_gt_bin_num, temp_mgt, pre_pdf_bin_tag1, pdf_bin_tag1);
                        ggl_fast_hist(data_info->mg_gt_bin, data_info->mg_gt_bin_num, temp_mgx, pre_pdf_bin_tag2, pdf_bin_tag2);
                        
                        data_info->worker_sub_chi_gt_count[chi_gt_pos_i + pdf_bin_tag1] +=1;
                        data_info->worker_sub_chi_gt_count[data_info->chi_gt_signal_block_len + chi_gt_pos_i + pdf_bin_tag2] +=1;
                        
                        pre_pdf_bin_tag1 = pdf_bin_tag1;
                        pre_pdf_bin_tag2 = pdf_bin_tag2;
                    }
#endif                  
                }
            }
        }
    }

    // add the count to the total count array, according to the jack id
    for(i=0; i<data_info->jack_num+1; i++)
    {   
        if(i == data_info->len_expo_jackid[len_expo_label]){continue;}

        for(j=0; j<2*data_info->chi_sigma_signal_block_len; j++)
        {
            data_info->worker_total_chi_sigma_count[i*data_info->chi_sigma_signal_block_len*2 + j] += data_info->worker_sub_chi_sigma_count[j];
        }
        for(j=0; j<2*data_info->chi_gt_signal_block_len; j++)
        {
            data_info->worker_total_chi_gt_count[i*data_info->chi_gt_signal_block_len*2 + j] += data_info->worker_sub_chi_gt_count[j];
        }
        for(j=0; j<data_info->sub_signal_count_len; j++)
        {
            data_info->worker_total_signal_count[i*data_info->sub_signal_count_len + j] += data_info->worker_sub_signal_count[j];
        }
    }

    ed = clock();
    char times[100];
    sprintf(times,"worker %d. %d Lens file. Finished in %.2f sec. %d Lenses, %d pairs",
            data_info->rank, len_expo_label, (ed-st)/CLOCKS_PER_SEC, data_info->len_data_row[len_expo_label], pair_count);
    std::cout<<times<<std::endl;
}

void ggl_collect_chi(ggl_data_info *data_info)
{
    int i, j;
    MPI_Status status;
    char set_name[50];

    if (data_info->rank > 0)
    {
        MPI_Send(data_info->worker_total_signal_count, data_info->total_signal_count_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD);
        sprintf(set_name,"/data");
        sprintf(data_info->ggl_result_path,"%s/result/chi_theta_cache_%d.hdf5", data_info->ggl_total_path, data_info->rank);

        write_h5(data_info->ggl_result_path, set_name, data_info->worker_total_signal_count,
                 data_info->jack_num+1, data_info->total_signal_count_len, true);
    }
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->worker_total_signal_count, data_info->total_signal_count_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
        
            for(j=0;j<data_info->total_signal_count_len;j++)
            {data_info->total_signal_count[j] += data_info->worker_total_signal_count[j];}
        }
        
        sprintf(set_name,"/data");
        sprintf(data_info->ggl_result_path,"%s/result/chi_signal_cache.hdf5", data_info->ggl_total_path);

        write_h5(data_info->ggl_result_path, set_name, data_info->total_signal_count,
                 data_info->jack_num+1,data_info->total_signal_count_len, true);
    }

#ifdef GGL_DELTA_SIGMA

    if (data_info->rank > 0)
    {
        MPI_Send(data_info->worker_total_chi_sigma_count, 2*data_info->total_chi_sigma_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD);
        sprintf(set_name,"/data");
        sprintf(data_info->ggl_result_path,"%s/result/chi_sigma_cache_%d.hdf5", data_info->ggl_total_path, data_info->rank);

        write_h5(data_info->ggl_result_path, set_name, data_info->worker_total_chi_sigma_count,
                 data_info->jack_num+1,2*data_info->total_chi_sigma_len, true);
    }
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->worker_total_chi_sigma_count, 2*data_info->total_chi_sigma_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
        
            for(j=0;j<2*data_info->total_chi_sigma_len;j++)
            {data_info->total_chi_sigma_count[j] += data_info->worker_total_chi_sigma_count[j];}
        }
        
        sprintf(set_name,"/data");
        sprintf(data_info->ggl_result_path,"%s/result/chi_sigma_cache.hdf5", data_info->ggl_total_path);

        write_h5(data_info->ggl_result_path, set_name, data_info->total_chi_sigma_count,
                 data_info->jack_num+1,2*data_info->total_chi_sigma_len, true);
    }
#endif

#ifdef GGL_GAMMA_T

    if (data_info->rank > 0)
    {
        MPI_Send(data_info->worker_total_chi_gt_count, 2*data_info->total_chi_gt_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD);
        sprintf(set_name,"/data");
        sprintf(data_info->ggl_result_path,"%s/result/chi_gt_cache_%d.hdf5", data_info->ggl_total_path, data_info->rank);

        write_h5(data_info->ggl_result_path, set_name, data_info->worker_total_chi_sigma_count,
                 data_info->jack_num+1,2*data_info->total_chi_sigma_len, true);
    }
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->worker_total_chi_gt_count, 2*data_info->total_chi_gt_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
        
            for(j=0;j<2*data_info->total_chi_gt_len;j++)
            {data_info->total_chi_gt_count[j] += data_info->worker_total_chi_gt_count[j];}
        }
        
        sprintf(set_name,"/data");
        sprintf(data_info->ggl_result_path,"%s/result/chi_gt_cache.hdf5", data_info->ggl_total_path);

        write_h5(data_info->ggl_result_path, set_name, data_info->total_chi_gt_count,
                 data_info->jack_num+1,2*data_info->total_chi_gt_len, true);
    }
#endif

}

void ggl_cal_signals(ggl_data_info * data_info)
{   
    sprintf(data_info->ggl_log_inform,"start calculate\n");
    std::cout<<data_info->ggl_log_inform;

    char set_name[50];
    int i, j, st_c, st_j;

    double *theta = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
    double *radius = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
    double *count = new double[(data_info->jack_num+1)*data_info->signal_pts_num];

    double *signals = new double[data_info->signal_pts_num];
    double *signals_err = new double[data_info->signal_pts_num];

    for(i=0; i<data_info->jack_num+1; i++)
    { 
        // theta
        st_c = i*data_info->signal_pts_num;
        st_j = st_c*3;
        for(j=0;j<data_info->signal_pts_num;j++)
        { theta[st_c + j] = data_info->total_signal_count[st_j + j];}

        // radius
        st_j = st_c*3 + data_info->signal_pts_num;
        for(j=0;j<data_info->signal_pts_num;j++)
        { radius[st_c + j] = data_info->total_signal_count[st_j + j];}

        // count
        st_j = st_c*3 + data_info->signal_pts_num*2;
        for(j=0;j<data_info->signal_pts_num;j++)
        { count[st_c + j] = data_info->total_signal_count[st_j + j];}           
    }

    sprintf(data_info->ggl_result_path,"%s/result/result.hdf5", data_info->ggl_total_path);
    sprintf(set_name,"/theta");
    write_h5(data_info->ggl_result_path, set_name, theta,
            data_info->jack_num+1,data_info->signal_pts_num, true);
    sprintf(set_name,"/radius");
    write_h5(data_info->ggl_result_path, set_name, radius,
            data_info->jack_num+1,data_info->signal_pts_num, false);
    sprintf(set_name,"/count");
    write_h5(data_info->ggl_result_path, set_name, count,
            data_info->jack_num+1,data_info->signal_pts_num, false);


#ifdef GGL_DELTA_SIGMA

    double *temp_sigma_count = new double[data_info->chi_sigma_signal_block_len];
    double *delta_sigma_t = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
    double *delta_sigma_x = new double[(data_info->jack_num+1)*data_info->signal_pts_num];

    for(i=0; i<data_info->jack_num+1; i++)
    {   
        st_j = i*data_info->signal_pts_num;

        // delta_sigma_t
        st_c = i*data_info->chi_sigma_signal_block_len*2;
        for(j=0; j<data_info->chi_sigma_signal_block_len; j++)
        { temp_sigma_count[j] = data_info->total_chi_sigma_count[st_c+j];}

        // show_arr(temp_count, 1, data_info->chi_signal_block_len);
        ggl_pdf_signals(temp_sigma_count, data_info->delta_sigma_guess, data_info->pdf_sigma_num, 
                        data_info->mg_sigma_bin_num, data_info->signal_pts_num, signals, signals_err);
        
        for(j=0; j<data_info->signal_pts_num; j++)
        {delta_sigma_t[st_j + j] = signals[j];}


        // delta_sigma_x
        st_c = i*data_info->chi_sigma_signal_block_len*2 + data_info->chi_sigma_signal_block_len;
        for(j=0;j<data_info->chi_sigma_signal_block_len;j++)
        { temp_sigma_count[j] = data_info->total_chi_sigma_count[st_c+j];}

        ggl_pdf_signals(temp_sigma_count, data_info->delta_sigma_guess, data_info->pdf_sigma_num, 
                        data_info->mg_sigma_bin_num, data_info->signal_pts_num, signals, signals_err);
        
        for(j=0; j<data_info->signal_pts_num; j++)
        {delta_sigma_x[st_j + j] = signals[j];}
    }

    sprintf(set_name,"/delta_t");
    write_h5(data_info->ggl_result_path, set_name, delta_sigma_t,
            data_info->jack_num+1,data_info->signal_pts_num, false);
    sprintf(set_name,"/delta_x");
    write_h5(data_info->ggl_result_path, set_name, delta_sigma_x,
            data_info->jack_num+1,data_info->signal_pts_num, false);

#endif

#ifdef GGL_GAMMA_T
    double *temp_gt_count = new double[data_info->chi_gt_signal_block_len];
    double *gt = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
    double *gx = new double[(data_info->jack_num+1)*data_info->signal_pts_num];

    for(i=0; i<data_info->jack_num+1; i++)
    {        
        // g_t
        st_c = i*data_info->chi_gt_signal_block_len*2;
        for(j=0; j<data_info->chi_gt_signal_block_len; j++)
        { temp_gt_count[j] = data_info->total_chi_gt_count[st_c+j];}

        ggl_pdf_signals(temp_gt_count, data_info->delta_sigma_guess, data_info->pdf_gt_num, 
                        data_info->mg_gt_bin_num, data_info->signal_pts_num, signals, signals_err);
        
        for(j=0; j<data_info->signal_pts_num; j++)
        {gt[st_j + j] = signals[j];}


        // g_x
        st_c = i*data_info->chi_gt_signal_block_len*2 + data_info->chi_gt_signal_block_len;
        for(j=0; j<data_info->chi_gt_signal_block_len; j++)
        { temp_gt_count[j] = data_info->total_chi_gt_count[st_c+j];}
        ggl_pdf_signals(temp_gt_count, data_info->delta_sigma_guess, data_info->pdf_gt_num, 
                        data_info->mg_gt_bin_num, data_info->signal_pts_num, signals, signals_err);
        
        for(j=0; j<data_info->signal_pts_num; j++)
        {gx[st_j + j] = signals[j];}
    }

    sprintf(set_name,"/g_t");
    write_h5(data_info->ggl_result_path, set_name, gt,
            data_info->jack_num+1,data_info->signal_pts_num, false);
    sprintf(set_name,"/g_x");
    write_h5(data_info->ggl_result_path, set_name, gx,
            data_info->jack_num+1,data_info->signal_pts_num, false);

#endif

    sprintf(data_info->ggl_log_inform,"finish calculate and save result.\n");
    std::cout<<data_info->ggl_log_inform;      
}

void ggl_pdf_signals(double *chi_count, double*pdf_signal_guess, int pdf_guess_num, int mg_bin_num, int signal_pts_num, double *signal, double *signal_err)
{
    int i, j, k, st;
    double chisq_i;
    double signal_i, signal_err_i;
    double *temp = new double[mg_bin_num];
    double *chisq = new double[pdf_guess_num];

    for(i=0; i<signal_pts_num; i++)
    {
        for(j=0; j<pdf_guess_num; j++)
        {   
            st = i*pdf_guess_num*mg_bin_num + j*mg_bin_num;
            for(k=0; k<mg_bin_num; k++)
            {
                temp[k] = chi_count[st+k];
            }
            cal_chisq_1d(temp, mg_bin_num, chisq_i);
            chisq[j] = chisq_i;
        }
        fit_shear(pdf_signal_guess, chisq, pdf_guess_num, signal_i, signal_err_i, chisq_i, 80);
        signal[i] = signal_i;
        signal_err[i] = signal_err_i;
    }
}
