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
    data_info->len_jackid_col = 5;
    
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
    data_info->src_com_dist_col = 10;
    
    if(data_info->jack_num <= 2){data_info->jack_num=0;}

    // no len/src data array exists in memory
    data_info->len_expo_read_tag = 0;
    data_info->src_expo_read_tag = 0;

    data_info->pos_inform_num = 4;

    data_info->crit_coeff = 1.6628952007121066*1.e6;

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
        data_info->chi_sigma_theta_block_len_sub, 
        data_info->signal_pts_num, 
        data_info->chi_sigma_theta_block_len,
        data_info->jack_num, data_info->chi_sigma_jack_block_len);
        std::cout<<data_info->ggl_log_inform;
        std::cout<<"2D hist: "<<data_info->hist2d_len<<" bins. "<<data_info->hist2d_total_len<<std::endl;
#endif 
#ifdef GGL_GAMMA_T
        sprintf(data_info->ggl_log_inform,"\nCalculate Tangential shear. G bins %d. Guess: %d\n", data_info->mg_gt_bin_num, data_info->pdf_gt_num);
        std::cout<<data_info->ggl_log_inform;
        show_arr(data_info->mg_gt_bin,1,data_info->mg_gt_bin_num+1);
        // std::cout<<std::endl;
        // show_arr(data_info->gt_guess,1,data_info->pdf_gt_num);
        // std::cout<<std::endl;

        sprintf(data_info->ggl_log_inform,"\nTangential shear Chi len\nof one pts: %d\nof %d pts: %d\nof %d jackknife: %d\n", 
        data_info->chi_g_theta_block_len_sub, 
        data_info->signal_pts_num, 
        data_info->chi_g_theta_block_len,
        data_info->jack_num, data_info->chi_g_jack_block_len);
        std::cout<<data_info->ggl_log_inform;
#endif  
        sprintf(data_info->ggl_log_inform,"\nSeparation bins: %d\n", data_info->sep_bin_num);
        std::cout<<data_info->ggl_log_inform;
        show_arr(data_info->separation_bin,1,data_info->sep_bin_num+1);
        std::cout<<std::endl;
        std::cout<<data_info->ggl_total_path<<std::endl;
        sprintf(data_info->ggl_log_inform,"\nThe background have z >= fore Z + %.3f\n", data_info->back_dz);
        std::cout<<data_info->ggl_log_inform;
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
         data_info->len_expo_name[i] = new char[100];
    }
    data_info->len_data_row = new int[data_info->len_expo_num];
    data_info->len_expo_jackid = new int[data_info->len_expo_num];
    data_info->src_expo_needed_tag = new int[data_info->src_expo_num];
    initialize_arr(data_info->src_expo_needed_tag, data_info->src_expo_num,-1);


    for(i=0; i<data_info->src_expo_num; i++)
    { 
        data_info->src_expo_path[i] = new char[450];
        data_info->src_expo_name[i] = new char[100];
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
    int i, j;

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



    // read hist2d bins
    sprintf(data_info->set_name,"/hist2d_mg_sigma_bin");
    read_h5_datasize(data_info->ggl_pdf_inform_path, data_info->set_name,data_info->hist2d_mg_sigma_bin_num);
    data_info->hist2d_mg_sigma_bin = new MY_FLOAT[data_info->hist2d_mg_sigma_bin_num];

    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->hist2d_mg_sigma_bin);
    data_info->hist2d_mg_sigma_bin_num = data_info->hist2d_mg_sigma_bin_num - 1;

    sprintf(data_info->set_name,"/hist2d_mn_sigma_bin");
    read_h5_datasize(data_info->ggl_pdf_inform_path, data_info->set_name,data_info->hist2d_mn_sigma_bin_num);
    data_info->hist2d_mn_sigma_bin = new MY_FLOAT[data_info->hist2d_mn_sigma_bin_num];

    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->hist2d_mn_sigma_bin);
    data_info->hist2d_mn_sigma_bin_num = data_info->hist2d_mn_sigma_bin_num - 1;

    data_info->hist2d_len = data_info->hist2d_mn_sigma_bin_num*data_info->hist2d_mg_sigma_bin_num;
    data_info->hist2d_total_len = data_info->hist2d_len*data_info->signal_pts_num;
    // there are "signal_pts_num" blocks, each one is a 2d array (1d actually in the memory)
    data_info->hist2d_count = new double[data_info->hist2d_total_len];
    initialize_arr(data_info->hist2d_count, data_info->hist2d_total_len, 0);

    data_info->hist2d_x = new double[data_info->hist2d_total_len];
    initialize_arr(data_info->hist2d_x, data_info->hist2d_total_len, 0);

    data_info->hist2d_y = new double[data_info->hist2d_total_len];
    initialize_arr(data_info->hist2d_y, data_info->hist2d_total_len, 0);



    // read delta sigma guesses
    sprintf(data_info->set_name,"/delta_sigma_guess");
    read_h5_datasize(data_info->ggl_pdf_inform_path, data_info->set_name,data_info->pdf_sigma_num);
    data_info->delta_sigma_guess = new double[data_info->pdf_sigma_num];
    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->delta_sigma_guess);

    data_info->chi_sigma_theta_block_len_sub = data_info->pdf_sigma_num*data_info->mg_sigma_bin_num;
    data_info->chi_sigma_theta_block_len = data_info->chi_sigma_theta_block_len_sub*data_info->signal_pts_num;
    data_info->chi_sigma_jack_block_len = data_info->chi_sigma_theta_block_len*(data_info->jack_num+1);

    data_info->worker_sub_chi_sigma_tan = new double[data_info->chi_sigma_theta_block_len];
    data_info->worker_total_chi_sigma_tan = new double[data_info->chi_sigma_jack_block_len];
    data_info->worker_sub_chi_sigma_cross = new double[data_info->chi_sigma_theta_block_len];
    data_info->worker_total_chi_sigma_cross = new double[data_info->chi_sigma_jack_block_len];

    initialize_arr(data_info->worker_total_chi_sigma_tan, data_info->chi_sigma_jack_block_len, 0);
    initialize_arr(data_info->worker_total_chi_sigma_cross, data_info->chi_sigma_jack_block_len, 0);

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

    data_info->chi_g_theta_block_len_sub = data_info->pdf_gt_num*data_info->mg_gt_bin_num;
    data_info->chi_g_theta_block_len = data_info->signal_pts_num*data_info->chi_g_theta_block_len_sub;
    data_info->chi_g_jack_block_len = data_info->chi_g_theta_block_len*(data_info->jack_num+1);

    data_info->worker_sub_chi_g_tan = new double[data_info->chi_g_theta_block_len];
    data_info->worker_total_chi_g_tan = new double[data_info->chi_g_jack_block_len];
    data_info->worker_sub_chi_g_cross = new double[data_info->chi_g_theta_block_len];
    data_info->worker_total_chi_g_cross = new double[data_info->chi_g_jack_block_len];

    initialize_arr(data_info->worker_total_chi_g_tan, data_info->chi_g_jack_block_len, 0);
    initialize_arr(data_info->worker_total_chi_g_cross, data_info->chi_g_jack_block_len, 0);

#endif


    if(data_info->rank == 0)
    {
        data_info->total_chi_sigma_tan = new double[data_info->chi_sigma_jack_block_len];
        data_info->total_chi_sigma_cross = new double[data_info->chi_sigma_jack_block_len];
        data_info->total_chi_g_tan = new double[data_info->chi_g_jack_block_len];
        data_info->total_chi_g_cross = new double[data_info->chi_g_jack_block_len];
        data_info->total_signal_count = new double[data_info->total_signal_count_len];

        initialize_arr(data_info->total_chi_sigma_tan, data_info->chi_sigma_jack_block_len, 0);
        initialize_arr(data_info->total_chi_sigma_cross, data_info->chi_sigma_jack_block_len, 0);
        initialize_arr(data_info->total_chi_g_tan, data_info->chi_g_jack_block_len, 0);
        initialize_arr(data_info->total_chi_g_cross, data_info->chi_g_jack_block_len, 0);
        initialize_arr(data_info->total_signal_count, data_info->total_signal_count_len, 0);

        data_info->hist2d_count_total = new double[data_info->hist2d_total_len];
        data_info->hist2d_x_total = new double[data_info->hist2d_total_len];
        data_info->hist2d_y_total = new double[data_info->hist2d_total_len];

        initialize_arr(data_info->hist2d_count_total, data_info->hist2d_total_len, 0);
        initialize_arr(data_info->hist2d_x_total, data_info->hist2d_total_len, 0);
        initialize_arr(data_info->hist2d_y_total, data_info->hist2d_total_len, 0);

    }

}


void ggl_find_src_needed(ggl_data_info *data_info, int len_expo_label)
{
    int ifg, ifg_row, bkg;
    MY_FLOAT max_sep_theta;
    MY_FLOAT dra, ddec, sep_theta;
    MY_FLOAT len_ra, len_dec, len_cos_dec, len_dist;

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
        len_dist = data_info->len_expo_data[ifg_row + data_info->len_com_dist_col]/(1+data_info->len_expo_data[ifg_row + data_info->len_z_col])
#else
        len_dist = data_info->len_expo_data[ifg_row + data_info->len_com_dist_col];
#endif     
        max_sep_theta = 1.2*data_info->separation_bin[data_info->sep_bin_num]/len_dist*180/Pi;  
        for(bkg=0; bkg<data_info->src_expo_num; bkg++)
        {
            dra = (data_info->src_pos_informs[bkg][0] - len_ra)*len_cos_dec;
            ddec = data_info->src_pos_informs[bkg][1] - len_dec;
            sep_theta = sqrt(dra*dra + ddec*ddec) - data_info->src_pos_informs[bkg][2];
            // std::cout<<data_info->src_pos_informs[bkg][0]<<" "<<data_info->src_pos_informs[bkg][1]<<" "<<len_ra<<" "<<len_dec<<std::endl;
            if(sep_theta <= max_sep_theta)
            {
                data_info->src_expo_needed_tag[bkg] = 1;
                // std::cout<<data_info->src_pos_informs[bkg][0]<<" "<<data_info->src_pos_informs[bkg][1]<<" "<<len_ra<<" "<<len_dec<<" "
                // <<max_sep_theta<<" "<<data_info->sep_bin_num<<" "<<data_info->separation_bin[data_info->sep_bin_num + 1]<<" "<<data_info->separation_bin[data_info->sep_bin_num]<<" "<<len_dist<<std::endl;
            }
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
    int i, j;
    int tag;
    if(val >= bins[pre_bin_tag])
    {   
        j=0;
        for(i=pre_bin_tag; i<bin_num; i++) 
        {
            if(val>=bins[i] and val<bins[i+1]){tag = i;break;}
        }
    }
    else
    {
        j=1;
        for(i=pre_bin_tag; i>0; i--)
        {
            if(val>=bins[i-1] and val<bins[i]){tag = i-1;break;}
        }
    }
    bin_tag = tag;
    if(bin_tag < 0)
    {
        std::cout<<"Outlier: val="<<val<<" bin_tag="<<tag<<" pre_bin="<<bins[pre_bin_tag]<<" bin[0]="<<bins[0]<<" bin[-1]="<<bins[bin_num]<<" "<<j<<std::endl;
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
    MY_FLOAT sep_dist, sep_theta, sep_theta_;
    MY_FLOAT sigma_crit, coeff;
    MY_FLOAT src_mg1, src_mg2, src_mn, src_mv, src_mu;
    MY_FLOAT src_mg1_rot, src_mg2_rot, src_mu_rot, src_mv_rot;
    MY_FLOAT temp_mgt, temp_mgx, temp_mnut, temp_mnux,temp_mnut_g, temp_mnux_g;

    int pre_pdf_bin_tag1, pdf_bin_tag1;
    int pre_pdf_bin_tag2, pdf_bin_tag2;
    int hist2d_mg_sigma_bin_mid, hist2d_mn_sigma_bin_mid, hist2d_total_tag;
    
    hist2d_mg_sigma_bin_mid = data_info->hist2d_mg_sigma_bin_num/2;
    hist2d_mn_sigma_bin_mid = data_info->hist2d_mn_sigma_bin_num/2;

    MY_FLOAT rotation_mat[6];
    int pair_count = 0;

    st = clock();

    initialize_arr(data_info->worker_sub_chi_sigma_tan, data_info->chi_sigma_theta_block_len, 0);
    initialize_arr(data_info->worker_sub_chi_sigma_cross,data_info->chi_sigma_theta_block_len, 0);
    initialize_arr(data_info->worker_sub_chi_g_tan, data_info->chi_g_theta_block_len, 0);
    initialize_arr(data_info->worker_sub_chi_g_cross, data_info->chi_g_theta_block_len, 0);
    initialize_arr(data_info->worker_sub_signal_count, data_info->sub_signal_count_len, 0);


    ggl_read_len_exp(data_info, len_expo_label);

    ggl_find_src_needed(data_info, len_expo_label);
    
    // std::cout<<"start... "<<std::endl;
    // i=0;
    // for(bkg=0; bkg<data_info->src_expo_num;bkg++)
    // {
    //     if(data_info->src_expo_needed_tag[bkg] == 1){i++;}
    // }
    // std::cout<<i<<" expos"<<std::endl;

    for(bkg=0; bkg<data_info->src_expo_num; bkg++)
    {   
        if(data_info->src_expo_needed_tag[bkg]< 1){continue;}

        ggl_read_src_exp(data_info, bkg);
        
        // std::cout<<data_info->rank<<" "<<bkg<<" start... "<<std::endl;

        for(ifg=0; ifg<data_info->len_data_row[len_expo_label]; ifg++)
        {    
            ifg_row = ifg*data_info->len_data_col;

            len_ra = data_info->len_expo_data[ifg_row + data_info->len_ra_col];
            len_dec = data_info->len_expo_data[ifg_row + data_info->len_dec_col];
            len_cos_dec = data_info->len_expo_data[ifg_row + data_info->len_cos_dec_col];

            len_z = data_info->len_expo_data[ifg_row + data_info->len_z_col];
            len_z_dz = len_z + data_info->back_dz;

            len_dist = data_info->len_expo_data[ifg_row + data_info->len_com_dist_col];

            // stacking in physical or comoving coordinate
#ifdef GGL_PROP_DIST_STACK
            coeff = data_info->crit_coeff/len_dist*(1+len_z);
#else
            coeff = data_info->crit_coeff/len_dist/(1+len_z);
#endif            
            for(ibkg=0; ibkg<data_info->src_data_row[bkg]; ibkg++)
            {   
                ibkg_row = ibkg*data_info->src_data_col;
                src_z = data_info->src_expo_data[ibkg_row + data_info->src_z_col];
                src_z_err = src_z - data_info->src_expo_data[ibkg_row + data_info->src_zerr_col];
                if(len_z >= src_z_err or len_z_dz >= src_z){continue;}

                src_ra = data_info->src_expo_data[ibkg_row + data_info->src_ra_col];
                src_dec = data_info->src_expo_data[ibkg_row + data_info->src_dec_col];
                src_dist = data_info->src_expo_data[ibkg_row + data_info->src_com_dist_col];

                sigma_crit = coeff*src_dist/(src_dist - len_dist);

                // dra = (len_ra - src_ra)*len_cos_dec;
                // ddec = len_dec - src_dec;
                // sep_theta = sqrt(dra*dra + ddec*ddec)*DEG2RAD;
                separation_angle_2(len_ra, len_dec, src_ra, src_dec, sep_theta);
                // sprintf
                // std::cout<<sep_theta_<<" "<<sep_theta<<std::endl;
#ifdef GGL_PROP_DIST_STACK
                sep_dist = sep_theta*data_info->len_expo_data[ifg_row + data_info->len_com_dist_col]/(1+len_dist);
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
                    // std::cout<<bkg<<" "<<sep_bin_tag<<std::endl;
                    pair_count ++;
                    // if(sep_bin_tag == 3){std::cout<<ifg<<" "<<ibkg<<std::endl;}
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
                    // std::cout<<"rot "<<len_ra<<" "<<len_dec<<" "<<len_cos_dec<<" "<<src_ra<<" "<<src_dec<<std::endl;
                    // show_arr(rotation_mat,1,6);
                    // std::cout<<"G "<<src_mg1<<" "<<src_mg2<<" "<<src_mn<<" "<<src_mu<<" "<<src_mv<<std::endl;
                    // std::cout<<"Gr "<<src_mg1_rot<<" "<<src_mg2_rot<<" "<<src_mn<<" "<<src_mu_rot<<std::endl;
                    // std::cout<<std::endl;
                    
                    temp_mnut = src_mn + src_mu_rot;
                    temp_mnux = src_mn - src_mu_rot;
                    // pdf
#ifdef GGL_GAMMA_T 

                    pre_pdf_bin_tag1 = 0;
                    pre_pdf_bin_tag2 = 0;

                    temp_mnut_g = temp_mnut;
                    temp_mnux_g = temp_mnux;

                    chi_gt_pos = sep_bin_tag*data_info->chi_g_theta_block_len_sub;
                    // std::cout<<chi_gt_pos<<" "<<sep_bin_tag<<" "<<data_info->chi_g_theta_block_len_sub<<std::endl;
                    for(i=0; i<data_info->pdf_gt_num; i++)
                    { 
                        // \gamma_t & \gamma_x
                        chi_gt_pos_i = chi_gt_pos + i*data_info->mg_gt_bin_num;
                        temp_mgt = src_mg1_rot - data_info->gt_guess[i]*temp_mnut_g;
                        temp_mgx = src_mg2_rot - data_info->gt_guess[i]*temp_mnux_g;
                        ggl_fast_hist(data_info->mg_gt_bin, data_info->mg_gt_bin_num, temp_mgt, pre_pdf_bin_tag1, pdf_bin_tag1);
                        ggl_fast_hist(data_info->mg_gt_bin, data_info->mg_gt_bin_num, temp_mgx, pre_pdf_bin_tag2, pdf_bin_tag2);
                        
                        data_info->worker_sub_chi_g_tan[chi_gt_pos_i + pdf_bin_tag1] +=1;
                        data_info->worker_sub_chi_g_cross[chi_gt_pos_i + pdf_bin_tag2] +=1;
                        
                        // std::cout<<sep_bin_tag<<" "<<i<<" "<<chi_gt_pos<<" "<<chi_gt_pos_i<<std::endl;
                        // std::cout<<temp_mgt<<" "<<src_mn<<" "<<src_mu_rot<<" "<<data_info->gt_guess[i]<<" "<<pre_pdf_bin_tag1<<" "<<pdf_bin_tag1<<std::endl;
                        // show_arr(data_info->mg_gt_bin,1,data_info->mg_gt_bin_num+1);
                        // std::cout<<std::endl;

                        pre_pdf_bin_tag1 = pdf_bin_tag1;
                        pre_pdf_bin_tag2 = pdf_bin_tag2;
                    }
#endif 

#ifdef GGL_DELTA_SIGMA

                    src_mg1_rot *= sigma_crit;
                    src_mg2_rot *= sigma_crit;

                    ggl_fast_hist(data_info->hist2d_mg_sigma_bin, data_info->hist2d_mg_sigma_bin_num, src_mg1_rot, hist2d_mg_sigma_bin_mid, pdf_bin_tag1);
                    ggl_fast_hist(data_info->hist2d_mn_sigma_bin, data_info->hist2d_mn_sigma_bin_num, temp_mnut, hist2d_mn_sigma_bin_mid, pdf_bin_tag2);
                    // // x: mgt*sigma_crit, y: N+U
                    hist2d_total_tag = sep_bin_tag*data_info->hist2d_len + pdf_bin_tag2*data_info->hist2d_mg_sigma_bin_num + pdf_bin_tag1;
                    
                    // std::cout<<data_info->rank<<" "<<hist2d_total_tag<<"("<<data_info->hist2d_total_len<<") "<<len_expo_label<<" "<<bkg<<std::endl;
                    // if(hist2d_total_tag > data_info->hist2d_total_len - 1 or hist2d_total_tag < 0)
                    // {
                    //     std::cout<<data_info->rank<<" "<<len_expo_label<<" "<<bkg<<" "<<pdf_bin_tag1<<" "<<pdf_bin_tag2<<" "<<hist2d_total_tag<<std::endl;
                    // }

                    data_info->hist2d_count[hist2d_total_tag] += 1;
                    data_info->hist2d_x[hist2d_total_tag] += src_mg1_rot;
                    data_info->hist2d_y[hist2d_total_tag] += temp_mnut;


                    chi_sigma_pos = sep_bin_tag*data_info->chi_sigma_theta_block_len_sub;

                    pre_pdf_bin_tag1 = 0;
                    pre_pdf_bin_tag2 = 0;

                    for(i=0; i<data_info->pdf_sigma_num; i++)
                    { 
                        // \Delta\Sigma(R) & \Delta\Sigma(R)_x
                        chi_sigma_pos_i = chi_sigma_pos + i*data_info->mg_sigma_bin_num;
                        temp_mgt = src_mg1_rot - data_info->delta_sigma_guess[i]*temp_mnut;
                        temp_mgx = src_mg2_rot - data_info->delta_sigma_guess[i]*temp_mnux;

                        // if(temp_mgt < data_info->mg_sigma_bin[0] or temp_mgt > data_info->mg_sigma_bin[data_info->mg_sigma_bin_num])
                        // {
                        //     std::cout<<"gt "<<temp_mgt<<" "<<src_mg1_rot<<" "<<data_info->delta_sigma_guess[i]<<" "<<temp_mnut<<std::endl;
                        // }

                        // if(temp_mgx < data_info->mg_sigma_bin[0] or temp_mgx > data_info->mg_sigma_bin[data_info->mg_sigma_bin_num])
                        // {
                        //     std::cout<<"gx "<<temp_mgx<<" "<<src_mg2_rot<<" "<<data_info->delta_sigma_guess[i]<<" "<<temp_mnux<<std::endl;
                        // }

                        ggl_fast_hist(data_info->mg_sigma_bin, data_info->mg_sigma_bin_num, temp_mgt, pre_pdf_bin_tag1, pdf_bin_tag1);
                        ggl_fast_hist(data_info->mg_sigma_bin, data_info->mg_sigma_bin_num, temp_mgx, pre_pdf_bin_tag2, pdf_bin_tag2);


                        // std::cout<<pre_pdf_bin_tag1<<" "<<pdf_bin_tag1<<" "<<data_info->mg_sigma_bin[pdf_bin_tag1]<<" "<<temp_mgt<<" "<<data_info->mg_sigma_bin[pdf_bin_tag1+1]<<std::endl;
                        // std::cout<<pre_pdf_bin_tag2<<" "<<pdf_bin_tag2<<" "<<data_info->mg_sigma_bin[pdf_bin_tag2]<<" "<<temp_mgx<<" "<<data_info->mg_sigma_bin[pdf_bin_tag2+1]<<std::endl;

                        data_info->worker_sub_chi_sigma_tan[chi_sigma_pos_i + pdf_bin_tag1] +=1;
                        data_info->worker_sub_chi_sigma_cross[chi_sigma_pos_i + pdf_bin_tag2] +=1;

                        //show_arr(data_info->mg_sigma_bin, 1, data_info->mg_sigma_bin_num+1);
                        // std::cout<<sigma_crit<<" "<<data_info->mg_sigma_bin[pdf_bin_tag1]<<" "<<temp_mgt<<" "<<data_info->mg_sigma_bin[pdf_bin_tag1]<<" "<<pre_pdf_bin_tag1<<" "<<pdf_bin_tag1<<std::endl;

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
        if(data_info->jack_num > 2)
        {if(i == data_info->len_expo_jackid[len_expo_label]){continue;}}
#ifdef GGL_DELTA_SIGMA
        for(j=0; j<data_info->chi_sigma_theta_block_len; j++)
        {
            data_info->worker_total_chi_sigma_tan[i*data_info->chi_sigma_theta_block_len + j] += data_info->worker_sub_chi_sigma_tan[j];
            data_info->worker_total_chi_sigma_cross[i*data_info->chi_sigma_theta_block_len + j] += data_info->worker_sub_chi_sigma_cross[j];
        }
#endif

#ifdef GGL_GAMMA_T
        for(j=0; j<data_info->chi_g_theta_block_len; j++)
        {
            data_info->worker_total_chi_g_tan[i*data_info->chi_g_theta_block_len + j] += data_info->worker_sub_chi_g_tan[j];
            data_info->worker_total_chi_g_cross[i*data_info->chi_g_theta_block_len + j] += data_info->worker_sub_chi_g_cross[j];
        }
#endif
        for(j=0; j<data_info->sub_signal_count_len; j++)
        {
            data_info->worker_total_signal_count[i*data_info->sub_signal_count_len + j] += data_info->worker_sub_signal_count[j];
        }
    }

    ed = clock();
    char times[100];
    sprintf(times,"worker %d. %dth Lens file, %dth jack. Finished in %.2f sec. %d Lenses, %d pairs",
            data_info->rank, len_expo_label, data_info->len_expo_jackid[len_expo_label], (ed-st)/CLOCKS_PER_SEC, data_info->len_data_row[len_expo_label], pair_count);
    std::cout<<times<<std::endl;
}

void ggl_cache(ggl_data_info *data_info)
{
    int i, j;
    MPI_Status status;
    char set_name[50];

    // sprintf(data_info->ggl_log_inform,"%d Start caching\n", data_info->rank);
    // if(data_info->rank >= 0) {std::cout<<data_info->ggl_log_inform;} 

    sprintf(set_name,"/data");
    if (data_info->rank > 0)
    {
        sprintf(data_info->ggl_result_path,"%s/result/chi_theta_cache_%d.hdf5", data_info->ggl_total_path, data_info->rank);
        write_h5(data_info->ggl_result_path, set_name, data_info->worker_total_signal_count, data_info->jack_num+1, data_info->sub_signal_count_len, true);
        // std::cout<<data_info->ggl_result_path<<std::endl;
    }
    else
    {        
        sprintf(data_info->ggl_result_path,"%s/result/chi_theta_cache.hdf5", data_info->ggl_total_path);
        write_h5(data_info->ggl_result_path, set_name, data_info->total_signal_count, data_info->jack_num+1,data_info->sub_signal_count_len, true);
        // std::cout<<data_info->ggl_result_path<<std::endl;
    }

#ifdef GGL_DELTA_SIGMA
    if (data_info->rank > 0)
    {
        sprintf(data_info->ggl_result_path,"%s/result/chi_sigma_cache_%d.hdf5", data_info->ggl_total_path, data_info->rank);

        sprintf(set_name,"/t");
        write_h5(data_info->ggl_result_path, set_name, data_info->worker_total_chi_sigma_tan,
                 data_info->jack_num+1, data_info->chi_sigma_theta_block_len, true);

        sprintf(set_name,"/x");
        write_h5(data_info->ggl_result_path, set_name, data_info->worker_total_chi_sigma_cross,
                 data_info->jack_num+1, data_info->chi_sigma_theta_block_len, false);
    }
    else
    {
        sprintf(data_info->ggl_result_path,"%s/result/chi_sigma_cache.hdf5", data_info->ggl_total_path);

        sprintf(set_name,"/t");
        write_h5(data_info->ggl_result_path, set_name, data_info->total_chi_sigma_tan,
                 data_info->jack_num+1, data_info->chi_sigma_theta_block_len, true);

        sprintf(set_name,"/x");
        write_h5(data_info->ggl_result_path, set_name, data_info->total_chi_sigma_cross,
                 data_info->jack_num+1, data_info->chi_sigma_theta_block_len, false);
    }
#endif

#ifdef GGL_GAMMA_T
if (data_info->rank > 0)
    {
        sprintf(data_info->ggl_result_path,"%s/result/chi_g_cache_%d.hdf5", data_info->ggl_total_path, data_info->rank);

        sprintf(set_name,"/t");
        write_h5(data_info->ggl_result_path, set_name, data_info->worker_total_chi_g_tan,
                 data_info->jack_num+1, data_info->chi_g_theta_block_len, true);
        
        sprintf(set_name,"/x");
        write_h5(data_info->ggl_result_path, set_name, data_info->worker_total_chi_g_cross,
                 data_info->jack_num+1, data_info->chi_g_theta_block_len, false);
    }
    else
    {               
        sprintf(data_info->ggl_result_path,"%s/result/chi_g_cache.hdf5", data_info->ggl_total_path);

        sprintf(set_name,"/t");
        write_h5(data_info->ggl_result_path, set_name, data_info->total_chi_g_tan,
                 data_info->jack_num+1, data_info->chi_g_theta_block_len, true);

        sprintf(set_name,"/x");
        write_h5(data_info->ggl_result_path, set_name, data_info->total_chi_g_cross,
                 data_info->jack_num+1, data_info->chi_g_theta_block_len, false);
    }
#endif

//     sprintf(data_info->ggl_log_inform,"%d Finish caching\n", data_info->rank);
//     if(data_info->rank >= 0) {std::cout<<data_info->ggl_log_inform;} 
}


void ggl_collect_chi(ggl_data_info *data_info)
{
    int i, j;
    MPI_Status status;

    // sprintf(data_info->ggl_log_inform,"%d Start collecting data\n", data_info->rank);
    // if(data_info->rank >= 0) {std::cout<<data_info->ggl_log_inform;} 

    if (data_info->rank > 0)
    {MPI_Send(data_info->worker_total_signal_count, data_info->total_signal_count_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD);}
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->worker_total_signal_count, data_info->total_signal_count_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
        
            for(j=0;j<data_info->total_signal_count_len;j++)
            {data_info->total_signal_count[j] += data_info->worker_total_signal_count[j];}
        }
    }

#ifdef GGL_DELTA_SIGMA
    if (data_info->rank > 0)
    {MPI_Send(data_info->worker_total_chi_sigma_tan, data_info->chi_sigma_jack_block_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD);}
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->worker_total_chi_sigma_tan, data_info->chi_sigma_jack_block_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            for(j=0;j<data_info->chi_sigma_jack_block_len;j++)
            {data_info->total_chi_sigma_tan[j] += data_info->worker_total_chi_sigma_tan[j];}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    if (data_info->rank > 0)
    {MPI_Send(data_info->worker_total_chi_sigma_cross, data_info->chi_sigma_jack_block_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD);}
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->worker_total_chi_sigma_cross, data_info->chi_sigma_jack_block_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            for(j=0;j<data_info->chi_sigma_jack_block_len;j++)
            {data_info->total_chi_sigma_cross[j] += data_info->worker_total_chi_sigma_cross[j];}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    if (data_info->rank > 0)
    {MPI_Send(data_info->hist2d_count, data_info->hist2d_total_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD);}
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->hist2d_count, data_info->hist2d_total_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);

            for(j=0;j<data_info->hist2d_total_len;j++)
            {data_info->hist2d_count_total[j] += data_info->hist2d_count[j];}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (data_info->rank > 0)
    {MPI_Send(data_info->hist2d_x, data_info->hist2d_total_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD);}
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->hist2d_x, data_info->hist2d_total_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);

            for(j=0;j<data_info->hist2d_total_len;j++)
            {data_info->hist2d_x_total[j] += data_info->hist2d_x[j];}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (data_info->rank > 0)
    {MPI_Send(data_info->hist2d_y, data_info->hist2d_total_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD);}
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->hist2d_y, data_info->hist2d_total_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);

            for(j=0;j<data_info->hist2d_total_len;j++)
            {data_info->hist2d_y_total[j] += data_info->hist2d_y[j];}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef GGL_GAMMA_T
    if (data_info->rank > 0)
    { MPI_Send(data_info->worker_total_chi_g_tan, data_info->chi_g_jack_block_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD);}
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->worker_total_chi_g_tan, data_info->chi_g_jack_block_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
        
            for(j=0;j<data_info->chi_g_jack_block_len;j++)
            {data_info->total_chi_g_tan[j] += data_info->worker_total_chi_g_tan[j];}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

if (data_info->rank > 0)
    {MPI_Send(data_info->worker_total_chi_g_cross, data_info->chi_g_jack_block_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD); }
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->worker_total_chi_g_cross, data_info->chi_g_jack_block_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
        
            for(j=0;j<data_info->chi_g_jack_block_len;j++)
            {data_info->total_chi_g_cross[j] += data_info->worker_total_chi_g_cross[j];}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // sprintf(data_info->ggl_log_inform,"%d Finish collecting data\n", data_info->rank);
    // if(data_info->rank >= 0) {std::cout<<data_info->ggl_log_inform;} 
}



void ggl_cal_signals(ggl_data_info * data_info)
{   
    sprintf(data_info->ggl_log_inform,"\n========================== start calculate ==========================\n");
    std::cout<<data_info->ggl_log_inform;

    char set_name[50];
    int i, j, st_c, st_j;

    double *theta = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
    double *radius = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
    double *count = new double[(data_info->jack_num+1)*data_info->signal_pts_num];

    double *signals = new double[data_info->signal_pts_num];
    double *signals_err = new double[data_info->signal_pts_num];

    double *chisq_all = new double[data_info->signal_pts_num*data_info->pdf_sigma_num];
    double *chisq_fit_coeff = new double[data_info->signal_pts_num*3];

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

        if(i == data_info->jack_num)
        {   
            std::cout<<"\n========================== Count ==========================\n";
            for(j=0;j<data_info->signal_pts_num;j++)
            {std::cout<<count[st_c+j]<<" ";}
            std::cout<<std::endl;
            std::cout<<"\n========================== Theta [armin] ==========================\n";
            for(j=0;j<data_info->signal_pts_num;j++)
            {std::cout<<theta[st_c+j]/count[st_c+j]*60<<" ";}
            std::cout<<std::endl;
            std::cout<<"\n========================== Radius [Mpc/h] ==========================\n";
            for(j=0;j<data_info->signal_pts_num;j++)
            {std::cout<<radius[st_c+j]/count[st_c+j]<<" ";}
            std::cout<<std::endl;
        }           
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

    sprintf(set_name,"/delta_t_guess");
    write_h5(data_info->ggl_result_path, set_name, data_info->delta_sigma_guess,
            data_info->pdf_sigma_num, 1, false);

    sprintf(set_name,"/g_t_guess");
    write_h5(data_info->ggl_result_path, set_name, data_info->gt_guess,
            data_info->pdf_gt_num, 1, false);
          
            
    sprintf(data_info->ggl_log_inform,"========================== Finish calculating theta & radius ==========================\n");
    std::cout<<data_info->ggl_log_inform;



#ifdef GGL_DELTA_SIGMA

    sprintf(data_info->ggl_log_inform,"\n========================== Start calculating GGL_DELTA_SIGMA ==========================\n");
    std::cout<<data_info->ggl_log_inform;

    double *temp_sigma = new double[data_info->chi_sigma_theta_block_len];
    double *delta_sigma_tan = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
    double *delta_sigma_tan_err = new double[(data_info->jack_num+1)*data_info->signal_pts_num];

    double *delta_sigma_cross = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
    double *delta_sigma_cross_err = new double[(data_info->jack_num+1)*data_info->signal_pts_num];


    for(i=0; i<data_info->jack_num+1; i++)
    {   
        st_j = i*data_info->signal_pts_num;

        // delta_sigma_t
        st_c = i*data_info->chi_sigma_theta_block_len;
        for(j=0; j<data_info->chi_sigma_theta_block_len; j++)
        { temp_sigma[j] = data_info->total_chi_sigma_tan[st_c+j];}

        // show_arr(temp_count, 1, data_info->chi_signal_block_len);
        ggl_pdf_signals(temp_sigma, data_info->delta_sigma_guess, data_info->pdf_sigma_num, 
                        data_info->mg_sigma_bin_num, data_info->signal_pts_num, signals, signals_err, chisq_all, chisq_fit_coeff);
        
        for(j=0; j<data_info->signal_pts_num; j++)
        {
            delta_sigma_tan[st_j + j] = signals[j];
            delta_sigma_tan_err[st_j + j] = signals_err[j];
        }

        // just for the total sample, 
        // save the chi^2 of the PDF_SYM process and the coefficients
        if(i == data_info->jack_num)
        {
            sprintf(set_name,"/Gt_pdf");
            write_h5(data_info->ggl_result_path, set_name, temp_sigma,
                    data_info->signal_pts_num,data_info->chi_sigma_theta_block_len_sub, false);

            sprintf(set_name,"/delta_t_chisq");
            write_h5(data_info->ggl_result_path, set_name, chisq_all,
                    data_info->signal_pts_num,data_info->pdf_sigma_num, false);

            sprintf(set_name,"/delta_t_chisq_coeff");
            write_h5(data_info->ggl_result_path, set_name, chisq_fit_coeff,
                    data_info->signal_pts_num, 3, false);
        }


        // delta_sigma_x
        for(j=0;j<data_info->chi_sigma_theta_block_len;j++)
        { temp_sigma[j] = data_info->total_chi_sigma_cross[st_c+j];}

        ggl_pdf_signals(temp_sigma, data_info->delta_sigma_guess, data_info->pdf_sigma_num, 
                        data_info->mg_sigma_bin_num, data_info->signal_pts_num, signals, signals_err,chisq_all, chisq_fit_coeff);
        
        for(j=0; j<data_info->signal_pts_num; j++)
        {
            delta_sigma_cross[st_j + j] = signals[j];
            delta_sigma_cross_err[st_j + j] = signals_err[j];
        }

        if(i == data_info->jack_num)
        {   
            std::cout<<"Delta Sigma_t\n";
            for(j=0;j<data_info->signal_pts_num;j++)
            {std::cout<<delta_sigma_tan[j]<<"("<<delta_sigma_tan_err[j]<<") ";}
            std::cout<<std::endl;

            std::cout<<"Delta Sigma_x\n";
            for(j=0;j<data_info->signal_pts_num;j++)
            {std::cout<<delta_sigma_cross[j]<<"("<<delta_sigma_cross_err[j]<<") ";}
            std::cout<<std::endl;
        } 
    }

    sprintf(set_name,"/delta_sigma_t");
    write_h5(data_info->ggl_result_path, set_name, delta_sigma_tan,
            data_info->jack_num+1,data_info->signal_pts_num, false);
    sprintf(set_name,"/delta_sigma_t_err");
    write_h5(data_info->ggl_result_path, set_name, delta_sigma_tan_err,
            data_info->jack_num+1,data_info->signal_pts_num, false);
    sprintf(set_name,"/delta_sigma_x");
    write_h5(data_info->ggl_result_path, set_name, delta_sigma_cross,
            data_info->jack_num+1,data_info->signal_pts_num, false);
    sprintf(set_name,"/delta_sigma_x_err");
    write_h5(data_info->ggl_result_path, set_name, delta_sigma_cross_err,
            data_info->jack_num+1,data_info->signal_pts_num, false);

    sprintf(set_name,"/mg_sigma_bin");
    write_h5(data_info->ggl_result_path, set_name, data_info->mg_sigma_bin,
            data_info->mg_sigma_bin_num+1, 1, false);
            
    // save the 2d hist for more calculations later
    for(i=0; i<data_info->signal_pts_num; i++)
    {   
        sprintf(set_name,"/delta_sigma_t_grid2d_count_%d", i);
        write_h5(data_info->ggl_result_path, set_name, &data_info->hist2d_count_total[i*data_info->hist2d_len], 
        data_info->hist2d_mn_sigma_bin_num, data_info->hist2d_mg_sigma_bin_num, false);

        sprintf(set_name,"/delta_sigma_t_grid2d_x_%d", i);
        write_h5(data_info->ggl_result_path, set_name, &data_info->hist2d_x_total[i*data_info->hist2d_len], 
        data_info->hist2d_mn_sigma_bin_num, data_info->hist2d_mg_sigma_bin_num, false);

        sprintf(set_name,"/delta_sigma_t_grid2d_y_%d", i);
        write_h5(data_info->ggl_result_path, set_name, &data_info->hist2d_y_total[i*data_info->hist2d_len], 
        data_info->hist2d_mn_sigma_bin_num, data_info->hist2d_mg_sigma_bin_num, false);
    }

    delete[] temp_sigma;
    delete[] delta_sigma_tan;
    delete[] delta_sigma_tan_err;
    delete[] delta_sigma_cross;
    delete[] delta_sigma_cross_err;

    sprintf(data_info->ggl_log_inform,"========================== Finish calculating GGL_DELTA_SIGMA ==========================\n");
    std::cout<<data_info->ggl_log_inform;
#endif

#ifdef GGL_GAMMA_T
    double *temp_g = new double[data_info->chi_g_theta_block_len];
    double *gt = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
    double *gt_err = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
    double *gx = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
    double *gx_err = new double[(data_info->jack_num+1)*data_info->signal_pts_num];


    sprintf(data_info->ggl_log_inform,"\n========================== Start calculating GGL_GAMMA_T ==========================\n");
    std::cout<<data_info->ggl_log_inform;
    for(i=0; i<data_info->jack_num+1; i++)
    {   
        st_j = i*data_info->signal_pts_num;

        // g_t
        st_c = i*data_info->chi_g_theta_block_len;

        for(j=0; j<data_info->chi_g_theta_block_len; j++)
        { temp_g[j] = data_info->total_chi_g_tan[st_c+j];}

        ggl_pdf_signals(temp_g, data_info->gt_guess, data_info->pdf_gt_num, 
                        data_info->mg_gt_bin_num, data_info->signal_pts_num, signals, signals_err, chisq_all, chisq_fit_coeff);

        for(j=0; j<data_info->signal_pts_num; j++)
        {
            gt[st_j + j] = signals[j]; 
            gt_err[st_j + j] = signals_err[j];
            // std::cout<<signals[j]<<" ";
        }
        // std::cout<<std::endl;
        // just for the total sample, 
        // save the chi^2 of the PDF_SYM process and the coefficients
        if(i == data_info->jack_num)
        {
            sprintf(set_name,"/g_t_chisq");
            write_h5(data_info->ggl_result_path, set_name, chisq_all,
                    data_info->signal_pts_num,data_info->pdf_gt_num, false);

            sprintf(set_name,"/g_t_chisq_coeff");
            write_h5(data_info->ggl_result_path, set_name, chisq_fit_coeff,
                    data_info->signal_pts_num, 3, false);
        }


        // g_x
        st_c = i*data_info->chi_g_theta_block_len;

        for(j=0; j<data_info->chi_g_theta_block_len; j++)
        { temp_g[j] = data_info->total_chi_g_cross[st_c+j];}

        ggl_pdf_signals(temp_g, data_info->gt_guess, data_info->pdf_gt_num, 
                        data_info->mg_gt_bin_num, data_info->signal_pts_num, signals, signals_err, chisq_all, chisq_fit_coeff);
        
        for(j=0; j<data_info->signal_pts_num; j++)
        {
            gx[st_j + j] = signals[j]; 
            gx_err[st_j + j] = signals_err[j]; 
        }

        if(i == data_info->jack_num)
        {   
            std::cout<<"gt\n";
            for(j=0;j<data_info->signal_pts_num;j++)
            {std::cout<<gt[j]<<"("<<gt_err[j]<<") ";}
            std::cout<<std::endl;
            std::cout<<"gx\n";
            for(j=0;j<data_info->signal_pts_num;j++)
            {std::cout<<gx[j]<<"("<<gx_err[j]<<") ";}
            std::cout<<std::endl;
        } 
    }

    sprintf(set_name,"/g_t");
    write_h5(data_info->ggl_result_path, set_name, gt,
            data_info->jack_num+1,data_info->signal_pts_num, false);
    sprintf(set_name,"/g_t_err");
    write_h5(data_info->ggl_result_path, set_name, gt_err,
            data_info->jack_num+1,data_info->signal_pts_num, false);
    sprintf(set_name,"/g_x");
    write_h5(data_info->ggl_result_path, set_name, gx,
            data_info->jack_num+1,data_info->signal_pts_num, false);
    sprintf(set_name,"/g_x_err");
    write_h5(data_info->ggl_result_path, set_name, gx_err,
            data_info->jack_num+1,data_info->signal_pts_num, false);

    sprintf(set_name,"/mg_gt_bin");
    write_h5(data_info->ggl_result_path, set_name, data_info->mg_gt_bin,
            data_info->mg_gt_bin_num+1, 1, false);

    sprintf(data_info->ggl_log_inform,"========================== Finish calculating GGL_GAMMA_T ==========================\n");
    std::cout<<data_info->ggl_log_inform;

    delete[] temp_g;
    delete[] gt;
    delete[] gt_err;
    delete[] gx;
    delete[] gx_err;

#endif

    sprintf(data_info->ggl_log_inform,"finish calculate and save result.\n");
    std::cout<<data_info->ggl_log_inform;  

    delete[] theta;
    delete[] radius;
    delete[] count;
    delete[] signals;
    delete[] signals_err;
    delete[] chisq_fit_coeff;
    delete[] chisq_all;
}

void ggl_pdf_signals(double *chi_count, double*pdf_signal_guess, int pdf_guess_num, int mg_bin_num, int signal_pts_num, double *signal, double *signal_err, double *chisq_all, double *chisq_fit_coeff)
{
    int i, j, k, st;
    double chisq_i;
    double signal_i, signal_err_i;
    double *temp = new double[mg_bin_num];
    double *chisq = new double[pdf_guess_num];
    double fit_coeff[3];

    for(i=0; i<signal_pts_num; i++)
    {
        for(j=0; j<pdf_guess_num; j++)
        {   
            st = i*pdf_guess_num*mg_bin_num + j*mg_bin_num;
            for(k=0; k<mg_bin_num; k++)
            {
                temp[k] = chi_count[st+k];
            }
            // show_arr(temp, 1, mg_bin_num);
            cal_chisq_1d(temp, mg_bin_num, chisq_i);
            chisq[j] = chisq_i;
            
            chisq_all[i*pdf_guess_num + j] = chisq_i;
        }
        // show_arr(chisq, 1, pdf_guess_num);
        fit_shear(pdf_signal_guess, chisq, pdf_guess_num, signal_i, signal_err_i, chisq_i, fit_coeff,1, 100);
        signal[i] = signal_i;
        signal_err[i] = signal_err_i;

        chisq_fit_coeff[3*i] = fit_coeff[0];
        chisq_fit_coeff[3*i + 1] = fit_coeff[1];
        chisq_fit_coeff[3*i + 2] = fit_coeff[2];

    }
    delete[] temp;
    delete[] chisq;
}
