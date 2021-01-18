#include<science_cal_lib.h>


/////////////////////////  GGL part  /////////////////////////////////////////
void ggl_initialize(ggl_data_info *data_info)
{   
    int i;

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

    data_info->back_dz = 0.05;

    data_info->pos_inform_num = 4;

    // when read a len exposure, it will calculate the separation between the 
    // len and src and labelled the src exposures.
    data_info->src_expo_needed_tag = new int[data_info->src_expo_num];
    initialize_arr(data_info->src_expo_needed_tag, data_info->src_expo_num,-1);
    

    // read G bins
    sprintf(data_info->ggl_pdf_inform_path,"%s/cata/pdf_inform.hdf5", data_info->ggl_total_path);
    sprintf(data_info->set_name,"/mg_bin");
    read_h5_datasize(data_info->ggl_pdf_inform_path, data_info->set_name,data_info->mg_bin_num);
    data_info->mg_bin = new MY_FLOAT[data_info->mg_bin_num];
    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->mg_bin);
    data_info->mg_bin_num = data_info->mg_bin_num - 1;


    // read signal guesses
    sprintf(data_info->set_name,"/delta_sigma_guess");
    read_h5_datasize(data_info->ggl_pdf_inform_path, data_info->set_name,data_info->pdf_guess_num);
    data_info->gt_guess = new double[data_info->pdf_guess_num];
    data_info->delta_sigma_guess = new double[data_info->pdf_guess_num];

    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->delta_sigma_guess);
    sprintf(data_info->set_name,"/g_guess");
    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->gt_guess);

    sprintf(data_info->set_name,"/separation_bin");
    read_h5_datasize(data_info->ggl_pdf_inform_path, data_info->set_name,data_info->sep_bin_num);
    data_info->separation_bin = new MY_FLOAT[data_info->sep_bin_num];
    read_h5(data_info->ggl_pdf_inform_path, data_info->set_name, data_info->separation_bin);
    data_info->signal_pts_num = data_info->sep_bin_num - 1;
    data_info->sep_bin_num = data_info->sep_bin_num - 1;



    data_info->chi_theta_block_len = data_info->pdf_guess_num*data_info->mg_bin_num;
    data_info->chi_signal_block_len = data_info->signal_pts_num*data_info->chi_theta_block_len;
    data_info->chi_jack_block_len = data_info->chi_signal_block_len*4 + data_info->signal_pts_num*3;
    data_info->total_chi_count_len = data_info->chi_jack_block_len*(data_info->jack_num+1);
    // all pair data will be stacked into the last block in total_chi_count for the estimation of signals
    data_info->worker_total_chi_count = new double[data_info->total_chi_count_len];
    data_info->worker_sub_chi_count = new double[data_info->chi_jack_block_len];
    if(data_info->rank == 0)
    {
        data_info->total_chi_count = new double[data_info->total_chi_count_len];
        initialize_arr(data_info->total_chi_count, data_info->total_chi_count_len, 0);
    }
    initialize_arr(data_info->worker_total_chi_count, data_info->total_chi_count_len, 0);
    initialize_arr(data_info->worker_sub_chi_count, data_info->chi_jack_block_len, 0);
    

    // read fore/background exposure informs
    sprintf(data_info->ggl_foreground_inform_path, "%s/foreground/foreground_list.dat", data_info->ggl_total_path);
    sprintf(data_info->ggl_background_inform_path, "%s/background/background_list.dat", data_info->ggl_total_path);
    line_count(data_info->ggl_foreground_inform_path, data_info->len_expo_num);
    line_count(data_info->ggl_background_inform_path, data_info->src_expo_num);

    // for task distribution
    data_info->len_expo_num_remain = data_info->len_expo_num;

    // data_info->len_nearest_dist = new MY_FLOAT[data_info->len_expo_num];


    for(i=0; i<data_info->len_expo_num; i++)
    {
         data_info->len_expo_path[i] = new char[450];
         data_info->len_expo_name[i] = new char[50];
    }
    data_info->len_data_row = new int[data_info->len_expo_num];
    data_info->len_expo_jackid = new int[data_info->len_expo_num];

    for(i=0; i<data_info->src_expo_num; i++)
    { 
        data_info->src_expo_path[i] = new char[450];
        data_info->src_expo_name[i] = new char[50];
        data_info->src_pos_informs[i] = new MY_FLOAT[data_info->pos_inform_num*data_info->src_expo_num];
    }
    data_info->src_data_row = new int[data_info->src_expo_num];

    ggl_read_len_list(data_info->ggl_foreground_inform_path, data_info);
    ggl_read_src_list(data_info->ggl_background_inform_path, data_info);

    if(data_info->rank == 0)
    {   
        char cal_inform[500];
        std::cout<<data_info->ggl_total_path<<std::endl;

        sprintf(cal_inform, "Foreground %d expos.  Background %d expos\n", data_info->len_expo_num, data_info->src_expo_num);
        std::cout<<cal_inform;

        sprintf(cal_inform,"G bins: %d.\n", data_info->mg_bin_num);
        std::cout<<cal_inform;
        show_arr(data_info->mg_bin,1,data_info->mg_bin_num+1);
        std::cout<<std::endl;

        sprintf(cal_inform,"PDF_guess_num: %d.\n", data_info->pdf_guess_num);
        std::cout<<cal_inform;
        show_arr(data_info->gt_guess,1,data_info->pdf_guess_num);
        std::cout<<std::endl;
        show_arr(data_info->delta_sigma_guess,1,data_info->pdf_guess_num);
        std::cout<<std::endl;

        sprintf(cal_inform,"Separation bins: %d", data_info->sep_bin_num);
        std::cout<<cal_inform;
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
        strs >> data_info->len_expo_name[line_count];
        strs >> data_info->len_data_row[line_count];
        strs >> data_info->len_nearest_dist[line_count];
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
        strs >> data_info->src_expo_name[line_count];
        strs >> data_info->src_data_row[line_count];
        strs >> data_info->src_pos_informs[data_info->pos_inform_num][0];
        strs >> data_info->src_pos_informs[data_info->pos_inform_num][1];
        strs >> data_info->src_pos_informs[data_info->pos_inform_num][2];
        strs >> data_info->src_pos_informs[data_info->pos_inform_num][3];

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
    
    size = data_info->len_data_col*data_info->len_data_row[len_expo_label];
    data_info->len_expo_data = new MY_FLOAT[size];
    read_h5(data_info->len_expo_path[len_expo_label], set_name, data_info->len_expo_data);
    
    data_info->len_expo_read_tag = 1;
}

void ggl_read_src_exp(ggl_data_info *data_info, int src_expo_label)
{
    if(data_info->src_expo_read_tag == 1)
    {delete[] data_info->src_expo_data;}

    char set_name[50];
    int size;
    sprintf(set_name,"/data");
    
    size = data_info->src_data_col*data_info->src_data_row[src_expo_label];
    data_info->src_expo_data = new MY_FLOAT[size];
    read_h5(data_info->src_expo_path[src_expo_label], set_name, data_info->src_expo_data);
    
    data_info->src_expo_read_tag = 1;
}

void ggl_find_src_needed(ggl_data_info *data_info, int len_expo_label)
{
    int ifg, ifg_row, bkg;
    MY_FLOAT max_sep_theta;
    MY_FLOAT dra, ddec, sep_theta;
    MY_FLOAT len_ra, len_dec, len_cos_dec;

    initialize_arr(data_info->src_expo_needed_tag, data_info->src_expo_num, -1);

    for(ifg=0; ifg<data_info->len_data_row[len_expo_label]; ifg++)
    {   
        ifg_row = ifg*data_info->len_data_col;
        len_ra = data_info->len_expo_data[ifg_row + data_info->len_ra_col];
        len_dec = data_info->len_expo_data[ifg_row + data_info->len_dec_col];
        len_cos_dec = data_info->len_expo_data[ifg_row + data_info->len_cos_dec_col];
        
        max_sep_theta = data_info->separation_bin[data_info->sep_bin_num+1]/
        data_info->len_expo_data[ifg_row + data_info->len_com_dist_col]*1.2;

        // if stack the signal in physical coordinate
#ifdef GGL_PROP_DIST_STACK
        max_sep_theta = 1.2*data_info->separation_bin[data_info->sep_bin_num+1]/
        data_info->len_expo_data[ifg_row + data_info->len_prop_dist_col];
#else
        max_sep_theta = 1.2*data_info->separation_bin[data_info->sep_bin_num+1]/
        data_info->len_expo_data[ifg_row + data_info->len_com_dist_col];
#endif

        for(bkg=0; bkg<data_info->src_expo_num; bkg++)
        {
            dra = (len_ra - data_info->src_pos_informs[bkg][0])*len_cos_dec;
            ddec = len_dec - data_info->src_pos_informs[bkg][1];
            sep_theta = (sqrt(dra*dra + ddec*ddec) - data_info->src_pos_informs[bkg][2])/180*Pi;

            if(sep_theta <= max_sep_theta){data_info->src_expo_needed_tag[bkg] = 1;}
        }
    }

}

void ggl_rotation_matrix(MY_FLOAT cent_ra, MY_FLOAT cent_dec, MY_FLOAT src_ra, MY_FLOAT src_dec, MY_FLOAT*rotation_matrix)
{
    MY_FLOAT sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta;
    MY_FLOAT dra, ddec, delta_radius;
    
    dra = src_ra - cent_ra;
    ddec = src_dec - cent_dec;
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
    int ir, sep_bin_tag, chi_pos, chi_pos_i;
    int i,j,k;
    double st, ed;

    int jack_id;

    MY_FLOAT len_ra, len_dec, src_ra, src_dec;
    MY_FLOAT len_z, len_dist, src_z, src_dist;
    MY_FLOAT dra, ddec, delta_radius;
    MY_FLOAT sep_dist, sep_theta;
    MY_FLOAT sigma_crit, coeff;
    MY_FLOAT src_mg1, src_mg2, src_mn, src_mv, src_mu;
    MY_FLOAT src_mg1_rot, src_mg2_rot, src_mu_rot, src_mv_rot;
    MY_FLOAT temp_mgt, temp_mgx, temp_mnut, temp_mnux,temp_mnut_g, temp_mnux_g;

    int pre_bin_tag_sig1, pre_bin_tag_sig2, bin_tag_sig1, bin_tag_sig2;
    int pre_bin_tag_shear1, pre_bin_tag_shear2, bin_tag_shear1, bin_tag_shear2;

    MY_FLOAT rotation_mat[6];


    st = clock();

    initialize_arr(data_info->worker_sub_chi_count, data_info->chi_jack_block_len, 0);
    
    ggl_read_len_exp(data_info, len_expo_label);
    
    for(bkg=0; bkg<data_info->src_expo_num; bkg++)
    {   
        if(data_info->src_expo_needed_tag[bkg]< 0){continue;}
        
        ggl_read_src_exp(data_info, bkg);

        for(ifg=0; ifg<data_info->len_data_row[len_expo_label]; ifg++)
        {    
            ifg_row = ifg*data_info->len_data_col;
            len_ra = data_info->len_expo_data[ifg_row + data_info->len_ra_col];
            len_dec = data_info->len_expo_data[ifg_row + data_info->len_dec_col];
            len_z = data_info->len_expo_data[ifg_row + data_info->len_z_col];
            len_dist = data_info->len_expo_data[ifg_row + data_info->len_prop_dist_col];

            // stacking in physical or comoving coordinate
#ifdef GGL_PROP_DIST_STACK
            coeff = 1./len_dist;
#else
            coeff = 1./len_dist/(1+len_z)/(1+len_z);
#endif            
            for(ibkg=0; ibkg<data_info->src_data_row[bkg]; ibkg++)
            {   
                ibkg_row = ibkg*data_info->src_data_col;
                src_z = data_info->src_expo_data[ibkg_row + data_info->src_z_col] - data_info->back_dz;
                if(len_z <= src_z){continue;}

                src_ra = data_info->src_expo_data[ibkg_row + data_info->src_ra_col];
                src_dec = data_info->src_expo_data[ibkg_row + data_info->src_dec_col];
                src_dist = data_info->src_expo_data[ibkg_row + data_info->src_prop_dist_col];

                sigma_crit = coeff*src_dist/(src_dist - len_dist);

                separation(len_ra, len_dec, src_ra, src_dec, sep_theta);

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
                
                if(sep_bin_tag > -1)
                {   

                    data_info->worker_sub_chi_count[data_info->chi_signal_block_len*4 + sep_bin_tag] += sep_theta;
                    data_info->worker_sub_chi_count[data_info->chi_signal_block_len*4 + data_info->signal_pts_num + sep_bin_tag] += sep_dist;
                    data_info->worker_sub_chi_count[data_info->chi_signal_block_len*4 + data_info->signal_pts_num + sep_bin_tag] += 1;

                    // rotation, sin_theta, cos_theta, sin_2theta, cos_2theta, sin_4theta, cos_4theta
                    ggl_rotation_matrix(len_ra,len_dec, src_ra, src_dec, rotation_mat);
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

                    temp_mnut_g = temp_mnut*sigma_crit;
                    temp_mnux_g = temp_mnux*sigma_crit;

                    // pdf
                    pre_bin_tag_sig1 = 0;
                    pre_bin_tag_sig2 = 0;
                    pre_bin_tag_shear1 = 0;
                    pre_bin_tag_shear2 = 0;

                    chi_pos = sep_bin_tag*data_info->chi_theta_block_len;

                    for(i=0; i<data_info->pdf_guess_num; i++)
                    {   
                        chi_pos_i = chi_pos + i*data_info->mg_bin_num;

                        // \Delta\Sigma(R) & \Delta\Sigma(R)_x
                        temp_mgt = src_mg1_rot - data_info->delta_sigma_guess[i]*temp_mnut;
                        temp_mgx = src_mg2_rot - data_info->delta_sigma_guess[i]*temp_mnux;
                        ggl_fast_hist(data_info->mg_bin, data_info->mg_bin_num, temp_mgt, pre_bin_tag_sig1, bin_tag_sig1);
                        ggl_fast_hist(data_info->mg_bin, data_info->mg_bin_num, temp_mgx, pre_bin_tag_sig2, bin_tag_sig2);
                        
                        data_info->worker_sub_chi_count[chi_pos_i + bin_tag_sig1] +=1;
                        data_info->worker_sub_chi_count[data_info->chi_signal_block_len + chi_pos_i + bin_tag_sig2] +=1;
                        
                        pre_bin_tag_sig1 = bin_tag_sig1;
                        pre_bin_tag_sig2 = bin_tag_sig2;

                        // \gamma_t & \gamma_x
                        temp_mgt = src_mg1_rot - data_info->gt_guess[i]*temp_mnut_g;
                        temp_mgx = src_mg2_rot - data_info->gt_guess[i]*temp_mnux_g;
                        ggl_fast_hist(data_info->mg_bin, data_info->mg_bin_num, temp_mgt, pre_bin_tag_shear1, bin_tag_shear1);
                        ggl_fast_hist(data_info->mg_bin, data_info->mg_bin_num, temp_mgx, pre_bin_tag_shear2, bin_tag_shear2);
                        
                        data_info->worker_sub_chi_count[data_info->chi_signal_block_len*2 + chi_pos_i + bin_tag_shear1] +=1;
                        data_info->worker_sub_chi_count[data_info->chi_signal_block_len*3 + chi_pos_i + bin_tag_shear2] +=1;
                        
                        pre_bin_tag_shear1 = bin_tag_shear1;
                        pre_bin_tag_shear2 = bin_tag_shear2;

                    }
                }
            }
        }
    }

    // add the count to the total count array, according to the jack id
    for(i=0; i<data_info->jack_num+1; i++)
    {   
        if(i == data_info->jack_id){continue;}

        for(j=0; j<data_info->chi_jack_block_len; j++)
        {
            data_info->worker_total_chi_count[i*data_info->chi_jack_block_len + j] += data_info->worker_sub_chi_count[j];
        }
    }

    ed = clock();
    char times[100];
    sprintf(times,"Finished in %.2f sec", (ed-st)/CLOCKS_PER_SEC);
    std::cout<<times<<std::endl;
}

void ggl_collect_chi(ggl_data_info *data_info)
{
    int i, j;
    MPI_Status status;

    if (data_info->rank > 0)
    {
        MPI_Send(data_info->worker_total_chi_count, data_info->total_chi_count_len, MPI_DOUBLE, 0, data_info->rank, MPI_COMM_WORLD);
    }
    else
    {
        for(i=1;i<data_info->numprocs;i++)
        {
            MPI_Recv(data_info->worker_total_chi_count, data_info->total_chi_count_len, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            
            for(j=0;j<data_info->total_chi_count_len;j++)
            {data_info->total_chi_count[j] += data_info->worker_total_chi_count[j];}
        }
    }   
}

void ggl_cal_signals(ggl_data_info * data_info)
{   
    if(data_info->rank == 0)
    {
        int i, j, st_c, st_j;

        double *temp_count = new double[data_info->chi_signal_block_len];
        double *signals = new double[data_info->signal_pts_num];
        double *signals_err = new double[data_info->signal_pts_num];


        double *delta_sigma_t = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
        double *delta_sigma_x = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
        double *gt = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
        double *gx = new double[(data_info->jack_num+1)*data_info->signal_pts_num];       

        double *theta = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
        double *radius = new double[(data_info->jack_num+1)*data_info->signal_pts_num];
        double *count = new double[(data_info->jack_num+1)*data_info->signal_pts_num];

        for(i=0; i<data_info->jack_id+1; i++)
        {   
            st_j = i*data_info->signal_pts_num;

            // delta_sigma_t
            st_c = i*data_info->chi_jack_block_len;
            for(j=0;j<data_info->chi_signal_block_len;j++)
            { temp_count[j] = data_info->total_chi_count[st_c+j];}
            ggl_pdf_signals(temp_count, data_info->delta_sigma_guess, data_info->pdf_guess_num, data_info->mg_bin_num, data_info->signal_pts_num, signals, signals_err);
            
            for(j=0; j<data_info->signal_pts_num; j++)
            {delta_sigma_t[st_j + j] = signals[j];}


            // delta_sigma_x
            st_c = i*data_info->chi_jack_block_len + data_info->chi_signal_block_len;
            for(j=0;j<data_info->chi_signal_block_len;j++)
            { temp_count[j] = data_info->total_chi_count[st_c+j];}
            ggl_pdf_signals(temp_count, data_info->delta_sigma_guess, data_info->pdf_guess_num, data_info->mg_bin_num, data_info->signal_pts_num, signals, signals_err);
            
            for(j=0; j<data_info->signal_pts_num; j++)
            {delta_sigma_x[st_j + j] = signals[j];}


            // g_t
            st_c = i*data_info->chi_jack_block_len + 2*data_info->chi_signal_block_len;
            for(j=0;j<data_info->chi_signal_block_len;j++)
            { temp_count[j] = data_info->total_chi_count[st_c+j];}
            ggl_pdf_signals(temp_count, data_info->delta_sigma_guess, data_info->pdf_guess_num, data_info->mg_bin_num, data_info->signal_pts_num, signals, signals_err);
            
            for(j=0; j<data_info->signal_pts_num; j++)
            {gt[st_j + j] = signals[j];}


            // g_x
            st_c = i*data_info->chi_jack_block_len + 3*data_info->chi_signal_block_len;
            for(j=0;j<data_info->chi_signal_block_len;j++)
            { temp_count[j] = data_info->total_chi_count[st_c+j];}
            ggl_pdf_signals(temp_count, data_info->delta_sigma_guess, data_info->pdf_guess_num, data_info->mg_bin_num, data_info->signal_pts_num, signals, signals_err);
            
            for(j=0; j<data_info->signal_pts_num; j++)
            {gx[st_j + j] = signals[j];}


            // theta
            st_c = i*data_info->chi_jack_block_len + 4*data_info->chi_signal_block_len;
            for(j=0;j<data_info->signal_pts_num;j++)
            { theta[st_j + j] = data_info->total_chi_count[st_c + j];}

            // radius
            st_c = i*data_info->chi_jack_block_len + 4*data_info->chi_signal_block_len + data_info->signal_pts_num;
            for(j=0;j<data_info->signal_pts_num;j++)
            { radius[st_j + j] = data_info->total_chi_count[st_c + j];}

            // count
            st_c = i*data_info->chi_jack_block_len + 4*data_info->chi_signal_block_len + 2*data_info->signal_pts_num;
            for(j=0;j<data_info->signal_pts_num;j++)
            { count[st_j + j] = data_info->total_chi_count[st_c + j];}           
        }

    }
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
