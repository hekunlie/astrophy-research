#include<GG_lensing_functions_expo_wise.h>

void ggl_initialize(ggl_data_info *data_info)
{
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


    // no len/src data array exists in memory
    data_info->len_expo_read_tag = 0;
    data_info->src_expo_read_tag = 0;


    data_info->back_dz = 0.05;

    data_info->pos_inform_num = 6;

    // when read a len exposure, it will calculate the separation between the 
    // len and src and labelled the src exposures.
    data_info->src_expo_needed_tag = new int[data_info->src_expo_num];
    initialize_arr(data_info->src_expo_needed_tag, data_info->src_expo_num,-1);

    data_info->gt_guess = new MY_FLOAT[data_info->pdf_guess_num];
    data_info->delta_sigma_guess = new MY_FLOAT[data_info->pdf_guess_num];
    data_info->mg_bin = new MY_FLOAT[data_info->mg_bin_num];

    data_info->sepration_bin = new MY_FLOAT[data_info->signal_pts_num+1];

    data_info->chi_theta_block_len = data_info->pdf_guess_num*data_info->mg_bin_num;
    data_info->chi_signal_block_len = data_info->signal_pts_num*data_info->chi_theta_block_len;
    data_info->chi_jack_block_len = data_info->chi_signal_block_len*4 + data_info->signal_pts_num*3;
    // all pair data will be stacked into the last block in total_chi_count for the estimation of signals
    data_info->worker_total_chi_count = new double[data_info->chi_jack_block_len*(data_info->jack_num+1)];
    data_info->worker_sub_chi_count = new double[data_info->chi_jack_block_len];

    initialize_arr(data_info->worker_total_chi_count, data_info->chi_jack_block_len*(data_info->jack_num+1), 0);
    initialize_arr(data_info->worker_sub_chi_count, data_info->chi_jack_block_len, 0);

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


void find_src_needed(ggl_data_info *data_info, int len_expo_label)
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
        data_info->len_expo_data[ifg_row + data_info->len_com_dist_col]*(1+data_info->len_expo_data[ifg_row + data_info->len_z_col]);
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
    MY_FLOAT len_z, len_com_dist, src_z, src_com_dist;
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
            len_com_dist = data_info->len_expo_data[ifg_row + data_info->len_com_dist_col];

            // stacking in physical or comoving coordinate
#ifdef GGL_PROP_DIST_STACK
            coeff = (1+len_z)/len_com_dist;
#else
            coeff = 1./len_com_dist/(1+len_z);
#endif            
            for(ibkg=0; ibkg<data_info->src_data_row[bkg]; ibkg++)
            {   
                ibkg_row = ibkg*data_info->src_data_col;
                src_z = data_info->src_expo_data[ibkg_row + data_info->src_z_col] - data_info->back_dz;
                if(len_z <= src_z){continue;}

                src_ra = data_info->src_expo_data[ibkg_row + data_info->src_ra_col];
                src_dec = data_info->src_expo_data[ibkg_row + data_info->src_dec_col];
                src_com_dist = data_info->src_expo_data[ibkg_row + data_info->src_com_dist_col];

                sigma_crit = coeff*src_com_dist/(src_com_dist - len_com_dist);

                separation(len_ra, len_dec, src_ra, src_dec, sep_theta);

#ifdef GGL_PROP_DIST_STACK
                sep_dist = sep_theta*data_info->len_expo_data[ifg_row + data_info->len_com_dist_col]/(1+len_z);
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

