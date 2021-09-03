#include<ctime>
#include<iostream>
#include<cmath>
// g++  c4py.cpp -shared -fPIC -o libc4py.so -std=c++11

extern "C"
{
    int locate(double x, double *bins, int bin_num)
    {
        int i,j,st,ed,mid;
        st = 0;
        ed = bin_num;
 
        while(true)
        {   
            if(ed - st <= 1){break;}

            mid = (st+ed)/2;
            if(x >= bins[mid]){st = mid;}
            else{ed = mid;}
        }
        return st;
    }
    
    int locate_f(float x, float *bins, int bin_num)
    {
        int i,j,st,ed,mid;
        st = 0;
        ed = bin_num;
 
        while(true)
        {   
            if(ed - st <= 1){break;}

            mid = (st+ed)/2;
            if(x >= bins[mid]){st = mid;}
            else{ed = mid;}
        }
        return st;
    }
    
    void hist1d_fq(float x, float *pdf_bin, int bin_num, int pre_bin_tag, int& bin_tag)
    {   
        int i, j;
        if(x >= pdf_bin[pre_bin_tag])
        {
            for(i=pre_bin_tag; i<bin_num; i++)
            {
                if(x >= pdf_bin[i] and x < pdf_bin[i+1]){bin_tag = i;break;}
                // std::cout<<x<<" "<<pdf_bin[i]<<" "<<pdf_bin[i+1]<<" "<<bin_tag<<" "<<i<<std::endl;
            }
        }
        else
        {
            for(i=pre_bin_tag; i>0; i--)
            {   
                j = i-1;
                if(x >= pdf_bin[j] and x < pdf_bin[i]){bin_tag = j;break;}
            }
        }
    }


    void hist1d_fast(double *x, int num, double *xbin,  int xbin_num, int *count)
    {   
        int i;

        int xst, xed, xmid;
        for(i=0;i<xbin_num;i++){count[i] = 0;}


        for(i=0; i<num; i++)
        {

            xst = 0;
            xed = xbin_num;
            while(xed - xst>1)
            {   
                xmid = (xst+xed)/2;
                if(x[i] >= xbin[xmid]){xst = xmid;}
                else{xed = xmid;}
            }
            count[xst] += 1;

        }
    }


    void hist2d(double *x, double *y, int num, double *xbin, double *ybin, int xbin_num, int ybin_num, int *count)//, double *grid_x, double *grid_y)
    {   
        int i, j, k;
        int xmid, ymid;
        
        int xtag, ytag;
        double t1, t2;

        t1 = clock();
        for(i=0; i<num; i++)
        {
            for(j=0;j<xbin_num;j++)
            {
                if(x[i]>=xbin[j] && x[i]<xbin[j+1])
                {
                    xtag = j;break;
                }
            }

            for(j=0;j<ybin_num;j++)
            {
                if(y[i]>=ybin[j] && y[i]<ybin[j+1])
                {
                    ytag = j;break;
                }
            }

            count[ytag*xbin_num + xtag] += 1;
            // grid_x[ytag*xbin_num + xtag] += x[i];
            // grid_y[ytag*xbin_num + xtag] += y[i];
        }
        // for(i=0;i=xbin_num*ybin_num;i++)
        // {
        //     grid_x[i] = grid_x[i]/count[i];
        //     grid_y[i] = grid_y[i]/count[i];
        // }
        t2 = clock();
        std::cout<<(t2-t1)/CLOCKS_PER_SEC<<std::endl;
    }

    void hist2d_fast(double *x, double *y, int num, double *xbin, double *ybin, int xbin_num, int ybin_num, int *count, double *grid_x, double *grid_y)
    {   
        int i, j, k;

        int xst, xed, xmid;
        int yst, yed, ymid;

        // double t1, t2;

        // t1 = clock();
        for(i=0; i<num; i++)
        {

            xst = 0;
            xed = xbin_num;
            while(xed - xst>1)
            {   
                xmid = (xst+xed)/2;
                if(x[i] >= xbin[xmid]){xst = xmid;}
                else{xed = xmid;}
            }
 
            
            yst = 0;
            yed = ybin_num;
            while(yed - yst > 1)
            {   
                ymid = (yst+yed)/2;
                if(y[i] >= ybin[ymid]){yst = ymid;}
                else{yed = ymid;}
            }

            
            k = yst*xbin_num + xst;
            count[k] += 1;
            grid_x[k] += x[i];
            grid_y[k] += y[i];
        }
        for(i=0;i<xbin_num*ybin_num;i++)
        {   
            if(count[i] > 0)
            {   
                grid_x[i] = grid_x[i]/count[i];
                grid_y[i] = grid_y[i]/count[i];
            }
            else
            {
                grid_x[i] = 0;
                grid_y[i] = 0;
            }
        }
        // t2 = clock();
        // std::cout<<(t2-t1)/CLOCKS_PER_SEC<<std::endl;
    }
    

    void find_overlap(float *ra, float *dec, int data_len, float*ra_bin, float *dec_bin, int ra_bin_num, int dec_bin_num, int *ra_dec_mask, int *labels)
    {
        int i, j, k;
        int ix, iy;
        float x, y;
        for(i=0; i< data_len; i++)
        {
            ix = locate_f(ra[i], ra_bin, ra_bin_num);
            iy = locate_f(dec[i], dec_bin, dec_bin_num);
            
            k = iy*ra_bin_num + ix;
            if(ra_dec_mask[k] > 0)
            {labels[i] = 1;}
        }
    }


    void get_chi2(double *xbin, int xbin_num, double *ybin, int ybin_num, int *count_xy, double ghat, double *Gbin, int Gbin_num)
    {
        // y = 1/ghat *(x - x_G) are the lines of the bins' boundaries
        // x_Gs are the value of Gbin

        int i,j,k, xy_tag;
        int *count_1d = new int[Gbin_num];
        double temp;
        double *xs = new double[ybin_num+1];
        double xs1, xs2;
        double xs1_G_k1, xs2_G_k1,xs1_G_k2, xs2_G_k2;
        // coordinates of the four coners of the grid in the 2d histogram
        double cx1, cx2, dx, cy1, cy2, dy, ds;
        double ix1, ix2, iy1, iy2;


        for(i=0; i<ybin_num+1; i++)
        {   
            // without x_G
            xs[i] = ghat*ybin[i];
        }
        if(ghat == 0)
        {
            // caution
            ;
        }
        else if(ghat > 0)
        {
            // ghat > 0 case
            for(i=0; i<ybin_num; i++)
            {   
                // the x on the lines without x_G,
                // which will be compared with the x of
                // grid to determine which bin this grid
                // belongs to.
                xs1 = xs[i];
                xs2 = xs[i+1];

                cy1 = ybin[i];
                cy2 = ybin[i+1];
                dy = ybin[i+1] - ybin[i];
                for(j=0; j<xbin_num; j++)
                {   
                    xy_tag = i*xbin_num + j;
                    if(count_xy[xy_tag] < 1){continue;}

                    cx1 = xbin[j];
                    cx2 = xbin[j+1];
                    dx = xbin[j+1]- xbin[j];

                    ds = dx*dy;

                    for(k= 0; k<Gbin_num; k++)
                    {
                        xs1_G_k1 = xs1 + Gbin[k];
                        xs1_G_k2 = xs1 + Gbin[k+1];

                        xs2_G_k1 = xs2 + Gbin[k];                    
                        xs2_G_k2 = xs2 + Gbin[k+1];

                        if(cx1 >= xs2_G_k1 and cx1 < xs2_G_k2)
                        {   
                            if(cx2 <= xs1_G_k2)
                            {
                                // the grid is in the k'th Gbin
                                count_1d[k] += count_xy[xy_tag];
                                break;
                            }
                            else
                            {
                                // part of the grid is in the k'th grid
                                // and part of it is in the (k+1)'th grid
                                if(cx2 <= xs2_G_k2)
                                {
                                    if(xs1_G_k2 >= cx1)
                                    {
                                        // only one corner in the (k+1)'th grid
                                        temp = (cx2 - xs1_G_k2)*((cx2 - Gbin[k+1])/ghat*0.5 - cy1)/ds;
                                    }
                                    else
                                    {  
                                        // two corners in the (k+1)'th grid
                                        temp = ((cx2 - Gbin[k+1])/ghat*0.5 - cy1 + (cx1 - Gbin[k+1])/ghat*0.5 - cy1)*dx*0.5/ds;
                                    }
                                }
                                else
                                {
                                    if(xs1_G_k2 >= cx1)
                                    {
                                        // only one corner in the (k+1)'th grid
                                        temp = (cx2 - xs1_G_k2 + cx2 - xs2_G_k2)*dy*0.5/ds;
                                    }
                                    else
                                    {  
                                        // two corners in the (k+1)'th grid
                                        temp = 1 - (xs2_G_k2 - cx1)*(cy2 - cx1 + Gbin[k+1])*0.5/ds;
                                    }
                                }
                                count_1d[k] += count_xy[xy_tag] * (1-temp);
                                count_1d[k+1] += count_xy[xy_tag] * temp;
                                break;
                            }
                        }                        
                    }
                }
            }
        }
        else
        {
            // ghat < 0 case
            for(i=0; i<ybin_num; i++)
            {   
                // the x on the lines without x_G,
                // which will be compared with the x of
                // grid to determine which bin this grid
                // belongs to.
                xs1 = xs[i];
                xs2 = xs[i+1];

                cy1 = ybin[i];
                cy2 = ybin[i+1];
                dy = ybin[i+1] - ybin[i];
                for(j=0; j<xbin_num; j++)
                {   
                    xy_tag = i*xbin_num + j;
                    if(count_xy[xy_tag] < 1){continue;}

                    cx1 = xbin[j];
                    cx2 = xbin[j+1];
                    dx = xbin[j+1]- xbin[j];

                    ds = dx*dy;

                    for(k= 0; k<Gbin_num; k++)
                    {
                        xs1_G_k1 = xs1 + Gbin[k];
                        xs1_G_k2 = xs1 + Gbin[k+1];

                        xs2_G_k1 = xs2 + Gbin[k];                    
                        xs2_G_k2 = xs2 + Gbin[k+1];

                        if(cx1 >= xs2_G_k1 and cx1 < xs2_G_k2)
                        {   
                            if(cx2 <= xs1_G_k2)
                            {
                                // the grid is in the k'th Gbin
                                count_1d[k] += count_xy[xy_tag];
                                break;
                            }
                            else
                            {
                                // part of the grid is in the k'th grid
                                // and part of it is in the (k+1)'th grid
                                if(cx2 <= xs2_G_k2)
                                {
                                    if(xs1_G_k2 >= cx1)
                                    {
                                        // only one corner in the (k+1)'th grid
                                        temp = (cx2 - xs1_G_k2)*((cx2 - Gbin[k+1])/ghat*0.5 - cy1)/ds;
                                    }
                                    else
                                    {  
                                        // two corners in the (k+1)'th grid
                                        temp = ((cx2 - Gbin[k+1])/ghat*0.5 - cy1 + (cx1 - Gbin[k+1])/ghat*0.5 - cy1)*dx*0.5/ds;
                                    }
                                }
                                else
                                {
                                    if(xs1_G_k2 >= cx1)
                                    {
                                        // only one corner in the (k+1)'th grid
                                        temp = (cx2 - xs1_G_k2 + cx2 - xs2_G_k2)*dy*0.5/ds;
                                    }
                                    else
                                    {  
                                        // two corners in the (k+1)'th grid
                                        temp = 1 - (xs2_G_k2 - cx1)*(cy2 - cx1 + Gbin[k+1])*0.5/ds;
                                    }
                                }
                                count_1d[k] += count_xy[xy_tag] * (1-temp);
                                count_1d[k+1] += count_xy[xy_tag] * temp;
                                break;
                            }
                        }                        
                    }
                }
            }
        }
 
        delete[] count_1d;
        delete[] xs;
    }
    
    void deblend(float *xc, float *yc, float *z, float *radius, int *chip_label, int data_len, float sep_radius_scale, float sep_z, int *blended_label)
    {
        // compare the sources in data_1 and data_2(close to data_1), and label the sources that are very close
        // to each other but with a large redshift difference
        int i,j,k;
        float xc_i, yc_i, z_i, rad_i;
        float dr, r2, dx, dy, dz;

        for(i=0; i<data_len; i++)
        {   
            xc_i = xc[i];
            yc_i = yc[i];
            z_i = z[i];
            rad_i = radius[i];           

            for(j=i+1; j<data_len; j++)
            {   
                if(chip_label[i] == chip_label[j])
                {
                    r2 = (rad_i + radius[j])*sep_radius_scale;
                    dx = fabs(xc_i - xc[j]) - r2;

                    if(dx <= 0)
                    {
                        dy = fabs(yc_i - yc[j]) - r2;
                        // dr = sqrt(dx*dx + dy*dy) - r2;
                        dz = fabs(z_i - z[j]);
                        if(dy <= 0 and dz >= sep_z)
                        {
                            blended_label[i] = -1;
                            blended_label[j] = -1;
                        }
                    }
                }
            }
        }
    }

    void deblend_i(float *xc, float *yc, float *z, int *chip_label, int data_len, float sep_pix, float sep_z, int *blended_label)
    {
        // compare the sources in data_1 and data_2(close to data_1), and label the sources that are very close
        // to each other but with a large redshift difference
        int i,j,k;
        float xc_i, yc_i, z_i, rad_i;
        float dr, r2, dx, dy, dz;

        for(i=0; i<data_len; i++)
        {   
            xc_i = xc[i];
            yc_i = yc[i];
            z_i = z[i];
        
            for(j=i+1; j<data_len; j++)
            {   
                if(chip_label[i] == chip_label[j])
                {
                    dx = fabs(xc_i - xc[j]);

                    if(dx <= sep_pix)
                    {
                        dy = fabs(yc_i - yc[j]);
                        // dr = sqrt(dx*dx + dy*dy);
                        dz = fabs(z_i - z[j]);

                        if(dy <= sep_pix and dz >= sep_z)
                        {
                            blended_label[i] = -1;
                            blended_label[j] = -1;
                            //std::cout<<xc_i<<" "<<xc[j]<<" "<<dx<<std::endl;
                        }
                    }
                }
            }
        }
    }

    void deblend_mutual(float *ra1, float *dec1, float *z1, int len1, 
                        float *ra2, float *dec2,float *z2, int len2, float sep_ang, float sep_z, int *blended_label)
    {
        // compare the sources in data_1 and data_2(close to data_1), and label the sources that are very close
        // to each other but with a large redshift difference
        int i,j,k;
        float dra, ddec, dz, dr;

        for(i=0; i<len1; i++)
        {
            if(blended_label[i] > 0)
            {
                for(j=0; j<len2; j++)
                {   
                    dra = ra2[j] - ra1[i];
                    ddec = dec2[j] - dec1[i];
                    dz = abs(z2[j] - z1[i]);
                    dr = sqrt(dra*dra + ddec*ddec);
                    
                    if(dr <= sep_ang and sep_z >= dz)
                    {blended_label[i] = -1;}
                }
            }
        }
    }

    void  cata_match(float *cata_ra1, float*cata_dec1, int cata_len1, 
                    float *cata_ra2, float*cata_dec2, int cata_len2, 
                    int *seq_labels, int *check_labels, float diff_ang_thresh)
    {  
        // match the source between these two catalogue
        // find the corresponding pairs in cata2, 
        // seq_labels label the row index in cata2 
        int i,j,k;
        float diff_rad, diff_z, diff_rad_min, diff_z_min;
        float diff_ra, diff_dec;
        int min_tag_1, min_tag_2;
        float ra1, dec1;

        for(i=0; i<cata_len1; i++)
        {   
            min_tag_1 = -1;
            // min_tag_2 = -2;
            diff_rad_min = 1;
            // diff_z_min = 1;
            ra1 = cata_ra1[i];
            dec1 = cata_dec1[i];

            for(j=0; j<cata_len2; j++)
            {   
                diff_ra = ra1 - cata_ra2[j];
                diff_dec = dec1 - cata_dec2[j];
                
                // diff_rad = abs(cata_ra1[i] - cata_ra2[j]) + abs(cata_dec1[i] - cata_dec2[j]);
                diff_rad = sqrt(diff_ra*diff_ra + diff_dec*diff_dec);

                if(diff_rad < diff_rad_min) 
                {
                    diff_rad_min = diff_rad;
                    min_tag_1 = j;
                }
            }

            if(diff_rad_min < diff_ang_thresh )
            {
                seq_labels[i] = min_tag_1;
                check_labels[i] += 1;
            }
        }
    }

    double search_shear_range_chi2(float *mg, float *mnu, float *count_weight, int data_num, float signal, 
                                    float *mg_bin, int mg_bin_num)
    {
        int i, bin_num2;
        int xst, xed, xmid;
        float mgh;
        double chisq, n1, n2;

        double *count = new double[mg_bin_num];

        for(i=0;i<mg_bin_num;i++){count[i] = 0;}

        for(i=0; i<data_num; i++)
        {
            mgh = mg[i] - signal*mnu[i];

            xst = 0;
            xed = mg_bin_num;
            while(xed - xst>1)
            {   
                xmid = (xst+xed)/2;
                if(mgh >= mg_bin[xmid]){xst = xmid;}
                else{xed = xmid;}
            }
            count[xst] += count_weight[i];
        }

        chisq = 0;
        bin_num2 = mg_bin_num/2;
        for(i=0;i<bin_num2; i++)
        {   
            n1 = count[i] - count[mg_bin_num-i-1];
            n2 = count[i] + count[mg_bin_num-i-1];
            chisq += n1*n1/n2;
        }
        chisq = chisq/2;

        delete[] count;
        return chisq;
    }

    void search_shear_range(float *mg, float *mnu, float *count_weight, int data_num, float *mg_bin, int mg_bin_num, 
                            float left_shear_guess, float right_shear_guess, int max_iters, float chi2_gap, 
                            float *fit_shear_range, float *fit_chi2, int fit_shear_num)
    {
        int i, j, change, iters;
        float mc, mcl, mcr, left, right;
        float fmc, fmcl, fmcr, temp;

        left = left_shear_guess;
        right = right_shear_guess;
        change = 1;
        iters = 0;
        while(true)
        {
            change = 0;
            mc = (left + right) / 2;
            mcl = left;
            mcr = right;

            fmc = search_shear_range_chi2(mg, mnu, count_weight, data_num, mc, mg_bin, mg_bin_num);
            fmcl = search_shear_range_chi2(mg, mnu, count_weight, data_num, mcl, mg_bin, mg_bin_num);
            fmcr = search_shear_range_chi2(mg, mnu, count_weight, data_num, mcr, mg_bin, mg_bin_num);

            temp = fmc + chi2_gap;

            if (fmcl > temp)
            { 
                left = mcl + (mc - mcl) / 3;
                change = 1;
            }
            if (fmcr > temp)
            {
                right = mcr - (mcr - mc) / 3;
                change = 1;
            }

            if(change == 0){break;}

            iters += 1;
            if (iters > max_iters){break;}
        }
        
        temp = (right - left)/(fit_shear_num-1);
        for(i=0;i<fit_shear_num;i++)
        {
            fit_shear_range[i] = left + i*temp;
            fit_chi2[i] = search_shear_range_chi2(mg, mnu, count_weight, data_num, fit_shear_range[i], mg_bin, mg_bin_num);
        }

    }


    double search_shear_range_chi2_corr(double *mg, double *mnu, int data_num, double *mg_corr, double *mnu_corr, int corr_num, double signal, 
                                    double *mg_bin, int mg_bin_num)
    {
        int i, bin_num2;
        int xst, xed, xmid;
        double mgh;
        double chisq, n1, n2, n1_corr, n2_corr;

        double *count = new double[mg_bin_num];
        double *count_corr = new double[mg_bin_num];

        for(i=0;i<mg_bin_num;i++){count[i] = 0;count_corr[i] = 0;}

        for(i=0; i<data_num; i++)
        {
            mgh = mg[i] - signal*mnu[i];

            xst = 0;
            xed = mg_bin_num;
            while(xed - xst>1)
            {   
                xmid = (xst+xed)/2;
                if(mgh >= mg_bin[xmid]){xst = xmid;}
                else{xed = xmid;}
            }
            count[xst] += 1;
        }

        for(i=0; i<corr_num; i++)
        {
            mgh = mg_corr[i] - signal*mnu_corr[i];

            xst = 0;
            xed = mg_bin_num;
            while(xed - xst>1)
            {   
                xmid = (xst+xed)/2;
                if(mgh >= mg_bin[xmid]){xst = xmid;}
                else{xed = xmid;}
            }
            count_corr[xst] += 1;
        }

        for(i=0;i<mg_bin_num;i++)
        {   
            if(count[i] >= count_corr[i])
            {count[i] = count[i] - count_corr[i];}
            else
            {count[i] = 0;}
        }

        chisq = 0;
        bin_num2 = mg_bin_num/2;
        for(i=0;i<bin_num2; i++)
        {   
            n2 = count[i] + count[mg_bin_num-i-1];
            if(n2 == 0){continue;}

            n1 = count[i] - count[mg_bin_num-i-1];
            chisq += n1*n1/n2;
        }

        chisq = chisq/2;

        delete[] count;
        delete[] count_corr;

        return chisq;
    }


    void search_shear_range_corr(double *mg, double *mnu, int data_num, double *mg_corr, double *mnu_corr, int corr_num, double *mg_bin, int mg_bin_num, 
                            double left_shear_guess, double right_shear_guess, int max_iters, double chi2_gap, 
                            double *fit_shear_range, double *fit_chi2, int fit_shear_num)
    {
        int i, j, change, iters;
        double mc, mcl, mcr, left, right;
        double fmc, fmcl, fmcr, temp;

        left = left_shear_guess;
        right = right_shear_guess;
        change = 1;
        iters = 0;
        while(true)
        {
            change = 0;
            mc = (left + right) / 2;
            mcl = left;
            mcr = right;

            fmc = search_shear_range_chi2_corr(mg, mnu, data_num, mg_corr, mnu_corr, corr_num, mc, mg_bin, mg_bin_num);
            fmcl = search_shear_range_chi2_corr(mg, mnu, data_num, mg_corr, mnu_corr, corr_num, mcl, mg_bin, mg_bin_num);
            fmcr = search_shear_range_chi2_corr(mg, mnu, data_num, mg_corr, mnu_corr, corr_num, mcr, mg_bin, mg_bin_num);

            temp = fmc + chi2_gap;

            if (fmcl > temp)
            { 
                left = mcl + (mc - mcl) / 3;
                change = 1;
            }
            if (fmcr > temp)
            {
                right = mcr - (mcr - mc) / 3;
                change = 1;
            }

            if(change == 0){break;}

            iters += 1;
            if (iters > max_iters){break;}
        }
        
        temp = (right - left)/(fit_shear_num-1);
        for(i=0;i<fit_shear_num;i++)
        {
            fit_shear_range[i] = left + i*temp;
            fit_chi2[i] = search_shear_range_chi2_corr(mg, mnu, data_num, mg_corr, mnu_corr, corr_num, fit_shear_range[i], mg_bin, mg_bin_num);
        }

    }
    

    void cal_chisq(float *mg, float*mu, float*count_weight, int data_num, float *pdf_bin, int bin_num, float *guess, int guess_num,  double*count)
    {
        int i, j, k;
        float temp;
        int pre_bin_tag, bin_tag, tag;


        for(i=0; i<data_num; i++)
        {
            pre_bin_tag = 0;
            for(j=0; j<guess_num; j++)
            {
                temp = mg[i] - guess[j]*mu[i];

                hist1d_fq(temp, pdf_bin, bin_num, pre_bin_tag, bin_tag);

                tag = j*bin_num + bin_tag;
                count[tag] += count_weight[i];

                pre_bin_tag = bin_tag;                
            }
        }

    }
}