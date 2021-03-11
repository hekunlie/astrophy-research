#include<ctime>
#include<iostream>
#include<cmath>

extern "C"
{
    void test1(double x, double y)
    {
        double c = x+y;
        // std::cout<<c<<std::endl;
    }

    // this one doesn't work, "segmentation fault"
    void test2(double x, double y, double &z)
    {
        z = x+y;
        // std::cout<<z<<std::endl;
    }
    
    double test3(double x, double y)
    {
        return x+y;
    }

    // 1-D array
    double test4(int a, double *data)
    {
        double b = 0;
        int i;
        for(i=0;i<a;i++) b+=data[i];
        return b;
    }

    // double test5(int a, double *data, int b, double *data2)
    // {
    //     double c = 0;
    //     int i;
    //     for(i=0;i<a;i++) c+=data[i];
    //     std::cout<<"Count "<<c<<std::endl;
    //     for(i=0;i<b;i++)
    //     {
    //         data2[i] = c;
    //         std::cout<<i<<" "<<c<<std::endl;
    //     }
    // }


    double test(double *data1, int len1, double *data2, int len2)
    {
        double accum = 0;
        int i;
        for(i=0; i<len1; i++) accum+=data1[i];

        for(i=0; i<len2; i++)data2[i] = accum;
    }


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
    
    double get_chi2(double *xbin, int xbin_num, double *ybin, int ybin_num, int *count_xy, double ghat, double *Gbin, int Gbin_num)
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

    void deblending_self(float *ra1, float *dec1, float *z1, int len1, float sep_deg, float sep_z, int *blended_label)
    {
        // compare the sources in data_1 and data_2(close to data_1), and label the sources that are very close
        // to each other but with a large redshift difference
        int i,j,k;
        float dra, ddec, dz, ddeg;

        for(i=0; i<len1; i++)
        {
            for(j=i; j<len1; j++)
            {   
                dra = ra1[j] - ra1[i];
                ddec = dec1[j] - dec1[i];
                dz = abs(z1[j] - z1[i]);
                ddeg = sqrt(dra*dra + ddec*ddec);
                
                if(ddeg <= sep_deg and sep_z >= dz)
                {blended_label[i] = 1;}
            }
        }
    }

    void deblending_mutual(float *ra1, float *dec1, float *z1, int len1, float *ra2, float *dec2,float *z2, int len2, float sep_deg, float sep_z, int *blended_label)
    {
        // compare the sources in data_1 and data_2(close to data_1), and label the sources that are very close
        // to each other but with a large redshift difference
        int i,j,k;
        float dra, ddec, dz, ddeg;

        for(i=0; i<len1; i++)
        {
            if(blended_label[i] < 1)
            {
                for(j=0; j<len2; j++)
                {   
                    dra = ra2[j] - ra1[i];
                    ddec = dec2[j] - dec1[i];
                    dz = abs(z2[j] - z1[i]);
                    ddeg = sqrt(dra*dra + ddec*ddec);
                    
                    if(ddeg <= sep_deg and sep_z >= dz)
                    {blended_label[i] = 1;}
                }
            }

        }
    }

}