#include<ctime>
#include<iostream>
extern "C"
{
    void locate(double x, double *bins, int bin_num, int &bin_tag)
    {
        int i,j,st,ed,mid;
        st = 0;
        ed = bin_num;
 
        while(ed-st >= 1)
        {   
            mid = bin_num/2;
            if(x >= bins[mid]){st = mid;}
            else{ed = mid;}
        }
        bin_tag = st;
    }

    void hist2d(double *x, double *y, int num, double *xbin, double *ybin, int xbin_num, int ybin_num, int *count, double *grid_x, double *grid_y)
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
            grid_x[ytag*xbin_num + xtag] += x[i];
            grid_y[ytag*xbin_num + xtag] += y[i];
        }
        for(i=0;i=xbin_num*ybin_num;i++)
        {
            grid_x[i] = grid_x[i]/count[i];
            grid_y[i] = grid_y[i]/count[i];
        }
        t2 = clock();
        std::cout<<(t2-t1)/CLOCKS_PER_SEC<<std::endl;
    }

    void hist2d_fast(double *x, double *y, int num, double *xbin, double *ybin, int xbin_num, int ybin_num, int *count, double *grid_x, double *grid_y)
    {   
        int i, j, k;
        int xmid, ymid;
        
        int xtag, ytag;
        double t1, t2;

        t1 = clock();
        for(i=0; i<num; i++)
        {
            locate(x[i], xbin, xbin_num, xtag);

            locate(y[i], ybin, ybin_num, ytag);

            count[ytag*xbin_num + xtag] += 1;
            grid_x[ytag*xbin_num + xtag] += x[i];
            grid_y[ytag*xbin_num + xtag] += y[i];
        }
        for(i=0;i=xbin_num*ybin_num;i++)
        {
            grid_x[i] = grid_x[i]/count[i];
            grid_y[i] = grid_y[i]/count[i];
        }
        t2 = clock();
        std::cout<<(t2-t1)/CLOCKS_PER_SEC<<std::endl;
    }

}