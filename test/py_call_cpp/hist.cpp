#include<ctime>
#include<iostream>
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

}