#include<FQlib.h>
#include<hk_iolib.h>

#define MY_FLOAT float

void hist_2d(MY_FLOAT x, MY_FLOAT*bins, int bin_num, int &ix)
{
    int i;
    for(i=0; i<bin_num; i++)
    {
        if(x > bins[i] and x <= bins[i+1]){ix=i;break;}
    }

}

void hist_2d_fast(MY_FLOAT x,  MY_FLOAT*bins, int bin_num, int bin_num2, int &ix)
{
    int i;
    if(x < 0)
    {
        for(i=0; i<bin_num2; i++)
        {
            if(x > bins[i] and x <= bins[i+1]){ix=i;break;}
        }
    }
    else
    {
        for(i=bin_num2; i<bin_num; i++)
        {
            if(x > bins[i] and x <= bins[i+1]){ix=i;break;}
        }
    }
}


void hist_2d_fast_new(MY_FLOAT x,  MY_FLOAT*bins, int bin_num, int bin_num1,int bin_num2, int bin_num3, int &ix)
{
    int i;
    if(x < 0)
    {
        if(x < bins[bin_num1])
        {
            for(i=0; i<bins[bin_num1]; i++)
            {if(x > bins[i] and x <= bins[i+1]){ix=i;break;}}
        }
        else
        {
            for(i=bin_num1; i<bin_num2; i++)
            {if(x > bins[i] and x <= bins[i+1]){ix=i;break;}}
        }
    }
    else
    {
        if(x < bins[bin_num3])
        {
            for(i=bin_num2; i<bin_num3; i++)
            {if(x > bins[i] and x <= bins[i+1]){ix=i;break;}}
        }
        else
        {
            for(i=bin_num3; i<bin_num; i++)
            { if(x > bins[i] and x <= bins[i+1]){ix=i;break;}}
        }
    }
}

int main()
{

    double st1, st2, st3, st4, st5, st6;
    int i, j, k;
    
    char data_path[400], set_name[50];
    int data_num, bin_num, bin_num1, bin_num2, bin_num3;
    data_num = 259661;
    MY_FLOAT *data = new MY_FLOAT[data_num]{};
    
    int loop_num = 2000;
    bin_num = 10;
    bin_num2 = bin_num/2;
    bin_num1 = bin_num2/2;
    bin_num3 = bin_num1 + bin_num2;
    MY_FLOAT *mg_bin = new MY_FLOAT[bin_num+1];

    MY_FLOAT *count1 = new MY_FLOAT[bin_num]{};
    MY_FLOAT *count2 = new MY_FLOAT[bin_num]{};
    MY_FLOAT *count3 = new MY_FLOAT[bin_num]{};

    sprintf(set_name,"/data");
    sprintf(data_path,"test.hdf5");
    read_h5(data_path, set_name, data);

    set_bin(data, data_num, mg_bin, bin_num,1000);
    show_arr(mg_bin, 1, bin_num+1);

    st1 = clock();
    // k = 0;
    // for(i=0;i<100000;i++)
    // {
    //     for(j=0;j<10000;j++)
    //     {
    //         k++;
    //         if(k == 20000)
    //         {k=0;}
    //     }
    // }
    // st2 = clock();
    // std::cout<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;

    // k = 0;
    // for(i=0;i<100000;i++)
    // {
    //     for(j=0;j<10000;j++)
    //     {
    //         k = k%20000;
    //         k++;
    //     }
    // }
    // st2 = clock();
    // std::cout<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
    for(i=0;i<loop_num;i++)
    {   
        for(j=0;j<data_num;j++)
        {hist_2d(data[j], mg_bin, bin_num, k); count1[k] += 1;}
    }
    st2 = clock();
    for(i=0;i<loop_num;i++)
    {   
        for(j=0;j<data_num;j++)
        {hist_2d_fast(data[j], mg_bin, bin_num,bin_num2, k);count2[k] += 1;}
    }
    st3 = clock();
    for(i=0;i<loop_num;i++)
    {   
        for(j=0;j<data_num;j++)
        {hist_2d_fast_new(data[j], mg_bin, bin_num,bin_num1,bin_num2,bin_num3, k);count3[k] += 1;}
    }
    st4 = clock();
    show_arr(count1, 1, bin_num);
    show_arr(count2, 1, bin_num);
    show_arr(count2, 1, bin_num);

    for(i=0; i<bin_num; i++)
    {
        std::cout<<count1[i] - count2[i]<<" "<<count2[i] - count3[i]<<" "<<count1[i] - count3[i]<<std::endl;
    }

    std::cout<<(st2-st1)/CLOCKS_PER_SEC<<" "<<(st3-st2)/CLOCKS_PER_SEC<<" "<<(st4-st3)/CLOCKS_PER_SEC<<std::endl;


    return 0;
}