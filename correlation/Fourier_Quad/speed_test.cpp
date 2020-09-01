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

void cal(MY_FLOAT a, MY_FLOAT b, int ia, int ib, MY_FLOAT c,MY_FLOAT d, int &ic,int &id)
{

    if (a>c){ic = ia;}
    else{ic = ia -1;}

    if (b>d){id = ib;}
    else{id = ib -1;}
    
}
void cal(MY_FLOAT* abcd, int ia, int ib, int &ic,int &id)
{

    if (abcd[0]>abcd[2]){ic = ia;}
    else{ic = ia -1;}

    if (abcd[1]>abcd[3]){id = ib;}
    else{id = ib -1;}
    
}

void cal(MY_FLOAT* abcd, int *i_ab, int &ic,int &id)
{

    if (abcd[0]>abcd[2]){ic = i_ab[0];}
    else{ic = i_ab[0] -1;}

    if (abcd[1]>abcd[3]){id = i_ab[1];}
    else{id = i_ab[1] -1;}
    
}
int main(int argc, char **argv)
{

    double st1, st2, st3, st4, st5, st6;
    double tt1, tt2, tt3, tt4;
    int i, j, k,m,n;
    
    char data_path[400], set_name[50];
    int data_num, bin_num, bin_num1, bin_num2, bin_num3;
    data_num = 10000;
    MY_FLOAT *data = new MY_FLOAT[data_num]{};
    register MY_FLOAT re_data[10000];

    int loop_num = atoi(argv[1]);
    int loop_num_1 = atoi(argv[2]);
    bin_num = 10;
    bin_num2 = bin_num/2;
    bin_num1 = bin_num2/2;
    bin_num3 = bin_num1 + bin_num2;
    MY_FLOAT *mg_bin = new MY_FLOAT[bin_num+1];

    MY_FLOAT *count1 = new MY_FLOAT[bin_num]{};
    MY_FLOAT *count2 = new MY_FLOAT[bin_num]{};
    MY_FLOAT *count3 = new MY_FLOAT[bin_num]{};

    MY_FLOAT a,b,c,d;
    int ia, ib, ic, id;
    MY_FLOAT *abcd = new MY_FLOAT[4]{1,2,3,4};
    MY_FLOAT abcd_[4] = {1,2,3,4};
    ia = 0;
    ib = 0;
    ic = 0;
    id = 0;
    int i_ab[2] = {0,0};

    // sprintf(set_name,"/data");
    // sprintf(data_path,"test.hdf5");
    // read_h5(data_path, set_name, data);

    // set_bin(data, data_num, mg_bin, bin_num,1000);
    // show_arr(mg_bin, 1, bin_num+1);

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
    // for(i=0;i<loop_num;i++)
    // {   
    //     for(j=0;j<data_num;j++)
    //     {hist_2d(data[j], mg_bin, bin_num, k); count1[k] += 1;}
    // }
    // st2 = clock();
    // for(i=0;i<loop_num;i++)
    // {   
    //     for(j=0;j<data_num;j++)
    //     {hist_2d_fast(data[j], mg_bin, bin_num,bin_num2, k);count2[k] += 1;}
    // }
    // st3 = clock();
    // for(i=0;i<loop_num;i++)
    // {   
    //     for(j=0;j<data_num;j++)
    //     {hist_2d_fast_new(data[j], mg_bin, bin_num,bin_num1,bin_num2,bin_num3, k);count3[k] += 1;}
    // }
    // st4 = clock();
    // show_arr(count1, 1, bin_num);
    // show_arr(count2, 1, bin_num);
    // show_arr(count2, 1, bin_num);

    // for(i=0; i<bin_num; i++)
    // {
    //     std::cout<<count1[i] - count2[i]<<" "<<count2[i] - count3[i]<<" "<<count1[i] - count3[i]<<std::endl;
    // }
    // for(i=0;i<loop_num_1;i++)
    // {   
    //     a = i;
    //     b = i;
    //     c = i;
    //     d = i;
    //     ia = 0;
    //     ib = 0;
    //     for(j=0;j<loop_num;j++)
    //     {
    //         cal(a,b,ia,ib,c,d,ic,id);
    //     }
    // }   

    // st2 = clock();
    // for(i=0;i<loop_num_1;i++)
    // { 
    //     a = i;
    //     b = i;
    //     c = i;
    //     d = i;
    //     ia = 0;
    //     ib = 0;
    //     for(j=0;j<loop_num;j++)
    //     {
    //         cal(abcd,ia,ib,ic, id);
    //     }
    // }  
    // st3 = clock();
    
    // for(i=0;i<loop_num_1;i++)
    // { 
    //     a = i;
    //     b = i;
    //     c = i;
    //     d = i;
    //     i_ab[0] = 0;
    //     i_ab[1] = 0;
    //     for(j=0;j<loop_num;j++)
    //     {
    //         cal(abcd, i_ab, ic, id);
    //     }
    // }
    // st4 = clock();
    // std::cout<<(st2-st1)/CLOCKS_PER_SEC<<" "<<(st3-st2)/CLOCKS_PER_SEC<<" "<<(st4-st3)/CLOCKS_PER_SEC<<std::endl;

    st2 = clock();
    k = 0;
    m = 0;
    // n = 1000;
    // for(i=0;i<1000;i++){re_data[i] = data[i];}
    for(i=0; i<loop_num_1; i++)
    {
        for(j=0;j<loop_num;j++)
        {
            a = re_data[k];
            if(k >= 10000)
            {   
                // if(n >= 10000){n=0;}
                // for(m=0;m<1000;m++){re_data[m] = data[m+n];}
                k = 0;
            }
            // n += 1000;
            k++;
        }
    }
    st3 = clock();
    k = 0;
    for(i=0; i<loop_num_1; i++)
    {
        for(j=0;j<loop_num;j++)
        {
            a = *(data+k);
            if(k >= data_num){k=0;}
            k++;
        }
    }
    st4 = clock();
    tt1 = (st2-st1)/CLOCKS_PER_SEC;
    tt2 = (st3-st2)/CLOCKS_PER_SEC;
    tt3 = (st4-st3)/CLOCKS_PER_SEC;

    std::cout<<tt1<<" "<<tt2<<" "<<tt3<<" "<<tt3/tt2<<std::endl;
    char log_inform[100];
    a = atof(argv[3]);
    sprintf(log_inform, "%g, %.0f",a,a);
    std::cout<<log_inform<<std::endl;
    return 0;
}