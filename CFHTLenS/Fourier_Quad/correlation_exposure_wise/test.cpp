#include<functions_expo_wise.h>
#include<hk_iolib.h>

struct data_test
{
    double *g1;
    double *g2;
    double *gg;
};

int main(int argc, char **argv)
{
    char file_path[400], hdf5_path[400];
    char set_name[20];
    
    sprintf(file_path,"test_data.hdf5");


    double a,b,c;
    double t1, t2, t3,t4,st1, st2, st3;
    int i, j, k;
    int num = 1000000;
    data_test dd_test;
    dd_test.g1 = new double[num]{};
    dd_test.g2 = new double[num]{};
    dd_test.gg = new double[2*num]{};



    int tag;
    t1 = clock();
    tag =0;
    for(i=0;i<10000;i++)
    {
        for(j=0;j<10000;j++)
        {   
            for(k=0; k<40; k++)
            {
                a = dd_test.g1[tag];
                b = dd_test.g2[tag];
                tag ++;
            }
            if(tag == num){tag =0;}
        }
    }
    t2 = clock();
    tag = 0;
    for(i=0;i<10000;i++)
    {
        for(j=0;j<10000;j++)
        {   
            for(k=0; k<80; k = k+2)
            {
                a = dd_test.gg[tag];
                b = dd_test.gg[tag+1];
                tag ++;
            }
            if(tag == num){tag =0;}
        }
    }
    t3 = clock();
    st1 = (t2-t1)/CLOCKS_PER_SEC;
    st2 = (t3-t2)/CLOCKS_PER_SEC;

    std::cout<<st1<<" "<<st2<<std::endl;
    return 0;
}