#include<cmath>
#include<ctime>
#include<FQlib.h>

int main(int argc, char ** argv)
{
    int seed = 123;
    int i,j,k;
    int num = 10000000;
    double s1,s2,s3;
    double *rand_val =new double[num]{};
    double *exp_val = new double[num]{};
    double *exp_val_fit = new double[num]{};

    double st1, st2;
    double rand_min, rand_max;
    rand_min = 0.01;
    rand_max = 10;

    gsl_initialize(seed);

    for(i=0;i<num;i++)
    {
        rand_uniform(rand_min, rand_max, rand_val[i]);
    }
    show_arr(rand_val,1,100);

    st1 = clock();
    for(i=0;i<num;i++)
    {
        s1 = rand_val[i] + rand_val[0];
    }
    st2 = clock();
    std::cout<<num<<" times (double) adding use "<<(st2-st1)/CLOCKS_PER_SEC<<" sec"<<std::endl;

    st1 = clock();
    for(i=0;i<num;i++)
    {
        s1 = rand_val[i] - rand_val[0];
    }
    st2 = clock();
    std::cout<<num<<" times (double) subtraction use "<<(st2-st1)/CLOCKS_PER_SEC<<" sec"<<std::endl;

    st1 = clock();
    for(i=0;i<num;i++)
    {
        s1 = rand_val[i] * rand_val[0];
    }
    st2 = clock();
    std::cout<<num<<" times (double) multiply use "<<(st2-st1)/CLOCKS_PER_SEC<<" sec"<<std::endl;

    st1 = clock();
    for(i=0;i<num;i++)
    {
        s1 = rand_val[i]/rand_val[0];
    }
    st2 = clock();
    std::cout<<num<<" times (double) division use "<<(st2-st1)/CLOCKS_PER_SEC<<" sec"<<std::endl;

    st1 = clock();
    for(i=0;i<num;i++)
    {
        s1 = pow(rand_val[i],2);
    }
    st2 = clock();
    std::cout<<num<<" times (double) pow^2 use "<<(st2-st1)/CLOCKS_PER_SEC<<" sec"<<std::endl;

    st1 = clock();
    for(i=0;i<num;i++)
    {
        s1 = pow(rand_val[i],4);
    }
    st2 = clock();
    std::cout<<num<<" times (double) pow^4 use "<<(st2-st1)/CLOCKS_PER_SEC<<" sec"<<std::endl;

    st1 = clock();
    for(i=0;i<num;i++)
    {
        exp_val[i] = exp(rand_val[i]);
    }
    st2 = clock();
    std::cout<<num<<" times (double) exp use "<<(st2-st1)/CLOCKS_PER_SEC<<" sec"<<std::endl;

    st1 = clock();
    for(i=0;i<num;i++)
    {
        j = int(rand_val[i]);
    }
    st2 = clock();
    std::cout<<num<<" times (double) to int use "<<(st2-st1)/CLOCKS_PER_SEC<<" sec"<<std::endl;


    int exp_bin_num = 200000;
    double *exp_x = new double[exp_bin_num];
    double *x = new double[exp_bin_num];
    double step;
    step = rand_max/(exp_bin_num-1);
    for(i=0;i<exp_bin_num;i++)
    {
        exp_x[i] = exp(step*i);
        x[i] = step*i;
    }

    char inform[200];
    st1 = clock();
    for(i=0;i<num;i++)
    {
        j = int(rand_val[i]/step);
        exp_val_fit[i] = exp_x[j];
        if(i<20)
        {   
            std::cout<<rand_val[i]<<" "<<exp(rand_val[i])<<" "<<x[j]<<" "<<exp_x[j]<<std::endl;
            std::cout<<rand_val[i]-x[j]<<" "<<exp(rand_val[i])-exp_x[j]<<std::endl;
            std::cout<<(rand_val[i]-x[j])/rand_val[i]<<" "<<(exp(rand_val[i])-exp_x[j])/exp(rand_val[i])<<std::endl;
            std::cout<<std::endl;

        }
    }
    st2 = clock();
    std::cout<<num<<" times (double) fit exp use "<<(st2-st1)/CLOCKS_PER_SEC<<" sec"<<std::endl;
    
    double diff_min, diff_max, diff;
    diff = 0;
    diff_min = 1000;
    diff_max = -10;
    int tag;
    for(i=0;i<num;i++)
    {
        s1 = fabs(exp_val_fit[i] - exp_val[i]);
        diff += s1;
        if (s1> diff_max )
        {
            diff_max = s1;
            tag = i;
        }
        if (s1< diff_min )
        {
            diff_min = s1;
        }
    }
    std::cout<<"Diff max: "<<diff_max<<" min: "<<diff_min<<" mean: "<<diff/num<<std::endl;
    std::cout<<rand_val[tag]<<" "<<exp_val_fit[i]<<" "<< exp_val[i]<<std::endl;
    gsl_free();
    return 0;
}