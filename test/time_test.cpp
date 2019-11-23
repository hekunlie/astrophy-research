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
    double st1, st2;

    gsl_initialize(seed);

    for(i=0;i<num;i++)
    {
        rand_uniform(0.5, 100, rand_val[i]);
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
        s1 = exp(rand_val[i]);
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
    gsl_free();

    return 0;
}