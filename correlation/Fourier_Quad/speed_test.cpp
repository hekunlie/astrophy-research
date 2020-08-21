#include<FQlib.h>

int main()
{
    double st1, st2;
    int i,j, k;
    st1 = clock();
    k = 0;
    for(i=0;i<100000;i++)
    {
        for(j=0;j<10000;j++)
        {
            k++;
            if(k == 20000)
            {k=0;}
        }
    }
    st2 = clock();
    std::cout<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;

    k = 0;
    for(i=0;i<100000;i++)
    {
        for(j=0;j<10000;j++)
        {
            k = k%20000;
            k++;
        }
    }
    st2 = clock();
    std::cout<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;

    return 0;
}