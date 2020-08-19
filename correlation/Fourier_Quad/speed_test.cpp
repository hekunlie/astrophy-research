#include<FQlib.h>

int main()
{
    double st1, st2;
    int i,j;
    st1 = clock();
    for(i=0;i<100000;i++)
    {
        for(j=0;j<10000;j++)
        {
            if(i<j){;}
            else{;}
        }
    }
    st2 = clock();
    std::cout<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;
    return 0;
}