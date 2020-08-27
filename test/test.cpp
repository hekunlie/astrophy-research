#include<iostream>
#include<algorithm>

int main()
{   
    double temp_x_, temp_y_;
    int ix_, iy_;
    double temp_x, temp_y;
    int ix, iy, bin_num,bin_num2;
    int i, im, im1, im2;
    bin_num = 10;
    bin_num2 = bin_num/2;

    ix_ = 8;
    iy_ = 4;
    
    temp_x = -1;
    temp_x_ = 0;
    if(temp_x >= temp_x_)
    {
        for(im=ix_; im<bin_num; im++)
        {std::cout<<im<<" "<<im+1<<std::endl;}
    }
    else
    {
        for(im=ix_; im>-1; im--)
        {std::cout<<im<<" "<<im+1<<std::endl;}
    }
    std::cin>>i;
    return 0;
}