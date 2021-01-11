#include<GG_lensing_functions_expo_wise.h>
#include<iostream>
#include<cmath>

int main(int argc, char *argv[])
{
    ggl_data_info data_info;
    data_info.jack_num = atoi(argv[1]);
    data_info.mg_bin_num = 10;
    data_info.pdf_guess_num = 100;
    data_info.signal_pts_num = 50;

    ggl_initialize(&data_info);

    std::cout<<data_info.jack_num<<" "<<data_info.pdf_guess_num<<" "
    <<data_info.signal_pts_num<<" "<<data_info.signal_chi_len<<" "<<data_info.chi_len<<std::endl;
    
    for( long i=0;i<10000000;i++)
    {atan2(tan(i*i*i*i*i*i),cos(i*i*i*i))* cos(i*i*i*i)*tan(i*i*i*i*i*i);}
    return 0;
}