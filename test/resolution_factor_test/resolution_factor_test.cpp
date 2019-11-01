#include <FQlib.h>

void get_quad_(const double *img, const int img_size, const double weight_sigma_sq, double &quad_size,int tag, char*inform)
{
    // calculate the gaussian-weighted quadrupole of a galaxy or PSF image 
    double temp_quad, temp_norm;
    int i,j,m;
    double ry_sq, r_sq;
    double cen, wei_coeff,wei_img;
    
    // the image center
    cen = img_size*0.5-0.5;

    wei_coeff = 0.5/weight_sigma_sq;

    temp_quad = 0;
    temp_norm = 0;
    for(i=0; i<img_size; i++)
    {   
        ry_sq = (i - cen)*(i-cen);
        m = i*img_size;
        for(j=0; j<img_size; j++)
        {
            r_sq = (j - cen)*(j-cen) + ry_sq;

            wei_img = exp(-r_sq*wei_coeff)*img[m+j];
   
            temp_quad += wei_img*r_sq;
            temp_norm += wei_img;
        }
    }
    tag = 1;
    if(temp_quad< 0 or temp_norm<=0)
    {
        sprintf(inform,"%.5f, %.5f",temp_quad, temp_norm);
        tag = -1;
    }
    
    quad_size = temp_quad/temp_norm;
}

int main(int argc, char **argv)
{
    char img_path[200], data_path[200], set_name[30];
    double *img, *stamp;
    double *sex_data, *data;
    int i,j,k,tag;
    int nx, ny, size;
    int sex_row, sex_col;

    double st1, st2;

    double gal_quad;
    double eff_radius_sq;

    int quad_tag;
    char inform[40];

    size = 64;
    nx = 100;
    ny = 100;

    img = new double[nx*ny*size*size];
    stamp = new double[size*size];

    sex_row = 10000000;
    sex_col = 4;
    sex_data = new double[sex_row*sex_col];

    data = new double[nx*ny]{};

    sprintf(img_path,"/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/pts_dimmer/1/gal_chip_0000.fits");
    read_fits(img_path, img);

    sprintf(data_path, "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/pts_dimmer/result/data/sex2_1.5/sex_1.hdf5");
    sprintf(set_name,"/data");
    read_h5(data_path, set_name, sex_data);
    st1 = clock();

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            tag = i*nx + j;
            segment(img, stamp, tag, size, nx, ny);
            // \pi * r^2 = pixel number
            eff_radius_sq = sex_data[tag*sex_col + 3]/Pi;           
    
            if (eff_radius_sq > 0)            
            {
                get_quad(stamp, size, eff_radius_sq, gal_quad);
                data[tag] = gal_quad;
            }

            if(tag < 20)
            {
                std::cout<<sex_data[tag*sex_col + 3]<<" "<<eff_radius_sq<<" "<<data[tag]<<std::endl;
            }

        }

    }

    sprintf(data_path,"quad_size.hdf5");
    sprintf(set_name,"/data");
    write_h5(data_path, set_name, data,nx*ny,1,true);
    st2 = clock();
    std::cout<<(st2-st1)/CLOCKS_PER_SEC<<std::endl;

    delete[] sex_data;
    delete[] data;
    delete[] img;
    delete[] stamp;

    return 0;
}