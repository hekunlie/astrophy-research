#include<FQlib.h>
#include<hk_mpi.h>
#include<hk_iolib.h>

void coord_rotation(const double*xy, const int pts_num, const double theta, double *xy_r)
{
    double sin_theta, cos_theta;
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    for (int i = 0; i < pts_num; i++)
    {
        xy_r[i] = cos_theta * xy[i] - sin_theta*xy[i + pts_num];
        xy_r[i + pts_num] = sin_theta * xy[i] + cos_theta*xy[i + pts_num];
        // std::cout<<xy[i]<<" "<<xy[i+pts_num]<<" "<<xy_r[i]<<" "<<xy_r[i+pts_num]<<std::endl;
    }
}


int main(int argc, char*argv[])
{   
    int pts_num;
    int i,j,k;
    int size, img_cent;
    int seed_pts;
    double theta;
    char img_path[200];

    pts_num = atoi(argv[1]);
    seed_pts = atoi(argv[2]);
    size = 50;
    img_cent = size*0.5 - 0.5;

    double *point = new double[2*pts_num]{};
    double *point_r = new double[2*pts_num]{};
    double *img = new double[size*size]{};
    

    gsl_initialize(seed_pts,1);

    create_points(point, pts_num, 10, 1, rng1);


    for(i=0;i<11;i++)
    {   
        theta = i*Pi/10*2;
        coord_rotation(point, pts_num, theta, point_r);

        convolve(img, point_r, 1, size, img_cent, pts_num, 0, 4, 0, 0, 2);

        sprintf(img_path, "!img_%d.fits", i);
        write_fits(img_path, img, size, size);
    }

    return 0;
}