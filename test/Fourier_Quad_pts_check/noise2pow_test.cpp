#include<FQlib.h>
#include<hk_iolib.h>

void noise2pow(double *pow_img, const int size, const double sigma, gsl_rng *gsl_rand_rng)
{
    int i,j, kx1, ky1, kx2, ky2;
    int cent;
    double noise;
    cent = size*0.5;

    for(i=1;i<cent;i++) // y
    {
        ky1 = (cent + i)*size;
        ky2 = (cent - i)*size;

        // ky = size/2
        rand_gauss(sigma, 0, noise, gsl_rand_rng);
        pow_img[cent*size + cent+i] += noise;
        pow_img[cent*size + cent- i] += noise;

        // kx = size/2
        rand_gauss(sigma, 0, noise, gsl_rand_rng);
        pow_img[ky1 + cent] += noise;
        pow_img[ky2 + cent] += noise;

        // ky = 0, kx = [1, size]
        rand_gauss(sigma, 0, noise, gsl_rand_rng);
        pow_img[cent + i] += noise;
        pow_img[cent - i] += noise;

        // kx=0, ky = [1, size]
        rand_gauss(sigma, 0, noise, gsl_rand_rng);
        pow_img[ky1] += noise;
        pow_img[ky2] += noise;

        for(j=1;j<cent;j++) // x
        {
            kx1 = cent + j;
            kx2 = cent - j;

            rand_gauss(sigma, 0, noise, gsl_rand_rng);
            pow_img[ky1 + kx1] += noise;
            pow_img[ky2 + kx2] += noise;
            
            rand_gauss(sigma, 0, noise, gsl_rand_rng);
            pow_img[ky1 + kx2] += noise;
            pow_img[ky2 + kx1] += noise;
        }
    }

    rand_gauss(sigma, 0, noise, gsl_rand_rng);
    pow_img[0] += noise; // [0,0]
    rand_gauss(sigma, 0, noise, gsl_rand_rng);
    pow_img[cent] += noise; // [0, size/2]
    rand_gauss(sigma, 0, noise, gsl_rand_rng);
    pow_img[cent*size] += noise; // [size/2, 0]
    rand_gauss(sigma, 0, noise, gsl_rand_rng);
    pow_img[cent*size + cent] += noise; // [size/2, size/2]
}

int main(int argc, char **argv)
{
    int i,j,k;
    int size;
    double *img;

    size = atoi(argv[1]);
    img = new double[size*size]{};

    gsl_initialize(123, 0);
    noise2pow(img, size, 10, rng0);
    gsl_free(0);

    char data_path[100];
    sprintf(data_path, "!pow_noise.fits");
    write_fits(data_path, img, size, size);
    delete[] img;


    return 0;
}