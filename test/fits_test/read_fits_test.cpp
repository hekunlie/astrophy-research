#include<FQlib.h>
#include<hk_iolib.h>

void read_fits_new(const char *filename, const int hdu_label, float *buffer_ptr, int move)
{
    // read the i'th tile from the compressed fits file

    fitsfile *fptr, *fptr_z;	/* pointer to the FITS file, defined in fitsio.h */
	int status = 0, nfound, anynull, hdutype, hdupos=-1,hdupos_ini;
	long naxes[2],fpixel = 1, nbuffer, npixels, totpix;
    int naxis=0, bitpix=0, hud_num=0;
	double datamin, datamax, nullval = 0;

    int ii, jj;

	fits_open_file(&fptr, filename, READONLY, &status);
    
    fits_get_num_hdus(fptr, &hud_num, &status);

    fits_get_hdu_num(fptr, &hdupos_ini);

    std::cout<<hud_num<<" "<<hdupos_ini<<std::endl;
    // fits_is_compressed_image(fptr, &status);
    char name[50];
    
    for(hdupos=hdupos_ini; hdupos < hud_num; hdupos++)  /* Main loop through each extension */
    {
        fits_movabs_hdu(fptr, hdupos, &hdutype, &status);

        // fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);
        /* read the NAXIS1 and NAXIS2 keyword to get image size */
        sprintf(name,"NAXIS1");
        fits_read_key(fptr, TINT,name, &ii, NULL, &status);
        sprintf(name,"NAXIS2");

        fits_read_key(fptr, TINT, name, &jj, NULL, &status);

        // fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);
        std::cout<<hdupos<<" "<<status<<" "<<" "<<ii<<" "<<jj<<std::endl;

    }

        
	// fits_read_keys_lng(fptr_z, "NAXIS", 1, 2, naxes, &nfound, &status);		/* read the NAXIS1 and NAXIS2 keyword to get image size */
    
	// npixels = naxes[0] * naxes[1];	/* number of pixels in the image, python_arr[naxes[1], naxes[0]] */
    // std::cout<<naxes[0]<<" "<<naxes[1]<<" "<<npixels<<std::endl;

    // // buffer_ptr = new float[npixels];
	// fits_read_img(fptr_z, TFLOAT, fpixel, npixels, &nullval, buffer_ptr, &anynull, &status);

	fits_close_file(fptr, &status);
}


int main(int argc, char const *argv[])
{
    float *img = new float[100];
    char file_name[500];
    char result_name[500];
    int hdu_label, move;

    strcpy(file_name, argv[1]);
    hdu_label = atoi(argv[2]);
    move = atoi(argv[3]);

    sprintf(result_name, "!/home/hklee/work/code_test/fits_test/test.fits");

    read_fits_new(file_name, hdu_label, img, move);
    // std::cout<<"Read image"<<std::endl;
    // show_arr(img, 10, 10);
    // write_fits(result_name, img, 4094, 2046);

    return 0;
}
