    The purpose of the code is to measure cosmic shear directly from the original images of individual exposures using the method called Fourier_Quad [Zhang, Zhang, Luo, ApJ, 834, 8(2017)]. Fourier_Quad has a number of advantages, including: 1) free of assumptions about the galaxy or PSF morphology; 2) rigorous removal of noise (background and Poisson) contamination; 3) extremely fast image processing speed (less than 1e-2 CPU*sec/Galaxy); 4) optimal shear statistics, enabling the statistical error of shear measurement to approach the theoretical lower limit; 5) shear results are immune to misidentification of stars as galaxies; 6) accurate to the second order in shear/convergence.
    The code provided in this package includes subroutines taking care of the following image processing procedures:

1) Background removal and its rms;
2) Defects identification (cosmic rays, bad pixels, etc.);
3) Astrometric calibration;
4) PSF reconstruction using bright stars;
5) Shear estimation with Fourier_Quad;
6) Shear accuracy test with field distortion.

To run the code, one has to provide a "source.list" file, which includes all the files required for shear measurement. For example, to compile and run the code, in the directory containing the fortran files, one can do the following:

mpif77 *.f -o main -mcmodel=medium -lcfitsio
mpiexec -n 31 ./main /mw/w1234/original/source.list

Number "31" here specifies the number of cores user tends to use for the MPI process. The "source.list" file should be prepared in the following way:

           1           1           1
/mw/w1234/original/w1m0m0/science/819719p_10.fits                               
/mw/w1234/original/gaia_cat/w1.cat                                              
/mw/w1234/w1m0m0/science/cat/w1m0m0.dat                                         
/mw/w1234/original/w1m0m0                                                       
           1           1           2
/mw/w1234/original/w1m0m0/science/819719p_11.fits                               
/mw/w1234/original/gaia_cat/w1.cat                                              
/mw/w1234/w1m0m0/science/cat/w1m0m0.dat                                         
/mw/w1234/original/w1m0m0                             

... ... ... ...

The three numbers in the first row of each data block correspond to the ID numbers of the field, exposure, and chip of each image. In cases when the  images are all independent (e.g., simulated images), user only needs to roll the field number, fixing both the exposure ID and the chip ID at one.

The second row specifies the path and the name of the image file; the third row gives the name of the standard catalog file (e.g., from GAIA) for astrometric calibration; the fourth row refers to the name of the source catalog (though the code itself can generate a source catalog); the fifth rwo refers to the parent path, in which all results are stored. The astrometrical parameters are stored in the subdirectory called "astrometry", and the shear catalog as well as related results are stored in subdirectory called "result", both of which can be created through shell scripts, such as those in "run_code.sh". 

Most of the key parameters of the code is specified in the "para.inc" file. The meanings of the parameters are provided below:

£¨npx£¬npy£©£º physical size of the 2D array that can contain the image from the fits file;

ASTROMETRY_READY (0,1,2): "0" means the code performs astrometric calibration; "1" means the code performs shear measurement based on the astrometric parameters that are already calculated; "2" means the astrometric parameters are trivial, i.e., only CRPIX, CRVAL, and CD's are used (usually for simulated data);

astrometry_search_radius: the maximal dispersion of displacement between the source position on the image and the position inferred by the standard catalog;

astrometry_shift_ratio: the maximal relative displacement between the source position on the image and the position inferred by the standard catalog (with respect to the size of the long side of the chip);

ext_cat (0,1): "0" means the code generates the source catalog by itself;
"1" means external source catalog is used for shear measurement;

ns: the size of the galaxy and star stamp;

chip_margin£ºthe distance (in pixel) that the edge of the source stamps must be away from the edge of the chip;

flag_thresh: The minimal distance (in pixel) that a valid source should be away from a defect or another neighboring source;

len_g: the side length £¨in terms of stamp size "ns") of the large fits file that contains source stamps;

len_s: the side length £¨in terms of stamp size "ns") of the large fits file that contains star stamps;

nstar_max: the maximal number of stars in a chip;

ngal_max: the maximal number of sources in a chip;

npara: the maximal number of parameters that one needs to store in the source catalog (including the shear estimators);

NMAX_CHIP: the maximal number of chips in an exposure;

npd: the number of free parameters in the polynomial function for astrometrical calibration;

blocksize: the side length of each subregion that is used for background evaluation; 

source_thresh: threshold for defining source pixels;

area_thresh£º minimal number of pixels that a source should contain;

nstar_min: minimal number of stars a chip should contain for performing PSF reconstruction;

psf_local(0,1): "0" means perform exposure-wise polynomial global PSF fitting; "1" means perform chip-wise polynomial local PSF fitting;

psf_order: the order of polynomial function for the global fitting;

nstar_thresh: the minimal number of stars that is required for global PSF fitting;

npo: the total number of free parameters the global PSF fitting contains;

npl: the total number of free parameters the local PSF fitting contains;
 
gal_smooth(0,1): "1" means smoothing the galaxy power stamp, "0" means not;

star_smooth(0,1): "1" means smoothing the star power stamp, "0" means not;

SNR_thresh_star: the minimal SNR of the stars that are used for PSF reconstruction;

saturation_thresh: the saturation level of the fits file. 
	
NMAX_EXPO: the maximal number of exposures per field;

NMAX_FIELD: the maximal number of fields;

ngfieldmax: the maximal number of galaxies per field (the images of each galaxy on different exposures are counted as different galaxies here).

