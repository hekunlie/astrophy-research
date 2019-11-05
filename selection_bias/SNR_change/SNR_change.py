import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
from astropy.io import fits
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
import galsim
import time
from subprocess import Popen

# generate galaxy and measure the SNR's

size = int(argv[1])
psf_r = float(argv[2])
e = float(argv[3])
btr = float(argv[4])
ra = float(argv[5])
seed = int(argv[6])
file_tag = int(argv[7])
# seed = 52405 # used in the paper

pts_source = 0

num = 11
pixel_scale = 0.187

if pts_source == 0:
    flux = numpy.array([tool_box.mag_to_flux(21.5), tool_box.mag_to_flux(22.2), tool_box.mag_to_flux(24.1)])
else:
    flux = numpy.array([tool_box.mag_to_flux(21.5), tool_box.mag_to_flux(22.5), tool_box.mag_to_flux(24.5)])

flux_num = len(flux)
noise_sig = 60

detect_thresh = 2

fq = Fourier_Quad(size, seed)
fq_p = Fourier_Quad(size, 1760)
# all the images are added by the same noise
noise = fq.draw_noise(0, noise_sig)

shear_beta = numpy.linspace(0, numpy.pi, num)
input_g = numpy.linspace(-0.06, 0.06, num)

total_path = "./imgs/"
img_path = total_path + "fits/%d/"%file_tag
os.mkdir(img_path)

if pts_source == 0:
    psf = galsim.Moffat(beta=3.5, fwhm=psf_r, flux=1.0, trunc=psf_r*3)
    psf_img = galsim.ImageD(size, size)
    psf.drawImage(image=psf_img, scale=pixel_scale)
    psf_img = psf_img.array
else:
    psf_img = fq.cre_psf(4, model="Moffat")

hdu = fits.PrimaryHDU(psf_img)
hdu.writeto(total_path+"psf.fits",overwrite=True)
idx = psf_img < psf_img.max()/2
psf_img[idx] = 0
idx = psf_img > 0
psf_img[idx] = 1
psf_radius = numpy.sqrt(psf_img.sum()/numpy.pi)
psf_quad = tool_box.get_quad(psf_img,size,psf_radius)[0]

criterion_num = 5
# the SNR_S, SNR_A, MAG, resolution factor, PK0 of the original galaxy
ori_data = numpy.zeros((flux_num*criterion_num, 4))
# flux_1: pk0, tra_SNR, 0, 0
# flux_2: pk0, tra_SNR, 0, 0
# flux_3: pk0, tra_SNR, 0, 0
# ...
# flux_1: SEX_SNR, tra_SNR, 0, 0
# flux_2: SEX_SNR, tra_SNR, 0, 0
# flux_3: SEX_SNR, tra_SNR, 0, 0
# ...
# flux_1: SNR_AUTO, tra_SNR, flux_auto, flux_err
# ...
# MAG, tra_SNR, 0, 0
# ...
# resolution factor,0,0
# ...

data = numpy.zeros((flux_num*criterion_num, num))
# each column:
# pk0 of flux_1
# pk0 of flux_2
# pk0 of flux_3
# ...
# SEX_SNR of flux_1
# SEX_SNR of flux_2
# SEX_SNR of flux_3
# ...
# SNR_AUTO
# ...
# SEX_MAG
# ...
# resolution factor
# ...

ext_data = numpy.zeros((3, num)) # for SNR_A

pool = []
pts_num = 100
rand_pts = fq_p.ran_pts(num=pts_num, radius=10, ellip=0.8)

print("Simulate")

for k in range(flux_num):
    # draw the pre-sheared galaxy
    # galsim galaxy
    if pts_source == 0:
        # rng = galsim.BaseDeviate(12300000)
        # knot = galsim.randwalk.RandomWalk(npoints=60, half_light_radius=ra, flux=1, rng=rng)
        bulge = galsim.Sersic(half_light_radius=ra, n=4, trunc=4.5 * ra)  # be careful
        disk = galsim.Sersic(scale_radius=ra, n=1, trunc=4.5 * ra)  # be careful
        gal = bulge * btr + disk * (1 - btr) #+ ktr*knot
        gal = gal.shear(e1=e, e2=0)#beta=0.*galsim.degrees)
        gal_f = gal.withFlux(flux[k])
        #
        # rng = galsim.BaseDeviate(12300000)
        # gal = galsim.randwalk.RandomWalk(npoints=200, half_light_radius=ra, flux=flux[k], rng=rng)
        # gal_f = gal.shear(e1=e, e2=0)

        gal_c = galsim.Convolve([gal_f, psf])
        img_0 = galsim.ImageD(size, size)
        gal_c.drawImage(image=img_0, scale=pixel_scale)
        ori_gal_img = img_0.array + noise
    # random walk galaxy
    else:
        img_0 = fq.convolve_psf(rand_pts, 4, flux[k]/pts_num, "Moffat")
        ori_gal_img = img_0 + noise

    path = img_path + "gal_pre_shear_%d.fits" %k
    hdu = fits.PrimaryHDU(ori_gal_img)
    hdu.writeto(path, overwrite=True)
    pool.append(ori_gal_img)

    # the sheared galaxy
    for i in range(num):
        if pts_source == 0:
            gal_s = gal_f.shear(g1=input_g[i], g2=0)#beta=shear_beta[i]*galsim.radians)
            gal_s_c = galsim.Convolve([gal_s, psf])
            img_s = galsim.ImageD(size, size)
            gal_s_c.drawImage(image=img_s, scale=pixel_scale)
            shear_gal_img = img_s.array + noise
        else:
            rand_pts_s = fq.shear(rand_pts, g1=input_g[i], g2=0)
            img_s = fq.convolve_psf(rand_pts_s, 4, flux[k]/pts_num, "Moffat")
            shear_gal_img = img_s + noise

        pool.append(shear_gal_img)
        path = img_path + "gal_shear_%d_%d.fits"%(k,i)
        hdu = fits.PrimaryHDU(shear_gal_img)
        hdu.writeto(path, overwrite=True)

big_img = fq.stack(pool, num+1)
hdu = fits.PrimaryHDU(big_img)
hdu.writeto(img_path+"total_img.fits",overwrite=True)


############################## run SEX ###############################
print("Run SEX")

for i in range(flux_num):

    # read the original galaxy
    ori_gal_path = img_path + "gal_pre_shear_%d.fits"%i
    ori_gal = fits.open(ori_gal_path)[0].data

    # the pk0
    ori_pow = fq.pow_spec(ori_gal)
    ori_pk0 = numpy.sqrt(ori_pow[int(size/2), int(size/2)])/size/noise_sig
    print(i)
    # detect the source
    detect_0 = tool_box.stamp_detector(ori_gal, size, size, 5, 5.5, detect_thresh*noise_sig)

    if detect_0:
        # the traditional SNR
        mask_0 = detect_0[0]
        snr_tradi_0 = numpy.sum(mask_0*ori_gal)/numpy.sqrt(mask_0.sum())/noise_sig

        # measure the SNR_S, SNR_A, MAG of original galaxy
        ori_gal_cat_name = img_path + "gal_pre_shear_%d.cat"%i
        cmd = "sex %s -CATALOG_NAME %s" % (ori_gal_path, ori_gal_cat_name)
        sex = Popen(cmd, shell=True)
        sex.wait()
        try:
            ori_gal_cata = numpy.loadtxt(ori_gal_cat_name)
            cat_row = ori_gal_cata.shape
            ori_sex_snr = ori_gal_cata[0]
            ori_flux_auto = ori_gal_cata[1]
            ori_flux_err = ori_gal_cata[2]
            ori_snr_auto = ori_flux_auto / ori_flux_err
            ori_sex_mag = ori_gal_cata[3]

            ori_radius = numpy.sqrt(ori_gal_cata[4]/numpy.pi)
            ori_gal_quad = tool_box.get_quad(ori_gal,size, ori_radius)[0]
            ori_r_factor = ori_gal_quad #1 - psf_quad/(ori_gal_quad + psf_quad)
            print("%d: "%i, ori_gal_cata)
        except:
            print("%d: Not found"%i)
            ori_sex_snr, ori_snr_auto, ori_sex_mag, ori_flux_auto, ori_flux_err, ori_r_factor = -99, -99, -99, -99, -99,-99

        ori_data[i, 0:2] = ori_pk0, snr_tradi_0
        ori_data[i + flux_num, 0:2] = ori_sex_snr, snr_tradi_0
        ori_data[i + int(flux_num*2)] = ori_snr_auto, snr_tradi_0, ori_flux_auto, ori_flux_err # save additional parameters
        ori_data[i + int(flux_num*3),0:2] = ori_sex_mag, snr_tradi_0
        ori_data[i + int(flux_num*4),0:2] = ori_r_factor, snr_tradi_0

        # measure each sheared galaxy
        for k in range(num):
            gals_path = img_path + 'gal_shear_%d_%d.fits' %(i, k)
            gal_img = fits.open(gals_path)[0].data

            fits_name = img_path + "gal_shear_%d_%d.fits" %(i, k)
            cat_name = img_path + "gal_shear_%d_%d.cat" %(i, k)
            cmd = "sex %s -CATALOG_NAME %s" % (fits_name, cat_name)
            sex = Popen(cmd, shell=True)
            sex.wait()

            try:
                cata_data = numpy.loadtxt(cat_name)
                sex_snr = cata_data[0]
                flux_auto = cata_data[1]
                flux_err = cata_data[2]
                snr_auto = flux_auto / flux_err
                sex_mag = cata_data[3]

                gal_radius = numpy.sqrt(cata_data[4] / numpy.pi)
                gal_quad = tool_box.get_quad(gal_img, size, gal_radius)[0]
                r_factor = gal_quad #1 - psf_quad / (gal_quad + psf_quad)

            except:
                sex_snr, snr_auto, sex_mag, flux_auto, flux_err,r_factor = -99,-99,-99,-99,-99,-99

            gal_pow = fq.pow_spec(gal_img)
            pk0 = numpy.sqrt(gal_pow[int(size / 2), int(size / 2)])/size/noise_sig

            data[i, k] = pk0
            data[i + flux_num, k] = sex_snr
            data[i + int(flux_num*2), k] = snr_auto
            ext_data[0,k] = snr_auto
            ext_data[1,k] = flux_auto
            ext_data[2,k] = flux_err
            data[i + int(flux_num*3), k] = sex_mag
            data[i + int(flux_num*4), k] = r_factor

result_path = img_path + "result.npz"
numpy.savez(result_path, input_g, ori_data, data, ext_data)

pic_path = total_path + "pic/change_%d.png"%file_tag
cmd = "python SNR_change_plot.py %d %d %s %s"%(5, flux_num, result_path, pic_path)
os.system(cmd)
