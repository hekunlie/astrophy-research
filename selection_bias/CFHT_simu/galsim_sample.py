from sys import path, argv
path.append("/home/hkli/work/mylib")
path.append('D:/GitHub/astrophy-research/mylib')
import numpy
from astropy.io import fits
import galsim


stamp_size = int(argv[1])
sersic_idx = float(argv[2])
flux = float(argv[3])
e1 = float(argv[4])
e2 = float(argv[5])

noise_sig = 60

ras = [0.2, 0.3, 0.4, 0.5, 0.6, 0.8]

imgs = numpy.zeros((2*stamp_size, len(ras)*stamp_size))

noise = numpy.random.normal(0,noise_sig, (stamp_size, stamp_size))

pixel_scale = 0.187
psf = galsim.Moffat(beta=3.5, fwhm=0.7, flux=1.0, trunc=1.4)
psf_img = galsim.ImageD(stamp_size, stamp_size)
psf.drawImage(image=psf_img, scale=pixel_scale)

for i in range(len(ras)):

    gal = galsim.Sersic(scale_radius=ras[i], n=sersic_idx, trunc=4.5 * ras[i],flux=1.0)
    gal_e = gal.shear(e1=e1, e2=e2)
    gal_f = gal_e.withFlux(flux)
    gal_c = galsim.Convolve([gal_f, psf])
    img = galsim.ImageD(stamp_size, stamp_size)
    gal_c.drawImage(image=img, scale=pixel_scale)

    imgs[:stamp_size,i*stamp_size:(i+1)*stamp_size] = img.array
    imgs[stamp_size:,i*stamp_size:(i+1)*stamp_size] = img.array + noise

hdu = fits.PrimaryHDU(imgs)
psf_path ='./img.fits'
hdu.writeto(psf_path, overwrite=True)