from Fourier_Quad import Fourier_Quad
import lsstetc
import numpy
import matplotlib.pyplot as plt
import tool_box

a = [1,2,3,4,5,6,7,8]

b = numpy.random.randint(10,20, 8)/20
plt.errorbar(a,a,b,capsize=2,c='coral')
plt.plot([1,8],[2,2],c='grey')
plt.show()
size = 52
num = 20000
prop = lsstetc.ETC(band='r', pixel_scale=0.2, stamp_size=size, nvisits=180)
mags = tool_box.mags_mock(num, 21, 26.5)
plt.hist(mags, 100)
plt.show()
plt.close()
noise_sig = prop.sigma_sky
print(noise_sig)
fq = Fourier_Quad(size, numpy.random.randint(0, 10000000, 1))
psf =fq.cre_psf(4, "Moffat")
plt.imshow(psf)
plt.show()
print("PSF:", numpy.sum(psf))
measures = numpy.zeros((num, 5))
rim = fq.border(3)
n = numpy.sum(rim)
for i in range(num):
    flux = prop.flux(mags[i]) / 45
    p = fq.ran_pos(45, 9)
    gal_final = fq.convolve_psf(p, 4, flux, 'Moffat') + fq.draw_noise(0, noise_sig)
    obj, oflux, signalsq, peak, flag = tool_box.stamp_detector(gal_final, noise_sig*2, size, size)
    snr = numpy.sqrt(signalsq) / noise_sig
    ori_snr = oflux/numpy.sqrt(len(obj))/noise_sig

    gpow = fq.pow_spec(gal_final)
    if numpy.max(gpow) == gpow[int(size/2), int(size/2)]:
        signal = gpow[int(size/2), int(size/2)]
    else:
        signal = numpy.sum(gpow[int(size/2-1):int(size/2+2), int(size/2-1):int(size/2+2)])/9
    noise_level = numpy.sum(rim*gpow)/n
    fsnr = numpy.sqrt(signal/noise_level)

    measures[i] = snr, ori_snr, peak/noise_sig, fsnr, oflux/noise_sig
    if i< 50:
        print(numpy.max(gal_final))
        # idx = gal_final < 2*noise_sig
        # gal_final[idx] = 0
        # plt.imshow(gal_final)
        # plt.show()

snr, ori_snr, peak, fsnr, oflux = measures[:,0], measures[:,1], measures[:,2], measures[:,3], measures[:,4]
print(numpy.max(snr), numpy.max(ori_snr), numpy.max(peak), numpy.max(fsnr), numpy.max(oflux))
plt.hist(snr[snr<40], 50)
plt.show()
plt.close()
plt.hist(ori_snr[ori_snr<40], 50)
plt.show()
plt.close()
plt.hist(peak[peak<30], 100)
plt.show()
plt.close()
plt.hist(fsnr[fsnr<20], 50)
plt.show()
plt.close()
plt.hist(oflux[oflux<100], 50)
plt.show()
plt.close()