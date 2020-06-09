from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
path.append("/home/hklee/work/mylib")
path.append("/home/hkli/work/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
from astropy.io import fits
import h5py
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

psf_type = "Moffat"
psf_flux = 1
psf_scale = 4
stamp_size = 44
seed = 141432

noise_sig = 60

pts_num = 40
max_radius = 7
gal_flux = 7000

g1, g2 = 0, 0

total_num = 10000

itemsize = MPI.DOUBLE.Get_size()

if rank == 0:
    # bytes for 10 double elements
    nbytes = 6*numprocs*total_num*5*itemsize
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
shear_est_buffer = numpy.ndarray(buffer=buf1, dtype='d', shape=(numprocs*total_num, 5*6))

fq = Fourier_Quad(stamp_size, seed)
fq_n = Fourier_Quad(stamp_size, seed + rank)

psf_img = fq.cre_psf(psf_scale, psf_flux, psf_type)
psf_pow = fq.pow_spec(psf_img)
fq_n.get_radius(psf_pow, 2)


shear_est_nf = numpy.zeros((total_num, 5))
shear_est_n = numpy.zeros((total_num, 5))
shear_est_ct = numpy.zeros((total_num, 5))
shear_est_ct_est = numpy.zeros((total_num, 5))
shear_est_pure_ct_est = numpy.zeros((total_num, 5))
shear_est_nr = numpy.zeros((total_num, 5))

dt = numpy.linspace(-numpy.pi, numpy.pi, total_num)

pts = fq.ran_pts(pts_num, max_radius)

gal_img_nf = fits.open("./gal_%d.fits"%rank)[0].data

# gal_img_nf = fq_n.convolve_psf(pts, psf_scale, gal_flux / pts_num, psf_type)
gal_pow_nf = fq_n.pow_spec(gal_img_nf)

for i in range(total_num):

    # noise free galaxy image
    # pts_r = fq_n.rotate(pts, dt[i])

    # noise image
    noise = fq_n.draw_noise(0, noise_sig)
    noise_pow = fq_n.pow_spec(noise)
    noise_new = fq_n.draw_noise(0, noise_sig)
    noise_new_pow = fq_n.pow_spec(noise_new)
    # noise power residual image
    noise_pow_residual = noise_pow - noise_new_pow
    # noise power cross term image
    noise_pow_ct = fq_n.pow_spec(noise + noise_new) - noise_pow - noise_new_pow

    # noisy galaxy image
    gal_img_n = gal_img_nf + noise
    gal_pow_n = fq_n.pow_spec(gal_img_n)
    gal_pow_dn = gal_pow_n - noise_new_pow

    # galaxy-noise power cross term image
    ct_pow = gal_pow_n - gal_pow_nf - noise_pow
    # galaxy-noise power cross term image est
    gal_img_new = gal_img_n + noise_new
    gal_img_new_pow = fq_n.pow_spec(gal_img_new)
    ct_pow_est = gal_img_new_pow - gal_pow_n - noise_new_pow
    # pure galaxy-noise power cross term image est
    pure_ct_pow_est = ct_pow_est - noise_pow_ct

    # estimate
    # noise free power
    shear_est_nf[i] = fq_n.shear_est(gal_pow_nf, psf_pow, F=True)

    # noisy power
    shear_est_n[i] = fq_n.shear_est(gal_pow_dn, psf_pow, F=True)

    # noise residual power
    shear_est_nr[i] = fq_n.shear_est(noise_pow_residual, psf_pow, F=True)

    # cross term
    shear_est_ct[i] = fq_n.shear_est(ct_pow, psf_pow, F=True)

    # est cross term
    shear_est_ct_est[i] = fq_n.shear_est(ct_pow_est, psf_pow, F=True)

    # pure est cross term
    shear_est_pure_ct_est[i] = fq_n.shear_est(pure_ct_pow_est, psf_pow, F=True)

    if rank == 0 and i == 0:
        img_buffer = numpy.row_stack((gal_img_nf,gal_img_n))
        hdu = fits.PrimaryHDU(img_buffer)
        hdu.writeto("./sample.fits", overwrite=True)


shear_est_buffer[rank*total_num:(rank+1)*total_num, :5] = shear_est_nf
shear_est_buffer[rank*total_num:(rank+1)*total_num, 5:10] = shear_est_nr
shear_est_buffer[rank*total_num:(rank+1)*total_num, 10:15] = shear_est_n
shear_est_buffer[rank*total_num:(rank+1)*total_num, 15:20] = shear_est_ct
shear_est_buffer[rank*total_num:(rank+1)*total_num, 20:25] = shear_est_ct_est
shear_est_buffer[rank*total_num:(rank+1)*total_num, 25:30] = shear_est_pure_ct_est
comm.Barrier()


if rank == 0:
    h5f = h5py.File("./sample_data.hdf5", "w")
    h5f["/noise_free"] = shear_est_buffer[:,:5]
    h5f["/noise_residual"] = shear_est_buffer[:,5:10]
    h5f["/noisy"] = shear_est_buffer[:,10:15]
    h5f["/cross_term"] = shear_est_buffer[:,15:20]
    h5f["/cross_term_est"] = shear_est_buffer[:,20:25]
    h5f["/pure_cross_term_est"] = shear_est_buffer[:,25:30]
    h5f.close()

comm.Barrier()
