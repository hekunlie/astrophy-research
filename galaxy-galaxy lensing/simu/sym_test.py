from sys import path,argv
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import gglensing_tool
import numpy
import h5py
import galsim
from Fourier_Quad import Fourier_Quad
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units
import time
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

mg_bin_num = int(argv[1])

data_path = "/mnt/perc/hklee/Galaxy_Galaxy_lensing_test/cata/background/single_shear_test"
data_type = argv[2]
shear_tag = 1
non_shear_tag = 0
num_s = 10000000

data_path1 = data_path + "/data/sheared_data_%d_%s.hdf5"%(shear_tag,data_type)
para_path1 = data_path + "/sheared_para_%d.hdf5"%(shear_tag)

data_path2 = data_path + "/data/non_sheared_data_%d_%s.hdf5"%(non_shear_tag,data_type)
para_path2 = data_path + "/non_sheared_para_%d.hdf5"%(non_shear_tag)

data_path3 = data_path + "/data/non_sheared_data_%d_%s.hdf5"%(non_shear_tag + 1, data_type)
para_path3 = data_path + "/non_sheared_para_%d.hdf5"%(non_shear_tag + 1,)

mg1 = numpy.zeros((int(2*num_s),),dtype=numpy.float32)
mnu1 = numpy.zeros((int(2*num_s),),dtype=numpy.float32)


# sheared galaxies
h5f = h5py.File(data_path1,"r")
mg1[:num_s] = h5f["/data"][()][:,0]
mnu1[:num_s] = h5f["/data"][()][:,2] + h5f["/data"][()][:,3]
h5f.close()
# non-sheared galaxies
h5f = h5py.File(data_path2,"r")
mg1[num_s:] = h5f["/data"][()][:,0]
mnu1[num_s:] = h5f["/data"][()][:,2] + h5f["/data"][()][:,3]
h5f.close()

h5f = h5py.File(data_path3,"r")
mg1_corr = h5f["/data"][()][:,0]
mnu1_corr = h5f["/data"][()][:,2] + h5f["/data"][()][:,3]
h5f.close()


G1_bin = tool_box.set_bin(mg1[:num_s], mg_bin_num, 1000)
G1_hist_bin = gglensing_tool.set_bin(mg1[:num_s], 2500, 1.001,"log")
NU1_hist_bin = gglensing_tool.set_bin(mnu1[:num_s], 2500, 1.001,"log")


dilute_case = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
dilute_num = len(dilute_case)
task_list = [i for i in range(dilute_num)]

task_list_sub = tool_box.alloc(task_list, numprocs)[rank]


itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes = 4*dilute_num*itemsize
else:
    nbytes = 0

# on rank 0 of comm, create the contiguous shared block
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)

result_min = numpy.ndarray(buffer=buf1, dtype='d', shape=(4, dilute_num))


comm.Barrier()

for tag in task_list_sub:
    num_non = int(num_s * dilute_case[tag])
    num_total = num_s + num_non
    num_corr = num_non

    img = Image_Plot(xpad=0.25, ypad=0.25)
    img.subplots(1, 2)
    t1 = time.time()
    gh_mix, gh_mix_sig, coeff_mix, asym_mix = gglensing_tool.find_shear_grid(mg1[:num_total], mnu1[:num_total], G1_bin,
                                                                             G1_hist_bin, NU1_hist_bin, chisq_gap=50,
                                                                             max_iters=50, fit_num=20, dg=0.002,
                                                                             ax=img.axs[0][0])[:4]
    t2 = time.time()

    chisq_min_mix = coeff_mix[0] - coeff_mix[1] ** 2 / 4 / coeff_mix[2]
    result_min[0,tag] = chisq_min_mix
    result_min[1,tag] = asym_mix

    print("%.5f(%.5f). asym: %.3e. %d source + %d contaminations. %.2f sec" % (
    gh_mix, gh_mix_sig, asym_mix, num_s, num_non, t2 - t1))

    # correction
    t1 = time.time()
    ghc, ghc_sig, coeffc, asymc = gglensing_tool.find_shear_grid_corr_new(mg1[:num_total],
                                                                          mnu1[:num_total],
                                                                          mg1_corr[:num_corr],
                                                                          mnu1_corr[:num_corr],
                                                                          G1_bin, G1_hist_bin, NU1_hist_bin,
                                                                          chisq_gap=50, max_iters=50, fit_num=20,
                                                                          dg=0.002, ax=img.axs[0][1])[:4]
    t2 = time.time()

    chisq_min_c = coeffc[0] - coeffc[1] ** 2 / 4 / coeffc[2]
    result_min[2,tag] = chisq_min_c
    result_min[3,tag] = asymc

    for m in range(2):
        img.axis_sci_ticklabel(0, m, 1)
    print("%.5f(%.5f). asym: %.3e. %d source + %d contaminations + %d corrections. %.2f sec" % (
    ghc, ghc_sig, asymc, num_s, num_non, num_corr, t2 - t1))

    img.save_img(data_path + "/result/%d/%s_%.2f_dilute_%d_mg_bins.png"%(mg_bin_num, data_type, dilute_case[tag], mg_bin_num))
    # img.show_img()
    img.close_img()

comm.Barrier()
if rank == 0:
    numpy.savez(data_path + "/result/%d/%s_min_change_%d.npz"%(mg_bin_num, data_type, mg_bin_num), result_min)
    img = Image_Plot(xpad=0.25)
    img.subplots(1, 2)
    img.axs[0][0].plot(dilute_case, result_min[0], label="$\chi^2_{min}$", marker="s")
    img.axs[0][0].plot(dilute_case, result_min[2], label="$\chi^2_{min}$_corr",marker="o")
    img.axs[0][1].plot(dilute_case, result_min[1], label="Asym", marker="s")
    img.axs[0][1].plot(dilute_case, result_min[3], label="Asym_corr",marker="o")

    img.axis_sci_ticklabel(0, 1, 0)
    img.set_label(0, 0, 0, "$\chi^2$")
    img.set_label(0, 0, 1, "dilution ratio")
    img.set_label(0, 1, 0, "Asym")
    img.set_label(0, 1, 1, "dilution ratio")
    img.axs[0][0].legend()
    img.axs[0][1].legend()
    # img.show_img()
    img.save_img(data_path + "/result/%d/%s_min_change_%d.png"%(mg_bin_num, data_type, mg_bin_num))
comm.Barrier()