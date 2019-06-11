import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import tool_box
from plot_tool import Image_Plot
import h5py
from sys import argv
import numpy
import time
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

fore_source_name = argv[1]

area_num = len(argv) - 2

if area_num <= 0:
    print("Need the area NO.")
    exit(0)

area_ids = [int(argv[i]) for i in range(2, len(argv))]
print(area_ids, area_num)

data_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/result/%s/"%fore_source_name
for tag, area_id in enumerate(area_ids):
    if rank == 0:
        if not os.path.exists(data_path + "pic/"):
            os.mkdir(data_path + "pic/")
        if not os.path.exists(data_path + "pic/w_%d"%area_id):
            os.mkdir(data_path + "pic/w_%d"%area_id)
        if not os.path.exists(data_path + "pic/total"):
            os.mkdir(data_path + "pic/total")
    comm.Barrier()

    h5f_path = data_path + "w_%d/radius_%d.hdf5"%(area_id, rank)
    h5f = h5py.File(h5f_path,"r")

    # the count in each bins
    # the tangential shear
    temp_chi_tan = h5f["/chi_tan"].value
    temp_chi_cross = h5f["/chi_cross"].value
    # the 'shear*critical_surface_density'
    temp_chi_crit_tan = h5f["/chi_crit_tan"].value
    temp_chi_crit_cross = h5f["/chi_crit_cross"].value

    if tag == 0:
        gh = h5f["/gh"].value[:,0]
        gh_crit = h5f["/gh_crit"].value[:,0]

        gh_num = gh.shape[0]
        gh_crit_num = gh_crit.shape[0]

        chi_tan = numpy.zeros_like(temp_chi_tan) + temp_chi_tan
        chi_cross = numpy.zeros_like(temp_chi_cross) + temp_chi_cross

        chi_crit_tan = numpy.zeros_like(temp_chi_crit_tan) + temp_chi_crit_tan
        chi_crit_cross = numpy.zeros_like(temp_chi_crit_cross) + temp_chi_crit_cross
    else:
        chi_tan += temp_chi_tan
        chi_cross += temp_chi_cross

        chi_crit_tan += temp_chi_crit_tan
        chi_crit_cross += temp_chi_crit_cross

    mg_bin_num = temp_chi_tan.shape[1]

half_bin = int(mg_bin_num/2)
inverse = range(half_bin-1, -1, -1)

chi_data = [chi_tan,chi_cross]
chi_crit_data = [chi_crit_tan, chi_crit_cross]

# the chi square
chisq = numpy.zeros((gh_num, 2))
chisq_crit = numpy.zeros((gh_crit_num, 2))

chisq_fit = []
chisq_crit_fit = []

itemsize = MPI.DOUBLE.Get_size()
data_row, data_col = cpus, 8
element_num = data_row*data_col
if rank == 0:
    nbytes = element_num*itemsize
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
final_signal = numpy.ndarray(buffer=buf1, dtype='d', shape=(data_row,data_col)) # array filled with zero
comm.Barrier()

for i in range(2):
    # shear
    chi_left = chi_data[i][:,:half_bin][:,inverse]
    chi_right = chi_data[i][:,half_bin:]
    dn = (chi_right - chi_left)**2
    dm = chi_right + chi_left
    chi = numpy.sum(dn/dm,axis=1)
    idx = chi == numpy.nan
    chi[idx] = -10
    chisq[:,i] = chi

    idx = chi >= 0
    chi_min = chi[idx].min()
    idx = chi <= chi_min + 40

    coeff = tool_box.fit_1d(gh[idx], chisq[:,i][idx], 2, "scipy")
    chisq_fit.append((chisq[:,i][idx], gh[idx],coeff))
    g_h = -coeff[1] / 2. / coeff[2]
    g_sig = 0.70710678118 / numpy.sqrt(coeff[2])
    final_signal[rank, i*2] = g_h
    final_signal[rank, i * 2 + 1] = g_sig

    # shear*critical_surface_density
    chi_left = chi_crit_data[i][:,:half_bin][:,inverse]
    chi_right = chi_crit_data[i][:,half_bin:]
    dn = (chi_right - chi_left)**2
    dm = chi_right + chi_left
    chi = numpy.sum(dn/dm,axis=1)
    idx = chi == numpy.nan
    chi[idx] = -10
    chisq_crit[:,i] = chi

    idx = chi >= 0
    chi_min = chi[idx].min()
    idx = chi <= chi_min + 40
    coeff = tool_box.fit_1d(gh_crit[idx], chisq_crit[:,i][idx], 2, "scipy")
    chisq_crit_fit.append((chisq_crit[:,i][idx], gh_crit[idx], coeff))
    g_h = -coeff[1] / 2. / coeff[2]*388.2833518
    g_sig = 0.70710678118 / numpy.sqrt(coeff[2])*388.2833518
    final_signal[rank, i*2 + 4] = g_h
    final_signal[rank, i * 2 + 5] = g_sig


img = Image_Plot()
img.create_subfig(2,2)

img.axs[0][0].scatter(gh, chisq[:,0],c="C1",s=3)
img.axs[0][0].scatter(chisq_fit[0][1], chisq_fit[0][0],c="C2",s=5)

coeff = chisq_fit[0][2]
fit_range = chisq_fit[0][1]
fx = coeff[0]+coeff[1]*fit_range+coeff[2]*fit_range**2
img.axs[0][0].plot(chisq_fit[0][1], fx, c="C2")


img.axs[0][1].scatter(gh, chisq[:,1],c="C1",s=3)
img.axs[0][1].scatter(chisq_fit[1][1], chisq_fit[1][0],c="C2",s=5)

coeff = chisq_fit[1][2]
fit_range = chisq_fit[1][1]
fx = coeff[0]+coeff[1]*fit_range+coeff[2]*fit_range**2
img.axs[0][1].plot(chisq_fit[1][1], fx, c="C2")


img.axs[1][0].scatter(gh_crit, chisq_crit[:,0],c="C1",s=3)
img.axs[1][0].scatter(chisq_crit_fit[0][1], chisq_crit_fit[0][0],c="C2",s=5)

coeff = chisq_crit_fit[0][2]
fit_range = chisq_crit_fit[0][1]
fx = coeff[0]+coeff[1]*fit_range+coeff[2]*fit_range**2
img.axs[1][0].plot(chisq_crit_fit[0][1], fx, c="C2")

img.axs[1][1].scatter(gh_crit, chisq_crit[:,1],c="C1",s=3)
img.axs[1][1].scatter(chisq_crit_fit[1][1], chisq_crit_fit[1][0],c="C2",s=5)

coeff = chisq_crit_fit[1][2]
fit_range = chisq_crit_fit[1][1]
fx = coeff[0]+coeff[1]*fit_range+coeff[2]*fit_range**2
img.axs[1][1].plot(chisq_crit_fit[1][1], fx, c="C2")

if area_num == 1:
    img.save_img(data_path + "pic/w_%d/%d.png"%(area_ids[0],rank))
else:
    img.save_img(data_path + "pic/total/%d.png" %rank)
img.close_img()
comm.Barrier()

if rank == 0:

    background_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/cata_result_ext_grid.hdf5"
    h5f = h5py.File(background_path,"r")
    radius_bin = h5f["/radius_bin"].value
    h5f.close()
    x = (radius_bin[:cpus] + radius_bin[1:])/2
    img = Image_Plot()
    img.create_subfig(1, 2)
    img.axs[0][0].errorbar(x, final_signal[:, 0],final_signal[:, 1], capsize=3, c="C1", label="T")
    img.axs[0][0].errorbar(x, final_signal[:, 2],final_signal[:, 3], capsize=3, c="C2", label="X")

    img.axs[0][1].errorbar(x, final_signal[:, 4],final_signal[:, 5], capsize=3, c="C1",  label="T")
    img.axs[0][1].errorbar(x, final_signal[:, 6],final_signal[:, 7], capsize=3, c="C2",  label="X")
    img.axs[0][0].legend()
    img.axs[0][1].legend()
    img.axs[0][0].set_xscale("log")
    img.axs[0][1].set_xscale("log")
    # img.axs[0][0].set_yscale("symlog")
    img.axs[0][1].set_yscale("symlog")
    img.tick_label(0, 0, 0,"$\gamma$")
    img.tick_label(0, 0, 1, "$R \quad Mpc\cdot h^{-1}$")
    img.tick_label(0, 1, 0, "$\Delta\Sigma$")
    img.tick_label(0, 1, 1, "$R \quad Mpc\cdot h^{-1}$")
    if area_num == 1:
        img.save_img(data_path + "pic/w_%d/signal.png" %area_ids[0])
    else:
        img.save_img(data_path + "pic/total/total_signal.png")