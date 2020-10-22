import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
from mpi4py import MPI
import tool_box
from Fourier_Quad import Fourier_Quad
from plot_tool import Image_Plot


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


total_path = argv[1]

src_path = total_path
pic_path = total_path + "/cat_inform/field_distortion/shear"

source_list_nm = "/cat_inform/nname_field_raw_avail.dat"

fields = []
with open(src_path + source_list_nm, "r") as f:
    conts = f.readlines()

for nm in conts:
    field_nm = nm.split("\n")[0]
    if field_nm not in fields:
        fields.append(field_nm)

if rank == 0:
    print("Totally ",len(fields)," fields")

sub_fields = tool_box.alloc(fields, cpus)[rank]

fq = Fourier_Quad(12,123)

bin_num = 20
gf_bin = numpy.linspace(-0.005, 0.005, bin_num+1)
gf_bin_c = (gf_bin[1:] + gf_bin[:-1])/2

result = numpy.zeros((2, bin_num))

for fnm in sub_fields:
    if os.path.exists(src_path + "/%s/result/%s_raw.hdf5" % (fnm, fnm)):
        h5f = h5py.File(src_path + "/%s/result/%s_raw.hdf5" % (fnm, fnm), "r")
        data = h5f["/data"][()]
        h5f.close()
        col_shift = 0
        # ichip = data[:, col_shift]
        # xc = data[:, col_shift + 2]
        # yc = data[:, col_shift + 3]

        nstar = data[:,col_shift+5]

        flux_alt = data[:, col_shift + 12]

        ra = data[:, col_shift + 13]
        dec = data[:, col_shift + 14]

        ra_cent = (ra.max() + ra.min())/2
        dec_cent = (dec.max() + dec.min())/2

        idx1 = nstar >= 12
        idx2 = flux_alt >= 3


        idxr = ra >= ra_cent
        idxl = ra < ra_cent
        idxd = dec < dec_cent
        idxu = dec >= dec_cent

        region = [idxu, idxd, idxl, idxr]
        labels = ["upper", "lower", "left", "right"]

        img = Image_Plot(xpad=0.25, ypad=0.1)
        img.subplots(1, 4)

        for i in range(4):
            m, n = divmod(i,2)

            result[:,:] = - 99

            gf1 = data[:, col_shift + 15][idx1&idx2&region[i]]
            # gf2 = data[:, col_shift + 16][idx1&idx2&region[i]]

            mg1 = data[:,col_shift + 17][idx1&idx2&region[i]]
            # mg2 = data[:,col_shift + 18][idx1&idx2&region[i]]
            mnu1 = data[:,col_shift + 19][idx1&idx2&region[i]] - data[:,col_shift + 20][idx1&idx2&region[i]]
            # mu = -data[:,col_shift + 20][idx1&idx2&region[i]]
            # mv = -data[:,col_shift + 21][idx1&idx2&region[i]]

            for j in range(bin_num):
                idxf1 = gf1 >= gf_bin[i]
                idxf2 = gf1 < gf_bin[i+1]

                idx = idxf1 & idxf2
                if idx.sum() > 1000:
                    gh, gh_sig = fq.find_shear(mg1[idx], mnu1[idx],8)[:2]
                    result[0,j] = gh
                    result[1,j] = gh_sig
            idxp = result[0] > - 99
            img.axs[m][n].errorbar(gf_bin_c[idxp], result[0][idxp], result[1][idxp],
                                   capsize=3, c="C1", label="%s g1"%labels[i], marker="o")
            xs = img.axs[m][n].set_xlim()
            ys = img.axs[m][n].set_ylim()
            max_x = max([abs(xs[0]), abs(xs[1])])
            max_y = max([abs(ys[0]), abs(ys[1])])
            max_range = max(max_x, max_y)

            img.axs[m][n].plot([-max_range, max_range], [-max_range, max_range], linestyle="--", c="grey")
            img.axs[m][n].legend(loc="upper left")
            img.axs[m][n].set_xlim(-max_x, max_x)
            img.axs[m][n].set_ylim(-max_y, max_y)

            img.set_label(m,n, 0, "$g_1$")
            img.set_label(m,n, 1, "True $g_1$")
        img.save_img(pic_path + "/%s.png"%fnm)
        img.close_img()