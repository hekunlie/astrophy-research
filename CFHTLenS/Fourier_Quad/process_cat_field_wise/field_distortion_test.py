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
import matplotlib.pyplot as plt



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


def plot_gf(ra, dec, gf1, gf2, pic_path, gf1_bin_num = 15, gf2_bin_num = 18,gf_bin_num = 18):
    gf = numpy.sqrt(gf1**2 + gf2**2)

    label_1 = numpy.zeros_like(gf1)
    label_2 = numpy.zeros_like(gf2)
    label_3 = numpy.zeros_like(gf)

    gf1_bin = numpy.linspace(gf1.min(), gf1.max(), gf1_bin_num+1)
    gf2_bin = numpy.linspace(gf2.min(), gf2.max(), gf2_bin_num+1)
    gf_bin = numpy.linspace(gf.min(), gf.max(), gf_bin_num+1)

    for i in range(gf1_bin_num):
        idx1 = gf1 >= gf1_bin[i]
        idx2 = gf1 < gf1_bin[i+1]
        idx = idx1 & idx2
        label_1[idx] = (gf1_bin[i]+ gf1_bin[i+1])/2

    for i in range(gf2_bin_num):
        idx1 = gf2 >= gf2_bin[i]
        idx2 = gf2 < gf2_bin[i+1]
        idx = idx1 & idx2
        label_2[idx] = (gf2_bin[i]+ gf2_bin[i+1])/2

    for i in range(gf_bin_num):
        idx1 = gf >= gf_bin[i]
        idx2 = gf < gf_bin[i+1]
        idx = idx1 & idx2
        label_3[idx] = (gf_bin[i]+ gf_bin[i+1])/2


    img = Image_Plot()
    img.subplots(1,3)

    norm = plt.Normalize(vmin=numpy.min(label_1), vmax=numpy.max(label_1))
    cmap = plt.get_cmap('YlOrRd',gf1_bin_num)
    cl = cmap(norm(label_1))
    fig = img.axs[0][0].scatter(ra, dec, color=cl, s=3)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    img.figure.colorbar(sm, ax=img.axs[0][0])


    norm = plt.Normalize(vmin=numpy.min(label_2), vmax=numpy.max(label_2))
    cmap = plt.get_cmap('YlOrRd',gf2_bin_num)
    cl = cmap(norm(label_2))
    fig = img.axs[0][1].scatter(ra, dec, color=cl, s=3)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    img.figure.colorbar(sm, ax=img.axs[0][1])


    norm = plt.Normalize(vmin=numpy.min(label_3), vmax=numpy.max(label_3))
    cmap = plt.get_cmap('YlOrRd',gf_bin_num)
    cl = cmap(norm(label_3))
    fig = img.axs[0][2].scatter(ra, dec, color=cl, s=3)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    img.figure.colorbar(sm, ax=img.axs[0][2])

    img.save_img(pic_path)
    # img.show_img()
    img.close_img()


total_path = argv[1]

src_path = total_path
pic_path = total_path + "/cat_inform/field_distortion/shear_stack"

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

# for fnm in sub_fields:
#     if os.path.exists(src_path + "/%s/result/%s_raw.hdf5" % (fnm, fnm)):
#         h5f = h5py.File(src_path + "/%s/result/%s_raw.hdf5" % (fnm, fnm), "r")
#         data = h5f["/data"][()]
#         h5f.close()
#         col_shift = 0
#         # ichip = data[:, col_shift]
#         # xc = data[:, col_shift + 2]
#         # yc = data[:, col_shift + 3]
#
#         nstar = data[:,col_shift+5]
#
#         flux_alt = data[:, col_shift + 12]
#
#         ra = data[:, col_shift + 13]
#         dec = data[:, col_shift + 14]
#
#         ra_cent = (ra.max() + ra.min())/2
#         dec_cent = (dec.max() + dec.min())/2
#
#         idx1 = nstar >= 12
#         idx2 = flux_alt >= 3
#         idx3 = numpy.abs(data[:, col_shift + 15]) <= 0.005
#         idx4 = numpy.abs(data[:, col_shift + 16]) <= 0.005
#
#         idxr = ra >= ra_cent
#         idxl = ra < ra_cent
#         idxd = dec < dec_cent
#         idxu = dec >= dec_cent
#
#         region = [idxu, idxd, idxl, idxr]
#         labels = ["upper", "lower", "left", "right"]
#
#         img = Image_Plot(xpad=0.25, ypad=0.1)
#         img.subplots(2, 2)
#
#         for i in range(4):
#             m, n = divmod(i,2)
#
#             result[:,:] = - 99
#
#             gf1 = data[:, col_shift + 15][idx1&idx2&idx3&idx4&region[i]]
#             # gf2 = data[:, col_shift + 16][idx1&idx2&region[i]]
#
#             mg1 = data[:,col_shift + 17][idx1&idx2&idx3&idx4&region[i]]
#             # mg2 = data[:,col_shift + 18][idx1&idx2&region[i]]
#             mnu1 = data[:,col_shift + 19][idx1&idx2&idx3&idx4&region[i]] - data[:,col_shift + 20][idx1&idx2&idx3&idx4&region[i]]
#             # mu = -data[:,col_shift + 20][idx1&idx2&region[i]]
#             # mv = -data[:,col_shift + 21][idx1&idx2&region[i]]
#
#             for j in range(bin_num):
#                 idxf1 = gf1 >= gf_bin[j]
#                 idxf2 = gf1 < gf_bin[j+1]
#
#                 idx = idxf1 & idxf2
#                 if idx.sum() > 1000:
#                     gh, gh_sig = fq.find_shear(mg1[idx], mnu1[idx],8)[:2]
#                     result[0,j] = gh
#                     result[1,j] = gh_sig
#             idxp = result[0] > - 99
#             img.axs[m][n].errorbar(gf_bin_c[idxp], result[0][idxp], result[1][idxp],
#                                    capsize=3, c="C1", label="%s g1"%labels[i], marker="o")
#             xs = img.axs[m][n].set_xlim()
#             ys = img.axs[m][n].set_ylim()
#             max_x = max([abs(xs[0]), abs(xs[1])])
#             max_y = max([abs(ys[0]), abs(ys[1])])
#             max_range = max(max_x, max_y)
#
#             img.axs[m][n].plot([-max_range, max_range], [-max_range, max_range], linestyle="--", c="grey")
#             img.axs[m][n].legend(loc="upper left")
#             img.axs[m][n].set_xlim(-max_x, max_x)
#             img.axs[m][n].set_ylim(-max_y, max_y)
#
#             img.set_label(m,n, 0, "$g_1$")
#             img.set_label(m,n, 1, "True $g_1$")
#         img.save_img(pic_path + "/%s.png"%fnm)
#         img.close_img()

sub_sample_num = divmod(len(sub_fields), 4)[0]
sub_fields_sep = tool_box.alloc(sub_fields, sub_sample_num)
col_shift = 0
for sub_tag in range(sub_sample_num):
    tag = 0
    for fnm in sub_fields_sep[sub_tag]:
        if os.path.exists(src_path + "/%s/result/%s_raw.hdf5" % (fnm, fnm)):
            h5f = h5py.File(src_path + "/%s/result/%s_raw.hdf5" % (fnm, fnm), "r")
            temp = h5f["/data"][()]
            h5f.close()

            if tag == 0:
                data = temp
            else:
                data = numpy.row_stack((data, temp))
            tag += 1

    if tag > 0:
        nstar = data[:,col_shift+5]

        flux_alt = data[:, col_shift + 12]

        ra = data[:, col_shift + 13]
        dec = data[:, col_shift + 14]

        ra_cent = (ra.max() + ra.min())/2
        dec_cent = (dec.max() + dec.min())/2

        idx1 = nstar >= 12
        idx2 = flux_alt >= 3
        idx3 = numpy.abs(data[:, col_shift + 15]) <= 0.005
        idx4 = numpy.abs(data[:, col_shift + 16]) <= 0.005
        idxs = idx1*idx2&idx3&idx4




        idxu = dec >= dec_cent
        idxd = dec < dec_cent
        idxl = ra < ra_cent
        idxr = ra >= ra_cent

        region = [idxs&idxu&idxl, idxs&idxu&idxr, idxs&idxd&idxl, idxs&idxd&idxr]
        labels = ["upper left", "upper right", "lower left", "lower right"]
        print(data.shape[0],idxs.sum(), [region[i].sum() for i in range(4)],sum([region[i].sum() for i in range(4)]))

        img = Image_Plot(xpad=0.25, ypad=0.1)
        img.subplots(2, 2)

        for i in range(4):
            m, n = divmod(i,2)

            result[:,:] = - 99

            gf1 = data[:, col_shift + 15][region[i]]
            # gf2 = data[:, col_shift + 16][idx1&idx2&region[i]]

            mg1 = data[:,col_shift + 17][region[i]]
            # mg2 = data[:,col_shift + 18][idx1&idx2&region[i]]
            mnu1 = data[:,col_shift + 19][region[i]] - data[:,col_shift + 20][region[i]]
            # mu = -data[:,col_shift + 20][idx1&idx2&region[i]]
            # mv = -data[:,col_shift + 21][idx1&idx2&region[i]]

            for j in range(bin_num):
                idxf1 = gf1 >= gf_bin[j]
                idxf2 = gf1 < gf_bin[j+1]

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
        img.save_img(pic_path + "/%s.png"%("_".join(sub_fields_sep[sub_tag])))
        img.close_img()