import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
from plot_tool import Image_Plot
from mpi4py import MPI
import tool_box
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

def plot_gf_celestial(ra, dec, gf1, gf2, gf, gf1_scale, gf2_scale, gf_scale,dot_size=1,pic_path=None):

    img = Image_Plot(xpad=0.2,ypad=0.1)
    img.subplots(1, 3)
    color_cm = 'bwr'

    norm = plt.Normalize(vmin=numpy.min(gf1_scale[0]), vmax=numpy.max(gf1_scale[1]))
    cmap = plt.get_cmap(color_cm)
    cl = cmap(norm(gf1))
    fig = img.axs[0][0].scatter(ra,dec, color=cl, s=dot_size)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    img.figure.colorbar(sm, ax=img.axs[0][0])

    norm = plt.Normalize(vmin=numpy.min(gf2_scale[0]), vmax=numpy.max(gf2_scale[1]))
    cmap = plt.get_cmap(color_cm)
    cl = cmap(norm(gf2))
    fig = img.axs[0][1].scatter(ra,dec, color=cl, s=dot_size)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    img.figure.colorbar(sm, ax=img.axs[0][1])

    norm = plt.Normalize(vmin=numpy.min(gf_scale[0]), vmax=numpy.max(gf_scale[1]))
    cmap = plt.get_cmap(color_cm)
    cl = cmap(norm(gf))
    fig = img.axs[0][2].scatter(ra,dec, color=cl, s=dot_size)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    img.figure.colorbar(sm, ax=img.axs[0][2])

    for i in range(3):
        img.axis_sci_ticklabel(0,i,0)
        img.axis_sci_ticklabel(0,i,1)
        img.set_label(0,i,0,"y")
        img.set_label(0,i,1,"x")
    if pic_path:
        img.save_img(pic_path)
    # img.show_img()
    img.close_img()

def plot_gf_pix(x, y, ichip, gf1, gf2, gf, gf1_scale, gf2_scale,gf_scale,dot_size=1,pic_path=None):

    chip_row, chip_col = numpy.divmod(ichip, 9)

    x_ccd = x + chip_col * 2112
    y_ccd = y + chip_row * 4644

    img = Image_Plot(xpad=0.2,ypad=0.1)
    img.subplots(1, 3)
    color_cm = 'bwr'

    norm = plt.Normalize(vmin=numpy.min(gf1_scale[0]), vmax=numpy.max(gf1_scale[1]))
    cmap = plt.get_cmap(color_cm)
    cl = cmap(norm(gf1))
    fig = img.axs[0][0].scatter(x_ccd, y_ccd, color=cl, s=dot_size)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    img.figure.colorbar(sm, ax=img.axs[0][0])

    norm = plt.Normalize(vmin=numpy.min(gf2_scale[0]), vmax=numpy.max(gf2_scale[1]))
    cmap = plt.get_cmap(color_cm)
    cl = cmap(norm(gf2))
    fig = img.axs[0][1].scatter(x_ccd, y_ccd, color=cl, s=dot_size)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    img.figure.colorbar(sm, ax=img.axs[0][1])

    norm = plt.Normalize(vmin=numpy.min(gf_scale[0]), vmax=numpy.max(gf_scale[1]))
    cmap = plt.get_cmap(color_cm)
    cl = cmap(norm(gf))
    fig = img.axs[0][2].scatter(x_ccd, y_ccd, color=cl, s=dot_size)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    img.figure.colorbar(sm, ax=img.axs[0][2])

    for i in range(3):
        img.axis_sci_ticklabel(0,i,0)
        img.axis_sci_ticklabel(0,i,1)
        img.set_label(0,i,0,"y")
        img.set_label(0,i,1,"x")
    if pic_path:
        img.save_img(pic_path)
    # img.show_img()
    img.close_img()

total_path = argv[1]

src_path = total_path + "/cat_hdf5"
pic_path = total_path + "/cat_inform/field_distortion"

source_list_nm = "exposure_avail.dat"

with open(total_path + "/cat_inform/"+source_list_nm, "r") as f:
    content = f.readlines()

all_exposures_path = []
all_exposures = []
for cc in content:
    all_exposures.append(cc.split("/")[-1].split("\n")[0])
    all_exposures_path.append(cc.split("\n")[0])

my_exposures_path = tool_box.alloc(all_exposures_path, cpus)[rank]
my_exposures = tool_box.alloc(all_exposures, cpus)[rank]


for tag, fnm in enumerate(my_exposures):
    # print(my_exposures_path[tag], fnm)

    h5f = h5py.File(my_exposures_path[tag], "r")
    data = h5f["/data"][()]
    h5f.close()

    col_shift = 16
    ichip = data[:, col_shift] - 1
    xc = data[:,col_shift+2]
    yc = data[:,col_shift+3]
    ra = data[:,col_shift+13]
    dec = data[:,col_shift+14]

    gf1 = data[:, col_shift + 15]
    gf2 = data[:, col_shift + 16]
    gf = numpy.sqrt(gf1 ** 2 + gf2 ** 2)

    gf1_scale = [gf1.min(), gf1.max()]
    gf2_scale = [gf2.min(), gf2.max()]
    gf_scale = [gf.min(), gf.max()]

    pic_nm = pic_path + "/%s.png"%fnm.split(".")[0]
    # plot_gf_celestial(ra, dec, gf1, gf2, gf, gf1_scale, gf2_scale, gf_scale, dot_size=1, pic_path=pic_nm)
    plot_gf_pix(xc,yc,ichip, gf1, gf2, gf, gf1_scale, gf2_scale, gf_scale, dot_size=1, pic_path=pic_nm)


