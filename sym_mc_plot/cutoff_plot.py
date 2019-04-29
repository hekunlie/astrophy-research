import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import plot_tool
from mpi4py import MPI
import tool_box
import numpy
from sys import argv
import h5py



# This program plots the figure of each cutoff

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

filter_name = ["sex2_4", "sex2_2", "sex2_1.5", "sex3_4", "sex3_2", "sex3_1.5","sex4_4", "sex4_2", "sex4_1.5"]
select_name = ["mag_auto", "snr_auto", 'sex_snr', "flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4", "flux2_ex5"]
source = argv[1]

ini_path = "%s/work/envs/envs.dat"%my_home
total_path = tool_box.config(ini_path, ['get'], [['selection_bias', "%s_path"%source, '1']])[0]

shear_path = total_path + "parameters/shear.npz"
g1 = numpy.load(shear_path)["arr_0"]
g2 = numpy.load(shear_path)["arr_1"]
# range of plot
g1_extr = [g1.min(), g1.max()]
g2_extr = [g2.min(), g2.max()]
g1_plt_s, g1_plt_e = g1.min() - (g1.max() - g1.min()) * 0.2, g1.max() + (g1.max() - g1.min()) * 0.2
g2_plt_s, g2_plt_e = g2.min() - (g2.max() - g2.min()) * 0.2, g2.max() + (g2.max() - g2.min()) * 0.2
plt_range = [g1_plt_s, g1_plt_e, g2_plt_s, g2_plt_e]

file_paths = []
for fnm in filter_name:
    for snm in select_name:
        file_path = total_path + "result/cuts/sym/%s/%s/"%(fnm, snm)
        file_paths.append(file_path)

my_files = tool_box.allot(file_paths, cpus)[rank]

for file_path in my_files:

    pdfs = os.listdir(file_path)
    for pdf_file in pdfs:
        if ".pdf" in pdf_file:
            os.remove(file_path + pdf_file)
    h5f = h5py.File(file_path + "total.hdf5","r")

    mc1 = h5f["/mc1"].value
    mc2 = h5f["/mc2"].value
    results = h5f["/shear"].value
    num = h5f["/num"].value
    cut_scale = h5f["/cut_scale"].value

    # the row labels are corresponding to the shears
    shear_num, cut_num = mc1.shape

    h5f.close()

    for i in range(cut_num):
        mc_title = ['0', '0', '0', '0']
        e1mc = (mc1[0, i]+1, mc1[1, i], mc1[2, i], mc1[3, i])
        e2mc = (mc2[0, i]+1, mc2[1, i], mc2[2, i], mc2[3, i])
        # m1/2 +- 2sig
        m_r = [[e1mc[0] - 1 - 2 * e1mc[1], e1mc[0] - 1 + 2 * e1mc[1]],
               [e2mc[0] - 1 - 2 * e2mc[1], e2mc[0] - 1 + 2 * e2mc[1]]]
        c_r = [[e1mc[2] - 2 * e1mc[3], e1mc[2] + 2 * e1mc[3]],
               [e2mc[2] - 2 * e2mc[3], e2mc[2] + 2 * e2mc[3]]]

        for ii in range(2):
            if tool_box.check_in(m_r[ii]):
                mc_title[ii] = ''
            else:
                mc_title[ii] = "_m" + str(ii+1)
            if tool_box.check_in(c_r[ii]):
                mc_title[ii + 2] = ''
            else:
                mc_title[ii + 2] = "_c" + str(ii+1)
        pic_mc = "".join(mc_title)

        pic_path = file_path + "/" + str(round(cut_scale[0,i], 4)) + pic_mc+".pdf"

        tool_box.mcplot(g1, results[:,i], results[:,i+cut_num], num[:,i],
                        g2, results[:,i+2*cut_num], results[:,i+3*cut_num], num[:,i],
                        e1mc, e2mc, str(round(cut_scale[0,i],4)), 'max', plt_range,pic_path)




