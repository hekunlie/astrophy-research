import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/' % my_home)
from subprocess import Popen
import numpy
import tool_box
import time

sources = ["dimmerm3"]
cpus = 70
num = 140
jobs = numpy.zeros((num, 1))

t1 = time.time()
for source in sources:

    envs_path = "%s/work/envs/envs.dat" % my_home
    get_contents = [['selection_bias', "%s_path_para" % source, '1']]
    para_path = tool_box.config(envs_path, ['get'], get_contents)[0]

    # the parameters
    para_contents = [["para", "noise_sig", 1], ["para","shear_num",1]]
    para_items = tool_box.config(para_path + "para.ini", ['get',"get"], para_contents)

    noise_sig = float(para_items[0])
    shear_num = int(para_items[1])

    while True:
        if jobs.sum() == num:
            break
        for i in range(num):
            if os.path.exists("%s/work/test/job/%s/finish_%d.dat"%(my_home, source, i)):
                jobs[i, 0] = 1

    if "pts" in source:
        max_radius = 8
    else:
        max_radius = 5.5

    filter_names = ["sex2_2", "sex3_2", "sex4_2",
                    "sex2_1.5", "sex3_1.5", "sex4_1.5"]
    gauss_filters = ["gauss_2.0_5x5", "gauss_3.0_5x5", "gauss_4.0_7x7",
                     "gauss_2.0_5x5", "gauss_3.0_5x5", "gauss_4.0_7x7"]

    for ii, filter_name in enumerate(filter_names):
        # change the filter
        with open("./default.sex_ori", "r") as f:
            contents = f.readlines()
        sig_level = float(filter_name.split("_")[1])*noise_sig
        contents[16] = 'DETECT_THRESH    %.2f          # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n'%sig_level
        contents[18] = 'ANALYSIS_THRESH  %.2f            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n' % sig_level
        contents[23] = "FILTER_NAME      %s.conv   # name of the file containing the filter\n"%gauss_filters[ii]
        with open("./default.sex", "w") as f:
            f.writelines(contents)

        # run the estimation code
        cmd = "mpirun -np %d python snr_est.py snr %s %s %.1f"%(cpus, source, filter_name, max_radius)
        a = Popen(cmd, shell=True)
        a.wait()

        # add the data in .cat files to the final catalogs
        cmd = "mpirun -np %d python snr_est.py add %s %s %.1f"%(shear_num,source, filter_name, max_radius)
        a = Popen(cmd, shell=True)
        a.wait()

        # check
        cmd = "mpirun -np %d python snr_est.py check %s %s %.1f"%(shear_num,source, filter_name, max_radius)
        a = Popen(cmd, shell=True)
        a.wait()
t2 = time.time()
print("SNR EST: ",sources, cpus, t2-t1)

