from subprocess import Popen
import time
import os
from sys import path,argv
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/fourier_quad/'%my_home)
import shutil
import shelve

with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        data_path = path.split("=")[1]
    elif "cfht_cut_path" in path:
        cut_path = path.split("=")[1]

source_list = data_path + "source.list"

curr_path = os.getcwd()+"/" #+ "/src_new_%s/"%argv[1]

test_file = curr_path + "shear_field_distortion_test.f"
with open(test_file, 'r') as f:
    conts = f.readlines()

# istar 59  iflux_alt 60
# stars = "                if (cat(istar).lt.%f) cycle"
# thres = "                if (cat(%s).lt.%f) cycle"

f = shelve.open(cut_path+"cut_dict")
cuts_list = f["dict"]
f.close()

fluxs = ["iflux_alt", "iflux2"]
cuts = [cuts_list["flux_alt"], cuts_list["flux2"]]

nstars = [int(os.path.basename(os.getcwd()).split("_")[-1])]
# nstars = [int(argv[1])]

for i in range(len(nstars)):
    stars = "                if (cat(istar).lt.%.1f) cycle\n"%nstars[i]
    star_path = cut_path + "fresh/%d/"%nstars[i]
    if not os.path.exists(star_path):
        os.mkdir(star_path)

    for j in range(len(fluxs)):
        flux_path = star_path + fluxs[j] + "/"
        if not os.path.exists(flux_path):
            os.mkdir(flux_path)

        for k in range(len(cuts[j])):
            cutoff_path = flux_path + "/%d/"%k
            if not os.path.exists(cutoff_path):
                os.mkdir(cutoff_path)

            thres = "                if (cat(%s).lt.%.4f) cycle\n"%(fluxs[j], cuts[j][k])
            conts[59], conts[60] = stars, thres
            with open(test_file, "w") as f:
                f.writelines(conts)

            main_path = curr_path + "main"
            if os.path.exists(main_path):
                os.remove(main_path)

            # compiling
            cmd = "mpif77 %s*.f -o %s -mcmodel=medium -lcfitsio"%(curr_path, main_path)
            a = Popen(cmd, shell=True)
            a.wait()
            print("Compiling")
            # run code
            cmd = "mpiexec -n 1 %s %s"%(main_path,source_list)
            a = Popen(cmd, shell=True)
            a.wait()

            all_files = os.listdir(curr_path)

            for target in all_files:
                if "full_shear_" in target:
                    old_name = curr_path + "/%s"%target
                    new_name = cutoff_path + target
                    shutil.move(old_name, new_name)

