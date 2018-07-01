from subprocess import Popen
import time
import os
from sys import path
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/fourier_quad/'%my_home)
import shutil
import shelve
from mpi4py import MPI
from sys import argv


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        data_path = path.split("=")[1]
    elif "cfht_cut_path" in path:
        cut_path = path.split("=")[1]

source_list = data_path + "source.list"

test_file = "shear_field_distortion_test.f"
with open(test_file, 'r') as f:
    conts = f.readlines()

f = shelve.open(cut_path+"cut_dict")
cuts_list = f["dict"]
f.close()
cuts_num = len(cuts_list)
fluxs = ["iflux_alt", "iflux2"]
cuts = [cuts_list["flux_alt"], cuts_list["flux2"]]
nstars = [12, 16, 20, 30]
# prepare the directory
if argv[1] == "mkdir":
    fresh_path = cut_path + "fresh/"
    if os.path.exists(fresh_path):
        shutil.rmtree(fresh_path)
    os.mkdir(fresh_path)
    for i in range(len(fluxs)):
        flux_path = fresh_path + "%s/"%fluxs[i]
        os.mkdir(flux_path)
        for ns in nstars:
            star_path = flux_path + "%d/"%ns
            for j in range(cuts_num):
                c_path = star_path + ""

f = shelve.open("code_file")
code_names = f['file']
f.close()
shear_fd_name = code_names[-1]
comp_name = ""
for name in code_names:
    if name != shear_fd_name:
        comp_name += "%s "%name


# istar 59  iflux_alt 60
# stars = "                if (cat(istar).lt.%f) cycle"
# thres = "                if (cat(%s).lt.%f) cycle"

if cpus != len(fluxs)*len(cuts):
    print("Number of task doesn't equal number threads!")
    exit()

crit_id, nstars_id = divmod(rank, len(nstars))

stars = "                if (cat(istar).lt.%.1f) cycle\n" % nstars[nstars_id]
star_path = cut_path + "fresh/%d/" % nstars[nstars_id]
if not os.path.exists(star_path):
    os.mkdir(star_path)

for i in range(cuts_num):
    flux_path = star_path + fluxs[crit_id] + "/"
    if not os.path.exists(flux_path):
        os.mkdir(flux_path)

    for k in range(len(cuts[j])):
        cutoff_path = flux_path + "/%d/"%k
        if not os.path.exists(cutoff_path):
            os.mkdir(cutoff_path)

        thres = "                if (cat(%s).lt.%.2f) cycle\n"%(fluxs[j], (cuts[j][k])**2)
        conts[59], conts[60] = stars, thres
        with open(test_file, "w") as f:
            f.writelines(conts)

        if os.path.exists("main"):
            os.remove("main")

        # compiling
        cmd = "mpif77 *.f -o main -mcmodel=medium -lcfitsio"
        a = Popen(cmd, shell=True)
        a.wait()
        print("Compiling")
        # run code
        cmd = "mpiexec -n 1 ./main %s"%source_list
        a = Popen(cmd, shell=True)
        a.wait()

        curr_path = os.getcwd()
        all_files = os.listdir(curr_path)
        for target in all_files:
            if "full_shear_" in target:
                old_name = curr_path + "/%s"%target
                new_name = cutoff_path + target
                shutil.move(old_name, new_name)


