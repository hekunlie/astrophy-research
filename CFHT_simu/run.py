from subprocess import Popen
import os
from sys import path
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
path.append('%s/work/fourier_quad/'%my_home)
import tool_box
import shutil


envs_path = "%s/work/envs/envs.dat" % my_home
get_contents = [['cfht', "cfht_path_cut", '1']]
gets = ["get" for i in range(len(get_contents))]
cut_path = tool_box.config(envs_path, gets, get_contents)[0]
print(cut_path)
cut = "flux_alt"
stars = [10, 12, 14, 16, 18]
areas = [1, 3, 5, 7, 9]
for i in stars:
    for j in areas:

        old_nm = cut_path + cut
        new_nm = cut_path + cut + "_s%d_a%d"%(i, j)
        if os.path.exists(old_nm):
            # remove the old directory
            shutil.rmtree(old_nm)
            # build new one
            os.mkdir(old_nm)
        else:
            os.mkdir(old_nm)

        cmd = "mpirun -np 50 python sym_mc_plot_cfht.py %s %d %d result_int"%(cut, i, j)
        a = Popen(cmd,shell=True)
        a.wait()

        os.rename(old_nm, new_nm)

