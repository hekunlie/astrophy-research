import matplotlib
matplotlib.use('Agg')
import os
from sys import path,argv
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/fourier_quad/'%my_home)
import tool_box
from Fourier_Quad import Fourier_Quad
import numpy
from mpi4py import MPI
import time
from astropy.io import fits
import copy
import warnings

warnings.filterwarnings("error")

# to find the binary on the each source chip. it will save the binary label for each source.
# '1' means binary
# if the command input is 'find', it will find the binaries
# if 'stack' is input, it will stack the existing binary label files.

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()

with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        data_path = path.split("=")[1]
    elif "cfht_res_path" in path:
        result_path = path.split("=")[1]
    elif "cfht_pic_path" in path:
        pic_path = path.split("=")[1]
    elif "cfht_field_path" in path:
        field_path = path.split("=")[1]

size = 48
fq = Fourier_Quad(size, 123)

nname_path = data_path + "nname.dat"
field_dict, fields = tool_box.field_dict(nname_path)
r_fields = tool_box.allot(fields,cpus)[rank]

# for the stacking process
count = 0

# the location of each galaxy is labeled by the field_label and exposure_label
# counting from the left, the first, third and fifth figure denotes "w_m(p)_(m)p_"
# the second and the fourth denotes "m" or "p" (1=m,0=p)
# the last two figure denote the chip NO.
for field in r_fields:
    expos = list(field_dict[field].keys())
    field_label = tool_box.cfht_label(field)
    for expo in expos:
        expo_label = int(expo.split("p")[0])
        chips = field_dict[field][expo]
        for chip in chips:
            chip_label = int(chip.split("_")[1].split(".")[0])
            shear_info_path = data_path + "%s/result/%s_shear.dat"%(field,chip)
            try:
                shear_data = numpy.loadtxt(shear_info_path, skiprows=1)
                lb_data = numpy.zeros((len(shear_data), 3))
                lb_data[:, 0] = field_label
                lb_data[:, 1] = expo_label
                lb_data[:, 2] = chip_label
                if count == 0:
                    stack_data = copy.deepcopy(lb_data)
                else:
                    stack_data = numpy.row_stack((stack_data, lb_data))
                count += 1
            except:
                print("Empty %s/%s"%(field, chip))

sp = stack_data.shape
recv_sp = comm.gather(sp, root=0)
if rank > 0:
    comm.Send(stack_data, dest=0, tag=rank)
else:
    for procs in range(1, cpus):
        recvs = numpy.empty(recv_sp[procs], dtype=numpy.float64)
        comm.Recv(recvs, source=procs, tag=procs)
        stack_data = numpy.row_stack((stack_data, recvs))
    final_path = result_path + "label.npz"

    numpy.savez(final_path, stack_data)
    print("Totally, %d galaxies"%(len(stack_data)))





