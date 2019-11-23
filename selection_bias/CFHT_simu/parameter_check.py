import platform
if platform.system() == 'Linux':
    import matplotlib
    matplotlib.use('Agg')
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/'%my_home)
import numpy
import matplotlib.pyplot as plt
import tool_box
from mpi4py import MPI
import h5py
import time

total_path, loops = argv[1], int(argv[2])

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

para_path = total_path + "/parameters"
para_ini_path = para_path + "/para.ini"
paras = tool_box.config(para_ini_path,["get",'get',"get",'get',"get",'get','get'],
                        [["para","total_num","0"],["para","stamp_size","0"],
                        ["para", "mag_s", "0"],["para","mag_e","0"],
                        ["para", "radius_s", "0"],["para","radius_e","0"],
                         ["para", "shear_num", "0"]])

chip_num = int(paras[0])
size = int(paras[1])
num_i = int(chip_num/loops)
mag_s, mag_e = float(paras[2]), float(paras[3])
radius_s, radius_e = float(paras[4]), float(paras[5])
shear_num = int(paras[6])
stamp_num = 10000

seed_ini = numpy.random.randint(1, 1000000, size=cpus)[rank]
rng_ini = numpy.random.RandomState(seed_ini)
seeds = rng_ini.randint(0, 100000, size=2000)

if rank == 0:
    print(chip_num*stamp_num, size, mag_s, mag_e, radius_s, radius_e)
comm.Barrier()

# there are many shear points,
# the threads will process their own points
m, n = divmod(shear_num, cpus)
shear_st = m*rank
shear_ed = m*(rank+1)
if rank == cpus - 1:
    shear_ed += n


for shear_id in range(shear_st, shear_ed):
    h5_path = para_path+'/para_%d.hdf5'%shear_id

    f = h5py.File(h5_path,"w")

    # magnitude & flux
    flux, mag = numpy.zeros((chip_num*stamp_num, )),numpy.zeros((chip_num*stamp_num, ))
    for i in range(loops):
        time.sleep(rank * 0.05)
        mag_i = tool_box.mag_generator(num_i*stamp_num, mag_s, mag_e)
        flux_i = tool_box.mag_to_flux(mag_i)
        sp, ep = i*num_i*stamp_num, (i+1)*num_i*stamp_num
        mag[sp: ep] = mag_i
        flux[sp: ep] = flux_i

    plt.figure(figsize=(8,6))
    f["/flux"] = flux
    f["/mag"] = mag
    f.close()

    pic = para_path + "/pic/mag_%d.png"%shear_id
    plt.title("MAG")
    plt.hist(mag, 100)
    plt.savefig(pic)


comm.Barrier()
