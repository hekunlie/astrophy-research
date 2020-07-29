import platform
if platform.system() == 'Linux':
    import matplotlib
    matplotlib.use('Agg')
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/'%my_home)
import numpy
import tool_box
from mpi4py import MPI
import h5py



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

para_path = argv[1]

total_num = 20000000

seed = 1234 + rank*1212

rng = numpy.random.RandomState(seed)

# disc_e_i = tool_box.ran_generator(tool_box.bulge_e_pdf, total_num, seed, 0, 0.7, 0, 2.1)[0]
disc_e_i = rng.uniform(0, 0.266, total_num)
# disc_e_i = 0.6

theta = rng.uniform(0, 2*numpy.pi, total_num)
disc_e1_i, disc_e2_i = disc_e_i*numpy.cos(theta), disc_e_i*numpy.sin(theta)

e1 = numpy.zeros((total_num,),dtype=numpy.float32)
e2 = numpy.zeros((total_num,),dtype=numpy.float32)
e1[:] = disc_e1_i
e2[:] = disc_e2_i

h5f = h5py.File(para_path+"/para_%d.hdf5"%rank,"w")
h5f["/e1"] = e1
h5f["/e2"] = e2
h5f.close()

comm.Barrier()
# checking
# if it shows a "m" or "c", run the code again with different seed
if rank == 0:
    h5f = h5py.File(para_path + "/shear.hdf5", "r")
    g1_t = h5f["/g1"][()]
    g2_t = h5f["/g2"][()]
    h5f.close()
    mean_e = numpy.zeros((2, cpus))
    sig_e = numpy.zeros((2, cpus))
    for i in range(cpus):
        h5f = h5py.File(para_path + "/para_%d.hdf5" %i, "r")
        e1 = h5f["/e1"][()]
        e2 = h5f["/e2"][()]
        h5f.close()
        mean_e[0, i] = e1.mean()
        mean_e[1, i] = e2.mean()
        sig_e[0, i] = e1.std()/numpy.sqrt(total_num)
        sig_e[1, i] = e2.std()/numpy.sqrt(total_num)
    mc1 = tool_box.data_fit(g1_t, mean_e[0], sig_e[0])
    mc2 = tool_box.data_fit(g2_t, mean_e[1], sig_e[1])
    print(mean_e[0])
    print(sig_e[0])
    print(mean_e[1])
    print(sig_e[1])
    print("m1*10^3: %.5f(%.5f), c1*10^4: %.6f(%.6f)"%(mc1[0]*1000, mc1[1]*1000, mc1[2]*10000, mc1[3]*10000))
    print("m2*10^3: %.5f(%.5f), c2*10^4: %.6f(%.6f)"%(mc2[0]*1000, mc2[1]*1000, mc2[2]*10000, mc2[3]*10000))
comm.Barrier()