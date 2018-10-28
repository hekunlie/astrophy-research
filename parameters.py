import platform
if platform.system() == 'Linux':
    import matplotlib
    matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/'%my_home)
import numpy
import matplotlib.pyplot as plt
import tool_box
import lsstetc
from mpi4py import MPI
import h5py
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

loops, source = int(argv[1]), argv[2]

envs_path = "%s/work/envs/envs.dat"%my_home
para_path = tool_box.config(envs_path,["get"],[["selection_bias","%s_path_para"%source,"0"]])[0]
para_ini_path = para_path+"para.ini"
paras = tool_box.config(para_ini_path,["get",'get',"get",'get',"get",'get'],
                        [["para","total_num","0"],["para","size","0"],
                        ["para", "mag_s", "0"],["para","mag_e","0"],
                        ["para", "radius_s", "0"],["para","radius_e","0"]])
num = int(paras[0])
size = int(paras[1])
num_i = int(num/loops)
mag_s, mag_e = float(paras[2]),float(paras[3])
radius_s, radius_e = float(paras[4]), float(paras[5])
if rank == 0:
    print(num, size, mag_s, mag_e, radius_s, radius_e)

prop = lsstetc.ETC(band='r', pixel_scale=0.2, stamp_size=size, nvisits=180)
path = para_path+'para_%d.hdf5'%rank
f = h5py.File(path,"w")

# def ran(start, end, num):
#     return numpy.random.uniform(start, end, num)
#
#
# g1 = numpy.append(numpy.append(ran(-0.02, 0, 4), ran(0, 0.021, 3)),numpy.append(ran(-0.02, 0, 3), ran(0, 0.021, 4)))
# g2 = numpy.append(numpy.append(ran(-0.02, 0, 3), ran(0, 0.021, 4)),numpy.append(ran(0, 0.021, 3), ran(-0.02, 0, 4)))
# # numpy.random.shuffle(g1)
# # numpy.random.shuffle(g2)
# plt.subplot(131)
# plt.scatter(g1,g2)
# plt.subplot(132)
# plt.scatter(g1,g1)
# plt.subplot(133)
# plt.scatter(g2,g2)
# plt.show()
# numpy.savez('E:/selection_bias/parameters/shear.npz', g1, g2)
# numpy.savetxt('E:/selection_bias/parameters/shear.dat',numpy.append(g1,g2))


# ellipticity
e1, e2, e, es = [], [], [], []
for i in range(loops):
    seed = rank * 43254 + int(numpy.random.randint(1, 1256542344, 1)[0])
    time.sleep(rank*0.1+0.5)
    e1_i, e2_i, e_i, es_i = tool_box.ellip_mock(num_i, seed)
    e1.extend(e1_i.tolist())
    e2.extend(e2_i.tolist())
    e.extend(e_i.tolist())
    es.extend(es_i.tolist())

e1 = numpy.array(e1)
e2 = numpy.array(e2)
e = numpy.array(e)
es = numpy.array(es)
print("Rank: %3d, mean(e1): %10.6f, std: %.4f, mean(e2): %10.6f, std: %.4f, max: %.5f, %.5f"
      %(rank, numpy.mean(e1)[0], numpy.std(e1)[0], numpy.mean(e2)[0], numpy.std(e2)[0], numpy.max(e1), numpy.max(e2)))
f["/e1"] = e1
f["/e2"] = e2
f["/e"] = e
f["/es"] = es

# magnitude & flux
flux = []
mag = []
for i in range(loops):
    time.sleep(rank * 0.1 + 0.5)
    mag_i = tool_box.mags_mock(num_i, mag_s, mag_e)
    mag.extend(mag_i.tolist())
    flux_i = prop.flux(mag_i).tolist()
    flux.extend(flux_i)
flux = numpy.array(flux)
mag = numpy.array(mag)
f["/flux"] = flux
f["/mag"] = mag

# galactic radius
radius = []
for i in range(loops):
    seed = rank * 43254 + int(numpy.random.randint(1, 1256542344, 1)[0])
    time.sleep(rank*0.1+0.5)
    rng = numpy.random.RandomState(seed=seed)
    radius_i = rng.uniform(radius_s, radius_e, num_i)
    radius.extend(radius_i.tolist())
radius = numpy.array(radius)
f["/radius"] = radius

# B/T ratio
btr = []
for i in range(loops):
    seed = rank * 43254 + int(numpy.random.randint(1, 1256542344, 1)[0])
    time.sleep(rank*0.1+0.5)
    rng = numpy.random.RandomState(seed=seed)
    btr_i = rng.normal(0, 0.1, 2*num_i)
    btr_i.shape = (len(btr_i), 1)
    idx1 = btr_i >= 0
    idx2 = btr_i < 1
    c_i = btr_i[idx1&idx2]
    c_i.shape = (len(c_i), 1)
    if len(c_i) > num_i:
        f_btr_i = c_i[:num_i]
    elif len(c_i) < num_i:
        while True:
            gap = num_i - len(c_i)
            if gap == 0:
                f_btr_i = c_i
                break
            plus = rng.normal(0, 0.1, 10*gap)
            plus = plus[plus >= 0][:gap]
            plus.shape = (len(plus), 1)
            c_i = numpy.row_stack((c_i, plus))
    else:
        f_btr_i = c_i
    btr.extend(f_btr_i.tolist())
f["/btr"] = btr
f.close()

pic = para_path + "/pic/ellip_%d.png"%rank
plt.figure(figsize=(16,16))
plt.subplot(221)
plt.title("e1")
plt.hist(e1, 100)
plt.subplot(222)
plt.title("e2")
plt.hist(e2, 100)
plt.subplot(223)
plt.title("e")
plt.hist(e, 100)
plt.subplot(224)
plt.title("es")
plt.hist(es, 100)
plt.savefig(pic)
plt.close()
#
pic = para_path + "/pic/fmrb_%d.png"%rank
plt.figure(figsize=(16,16))
plt.subplot(221)
plt.title("MAG")
plt.hist(mag, 100)
plt.subplot(222)
plt.title("FLUX")
plt.hist(flux, 100)
plt.subplot(223)
plt.title("RADIUS")
plt.hist(radius, 100)
plt.subplot(224)
plt.title("BTR")
plt.hist(btr, 100)
plt.savefig(pic)
plt.close()


