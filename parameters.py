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

loops, source = int(argv[1]), argv[2]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

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
h5_path = para_path+'para_%d.hdf5'%rank
para_logs_path = para_path + "logs/logs_%d.dat"%rank
log_inform = "RANK: %d, LOOP: %d, TOTAL NUM: %d, NUM in LOOP: %d, " \
             "SIZEL %d, MAG: %d ~ %d, RADIUS: %.2f ~ %.2f\n"%(rank, loops, num, num_i, size, mag_s, mag_e, radius_s, radius_e)

tool_box.write_log(log_path=para_logs_path, content=log_inform, way='direct')
f = h5py.File(h5_path,"w")


# g1 = numpy.append(numpy.append(ran(-0.02, 0, 4), ran(0, 0.021, 3)),numpy.append(ran(-0.02, 0, 3), ran(0, 0.021, 4)))
# g2 = numpy.append(numpy.append(ran(-0.02, 0, 3), ran(0, 0.021, 4)),numpy.append(ran(0, 0.021, 3), ran(-0.02, 0, 4)))
# numpy.random.shuffle(g1)
# numpy.random.shuffle(g2)
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
e1, e2, e = numpy.zeros((num, 1)),numpy.zeros((num, 1)),numpy.zeros((num, 1))
for i in range(loops):
    seed = rank * 43254 + int(numpy.random.randint(1, 12565, 1)[0])

    time.sleep(rank*0.05)
    e1_i, e2_i, e_i = tool_box.ellip_mock(num_i, seed)
    e1[i*num_i: (i+1)*num_i, 0] = e1_i
    e2[i*num_i: (i+1)*num_i, 0] = e2_i
    e[i*num_i: (i+1)*num_i, 0] = e_i

    log_inform = "ELLIP loop: %d, seed: %d, mean(e1): %.4f, std(e1): %.4f, mean(e2): %.4f, std(e2): %.4f, " \
                 "max(e1): %.4f, max(e2): %.4f\n"\
                 %(i, seed, e1_i.mean(), e1_i.std(), e2_i.mean(), e2_i.std(), e1_i.max(), e2_i.max())
    tool_box.write_log(log_path=para_logs_path, content=log_inform, way="direct")

print("Rank: %3d, mean(e1): %10.6f, std: %.4f, mean(e2): %10.6f, std: %.4f, max: %.5f, %.5f"
      %(rank, numpy.mean(e1), numpy.std(e1), numpy.mean(e2), numpy.std(e2), numpy.max(e1), numpy.max(e2)))
f["/e1"] = e1
f["/e2"] = e2
f["/e"] = e


# magnitude & flux
flux, mag = numpy.zeros((num, 1)),numpy.zeros((num, 1))
for i in range(loops):
    time.sleep(rank * 0.05)
    mag_i = tool_box.mags_mock(num_i, mag_s, mag_e)
    flux_i = prop.flux(mag_i)
    mag[i*num_i: (i+1)*num_i, 0] = mag_i
    flux[i*num_i: (i+1)*num_i, 0] = flux_i

    log_inform = "MAGNITUDE loop: %d, min(mag): %.2f, max(mag): %.2f\n"%(i, mag_i.min(), mag_i.max())
    tool_box.write_log(para_logs_path, log_inform, "direct")

f["/flux"] = flux
f["/mag"] = mag

# galactic radius
radius = numpy.zeros((num, 1))
for i in range(loops):
    seed = rank * 43254 + int(numpy.random.randint(1, 125654, 1)[0])
    time.sleep(rank*0.05)
    rng = numpy.random.RandomState(seed=seed)
    radius_i = rng.uniform(radius_s, radius_e, num_i)
    radius[i*num_i: (i+1)*num_i, 0] = radius_i

    log_inform = "RADIUS loop: %d, min: %.2f, max: %.2f\n"%(i, radius_i.min(), radius_i.max())
    tool_box.write_log(para_logs_path, log_inform, "direct")
f["/radius"] = radius

# B/T ratio
btr = numpy.zeros((num, 1))
for i in range(loops):
    seed = rank * 43254 + int(numpy.random.randint(1, 125654, 1)[0])
    time.sleep(rank*0.05)
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
    f_btr_i.shape = (num_i, 1)
    btr[i*num_i: (i+1)*num_i] = f_btr_i

    log_inform = "B/T loop: %d, seed: %d, min: %.2f, max: %.2f\n"%(i, seed, f_btr_i.min(), f_btr_i.max())
    tool_box.write_log(para_logs_path, log_inform, "direct")

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


