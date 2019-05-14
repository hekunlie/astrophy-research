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

sect, source,loops = argv[1], argv[2], int(argv[3])

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

envs_path = "%s/work/envs/envs.dat"%my_home
para_path = tool_box.config(envs_path,["get"],[["%s"%sect,"%s_path_para"%source,"0"]])[0]

para_ini_path = para_path+"para.ini"
paras = tool_box.config(para_ini_path,["get",'get',"get",'get',"get",'get'],
                        [["para","total_num","0"],["para","stamp_size","0"],
                        ["para", "mag_s", "0"],["para","mag_e","0"],
                        ["para", "radius_s", "0"],["para","radius_e","0"]])

num = int(paras[0])
size = int(paras[1])
num_i = int(num/loops)
mag_s, mag_e = float(paras[2]),float(paras[3])
radius_s, radius_e = float(paras[4]), float(paras[5])
stamp_num = 10000
if rank == 0:
    print(num*stamp_num, size, mag_s, mag_e, radius_s, radius_e)
comm.Barrier()
pic = para_path + "/pic/ellip_%d.png"%rank
plt.figure(figsize=(16,16))

h5_path = para_path+'para_%d.hdf5'%rank
para_logs_path = para_path + "logs/logs_%d.dat"%rank
logger = tool_box.get_logger(para_logs_path)
log_inform = "RANK: %d, LOOP: %d, TOTAL NUM: %d, NUM in LOOP: %d, " \
             "SIZEL %d, MAG: %f ~ %f, RADIUS: %.2f ~ %.2f\n"%(rank, loops, num, num_i, size, mag_s, mag_e, radius_s, radius_e)
logger.info(log_inform)
f = h5py.File(h5_path,"w")


itemsize = MPI.DOUBLE.Get_size()
element_num = cpus
if rank == 0:
    # bytes for 10 double elements
    nbytes = element_num*itemsize
else:
    nbytes = 0

# on rank 0 of comm, create the contiguous shared block
win1 = MPI.Win.Allocate_shared(nbytes*2, itemsize, comm=comm)
win2 = MPI.Win.Allocate_shared(nbytes*2, itemsize, comm=comm)
# create a numpy array whose data points to the shared block
# buf is the block's address in the memory
buf1, itemsize = win1.Shared_query(0)
buf2, itemsize = win2.Shared_query(0)

# create a numpy array from buf
# code can run successfully without the following step
# buf = np.array(buf, dtype='float64', copy=False) # may be redundant
# "d" means double = 'float64'
mean_e = numpy.ndarray(buffer=buf1, dtype='d', shape=(element_num,2)) # array filled with zero
sig_e = numpy.ndarray(buffer=buf2, dtype='d', shape=(element_num,2))


seed_ini = numpy.random.randint(1, 100000, size=cpus)[rank]
rng_ini = numpy.random.RandomState(seed_ini)
seeds = rng_ini.randint(0, 100000, size=2000)

# ellipticity
e1, e2, e = numpy.zeros((num*stamp_num, 1)), numpy.zeros((num*stamp_num, 1)),numpy.zeros((num*stamp_num, 1))
gal_type = numpy.zeros((num*stamp_num, 1))
disc_frac = 0.9
bulge_frac = 0.1
disc_num_i = int(disc_frac*num_i)*stamp_num
disc_num = int(disc_frac*num)*stamp_num
bulge_num_i = int(bulge_frac*num_i)*stamp_num
bulge_num = int(bulge_frac*num)*stamp_num
# disc-dominated galaxies

pesb = numpy.linspace(0, 0.8, int(0.8 / 0.0000001) + 1)

counts = 0
if disc_frac > 0:

    pef = tool_box.disc_e_pdf([pesb])
    pe = pef / pef.sum()

    for i in range(loops):

        rng = numpy.random.RandomState(seeds[counts])
        counts += 1

        disc_e_i = rng.choice(pesb, disc_num_i, p=pe)

        # disc_e_i = tool_box.ran_generator(tool_box.disc_e_pdf, disc_num_i, seed, 0, 0.804, 0, 2.1)[0]

        # disc_e_i = rng.uniform(0, 0.9, disc_num_i)

        theta = rng.uniform(0, 2*numpy.pi, disc_num_i)
        disc_e1_i, disc_e2_i = disc_e_i*numpy.cos(theta), disc_e_i*numpy.sin(theta)

        sp, ep = i*disc_num_i, (i + 1)*disc_num_i
        e1[sp:ep, 0] = disc_e1_i
        e2[sp:ep, 0] = disc_e2_i
        e[sp:ep, 0] = disc_e_i

        log_inform = "Disc loop: %d, %d: %d, seed: %d, mean(e1): %.4f, std(e1): %.4f, mean(e2): %.4f, std(e2): %.4f, " \
                     "max(e1): %.4f, max(e2): %.4f\n"\
                     %(i, sp, ep, seeds[counts], e1[sp:ep, 0].mean(), e1[sp:ep, 0].std(), e2[sp:ep, 0].mean(), e2[sp:ep, 0].std(),
                       e1[sp:ep, 0].max(), e2[sp:ep, 0].max())
        logger.info(log_inform)
    plt.subplot(331)
    plt.hist(e1[:disc_num], 100)
    plt.title("Disc e1")
    plt.subplot(332)
    plt.hist(e2[:disc_num], 100)
    plt.title("Disc e2")
    plt.subplot(333)
    plt.hist(e[:disc_num], 100)
    plt.title("Disc e")

# bulge-dominated galaxies
if bulge_num > 0:

    pef = tool_box.bulge_e_pdf([pesb])
    pe = pef / pef.sum()

    for i in range(loops):

        rng = numpy.random.RandomState(seeds[counts])
        counts += 1

        bulge_e_i = rng.choice(pesb, bulge_num_i, p=pe)

        # bulge_e_i = tool_box.ran_generator(tool_box.bulge_e_pdf, bulge_num_i, seed, 0, 0.804, 0, 2.7)[0]

        # bulge_e_i = rng.uniform(0, 0.9, bulge_num_i)

        theta = rng.uniform(0, 2*numpy.pi, bulge_num_i)
        bulge_e1_i, bulge_e2_i = bulge_e_i*numpy.cos(theta), bulge_e_i*numpy.sin(theta)

        sp, ep = disc_num + i*bulge_num_i, disc_num + (i+1)*bulge_num_i
        e1[sp:ep, 0] = bulge_e1_i
        e2[sp:ep, 0] = bulge_e2_i
        e[sp:ep, 0] = bulge_e_i
        gal_type[sp:ep, 0] = 1

        log_inform = "Bulge loop: %d, %d: %d, seed: %d, mean(e1): %.4f, std(e1): %.4f, mean(e2): %.4f, std(e2): %.4f, " \
                     "max(e1): %.4f, max(e2): %.4f\n"\
                     %(i, sp, ep, seeds[counts], e1[sp:ep, 0].mean(), e1[sp:ep, 0].std(), e2[sp:ep, 0].mean(), e2[sp:ep, 0].std(),
                       e1[sp:ep, 0].max(), e2[sp:ep, 0].max())
        logger.info(log_inform)

    plt.subplot(334)
    plt.hist(e1[disc_num:], 100)
    plt.title("Bulge e1")
    plt.subplot(335)
    plt.hist(e2[disc_num:], 100)
    plt.title("Bulge e2")
    plt.subplot(336)
    plt.hist(e[disc_num:], 100)
    plt.title("Bulge e")
for i in range(cpus):
    if i == rank:
        print("Rank: %3d, mean(e1): %10.6f, std: %.4f, mean(e2): %10.6f, std: %.4f, max: %.5f, %.5f"
              %(rank, numpy.mean(e1), numpy.std(e1), numpy.mean(e2), numpy.std(e2), numpy.max(e1), numpy.max(e2)))
        comm.Barrier()
f["/e1"] = e1
f["/e2"] = e2
f["/e"] = e
f["/type"] = gal_type

mean_e[rank,0] = e1.mean()
mean_e[rank,1] = e2.mean()
sig_e[rank, 0] = e1.std()
sig_e[rank, 1] = e2.std()

# magnitude & flux
flux, mag = numpy.zeros((num*stamp_num, 1)),numpy.zeros((num*stamp_num, 1))
for i in range(loops):
    time.sleep(rank * 0.05)
    mag_i = tool_box.mag_generator(num_i*stamp_num, mag_s, mag_e)
    flux_i = tool_box.mag_to_flux(mag_i)
    sp, ep = i*num_i*stamp_num, (i+1)*num_i*stamp_num
    mag[sp: ep, 0] = mag_i
    flux[sp: ep, 0] = flux_i

    log_inform = "MAGNITUDE loop: %d, min(mag): %.2f, max(mag): %.2f\n"%(i, mag_i.min(), mag_i.max())
    logger.info(log_inform)

f["/flux"] = flux
f["/mag"] = mag
plt.subplot(337)
plt.title("MAG")
plt.hist(mag, 100)

# galactic radius
radius = numpy.zeros((num*stamp_num, 1))
for i in range(loops):

    time.sleep(rank*0.05)
    rng = numpy.random.RandomState(seeds[counts])
    counts += 1

    # radius_i = rng.uniform(radius_s, radius_e, num_i*stamp_num)

    sp, ep = i * num_i * stamp_num, (i + 1) * num_i * stamp_num
    radius_i = tool_box.radii_from_mags(mag[sp: ep, 0], radius_s, radius_e)

    radius[sp: ep, 0] = radius_i

    log_inform = "RADIUS loop: %d, min: %.2f, max: %.2f\n"%(i, radius_i.min(), radius_i.max())
    logger.info(log_inform)
f["/radius"] = radius
plt.subplot(338)
plt.title("Radius")
plt.hist(radius, 100)

# B/T ratio
btr = numpy.zeros((num*stamp_num, 1))
for i in range(loops):

    f_btr_i = tool_box.ran_generator(tool_box.bulge_frac_pdf, num_i*stamp_num, seeds[counts], 0, 1, 0, 4.1)[0]
    counts += 1

    sp, ep = i*num_i*stamp_num, (i + 1)*num_i*stamp_num
    btr[sp: ep, 0] = f_btr_i

    log_inform = "B/T loop: %d, seed: %d, min: %.2f, max: %.2f\n"%(i, seeds[counts], f_btr_i.min(), f_btr_i.max())
    logger.info(log_inform)
plt.subplot(339)
plt.title("B/T")
plt.hist(btr, 100)
plt.savefig(pic)
plt.close()

f["/btr"] = btr
f.close()

comm.Barrier()
if rank == 0:
    shears = numpy.loadtxt(para_path + "shear.dat")
    g1_t = shears[:cpus]
    g2_t = shears[cpus:2*cpus]

    mc1 = tool_box.data_fit(g1_t, mean_e[:,0], sig_e[:,0]/numpy.sqrt(num*stamp_num))
    mc2 = tool_box.data_fit(g2_t, mean_e[:,1], sig_e[:,1]/numpy.sqrt(num*stamp_num))
    print(mean_e[:,0])
    print(sig_e[:,0])
    print(mean_e[:,1])
    print(sig_e[:,1])
    print("m1: %.5f(%.5f), c1: %.6f(%.6f)"%(mc1[0], mc1[1], mc1[2], mc1[3]))
    print("m2: %.5f(%.5f), c2: %.6f(%.6f)"%(mc2[0], mc2[1], mc2[2], mc2[3]))
comm.Barrier()