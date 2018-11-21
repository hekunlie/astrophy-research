import platform
if platform.system() == 'Linux':
    import matplotlib
    matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/'%my_home)
path.append('E:/Github/astrophy-research/my_lib/')
import numpy
import tool_box
from mpi4py import MPI
import time
import h5py
from astropy.io import fits
import matplotlib.pyplot as plt

cmd = argv[1]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()
# cpus = 1
# rank = 0
log_path = "./log_%d.dat"%rank
logger = tool_box.get_logger(log_path)

nm_path = "/mw/w1234/original/nname.dat"
all_expos, all_fields = tool_box.field_dict(nm_path)

fields = tool_box.allot(all_fields, cpus)[rank]

chip_data_path = "/mw/w1234/original/"
result_path = "/lmc/cfht/para_fit/"


my, mx = numpy.mgrid[0:4644, 0:2112]
myf = my.flatten()
mxf = mx.flatten()
tags = numpy.arange(0, len(myf))
if cmd == "files":
    if rank == 0:
        for field_name in all_fields:
            pic_field_path = result_path + "pic/" + '%s/'%field_name
            if not os.path.exists(pic_field_path):
                os.mkdir(pic_field_path)
            data_field_path = result_path + "data/" + '%s/'%field_name
            if not os.path.exists(data_field_path):
                os.mkdir(data_field_path)
            # for expo_name in all_expos[field_name].keys():
            #     expo_path = field_path + "%s/"%expo_name
            #     if not os.path.exists(expo_path):
            #         os.mkdir(expo_path)

if cmd == "fit":
    t1 = time.time()
    final_data_list = []
    all_names = []
    fails = []
    count = 0
    bad_pix, saturated = 0.1, 50000
    for field_name in fields:
        for expo_name in all_expos[field_name].keys():
            # print(field_name, expo_name)
            ts = time.time()
            logger.info("RANK: %d %s -- %s starts..."%(rank, field_name, expo_name))
            # sigma, mean, amplitude, a1, a2 (x), a3 (y)
            noise_data = numpy.zeros((36, 6)) - 1

            pic_path = result_path + 'pic/%s/%s.png' % (field_name, expo_name)
            fig = plt.figure(figsize=(36, 16))
            for i in range(36):
                ax = fig.add_subplot(4, 9, i+1)
                chip_path = chip_data_path + "%s/science/%s_%d.fits"%(field_name, expo_name, i+1)
                if os.path.exists(chip_data_path):
                    # try:
                    img = fits.open(chip_path)[0].data
                    back_grd_info = tool_box.fit_background(img, 200000, "flat", bad_pix, saturated,myf,mxf,tags)[0][0]
                    noise_data[i, 3:6] = back_grd_info.reshape(1, 3)
                    img_zero = img - back_grd_info[0] - back_grd_info[1]*mx - back_grd_info[2]*my
                    noise_info = tool_box.fit_background(image=img_zero, pix_num=200000, function="gauss",
                                                         pix_lb=bad_pix-back_grd_info[0,0], pix_ub=saturated,
                                                         my=myf, mx=mxf, seqs=tags, ax=ax)[0]
                    noise_data[i,0:3] = noise_info[0][1], noise_info[0][0], noise_info[1][0]
                    # except:
                    #     print("RANK %02d:  CHIP %s/%s_%d.fits fitting failed !!" % (rank, field_name, expo_name, i))

                else:
                    print("RANK %02d:  CHIP %s_%d.fits does not exist!!" % (rank, expo_name, i))
            plt.subplots_adjust(wspace=0)
            plt.savefig(pic_path)
            plt.close()

            numpy.savez(result_path+'data/%s/%s.npz' % (field_name, expo_name), noise_data)

            if count == 0:
                final_data = noise_data.copy()
            else:
                final_data = numpy.row_stack((final_data, noise_data))
            te = time.time()
            logger.info("RANK: %d %s -- %s end, time: %.2f" % (rank, field_name, expo_name, te - ts))
            count += 1
            all_names.append("/%s/%s"%(field_name, expo_name))
    final_data_list.append(final_data)
    sp = final_data.shape
    sps = comm.gather(sp, root=0)
    all_names = comm.gather(all_names, root=0)

    if rank > 0:
        comm.Send(final_data, dest=0, tag=rank)
    else:
        for procs in range(1, cpus):
            recvs = numpy.empty(sps[procs], dtype=numpy.float64)
            comm.Recv(recvs, source=procs, tag=procs)
            final_data_list.append(recvs)

        h5file = h5py.File(result_path + "data/sigma.hdf5", "w")
        for i in range(len(all_names)):
            if i == 0:
                all_data = final_data_list[i].copy()
            else:
                all_data = numpy.row_stack((final_data_list[i], all_data))

            for tag, name in enumerate(all_names[i]):
                h5file[name] = final_data_list[i][tag*36: (tag+1)*36]
        h5file["/all"] = all_data
        h5file.close()

    t2 = time.time()
    if rank == 0:
        print(t2 - t1)
