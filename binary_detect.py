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
count_0 = 1
f = []
for field in r_fields:
    expos = list(field_dict[field].keys())
    for expo in expos:

        if argv[1] == "find":
            count = 0
            chips = field_dict[field][expo]
            for chip in chips:
                chip_path = data_path + "%s/stamps/%s_source.fits"%(field,chip)
                chip_info_path = data_path + "%s/stamps/%s_source_info.dat"%(field,chip)
                shear_info_path = data_path + "%s/result/%s_shear.dat"%(field,chip)
                try:
                    shear_data = numpy.loadtxt(shear_info_path, skiprows=1)
                    sigs = numpy.loadtxt(chip_info_path, skiprows=1)[:, 3]
                    f = fits.open(chip_path)
                    img = f[0].data
                    pool = fq.segment(img)
                    binary_tag = numpy.zeros((len(pool),1))
                    binary_pool = []
                    binary_source_pool = []
                    f.close()

                    for i in range(len(pool)):
                        objs, peaks = tool_box.find_binary(pool[i], size, size, sigs[i])[0:2]
                        num = len(objs)
                        if num > 1 or (num == 1 and len(peaks[0]) > 1):
                            arr = numpy.zeros((size,size))
                            for j in range(num):
                                for p in objs[j]:
                                    arr[p[0]-2, p[1]-2] = 1
                                for p in peaks[j]:
                                    arr[p[0]-2, p[1]-2] = 3
                            binary_pool.append(arr)
                            binary_source_pool.append(pool[i])
                            binary_tag[i, 0] = 1

                    if count == 0:
                        binary_data = copy.deepcopy(binary_tag)
                    else:
                        binary_data = numpy.row_stack((binary_data, binary_tag))
                    count += 1

                    if len(binary_pool) > 0:
                        binary_img = fq.stack(binary_pool, 15)
                        binary_img_path = field_path + "%s/%s/bi_%s.fits"%(field, expo, chip)
                        hdu = fits.PrimaryHDU(binary_img)
                        hdu.writeto(binary_img_path, overwrite=True)

                        binary_img = fq.stack(binary_source_pool, 15)
                        binary_img_path = field_path + "%s/%s/bi_%s_source.fits"%(field, expo, chip)
                        hdu = fits.PrimaryHDU(binary_img)
                        hdu.writeto(binary_img_path, overwrite=True)
                except:
                    print("Empty %s/%s"%(field, chip))

            binary_tag_path = field_path + "%s/bi_%s.npz" % (field, expo)
            numpy.savez(binary_tag_path, binary_data)

        elif argv[1] == 'stack':
            bi_path = field_path + "%s/bi_%s.npz"%(field, expo)
            shear_info_path = field_path + "%s/%s.npz"%(field, expo)

            data = numpy.load(bi_path)["arr_0"]
            shear_cat = numpy.load(shear_info_path)["arr_0"]
            if data.shape[0] != shear_cat.shape[0]:
                print(field, expo,data.shape,shear_cat.shape)
                if field not in f:
                    f.append(field)

            # try:
            #     data = numpy.load(bi_path)["arr_0"]
            #     shear_cat = numpy.load(shear_info_path)["arr_0"]
            #     if data.shape[0] == shear_cat.shape[0]:
            #         if count_0 == 0:
            #             stack_data = data
            #         else:
            #             stack_data = numpy.row_stack((stack_data, data))
            #         count_0 += 1
            #     else:
            #         print(field,"The binary catalog doesn't match the shear catalog!")
            #         exit()
            # except:
            #     print(field,expo,"file doesn't exist!")
        else:
            print("Wrong input!")
            exit()
print(len(field))
# if argv[1] == "stack":
#     sp = stack_data.shape
#     recv_sp = comm.gather(sp, root=0)
#     if rank > 0:
#         comm.Send(stack_data, dest=0, tag=rank)
#     else:
#         for procs in range(1, cpus):
#             recvs = numpy.empty(recv_sp[procs], dtype=numpy.float64)
#             comm.Recv(recvs, source=procs, tag=procs)
#             stack_data = numpy.row_stack((stack_data, recvs))
#         final_path = result_path + "binary_label.npz"
#         numpy.savez(final_path, stack_data)




