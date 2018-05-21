from sys import path
path.append('/home/hkli/work/fourier_quad')
import numpy
from subprocess import Popen
import os
from mpi4py import MPI
from astropy.io import fits
import time
import warnings
import tool_box
import shutil


warnings.filterwarnings("error")


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

t1 = time.time()
size = 48

with open("/home/hkli/work/envs/envs.dat", "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        data_path = path.split("=")[1]
    elif "cfht_field_path" in path:
        field_path = path.split("=")[1]

cfht_dict, fields = tool_box.field_dict(data_path + "nname.dat")
field_pool = tool_box.allot(fields, cpus)[rank]

empty_chip = []
snr_expos = []
for field in field_pool:
    expos = list(cfht_dict[field].keys())
    for expo in expos:
        chip_num = 0
        expo_path = field_path + "%s/%s/"%(field, expo)

        # initialize
        if os.path.exists(expo_path):
            shutil.rmtree(expo_path)
        os.mkdir(expo_path)

        for chip in list(cfht_dict[field][expo]):
            cat_name = data_path + "%s/stamps/%s_source_info.dat"%(field, chip)
            try:
                source_info = numpy.loadtxt(cat_name, skiprows=1)
                img_path = data_path + "%s/stamps/%s_source.fits" %(field, chip)
                sex_cat = expo_path + "%s.cat" % chip
                cmd = "sex %s -CATALOG_NAME %s" % (img_path, sex_cat)
                a = Popen(cmd, shell=True)
                a.wait()

                arr = numpy.loadtxt(sex_cat)
                a_snr = arr[:, 0]
                a_mag = arr[:, 1]
                a_area = arr[:, 2]
                a_x = arr[:, 3] - 1
                a_y = arr[:, 4] - 1

                fits_img = fits.open(img_path)
                h, w = fits_img[0].data.shape
                fits_img.close()
                row, col = int(h / size), int(w / size)
                sex_data = numpy.zeros((len(source_info), 3))
                for m in range(len(arr)):
                    xx = a_x[m]
                    yy = a_y[m]
                    snr = a_snr[m]
                    area = a_area[m]
                    mag = a_mag[m]
                    mx, modx = divmod(xx, size)
                    my, mody = divmod(yy, size)
                    tag = int(col * my + mx)

                    distance = numpy.sqrt((modx - size / 2) ** 2 + (mody - size / 2) ** 2)
                    if area >= 5 and distance <= 5:
                        if area > sex_data[tag, 1]:
                            sex_data[tag, 0] = snr
                            sex_data[tag, 1] = area
                            sex_data[tag, 2] = mag
                if chip_num == 0:
                    stack_data = sex_data
                else:
                    stack_data = numpy.row_stack((stack_data, sex_data))
                chip_num += 1
            except:
                excluded = "%s/%s\n"%(field, chip)
                print("Empty :%s"%excluded)
                empty_chip.append(excluded)

        if len(stack_data) > 1:
            sex_npz = field_path + "%s/%s_sex.npz"%(field, expo)
            numpy.savez(sex_npz, stack_data)
            snr_expos.append(stack_data)

for i in range(len(snr_expos)):
    if i == 0:
        data = snr_expos[i]
    else:
        data = numpy.row_stack((data, snr_expos[i]))

recv_emp_chip = comm.gather(empty_chip, root=0)

sp = data.shape
recv_sp = comm.gather(sp, root=0)

if rank > 0:
    comm.Send(data, dest=0, tag=rank)
else:
    for procs in range(1, cpus):
        recvs = numpy.empty(recv_sp[procs], dtype=numpy.float64)
        comm.Recv(recvs, source=procs, tag=procs)
        data = numpy.row_stack((data, recvs))
    snr_final_path = field_path + "sex_snr.npz"
    numpy.savez(snr_final_path, data)

t2 = time.time()

if rank == 0:
    gather_empty = []
    for chips in recv_emp_chip:
        gather_empty.extend(chips)
    ex_chip_path = field_path + "exclude_chips.dat"
    if not os.path.exists(ex_chip_path):
        with open(ex_chip_path, "w") as f:
            f.writelines(gather_empty)

    print((t2 - t1) / 3600)