import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append("%s/work/mylib/"% my_home)
import h5py
import numpy
import numpy.ctypeslib as ctl
from mpi4py import MPI
import tool_box
import ctypes
import time


lib = ctypes.cdll.LoadLibrary("libhist.so")

deblending_self = lib.deblending_self
deblending_mutual = lib.deblending_mutual


deblending_self.restype = None
deblending_mutual.restype = None
deblending_self.argtypes = [ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                            ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                            ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                            ctypes.c_int,
                            ctypes.c_float,
                            ctypes.c_float,
                            ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous')]

deblending_mutual.argtypes = [ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctypes.c_int,
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctypes.c_int,
                                ctypes.c_float,
                                ctypes.c_float,
                                ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous')]


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

parent_path = argv[1]

# <= sep_deg & >= dz
sep_deg = 0.0013
dz = 0.5

ra_col = 0
dec_col = 1
z_col = 10


file_name = []
pos = []
with open(parent_path + "/cat_inform/exposure_avail_inform.dat","r") as f:
    cc = f.readlines()
for tag, row in enumerate(cc):
    if tag > 0:
        iterms = row.split("\n")[0].split("\t")
        file_name.append(iterms[0])
        temp = []
        for i in range(1, len(iterms)-1):
            temp.append(float(iterms[i]))
pos_arr = numpy.array(pos)
file_num = pos_arr.shape[0]

task_label = [i for i in range(file_num)]

task_label_sub = tool_box.alloc(task_label, numprocs)[rank]

for i in task_label_sub:

    t1 = time.time()

    h5f = h5py.File(parent_path + "/cat_hdf5/%s.hdf5"%file_name[i],"r")
    data = h5f["/data"][()]
    h5f.close()

    ra1 = numpy.ascontiguousarray(data[:,ra_col], dtype=numpy.float32)
    dec1 = numpy.ascontiguousarray(data[:,dec_col], dtype=numpy.float32)
    z1 = numpy.ascontiguousarray(data[:,z_col], dtype=numpy.float32)

    ra1_cent = ra1.mean()
    dec1_cent = dec1.mean()
    src_num1 = z1.shape[0]

    blended_label = numpy.zeros_like(z1, dtype=numpy.intc)

    deblending_self(ra1, dec1, z1, src_num1, sep_deg, dz, blended_label)

    for j in range(file_num):
        if i != j:
            diff = numpy.abs(pos_arr[j] - pos_arr[i])
            if diff.max() < 1:
                pass
            else:
                pass
            h5f = h5py.File(files_sub[j], "r")
            data = h5f["/data"][()]
            h5f.close()

            ra2 = numpy.ascontiguousarray(data[:, ra_col], dtype=numpy.float32)
            dec2 = numpy.ascontiguousarray(data[:, dec_col], dtype=numpy.float32)
            z2 = numpy.ascontiguousarray(data[:, z_col], dtype=numpy.float32)

            ra2_cent = ra2.mean()
            dec2_cent = dec2.mean()
            src_num2 = z1.shape[0]

            deblending_mutual(ra1, dec1, z1, src_num1, ra2, dec2, z2, src_num2,sep_deg, dz, blended_label)

    label_path = files_sub[i].replace("cat_hdf5", "blended_label")
    h5f = h5py.File(label_path,"w")
    h5f["/blended_label"] = blended_label
    h5f.close()

    t2 = time.time()
    print("rank %d. %.2f %s"%(rank, t2-t1, files_sub[i]))

comm.Barrier()