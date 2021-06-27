import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
from mpi4py import MPI
import tool_box
import FQlib
import time




comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

data_path = argv[1]

if rank == 0:
    h5f = h5py.File(data_path + "/bias_test.hdf5", "w")
    h5f.close()
comm.Barrier()

for i in range(4):

    h5f = h5py.File(data_path + "/group_predict_%d.hdf5"%i,"r")
    label = h5f["/data"][()]
    h5f.close()

    h5f = h5py.File(data_path + "/stack_data_%d.hdf5"%i, "r")
    data = h5f["/data"][()][:,:4]
    h5f.close()

    label_min, label_max = label.min(), label.max()

    label_list = [i for i in range(label_min, label_max+1)]

    label_list_sub = tool_box.alloc(label_list, cpus)[rank]

    results = numpy.zeros((len(label_list_sub), 5))
    sp = (len(label_list_sub), 5)

    for tag, ib in enumerate(label_list_sub):
        idx = label == ib

        data_sub = data[idx]

        mg1, mg2, mn, mu = data_sub[:,0],data_sub[:,1],data_sub[:,2],data_sub[:,3]

        gh1, gh1_sig = FQlib.find_shear_cpp(data_sub[:,0], data_sub[:,2] + data_sub[:,3],20)[:2]
        gh2, gh2_sig = FQlib.find_shear_cpp(data_sub[:,1], data_sub[:,2] - data_sub[:,3],20)[:2]

        results[tag] = ib, gh1, gh1_sig, gh2, gh2_sig

    comm.Barrier()

    total_sp = comm.gather(sp,root=0)
    comm.Barrier()
    if rank > 0:
        # !!!! remember the data type, MPI.DOUBLE, MPI.FLOAT, ...
        # or it will raise an error, Keyerror
        comm.Send([results, MPI.DOUBLE], dest=0, tag=rank)
    else:
        # receive the data from other CPUs
        # !!!! the start points is 1 in range() not 0
        for procs in range(1, cpus):
            # prepare a buffer for the data, the shape must be the same
            # with that of what the other CPUs send, you have collected them in 'data_sps'
            recvs = numpy.empty(total_sp[procs], dtype=numpy.double)
            # receive it using the buffer,
            comm.Recv(recvs, source=procs, tag=procs)
            # then do whatever you want ...
            results = numpy.row_stack((results, recvs))

        h5f = h5py.File(data_path + "/bias_test.hdf5","r+")
        print(results)
        h5f["/%d"%i] = results
        h5f.close()
    comm.Barrier()


