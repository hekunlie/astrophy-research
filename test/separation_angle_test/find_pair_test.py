import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import time
import tool_box
import h5py
import numpy
from mpi4py import MPI
from astropy.coordinates import SkyCoord
from astropy import units



area_id = "/w_" +argv[1]
radius_id = int(argv[2])
cmd = argv[3]

# read radius bin
h5f = h5py.File("./data/result/radius_bin.hdf5", "r")
radius_bin = h5f["/radius_bin"].value[:, 0]
h5f.close()

if cmd == "find":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()
    h = 0.7
    C_0_hat = 2.99792458
    H_0 = 100*h
    coeff = 1000*C_0_hat

    # read foreground
    h5f = h5py.File("./data/%s_sub.hdf5"%area_id,"r")
    RA_f = h5f["/RA"].value
    DEC_f = h5f["/DEC"].value
    Z_f = h5f["/Z"].value
    DIST_f = h5f["/DISTANCE"].value[:,0]
    h5f.close()
    fore_num = Z_f.shape[0]

    tasks = [i for i in range(fore_num)]
    my_task = tool_box.allot(tasks, cpus)[rank]
    print("Rank: %d. My: %d ~ %d"%(rank, min(my_task), max(my_task)))

    # read background
    # area_id = "/w_1"
    background_path = "/mnt/perc/hklee/CFHT/gg_lensing/data/cfht_cata_grid.hdf5"
    h5f = h5py.File(background_path,"r")
    Z = h5f["%s/Z"%area_id].value[:,0]

    back_num = Z.shape[0]
    back_label = numpy.arange(back_num)

    BACK_DATA = numpy.zeros((Z.shape[0], 8))

    BACK_DATA[:,0] = back_label
    BACK_DATA[:,1] = h5f["%s/RA"%area_id].value[:,0]
    BACK_DATA[:,2] = h5f["%s/DEC"%area_id].value[:,0]
    BACK_DATA[:,3] = h5f["%s/Z"%area_id].value[:,0]
    BACK_DATA[:,4] = h5f["%s/DISTANCE"%area_id].value[:,0]
    h5f.close()

    backgal_pos = SkyCoord(ra=BACK_DATA[:,1] * units.deg, dec=BACK_DATA[:,2] * units.deg, frame="fk5")

    itemsize = MPI.INT.Get_size()
    if rank == 0:
        # bytes for 10 double elements
        nbytes = fore_num*itemsize
    else:
        nbytes = 0
    win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
    buf1, itemsize = win1.Shared_query(0)
    pair_count = numpy.ndarray(buffer=buf1, dtype=numpy.intc, shape=(fore_num,)) # array filled with zero

    stack_count = 0
    for ig in my_task:

        fore_pos = SkyCoord(ra=RA_f[ig] * units.deg, dec=DEC_f[ig] * units.deg, frame="fk5")
        sep = fore_pos.separation(backgal_pos).radian

        idx_z = Z >= Z_f[ig] + 0.1
        # the transverse distance
        dist_tran = DIST_f[ig]*coeff*sep
        idx1 = dist_tran >= radius_bin[radius_id]
        idx2 = dist_tran < radius_bin[radius_id + 1]
        idx_t = idx1 & idx2 & idx_z
        pair_num = idx_t.sum()
        pair_count[ig] = pair_num

        if pair_num > 0:
            # foregal id
            stack_data = numpy.ones((pair_num,8))*ig
            # backgal id
            stack_data[:, 1] = BACK_DATA[idx_t][:,0]
            # separation radian
            stack_data[:, 2] = sep[idx_t]
            # RA & DEC
            stack_data[:, 3] = BACK_DATA[idx_t][:,1]
            stack_data[:, 4] = BACK_DATA[idx_t][:,2]
            # Z
            stack_data[:, 5] = BACK_DATA[idx_t][:,3]
            # \Sigma_crit
            stack_data[:, 6] = BACK_DATA[idx_t][:,4]/DIST_f[ig]/(BACK_DATA[idx_t][:,4] - DIST_f[ig])/(1+Z_f[ig])
            # transverse distance
            stack_data[:, 7] = dist_tran[idx_t]

            if stack_count == 0:
                pair_data = stack_data
            else:
                pair_data = numpy.row_stack((pair_data, stack_data))
            stack_count += 1

    # collect data from each cpu
    if stack_count > 1:
        data_sp = pair_data.shape
    else:
        pair_data = numpy.zeros((1, 1))
        data_sp = pair_data.shape

    data_sps = comm.gather(data_sp, root=0)


    if rank > 0:
        comm.Send([pair_data, MPI.DOUBLE], dest=0, tag=rank)
    else:
        if data_sp[0] > 1 and data_sp[1] > 1:
            stack_pool = [pair_data]
        else:
            stack_pool = []

        for procs in range(1, cpus):
            recvs = numpy.empty(data_sps[procs], dtype=numpy.double)
            comm.Recv(recvs, source=procs, tag=procs)
            if data_sps[procs][0] > 1 and data_sps[procs][1] > 1:
                stack_pool.append(recvs)

        for stack_count, arr in enumerate(stack_pool):
            if stack_count == 0:
                recv_buffer = arr
            else:
                recv_buffer = numpy.row_stack((recv_buffer, arr))

        h5f = h5py.File("./data/result/radius_%d_py.hdf5"%radius_id,"w")
        h5f["/data"] = recv_buffer
        h5f["/pair_count"] = pair_count
        h5f.close()
        print("%d galaxies have been found in radius bin [%.4f, %.4f]"%(pair_count.sum(), radius_bin[radius_id], radius_bin[radius_id+1]))
        idx = pair_count > 0
        print(pair_count[idx])

if cmd == "compare":
    # read foregal number
    h5f = h5py.File("./data/%s_sub.hdf5"%area_id,"r")
    Z_f = h5f["/Z"].value
    h5f.close()
    fore_num = Z_f.shape[0]

    # result of Py
    h5f = h5py.File("./data/result/radius_%d_py.hdf5" % radius_id, "r")
    py_data = h5f["/data"].value
    print("Tan R: [%f, %f]"%(py_data[:,7].min(), py_data[:,7].max()))
    py_pair_count = h5f["/pair_count"].value
    h5f.close()
    print("Py %d galaxies have been found in radius bin [%.5f, %.5f]" % (py_pair_count.sum(), radius_bin[radius_id], radius_bin[radius_id + 1]))
    idx_py = py_pair_count > 0
    # print(py_pair_count[idx_py])

    # result of CPP
    h5f = h5py.File("./data/result/radius_%d.hdf5" % radius_id, "r")
    cpp_data = h5f["/pair_data"].value
    cpp_pair_count = h5f["/pair_count"].value[-fore_num:]
    h5f.close()
    print("Tan R: [%f, %f]"%(cpp_data[:,-1].min(), cpp_data[:,-1].max()))
    print("CPP %d galaxies have been found in radius bin [%.5f, %.5f]" % (cpp_pair_count.sum(), radius_bin[radius_id], radius_bin[radius_id + 1]))
    idx_cpp = cpp_pair_count > 0
    # print(cpp_pair_count[idx_cpp])

    print("Number difference: ",numpy.sum(cpp_pair_count[idx_cpp] - py_pair_count[idx_py]))

    print("Pair data difference: ", numpy.sum(cpp_data[:,:6] - py_data[:,:6]))