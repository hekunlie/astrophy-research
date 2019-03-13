import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import numpy
from mpi4py import MPI
import time
import tool_box
import warnings
import h5py


warnings.filterwarnings("error")

# stack the CFHT shear catalog files into one big .hdf5 file.

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()

envs_path = "%s/work/envs/envs.dat" % my_home
get_contents = [['cfht', "cfht_path_catalog", '1'], ['cfht', "cfht_path_result", '1']]
path_items = tool_box.config(envs_path, ['get', 'get'], get_contents)
total_cata_path, result_path = path_items

field_dict, fields = tool_box.field_dict(total_cata_path + "cfht_cata/nname.dat")

area_num = 4
if rank == 0:
    h5f = h5py.File(total_cata_path + "cfht_cata/cata.hdf5", "w")
    h5f.close()
    gal_count = 0

num_in_field = []

for area_id in range(1, 1+area_num):
    field_paths = []
    field_count = 0
    for field in fields:
        if "w%d"%area_id in field:
            field_path = total_cata_path + "cfht_cata/%s.dat" %field
            field_paths.append(field_path)
            field_count += 1
    sub_fields = tool_box.allot(field_paths, cpus)[rank]

    fcount = 0
    data = None
    for field_name in sub_fields:
            t11 = time.time()
            try:
                temp = numpy.loadtxt(field_name)
                if fcount == 0:
                    data = temp.copy()
                else:
                    data = numpy.row_stack((data, temp))
                fcount += 1
            except:
                print(rank, " can't find ", field_name)
            t12 = time.time()

    if data is not None:
        data_sp = data.shape
    else:
        data = numpy.array([[1]])
        data_sp = data.shape
    num_in_field.append(data_sp)

    data_sps = comm.gather(data_sp, root=0)

    if rank > 0:
        comm.Send(data, dest=0, tag=rank)
    else:
        if data is not None:
            data_pool = [data]
        else:
            data_pool = []

        for procs in range(1, cpus):
            recvs = numpy.empty(data_sps[procs], dtype=numpy.float64)
            comm.Recv(recvs, source=procs, tag=procs)
            if data_sps[procs][0] > 1 and data_sps[procs][1] > 1:
                data_pool.append(recvs)

        for i in range(len(data_pool)):
            if i == 0:
                data = data_pool[i]
            else:
                data = numpy.row_stack((data, data_pool[i]))

        h5f = h5py.File(total_cata_path + "cfht_cata/cata.hdf5", "r+")
        h5f['/w_%d'%area_id] = data
        h5f.close()
        gal_count += len(data)
        print("Totally, %d galaxies are detected in W_%d"%(len(data), area_id))
print(rank, num_in_field)
te = time.clock()
if rank == 0:
    print(te-ts, gal_count)