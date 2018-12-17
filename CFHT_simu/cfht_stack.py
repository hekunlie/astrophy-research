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

field_dict, fields = tool_box.field_dict(total_cata_path + "nname.dat")

# filter_path = result_path + "filtered.dat"
# filter_exist = os.path.exists(filter_path)
#
# if filter_exist:
#     if rank == 0:
#         print("Stack filtered data")
#     with open(filter_path, "r") as f:
#         fields = f.readlines()
# else:
#     if rank == 0:
#         print("Stack unfiltered data")
#
# area_paths = []
# field_count = 0
# for area in fields:
#     area = area.split("\n")[0]
#     if filter_exist:
#         area_path = result_path + "field/%s/%s_shear.dat.npz" %(area,area)
#         field_count += 1
#     else:
#         area_path = total_path + area + "/result/%s_shear.dat"%area
#         field_count += 1
#     area_paths.append(area_path)

area_paths = []
field_count = 0
for field in fields:
    field_path = total_cata_path + "%s/result/%s_shear.dat" %(field, field)
    field_count += 1
sub_fields = tool_box.allot(fields, cpus)[rank]

fcount = 0
for field_name in sub_fields:
        sub_filed_path = total_cata_path + "%s/result/%s_shear.dat" %(field_name, field_name)
        t11 = time.time()
        try:
            temp = numpy.loadtxt(sub_filed_path)
            if fcount == 0:
                data = temp.copy()
            else:
                data = numpy.row_stack((data, temp))
            fcount += 1
        except:
            print(rank," can't find ", field_name)
        t12 = time.time()
        if rank == 0:
            print(t12- t11)

recv_fcount = comm.gather(fcount, root=0)
sp = data.shape
recv_sp = comm.gather(sp, root=0)

if rank > 0:
    comm.Send(data, dest=0, tag=rank)
else:
    for procs in range(1, cpus):
        recvs = numpy.empty(recv_sp[procs], dtype=numpy.float64)
        comm.Recv(recvs, source=procs, tag=procs)
        data = numpy.row_stack((data, recvs))
    fh5 = h5py.File(result_path+"cata.hdf5","w")
    fh5['/data'] = data
    fh5.close()
    print("Totally, %d galaxies are detected in %d (%d) fields"%(len(data), sum(recv_fcount), field_count))

te = time.clock()
if rank == 0:
    print(te-ts)