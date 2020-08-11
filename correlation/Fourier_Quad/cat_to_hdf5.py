import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
from mpi4py import MPI
import tool_box
import warnings

warnings.filterwarnings('error')


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

mode = argv[1]

total_path = "/mnt/ddnfs/data_users/hkli/CFHT/catalog/fourier_cata"

area_num = 4

if mode == "prepare":
    if rank == 0:
        original_path = total_path + "/original_cata"
        files = os.listdir(original_path)
        for i in range(1,1+area_num):
            for fn in files:
                if "w%d"%i in fn:
                    pass


    comm.Barrier()
    fields, field_name = tool_box.field_dict("nname_all.dat")
    if rank == 0:
        print(len(field_name))

    field_name_sub = tool_box.alloc(field_name, cpus)[rank]

    field_avail_sub = []

    for fns in field_name_sub:

        # read the the field data
        field_src_path = total_path + "/original_cata/%s"%fns

        try:
            fdat = numpy.loadtxt(field_src_path + "/result/%s.cat"%fns,dtype=numpy.float32)

            field_dst_path = total_path + "/hdf5_cata/%s" % fns
            if not os.path.exists(field_dst_path):
                os.makedirs(field_dst_path)

            h5f = h5py.File(field_dst_path + "/%s.hdf5"%fns, "w")
            h5f["/field"] = fdat
            h5f["/field"].attrs["gal_num"] = fdat.shape[0]

            expos = list(fields[fns].keys())
            expos_num = len(expos)

            # read the exposure files
            expo_label = 0
            for exp_nm in expos:
                try:
                    edat = numpy.loadtxt(field_src_path + "/result/%s_all.cat" % exp_nm, dtype=numpy.float32)
                    h5f["/expo_%d"%expo_label] = edat
                    h5f["/expo_%d"%expo_label].attrs["exposure_name"] = exp_nm
                    h5f["/expo_%d"%expo_label].attrs["gal_num"] = edat.shape[0]
                    expo_label += 1
                except:
                    print("%d %s-%s empty!" % (rank, fns, exp_nm))

            # how many exposures in this field
            h5f["/expos_num"] = numpy.array([expo_label], dtype=numpy.intc)

            h5f.close()

            field_avail_sub.append(fns+"\n")

        except:
            print("%d %s empty!"%(rank, fns))

    field_collection = comm.gather(field_avail_sub, root=0)
    if rank == 0:
        field_avail = []
        for fsb in field_collection:
            field_avail.extend(fsb)
        print(len(field_avail))
        with open(total_path + "/hdf5_cata/nname_avail.dat","w") as f:
            f.writelines(field_avail)

else:
    # collection
    result_path = total_path + "/hdf5_cata/total.hdf5"

    fields, field_name = tool_box.field_dict(total_path + "/hdf5_cata/nname_avail.dat")
    if rank == 0:
        print(len(field_name))
        h5f = h5py.File(result_path, "w")
        h5f.close()

    area_num = 4

    for i in range(area_num):
        sub_area_list = []
        for fns in field_name:
            if "w%d" % (i + 1) in fns:
                sub_area_list.append(fns)
        my_sub_area_list = tool_box.alloc(sub_area_list,cpus)[rank]
        # print(rank, i, len(my_sub_area_list))

        if len(my_sub_area_list) > 0:
            for tag, fns in enumerate(my_sub_area_list):

                h5f = h5py.File(total_path + "/hdf5_cata/%s/%s.hdf5" %(fns, fns),"r")
                temp = h5f["/field"][()]

                if tag == 0:
                    data = temp
                else:
                    data = numpy.row_stack((data, temp))
                h5f.close()

            sp = data.shape
        else:
            sp = (0,0)

        sp_total = comm.gather(sp, root=0)
        # print(i,rank, data.shape, data.dtype, sp, data[0,:5])
        comm.Barrier()

        if rank > 0 and sp[0] > 0:
            comm.Send([data,MPI.FLOAT], dest=0, tag=rank)
        else:
            for ir in range(1, cpus):
                if sp_total[ir][0] > 0:
                    recv_buf = numpy.empty(sp_total[ir],dtype=numpy.float32)
                    comm.Recv(recv_buf,source=ir, tag=ir)
                    data = numpy.row_stack((data, recv_buf))

            h5f = h5py.File(result_path, "r+")
            h5f["/w%d"%i] = data
            h5f.close()

        comm.Barrier()