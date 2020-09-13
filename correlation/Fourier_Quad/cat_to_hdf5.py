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

mode = argv[2]

area_num = 4
chip_num = 32

total_path = argv[1]

if mode == "cata_name":
    if rank == 0:

        files_nm = os.listdir(total_path)
        fields = []
        all_files = []
        # this is how they name the fields
        xlbs = ["p%d" % abs(i) if i > 0 else "m%d" % abs(i) for i in range(-5, 6)]
        ylbs = ["m%d" % abs(i) if i <= 0 else "p%d" % abs(i) for i in range(-5, 6)]
        print(xlbs)
        print(ylbs)

        for j in range(area_num):
            for ylb in ylbs:
                for xlb in xlbs:
                    field_nm = "w%d%s%s" % (j + 1, xlb, ylb)
                    if field_nm in files_nm:
                        print(field_nm)
                        fields.append(field_nm + "\n")

                        all_files.append(field_nm + "\n")

                        sub_expos = []

                        # find the exposures
                        sub_files = os.listdir(total_path + "/%s/result/"%field_nm)
                        for nm in sub_files:
                            if "_all.cat" in nm:
                                sub_expos.append(int(nm.split("p")[0]))
                        # sort the exposure labels from small to large
                        sub_expos = sorted(sub_expos)
                        # check the existence of chips
                        for expo in sub_expos:
                            for i in range(1, chip_num + 1):
                                chip_nm = "%dp_%d_shear.dat" % (expo, i)
                                if os.path.exists(total_path + "/%s/result/%s" % (field_nm, chip_nm)):
                                    all_files.append(chip_nm + "\n")
                                    print(chip_nm)

        with open(total_path + "/nname_all.dat", "w") as f:
            f.writelines(all_files)
        with open(total_path + "/nname.dat", "w") as f:
            f.writelines(fields)
        print(len(fields))

elif mode == "hdf5_cata":
    # convert the .dat to .hdf5

    fields, field_name = tool_box.field_dict(total_path + "/nname_all.dat")
    if rank == 0:
        print(len(field_name))

    field_name_sub = tool_box.alloc(field_name, cpus,"seq")[rank]

    field_avail_sub = []
    field_all_avail_sub = []

    for fns in field_name_sub:
        # read the the field data
        field_src_path = total_path + "/%s/result"%fns
        buffer = []
        try:
            fdat = numpy.loadtxt(field_src_path + "/%s.cat"%fns, dtype=numpy.float32)

            field_dst_path = field_src_path + "/%s.hdf5" % fns

            h5f = h5py.File(field_dst_path, "w")
            h5f["/field"] = fdat
            h5f["/field"].attrs["gal_num"] = fdat.shape[0]

            expos = list(fields[fns].keys())
            expos_num = len(expos)

            # read the exposure files
            expo_label = 0
            for exp_nm in expos:
                try:
                    edat = numpy.loadtxt(field_src_path + "/%s_all.cat" % exp_nm, dtype=numpy.float32)
                    h5f["/expo_%d"%expo_label] = edat
                    h5f["/expo_%d"%expo_label].attrs["exposure_name"] = exp_nm
                    h5f["/expo_%d"%expo_label].attrs["gal_num"] = edat.shape[0]
                    expo_label += 1
                    h5f_expo = h5py.File(field_src_path + "/%s_all.hdf5" % exp_nm,"w")
                    h5f_expo["/data"] = edat
                    h5f_expo.close()
                    buffer.append("%s_all.cat\n" % exp_nm)
                except:
                    print("%d %s-%s empty!" % (rank, fns, exp_nm))

            # how many exposures in this field
            h5f["/expos_num"] = numpy.array([expo_label], dtype=numpy.intc)

            h5f.close()

            field_avail_sub.append(fns+"\n")
            field_all_avail_sub.append(fns + "\n")
            if len(buffer)> 0:
                field_all_avail_sub.extend(buffer)

        except:
            print("%d %s/%s.cat empty!"%(rank, field_src_path, fns),os.path.exists(field_src_path + "/%s.cat"%fns))

    field_collection = comm.gather(field_avail_sub, root=0)
    field_all_collection = comm.gather(field_all_avail_sub, root=0)
    if rank == 0:
        field_avail = []
        for fsb in field_collection:
            field_avail.extend(fsb)
        print(len(field_avail))
        with open(total_path + "/nname_avail.dat","w") as f:
            f.writelines(field_avail)

        field_all_avail = []
        for fsb in field_all_collection:
            field_all_avail.extend(fsb)
        print(len(field_all_avail))
        with open(total_path + "/nname_all_avail.dat","w") as f:
            f.writelines(field_all_avail)

else:
    # collection
    result_path = total_path + "/total.hdf5"

    fields, field_name = tool_box.field_dict(total_path + "/nname_avail.dat")
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

                h5f = h5py.File(total_path + "/%s/result/%s.hdf5" %(fns, fns),"r")
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