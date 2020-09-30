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
chip_num = 36

total_path = argv[1]

if mode == "cata_name":
    if rank == 0:

        with open("anomaly_expo.dat", "r") as f:
            contents = f.readlines()
        anomaly_expo = []
        anomaly_val = []
        for cc in contents:
            val, field_nm, expo = cc.rstrip().split()[1:4]
            anomaly_expo.append(int(expo.split("p")[0]))
            anomaly_val.append(float(val))
            # print(anomaly_expo[-1],anomaly_val[-1])

        files_nm = os.listdir(total_path)
        fields = []
        all_files = []
        # # this is how they name the fields
        # xlbs = ["p%d" % abs(i) if i > 0 else "m%d" % abs(i) for i in range(-5, 6)]
        # ylbs = ["m%d" % abs(i) if i <= 0 else "p%d" % abs(i) for i in range(-5, 6)]
        # print(xlbs)
        # print(ylbs)

        for j in range(area_num):
            for fnm in files_nm:
                if "w%d"%(j+1) in fnm:
                    # print(field_nm)

                    sub_expos = []
                    # find the exposures
                    sub_files = os.listdir(total_path + "/%s/result/" % fnm)
                    for nm in sub_files:
                        if "_shear" in nm:
                            # check the anomaly of each exposure
                            expo = int(nm.split("p")[0])
                            if expo not in sub_expos:
                                sub_expos.append(expo)

                    if len(sub_expos) > 0:
                        fields.append(fnm + "\n")
                        all_files.append(fnm + "\n")

                    # sort the exposure labels from small to large
                    sub_expos = sorted(sub_expos)
                    # check the existence of chips
                    for expo in sub_expos:
                        for i in range(1, chip_num + 1):
                            chip_nm = "%dp_%d_shear.dat" % (expo, i)
                            chip_path = total_path + "/%s/result/%s" % (fnm, chip_nm)
                            if os.path.exists(chip_path):
                                all_files.append(chip_nm + "\n")
                                # print(chip_nm)
                            else:
                                print("Can't find %s" % chip_path)
            # for ylb in ylbs:
            #     for xlb in xlbs:
            #         field_nm = "w%d%s%s" % (j + 1, xlb, ylb)
            #         if field_nm in files_nm:
            #             # print(field_nm)
            #             fields.append(field_nm + "\n")
            #
            #             all_files.append(field_nm + "\n")
            #
            #             sub_expos = []
            #
            #             # find the exposures
            #             sub_files = os.listdir(total_path + "/%s/result/"%field_nm)
            #             for nm in sub_files:
            #                 if "_all.cat" in nm:
            #                     # check the anomaly of each exposure
            #                     expo = int(nm.split("p")[0])
            #                     expo_lb = anomaly_expo.index(expo)
            #                     if anomaly_val[expo_lb] > 0.01 or anomaly_val[expo_lb] < 0:
            #                         print("Anomaly %d %.6f"%(expo,anomaly_val[expo_lb]))
            #                         continue
            #                     else:
            #                         sub_expos.append(expo)
            #
            #             # sort the exposure labels from small to large
            #             sub_expos = sorted(sub_expos)
            #             # check the existence of chips
            #             for expo in sub_expos:
            #                 for i in range(1, chip_num + 1):
            #                     chip_nm = "%dp_%d_shear.dat" % (expo, i)
            #                     chip_path = total_path + "/%s/result/%s" % (field_nm, chip_nm)
            #                     if os.path.exists(chip_path):
            #                         all_files.append(chip_nm + "\n")
            #                         # print(chip_nm)
            #                     else:
            #                         print("Can't find %s"%chip_path)

        with open(total_path + "/nname_field_chips.dat", "w") as f:
            f.writelines(all_files)
        with open(total_path + "/nname_field.dat", "w") as f:
            f.writelines(fields)
        print(len(fields))

elif mode == "hdf5_cata":
    # convert the .dat to .hdf5

    fields, field_name = tool_box.field_dict(total_path + "/nname_field_chips.dat")
    if rank == 0:
        print(len(field_name))

    field_name_sub = tool_box.alloc(field_name, cpus, "seq")[rank]

    fields_sub_avail_sub = []
    fields_sub_raw_avail_sub = []

    for fns in field_name_sub:
        # read the the field data
        field_src_path = total_path + "/%s/result"%fns

        try:
            expos = list(fields[fns].keys())
            expos_num = len(expos)

            # read the exposure files
            expo_label = 0

            for exp_nm in expos:
                expo_src_path = field_src_path + "/%s_all.cat" % exp_nm
                expo_h5_path = field_src_path + "/%s_all.hdf5" % exp_nm
                try:
                    edat = numpy.loadtxt(expo_src_path, dtype=numpy.float32)

                    # Nan check
                    idx = numpy.isnan(edat)
                    if idx.sum() > 0:
                        print("Find Nan in ",expo_src_path)
                    expo_label += 1

                    h5f_expo = h5py.File(expo_h5_path,"w")
                    h5f_expo["/data"] = edat
                    h5f_expo.close()
                    fields_sub_avail_sub.append(expo_h5_path+"\n")
                except:
                    if os.path.exists(expo_src_path):
                        print("%d Failed in reading %s %d Bytes !" % (rank, expo_src_path, os.path.getsize(expo_src_path)))
                    else:
                        print("%d can't find %s!" % (rank, expo_src_path))

                # read & stack the chip data
                chip_label = 0
                for iexp in range(1, chip_num+1):
                    chip_src_path = field_src_path + "/%s_%d_shear.dat" %(exp_nm, iexp)
                    chip_h5_path = field_src_path + "/%s_%d_shear.hdf5" %(exp_nm, iexp)
                    try:
                        chip_data_raw = numpy.loadtxt(chip_src_path, dtype=numpy.float32, skiprows=1)
                        h5f_chip = h5py.File(chip_h5_path,"w")
                        h5f_chip["/data"] = chip_data_raw
                        h5f_chip.close()

                        row, col = chip_data_raw.shape
                        # Nan check
                        idx = numpy.isnan(chip_data_raw)
                        if idx.sum() > 0:
                            print("Find Nan in ", chip_src_path)

                        chip_data_raw_ic = numpy.zeros((row, col+1))
                        chip_data_raw_ic[:,0] = iexp - 1
                        chip_data_raw_ic[:,1:] = chip_data_raw

                        if chip_label == 0:
                            stack_chip_data_raw = chip_data_raw_ic
                        else:
                            stack_chip_data_raw = numpy.row_stack((stack_chip_data_raw, chip_data_raw_ic))
                        chip_label += 1

                    except:
                        if os.path.exists(chip_src_path):
                            print("Failed in read chip %s %d Bytes"%(chip_src_path, os.path.getsize(chip_src_path)))
                        else:
                            print("Failed in read chip %s 0 Bytes"%(chip_src_path))

                if chip_label > 0:
                    h5f_chip_stack = h5py.File(field_src_path + "/%s_all_raw.hdf5" % exp_nm, "w")
                    h5f_chip_stack["/data"] = stack_chip_data_raw
                    h5f_chip_stack.close()
                    fields_sub_raw_avail_sub.append(field_src_path + "/%s_all_raw.hdf5\n" % exp_nm)

        except:
            src_cat_path = field_src_path + "/%s.cat"%fns
            if os.path.exists(src_cat_path):
                print("%d %s empty! %d Bytes"%(rank, src_cat_path,os.path.getsize(src_cat_path)), os.path.exists(src_cat_path))
            else:
                print("%d %s empty! 0 Bytes" % (rank, src_cat_path), os.path.exists(src_cat_path))

    field_collection = comm.gather(fields_sub_avail_sub, root=0)
    field_raw_collection = comm.gather(fields_sub_raw_avail_sub, root=0)
    comm.Barrier()
    if rank == 0:
        field_avail = []
        for fsb in field_collection:
            field_avail.extend(fsb)
        print("Totally: ",len(field_avail), " exposures")
        with open(total_path + "/nname_field_expo_avail.dat","w") as f:
            f.writelines(field_avail)

        field_raw_avail = []
        for fsb in field_raw_collection:
            field_raw_avail.extend(fsb)
        print("Totally: ",len(field_raw_avail), " raw exposures")
        with open(total_path + "/nname_field_raw_expo_avail.dat","w") as f:
            f.writelines(field_raw_avail)


else:
    # collection

    source_list_nm = argv[2]
    result_nm = argv[3]

    expos = []
    with open(total_path + "/"+source_list_nm, "r") as f:
        conts = f.readlines()
    for nm in conts:
        expos.append(nm.split("\n")[0])

    if rank == 0:
        print(len(expos)," exposures")

    my_sub_area_list = tool_box.alloc(expos,cpus)[rank]
    # print(rank, i, len(my_sub_area_list))

    if len(my_sub_area_list) > 0:
        for tag, expo_path in enumerate(my_sub_area_list):

            h5f = h5py.File(expo_path,"r")
            temp = h5f["/data"][()]

            # Nan check
            idx = numpy.isnan(temp)
            if idx.sum() > 0:
                print("Find Nan ", expo_path)
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

        # Nan check
        idx = numpy.isnan(data)
        if idx.sum() > 0:
            print("Find Nan in final data")

        idx_v = numpy.invert(idx)
        h5f = h5py.File(total_path + "/%s"%result_nm, "w")
        h5f["/data"] = data
        h5f.close()



    comm.Barrier()