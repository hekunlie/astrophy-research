from sys import path,argv
path.append("/home/hklee/work/mylib")
import h5py
import numpy
import os


parent_path = argv[1]
# file_num = int(argv[2])
data_len_i = int(argv[2])

names = ["sheared","non_sheared"]
file_num = [50, 10]
data_type = ["noise_free.hdf5", "noisy_cpp.hdf5"]

for nm_tag, nm in enumerate(names):

    total_data_len = int(file_num[nm_tag] * data_len_i)
    stack_data = numpy.zeros((total_data_len, 5), dtype=numpy.float32)

    for dt in data_type:
        stack_tag = 0
        for i in range(file_num[nm_tag]):
            file_path_nf = parent_path + "/data/" + "%s_data_%d_%s"%(nm, i, dt)
            file_path_ny = parent_path + "/data/" + "%s_data_%d_%s"%(nm, i, dt)

            if not os.path.exists(file_path_nf):
                print("No ", "%s_%d_%s"%(nm, i, dt))
                break

            st, ed = int(i * data_len_i), int((i + 1) * data_len_i)

            h5f = h5py.File(file_path_nf, "r")
            stack_data[st:ed] = h5f["/data"][()]
            h5f.close()
            stack_tag += 1

        if stack_tag > 0:
            stack_path = parent_path + "/data/stack_%s_%s"%(nm, dt)
            h5f = h5py.File(stack_path, "w")
            h5f["/data"] = stack_data
            h5f.close()
            print("Stack data ",stack_path)




for nm_tag, nm in enumerate(names):

    total_data_len = int(file_num[nm_tag] * data_len_i)

    src_ra = numpy.zeros((total_data_len,), dtype=numpy.float32)
    src_dec = numpy.zeros((total_data_len,), dtype=numpy.float32)
    src_z = numpy.zeros((total_data_len,), dtype=numpy.float32)
    src_g1 = numpy.zeros((total_data_len,), dtype=numpy.float32)
    src_g2 = numpy.zeros((total_data_len,), dtype=numpy.float32)
    stack_tag = 0
    for i in range(file_num[nm_tag]):
        param_path = parent_path + "/params/%s_para_%d.hdf5"% (nm, i)
        if not os.path.exists(param_path):
            print("No ", "%s_para_%d.hdf5"% (nm, i))
            break

        st, ed = int(i * data_len_i), int((i + 1) * data_len_i)

        h5f = h5py.File(param_path, "r")
        src_ra[st:ed] = h5f["/ra"][()] / 3600
        src_dec[st:ed] = h5f["/dec"][()] / 3600
        src_z[st:ed] = h5f["/z"][()]
        src_g1[st:ed] = h5f["/gamma2"][()]
        src_g2[st:ed] = h5f["/gamma1"][()]
        h5f.close()
        stack_tag += 1

    if stack_tag > 0:
        stack_param_path = parent_path + "/params/stack_%s_para.hdf5"%nm
        h5f = h5py.File(stack_param_path, "w")
        h5f["/ra"] = src_ra
        h5f["/dec"] = src_dec
        h5f["/z"] = src_z
        h5f["/g1"] = src_g1
        h5f["/g2"] = src_g2
        h5f.close()
        print("Stack param, ", stack_param_path)