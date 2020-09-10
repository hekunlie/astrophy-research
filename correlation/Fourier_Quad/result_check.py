import h5py

data_path = "D:"
mg_bin_num = 10
zbin_num = 6
theta_bin_num = 5
chi_guess_num = 50

h5f = h5py.File(data_path + "/w1m0m0-819724p_num_count.hdf5","r")
gtt_1 =h5f["/t"][()]
gxx_1 =h5f["/t"][()]
theta_1 = h5f["/theta"][()]
h5f.close()

h5f = h5py.File(data_path + "/w1m0m0-819724p_num_count_zstack.hdf5","r")
gtt_2 =h5f["/tt"][()]
gxx_2 =h5f["/xx"][()]
theta_2 = h5f["/theta"][()]
h5f.close()
print(gtt_1.shape, gtt_2.shape,theta_1.shape, theta_2.shape)

t = 0
for i in range(zbin_num):
    for j in range(i, zbin_num):
        ij = int(i*zbin_num + j)
        ji = int(j*zbin_num + i)
        print("i, j",i, j, ij)
        print("j, i",j, i, ji)
        print(t)

        theta_block_len = mg_bin_num*chi_guess_num*theta_bin_num
        print(theta_block_len)
        a = gtt_1[int(ij*theta_block_len):int((ij+1)*theta_block_len)]
        b = gtt_1[int(ji*theta_block_len):int((ji+1)*theta_block_len)]
        c = gtt_2[int(t*theta_block_len):int((t+1)*theta_block_len)]
        print("Diff: ")
        if i==j:
            diff = c - b
            print(diff.min(), diff.max())
            diff = c - a
            print(diff.min(), diff.max())
        else:
            diff = c - b - a
            print(diff.min(), diff.max())

        print("Diff: ")
        a = theta_1[int(ij)]
        b = theta_1[int(ji)]
        c = theta_2[int(t)]
        if i == j:
            diff = c - b
            print(diff.min(), diff.max())
            diff = c - a
            print(diff.min(), diff.max())
        else:
            diff = c - b - a
            print(diff.min(), diff.max())
        t += 1
        print("\n")