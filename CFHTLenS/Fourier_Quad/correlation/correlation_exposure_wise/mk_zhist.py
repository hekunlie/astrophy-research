import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import correlation_function_tool as cf_tool
from plot_tool import Image_Plot
import h5py
import numpy

argc = len(argv)

src_data = argv[1]
dst_path = argv[2]

# z bins not used in the MCMC
throw_bins = []
for i in range(3, argc):
    throw_bins.append(int(argv[i]))
print("Throw: ",throw_bins)

z_col, ze_col = 8,9

h5f = h5py.File(src_data,"r")
data = h5f["/data"][()]
h5f.close()

redshift_bin = numpy.array([0.2, 0.39, 0.58, 0.72, 0.86, 1.02, 1.3],dtype=numpy.float32)
bin_num = len(redshift_bin) - 1
redshift_bin_resort = numpy.zeros((bin_num,), dtype=numpy.intc)

nz_bin_num = 340
zehist, zebin, zebin_cent = cf_tool.get_nz(data[:, z_col], redshift_bin, data[:, ze_col], nz_bin_num, (0,3), False)

zehist_resort = numpy.zeros((bin_num-len(throw_bins),nz_bin_num-1))
tag = 0

img = Image_Plot(xpad=0.2,ypad=0.2)
img.subplots(1,1)

for i in range(bin_num):
    alpha = 0.3
    if i not in throw_bins:
        zehist_resort[tag] = zehist[i]
        redshift_bin_resort[i] = 1
        tag += 1
        alpha = 1
    img.axs[0][0].plot(zebin_cent, zehist[i], alpha=alpha,label="[%.2f~%.2f]" % (redshift_bin[i], redshift_bin[i + 1]))
img.save_img("%s/zhist.png"%dst_path)

print(tag," bins used")
print("used bins: ",redshift_bin_resort)


h5f = h5py.File("%s/zhist.hdf5"%dst_path,"w")
h5f["/zhist_ori"] = zehist
h5f["/zhist"] = zehist_resort
h5f["/uesd_bin_labels"] = redshift_bin_resort
h5f["/zbin"] = zebin
h5f["/zbin_cent"] = zebin_cent
h5f.close()