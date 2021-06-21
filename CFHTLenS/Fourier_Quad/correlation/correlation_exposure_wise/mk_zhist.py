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


z_col, ze_col = 8,9
for i in range(4):
    h5f = h5py.File(src_data + "/stack_data_%d.hdf5"%i,"r")
    temp = h5f["/data"][()]
    h5f.close()
    if i == 0:
        data = temp
    else:
        data = numpy.row_stack((data, temp))

redshift_bin = numpy.array([0.2, 0.39, 0.58, 0.72, 0.86, 1.02, 1.3], dtype=numpy.float32)
# redshift_bin = numpy.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4],dtype=numpy.float32)
bin_num = len(redshift_bin) - 1


nz_bin_num = 200
zehist, zebin, zebin_cent = cf_tool.get_nz(data[:, z_col], redshift_bin, data[:, ze_col], nz_bin_num, (0,3), False)


tag = 0

img = Image_Plot(xpad=0.2,ypad=0.2)
img.subplots(1,1)

for i in range(bin_num):
    alpha = 1
    img.axs[0][0].plot(zebin_cent, zehist[i], alpha=alpha,label="[%.2f~%.2f]" % (redshift_bin[i], redshift_bin[i + 1]))
img.axs[0][0].legend()
img.save_img("%s/zhist.png"%dst_path)


h5f = h5py.File("%s/zhist.hdf5"%dst_path,"w")
h5f["/zhist"] = zehist
h5f["/zbin"] = zebin
h5f["/zbin_cent"] = zebin_cent
h5f.close()