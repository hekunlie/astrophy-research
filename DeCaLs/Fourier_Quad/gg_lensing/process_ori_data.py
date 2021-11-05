import os
my_home = os.popen("echo $HK_MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import shutil
import h5py
import numpy
from mpi4py import MPI
import hk_tool_box
from astropy.io import fits
import hk_healpy_tool
from sklearn.cluster import KMeans



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

deg2rad = numpy.pi/180

cata_path = argv[1]
with open(cata_path + "/file_list", "r") as f:
    ff = f.readlines()
ff_sub = hk_tool_box.alloc(ff, cpus)[rank]

NSIDE = 256

cent_num = 200
for fnm in ff_sub:
    fnm = fnm.split("\n")[0]

    for tag, ff in enumerate(["gal_jkf","rand_jkf"]):

        data = fits.open(cata_path + "/%s/%s.fits" % (fnm, ff))[1].data
        col_nm = data.dtype
        src_num = data.shape[0]
        print(fnm, src_num)

        # if tag == 0:
        #     src_num_0 = src_num
        #     zp = data["Z_PHOT_MEDIAN"]
        #
        # sub_data = numpy.zeros((src_num_0, 5), dtype=numpy.float32)
        # sub_data[:, 1] = data["RA"][:src_num_0]
        # sub_data[:, 2] = data["DEC"][:src_num_0]
        # sub_data[:, 3] = zp[:src_num_0]

        # pix_ra_dec, pixel_count, src_pix_label = hk_healpy_tool.get_healpy_pix(sub_data[:, 1],sub_data[:, 2], NSIDE)[:3]
        #
        # h5f = h5py.File(cata_path + "/%s/%s.hdf5" % (fnm, ff), "w")
        # h5f["/data"] = sub_data
        # h5f["/pix_ra_dec"] = pix_ra_dec
        # h5f["/pix_count"] = pixel_count
        # h5f["/data_pix_label"] = src_pix_label
        # h5f.close()

        h5f = h5py.File(cata_path + "/%s/%s.hdf5" % (fnm, ff), "r")
        pix_ra_dec = h5f["/pix_ra_dec"][()]
        pixel_count = h5f["/pix_count"][()]
        src_pix_label = h5f["/data_pix_label"][()]
        sub_data = h5f["/data"][()]
        h5f.close()

        rs = numpy.random.randint(1, 100000)
        eff_pix_num = len(pixel_count)
        if src_num <= int(eff_pix_num*2):
            group_label = KMeans(n_clusters=cent_num, random_state=rs).fit_predict(sub_data[:,[1,2]])
        else:
            group_label, group_label_pixel = hk_healpy_tool.kmeans_pix(pix_ra_dec, pixel_count, src_pix_label, eff_pix_num, cent_num, rs)

        h5f = h5py.File(cata_path + "/%s/%s.hdf5" % (fnm, ff), "r+")
        try:
            del h5f["/jkf_label"]
            print("Find /jkf_label")
        except:
            pass
        h5f["/jkf_label"] = group_label
        h5f.close()


comm.Barrier()