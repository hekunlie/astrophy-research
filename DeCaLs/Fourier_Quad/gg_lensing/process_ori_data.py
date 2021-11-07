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
import hk_c4py


def find_overlap(ra, dec):
    npz = numpy.load("/home/hklee/work/catalog/DECALS_cat_hist.npz")
    total_hist = npz["arr_0"]
    # print(total_hist.shape)
    ra_bin, dec_bin = npz["arr_1"], npz["arr_2"]
    # print(ra_bin.max(), dec_bin.max())

    ra_bin_num, dec_bin_num = ra_bin.shape[0] - 1, dec_bin.shape[0] - 1
    total_hist_num = int(total_hist.shape[0] / dec_bin_num)
    # print(ra_bin_num, dec_bin_num)

    decals_mask = numpy.zeros((dec_bin_num, ra_bin_num))
    for i in range(total_hist_num):
        st, ed = int(i * dec_bin_num), int((i + 1) * dec_bin_num)
        decals_mask += total_hist[st:ed]

    labels, idx = hk_c4py.find_overlap_mask(ra, dec, decals_mask, ra_bin, dec_bin)
    return labels, idx

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

    for tag, ff in enumerate(["gal_jkf"]):

        data = fits.open(cata_path + "/%s/%s.fits" % (fnm, ff))[1].data
        col_nm = data.dtype
        src_num = data.shape[0]

        if tag == 0:
            src_num_0 = src_num
            zp = data["Z_PHOT_MEDIAN"]

        sub_data = numpy.zeros((src_num_0, 5), dtype=numpy.float32)
        sub_data[:, 1] = data["RA"][:src_num_0]
        sub_data[:, 2] = data["DEC"][:src_num_0]
        sub_data[:, 3] = zp[:src_num_0]

        idx = find_overlap(sub_data[:, 1],sub_data[:, 2])[1]

        print(fnm, src_num, idx.sum(), idx.sum()/src_num)

        data_overlap = sub_data[idx]
        pix_ra_dec, pixel_count, src_pix_label = hk_healpy_tool.get_healpy_pix(data_overlap[:, 1],data_overlap[:, 2], NSIDE)[:3]
        #
        h5f = h5py.File(cata_path + "/%s/%s.hdf5" % (fnm, ff), "w")
        h5f["/data"] = data_overlap
        h5f["/data_ori"] = sub_data
        h5f["/pix_ra_dec"] = pix_ra_dec
        h5f["/pix_count"] = pixel_count
        h5f["/data_pix_label"] = src_pix_label
        h5f.close()

        # h5f = h5py.File(cata_path + "/%s/%s.hdf5" % (fnm, ff), "r")
        # pix_ra_dec = h5f["/pix_ra_dec"][()]
        # pixel_count = h5f["/pix_count"][()]
        # src_pix_label = h5f["/data_pix_label"][()]
        # sub_data = h5f["/data"][()]
        # h5f.close()

        rs = numpy.random.randint(1, 100000)
        eff_pix_num = len(pixel_count)
        for cent_num in [120, 200]:
            if src_num <= int(eff_pix_num*2):
                group_label = KMeans(n_clusters=cent_num, random_state=rs).fit_predict(sub_data[:,[1,2]])
            else:
                group_label, group_label_pixel = hk_healpy_tool.kmeans_pix(pix_ra_dec, pixel_count, src_pix_label, eff_pix_num, cent_num, rs)

            h5f = h5py.File(cata_path + "/%s/%s.hdf5" % (fnm, ff), "r+")
            try:
                del h5f["/jkf_label_%d"%cent_num]
                print("Find /jkf_label")
            except:
                pass
            h5f["/jkf_label_%d"%cent_num] = group_label
            h5f.close()


comm.Barrier()