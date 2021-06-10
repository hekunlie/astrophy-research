import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append("%s/work/mylib/"% my_home)
import h5py
import numpy
from plot_tool import Image_Plot
from mpi4py import MPI
import tool_box
import warnings

warnings.filterwarnings('error')


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

mode = argv[1]

ori_cat_chara = "_all_alt.cat"

data_band = ["g", "r", "z"]

if mode == "cata_name":

    total_path = argv[2]
    ori_cat_path = argv[3]

    if rank == 0:

        # with open("anomaly_expo.dat", "r") as f:
        #     contents = f.readlines()
        # anomaly_expo = []
        # anomaly_val = []
        # for cc in contents:
        #     val, field_nm, expo = cc.rstrip().split()[1:4]
        #     anomaly_expo.append(int(expo.split("p")[0]))
        #     anomaly_val.append(float(val))
        #     # print(anomaly_expo[-1],anomaly_val[-1])
        all_files = []
        for band in data_band:
            files_nm = os.listdir(ori_cat_path + "/%s"%(band))
            band_files = []

            for fnm in files_nm:
                if ori_cat_chara in fnm:
                    band_files.append("%s\n"%fnm.split(".")[0])
                    all_files.append("%s %s\n"%(fnm.split(".")[0],band))

            with open(total_path + "/cat_inform/exposure_name_%s_band.dat"%band, "w") as f:
                f.writelines(band_files)

            print(len(band_files), " exposures in %s band\n%s/%s"%(band,ori_cat_path,band))

        with open(total_path + "/cat_inform/exposure_name_all_band.dat", "w") as f:
            f.writelines(all_files)
        print(len(all_files), " exposures")

if mode == "hdf5_cata":
    # convert the .dat to .hdf5

    total_path = argv[2]
    ori_cat_path = argv[3]

    exposures_candidates = []
    exposures_candidates_band = []
    with open(total_path + "/cat_inform/exposure_name_all_band.dat", "r") as f:
        contents = f.readlines()
    for expo_name in contents:
        expo_nm, band = expo_name.split("\n")[0].split()
        exposures_candidates.append(expo_nm)
        exposures_candidates_band.append(band)

    exposures_candidates_sub = tool_box.alloc(exposures_candidates, numprocs, "seq")[rank]
    exposures_candidates_band_sub = tool_box.alloc(exposures_candidates_band, numprocs, "seq")[rank]

    exposures_candidates_avail_sub = []
    exception_sub = []

    for tag, fns in enumerate(exposures_candidates_sub):
        # read the the field data
        iband = exposures_candidates_band_sub[tag]
        expo_src_path = ori_cat_path + "/%s/%s.cat" % (iband,fns)
        expo_h5_path = total_path + "/cat_hdf5/%s/%s.hdf5" % (iband,fns)
        try:
            src_data = numpy.loadtxt(expo_src_path, dtype=numpy.float32)
            row, col = src_data.shape
        except:
            if os.path.exists(expo_src_path):
                log_inform = "%d Failed in reading %s %d Bytes !\n" % (
                    rank, expo_src_path, os.path.getsize(expo_src_path))
            else:
                log_inform = "%d can't find %s!\n" % (rank, expo_src_path)
            exception_sub.append(log_inform)
            row, col = 0, 0

        if row > 0:
            # Nan check
            idx = numpy.isnan(src_data)
            if idx.sum() > 0:
                print("Find Nan in ", expo_src_path)

            h5f_expo = h5py.File(expo_h5_path, "w")
            h5f_expo["/data"] = src_data
            h5f_expo.close()

            exposures_candidates_avail_sub.append(expo_h5_path + " %s"%exposures_candidates_band_sub[tag])

    exposures_candidates_avail = comm.gather(exposures_candidates_avail_sub, root=0)
    exception_collection = comm.gather(exception_sub, root=0)

    comm.Barrier()

    if rank == 0:
        exception_all = []
        for ec in exception_collection:
            exception_all.extend(ec)

        exposures_collect = []
        for fsb in exposures_candidates_avail:
            exposures_collect.extend(fsb)

        exposures_avail_all = []
        exposures_avail_band = [[] for i in range(len(data_band))]

        for fsb in exposures_collect:
            expo_path, band = fsb.split()
            exposures_avail_all.append(expo_path + "\n")

            for i in range(len(data_band)):
                if band == data_band[i]:
                    exposures_avail_band[i].append(expo_path + "\n")

        with open(total_path + "/cat_inform/exposure_avail_all.dat","w") as f:
            f.writelines(exposures_avail_all)

        for tag, band in enumerate(data_band):
            with open(total_path + "/cat_inform/exposure_avail_%s_band.dat"%band, "w") as f:
                f.writelines(exposures_avail_band[tag])

            exception_all.append("%d available exposures in %s band\n"%(len(exposures_avail_band[tag]), band))

        exception_all.append("Totally: %d available exposures\n"%len(exposures_avail_all))
        with open("log.dat", "w") as f:
            f.writelines(exception_all)


if mode == "hist":

    total_path = argv[2]

    band_num = len(data_band)

    zbin_num = 100
    ra_bin_num, dec_bin_num = 2000, 1000

    z_hist_bin = numpy.linspace(0, 4, zbin_num+1)
    ra_hist_bin = numpy.linspace(0, 360, ra_bin_num+1)
    dec_hist_bin = numpy.linspace(-90, 90, dec_bin_num+1)

    # source position hist
    sub_pos_hist = numpy.zeros((int(band_num*dec_bin_num), ra_bin_num))
    # for the source who has spectral z position hist
    sub_z_pos_hist = numpy.zeros((int(band_num*dec_bin_num), ra_bin_num))
    # the z hist
    sub_zhist = numpy.zeros((int(band_num*2), zbin_num))


    itemsize = MPI.DOUBLE.Get_size()

    pos_element_num = int(band_num*ra_bin_num*dec_bin_num)
    z_element_num = int(band_num*2*zbin_num)

    if rank == 0:
        pos_nbytes = pos_element_num * itemsize
        z_nbytes = z_element_num * itemsize
    else:
        pos_nbytes = 0
        z_nbytes = 0

    win1 = MPI.Win.Allocate_shared(pos_nbytes, itemsize, comm=comm)
    win2 = MPI.Win.Allocate_shared(pos_nbytes, itemsize, comm=comm)
    win3 = MPI.Win.Allocate_shared(z_nbytes, itemsize, comm=comm)

    buf1, itemsize = win1.Shared_query(0)
    buf2, itemsize = win2.Shared_query(0)
    buf3, itemsize = win3.Shared_query(0)

    total_pos_hist = numpy.ndarray(buffer=buf1, dtype='d', shape=(int(band_num*dec_bin_num), ra_bin_num))
    total_z_pos_hist = numpy.ndarray(buffer=buf2, dtype='d', shape=(int(band_num*dec_bin_num), ra_bin_num))
    total_zhist = numpy.ndarray(buffer=buf3, dtype='d', shape=(int(band_num*2), zbin_num))

    comm.Barrier()

    if rank == 0:
        total_pos_hist[:,:] = 0
        total_z_pos_hist[:,:] = 0
        total_zhist[:,:] = 0
    comm.Barrier()


    for tag, band in enumerate(data_band):
        with open(total_path + "/cat_inform/exposure_avail_%s_band.dat"%band, "r") as f:
            file_path = f.readlines()
        file_path_sub = tool_box.alloc(file_path, numprocs)[rank]

        for fps in file_path_sub:
            src = fps.split("\n")[0]

            h5f = h5py.File(src, "r")
            data = h5f["/data"][()]
            h5f.close()

            ra, dec, zp, zs = data[:,0], data[:,1], data[:,16], data[:,17]

            pos_num = numpy.histogram2d(dec, ra, [dec_hist_bin, ra_hist_bin])[0]
            zp_num = numpy.histogram(zp[zp >0], z_hist_bin)[0]
            zs_num = numpy.histogram(zs[zs>0], z_hist_bin)[0]
            idx = zs > 0
            z_pos_num = numpy.histogram2d(dec[idx], ra[idx], [dec_hist_bin, ra_hist_bin])[0]

            st1, st2 = int(tag*dec_bin_num), int(tag*2)
            ed = int((tag+1)*dec_bin_num)

            sub_pos_hist[st1:ed] += pos_num
            sub_z_pos_hist[st1:ed] += z_pos_num
            sub_zhist[st2] += zp_num
            sub_zhist[st2+1] += zs_num

    comm.Barrier()

    for i in range(numprocs):
        total_pos_hist += sub_pos_hist
        total_z_pos_hist += sub_z_pos_hist
        total_zhist += sub_zhist
        comm.Barrier()

    comm.Barrier()
    if rank == 0:

        h5f = h5py.File("./data_hist/hist.hdf5","w")
        h5f["/pos_hist"] = total_pos_hist
        h5f["/z_pos_hist"] = total_z_pos_hist
        h5f["/z_hist"] = total_zhist

        h5f["/ra_bin"] = ra_hist_bin
        h5f["/dec_bin"] = dec_hist_bin
        h5f["/z_bin"] = z_hist_bin
        h5f.close()

        idx = total_pos_hist < 1
        total_pos_hist[idx] = numpy.nan
        idx = total_z_pos_hist < 1
        total_z_pos_hist[idx] = numpy.nan

        for i in range(band_num):
            st, ed = int(i * dec_bin_num), int((i + 1) * dec_bin_num)

            img = Image_Plot(xpad=0.25, ypad=0.25)
            img.subplots(1,3)

            img.axs[0][0].imshow(numpy.flip(total_pos_hist[st:ed],axis=0))
            img.axs[0][1].imshow(numpy.flip(total_z_pos_hist[st:ed],axis=0))

            for j in range(2):
                pos_dec = [0, int(dec_bin_num/4), int(dec_bin_num/2), int(dec_bin_num*3/4), int(dec_bin_num-1)]
                img.set_ticklabel_str(0, j, 0, pos_dec, ["90","45", "0", "-45", "-90"])

                pos_ra = [0,int(ra_bin_num/4), int(ra_bin_num/2), int(ra_bin_num*3/4), int(ra_bin_num-1)]
                img.set_ticklabel_str(0, j, 1, pos_ra, ["0", "90", "180", "270", "360"])

                img.set_label(0,j,0,"Dec [Deg]")
                img.set_label(0,j,1,"RA [Deg]")
            img.axs[0][0].set_title("source hist")
            img.axs[0][1].set_title("source with spectral z hist")


            zpts = (z_hist_bin[1:] + z_hist_bin[:-1])/2

            st = int(i*2)

            img.axs[0][2].plot(zpts, total_zhist[st]/total_zhist[st].sum(), label="Photo Z")
            img.axs[0][2].plot(zpts, total_zhist[st+1]/total_zhist[st+1].sum(), label="Spec Z")
            img.axs[0][2].legend()

            img.set_label(0,2,0,"P(z)dz")
            img.set_label(0,2,1,"Z")
            img.save_img("./data_hist/%s_band.pdf" % data_band[i])
            img.close_img()
    comm.Barrier()