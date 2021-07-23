import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append("%s/work/mylib/"% my_home)
import h5py
import numpy
from mpi4py import MPI
import tool_box
import warnings

warnings.filterwarnings('error')


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

total_path = argv[1]
mode = argv[2]

if mode == "cata_name":
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

        files_nm = os.listdir(total_path + "/cat_ori")
        all_files = []

        for fnm in files_nm:
            if "_all.cat" in fnm:
                all_files.append("%s\n"%fnm.split(".")[0])

        with open(total_path + "/cat_inform/exposure_name.dat", "w") as f:
            f.writelines(all_files)

        print(len(all_files), " exposures")

elif mode == "hdf5_cata":
    # convert the .dat to .hdf5

    exposures_candidates = []
    with open(total_path + "/cat_inform/exposure_name.dat", "r") as f:
        contents = f.readlines()
    for expo_name in contents:
        exposures_candidates.append(expo_name.split("\n")[0])

    exposures_candidates_sub = tool_box.alloc(exposures_candidates, numprocs, "seq")[rank]

    exposures_candidates_avail_sub = []
    exception_sub = []

    for fns in exposures_candidates_sub:
        # read the the field data
        expo_src_path = total_path + "/cat_ori/%s.cat" % fns
        expo_h5_path = total_path + "/cat_hdf5/%s.hdf5" % fns
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

        if row > 0 and col > 1:
            # Nan check
            idx = numpy.isnan(src_data)
            if idx.sum() > 0:
                print("Find Nan in ", expo_src_path)

            h5f_expo = h5py.File(expo_h5_path, "w")
            h5f_expo["/data"] = src_data
            h5f_expo.close()

            exposures_candidates_avail_sub.append(expo_h5_path + "\n")

    exposures_candidates_avail = comm.gather(exposures_candidates_avail_sub, root=0)
    exception_collection = comm.gather(exception_sub, root=0)

    comm.Barrier()
    if rank == 0:
        exception_all = []
        for ec in exception_collection:
            exception_all.extend(ec)

        exposures_avail = []
        for fsb in exposures_candidates_avail:
            exposures_avail.extend(fsb)
        with open(total_path + "/cat_inform/exposure_avail.dat","w") as f:
            f.writelines(exposures_avail)

        exception_all.append("Totally: %d/%d available exposures\n"
                             "The expectation of Z calculated from the P(z) and normalized P(z), "
                             "and the ODDS have been added to the hdf5 cat\n"%(len(exposures_avail),len(exposures_candidates)))
        with open("log.dat","w") as f:
            f.writelines(exception_all)


if mode == "hist":

    zbin_num = 100
    ra_bin_num, dec_bin_num = 500, 250

    z_hist_bin = numpy.linspace(0, 4, zbin_num+1)
    ra_hist_bin = numpy.linspace(0, 360, ra_bin_num+1)
    dec_hist_bin = numpy.linspace(-90, 90, dec_bin_num+1)

    # source position hist
    sub_pos_hist = numpy.zeros((dec_bin_num, ra_bin_num))
    # for the source who has spectral z position hist
    sub_z_pos_hist = numpy.zeros((dec_bin_num, ra_bin_num))
    # the z hist
    sub_zhist = numpy.zeros((2, zbin_num))


    itemsize = MPI.DOUBLE.Get_size()

    pos_element_num = int(ra_bin_num*dec_bin_num)
    z_element_num = int(2*zbin_num)

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
        if i == rank:
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

elif mode == "field_expo":

    if rank == 0:
        field_expos = []
        with open(total_path + "/cat_inform/field_expo.dat", "r") as f:
            cc = f.readlines()
        expo_num = 0
        for c in cc:
            nms = c.split("\n")[0].split("\t")
            temp = [nms[0]]
            for expo_name in nms[1:]:
                expo_path = total_path + "/cat_hdf5/%s_all.hdf5"%expo_name
                if os.path.exists(expo_path):
                    temp.append(expo_name)
                    expo_num += 1
                else:
                    print("Not found %s"%expo_path)
            temp_str = "\t".join(temp) + "\n"
            field_expos.append(temp_str)
        print(expo_num)
        with open(total_path + "/cat_inform/field_expo_avail.dat", "w") as f:
            f.writelines(field_expos)
    comm.Barrier()


elif mode == "position":
    expo_files = []
    with open(total_path + "/cat_inform/field_expo_avail.dat", "r") as f:
        cc = f.readlines()
    field_num = len(cc)

    max_expo_num = 0
    for c in cc:
        nms = c.split("\n")[0].split("\t")
        expo_num = len(nms) - 1
        if expo_num > max_expo_num:
            max_expo_num = expo_num
        expo_files.append(nms[1:])

    expo_pos = numpy.zeros((field_num, max_expo_num, 8), dtype=numpy.float32)

    field_label = [i for i in range(field_num)]
    field_sub = tool_box.alloc(field_label, numprocs,"seq")[rank]

    for i in field_sub:
        for j, expo_name in enumerate(expo_files[i]):
            file_path = total_path + "/cat_hdf5/%s_all.hdf5"%expo_name

            h5f = h5py.File(file_path, "r")
            data = h5f["/data"][()]
            h5f.close()

            ra, dec = data[:,0], data[:,1]
            ra_min, ra_max = ra.min(), ra.max()
            dec_min, dec_max = dec.min(), dec.max()

            dra = (ra_max - ra_min)/2
            ddec = (dec_max - dec_min)/2

            ra_cent = (ra_max + ra_min)/2
            dec_cent = (dec_max + dec_min)/2

            expo_pos[i][j] = ra_cent, dec_cent, ra_min, ra_max, dra, dec_min, dec_max, ddec

    comm.Barrier()

    if rank > 0:
        comm.Send([expo_pos, MPI.FLOAT], dest=0, tag=rank)
    else:
        for procs in range(1, numprocs):
            recvs = numpy.empty((field_num, max_expo_num, 8), dtype=numpy.float32)
            comm.Recv(recvs, source=procs, tag=procs)
            expo_pos += recvs

        h5f = h5py.File(total_path + "/cat_inform/expo_pos.hdf5", "w")
        h5f["/data"] = expo_pos
        h5f.close()

    comm.Barrier()


else:
    # collection

    source_list_nm = argv[2]
    result_nm = argv[3]

    expos = []
    gal_num = []
    with open(total_path + "/"+source_list_nm, "r") as f:
        conts = f.readlines()
    for nm in conts:
        informs = nm.split("\n")
        expos.append(informs[0])
        gal_num.append(informs[2])

    if rank == 0:
        print(len(expos)," exposures")

    my_sub_area_list = tool_box.alloc(expos,numprocs)[rank]
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
        for ir in range(1, numprocs):
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
        print("Totally %d (%d) galaxies"%(data.shape[0], sum(gal_num)))


    comm.Barrier()