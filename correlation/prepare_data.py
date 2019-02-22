import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
from Fourier_Quad import Fourier_Quad
import tool_box
import h5py
from mpi4py import MPI
from sys import argv
import matplotlib.pyplot as plt
import numpy
import time


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cmd = argv[1]
cmds = ["collect","grid"]
if cmd not in cmds:
    if rank == 0:
        print("parameter must be one of ", cmds)
    exit()

area_num = 4

result_source = "result_int"
sect = "cfht"
envs_path = "%s/work/envs/envs.dat"%my_home

gets_item = [["cfht", "cfht_path_catalog", "0"], ["cfht", "cfht_path_data", "0"]]
path_items = tool_box.config(envs_path, ["get", "get"], gets_item)

cata_path, data_path = path_items

# data collection
if cmd == "collect":
    t1 = time.time()
    dicts, fields = tool_box.field_dict(cata_path+"nname.dat")

    if rank == 0:
        h5f_path = data_path + "cata_%s.hdf5" % result_source
        h5f = h5py.File(h5f_path, "w")
        h5f.close()
    # loop the areas,
    # the catalog of each area will be stored in "w_i"
    for i in range(area_num):
        if rank == 0:
            field_tar = []
            for field in fields:
                if "w%d"%(i+1) in field:
                    field_tar.append(field)
            # the allot() returns a list of lists
            field_pool = tool_box.allot(field_tar, cpus)
        else:
            field_pool = None
        # distribution the files
        field_pool = comm.scatter(field_pool, root=0)

        field_count = 0
        for field_name in field_pool:
            field_cat_path = cata_path + "%s/%s/%s_shear.dat"%(field_name, result_source, field_name)
            if os.path.exists(field_cat_path):
                try:
                    cat_arr = numpy.loadtxt(field_cat_path)
                    if field_count == 0:
                        cata_data = cat_arr
                    else:
                        cata_data = numpy.row_stack((cata_data, cat_arr))
                    field_count += 1
                except:
                    print(rank, "%s.dat doesn't exist"%field_name)
        data_sp = cata_data.shape
        data_sps = comm.gather(data_sp, root=0)
        npz_name = data_path+"rank_%d_%d.npz"%(i,rank)
        numpy.savez(npz_name, cata_data)

        if rank == 0:
            data_sps = numpy.array(data_sps)
            rows, cols = data_sps[:, 0], data_sps[:, 1]
            displ = []
            count = rows*cols
            for j in range(cpus):
                displ.append(count[0:j].sum())
            count = count.tolist()
            recv_buffer = numpy.empty((rows.sum(), cols[0]))
        else:
            count = None
            displ = None
            recv_buffer = None
        count = comm.bcast(count, root=0)
        displ = comm.bcast(displ, root=0)
        comm.Gatherv(cata_data, [recv_buffer, count, displ, MPI.DOUBLE], root=0)
        if rank == 0:
            h5f = h5py.File(h5f_path)
            h5f["/w_%d"%i] = recv_buffer
            h5f.close()

    t2 = time.time()
    if rank == 0:
        print(t2-t1)

# make grid
if cmd == "grid":

    fq = Fourier_Quad(12, 1123)

    t1 = time.time()

    gets_item = [["fresh_para_idx", "nstar", "0"], ["fresh_para_idx", "flux_alt", "0"],
                 ["fresh_para_idx", "ra", "0"], ["fresh_para_idx", "dec", "0"],
                 ["fresh_para_idx", "gf1", "0"], ["fresh_para_idx", "gf2", "0"],
                 ["fresh_para_idx", "g1", "0"], ["fresh_para_idx", "g2", "0"],
                 ["fresh_para_idx", "de", "0"], ["fresh_para_idx", "h1", "0"],
                 ["fresh_para_idx", "h2", "0"]]
    gets = ["get" for i in range(len(gets_item))]
    para_items = tool_box.config(envs_path, gets, gets_item)

    nstar_lb = int(para_items[0])
    flux_alt_lb = int(para_items[1])

    ra_lb = int(para_items[2])
    dec_lb = int(para_items[3])

    field_g1_lb = int(para_items[4])
    field_g2_lb = int(para_items[5])

    mg1_lb = int(para_items[6])
    mg2_lb = int(para_items[7])
    mn_lb = int(para_items[8])
    mu_lb = int(para_items[9])
    mv_lb = int(para_items[10])

    flux_alt_thresh = 4
    nstar_thresh = 12
    field_g1_bound = 0.005
    field_g2_bound = 0.0075

    data_col = 7
    mg_bin_num = 8
    block_scale = 4  # arcsec
    margin = 0.1 * block_scale

    radius_scale_bin_num = 21
    radius_scale = numpy.zeros((radius_scale_bin_num+1,), dtype=numpy.double)
    for i in range(1,radius_scale_bin_num+1):
        radius_scale[i] = 1+10 ** (0.1 * i)

    g_hat_bin_num = 60
    g_hat_bin = numpy.linspace(-0.0011, 0.0011, g_hat_bin_num)

    cf_cata_data_path = data_path + "cf_cata_%s.hdf5"%result_source

    if rank == 0:
        h5f = h5py.File(cf_cata_data_path,"w")

        # the radius bin for correlation function calculation
        h5f["/radius_bin"] = radius_scale
        h5f["/radius_bin"].attrs["shape"] = radius_scale_bin_num+1

        # the guess of gg correlation for the algorithm
        h5f["/g_hat_bin"] = g_hat_bin
        h5f["/g_hat_bin"].attrs["shape"] = g_hat_bin_num

        # the number of the sky area
        h5f.create_group("/grid")
        h5f["/grid"].attrs["max_area"] = area_num

        h5f.close()

    comm.Barrier()

    cata_data_path = data_path + "cata_%s.hdf5" % result_source

    if rank < area_num:

        h5f = h5py.File(cata_data_path, "r")
        ori_cat_data = h5f["/w_%d"%rank].value
        h5f.close()

        # cut off
        flux_alt_idx = ori_cat_data[:, flux_alt_lb] >= flux_alt_thresh
        nstar_idx = ori_cat_data[:, nstar_lb] >= nstar_thresh

        fg1 = numpy.abs(ori_cat_data[:, field_g1_lb])
        fg2 = numpy.abs(ori_cat_data[:, field_g2_lb])

        fg1_idx = fg1 <= field_g1_bound
        fg2_idx = fg2 <= field_g2_bound

        # the esitmators
        mg1 = ori_cat_data[:, mg1_lb][flux_alt_idx & nstar_idx & fg1_idx & fg2_idx]
        mg2 = ori_cat_data[:, mg2_lb][flux_alt_idx & nstar_idx & fg1_idx & fg2_idx]
        mn = ori_cat_data[:, mn_lb][flux_alt_idx & nstar_idx & fg1_idx & fg2_idx]
        mu = ori_cat_data[:, mu_lb][flux_alt_idx & nstar_idx & fg1_idx & fg2_idx]
        mv = ori_cat_data[:, mv_lb][flux_alt_idx & nstar_idx & fg1_idx & fg2_idx]
        # !!! be careful with the sign of mn, mu, mv

        # ra dec mg1 mg2 mn mu mv
        buffer = numpy.zeros((len(mg1), data_col))

        # correct the field distortion & the "c" and "m" !!!!
        field_g1 = ori_cat_data[:, field_g1_lb][flux_alt_idx & nstar_idx & fg1_idx & fg2_idx]
        field_g2 = ori_cat_data[:, field_g2_lb][flux_alt_idx & nstar_idx & fg1_idx & fg2_idx]
        buffer[:, 2] = mg1 - field_g1 * (mn + mu) - field_g2 * mv
        buffer[:, 3] = mg2 - field_g2 * (mn - mu) - field_g1 * mv
        buffer[:, 4] = mn
        buffer[:, 5] = mu
        buffer[:, 6] = mv

        # the bin for mg's
        mg1_bin = fq.set_bin(buffer[:, 2], mg_bin_num)
        mg2_bin = fq.set_bin(buffer[:, 3], mg_bin_num)

        # the Ra & Dec, convert degree to arcmin
        ra = ori_cat_data[:, ra_lb][flux_alt_idx & nstar_idx & fg1_idx & fg2_idx] * 60
        dec = ori_cat_data[:, dec_lb][flux_alt_idx & nstar_idx & fg1_idx & fg2_idx] * 60
        buffer[:, 0] = ra
        buffer[:, 1] = dec

        # the grid: x-axis--ra, y-axis--dec
        ra_min, ra_max = ra.min() - margin, ra.max() + margin
        dec_min, dec_max = dec.min() - margin, dec.max() + margin

        rows = int( (dec_max - dec_min)/block_scale + 1 )
        cols = int( (ra_max - ra_min)/block_scale + 1 )
        grid_num = rows*cols
        grid_shape = numpy.array([rows, cols], dtype=numpy.intc)

        num_in_block = []
        block_list = []
        block_names = []
        boundy = []
        boundx = []

        # find each block
        fig = plt.figure(figsize=(int(cols*1.0/rows*14), 14))
        ax = fig.add_subplot(111)
        for row in range(rows):

            y1 = dec_min + row * block_scale
            y2 = dec_min + (row + 1) * block_scale
            ir1 = dec >= y1
            ir2 = dec < y2

            for col in range(cols):

                x1 = ra_min + col*block_scale
                x2 = ra_min + (col+1)*block_scale
                ic1 = ra >= x1
                ic2 = ra < x2

                block = buffer[ir1&ir2&ic1&ic2]
                block_sp = block.shape

                block_list.append(block)
                num_in_block.append(block_sp[0])
                block_names.append("/grid/w_%d/row_%d/col_%d"%(rank, row, col))

                boundy.append([y1, y1, y2, y2])
                boundx.append([x1, x2, x1, x2])

                # scatters for checking
                if block_sp[0] > 0:
                    # x: Ra   y: Dec
                    ax.scatter(block[:,0], block[:,1], s=0.3)
        print(rank, ra_min, ra_max, dec_min, dec_max)
        for i in range(rows + 1):
            # ax.plot([x1, x2..],[y1, y2,..])
            # the horizontal line
            ax.plot([ra_min, ra_min+cols * block_scale], [dec_min+i*block_scale, dec_min+i*block_scale], c="black", linewidth=1)
        for j in range(cols + 1):
            # the vertical line
            ax.plot([ra_min+j*block_scale, ra_min+j*block_scale], [dec_min, dec_min+rows * block_scale], c="black", linewidth=1)
        ax.set_ylabel("DEC.")
        ax.set_xlabel("R.A.")
        pic_nm = data_path + "w%d_ra_dec.png" %rank
        plt.savefig(pic_nm)
        plt.close()

        max_sp = numpy.array([max(num_in_block), data_col],dtype=numpy.intc)

        num_in_block = numpy.array(num_in_block, dtype=numpy.intc)
        block_start = numpy.zeros((grid_num,), dtype=numpy.intc)
        block_end = numpy.zeros((grid_num,), dtype=numpy.intc)
        # for the memory saving
        block_start_s = numpy.zeros((grid_num,), dtype=numpy.intc)
        block_end_s = numpy.zeros((grid_num,), dtype=numpy.intc)

        boundy = numpy.array(boundy, dtype=numpy.double)
        boundx = numpy.array(boundx, dtype=numpy.double)

        final_data = numpy.zeros((max_sp[0]*grid_num, data_col),dtype=numpy.double)
        count = 0
        for i in range(grid_num):
            block_start[i] = i*max_sp[0]
            block_end[i] = (i+1)*max_sp[0]
            # for the memory saving
            if i == 0:
                block_start_s[i] = 0
            else:
                block_start_s[i] = sum(num_in_block[:i])
            block_end_s[i] = block_start_s[i] + num_in_block[i]

            if num_in_block[i] != 0:
                final_data[i*max_sp[0]:i*max_sp[0]+num_in_block[i]] = block_list[i]
                # for the memory saving
                if count == 0:
                    final_data_s = block_list[i]
                else:
                    final_data_s = numpy.row_stack((final_data_s, block_list[i]))
                count += 1

        final_ra = final_data[:, 0]
        final_dec = final_data[:, 1]
        final_mg1 = final_data[:, 2]
        final_mg2 = final_data[:, 3]
        final_mn = final_data[:, 4]
        final_mu = final_data[:, 5]
        final_mv = final_data[:, 6]
        num = ra.shape

        # for the memory saving
        final_ra_s = final_data_s[:, 0]
        final_dec_s = final_data_s[:, 1]
        final_mg1_s = final_data_s[:, 2]
        final_mg2_s = final_data_s[:, 3]
        final_mn_s = final_data_s[:, 4]
        final_mu_s = final_data_s[:, 5]
        final_mv_s = final_data_s[:, 6]
        compact_len = len(final_ra_s)

        for i in range(area_num):
            if i == rank:
                h5f = h5py.File(cf_cata_data_path)

                h5f.create_group("/grid/w_%d"%i)
                h5f["/grid/w_%d"%i].attrs["grid_shape"] = grid_shape
                h5f["/grid/w_%d"%i].attrs["max_block_size"] = max_sp
                h5f["/grid/w_%d"%i].attrs["block_scale"] = block_scale

                h5f["/grid/w_%d/boundy" % i] = boundy
                h5f["/grid/w_%d/boundy" % i].attrs["shape"] = numpy.array([grid_num, 4], dtype=numpy.intc)
                h5f["/grid/w_%d/boundx" % i] = boundx
                h5f["/grid/w_%d/boundx" % i].attrs["shape"] = numpy.array([grid_num, 4], dtype=numpy.intc)

                h5f["/grid/w_%d/num_in_block" % i] = num_in_block
                h5f["/grid/w_%d/num_in_block" % i].attrs["shape"] = numpy.array([grid_num], dtype=numpy.intc)

                h5f["/grid/w_%d/block_start" % i] = block_start
                h5f["/grid/w_%d/block_start" % i].attrs["shape"] = numpy.array([grid_num], dtype=numpy.intc)
                h5f["/grid/w_%d/block_end" % i] = block_end
                h5f["/grid/w_%d/block_end" % i].attrs["shape"] = numpy.array([grid_num], dtype=numpy.intc)
                # for the memory saving
                h5f["/grid/w_%d/block_start_s" % i] = block_start_s
                h5f["/grid/w_%d/block_start_s" % i].attrs["shape"] = numpy.array([compact_len], dtype=numpy.intc)
                h5f["/grid/w_%d/block_end_s" % i] = block_end_s
                h5f["/grid/w_%d/block_end_s" % i].attrs["shape"] = numpy.array([compact_len], dtype=numpy.intc)

                h5f["/grid/w_%d/G1_bin" % i] = mg1_bin
                h5f["/grid/w_%d/G1_bin" % i].attrs["shape"] = numpy.array([mg_bin_num], dtype=numpy.intc)
                h5f["/grid/w_%d/G2_bin" % i] = mg2_bin
                h5f["/grid/w_%d/G2_bin" % i].attrs["shape"] = numpy.array([mg_bin_num], dtype=numpy.intc)

                h5f["/grid/w_%d/data/RA" % i] = final_ra
                h5f["/grid/w_%d/data/RA" % i].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                h5f["/grid/w_%d/data/DEC" % i] = final_dec
                h5f["/grid/w_%d/data/DEC" % i].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                h5f["/grid/w_%d/data/G1" % i] = final_mg1
                h5f["/grid/w_%d/data/G1" % i].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                h5f["/grid/w_%d/data/G2" % i] = final_mg2
                h5f["/grid/w_%d/data/G2" % i].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                h5f["/grid/w_%d/data/N" % i] = final_mn
                h5f["/grid/w_%d/data/N" % i].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                h5f["/grid/w_%d/data/U" % i] = final_mu
                h5f["/grid/w_%d/data/U" % i].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                h5f["/grid/w_%d/data/V" % i] = final_mv
                h5f["/grid/w_%d/data/V" % i].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                # for the memory saving
                h5f["/grid/w_%d/data/RA_s" % i] = final_ra_s
                h5f["/grid/w_%d/data/RA_s" % i].attrs["shape"] = numpy.array([compact_len], dtype=numpy.intc)

                h5f["/grid/w_%d/data/DEC_s" % i] = final_dec_s
                h5f["/grid/w_%d/data/DEC_s" % i].attrs["shape"] = numpy.array([compact_len], dtype=numpy.intc)

                h5f["/grid/w_%d/data/G1_s" % i] = final_mg1_s
                h5f["/grid/w_%d/data/G1_s" % i].attrs["shape"] = numpy.array([compact_len], dtype=numpy.intc)

                h5f["/grid/w_%d/data/G2_s" % i] = final_mg2_s
                h5f["/grid/w_%d/data/G2_s" % i].attrs["shape"] = numpy.array([compact_len], dtype=numpy.intc)

                h5f["/grid/w_%d/data/N_s" % i] = final_mn_s
                h5f["/grid/w_%d/data/N_s" % i].attrs["shape"] = numpy.array([compact_len], dtype=numpy.intc)

                h5f["/grid/w_%d/data/U_s" % i] = final_mu_s
                h5f["/grid/w_%d/data/U_s" % i].attrs["shape"] = numpy.array([compact_len], dtype=numpy.intc)

                h5f["/grid/w_%d/data/V_s" % i] = final_mv_s
                h5f["/grid/w_%d/data/V_s" % i].attrs["shape"] = numpy.array([compact_len], dtype=numpy.intc)

                for ig in range(grid_num):
                    h5f.create_group(block_names[ig])
                    if num_in_block[ig] != 0:

                        num = max_sp[0]
                        h5f[block_names[ig] + "/RA"] = block_list[ig][:, 0]
                        h5f[block_names[ig] + "/RA"].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                        h5f[block_names[ig] + "/DEC"] = block_list[ig][:, 1]
                        h5f[block_names[ig] + "/DEC"].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                        h5f[block_names[ig] + "/G1"] = block_list[ig][:, 2]
                        h5f[block_names[ig] + "/G1"].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                        h5f[block_names[ig] + "/G2"] = block_list[ig][:, 3]
                        h5f[block_names[ig] + "/G2"].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                        h5f[block_names[ig] + "/N"] = block_list[ig][:, 4]
                        h5f[block_names[ig] + "/N"].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                        h5f[block_names[ig] + "/U"] = block_list[ig][:, 5]
                        h5f[block_names[ig] + "/U"].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                        h5f[block_names[ig] + "/V"] = block_list[ig][:, 6]
                        h5f[block_names[ig] + "/V"].attrs["shape"] = numpy.array([num], dtype=numpy.intc)

                h5f.close()
            comm.Barrier()

    t2 = time.time()
    if rank == 0:
        print(t2-t1)









