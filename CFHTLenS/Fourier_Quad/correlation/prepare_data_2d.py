import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
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

gets_item = [["fresh_para_idx", "g1", "0"], ["fresh_para_idx", "g2", "0"],
             ["fresh_para_idx", "de", "0"], ["fresh_para_idx", "h1", "0"],
             ["fresh_para_idx", "h2", "0"]]
gets = ["get" for i in range(len(gets_item))]
para_items = tool_box.config(envs_path, gets, gets_item)

mg1_lb = int(para_items[0])
mg2_lb = int(para_items[1])
mn_lb = int(para_items[2])
mu_lb = int(para_items[3])
mv_lb = int(para_items[4])

# data collection
if cmd == "collect":

    t1 = time.time()
    dicts, fields = tool_box.field_dict(cata_path+"nname.dat")
    h5f_path = data_path + "cata_%s.hdf5" % result_source

    if rank == 0:
        h5f = h5py.File(h5f_path, "w")
        h5f.close()
        buffer = []

    # loop the areas,
    # the catalog of each area will be stored in "w_i"
    num = []
    for area_id in range(1,area_num+1):

        read_err = 0

        # distribution the files
        field_tar = []
        for field in fields:
            if "w%d" % area_id in field:
                field_tar.append(field)
        # make sure that each thread gets 1 field to read at least
        # if len(field_tar) < cpus:
        #     print("Too many threads. %d at most"%(len(field_tar)))
        #     exit()
        field_pool = tool_box.allot(field_tar, cpus)[rank]

        # check
        anothers = ["w%d"%i for i in range(1,area_num+1) if i != area_id]
        for field_name in field_pool:
            for ext_field in anothers:
                if ext_field in field_name:
                    print("WRONG",rank, area_id, field_name)
                    exit()

        field_count = 0
        cata_data = None
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
        # in case of that some thread gets nothing from the file (non-existed)
        if cata_data is not None:
            data_sp = cata_data.shape
        else:
            cata_data = numpy.array([[1]])
            data_sp = cata_data.shape

        data_sps = comm.gather(data_sp, root=0)
        num.append(data_sp)
        # npz_name = data_path+"rank_%d_%d.npz"%(i, rank)
        # numpy.savez(npz_name, cata_data)

        if rank > 0:
            comm.Send(cata_data, dest=0, tag=rank)
        else:
            if data_sp[0] > 1 and data_sp[1] > 1:
                stack_pool = [cata_data]
            else:
                stack_pool = []

            for procs in range(1, cpus):
                recvs = numpy.empty(data_sps[procs], dtype=numpy.double)
                comm.Recv(recvs, source=procs, tag=procs)
                if data_sps[procs][0] > 1 and data_sps[procs][1] > 1:
                    stack_pool.append(recvs)

            for stack_count, arr in enumerate(stack_pool):
                if stack_count == 0:
                    recv_buffer = arr
                else:
                    recv_buffer = numpy.row_stack((recv_buffer, arr))

        # not safe because some threads will not get data
        # if rank == 0:
        #     data_sps = numpy.array(data_sps)
        #     rows, cols = data_sps[:, 0], data_sps[:, 1]
        #     displ = []
        #     count = rows*cols
        #     for j in range(cpus):
        #         displ.append(count[0:j].sum())
        #     count = count.tolist()
        #     recv_buffer = numpy.empty((rows.sum(), cols[0]))
        # else:
        #     count = None
        #     displ = None
        #     recv_buffer = None
        # count = comm.bcast(count, root=0)
        # displ = comm.bcast(displ, root=0)
        # comm.Gatherv(cata_data, [recv_buffer, count, displ, MPI.DOUBLE], root=0)

            h5f = h5py.File(h5f_path)

            # each shear estimator is 2 times the one in the paper,
            # Zhang et al. 2017 ApJ 834,
            # G1 = 2*G1_paper, G2 = 2*G2_paper
            # And N = -2N_paper, U = -2U_paper, V = -2V_paper
            recv_buffer[:,mg1_lb] = recv_buffer[:,mg1_lb] / 2
            recv_buffer[:,mg2_lb] = recv_buffer[:,mg2_lb] / 2
            recv_buffer[:,mn_lb] = recv_buffer[:,mn_lb] / 2
            recv_buffer[:,mu_lb] = -recv_buffer[:,mu_lb] / 2
            recv_buffer[:,mv_lb] = -recv_buffer[:,mv_lb] / 2

            h5f["/w_%d" % area_id] = recv_buffer
            h5f.close()
            print("RECV: ", recv_buffer.shape)
    print(rank, num)
    # stack the sub-catalogs from each area
    if rank == 0:
        fig = plt.figure(figsize=(14, 14))
        h5f = h5py.File(h5f_path)
        for area_id in range(1,area_num+1):
            temp = h5f["/w_%d" % area_id].value

            print(area_id, temp.shape,temp[:,12].min()*60, temp[:,12].max()*60,
                  temp[:,13].min()*60, temp[:,13].max()*60)

            ax = fig.add_subplot(220 + area_id)
            ax.scatter(temp[:,12]*60, temp[:,13]*60, s=0.3)
            if area_id == 1:
                data = temp
            else:
                data = numpy.row_stack((data, temp))
        plt.savefig(data_path + "W_ra_dec.png")
        plt.close()
        h5f["/total"] = data
        h5f.close()
    t2 = time.time()
    if rank == 0:
        print(t2-t1)

# make grid
if cmd == "grid":

    fq = Fourier_Quad(12, 1123)

    gets_item = [["fresh_para_idx", "nstar", "0"], ["fresh_para_idx", "flux_alt", "0"],
                 ["fresh_para_idx", "ra", "0"], ["fresh_para_idx", "dec", "0"],
                 ["fresh_para_idx", "gf1", "0"], ["fresh_para_idx", "gf2", "0"],
                 ["fresh_para_idx", "g1", "0"], ["fresh_para_idx", "g2", "0"],
                 ["fresh_para_idx", "de", "0"], ["fresh_para_idx", "h1", "0"],
                 ["fresh_para_idx", "h2", "0"], ["fresh_para_idx", "total_area", "0"]]
    gets = ["get" for i in range(len(gets_item))]
    para_items = tool_box.config(envs_path, gets, gets_item)

    logger = tool_box.get_logger("./logs/%d_logs.dat"%rank)

    nstar_lb = int(para_items[0])
    flux_alt_lb = int(para_items[1])
    total_area_lb = int(para_items[11])

    ra_lb = int(para_items[2])
    dec_lb = int(para_items[3])

    field_g1_lb = int(para_items[4])
    field_g2_lb = int(para_items[5])

    mg1_lb = int(para_items[6])
    mg2_lb = int(para_items[7])
    mn_lb = int(para_items[8])
    mu_lb = int(para_items[9])
    mv_lb = int(para_items[10])

    flux_alt_thresh = 4.79
    nstar_thresh = 14
    total_area_thresh = 7
    field_g1_bound = 0.005
    field_g2_bound = 0.0075

    data_col = 7
    mg_bin_num = 8
    block_scale = 4  # arcsec
    margin = 0.1 * block_scale

    cf_cata_data_path = data_path + "cf_cata_%s.hdf5" % result_source
    logger.info("Start....")
    for area_id in range(1,area_num+1):

        logger.info("Start area: %d"%area_id)

        t1 = time.time()

        logger.info("Prepare the data")

        # open the original data
        cata_data_path = data_path + "cata_%s.hdf5" % result_source
        h5f = h5py.File(cata_data_path, "r")
        ori_cat_data = h5f["/w_%d"%area_id].value
        h5f.close()

        # cut off
        flux_alt_idx = ori_cat_data[:, flux_alt_lb] >= flux_alt_thresh
        nstar_idx = ori_cat_data[:, nstar_lb] >= nstar_thresh
        total_area_idx = ori_cat_data[:, total_area_lb] >= total_area_thresh

        fg1 = numpy.abs(ori_cat_data[:, field_g1_lb])
        fg2 = numpy.abs(ori_cat_data[:, field_g2_lb])

        fg1_idx = fg1 <= field_g1_bound
        fg2_idx = fg2 <= field_g2_bound
        cut_idx =flux_alt_idx & nstar_idx & total_area_idx & fg1_idx & fg2_idx
        # the esitmators
        mg1 = ori_cat_data[:, mg1_lb][cut_idx]
        mg2 = ori_cat_data[:, mg2_lb][cut_idx]
        mn = ori_cat_data[:, mn_lb][cut_idx]
        mu = ori_cat_data[:, mu_lb][cut_idx]
        mv = ori_cat_data[:, mv_lb][cut_idx]
        # !!! be careful with the sign of mn, mu, mv

        # ra dec mg1 mg2 mn mu mv
        data_row = len(mg1)
        buffer = numpy.zeros((data_row, data_col), dtype=numpy.double)

        # correct the field distortion & the "c" and "m" !!!!
        field_g1 = ori_cat_data[:, field_g1_lb][cut_idx]
        field_g2 = ori_cat_data[:, field_g2_lb][cut_idx]

        # the Ra & Dec, convert degree to arcmin
        ra = ori_cat_data[:, ra_lb][cut_idx] * 60
        dec = ori_cat_data[:, dec_lb][cut_idx]* 60

        buffer[:,0] = ra
        buffer[:,1] = dec
        # each shear estimator is 2 times the one in the paper,
        # Zhang et al. 2017 ApJ 834,
        # G1 = 2*G1, G2 = 2*G2
        # And N = -2N, U = -2U, V = -2V
        # the sign has been corrected in the "collect" process
        # correction: flux_alt>= 4.79, nstar>= 14, area>=7
        # c2 ~ 0.000292
        c2_correction = 0.000292
        buffer[:,2] = mg1
        buffer[:,3] = mg2 - c2_correction*(mn - mu)
        buffer[:,4] = mn
        buffer[:,5] = mu
        buffer[:,6] = mv

        # the grid: x-axis--ra, y-axis--dec
        ra_min, ra_max = ra.min() - margin, ra.max() + margin
        dec_min, dec_max = dec.min() - margin, dec.max() + margin
        logger.info("R.A.: %.4f ~ %.4f  DEC: %.4f ~ %.4f" % (ra_min, ra_max, dec_min, dec_max))
        # the bin for mg's
        mg1_bin = fq.set_bin(buffer[:, 2], mg_bin_num)
        mg2_bin = fq.set_bin(buffer[:, 3], mg_bin_num)

        # the grid
        grid_rows = int((dec_max - dec_min) / block_scale + 1)
        grid_cols = int((ra_max - ra_min) / block_scale + 1)
        grid_num = grid_rows * grid_cols

        logger.info("Grid: %d x %d"%(grid_rows, grid_cols))

        boundy = numpy.zeros((grid_num, 4), dtype=numpy.double)
        boundx = numpy.zeros((grid_num, 4), dtype=numpy.double)
        for row in range(grid_rows):
            for col in range(grid_cols):

                grid_id = row*grid_cols + col

                x1 = ra_min + col * block_scale
                x2 = ra_min + (col + 1) * block_scale
                y1 = dec_min + row * block_scale
                y2 = dec_min + (row + 1) * block_scale

                boundy[grid_id] = y1, y1, y2, y2
                boundx[grid_id] = x1, x2, x1, x2
        if rank == 0:

            if area_id == 1:
                # each time the program starts the file will be truncated
                h5f = h5py.File(cf_cata_data_path, "w")
                # the radius bin for correlation function calculation
                radius_scale_bin_num = 19
                radius_scale = numpy.zeros((radius_scale_bin_num + 1,), dtype=numpy.double)
                for i in range(1, radius_scale_bin_num + 1):
                    radius_scale[i] = 1 + 10 ** (0.1 * i)
                h5f["/radius_bin"] = radius_scale
                h5f["/radius_bin"].attrs["shape"] = numpy.array([radius_scale_bin_num + 1, 1], dtype=numpy.intc)

                # the guess of gg correlation for the algorithm
                g_hat_bin_num = 100
                g_hat_bin = numpy.linspace(10**(-7), 0.00014, g_hat_bin_num)
                h5f["/g_hat_bin"] = - g_hat_bin
                h5f["/g_hat_bin"].attrs["shape"] = numpy.array([g_hat_bin_num, 1], dtype=numpy.intc)

                # the number of the sky area
                h5f.create_group("/grid")
                h5f["/grid"].attrs["max_area"] = numpy.array([area_num, 1], dtype=numpy.intc)
            else:
                h5f = h5py.File(cf_cata_data_path, "r+")

            logger.info("Open the hdf5: %s" % cf_cata_data_path)
            h5f.create_group("/grid/w_%d"%area_id)

            h5f["/grid/w_%d" % area_id].attrs["grid_shape"] = numpy.array([grid_rows, grid_cols], dtype=numpy.intc)
            h5f["/grid/w_%d" % area_id].attrs["block_scale"] = numpy.array([block_scale, 1], dtype=numpy.double)

            h5f["/grid/w_%d/boundy"%area_id] = boundy
            h5f["/grid/w_%d/boundy"%area_id].attrs["shape"] = numpy.array([grid_num, 4], dtype=numpy.intc)
            h5f["/grid/w_%d/boundx"%area_id] = boundx
            h5f["/grid/w_%d/boundx"%area_id].attrs["shape"] = numpy.array([grid_num, 4], dtype=numpy.intc)

            h5f["/grid/w_%d/G1_bin"%area_id] = mg1_bin
            h5f["/grid/w_%d/G1_bin"%area_id].attrs["shape"] = numpy.array([mg_bin_num+1, 1], dtype=numpy.intc)
            h5f["/grid/w_%d/G2_bin"%area_id] = mg2_bin
            h5f["/grid/w_%d/G2_bin"%area_id].attrs["shape"] = numpy.array([mg_bin_num+1, 1], dtype=numpy.intc)

            h5f.close()
        logger.info("Close the hdf5: %s"%cf_cata_data_path)

        grid_num = grid_rows*grid_cols

        logger.info("Create the shared buffer, %d x %d"%(data_row, data_col))

        # then, distribute the grids to threads
        tasks_list = [i for i in range(grid_num)]
        # each thread gets its parts,
        # the list of grid id (row_i*col_j), e.g. [ 10,11,12,17, ...]
        # the sequence of the blocks will be preserved
        my_tasks = tool_box.allot(tasks_list, cpus, "seq")[rank]
        if cpus > grid_num:
            print("Too many cpus (%d) > grids (%d)"%(cpus, grid_num))
            exit()

        # the sub-list of a total
        # they will merged into the final one
        sub_block_list = []
        sub_num_in_block = []
        sub_block_names = []

        logger.info("Rank: %d start block searching. %d blocks (%d ~ %d)"%(rank, len(my_tasks), my_tasks[0], my_tasks[-1]))
        # find each block
        for grid_id in my_tasks:

            row, col = divmod(grid_id, grid_cols)

            x1 = ra_min + col * block_scale
            x2 = ra_min + (col + 1) * block_scale
            ic1 = buffer[:, 0] >= x1
            ic2 = buffer[:, 0] < x2

            y1 = dec_min + row * block_scale
            y2 = dec_min + (row + 1) * block_scale
            ir1 = buffer[:, 1] >= y1
            ir2 = buffer[:, 1] < y2

            block = buffer[ir1 & ir2 & ic1 & ic2]
            block_sp = block.shape

            sub_block_list.append(block)
            sub_num_in_block.append(block_sp[0])
            sub_block_names.append("/grid/w_%d/row_%d/col_%d" % (area_id, row, col))
        logger.info("Rank: %d finish block searching." %rank)

        # each thread write the (sub) blocks to the hdf5 file
        for rank_id in range(cpus):
            if rank_id == rank:
                h5f = h5py.File(cf_cata_data_path,"r+")
                esti_nm = ["/RA", "/DEC", "/G1", "/G2", "/N", "/U", "/V"]
                for ir in range(len(my_tasks)):
                    for inm in range(data_col):
                        set_nm = sub_block_names[ir]+esti_nm[inm]
                        # some blocks contains nothing
                        if sub_num_in_block[ir] > 0:
                            h5f[set_nm] = sub_block_list[ir][:, inm]
                        else:
                            h5f.create_group(set_nm)
                        h5f[set_nm].attrs["shape"] = numpy.array([sub_num_in_block[ir], 1],dtype=numpy.intc)
                h5f.close()
            comm.Barrier()
        logger.info("Rank: %d write the blocks into hdf5 %s." % (rank, cf_cata_data_path))

        # stack the find blocks into a entirety,
        # and it will be merged into the final catalog
        sub_stack_block_size = sum(sub_num_in_block)

        # in case of that some threads may get nothing
        if sub_stack_block_size > 0:
            sub_sp = (sub_stack_block_size, data_col)
            sub_block_stack = numpy.zeros(sub_sp, dtype=numpy.double)
            # for checking
            start_end_pool = []
            for ir in range(len(my_tasks)):
                if ir == 0:
                    std_idx = 0
                else:
                    std_idx = sum(sub_num_in_block[:ir])
                if sub_num_in_block[ir] > 0:
                    sub_block_stack[std_idx: std_idx+sub_num_in_block[ir]] = sub_block_list[ir]
                    std_idx_ = (std_idx, std_idx+sub_num_in_block[ir])
                    if std_idx_ in start_end_pool:
                        # the start & end idx of a block in the stacked array is unique
                        # if not, something goes wrong
                        print("Rank: %d Area: %d. Something went wrong in stack process"%(rank, area_id))
                        exit()
                    else:
                        start_end_pool.append(std_idx_)
        else:
            # if the thread gets nothing
            sub_sp = (1, 1)
            sub_block_stack = numpy.zeros(sub_sp, dtype=numpy.double)

        # send the sub_stack_block to rank 0
        # the sequence of the sub-list follows the the rank sequence
        stack_block_sps = comm.gather(sub_sp, root=0)
        sub_num_in_block_g = comm.gather(sub_num_in_block, root=0)
        logger.info("Rank: %d start to send stacked block"%rank)

        if rank > 0:
            comm.Send(sub_block_stack, dest=0, tag=rank)
        else:
            # rank 0 creates a buffer
            if sub_stack_block_size > 0:
                data_pool = [sub_block_stack]
            else:
                data_pool = []

            # receive the data from the others and put them into the buffer
            for procs in range(1, cpus):
                recvs = numpy.empty(stack_block_sps[procs], dtype=numpy.double)
                comm.Recv(recvs, source=procs, tag=procs)
                if stack_block_sps[procs][1] > 1:
                    data_pool.append(recvs)

            # stack the data in the buffer
            for tag, data_arr in enumerate(data_pool):
                if tag == 0:
                    final_data = data_arr
                else:
                    final_data = numpy.row_stack((final_data, data_arr))
                # all the blocks have been stacked into "final_data"

            logger.info("Rank: %d receive all the stacked block" % rank)
            # the sub_num_in_block_g is a list of lists,
            # a list of "sub_num_in_block"(list) from each thread
            num_in_block = []
            for sub_list in sub_num_in_block_g:
                num_in_block.extend(sub_list)
            max_block_size = numpy.array([max(num_in_block), data_col], dtype=numpy.intc)

            num_in_block = numpy.array(num_in_block, dtype=numpy.intc)
            block_start = numpy.zeros((grid_num,), dtype=numpy.intc)
            block_end = numpy.zeros((grid_num,), dtype=numpy.intc)
            # calculate the start & end row of each block
            for ig in range(grid_num):
                if ig == 0:
                    block_start[ig] = 0
                else:
                    block_start[ig] = sum(num_in_block[:ig])
                block_end[ig] = block_start[ig] + num_in_block[ig]

            logger.info("Rank: %d Plot the grid" % rank)
            fig = plt.figure(figsize=(int(grid_cols*1.0/grid_rows*14),14))
            ax = fig.add_subplot(111)
            # plot the grid lines
            for i in range(grid_rows + 1):
                # ax.plot([x1, x2..],[y1, y2,..])
                # the horizontal line
                ax.plot([ra_min, ra_min + grid_cols * block_scale],
                        [dec_min + i * block_scale, dec_min + i * block_scale], c="black",linestyle="--", linewidth=0.8)
            for j in range(grid_cols + 1):
                # the vertical line
                ax.plot([ra_min + j * block_scale, ra_min + j * block_scale],
                        [dec_min, dec_min + grid_rows * block_scale], c="black",linestyle="--", linewidth=0.8)

            for ig in range(grid_num):
                # x: Ra   y: Dec
                if num_in_block[ig] > 0:
                    ax.scatter(final_data[block_start[ig]:block_end[ig], 0],
                               final_data[block_start[ig]:block_end[ig], 1], s=0.3)
            ax.set_ylabel("DEC.")
            ax.set_xlabel("R.A.")
            pic_nm = data_path + "w%d_ra_dec.png"%area_id
            plt.savefig(pic_nm)
            plt.close()

            logger.info("Rank: %d start to write data to hdf5 %s" % (rank, cf_cata_data_path))
            f_row, f_col = final_data.shape

            h5f = h5py.File(cf_cata_data_path,"r+")

            h5f["/grid/w_%d"%area_id].attrs["max_block_size"] = numpy.array([num_in_block.max(),data_col],dtype=numpy.intc)

            h5f["/grid/w_%d/num_in_block" % area_id] = num_in_block
            h5f["/grid/w_%d/num_in_block" % area_id].attrs["shape"] = numpy.array([grid_num, 1], dtype=numpy.intc)

            h5f["/grid/w_%d/block_start" % area_id] = block_start
            h5f["/grid/w_%d/block_start" % area_id].attrs["shape"] = numpy.array([grid_num, 1], dtype=numpy.intc)

            h5f["/grid/w_%d/block_end" % area_id] = block_end
            h5f["/grid/w_%d/block_end" % area_id].attrs["shape"] = numpy.array([grid_num, 1], dtype=numpy.intc)

            esti_nm = ["/RA", "/DEC", "/G1", "/G2", "/N", "/U", "/V"]
            for ie, enm in enumerate(esti_nm):
                set_name = "/grid/w_%d/data/%s" % (area_id, enm)
                h5f[set_name] = final_data[:, ie]
                h5f[set_name].attrs["shape"] = numpy.array([f_row, 1], dtype=numpy.intc)
            h5f.close()

        logger.info("Rank: %d finish" % rank)
        comm.Barrier()
        t2 = time.time()
        if rank == 0:
            print("AREA: %d, "%area_id, t2-t1)








