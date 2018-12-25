import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import tool_box
import h5py
from mpi4py import MPI
from sys import argv
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

grid_scale = numpy.array([5., 10., 25., 60.]) # arcsec
corre_scale = [10**(0.2*i) for i in range(11)]

cpu_block_size = int(cpus/area_num)



envs_path = "%s/work/envs/envs.dat"%my_home
gets_item = [["cfht","%s_path_para","0"]]
para_path = tool_box.config(envs_path,["get"],[])[0]

!!!
data_path = config.get("cfht", "data_path_in")
res_path = config.get("cfht", "result_path")
flux_alt = int(config.get("fresh_para_idx", "flux_alt"))
nstar = int(config.get("fresh_para_idx", "nstar"))
flux_alt_thresh = 3.
nstar_thresh = 12

# cpus_0 read the fields catalog and distribute them to related cpus
# cpu_0 & block_0, cpu_1 & block_1 ...
# cpu_i & block_i will be in charge of area_i

# data collection
if cmd == "collect":
    dicts, fields = tool_box.field_dict(data_path+"nname.dat")
    for i in range(area_num):
        if rank == 0:
            field_tar = []
            for field in fields:
                if "w%d"%(i+1) in field:
                    field_tar.append(field)
            field_pool = tool_box.allot(field_tar, cpus)
        else:
            field_pool = None

        field_pool = comm.scatter(field_pool, root=0)

        field_count = 0
        for field_name in field_pool:
            field_cat_path = data_path + "%s/result/%s_shear.dat"%(field_name, field_name)
            if os.path.exists(field_cat_path):
                try:
                    cat_arr = numpy.loadtxt(field_cat_path)
                    if field_count == 0:
                        cat_data = cat_arr
                    else:
                        cat_data = numpy.row_stack((cat_data, cat_arr))
                    field_count += 1
                except:
                    print(rank, "%s.dat doesn't exist"%field_name)
        data_sp = cat_data.shape
        data_sps = comm.gather(data_sp, root=0)

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
        comm.Gatherv(cat_data, [recv_buffer, count, displ, MPI.DOUBLE], root=0)
        if rank == 0:
            final_data_path = res_path + "w%d.npz"%i
            numpy.savez(final_data_path, recv_buffer)

# make grid
if cmd == "grid":
    if rank < area_num:
        t1 = time.time()
        # test
        grid_h5_path = data_path + "w%d_grid.hdf5" % rank
        cat_path = data_path + "w%d.npz"%rank
        ori_cat_data = numpy.load(cat_path)["arr_0"]

        # cut off
        flux_alt_idx = ori_cat_data[:,flux_alt] >= flux_alt_thresh
        nstar_idx = ori_cat_data[:,nstar] >= nstar_thresh
        cat_data = ori_cat_data[flux_alt_idx&nstar_idx][:,12:21]

        # convert to arcmin
        ra, dec = cat_data[:,0]*60, cat_data[:,1]*60
        ra_min, ra_max, dec_min, dec_max = ra.min(), ra.max(), dec.min(), dec.max()

        plt.scatter(dec, ra)
        pic_nm = data_path + "w%d_ra_dec.png" %rank
        plt.savefig(pic_nm)
        plt.close()

        # the grid indices are labeled by the group name under the specific grid scale,
        # e.g. /0.5/1/2/block (/scale/col/row/block),
        # the attribute of block shows the shape of the data array in the grid
        f = h5py.File(grid_h5_path, "w")

        for i, scale in enumerate(grid_scale):
            # the grid: x-axis--ra, y-axis--dec
            rows, cols = int((dec_max - dec_min)/scale+1), int((ra_max - ra_min)/scale+1)
            test_arr = numpy.zeros((5 * rows, 5 * cols))
            for col in range(cols):
                for row in range(rows):

                    ir1 = ra >= ra_min + row*scale
                    ir2 = ra < ra_min + (row+1)*scale
                    ic1 = dec >= dec_min + col*scale
                    ic2 = dec < dec_min + (col+1)*scale
                    block = cat_data[ir1&ir2&ic1&ic2]
                    block_sp = block.shape

                    if block_sp[0] > 0:
                        test_arr[row*5:(row+1)*5, col*5:(col+1)*5] = block_sp[0]

                    group_name = "/%d/%d/%d"%(scale, col, row)
                    f.create_group(group_name)
                    data_set_name = group_name+"/block"
                    f.create_dataset(name=data_set_name, data=block)
                    f[data_set_name].attrs["info"] = [block_sp[0], block_sp[1],
                                                      ra_min + row*scale, ra_min + (row+1)*scale,
                                                      dec_min + col*scale, dec_min + (col+1)*scale]

            plt.imshow(test_arr)
            plt.colorbar()
            pic_nm = data_path + "w%d_grid_%.f.png" %(rank,scale)
            plt.savefig(pic_nm)
            plt.close()
        f.close()
        t2 = time.time()
        if rank == 0:
            print(t2-t1)









