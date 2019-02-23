import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/'%my_home)
import numpy
import time
import tool_box
import h5py
from mpi4py import MPI

# divide the data into many sub-sets for memory saving
# two cpus is needed

result_source = argv[1]

ts = time.clock()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()
if cpus > 2:
    print("Only two cpus is needed")
    exit()

fresh = "fresh_para_idx"
envs_path = "%s/work/envs/envs.dat" % my_home
get_contents = [['cfht', "cfht_path_result", '1'],['cfht', "cfht_path_data", '1'],
                ["%s"%fresh, "nstar", "1"],["%s"%fresh, "total_area", "1"],
                ["%s"%fresh, "flux2", '1'],["%s"%fresh, "flux_alt", '1'],
                ["%s"%fresh, "gf1", '1'], ["%s"%fresh, "gf2", '1'],
                ["%s"%fresh, "g1", '1'], ["%s"%fresh, "g2", '1'], ["%s"%fresh, "de", '1'],
                ["%s"%fresh, "h1", '1'], ["%s"%fresh, "h2", '1'], ]
path_items = tool_box.config(envs_path, ["get" for i in range(len(get_contents))], get_contents)
result_path, data_path = path_items[0], path_items[1]

# column labels
nstar_lb = int(path_items[2])
area_lb = int(path_items[3])
flux2_lb = int(path_items[4])
flux_alt_lb = int(path_items[5])
gf1_lb = int(path_items[6])
gf2_lb = int(path_items[7])

g1_lb = int(path_items[8])
g2_lb = int(path_items[9])
de_lb = int(path_items[10])
h1_lb = int(path_items[11])
h2_lb = int(path_items[12])


para_path = "./para.ini"
get_contents = [['field', "g1_num", '1'], ['field', "g2_num", '1'], ['field', "g1_s", '1'],
                ['field', "g1_e", '1'], ['field', "g2_s", '1'], ['field', "g2_e", '1']]
gets = ["get" for i in range(len(get_contents))]
path_items = tool_box.config(para_path, gets, get_contents)

# parameters for segment
g1_num = int(path_items[0])
g2_num = int(path_items[1])
g1_s = float(path_items[2])
g1_e = float(path_items[3])
g2_s = float(path_items[4])
g2_e = float(path_items[5])

# read the original catalog
cata_data_path = data_path + "cata_%s.hdf5"%result_source
h5f = h5py.File(cata_data_path, "r")
cata_data = h5f["/total"].value
h5f.close()

fg1 = numpy.linspace(g1_s, g1_e, g1_num)
fg2 = numpy.linspace(g2_s, g2_e, g2_num)
dfg1 = (fg1[1] - fg1[0])/2
dfg2 = (fg2[1] - fg2[0])/2
fgs = [fg1, fg2]
dfgs = [dfg1, dfg2]
gnums = [g1_num, g2_num]
gflbs = [gf1_lb, gf2_lb]

# nstar total_area flux2 flux_alt gf1 gf2 g1 g2 de h1 h2
# 0     1          2     3         4   5   6  7  8  9  10
ch = [nstar_lb, area_lb, flux2_lb, flux_alt_lb, gf1_lb, gf2_lb, g1_lb, g2_lb, de_lb, h1_lb, h2_lb]

# two files for g1 & g2 respectively
seg_data_path = data_path + "cata_%s_g%d.hdf5"%(result_source, rank+1)
h5f = h5py.File(seg_data_path, "w")
# two cpus
for i in range(gnums[rank]):
    t1 = time.time()
    idx1 = cata_data[:, gflbs[rank]] >= fgs[rank][i] - dfgs[rank]
    idx2 = cata_data[:, gflbs[rank]] < fgs[rank][i] + dfgs[rank]
    sub_data = cata_data[idx1&idx2][:, ch]
    set_name = "/fg%d_%d"%(rank+1, i)
    t2 = time.time()
    h5f[set_name] = sub_data
    t3 = time.time()
    print("Field g%d: %.6f ~ %.6f, galaxy: %d, t1: %.2f, t2: %.2f"
          %(rank+1, fgs[rank][i] - dfgs[rank], fgs[rank][i] + dfgs[rank], len(sub_data), t2-t1, t3-t2))
h5f.close()
te = time.clock()
print("TIME: %.2f"%(te-ts))




