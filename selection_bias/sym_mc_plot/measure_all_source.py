import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append('D:/GitHub/astrophy-research/mylib')
import numpy
import h5py
from Fourier_Quad import Fourier_Quad
from mpi4py import MPI
import tool_box
from plot_tool import Image_Plot


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# calculate the shear of the whole sample
source_nm = argv[1]
total_path = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s/"%source_nm

shear = numpy.load(total_path + "parameters/shear.npz")
try:
    g1 = shear["arr_0"][:,0]
    g2 = shear["arr_1"][:,0]
except:
    g1 = shear["arr_0"]
    g2 = shear["arr_1"]


itemsize = MPI.DOUBLE.Get_size()
element_num = g1.shape[0]
if rank == 0:
    # bytes for 10 double elements
    nbytes = element_num*itemsize*4
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
result = numpy.ndarray(buffer=buf1, dtype='d', shape=(4, element_num)) # array filled with zero

# read data
data_path = total_path + "result/data/data_%d.hdf5"%rank
h5f = h5py.File(data_path, "r")
data = h5f["/data"].value
h5f.close()

mg1 = data[:,2]
mg2 = data[:,3]
mnu1 = data[:,4] + data[:,5]
mnu2 = data[:,4] - data[:,5]

fq = Fourier_Quad(10,123)
gh1, gh1_sig = fq.find_shear(mg1, mnu1, 8)[:2]
gh2, gh2_sig = fq.find_shear(mg2, mnu2, 8)[:2]

result[0, rank] = gh1
result[1, rank] = gh1_sig
result[2, rank] = gh2
result[3, rank] = gh2_sig

comm.Barrier()
if rank == 0:
    mc1 = numpy.array(tool_box.data_fit(g1,result[0], result[1]))
    mc2 = numpy.array(tool_box.data_fit(g2,result[2], result[3]))

    mc1[0] = mc1[0] - 1
    mc2[0] = mc2[0] - 1

    result_path = total_path + "result/data/shear_result.hdf5"
    h5f = h5py.File(result_path,"w")
    h5f["/mc1"] = mc1
    h5f["/mc2"] = mc2
    h5f["/shear"] = result
    h5f.close()

    print(mc1)
    print(mc2)
    
    img = Image_Plot()
    img.subplots(1,1)
    img.axs[0][0].errorbar(g1, result[0], result[1],label="g1: $10^2$m1:%.2f(%.2f),$10^4$c1:%.2f(%.2f)"
                                                         %(mc1[0]*100,mc1[1]*100,mc1[2]*10000,mc1[3]*10000))
    img.axs[0][0].errorbar(g2, result[2], result[3],label="$g2: 10^2$m1:%.2f(%.2f),$10^4$c1:%.2f(%.2f)"
                                                         %(mc2[0]*100,mc2[1]*100,mc2[2]*10000,mc2[3]*10000))
    img.axs[0][0].legend()
    img.save_img("%s.png"%source_nm)


