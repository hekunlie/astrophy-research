from subprocess import Popen
from sys import argv

nargv = len(argv)
process = []
for i in range(nargv):
    if 0 < i < 4:
        c = int(float(argv[i]))
        process.append(c)

bin_num = 8
ch = "fsnr"
ch_thres = 0
g2num = 21
for i in process:
    if i == 1:
        print("stack data")
        cmd = "mpirun -np 40 python cfht_stack.py"
        a = Popen(cmd, shell=True)
        a.wait()
    if i==2:
        print("fit g")
        cmd = "mpirun -np %d python process_data.py 0 %d %s %.2f"%(g2num, bin_num, ch, ch_thres)
        a = Popen(cmd, shell=True)
        a.wait()
    if i==3:
        print("filter data")
        cache = "0_%s_%.2f_final_cache.npz"%(ch,ch_thres)
        cmd = "python cfht_data_filter.py %s %d %d %d %.2f"%(cache, g2num-4, g2num, bin_num, ch_thres)
        a = Popen(cmd, shell=True)
        a.wait()
