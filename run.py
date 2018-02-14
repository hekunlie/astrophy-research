from subprocess import Popen
from sys import argv

if int(argv[1]) == 1:
    comm = ['mpiexec -n 14 /home/hklee/work/c/simu',
            'mpirun -np 14 python selection_bias_proce_dat.py flux 0 1000000 0 sym']
else:
    comm = ['mpirun -np 14 python sym_mc_plot.py peak', 'mpirun -np 14 python sym_mc_plot.py flux',
            'mpirun -np 14 python sym_mc_plot.py fsnr', 'mpirun -np 14 python sym_mc_plot.py fsnr1',
            'mpirun -np 14 python sym_mc_plot.py fsnr4', 'mpirun -np 14 python sym_mc_plot.py fsnr9']#,
           # 'mpirun -np 14 python sym_mc_plot.py osnr', 'mpirun -np 14 python sym_mc_plot.py snr']
for cmd in comm:
    a = Popen(cmd, shell=True)
    a.wait()