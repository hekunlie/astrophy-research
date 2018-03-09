from subprocess import Popen
from sys import argv

if int(argv[1]) == 1:
    comm = ['mpirun -np 14 python selection_bias_images.py','mpiexec -n 14 /home/hklee/work/c/measure',
            'mpirun -np 14 python selection_bias_proce_dat.py flux 0 1000000 0 sym']
else:
    comm = ['mpirun -np 14 python sym_mc_plot.py peak', 'mpirun -np 14 python sym_mc_plot.py flux',
            'mpirun -np 14 python sym_mc_plot.py fsnr_f', 'mpirun -np 14 python sym_mc_plot.py fsnr_c',
            'mpirun -np 14 python sym_mc_plot.py sesnr', 'mpirun -np 14 python sym_mc_plot.py osnr',
            'mpirun -np 14 python sym_mc_plot.py snr']
for cmd in comm:
    a = Popen(cmd, shell=True)
    a.wait()
