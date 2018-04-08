from subprocess import Popen
from sys import argv

if int(argv[1]) == 1:
    comm = ['mpirun -np 70 python2 selection_bias_images_h.py','mpiexec -n 70 /home/hkli/work/c/measure',
            'mpirun -np 14 python selection_bias_proce_dat_h.py flux 0 1000000 0 sym 5']
else:
    comm = ['mpirun -np 14 python sym_mc_plot_h.py fsnr_c', 'mpirun -np 14 python sym_mc_plot_h.py fsnr_c_m',
            'mpirun -np 14 python sym_mc_plot_h.py fsnr_f','mpirun -np 14 python sym_mc_plot_h.py fsnr_f_m']
    # comm = ['mpirun -np 14 python sym_mc_plot_h.py peak', 'mpirun -np 14 python sym_mc_plot_h.py flux',
    #         'mpirun -np 14 python sym_mc_plot_h.py sex25','mpirun -np 14 python sym_mc_plot_h.py fsnr_c',
    #         'mpirun -np 14 python sym_mc_plot_h.py sex37', 'mpirun -np 14 python sym_mc_plot_h.py osnr',
    #         'mpirun -np 14 python sym_mc_plot_h.py snr', 'mpirun -np 14 python sym_mc_plot_h.py sex59']

for cmd in comm:
    a = Popen(cmd, shell=True)
    a.wait()
