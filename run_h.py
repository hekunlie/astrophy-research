from subprocess import Popen
from sys import argv

# if int(argv[1]) == 1:
comm = ['mpirun -np 70 python2 selection_bias_images_h.py', 'mpiexec -n 70 /home/hkli/work/c/measure_m2',
        'mpirun -np 14 python selection_bias_proce_dat_h.py flux 0 10000000 0 sym 5']
# else:
#     comm = ['mpirun -np 14 python sym_mc_plot_h.py fsnr_c', 'mpirun -np 14 python sym_mc_plot_h.py fsnr_c_m',
#             'mpirun -np 14 python sym_mc_plot_h.py fsnr_f','mpirun -np 14 python sym_mc_plot_h.py fsnr_f_m']
# comm = ['mpirun -np 14 python sym_mc_plot_h.py peak', 'mpirun -np 14 python sym_mc_plot_h.py flux2',
#         'mpirun -np 14 python sym_mc_plot_h.py flux_alt','mpirun -np 14 python sym_mc_plot_h.py area',
#         'mpirun -np 14 python sym_mc_plot_h.py harea', 'mpirun -np 14 python sym_mc_plot_h.py flux',
#         'mpirun -np 14 python sym_mc_plot_h.py hflux', 'mpirun -np 14 python sym_mc_plot_h.py sex',
#         'mpirun -np 14 python sym_mc_plot_h.py snr']

for cmd in comm:
    a = Popen(cmd, shell=True)
    a.wait()


# for i in range(19):
#     cmd = 'mpirun -np 25 python process_data.py %d 20 peak 10'%i
#     a = Popen(cmd, shell=True)
#     a.wait()

# for i in range(19):
#     cmd = 'mpirun -np 35 python process_data.py %d 20 peak 10'%i
#     a = Popen(cmd, shell=True)
#     a.wait()
#
#
# for i in range(19):
#     cmd = 'mpirun -np 25 python process_data.py %d 20 fsnr 5'%i
#     a = Popen(cmd, shell=True)
#     a.wait()
#
# for i in range(19):
#     cmd = 'mpirun -np 35 python process_data.py %d 20 fsnr 5'%i
#     a = Popen(cmd, shell=True)
#     a.wait()