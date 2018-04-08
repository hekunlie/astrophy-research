from subprocess import Popen


# comm = ['mpirun -np 70 python2 selection_bias_images_h.py','mpiexec -n 70 /home/hkli/work/c/measure',
#         'mpirun -np 14 python selection_bias_proce_dat_h.py flux 0 1000000 0 sym 5',
#         'mpirun -np 14 python sym_mc_plot_h_1.py peak', 'mpirun -np 14 python sym_mc_plot_h_1.py flux',
#         'mpirun -np 14 python sym_mc_plot_h_1.py fsnr_f', 'mpirun -np 14 python sym_mc_plot_h_1.py fsnr_c',
#         'mpirun -np 14 python sym_mc_plot_h_1.py sesnr', 'mpirun -np 14 python sym_mc_plot_h_1.py osnr',
#         'mpirun -np 14 python sym_mc_plot_h_1.py snr']
comm = ['mpirun -np 14 python sym_mc_plot_h_1.py sex25', 'mpirun -np 14 python sym_mc_plot_h_1.py sex59',
        'mpirun -np 14 python sym_mc_plot_h_1.py sex37',]
for cmd in comm:
    a = Popen(cmd, shell=True)
    a.wait()