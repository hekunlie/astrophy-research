from subprocess import Popen


# comm = ['mpirun -np 70 python2 selection_bias_images_h.py','mpiexec -n 70 /home/hkli/work/c/measure',
#         'mpirun -np 14 python selection_bias_proce_dat_h.py flux 0 1000000 0 sym 5',
#         'mpirun -np 14 python sym_mc_plot_h_1.py peak', 'mpirun -np 14 python sym_mc_plot_h_1.py flux',
#         'mpirun -np 14 python sym_mc_plot_h_1.py fsnr_f', 'mpirun -np 14 python sym_mc_plot_h_1.py fsnr_c',
#         'mpirun -np 14 python sym_mc_plot_h_1.py sesnr', 'mpirun -np 14 python sym_mc_plot_h_1.py osnr',
#         'mpirun -np 14 python sym_mc_plot_h_1.py snr']


sex_cut = [20, 30, 45, 60]
for j in range(4):
    for i in range(8):
        cmd = 'mpirun -np 14 python selection_bias_proce_dat_h.py sex %f 1000000 0 sym 5 %d'%(sex_cut[j], i)
        a = Popen(cmd, shell=True)
        a.wait()

snr_cut = [16., 25., 40., 70.]
for j in range(4):
    for i in range(8):
        cmd = 'mpirun -np 14 python selection_bias_proce_dat_h.py snr %f 1000000 0 sym 5 %d'%(snr_cut[j], i)
        a = Popen(cmd, shell=True)
        a.wait()
fsnr_cut = [3.5, 6., 10., 15.]
for j in range(4):
    for i in range(8):
        cmd = 'mpirun -np 14 python selection_bias_proce_dat_h.py fsnr %f 1000000 0 sym 5 %d'%(fsnr_cut[j], i)
        a = Popen(cmd, shell=True)
        a.wait()