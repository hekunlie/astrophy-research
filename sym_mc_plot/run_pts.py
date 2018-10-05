from subprocess import Popen

file_num = 5
comm = ['mpirun -np 14 python sym_mc_plot_pts.py sex_snr %d'%file_num,
        'mpirun -np 14 python sym_mc_plot_pts.py flux2 %d'%file_num,
        'mpirun -np 14 python sym_mc_plot_pts.py mag_auto %d'%file_num]


for cmd in comm:
    a = Popen(cmd, shell=True)
    a.wait()

