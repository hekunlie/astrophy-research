from subprocess import Popen
from sys import argv

# if int(argv[1]) == 1:
#     comm = ['mpirun -np 14 python selection_bias_images.py','mpiexec -n 14 /home/hklee/work/c/measure',
#             'mpirun -np 14 python selection_bias_proce_dat.py flux 0 1000000 0 sym']
# else:
comm = ['mpirun -np 80 python sym_mc_plot_cfht.py flux2', 'mpirun -np 80 python sym_mc_plot_cfht.py flux_alt']
for cmd in comm:
    a = Popen(cmd, shell=True)
    a.wait()
