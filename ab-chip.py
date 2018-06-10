import matplotlib
matplotlib.use('Agg')
import os
from sys import path
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/fourier_quad/'%my_home)
import tool_box
import numpy
from mpi4py import MPI
import time
from astropy.io import fits
import warnings
import logging

warnings.filterwarnings("error")

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()

with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        total_path = path.split("=")[1]
    elif "cfht_res_path" in path:
        result_path = path.split("=")[1]
    elif "cfht_pic_path" in path:
        pic_path = path.split("=")[1]
    elif "cfht_field_path" in path:
        field_path = path.split("=")[1]
log_path = "/home/hkli/work/logs/"
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logfile = log_path + '%d_log.dat' %rank

lf = logging.FileHandler(logfile, 'w')
form = logging.Formatter('%(asctime)s - %(message)s')
lf.setFormatter(form)
logger.addHandler(lf)

nname_path = total_path + "nname.dat"
field_dict, fields = tool_box.field_dict(nname_path)
r_fields = tool_box.allot(fields,cpus)[rank]

data_list = []
data_list_2 = []
step_list = []

for field in r_fields:
    expos = list(field_dict[field].keys())
    field_label = tool_box.cfht_label(field)
    for expo in expos:
        expo_label = int(expo.split("p")[0])
        chips = field_dict[field][expo]
        for chip in chips:
            chip_label = int(chip.split("_")[1].split(".")[0])
            chip_path = total_path + "%s/science/%s.fits"%(field, chip)

            f = fits.open(chip_path)
            img = f[0].data
            f.close()

            step = tool_box.find_step(img, 140)[1]
            if step == 1:
                step_list.extend([field_label, expo_label, chip_label, step])
            logger.info("%s--%s--%s: step: %d"%(field, expo, chip, step))
            # res_1 = tool_box.fit_backgroud(img, 2, 1, 10000)
            # slope_1 = [res_1[0][0, 0], res_1[1][0, 0]]
            # logger.info("%s--%s--%s: backgroud: %f, %f" % (field, expo, chip, slope_1[0],slope_1[1]))
            # slope_min_1, slope_max_1 = numpy.min(slope_1), numpy.max(slope_1)
            # a = (slope_max_1 - slope_min_1) / slope_max_1
            #
            # res_2 = tool_box.fit_backgroud(img, 1, 2, 10000)
            # slope_2 = [res_2[0][0, 0], res_2[1][0, 0]]
            # logger.info("%s--%s--%s: backgroud: %f, %f" % (field, expo, chip, slope_2[0], slope_2[1]))
            # slope_min_2, slope_max_2 = numpy.min(slope_2), numpy.max(slope_2)
            # b = (slope_max_2 - slope_min_2) / slope_max_2
            #
            # if a > 0.09 or b > 0.09:
            #     data_list.extend([field, expo_label, chip_label])
            #     data_list_2.extend([field_label, expo_label, chip_label])

# data = comm.gather(data_list, root=0)
# data_2 = comm.gather(data_list_2, root=0)
data_3 = comm.gather(step_list, root=0)

if rank == 0:

    # all_tags = []
    # for t in data_2:
    #     all_tags.extend(t)
    # data_tags = numpy.array(all_tags)
    # numpy.savez(result_path + "asterism.npz", data_tags)

    bi_path = result_path + "binary_label.npz"
    bi_data = numpy.load(bi_path)["arr_0"]
    f_lb = bi_data[:, 1]
    e_lb = bi_data[:, 2]
    c_lb = bi_data[:, 3]
    # stellar = numpy.zeros_like(c_lb)

    # labels = []
    # for lb in data:
    #     labels.extend(lb)
    # for tag, f in enumerate(labels):
    #     if "w" in f:
    #         field = tool_box.cfht_label(f)
    #         expo = labels[tag+1]
    #         chip = labels[tag+2]
    #         idx_f = f_lb == field
    #         idx_e = e_lb == expo
    #         idx_c = c_lb == chip
    #         try:
    #             stellar[idx_f&idx_e&idx_c] = 1
    #         except:
    #             print(field, expo, chip)
    steps = numpy.zeros_like(c_lb)
    all_steps = []
    for s in data_3:
        all_steps.extend(s)
    numpy.savez(result_path + "step_cache.npz", numpy.array(all_steps))
    for tag in range(0, len(all_steps), 4):
        field = all_steps[tag]
        expo = all_steps[tag+1]
        chip = all_steps[tag+2]
        step = all_steps[tag+3]
        idx_f = f_lb == field
        idx_e = e_lb == expo
        idx_c = c_lb == chip
        try:
            steps[idx_f&idx_e&idx_c] = step
        except:
            print(field, expo, chip)

    final_data = numpy.column_stack((bi_data,steps))
    numpy.savez(result_path+"final_ab_label.npz", final_data)







