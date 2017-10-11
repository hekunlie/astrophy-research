import Fourier_Quad
import pandas
from multiprocessing import Pool, Manager
import numpy
import time
import copy

def task_distri(target_list, cpu_num):
    # it will divide the target_list into some piece (small lists in a diction) of which the number depends on the specific cpu core number
    # target_list must be a python list
    m, n = divmod(len(target_list), cpu_num)
    distri_pool = {}
    if m == 0 and n != 0:
        for i in range(n):
            distri_pool[i] = map(int, str(target_list[i]))
    elif m != 0:
        for i in range(cpu_num):
            distri_pool[i] = target_list[i * m:(i + 1) * m]
        if n != 0:
            for i in range(n):
                distri_pool[i].append(target_list[-i - 1])
    else:
        print("Caution! Something goes wrong!!!")
    return distri_pool

def list_add(target_list, files_paths):
    # this function is designed for 'function 'classify'
    # target_list is the target list that this function will put data array into
    # files_paths is a list of paths of excel files
    for i in range(len(files_paths)):
        data = pandas.read_excel(files_paths[i]).values
        if i==0:
            temp_data = data
        else:
            temp_data = numpy.row_stack((temp_data, data))
    target_list.append(temp_data)
    target_list.reverse()

def classify(files_path_list, cpu_num):
    # the data will be assembled into same big arrays of which the number equals to cpu_num
    paths_distri = task_distri(files_path_list, cpu_num)
    final_data_list = Manager().list()
    p = Pool()
    cl_ts = time.clock()
    for i in range(cpu_num):
        p.apply_async(list_add, args=(final_data_list, paths_distri[i],))
    p.close()
    p.join()
    cl_te = time.clock()
    return final_data_list,cl_te-cl_ts

def detect(mask, ini_y, ini_x, signal, signal_val, img_size):
    if mask[ini_y, ini_x] > 0:
        signal.append((ini_y, ini_x))
        signal_val.append(mask[ini_y, ini_x])
        mask[ini_y, ini_x] = 0
        for cor in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            if ini_y + cor[0] < img_size and ini_x + cor[1] < img_size and mask[ini_y + cor[0], ini_x + cor[1]] > 0:
                detect(mask, ini_y + cor[0], ini_x + cor[1], signal, signal_val, img_size)
    return signal, signal_val

def stamp_detector(image, thres, size, ra=8, obj='stamp'):
    # get the source object
    image_c = copy.copy(image)
    img_idx = image_c < thres
    image_c[img_idx] = 0.

    center = image_c[int(size/2-ra):int(size/2+ra), int(size/2-ra):int(size/2+ra)]
    maxi = numpy.max(center)
    y, x = numpy.where(center > 0)
    y_sour = y + int(size/2-ra)
    x_sour = x + int(size / 2 - ra)

    final_obj = []
    final_flux = []
    for m in range(len(x_sour)):
        sour_pool = []
        sour_flux = []
        yp = y_sour[m]
        xp = x_sour[m]

        if image_c[yp, xp] > 0:
            sour_pool, flux = detect(image_c, yp, xp, sour_pool, sour_flux, size)

        if obj == 'stamp':
            if len(sour_pool) > len(final_obj):
                final_obj = sour_pool
                final_flux = sour_flux
        else:
            if len(sour_pool) > 5:
                final_obj.append(sour_pool)
                final_flux.append(sour_flux)

    return numpy.sqrt(len(final_obj) / numpy.pi), final_obj, numpy.sum(final_flux), maxi

