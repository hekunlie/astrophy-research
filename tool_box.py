from __future__ import division
import matplotlib
matplotlib.use('Agg')
from multiprocessing import Pool, Manager
import numpy
import copy
import matplotlib.pyplot as plt

def task_distri(target_list, cpu_num):
    # it will divide the target_list into some piece (small lists in a diction)
    # of which the number depends on the specific cpu core number
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
        data = numpy.load(files_paths[i])["arr_0"]
        # data = pandas.read_excel(files_paths[i]).values
        if i == 0:
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
    for i in range(cpu_num):
        p.apply_async(list_add, args=(final_data_list, paths_distri[i],))
    p.close()
    p.join()
    return final_data_list

def detect(mask, ini_y, ini_x, signal, signal_val, y_size, x_size):
    if mask[ini_y, ini_x] > 0:
        signal.append((ini_y, ini_x))
        signal_val.append(mask[ini_y, ini_x])
        mask[ini_y, ini_x] = 0
        for cor in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            if -1 < ini_y + cor[0] < y_size and -1 < ini_x + cor[1] < x_size and mask[ini_y + cor[0], ini_x + cor[1]] > 0:
                detect(mask, ini_y + cor[0], ini_x + cor[1], signal, signal_val, y_size, x_size)
    return signal, signal_val

def stamp_detector(image, thres, y_size, x_size, ra=10):
    # get the source object
    image_c = copy.copy(image)
    img_idx = image_c < thres
    image_c[img_idx] = 0.

    center = image_c[int(y_size/2-ra):int(y_size/2+ra), int(x_size/2-ra):int(x_size/2+ra)]
    maxi = numpy.max(center)
    y, x = numpy.where(center > 0)
    y_sour = y + int(y_size/2 - ra)
    x_sour = x + int(x_size/2 - ra)

    final_obj = []
    final_flux = []
    for m in range(len(x_sour)):
        sour_pool = []
        sour_flux = []
        yp = y_sour[m]
        xp = x_sour[m]

        if image_c[yp, xp] > 0:
            sour_pool, sour_flux = detect(image_c, yp, xp, sour_pool, sour_flux, y_size, x_size)

        if len(sour_pool) > len(final_obj):
            final_obj = sour_pool
            final_flux = sour_flux
        elif len(sour_pool) == len(final_obj):
            if numpy.max(final_flux) < numpy.max(sour_flux):
                final_obj = sour_pool
                final_flux = sour_flux

    return final_obj, numpy.sum(final_flux), numpy.sum((numpy.array(final_flux))**2), numpy.max(final_flux)

def source_detector(mask, ysize, xsize):
    # get the source object
    objects = []
    p = numpy.where(mask > 0)
    xp, yp = p[0], p[1]
    for j in range(len(xp)):
        if mask[xp[j], yp[j]] > 0:
            cache = [(xp[j], yp[j])]
            mask[xp[j], yp[j]] = 0
            num = 2
            num0 = 1
            while True:
                num_new = num0 - num
                if num == num0:
                    break
                num0 = len(cache)
                p_new = []
                for k in range(num_new, 0):
                    xy = cache[k]
                    for coord in ((xy[0] + 1, xy[1]), (xy[0] - 1, xy[1]), (xy[0], xy[1] - 1), (xy[0], xy[1] + 1)):
                        if -1 < coord[0] < ysize and -1 < coord[1] < xsize and mask[coord[0], coord[1]] > 0:
                            p_new.append(coord)
                            mask[coord[0], coord[1]] = 0
                cache.extend(p_new)
                num = len(cache)
            if num > 5:
                objects.append(cache)
    return objects

def get_hlr(image, scale, size,):
    # get the source object, to avoid the overflow of the stack
    mask = copy.copy(image)
    maxi = numpy.max(image[int(size/2-6):int(size/2+6), int(size/2-6):int(size/2+6)])
    y, x = numpy.where(mask == maxi)
    idx = mask < maxi / scale
    mask[idx] = 0.
    flux = maxi
    cache = [(x, y)]
    mask[y, x] = 0
    num = 2
    num0 = 1
    while True:
        num_new = num0 - num
        if num == num0:
            break
        num0 = len(cache)
        p_new = []
        for k in range(num_new, 0):
            xy = cache[k]
            for coord in ((xy[0] + 1, xy[1]), (xy[0] - 1, xy[1]), (xy[0], xy[1] - 1), (xy[0], xy[1] + 1)):
                if -1 < coord[0] < size and -1 < coord[1] < size and mask[coord[0], coord[1]] > 0:
                    p_new.append(coord)
                    flux += mask[coord[0], coord[1]]
                    mask[coord[0], coord[1]] = 0
        cache.extend(p_new)
        num = len(cache)

    return numpy.sqrt(len(cache)/numpy.pi), flux

def data_fit(x_data, y_data, y_err):
    # Y = A*X ,   y = m*x+c
    # Y = [y1,y2,y3,...].T  the measured data
    # A = [[1,1,1,1,...]
    #         [x1,x2,x3..]].T
    # X = [c,m].T
    # C = diag[sig1^2, sig2^2, sig3^2, .....]
    # the inverse of C is used as weight of data
    # X = [A.T*C^(-1)*A]^(-1) * [A.T*C^(-1) *Y]

    A1 = numpy.column_stack((numpy.ones_like(x_data.T), x_data.T))
    Y1 = y_data.T
    C1 = numpy.diag(y_err ** 2)
    L1 = numpy.linalg.inv(numpy.dot(numpy.dot(A1.T, numpy.linalg.inv(C1)), A1))
    R1 = numpy.dot(numpy.dot(A1.T, numpy.linalg.inv(C1)), Y1)
    sig_m1 = numpy.sqrt(L1[1, 1])
    sig_c1 = numpy.sqrt(L1[0, 0])
    mc = numpy.dot(L1, R1)

    # A1 = numpy.column_stack((numpy.ones_like(fg1.T), fg1.T))
    # Y1 = res_arr1[i].T
    # C1 = numpy.diag(res_arr1[i+4]**2)
    #
    # A2 = numpy.column_stack((numpy.ones_like(fg2.T), fg2.T))
    # Y2 = res_arr2[i].T
    # C2 = numpy.diag(res_arr2[i+4]**2)
    #
    # L1 = numpy.linalg.inv(numpy.dot(numpy.dot(A1.T, numpy.linalg.inv(C1)), A1))
    # R1 = numpy.dot(numpy.dot(A1.T, numpy.linalg.inv(C1)), Y1)
    # L2 = numpy.linalg.inv(numpy.dot(numpy.dot(A2.T, numpy.linalg.inv(C2)), A2))
    # R2 = numpy.dot(numpy.dot(A2.T, numpy.linalg.inv(C2)), Y2)
    #
    # sig_m1 = numpy.sqrt(L1[1, 1])
    # sig_c1 = numpy.sqrt(L1[0, 0])
    # sig_m2 = numpy.sqrt(L2[1, 1])
    # sig_c2 = numpy.sqrt(L2[0, 0])
    # e1mc = numpy.dot(L1, R1)
    # e2mc = numpy.dot(L2, R2)
    return mc[1], sig_m1, mc[0], sig_c1

def mcplot(x1_data, y1_data, x2_data, y2_data, e1mc, e2mc, cut_start, cut_end, path=None):
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(121)
    ax.errorbar(x1_data, y1_data[0], y1_data[1], ecolor='black', elinewidth='1', fmt='none', capsize=2)
    ax.plot(x1_data, e1mc[0] * x1_data + e1mc[2], color='red')
    ax.plot([-0.1, 0.1], [-0.1, 0.1], label='y=x', color='blue')
    ax.scatter(x1_data, y1_data[0], c='black')
    for j in range(len(x1_data)):
        ax.text(x1_data[j], y1_data[0,j], str(round(y1_data[2,j] / 1000, 1)) + "K", color="red")
    ax.text(0.1, 0.85, 'm=' + str(round(e1mc[0] - 1, 6)) + '$\pm$' + str(round(e1mc[1], 6)), color='green', ha='left',
            va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.8, 'c=' + str(round(e1mc[2], 6)) + '$\pm$' + str(round(e1mc[3], 6)), color='green', ha='left',
            va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.75, "[ " + cut_start + ", " + cut_end + "]", color='green', ha='left', va='center', transform=ax.transAxes,
            fontsize=20)
    plt.xlabel('True  g1', fontsize=20)
    plt.ylabel('Est  g1', fontsize=20)
    plt.legend(fontsize=15)
    plt.ylim(-0.04, 0.04)
    plt.xlim(-0.04, 0.04)

    # plot g2 line
    ax = fig.add_subplot(122)
    ax.errorbar(x2_data, y2_data[0], y2_data[1], ecolor='black', elinewidth='1', fmt='none', capsize=2)
    ax.plot(x2_data, e2mc[0] * x2_data + e2mc[2], color='red')
    ax.plot([-0.1, 0.1], [-0.1, 0.1], label='y=x', color='blue')
    ax.scatter(x2_data, y2_data[0], c='black')
    for j in range(len(x2_data)):
        ax.text(x2_data[j], y2_data[0, j], str(round(y2_data[2, j] / 1000, 1)) + "K", color="red")
    ax.text(0.1, 0.85, 'm=' + str(round(e2mc[0] - 1, 6)) + '$\pm$' + str(round(e2mc[1], 6)), color='green', ha='left',
            va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.8, 'c=' + str(round(e2mc[2], 6)) + '$\pm$' + str(round(e2mc[3], 6)), color='green', ha='left',
            va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.75, "[ " + cut_start + ", " + cut_end + "]", color='green', ha='left', va='center', transform=ax.transAxes,
            fontsize=20)
    plt.xlabel('True  g2', fontsize=20)
    plt.ylabel('Est  g2', fontsize=20)
    plt.legend(fontsize=15)
    plt.ylim(-0.04, 0.04)
    plt.xlim(-0.04, 0.04)
    if path is not None:
        plt.savefig(path)
    plt.close()