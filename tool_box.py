import matplotlib
matplotlib.use('Agg')
from multiprocessing import Pool, Manager
import numpy
import copy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

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
        data = numpy.loadtxt(files_paths[i])
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

def detect(mask, ini_y, ini_x, signal, signal_val, y_size, x_size, sflag):
    if mask[ini_y, ini_x] > 0:
        signal.append((ini_y, ini_x))
        signal_val.append(mask[ini_y, ini_x])
        mask[ini_y, ini_x] = 0
        for cor in ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)):
            if -1 < ini_y + cor[0] < y_size and -1 < ini_x + cor[1] < x_size and mask[ini_y + cor[0], ini_x + cor[1]]>0:
                if ini_y + cor[0] == 0 or ini_y + cor[0] == y_size-1 or ini_x + cor[1] == 0 or ini_x + cor[1]==x_size-1:
                    sflag = 1
                detect(mask, ini_y + cor[0], ini_x + cor[1], signal, signal_val, y_size, x_size, sflag)
    return signal, signal_val, sflag

def stamp_detector(image, thres, y_size, x_size, ra=10):
    # get the source object
    image_c = copy.copy(image)
    img_idx = image_c < thres
    image_c[img_idx] = 0.

    center = image_c[int(y_size/2-ra):int(y_size/2+ra), int(x_size/2-ra):int(x_size/2+ra)]
    y, x = numpy.where(center > 0)
    y_sour = y + int(y_size/2 - ra)
    x_sour = x + int(x_size/2 - ra)

    final_obj = []
    final_flux = []
    flag = 0
    for m in range(len(x_sour)):
        sour_pool = []
        sour_flux = []
        yp = y_sour[m]
        xp = x_sour[m]

        if image_c[yp, xp] > 0:
            sour_pool, sour_flux, flag = detect(image_c, yp, xp, sour_pool, sour_flux, y_size, x_size, flag)

        if len(sour_pool) > len(final_obj):
            final_obj = sour_pool
            final_flux = sour_flux
        elif len(sour_pool) == len(final_obj):
            if numpy.max(final_flux) < numpy.max(sour_flux):
                final_obj = sour_pool
                final_flux = sour_flux
    if len(final_flux) == 0:
        peak = 0
    else:
        peak = numpy.max(final_flux)
    return final_obj, numpy.sum(final_flux), numpy.sum((numpy.array(final_flux))**2), peak, flag

def source_detector(mask, ysize, xsize):
    # get the source object
    objects = []
    p = numpy.where(mask > 0)
    yp, xp = p[0], p[1]
    for j in range(len(xp)):
        if mask[yp[j], xp[j]] > 0:
            cache = [(yp[j], xp[j])]
            mask[yp[j], xp[j]] = 0
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

def get_hlr(image, scale, size, thres=None):
    # get the source object, to avoid the overflow of the stack
    mask = copy.copy(image)
    maxi = numpy.max(image[int(size/2-6):int(size/2+6), int(size/2-6):int(size/2+6)])
    y, x = numpy.where(mask == maxi)
    if thres:
        idx = mask < thres
        mask[idx] = 0.
    else:
        idx = mask < maxi / scale
        mask[idx] = 0.
    flux = maxi
    flux_sq = maxi**2

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
                    flux_sq += mask[coord[0], coord[1]] ** 2

                    mask[coord[0], coord[1]] = 0
        cache.extend(p_new)
        num = len(cache)

    return numpy.sqrt(len(cache)/numpy.pi), flux, flux_sq, num, maxi


def gauss_fit(data, bin_num):
    # fit the Gaussian distribution (amplitude, sigma, mu)
    def fun(x, ampli, sig, mu):
        return ampli * numpy.exp(-(x - mu) ** 2 / 2 / sig ** 2)

    num, bins = numpy.histogram(data, bin_num)
    num = num/numpy.sum(num)
    mp = numpy.where(num == numpy.max(num))[0][0]
    # if mp+1 < bin_num/3:
    #     bins = bins[:mp*2]
    #     num = num[:mp*2]
    # elif bin_num/3 <= mp+1 < bin_num*2./3:
    #     bins = bins[int(mp+1-bin_num/3):int(mp+1-bin_num/3)]
    #     num = num[int(mp+1-bin_num/3):int(mp+1-bin_num/3)]
    # else:
    #     tag = bin_num - mp - 1
    #     bins = bins[mp-tag:mp+tag]
    #     num = num[mp-tag:mp+tag]
    # print(len(num))

    coeff, coerr = curve_fit(f=fun, xdata=bins[1:], ydata=num)
    # the fitted sigma can be negative
    return coeff, coerr, bins, num


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
    return mc[1], sig_m1, mc[0], sig_c1

def mcplot(x1_data, y1_data, x2_data, y2_data, e1mc, e2mc, cut_start, cut_end, xylim, path=None):
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(121)
    ax.errorbar(x1_data, y1_data[0], y1_data[1], ecolor='black', elinewidth=1, fmt='none', capsize=2)
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
    plt.ylim(xylim[0], xylim[1])
    plt.xlim(xylim[0], xylim[1])

    # plot g2 line
    ax = fig.add_subplot(122)
    ax.errorbar(x2_data, y2_data[0], y2_data[1], ecolor='black', elinewidth=1, fmt='none', capsize=2)
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
    plt.ylim(xylim[2], xylim[3])
    plt.xlim(xylim[2], xylim[3])
    plt.show()
    if path is not None:
        plt.savefig(path)
    plt.close()

def mags_mock(num, mag_min, mag_max):
    m = numpy.linspace(mag_min, mag_max, 1000000)
    pm = 10**(23.04187527*numpy.log10(m) - 32.50618926)
    pm = pm/numpy.sum(pm)
    new_pdf = numpy.random.choice(m, num, p=pm)
    return new_pdf


def ellip_mock(num, seed=123400, figout=None):
    """
    Generate random ellipticity for a given number
    following the distribution of the ellipticity
    of bulge-dominated galaxies

    See Miller et al, 2013, MNRAS
    """
    numpy.random.RandomState(seed)
    b, c = 2.368, 6.691
    # probability
    pe = lambda e: 27.7478 * e * numpy.exp(-b * e - c * e * e)

    emin, emax = 0.0, 0.6
    pmin, pmax = 0.0, 2.64456

    es = numpy.linspace(emin, emax, int((emax - emin) / 0.000001) + 1)
    pe_base = pe(es)
    # normalize
    pe_base = pe_base / numpy.sum(pe_base)
    rbe = numpy.random.choice(es, num, p=pe_base)
    return rbe

def check_in(interval):
    if interval[0] <= 0 <= interval[1]:
        return True
    else:
        return False

def f(x, a, b, c, d, e, f):
    return a * x[0] ** 2 + b * x[0] * x[1] + c * x[1] ** 2 + d * x[0] + e * x[1] + f

def smooth(image,size):
    my, mx = numpy.mgrid[-2:3, -2:3]
    x, y = mx.reshape((1, 25)), my.reshape((1, 25))
    cen = int((size * size + size) / 2)
    fit_img = numpy.zeros_like(image)
    for i in range(size):
        for j in range(size):
            arr = numpy.zeros((size, size))
            arr[int(size / 2), int(size / 2)] = 0.5
            pos = []
            tag = 0
            pk = 0
            z = []
            x_d = []
            for m in range(-2, 3):
                p = (i + m + size) % size
                for n in range(-2, 3):
                    q = (j + n + size) % size

                    if tag not in [0, 4, 20, 24]:  # abs(m) != 2 or abs(n) != 2:
                        if p * size + q != cen:
                            pos.append((p, q))
                            z.append(image[p, q])
                            x_d.append((x[0, tag], y[0, tag]))
                        else:
                            pk = tag
                    tag += 1

            x_d = numpy.array(x_d).T
            a1, a2 = curve_fit(f, x_d, z)
            fit_img[i, j] = a1[5]
    return fit_img

def set_bin(data, bin_num):
    temp_data = numpy.sort(data[data>0])[:int(len(data[data>0])*0.99)]
    bin_size = len(temp_data)/bin_num*2
    bins = numpy.array([temp_data[int(i*bin_size)] for i in range(1, int(bin_num / 2))])
    bins = numpy.sort(numpy.append(numpy.append(-bins, [0.]), bins))
    bound = numpy.max(numpy.abs(data)) * 10000.
    bins = numpy.append(-bound, numpy.append(bins, bound))
    return bins


def field_dict(expo_file):
    # to build a dictionary that contains the exposures as {"field": {"exposure":[expo1, expo2..]}....}
    with open(expo_file, "r") as f:
        contents = f.readlines()
    file_dict = {}
    fields = [] # to avoid the disorder of the field sequence
    for c in contents:
        c = c.split("\n")[0]
        if "w" in c:
            field = c
            file_dict.setdefault(field, {})
            fields.append(field)
        else:
            chip = c.split(".fits")[0]
            expo = c.split("_")[0]
            if expo in file_dict[field].keys():
                file_dict[field][expo].append(chip)
            else:
                file_dict[field].setdefault(expo, [chip])
    return file_dict, fields


def allot(allot_list, fractions):
    num = len(allot_list)
    m, n = divmod(num, fractions)
    pool = []
    for i in range(fractions):
        temp = []
        if i == fractions-1:
            more = n
        else:
            more = 0
        for j in range(m*i, m*(i+1)+more):
            temp.append(allot_list[j])
        pool.append(temp)
    return pool

def file_name(path):
    ex = os.path.exists(path)
    ori_name = os.path.basename(path)
    ab_path = path.split(ori_name)[0]
    i = 1
    name = ori_name.split(".")
    n = len(name)
    if n == 1:
        a = name[0]
        b = ""
    else:
        a = ".".join(name[:n-1])
        b = ".".join(("",name[-1]))
    while ex:
        path = "".join((ab_path, "".join((a, "_%s"%str(i), b))))
        ex = os.path.exists(path)
        i += 1
    return path