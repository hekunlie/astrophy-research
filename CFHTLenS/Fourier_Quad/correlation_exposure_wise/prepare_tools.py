import numpy


def set_min_bin(x_min, x_max, bin_width, dx=1):
    bin_num = 0
    xbin = []
    while True:
        x_ = x_min - dx + bin_num * bin_width
        xbin.append(x_)
        bin_num += 1
        if x_ > x_max:
            break
    bin_num = len(xbin) - 1
    return numpy.array(xbin), bin_num


def get_bin_label(data, bins, bins_num):
    bins_label = numpy.zeros_like(data, dtype=numpy.intc)
    for i in range(bins_num):
        idx1 = data >= bins[i]
        idx2 = data < bins[i + 1]
        idx = idx1 & idx2
        bins_label[idx] = i
    return bins_label
