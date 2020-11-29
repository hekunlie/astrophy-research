import numpy


def set_min_bin(x_min, x_max, bin_width, dx=1):
    # set up bins
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
    # return the bin labels of each data point
    bins_label = numpy.zeros_like(data, dtype=numpy.intc)
    for i in range(bins_num):
        idx1 = data >= bins[i]
        idx2 = data < bins[i + 1]
        idx = idx1 & idx2
        bins_label[idx] = i
    return bins_label

def even_area(num_in_area, area_num, total_ncent):
    total_num = num_in_area.sum()
    raw_ncent = num_in_area/total_num*total_ncent
    # print(raw_ncent)
    ncent = numpy.zeros((area_num,), dtype=numpy.intc)
    ncent_bk = numpy.zeros((area_num,), dtype=numpy.intc)
    for i in range(area_num):
        sub_num = int(raw_ncent[i])
        if sub_num == 0:
            ncent[i] = 1
            ncent_bk[i] = 1
        else:
            ncent[i] = sub_num
            ncent_bk[i] = sub_num

    if ncent.sum() > total_ncent:
        idx = ncent == ncent.max()
        ncent[idx] = ncent[idx] - 1
    else:
        num_in_area_sort = numpy.sort(num_in_area)

        diff = total_ncent - ncent.sum()
        # print(ncent, ncent.sum(), diff)

        for i in range(area_num-1,-1,-1):
            if diff == 0:
                break
            for j in range(area_num):
                if numpy.abs(num_in_area[j] - num_in_area_sort[i])<0.01:
                    ncent[j] += 1
                    diff -= 1
                    break
    if ncent.sum() != total_ncent:
        print(ncent,ncent.sum())
        print(ncent_bk,ncent_bk.sum())
        print(ncent - ncent_bk)
    return ncent,num_in_area/ncent