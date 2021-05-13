import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
import numpy
import ctypes
import numpy.ctypeslib as ctl
libc4py = ctypes.cdll.LoadLibrary("%s/work/mylib/libc4py.so"%my_home)
import tool_box



deblend_c = libc4py.deblend
deblend_c.restype = None
deblend_c.argtypes = [ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                     ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                     ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                     ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                     ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous'),
                     ctypes.c_int,
                     ctypes.c_float,
                     ctypes.c_float,
                     ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous')]


def deblend(xc, yc, z, radius, chip_labels, sep_radius_scale, sep_z):
    xc = numpy.ascontiguousarray(xc, dtype=numpy.float32)
    yc = numpy.ascontiguousarray(yc, dtype=numpy.float32)
    z = numpy.ascontiguousarray(z, dtype=numpy.float32)
    radius = numpy.ascontiguousarray(radius, dtype=numpy.float32)
    chip_labels = numpy.ascontiguousarray(chip_labels, dtype=numpy.intc)

    data_len = xc.shape[0]
    labels = numpy.ones((data_len,), dtype=numpy.intc)

    deblend_c(xc, yc, z, radius, chip_labels, data_len, sep_radius_scale, sep_z, labels)

    return labels


deblend_i_c = libc4py.deblend_i
deblend_i_c.restype = None
deblend_i_c.argtypes = [ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                         ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                         ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                         ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous'),
                         ctypes.c_int,
                         ctypes.c_float,
                         ctypes.c_float,
                         ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous')]


def deblend_i(xc, yc, z, chip_labels, sep_pix, sep_z):
    xc = numpy.ascontiguousarray(xc, dtype=numpy.float32)
    yc = numpy.ascontiguousarray(yc, dtype=numpy.float32)
    z = numpy.ascontiguousarray(z, dtype=numpy.float32)
    chip_labels = numpy.ascontiguousarray(chip_labels, dtype=numpy.intc)

    data_len = xc.shape[0]
    labels = numpy.ones((data_len,), dtype=numpy.intc)

    deblend_i_c(xc, yc, z, chip_labels, data_len, sep_pix, sep_z, labels)

    return labels


deblend_mutual_c = libc4py.deblend_mutual
deblend_mutual_c.restype = None
deblend_mutual_c.argtypes = [ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctypes.c_int,
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctypes.c_int,
                                ctypes.c_float,
                                ctypes.c_float,
                                ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous')]

def deblend_mutual(ra_i, dec_i, z_i, ra_j, dec_j, z_j, sep_ang, sep_z, labels_i):
    ra_i = numpy.ascontiguousarray(ra_i, dtype=numpy.float32)
    dec_i = numpy.ascontiguousarray(dec_i, dtype=numpy.float32)
    z_i = numpy.ascontiguousarray(z_i, dtype=numpy.float32)

    ra_j = numpy.ascontiguousarray(ra_j, dtype=numpy.float32)
    dec_j = numpy.ascontiguousarray(dec_j, dtype=numpy.float32)
    z_j = numpy.ascontiguousarray(z_j, dtype=numpy.float32)

    labels_i = numpy.ascontiguousarray(labels_i, dtype=numpy.intc)

    data_len_i = ra_i.shape[0]
    data_len_j = ra_j.shape[0]

    deblend_mutual_c(ra_i, dec_i, z_i, data_len_i, ra_j, dec_j, z_j, data_len_j, sep_ang, sep_z, labels_i)



cata_match_c = libc4py.cata_match
cata_match_c.restype = None
cata_match_c.argtypes = [ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctypes.c_int,
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                                ctypes.c_int,
                                ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous'),
                                ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous'),
                                ctypes.c_float]

def cata_match(cata_ra1, cata_dec1, cata_ra2, cata_dec2,  diff_ang_thresh):
    cata_len1 = cata_ra1.shape[0]
    cata_len2 = cata_ra2.shape[0]

    cata_ra1 = numpy.ascontiguousarray(cata_ra1, dtype=numpy.float32)
    cata_dec1 = numpy.ascontiguousarray(cata_dec1, dtype=numpy.float32)

    cata_ra2 = numpy.ascontiguousarray(cata_ra2, dtype=numpy.float32)
    cata_dec2 = numpy.ascontiguousarray(cata_dec2, dtype=numpy.float32)

    seq_labels = numpy.zeros((cata_len1, ), dtype=numpy.intc) - 1
    check_labels = numpy.zeros((cata_len1, ), dtype=numpy.intc)

    cata_match_c(cata_ra1, cata_dec1, cata_len1, cata_ra2, cata_dec2, cata_len2,
                 seq_labels, check_labels, diff_ang_thresh)

    idx = seq_labels < 0
    miss_num = idx.sum()
    if idx.sum() > 0:
        print("%d sources don't have pairs!"%miss_num)
    idx = check_labels > 1
    multi_pair = idx.sum()
    if multi_pair > 0:
        print("%d sources have more than one pair!" % multi_pair)

    return seq_labels, check_labels


find_overlap_c = libc4py.find_overlap
find_overlap_c.restype = None
find_overlap_c.argtypes = [ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                            ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                            ctypes.c_int,
                            ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                            ctl.ndpointer(numpy.float32, flags='aligned, c_contiguous'),
                            ctypes.c_int,
                            ctypes.c_int,
                            ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous'),
                            ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous')]


def find_overlap_mask(ra, dec, mask, ra_bin, dec_bin, extend_step=0):
    """
    :param ra: the source ra to be selected
    :param dec: the source dec to be selected
    :param mask: the mask of the target survey, pixels must have values >0 or 0
    :param ra_bin: the bins to make the mask
    :param dec_bin: the bins to make the mask
    :param extend_step: how many pixels beyond the edge of the mask to be added
    :return:
    """

    mask_c = numpy.zeros_like(mask) + mask

    data_len = ra.shape[0]
    ra_bin_num = ra_bin.shape[0] - 1
    dec_bin_num = dec_bin.shape[0] - 1

    labels = numpy.zeros_like(ra, dtype=numpy.intc)

    ra = numpy.ascontiguousarray(ra, dtype=numpy.float32)
    dec = numpy.ascontiguousarray(dec, dtype=numpy.float32)

    ra_bin = numpy.ascontiguousarray(ra_bin, dtype=numpy.float32)
    dec_bin = numpy.ascontiguousarray(dec_bin, dtype=numpy.float32)

    source_area = []
    for iy in range(dec_bin_num):
        for ix in range(ra_bin_num):
            if mask_c[iy, ix] > 0:
                source_area.append([iy, ix])

    if extend_step > 0:
        tool_box.edge_extend(mask_c, dec_bin_num, ra_bin_num, source_area, extend_step)

    mask_c = mask_c.astype(dtype=numpy.intc).flatten()
    mask_c = numpy.ascontiguousarray(mask_c, dtype=numpy.intc)

    find_overlap_c(ra, dec, data_len, ra_bin, dec_bin, ra_bin_num, dec_bin_num, mask_c, labels)
    idx_n = labels > 0

    return labels, idx_n


