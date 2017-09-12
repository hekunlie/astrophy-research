from sys import path
path.append('/home/hklee/work/fourier_quad/')
import numpy
#import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
#import lsstetc
import time
import os
from multiprocessing import Pool
import pandas
import tool_box

def simu(gal_paths_list, info_paths_list, process_id):
    pass


if __name__=='__main__':
    CPU_num = 20
    chip_num = 250

    for i in range(10):
        files_path = '/lmc/selection_bias/%d/'%i
        if not os.path.isdir(files_path):
            os.mkdir(files_path)

    chip_paths_pool = ['/lmc/selection_bias/%d/gal_chip_%d.fits'%(i, j) for i in range(10) for j in range(chip_num)]
    info_path_pool = ['/lmc/selection_bias/%d/gal_info_%d.xlsx'%(i, j) for i in range(10) for j in range(chip_num)]
    chip_paths_list = tool_box.task_distri(chip_paths_pool, 20)
    info_path_list = tool_box.task_distri(info_path_pool, 20)

    arr = numpy.load('/lmc/selection_bias/shear.npz')
    shear1 = arr['arr_0']
    shear2 = arr['arr_1']
    mags = numpy.load('/home/hklee/work/selection_bias/parameters/lsstmagsims.npz')['arr_0']
    gal_rad = numpy.load('/home/hklee/work/selection_bias/parameters/gal_radius.npz')['arr_0']
    sersic_rad = numpy.load('/home/hklee/work/selection_bias/parameters/sersic_rd.npz')
    sersic_bulge = sersic_rad['arr_0']
    sersic_disk = sersic_rad['arr_1']
    e1e2 = numpy.load('/home/hklee/work/selection_bias/parameters/e1e2.npz')
    ie1 = e1e2['arr_0']
    ie2 = e1e2['arr_1']

    p = Pool()
    t1 = time.time()
    for m in range(CPU_num):
        p.apply_async(simu, args=(chip_paths_list, info_path_list, m,))
    p.close()
    p.join()
    # simu(shear1[0], shear2[0], ie1, ie2, gal_ra, sersic_rd, mags, 0)
    t2 = time.time()
    print('Time consuming: %.2f') % (t2 - t1)
    os.system('python selection_bias_est.py')