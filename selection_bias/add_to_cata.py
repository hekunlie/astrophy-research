import pandas as pd
import os
import time
from multiprocessing import Pool
import numpy

def cat_add(path, chip_num, stampsize, id_tag, fil_type):
    for k in range(chip_num):
        kk = str(k).zfill(2)
        print('Process %d: chip_%s start>>>>'%(id_tag,kk))

        # the catalog in which the SNR and some other parameters are stored
        gal_info_path = path + str(id_tag) + '/gal_info_' + kk + '.xlsx'
        gal_info = pd.read_excel(gal_info_path)

        # the catalog which produced by sextractor
        cata_path  =  path + str(id_tag) + '/gal_chip_' + kk + '.fits.cat'
        cata_data = numpy.loadtxt(cata_path)

        # the shear data catalog
        data_path = path + '/result/data/' + str(id_tag) + '_chip_' + kk + '.xlsx'

        x = cata_data[:,2]
        y = cata_data[:,3]
        snr = cata_data[:,0]
        sex_data = numpy.zeros((10000,3))
        for i in range(len(snr)):
            xx = x[i]
            yy = y[i]
            mx, modx = divmod(xx, stampsize)
            my, mody = divmod(yy, stampsize)
            tag = int(100*my+mx)
            sex_data[tag] = snr[i], xx, yy

        # add the SNR data to gal_info catalog
        gal_info[fil_type] = sex_data[:, 0]
        gal_info[fil_type + '_x'] = sex_data[:, 1]
        gal_info[fil_type + '_y'] = sex_data[:, 2]
        gal_info.to_excel(gal_info_path)

        # add the SNR data to shear data catalog
        if os.path.exists(data_path):
            data = pd.read_excel(data_path)
            data[fil_type] = sex_data[:, 0]
            data.to_excel(data_path)

        print('Process %d: chip_%s complete<<<<'%(id_tag, kk))

if __name__=='__main__':
    filter_type = 'gauss_5_9'
    chip_num = 100
    size = 80
    head = '/lmc/selection_bias/'
    p = Pool()
    t1 =time.time()
    for i in range(10):
        p.apply_async(cat_add,args=(head, chip_num, size, i, filter_type))
    p.close()
    p.join()
    #cat_add(paths[1],40,50)
    t2=time.time()
    print('Complete in %.2f'%(t2-t1))


