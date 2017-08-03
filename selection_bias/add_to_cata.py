import pandas as pd
import os
import time
from multiprocessing import Pool
import numpy

def cat_add(path,chip_num,stampsize,id_tag,type):
    for k in range(chip_num):
        kk = str(k).zfill(2)
        print('Process %d: chip_%s start>>>>'%(id_tag,kk))
        gal_info_path = path + 'gal_info_' + kk + '.xlsx'
        cata_path  =  path + 'gal_chip_' + kk + '.fits.cat'
        gal_info = pd.read_excel(gal_info_path)
        cata_data = numpy.loadtxt(cata_path)
        x = cata_data[:,2]
        y = cata_data[:,3]
        snr = cata_data[:,0]
        sex_data = numpy.zeros((10000,3))
        for i in range(len(snr)):
            xx = x[i]
            yy = y[i]
            mx,modx = divmod(xx,stampsize)
            my,mody = divmod(yy,stampsize)
            tag = int(100*my+mx)
            sex_data[tag] = snr[i],xx,yy
        col = [type] + ['sex_x','sex_y']
        df = pd.DataFrame(data=sex_data,columns=col)
        df_merge = pd.concat([gal_info,df],axis=1)
        df_merge.to_excel(gal_info_path)
        print('Process %d: chip_%s complete<<<<'%(id_tag,kk))


if __name__=='__main__':
    type = 'tophat'
    chip_num = 100
    size = 80
    head = '/lmc/selection_bias/'
    paths = [ head+str(i)+'/' for i in range(10)]
    p = Pool()
    t1 =time.time()
    for i in range(len(paths)):
        p.apply_async(cat_add,args=(paths[i],chip_num,size,i,type))
    p.close()
    p.join()
    #cat_add(paths[1],40,50)
    t2=time.time()
    print('Complete in %.2f'%(t2-t1))

