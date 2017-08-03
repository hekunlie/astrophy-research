import pandas as pd
import os
import time
from multiprocessing import Pool
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def check(path,chip_num,stampsize,tag):
    for k in range(chip_num):
        kk = str(k).zfill(2)
        cata_path  =  path + 'gal_chip_' + kk + '.fits.cat'
        fig = '/home/hklee/test/'+str(tag)+'_'+kk+'.png'
        cata_data = numpy.loadtxt(cata_path)
        x = cata_data[:,2]
        y = cata_data[:,3]
        arr = numpy.zeros((100,100))
        for i in range(len(x)):
            xx = x[i]
            yy = y[i]
            mx,modx = divmod(xx,stampsize)
            my,mody = divmod(yy,stampsize)
            arr[int(my),int(mx)]+=10
        plt.figure()
        plt.imshow(arr)
        plt.colorbar()
        plt.savefig(fig)
        plt.close()


if __name__=='__main__':
    chip_num = 40
    size = 50
    head = '/lmc/selection_bias/'
    paths = [ head+str(i)+'/' for i in range(10)]
    p = Pool()
    t1 =time.time()
    for i in range(len(paths)):
        p.apply_async(check,args=(paths[i],chip_num,size,i))
    p.close()
    p.join()
    #check(paths[0],chip_num,size,0)
    t2=time.time()
    print('Complete in %.2f'%(t2-t1))

