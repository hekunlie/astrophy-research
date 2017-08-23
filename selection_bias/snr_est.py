from subprocess import Popen
import os
import time
from multiprocessing import Pool


def sextract(path):
    files = os.listdir(path)
    gal_chips = []
    for file in files:
        if 'gal_chip' in file:
            gal_chips.append(path+file)
    for i in range(len(gal_chips)):
        comm = 'sex %s -CATALOG_NAME %s.cat '%(gal_chips[i], gal_chips[i])
        a = Popen(comm, shell=True)
        a.wait()

if __name__=='__main__':
    head = '/lmc/selection_bias/'
    paths = [ head+str(i)+'/' for i in range(10)]
    p = Pool()
    t1 = time.time()
    for i in range(len(paths)):
        p.apply_async(sextract,args=(paths[i],))
    p.close()
    p.join()
    t2=time.time()
    print('Complete in %.2f'%(t2-t1))
    os.system('sudo python cat_check.py')

