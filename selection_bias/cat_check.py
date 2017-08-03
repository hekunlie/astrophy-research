from subprocess import Popen
import os
import time
from multiprocessing import Pool


def check_cat(path):
    files = os.listdir(path)
    gal_chips = []
    n = 0
    exist_n = 0
    for file in files:
        if 'gal_chip' in file and 'cat' not in file:
            gal_chips.append(path+file)
    for exit in gal_chips:
        cat = exit+'.cat'
        if not os.path.exists(cat):
            n+=1
            print(cat)
        else:
            exist_n+=1
    if n==0:
        print('ALL exist,%d'%exist_n)

if __name__=='__main__':
    head = '/lmc/selection_bias/'
    paths = [ head+str(i)+'/' for i in range(10)]
    p = Pool()
    t1 =time.time()
    for i in range(len(paths)):
        p.apply_async(check_cat, args=(paths[i],))
    p.close()
    p.join()
    t2=time.time()
    print('Complete in %.2f'%(t2-t1))

