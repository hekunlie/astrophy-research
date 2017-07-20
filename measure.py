#import galsim
import os
from astropy.io import fits
import time
from multiprocessing import Pool
from Fourier_Quad import  Fourier_Quad
from sys import argv
import numpy
import pandas
import tool_box

def measure(path_list,tag,area):
    ahead = '/lmc/'+area+'/'
    res_ahead = '/home/hklee/result/'+area+'/'
    stampsize = 48
    for list_num in range(len(path_list)):
        t1=time.time()
        path   = path_list[list_num]
        location,number = path.split('/')
        print ("Process %d: %s_%s starts..."%(tag,location,number))

        shear_path = ahead+location +'/step2/'
        res_path = res_ahead+location+ "_exposure_%s.xlsx"%number
        data_col = ["KSB_e1","BJ_e1","RG_e1","FQ_G1","FG_N","fg1", "KSB_e2","BJ_e2","RG_e2","FQ_G2","FG_N","fg2","FQ_U","FQ_V"]
        p = 0
        for k in range(1,37):
            kk = str(k).zfill(2)
            gal_img_path   = ahead+location+'/step1/'+'gal_%s_%s.fits'%(number,kk)
            gal_data_path  = ahead+location+'/step1/'+'gal_info%s_%s.dat'%(number,kk)
            star_img_path  = ahead+location+'/step1/'+'star_%s_%s.fits'%(number,kk)
            star_data_path = ahead+location+'/step1/'+'star_info%s_%s.dat'%(number,kk)
            shear_data_path= shear_path+"shear_info%s_%s.dat"%(number,kk)
            noise_path = ahead+location+'/step1/'+'noise'+'%s_%s.fits'%(number,kk)
            star_noise_path = ahead+location+'/step1/'+'star_noise_'+'%s_%s.fits'%(number,kk)

            if os.path.getsize(gal_data_path)/1024. < 20 or os.path.getsize(shear_data_path)/1024. < 20 or os.path.getsize(star_data_path)<2100:
                print ('Process %d: skipped chip %s'%(tag,kk))
            else:
                star_data = numpy.loadtxt(star_data_path, skiprows=1)[:, 1:3]
                gal_stamps = fits.open(gal_img_path)[0].data
                gal_pool   = Fourier_Quad().divide_stamps(gal_stamps,stampsize)
                gal_data   = numpy.loadtxt(gal_data_path,skiprows=1)[:,17:20]

                star_stamps= fits.open(star_img_path)[0].data
                star_noise = fits.open(star_noise_path)[0].data

                noise_stamps = fits.open(noise_path)[0].data
                noise_pool   = Fourier_Quad().divide_stamps(noise_stamps,stampsize)

                shear_data = numpy.loadtxt(shear_data_path,skiprows=1)[:,31:33]
                ax,by,c    = Fourier_Quad().fit(star_stamps,star_noise,star_data,stampsize,mode=2)
                galnum = len(gal_pool)

                for i in range(galnum):
                    if gal_data[i,2]>=10.:
                        # index = kk+"_"+str(i).zfill(galnum_len)
                        # gal_index.append(index)
                        gal = gal_pool[i]
                        noise = noise_pool[i]
                        gal_x= gal_data[i,0]
                        gal_y= gal_data[i,1]
                        psf = gal_x*ax+gal_y*by+c

                        # gal = galsim.Image(galo)
                        # psf = galsim.Image(psfo)

                        # res_k = galsim.hsm.EstimateShear(gal,psf,shear_est='KSB',strict=False)
                        # res_k.corrected_g1
                        # res_b = galsim.hsm.EstimateShear(gal,psf,shear_est='BJ',strict=False)
                        # res_b.corrected_e1
                        # res_r = galsim.hsm.EstimateShear(gal,psf,shear_est='REGAUSS',strict=False)
                        # res_r.corrected_e1

                        G1,G2,N,U,V= Fourier_Quad().shear_est(gal, psf, stampsize, noise)
                        ith_row = numpy.array([0, 0, 0, G1, N, shear_data[i,0], 0, 0, 0, G2, N, shear_data[i,1], U, V ])
                        if p==0 and k==1:
                            data_matrix =ith_row
                        else:
                            data_matrix = numpy.row_stack((data_matrix,ith_row))
                        p=1
        df = pandas.DataFrame(data_matrix,columns=data_col )
        df.columns.name='Chip&NO'
        df.to_excel(res_path)
        t2=time.time()

        print ("Process %d: %s_%s done within %.2f sec (%d/%d)"%(tag,location,number,t2-t1,list_num+1,len(path_list)))


if __name__=="__main__":
    area ,ncpu= argv[1:3]
    corenum =int(ncpu)
    chipsnum = 36
    paths   = []
    paths_pool = {}
    data    = open('/lmc/'+area+'/nname.dat')
    print( "open nname.data")
    datalen = len(data.readlines())
    data.seek(0)

    for i in range(datalen):
        content = data.readline().split()[0]
        if 'w' in content:
            head = content
        else:
            path =  head + '/'+content
            paths.append(path)
    data.seek(0)
    data.close()
    print( "get all paths")
    paths_pool = tool_box.task_distri(paths,corenum)
    # m,n = divmod(len(paths),corenum)
    # if m==0 and n!=0:
    #     for i in range(n):
    #         paths_pool[i]=map(int,str(paths[i]))
    # elif m!=0:
    #     for i in range(corenum):
    #         paths_pool[i]=paths[i*m:(i+1)*m]
    #     if n!=0:
    #         for i in range(n):
    #             paths_pool[i].append(paths[-i-1])
    # else:
    #     print( "Caution! Something goes wrong!!!")
    # print ("all paths have been distributed")
    
    print( "Progress starts...")
    p = Pool()
    ts=time.time()
    for i in range(corenum):
        p.apply_async(measure,args=(paths_pool[i],i,area))
    p.close()
    p.join()
   #measure(paths_pool[5],5,area)
    te=time.time()
    print("Progress completes consuming %.3f hours." % ((te - ts) / 3600.))
    os.system('python process_data.py')



    
    
    
