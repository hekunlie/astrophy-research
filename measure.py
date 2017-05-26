#import galsim
import os
from astropy.io import fits
import time
from multiprocessing import Pool
from Fourier_Quad import *

def measure(path_list,tag):
    ahead = '/lmc/w1/'
    res_ahead = '/home/hklee/result/w1/'

    for list_num in range(len(path_list)):
        t1=time.time()
        stampsize = 48
        my, mx = numpy.mgrid[0:stampsize, 0:stampsize]
        path   = path_list[list_num]
        location,number = path.split('/')
        print ("Process %d: %s_%s starts..."%(tag,location,number))

        shear_path = ahead+location +'/step2/'
        res_path = res_ahead+location+ "_exposure_%s.txt"%number
        res_data = open(res_path,"w+")

        res_data.writelines("chip"+'\t'+"KSB_e1"+"\t"+"BJ_e1"+"\t"+"RG_e1"+"\t"+"FQ_G1"+"\t"+"FG_N"+"\t"+"fg1"+"\t"
                            +"KSB_e2"+"\t"+"BJ_e2"+"\t"+"RG_e2"+"\t"+"FQ_G2"+"\t"+"FG_N"+"\t"+"fg2"+"\t"+"FQ_U"+"\t"+"FQ_V"+"\n")
        
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
                ax,by,c    = Fourier_Quad().fit(star_stamps,star_noise,star_data,stampsize)
                galnum = len(gal_pool)

                for i in range(galnum):

                    if gal_data[i,2]>=10.:
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

                        beta   = Fourier_Quad().get_radius(psf,2.)
                        w_beta = Fourier_Quad().wbeta(beta, stampsize, mx, my)
                        G1,G2,N,U,V= Fourier_Quad().shear_est(gal, w_beta, psf, stampsize, mx, my,noise,N=True)

                        res_data.writelines(kk+'\t'+str(0)+"\t"+str(0)+"\t"+str(0)+"\t"+str(G1)+"\t"+str(N)+"\t"
                                        +str(shear_data[i,0])+"\t"+str(0)+"\t"+str(0)+"\t"+str(0)
                                        +"\t"+str(G2)+"\t"+str(N)+"\t"+str(shear_data[i,1])+"\t"+str(U)+"\t"+str(V)+'\n')


        res_data.close()
        t2=time.time()

        print ("Process %d : (%d/%d) %s_%s done within %f sec."%(tag,list_num+1,len(path_list),location,number,t2-t1),)


if __name__=="__main__":
    corenum =6
    chipsnum = 36
    paths   = []
    paths_pool = {}
    data    = open('/lmc/w1/nname.dat')
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
    
    m,n = divmod(len(paths),corenum)
    if m==0 and n!=0:
        for i in range(n):
            paths_pool[i]=map(int,str(paths[i]))
    elif m!=0:
        for i in range(corenum):
            paths_pool[i]=paths[i*m:(i+1)*m]
        if n!=0:
            for i in range(n):
                paths_pool[i].append(paths[-i-1])
    else:
        print( "Caution! Something goes wrong!!!")
    print ("all paths have been distributed")
    
    print( "Progress starts...")
    p = Pool()
    ts=time.time()
    for i in range(corenum):
        p.apply_async(measure,args=(paths_pool[i],i,))
    p.close()
    p.join()
    te=time.time()
    print ("Progress completes consuming %.3f hours."%((te-ts)/3600.))


    
    
    
