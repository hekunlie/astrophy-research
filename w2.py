import sys
import math
import numpy
import galsim
import os
from numpy import fft
import astropy.io
from astropy.io import fits
import time
import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.optimize import least_squares
from scipy.optimize import fmin_cg
import gc
###############################################
class Fourier_Quad:                

    def pow_spec(self, image):
        image_ps = fft.fftshift((numpy.abs(fft.fft2(image)))**2)        
        return image_ps
    
    def shear_est(self, gal, wbeta, psf, x, mx, my):
        gal_ps = self.pow_spec(gal)         
        mn1    = (-0.5)*((mx-0.5*x)**2 - (my-0.5*x)**2)
        mn2    = (-mx+0.5*x)*(my-0.5*x)
        mn3    = (mx-0.5*x)**2+(my-0.5*x)**2-0.5*wbeta[1]**2*((mx-0.5*x)**2+(my-0.5*x)**2)**2
        dg1    = (mn1*wbeta[0]/psf*gal_ps).sum()
        dg2    = (mn2*wbeta[0]/psf*gal_ps).sum()
        dn     = (mn3*wbeta[0]/psf*gal_ps).sum()
        return dg1,dg2,dn
    
    def wbeta(self, beta, x, mx, my):
        w_temp  = numpy.exp(-((mx-0.5*x)**2+(my-0.5*x)**2)/2./beta**2)        
        wk_beta = self.pow_spec(w_temp)
        beta1   = beta*2.0*math.pi/x        
        return wk_beta,beta1
        
    def ran_pos(self,num,imagesize):
        position   = numpy.matrix(numpy.zeros((2,num)))
        position[0]= numpy.matrix(numpy.random.normal(loc=0., scale=imagesize/12., size=num))
        position[1]= numpy.matrix(numpy.random.normal(loc=0., scale=imagesize/12., size=num))
        return position
        
    def rotate(self, pos, theta):
        rot_matrix = numpy.matrix([[math.cos(theta),math.sin(theta)],[-math.sin(theta),math.cos(theta)]])
        rot_pos    = rot_matrix*pos
        return rot_pos
    
    def shear(self, pos, g1, g2):
        shear_matrix = numpy.matrix(([(1+g1)/(1-g1**2-g2**2),g2/(1-g1**2-g2**2)],[g2/(1-g1**2-g2**2),(1-g1)/(1-g1**2-g2**2)]))
        shear_pos    = shear_matrix*pos
        return shear_pos
        
    def image_arr(self, pos, psf_scale,imagesize, mx, my,psf="GAUSS"):
        x   = pos.shape[1]
        pos = numpy.array(pos+imagesize/2)
        arr = numpy.zeros((imagesize, imagesize))
        
        if psf is 'GAUSS':
            for i in range(x):
                arr += numpy.exp(-((mx-pos[0,i])**2+(my-pos[1,i])**2)/2./psf_scale**2)
            return arr
        elif psf is "Moffat":
            for l in range(x):
                hstep = 3*psf_scale-numpy.sqrt((mx-pos[0,l])**2+(my-pos[1,l])**2)
                idx   = hstep < 0.
                hstep[idx] = 0.
                idx   = hstep!= 0.
                hstep[idx] = 1.            
                arr += (1+((mx-pos[0,l])**2+(my-pos[1,l])**2)/psf_scale**2)**(-3.5)*hstep
            return arr
            
    def cre_psf(self, psf_scale, imagesize, psf="GAUSS",x=0,y=0):
        xx = numpy.linspace(0,imagesize-1,imagesize)
        mx, my = numpy.meshgrid(xx,xx)
        if psf is 'GAUSS':
            arr = numpy.exp(-((mx-imagesize/2.+x)**2+(my-imagesize/2.+y)**2)/2./psf_scale**2)
            return arr
        
        if psf is 'Moffat':
            hstep = 3*psf_scale-numpy.sqrt((mx-imagesize/2.+x)**2+(my-imagesize/2.+y)**2)
            idx   = hstep < 0.
            hstep[idx] = 0.
            idx   = hstep!= 0.
            hstep[idx] = 1.
            arr   = (1+((mx-imagesize/2.+x)**2+(my-imagesize/2.+y)**2)/psf_scale**2)**(-3.5)*hstep
            return arr

    def gethlr(self, image,disp=1):# the sum of the image array must be positive!!! or it will fail
        size = image.shape[0]
        ra   = size/2*10
        half = 0.
        maxi = numpy.max(image)
        flux = numpy.sum(image)        
        y,x  = numpy.where(image==maxi)
        y    = int(y)
        x    = int(x)
        yy,xx= numpy.mgrid[0:size,0:size]        
        xx  = numpy.abs(xx-x)
        yy  = numpy.abs(yy-y)
        if flux >0:
            for r in range(10,ra,2):
                if half<flux/2.:                    
                    cir = xx**2+yy**2-(r/10.)**2
                    idx=cir<=0.
                    cir[idx]=1.
                    idx=cir>1.
                    cir[idx]=0.
                    half = numpy.sum(cir*image)
                    hlr  = r/10.
                else:
                    break
        else:
            if disp==1:                
                print "Failed! The sum of image array is non-positive!"
            hlr = 2.        
        return hlr

    def get_hlr(arr):
        maxi = numpy.max(arr)
        y, x = numpy.where(arr == maxi)
        idx = arr < maxi / 2
        arr[idx] = 0.
        idx = arr > 0.
        arr[idx] = 1.
        pool = []
        pool.append((int(y), int(x)))

        def check(x):  # x is a two components  tuple -like input
            if arr[x[0] - 1, x[1]] == 1 and (x[0] - 1, x[1]) not in pool:
                pool.append((x[0] - 1, x[1]))
            if arr[x[0] + 1, x[1]] == 1 and (x[0] + 1, x[1]) not in pool:
                pool.append((x[0] + 1, x[1]))
            if arr[x[0], x[1] + 1] == 1 and (x[0], x[1] + 1) not in pool:
                pool.append((x[0], x[1] + 1))
            if arr[x[0], x[1] - 1] == 1 and (x[0], x[1] - 1) not in pool:
                pool.append((x[0], x[1] - 1))
            return len(pool)

        while True:
            for cor in pool:
                num0 = len(pool)
                num1 = check(cor)
            if num0 == num1:
                break
        return numpy.sqrt(len(pool) / numpy.pi)

    def move(self, image,x,y):
        imagesize = image.shape[0]
        cent      = numpy.where(image==numpy.max(image))   
        dx        = int(numpy.max(cent[1])-x)
        dy        = int(numpy.max(cent[0])-y)
        if dy > 0:
            if dx > 0:
                arr_y = image[0:dy,0:imagesize]
                arr   = image[dy:imagesize,0:imagesize]
                arr   = numpy.row_stack((arr,arr_y))
                arr_x = arr[0:imagesize,0:dx]
                arr   = arr[0:imagesize,dx:imagesize]
                arr   = numpy.column_stack((arr,arr_x))
                return arr
            elif dx == 0:
                arr_y = image[0:dy,0:imagesize]
                arr   = image[dy:imagesize,0:imagesize]
                arr   = numpy.row_stack((arr,arr_y))
                return arr
            else:
                arr_y = image[0:dy,0:imagesize]
                arr   = image[dy:imagesize,0:imagesize]
                arr   = numpy.row_stack((arr,arr_y))
                arr_x = arr[0:imagesize,0:imagesize+dx]
                arr   = arr[0:imagesize,imagesize+dx:imagesize]
                arr   = numpy.column_stack((arr,arr_x))
                return arr
        elif dy==0:
            if dx >0:
                arr_x = image[0:imagesize,0:dx]
                arr   = image[0:imagesize,dx:imagesize]
                arr   = numpy.column_stack((arr,arr_x))
                return arr
            elif dx==0:
                return image
            else:
                arr    = image[0:imagesize,0:imagesize+dx]
                arr_x  = image[0:imagesize,imagesize+dx:imagesize]
                arr    = numpy.column_stack((arr_x,arr))
                return arr
        elif dy <0:
            if dx >0:
                arr_y   = image[imagesize+dy:imagesize,0:imagesize]
                arr     = image[0:imagesize+dy,0:imagesize]
                arr     = numpy.row_stack((arr_y,arr))
                arr_x   = arr[0:imagesize,0:dx]
                arr     = arr[0:imagesize,dx:imagesize]
                arr     = numpy.column_stack((arr,arr_x))
                return arr
            elif dx==0:
                arr_y   = image[imagesize+dy:imagesize,0:imagesize]
                arr     = image[0:imagesize+dy,0:imagesize]
                arr     = numpy.row_stack((arr_y,arr))
                return arr
            else:
                arr_y   = image[imagesize+dy:imagesize,0:imagesize]
                arr     = image[0:imagesize+dy,0:imagesize]
                arr     = numpy.row_stack((arr_y,arr))
                arr_x   = arr[0:imagesize,0:imagesize+dx]
                arr     = arr[0:imagesize,imagesize+dx:imagesize]
                arr     = numpy.column_stack((arr,arr_x))
                return arr
            
    def get_centroid(self,image):
        m0  = numpy.sum(image)
        y,x = numpy.mgrid[0:image.shape[0],0:image.shape[1]]
        x0  = numpy.sum(x*image)/m0
        y0  = numpy.sum(y*image)/m0
        return x0,y0
    
    def psf_align(self,image):
        imagesize= image.shape[0]
        arr      = self.move(image,0,0)
        image_f  = fft.fft2(arr)
        yy,xx    = numpy.mgrid[0:48,0:48]
        xx       = numpy.mod(xx+24,48)-24
        yy       = numpy.mod(yy+24,48)-24
        fk       = numpy.abs(image_f)**2
        line     = numpy.sort(numpy.array(fk.flat))
        idx = fk < 4*line[imagesize**2/2]
        fk[idx]  = 0
        weight   = fk/line[imagesize**2/2]
        kx   = xx*2*numpy.pi/imagesize
        ky   = yy*2*numpy.pi/imagesize        
        def pha(p):
            x, y = p
            return numpy.sum((numpy.angle(image_f*numpy.exp(-1.0j*(kx*x+ky*y))))**2*weight)
    
        res = fmin_cg(pha,[0,0],disp=False)        
        inve = fft.fftshift(numpy.real(fft.ifft2(image_f*numpy.exp(-1.0j*(kx*res[0]+ky*res[1])))))
        return inve
        
    def divide_stamps(self,image,stampsize):
        shape = image.shape
        y = shape[0]/stampsize
        x = shape[1]/stampsize    
        star=[image[iy*stampsize:(iy+1)*stampsize,ix*stampsize:(ix+1)*stampsize] for iy in range(y) for ix in range(x)]
        for i in range(x):
            if numpy.sum(star[-1])==0:
                star.pop()
        return star  # a list of psfs      


    
    def fitting(self,star_stamps,star_data,stampsize):
        psf_fit = self.divide_stamps(star_stamps,stampsize=stampsize)
        y,x = numpy.mgrid[0:48,0:48]
        psfn= len(psf_fit)
        data_list = [(self.psf_align(psf_fit[i]/numpy.sum(psf_fit[i])), 
                      star_data[i,0]-int(numpy.where(psf_fit[i]==numpy.max(psf_fit[i]))[1])+x, 
                      star_data[i,1]-int(numpy.where(psf_fit[i]==numpy.max(psf_fit[i]))[0])+y) for i in range(psfn)]   #[..[],..,([psf],[x],[y]),..]

        arr = numpy.array([[(data_list[i][0][m,n],
                             data_list[i][1][m,n],
                             data_list[i][2][m,n]) for i in range(psfn)] for m in range(stampsize) for n in range(stampsize)])
        #[[],...,[the array(3 columns) of the (m,n) element of all psf and the corespongding coordinates x and y]...]
        def func(alpha,x,y,z):
            return alpha[0]*x+alpha[1]*y+alpha[2]-z
        x0=numpy.array([0,0,0])       
        para_list=[(least_squares(func,x0,args=(arr[i][:,1],arr[i][:,2],arr[i][:,0])).x[0],
                    least_squares(func,x0,args=(arr[i][:,1],arr[i][:,2],arr[i][:,0])).x[1],
                    least_squares(func,x0,args=(arr[i][:,1],arr[i][:,2],arr[i][:,0])).x[2]) for i in range(stampsize**2)]
        para_arr = numpy.array(para_list)
        para_x=para_arr[:,0].reshape((stampsize,stampsize))
        para_y=para_arr[:,1].reshape((stampsize,stampsize))
        para_c=para_arr[:,2].reshape((stampsize,stampsize))
        
        return para_x, para_y, para_c

    def fit(self,star_stamps,star_data,stampsize):
        psf_fit_pool = self.divide_stamps(star_stamps,stampsize)
        x = star_data[:,0]
        y = star_data[:,1]
        sxx = numpy.sum(x*x)
        syy = numpy.sum(y*y)
        sxy = numpy.sum(x*y)
        sx  = numpy.sum(x)
        sy  = numpy.sum(y)
        sz  = 0
        szx = 0
        szy = 0
        for i in range(len(psf_fit_pool)):
            psf = self.psf_align(psf_fit_pool[i]/numpy.sum(psf_fit_pool[i]))
            sz += psf
            szx+= psf*x[i]
            szy+= psf*y[i]
        a = numpy.zeros((stampsize,stampsize))
        b = numpy.zeros((stampsize,stampsize))
        c = numpy.zeros((stampsize,stampsize))
        co_matr = numpy.array([[sxx,sxy,sx],[sxy,syy,sy],[sx,sy,1]])
        for m in range(stampsize):
            for n in range(stampsize):
               re = numpy.linalg.solve(co_matr,numpy.array([szx[m,n],szy[m,n],sz[m,n]]))
               a[m,n]=re[0]
               b[m,n]=re[1]
               c[m,n]=re[2]
        return a,b,c        

###########################################################
   
    
def measure(path_list):

    for list_num in range(len(path_list)):
        t1=time.time()
        stampsize = 48
        path   = path_list[list_num]
        location = path[-19:-13]
        number = path[-6::]
        print "Shear measurement of exposure %s in area %s is starting..."%(number,location)
        
        shear_path = '/run/media/lihekun/My Passport/w2/'+location +'/step2/'
        res_path = '/run/media/lihekun/My Passport/result/'+location+ '_exposure_%s.txt'%number
        res_data = open(res_path,'w+')
        res_data.writelines("KSB_e1"+"\t"+"BJ_e1"+"\t"+"RG_e1"+"\t"+"FQ_G1"+"\t"+"FG_N"+"\t"+"fg1"+"\t"
                            +"KSB_e2"+"\t"+"BJ_e2"+"\t"+"RG_e2"+"\t"+"FQ_G2"+"\t"+"FG_N"+"\t"+"fg2"+"\n")
        
        for k in range(1,chipsnum+1):
            try:
                kk = str(k).zfill(2)
                gal_img_path   = path.strip(number)+'gal_%s_%s.fits'%(number,kk)
                gal_data_path  = path.strip(number)+'gal_info%s_%s.dat'%(number,kk)
                star_img_path  = path.strip(number)+'star_%s_%s.fits'%(number,kk)
                star_data_path = path.strip(number)+'star_info%s_%s.dat'%(number,kk)
                shear_data_path= shear_path+"shear_info%s_%s.dat"%(number,kk)
            
                gal_stamps = fits.open(gal_img_path)[0].data
                gal_pool   = Fourier_Quad().divide_stamps(gal_stamps,stampsize)
                gal_data   = numpy.loadtxt(gal_data_path,skiprows=1)[:,17:19]
                star_stamps= fits.open(star_img_path)[0].data
                star_data  = numpy.loadtxt(star_data_path,skiprows=1)[:,1:3]
                shear_data = numpy.loadtxt(shear_data_path,skiprows=1)[:,31:33]
                ax,by,c    = Fourier_Quad().fit(star_stamps,star_data,stampsize)
                galnum = len(gal_pool)
                if len(shear_data)==0:
                    shear_data=numpy.zeros((galnum,2))
            except:
                print "Something went wrong in the exposure %s of location %"%(number,location)
            else:
                for i in range(galnum):

                    galo = gal_pool[i]
                    gal_x= gal_data[i,0]
                    gal_y= gal_data[i,1]
                    psfo  = (gal_x-stampsize/2)*ax+(gal_y-stampsize/2)*by+c

                    if numpy.sum(galo[46:48])==0:
                        galo = galo[0:32,0:32]
                        psfo = psfo[8:40,8:40]
                
                    gal_f = galo
                    psf_f = psfo
                    gal = galsim.Image(galo)
                    psf = galsim.Image(psfo)
                
                    res_k = galsim.hsm.EstimateShear(gal,psf,shear_est='KSB',strict=False)

                    res_b = galsim.hsm.EstimateShear(gal,psf,shear_est='BJ',strict=False)
                
                    res_r = galsim.hsm.EstimateShear(gal,psf,shear_est='REGAUSS',strict=False)
                
                    image_size = gal_f.shape[0]
                    beta   = 1.4*Fourier_Quad().gethlr(psf_f,disp=False)
                    xx     = numpy.linspace(0,image_size-1,image_size)
                    mx,my  = numpy.meshgrid(xx,xx)
                    w_beta = Fourier_Quad().wbeta(beta, image_size, mx, my)
                    psf_ps = Fourier_Quad().pow_spec(psf_f)
                    G1,G2,N= Fourier_Quad().shear_est(gal_f, w_beta, psf_ps, image_size, mx, my)

                    res_data.writelines(str(res_k.corrected_g1)+"\t"+str(res_b.corrected_e1)+"\t"+str(res_r.corrected_e1)+"\t"+str(G1)+"\t"+str(N)+"\t"
                                    +str(shear_data[i,0])+"\t"+str(res_k.corrected_g2)+"\t"+str(res_b.corrected_e2)+"\t"+str(res_b.corrected_e2)
                                    +"\t"+str(G2)+"\t"+str(N)+"\t"+str(shear_data[i,1])+'\n')        
        res_data.close()
        t2=time.time()
        print "Progress %d/%d:exposure %s in %s area has been done within %f sec."%(list_num+1,len(path_list),number,location,t2-t1)
        #gc.collect()
    
if __name__=="__main__":
    corenum =3  
    chipsnum = 36
    paths   = []
    paths_pool = {}
    data    = open('nname.dat')
    print "open nname.data"
    datalen = len(data.readlines())
    data.seek(0)

    for i in range(datalen):
        content = data.readline().split()[0][0:6]
        if 'w' in content:
            head = content
        else:
            path = '/run/media/lihekun/My Passport/w2/' + head + '/step1/' + content  
            paths.append(path)
    data.seek(0)
    data.close()
    print "get all paths"
    
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
        print "Caution! Something goes wrong!!!"
    print "all paths have been distributed"
    
    print "Progress starts..."
    p = Pool()
    ts=time.time()
    for i in range(corenum):
        p.apply_async(measure,args=(paths_pool[i],))
    p.close()
    p.join()
    te=time.time()
    print "Progress completes consuming %.3f hours."%((te-ts)/3600.)
    
    
    
    
    
    
    
