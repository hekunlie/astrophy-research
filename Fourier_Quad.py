import numpy
from numpy import fft
from scipy.optimize import fmin_cg
from scipy import ndimage, signal
import copy
import matplotlib.pyplot as plt
import random
from scipy import optimize
class Fourier_Quad:

    def pow_spec(self, image):
        image_ps = fft.fftshift((numpy.abs(fft.fft2(image)))**2)
        return image_ps

    def shear_est(self, gal_image, psf_image, imagesize, background_noise=None, F=False):
        x = imagesize
        my, mx = numpy.mgrid[0:x, 0:x]
        gal_ps = self.pow_spec(gal_image)

        if background_noise is not None: # to deduct the noise
            d=1
            nbg = self.pow_spec(background_noise)
            rim = self.border(d,x)
            n    = numpy.sum(rim)
            gal_pnoise = numpy.sum(gal_ps*rim)/n               #the Possion noise of galaxy image
            nbg_pnoise =numpy.sum(nbg*rim)/n                   #the  Possion noise of background noise image
            gal_ps = gal_ps  - nbg + nbg_pnoise- gal_pnoise

        if F == True:
            psf_ps = psf_image
        else:
            psf_ps = self.pow_spec(psf_image)

        hlr = self.get_radius_new(psf_ps,2.)
        wb,beta = self.wbeta(hlr,x)
        maxi = numpy.max(wb)
        idx = wb < maxi / 100000.
        wb[idx] = 0.
        maxi = numpy.max(psf_ps)
        idx = psf_ps < maxi / 100000.
        psf_ps[idx] = 1.

        tk  = wb/ psf_ps * gal_ps
        alpha =2.*numpy.pi/x
        kx = mx-0.5*x
        ky = my-0.5*x
        mn1 = (-0.5)*(kx**2 - ky**2)
        mn2 = -kx*ky
        mn3 = kx**2 + ky**2 - 0.5*beta**2*( kx**2 + ky**2 )**2
        mn4 = kx**4 - 6*kx**2*ky**2 + ky**4
        mn5 = kx**3*ky - kx*ky**3
        g1 = numpy.sum(mn1 * tk)*(alpha**4)
        g2 = numpy.sum(mn2 * tk)*(alpha**4)
        n  = numpy.sum(mn3 * tk)*(alpha**4)
        u  = numpy.sum(mn4 * tk)*(-0.5*beta**2)*(alpha**4)
        v  = numpy.sum(mn5 * tk)*(-2.*beta**2)*(alpha**4)
        return g1, g2, n, u, v

    def wbeta(self, beta, imagesize):
        my, mx = numpy.mgrid[0:imagesize, 0:imagesize]
        sigma = beta/numpy.sqrt(2)
        w_temp = numpy.exp(-((mx-0.5*imagesize)**2+(my-0.5*imagesize)**2)/2./sigma**2)
        beta = 1./beta
        return w_temp, beta

    def ran_pos(self, num, imagesize):
        position = numpy.matrix(numpy.zeros((2, num)))
        position[0] = numpy.matrix(numpy.random.normal(loc=0., scale=imagesize / 12., size=num))
        position[1] = numpy.matrix(numpy.random.normal(loc=0., scale=imagesize / 12., size=num))
        return position

    def rotate(self, pos, theta):
        rot_matrix = numpy.matrix([[numpy.cos(theta), numpy.sin(theta)], [-numpy.sin(theta), numpy.cos(theta)]])
        rot_pos = rot_matrix * pos
        return rot_pos

    def shear(self, pos, g1, g2):
        shear_matrix = numpy.matrix(([(1+g1)/(1-g1**2-g2**2), g2/(1-g1**2-g2**2)],
                                     [g2/(1-g1**2-g2**2), (1-g1)/(1-g1**2-g2**2)]))
        shear_pos = shear_matrix * pos
        return shear_pos

    def convolve_psf(self, pos, psf_scale, imagesize, psf="GAUSS"):
        x = pos.shape[1]
        my,mx = numpy.mgrid[0:imagesize,0:imagesize]
        pos = numpy.array(pos+imagesize / 2.)
        arr = numpy.zeros((imagesize, imagesize))

        if psf is 'GAUSS':
            for i in range(x):
                arr += numpy.exp(-((mx-pos[0,i])**2+(my-pos[1,i])** 2)/2./psf_scale**2)
            return arr
        elif psf is "Moffat":
            for l in range(x):
                hstep = 3 * psf_scale-numpy.sqrt((mx-pos[0,l])**2+(my-pos[1,l])**2)
                idx = hstep < 0.
                hstep[idx] = 0.
                idx = hstep != 0.
                hstep[idx] = 1.
                arr += (1+((mx-pos[0,l])**2+(my-pos[1,l])**2)/psf_scale**2)**(-3.5)*hstep
            return arr

    def cre_psf(self, psf_scale, imagesize, model="GAUSS", x=0, y=0):
        xx = numpy.linspace(0, imagesize - 1, imagesize)
        mx, my = numpy.meshgrid(xx, xx)
        if model is 'GAUSS':
            arr = numpy.exp(-((mx -imagesize/2.+x)**2+(my-imagesize/2.+y)**2)/2./psf_scale**2)
            return arr

        if model is 'Moffat':
            hstep = 3*psf_scale-numpy.sqrt((mx-imagesize/2.+x)**2+(my-imagesize/2.+y)**2)
            idx = hstep < 0.
            hstep[idx] = 0.
            idx = hstep != 0.
            hstep[idx] = 1.
            arr = (1+((mx-imagesize/2.+x)**2+(my-imagesize/2.+y)**2)/psf_scale**2)**(-3.5)*hstep
            return arr

    def get_radius(self, image, scale):
        # get the radius of the flux descends to the maximum/scale
        radi_arr = copy.copy(image)
        maxi = numpy.max(radi_arr)
        y, x = numpy.where(radi_arr == maxi)
        idx = radi_arr < maxi / scale
        radi_arr[idx] = 0.
        idx = radi_arr > 0.
        radi_arr[idx] = 1.
        half_radi_pool = []
        half_radi_pool.append((int(y[0]), int(x[0])))

        def check(x):  # x is a two components  tuple -like input
            if (x[0] - 1, x[1]) not in half_radi_pool and radi_arr[x[0] - 1, x[1]] == 1:
                half_radi_pool.append((x[0] - 1, x[1]))
            if (x[0] + 1, x[1]) not in half_radi_pool and radi_arr[x[0] + 1, x[1]] == 1:
                half_radi_pool.append((x[0] + 1, x[1]))
            if (x[0], x[1] + 1) not in half_radi_pool and radi_arr[x[0], x[1] + 1] == 1:
                half_radi_pool.append((x[0], x[1] + 1))
            if (x[0], x[1] - 1) not in half_radi_pool and radi_arr[x[0], x[1] - 1] == 1:
                half_radi_pool.append((x[0], x[1] - 1))
            return len(half_radi_pool)

        while True:
            for cor in half_radi_pool:
                num0 = len(half_radi_pool)
                num1 = check(cor)
            if num0 == num1:
                break
        return numpy.sqrt(len(half_radi_pool) / numpy.pi)

    def get_radius_new(self, image, scale):
        # get the radius of the flux descends to the maximum/scale
        radi_arr = copy.copy(image)
        maxi = numpy.max(radi_arr)
        y, x = numpy.where(radi_arr == maxi)
        idx = radi_arr < maxi / scale
        radi_arr[idx] = 0.
        half_radi_pool = []

        def detect(mask, ini_y, ini_x, signal):
            if mask[ini_y, ini_x] > 0:
                signal.append((ini_y, ini_x))
                mask[ini_y, ini_x] = 0
                for cor in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                    if mask[ini_y + cor[0], ini_x + cor[1]] > 0:
                        detect(mask, ini_y + cor[0], ini_x + cor[1], signal)
            return signal
        half_radi_pool = detect(radi_arr,y[0],x[0],half_radi_pool)

        return numpy.sqrt(len(half_radi_pool) / numpy.pi)

    def move(self, image, x, y):
        imagesize = image.shape[0]
        cent = numpy.where(image == numpy.max(image))
        dx = int(numpy.max(cent[1]) - x)
        dy = int(numpy.max(cent[0]) - y)
        if dy > 0:
            if dx > 0:
                arr_y = image[0:dy, 0:imagesize]
                arr = image[dy:imagesize, 0:imagesize]
                arr = numpy.row_stack((arr, arr_y))
                arr_x = arr[0:imagesize, 0:dx]
                arr = arr[0:imagesize, dx:imagesize]
                arr = numpy.column_stack((arr, arr_x))
                return arr
            elif dx == 0:
                arr_y = image[0:dy, 0:imagesize]
                arr = image[dy:imagesize, 0:imagesize]
                arr = numpy.row_stack((arr, arr_y))
                return arr
            else:
                arr_y = image[0:dy, 0:imagesize]
                arr = image[dy:imagesize, 0:imagesize]
                arr = numpy.row_stack((arr, arr_y))
                arr_x = arr[0:imagesize, 0:imagesize + dx]
                arr = arr[0:imagesize, imagesize + dx:imagesize]
                arr = numpy.column_stack((arr, arr_x))
                return arr
        elif dy == 0:
            if dx > 0:
                arr_x = image[0:imagesize, 0:dx]
                arr = image[0:imagesize, dx:imagesize]
                arr = numpy.column_stack((arr, arr_x))
                return arr
            elif dx == 0:
                return image
            else:
                arr = image[0:imagesize, 0:imagesize + dx]
                arr_x = image[0:imagesize, imagesize + dx:imagesize]
                arr = numpy.column_stack((arr_x, arr))
                return arr
        elif dy < 0:
            if dx > 0:
                arr_y = image[imagesize + dy:imagesize, 0:imagesize]
                arr = image[0:imagesize + dy, 0:imagesize]
                arr = numpy.row_stack((arr_y, arr))
                arr_x = arr[0:imagesize, 0:dx]
                arr = arr[0:imagesize, dx:imagesize]
                arr = numpy.column_stack((arr, arr_x))
                return arr
            elif dx == 0:
                arr_y = image[imagesize + dy:imagesize, 0:imagesize]
                arr = image[0:imagesize + dy, 0:imagesize]
                arr = numpy.row_stack((arr_y, arr))
                return arr
            else:
                arr_y = image[imagesize + dy:imagesize, 0:imagesize]
                arr = image[0:imagesize + dy, 0:imagesize]
                arr = numpy.row_stack((arr_y, arr))
                arr_x = arr[0:imagesize, 0:imagesize + dx]
                arr = arr[0:imagesize, imagesize + dx:imagesize]
                arr = numpy.column_stack((arr, arr_x))
                return arr

    def get_centroid(self, image,filt=False,radius=2):
        y0,x0 = self.mfpoly(image)
        y, x = numpy.mgrid[0:image.shape[0], 0:image.shape[1]]
        # m0 = numpy.sum(image)
        # mx  = numpy.sum(x * image)
        # my  = numpy.sum(y * image)
        # x0 = mx / m0
        # y0 = my / m0
        #yc,xc = numpy.where(image==numpy.max(image))
        y = y - y0
        x = x - x0
        if filt == True:
            gaus_filt = numpy.exp(-(y**2+x**2)/2/numpy.pi/radius**2)
            image = image *gaus_filt
        mxx = numpy.sum(x*x*image)
        mxy = numpy.sum(x * y * image)
        myy = numpy.sum(y* y * image)
        e1 = (mxx-myy)/(mxx+myy)
        e2 = 2*mxy/(mxx+myy)
        return y0, x0,e1,e2

    def psf_align(self, image):
        imagesize = image.shape[0]
        arr = self.move(image, 0, 0)
        image_f = fft.fft2(arr)
        yy, xx = numpy.mgrid[0:48, 0:48]
        xx = numpy.mod(xx + 24, 48) - 24
        yy = numpy.mod(yy + 24, 48) - 24
        fk = numpy.abs(image_f) ** 2
        line = numpy.sort(numpy.array(fk.flat))
        idx = fk < 4 * line[int(imagesize ** 2/ 2)]
        fk[idx] = 0
        weight = fk / line[int(imagesize ** 2/ 2)]
        kx = xx * 2 * numpy.pi / imagesize
        ky = yy * 2 * numpy.pi / imagesize

        def pha(p):
            x, y = p
            return numpy.sum((numpy.angle(image_f*numpy.exp(-1.0j*(kx*x+ky*y))))**2*weight)

        res = fmin_cg(pha, [0, 0], disp=False)
        inve = fft.fftshift(numpy.real(fft.ifft2(image_f*numpy.exp(-1.0j*(kx*res[0]+ky*res[1])))))
        return inve

    def gaussfilter(self, psfimage):
        x, y = numpy.mgrid[0:3, 0:3]
        xc, yc = 1.0, 1.0
        w = 1.3
        ker = (1.0/2.0/numpy.pi/w/w)*numpy.exp(-0.5*((x-xc)**2+(y-yc)**2)/2.0/w/w)
        imcov = signal.convolve(psfimage, ker, mode='same')
        return imcov

    def mfpoly(self,psf):
        p = numpy.where(psf == numpy.max(psf))
        yr, xr = p[0][0], p[1][0]
        yp, xp = numpy.mgrid[yr - 1:yr + 2, xr - 1:xr + 2]
        patch = psf[yr - 1:yr + 2, xr - 1:xr + 2]
        zz = patch.reshape(9)
        xx = xp.reshape(9)
        yy = yp.reshape(9)
        xy = xx * yy
        x2 = xx * xx
        y2 = yy * yy
        A = numpy.array([numpy.ones_like(zz), xx, yy, x2, xy, y2]).T
        cov = numpy.linalg.inv(numpy.dot(A.T, A))
        a, b, c, d, e, f = numpy.dot(cov, numpy.dot(A.T, zz))
        coeffs = numpy.array([[2.0 * d, e], [e, 2.0 * f]])
        mult = numpy.array([-b, -c])
        xc, yc = numpy.dot(numpy.linalg.inv(coeffs), mult)
        return yc, xc

    def divide_stamps(self, image, stampsize):
        shape = image.shape
        y = int(shape[0] / stampsize)
        x = int(shape[1] / stampsize)
        star = [image[iy * stampsize:(iy + 1) * stampsize, ix * stampsize:(ix + 1) * stampsize]
                for iy in range(y) for ix in range(x)]
        for i in range(x):
            if numpy.sum(star[-1]) == 0:
                star.pop()
        return star  # a list of psfs

    def image_stack(self, image_list, stampsize, columns):
        # the inverse operation of divide_stamps
        num = len(image_list)
        row_num, c = divmod(num, columns)
        if c != 0:
            row_num += 1
        arr = numpy.zeros((stampsize, row_num * columns * stampsize))
        for i in range(num):
            arr[0:stampsize, i*stampsize:(i+1)*stampsize]=image_list[i]
        arr0 = arr[0:stampsize, 0:columns * stampsize]
        for i in range(1, row_num):
            arr0 = numpy.row_stack((arr0,arr[0:stampsize,i*columns*stampsize:(i+1)*columns*stampsize]))
        return arr0

    def fit(self, star_stamp, noise_stamp, star_data, stampsize,mode=1):
        psf_pool = self.divide_stamps(star_stamp, stampsize)
        noise_pool = self.divide_stamps(noise_stamp, stampsize)
        x = star_data[:, 0]
        y = star_data[:, 1]
        sxx = numpy.sum(x * x)
        syy = numpy.sum(y * y)
        sxy = numpy.sum(x * y)
        sx = numpy.sum(x)
        sy = numpy.sum(y)
        sz = 0.
        szx = 0.
        szy = 0.
        d=1
        rim = self.border(d,stampsize)
        n = numpy.sum(rim)
        for i in range(len(psf_pool)):
            if mode==1:
                pmax = numpy.max(psf_pool[i])
                arr = psf_pool[i]/pmax#[p[0][0]-1:p[0][0]+1,p[1][0]-1:p[1][0]+1])
                conv = self.gaussfilter(arr)
                dy,dx = self.mfpoly(conv)
                psf = ndimage.shift(arr,(24-dy,24-dx),mode='reflect')
            elif mode==2:
                psf = self.pow_spec(psf_pool[i])
                noise = self.pow_spec(noise_pool[i])
                # psf_pnoise = numpy.sum(rim*psf)/n
                # noise_pnoise = numpy.sum(rim*noise)/n
                psf =psf - noise# -psf_pnoise+noise_pnoise
                pmax = (psf[23,24]+psf[24,23]+psf[24,25]+psf[25,24])/4
                #pmax = numpy.sum(psf)
                psf = psf / pmax
                psf[24,24]=1
            else:
                psf = self.psf_align(psf_pool[i])
                psf = psf/numpy.max(psf)
            sz += psf
            szx += psf * x[i]
            szy += psf * y[i]
        a = numpy.zeros((stampsize, stampsize))
        b = numpy.zeros((stampsize, stampsize))
        c = numpy.zeros((stampsize, stampsize))
        co_matr = numpy.array([[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, len(x)]])
        for m in range(stampsize):
            for n in range(stampsize):
                re = numpy.linalg.solve(co_matr, numpy.array([szx[m, n], szy[m, n], sz[m, n]]))
                a[m, n] = re[0]
                b[m, n] = re[1]
                c[m, n] = re[2]
        return a, b, c

    def border(self,edge,size):
        if edge>=size/2. :
            print("Edge must be smaller than half of  the size!")
        else:
            arr = numpy.ones((size,size))
            arr[edge:size-edge,edge:size-edge] = 0.
            return arr

    def set_bins(self,data_array,bin_num,sample=None): # checked 2017-7-9!!!
        # The input must be one dimensional array.
        if sample is not None:
            temp = numpy.sort(data_array)[sample:-sample]
            #temp = numpy.random.choice(data_array,sample,replace=False)
        else:
            temp = data_array
        dat_ma = numpy.max(numpy.abs(temp))
        bins = numpy.linspace(-dat_ma, dat_ma, bin_num + 1)
        bound = numpy.max(numpy.abs(data_array))*1.5
        bins_r = bins[:-1]
        bin_size = bins[1] - bins[0]
        # Because of the bins set up which bases on a sample of the original data,
        # the bins should be extended to include some data points that are out of the bounds.
        bins = numpy.append(-bound,numpy.append(bins[1:-1],bound))
        num_in_bin = numpy.histogram(data_array,bins)[0]
        return bins_r, num_in_bin, bin_size

    def G_bin(self, g, n, u, g_h, mode, bin_num,sample):#checked 2017-7-9!!!
        # mode 1 is for g1
        # mode 2 is for g2
        inverse = range(int(bin_num / 2 - 1), -1, -1)
        if mode==1:
            G_h = g - (n+u)*g_h
        else:
            G_h = g - (n-u)*g_h
        num = self.set_bins(G_h, bin_num,sample=sample )[1]
        n1 = num[0:int(bin_num / 2)]
        n2 = num[int(bin_num / 2):][inverse]
        return numpy.sum((n1 - n2)**2 / (n1 + n2))*0.5


    def fmin_g(self, g, n, u, mode, bin_num, left=-0.05, right=0.05, method=2,sample=100): #checked 2017-7-9!!!
        # model 1 for  g1
        # model 2 for g2
        if method==1:
            def func(g_g):
                return self.G_bin(g, n,u,g_g, mode, bin_num, sample=sample)
            g_h = optimize.fmin(func, [0.], xtol=1.e-8, ftol=1.e-8,maxfun=800, disp=0)[0]
        else:
            iters = 0
            same = 0
            while True:
                templ = left
                tempr = right
                m1 = (left+right)/2.
                m2 = (m1+left)/2.
                m3 = (m1+right)/2.
                fL = self.G_bin(g, n, u, left, mode, bin_num, sample=sample)
                fR = self.G_bin(g, n, u, right, mode, bin_num, sample=sample)
                fm1 = self.G_bin(g, n, u, m1, mode, bin_num, sample=sample)
                fm2 = self.G_bin(g, n, u, m2, mode, bin_num, sample=sample)
                fm3 = self.G_bin(g, n, u, m3, mode, bin_num, sample=sample)
                # print(fL,fm2,fm1,fm3,fR)
                # print(left,m2,m1,m3,right)
                # plt.scatter([left,m2,m1,m3,right],[fL,fm2,fm1,fm3,fR])
                # plt.show()
                if fL>max(fm1,fm2,fm3) and fR>max(fm1,fm2,fm3):
                    if fm1==fm2:
                        left = m2
                        right = m1
                    elif fm1==fm3:
                        left = m1
                        right = m3
                    elif fm2==fm3:
                        left = m2
                        right = m3
                    elif fm1<fm2 and fm1<fm3:
                        left = m2
                        right = m3
                    elif fm2<fm1 and fm2 <fm3:
                        right = m1
                    elif fm3<fm1 and fm3<fm2:
                        left = m1
                elif fR>fm3>=fL:
                    if fL==fm3:
                        right = m3
                    elif fL==fm1:
                        right=m1
                    elif fL==fm2:
                        right = m2
                    elif fm1==fm2:
                        right = right
                    elif fm1<fm2 and fm1<fL:
                        left = m2
                        right = m3
                    elif fm2<fL and fm2<fm1:
                        right = m1
                    elif fL<fm1 and fL <fm2:
                        right = m2
                elif fL>fm2>=fR:
                    if fR==fm2:
                        left = m2
                    elif fR==fm1:
                        left = m1
                    elif fR==fm3:
                        left = m3
                    elif fm1<fR and fm1 <fm3:
                        left = m2
                        right = m3
                    elif fm3<fm1 and fm3<fR:
                        left = m1
                    elif fR < fm1 and fR <fm3:
                        left = m3

                if abs(left-right)<1.e-5:
                    g_h = (left+right)/2.
                    break
                iters+=1
                if left==templ and right==tempr:
                    same+=1
                if iters>10 and same>3 or iters>15:
                    g_h = (left+right)/2.
                    break
                #print(left,right,abs(left-right))
        # fitting
        g_range = numpy.linspace(g_h-0.01,g_h+0.01,11)
        xi2 = numpy.array([self.G_bin(g,n,u,g_hat,mode,bin_num,sample=sample) for g_hat in g_range])
        gg4 = numpy.sum(g_range ** 4)
        gg3 = numpy.sum(g_range ** 3)
        gg2 = numpy.sum(g_range ** 2)
        gg1 = numpy.sum(g_range)
        xigg2 = numpy.sum(xi2 * (g_range ** 2))
        xigg1 = numpy.sum(xi2 * g_range)
        xigg0 = numpy.sum(xi2)
        cov = numpy.linalg.inv(numpy.array([[gg4, gg3, gg2], [gg3, gg2, gg1], [gg2, gg1, len(g_range)]]))
        paras = numpy.dot(cov, numpy.array([xigg2, xigg1, xigg0]))
        g_sig = numpy.sqrt(1/2./paras[0])
        # x = numpy.linspace(g_h-0.005,g_h+0.005,50)
        # plt.scatter(g_range,xi2,linewidths=1)
        # plt.plot(x,paras[0]*x**2+paras[1]*x+paras[2])
        # plt.show()
        return g_h,g_sig

    def ellip_plot(self, ellip, coordi, lent, width, title, mode=1,path=None,show=True):
        e1 = ellip[:, 0]
        e2 = ellip[:, 1]
        e = numpy.sqrt(e1 ** 2 + e2 ** 2)
        scale = numpy.mean(1 / e)
        x = coordi[:, 0]
        y = coordi[:, 1]
        if mode== 1:
            dx = lent * e1 / e / 2
            dy = lent * e2 / e / 2
        else:
            dx = scale * lent * e1 / e / 2
            dy = scale * lent * e2 / e / 2
        x1 = x + dx
        x2 = x - dx
        y1 = y + dy
        y2 = y - dy

        norm = plt.Normalize(vmin=numpy.min(e), vmax=numpy.max(e))
        cmap = plt.get_cmap('YlOrRd')
        fig = plt.figure(figsize=(20,10))
        plt.axes().set_aspect('equal', 'datalim')
        for i in range(len(x)):
            cl = cmap(norm(e[i]))
            plt.plot([y1[i], y2[i]], [x1[i], x2[i]], color=cl, linewidth=width)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        plt.colorbar(sm)
        plt.title(title,fontsize = 18)
        if path is not None:
            plt.savefig(path)
        if show is True:
            plt.show()
