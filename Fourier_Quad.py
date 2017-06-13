import numpy
from numpy import fft
from scipy.optimize import least_squares
from scipy.optimize import fmin_cg
from scipy import ndimage, signal
import copy


class Fourier_Quad:

    def pow_spec(self, image):
        image_ps = fft.fftshift((numpy.abs(fft.fft2(image)))**2)
        return image_ps

    def shear_est(self, gal, psf, imagesize, background_noise=None, F=True):
        x = imagesize#int(gal.shape[0]) #resolution of the image
        my, mx = numpy.mgrid[0:x, 0:x]
        gal_ps = self.pow_spec(gal)

        if background_noise is not None: # to deduct the noise
            d=1
            nbg = self.pow_spec(background_noise)
            rim = self.border(d,x)
            n    = numpy.sum(rim)
            gal_pnoise = numpy.sum(gal_ps*rim)/n               #the Possion noise of galaxy image
            nbg_pnoise =numpy.sum(nbg*rim)/n                   #the  Possion noise of background noise image
            gal_ps = (gal_ps - gal_pnoise) - (nbg - nbg_pnoise)

        if F == True:
            psf = psf
        else:
            psf = self.pow_spec(psf)

        hlr = self.get_radius(psf,2.)
        wb,beta = self.wbeta(hlr,x)
        maxi = numpy.max(wb)
        idx = wb < maxi / 100000.
        wb[idx] = 0.
        maxi = numpy.max(psf)
        idx = psf < maxi / 100000.
        psf[idx] = 1.

        tk  = wb/ psf * gal_ps
        alpha = 2.*numpy.pi/x
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

    def cre_psf(self, psf_scale, imagesize, psf="GAUSS", x=0, y=0):
        xx = numpy.linspace(0, imagesize - 1, imagesize)
        mx, my = numpy.meshgrid(xx, xx)
        if psf is 'GAUSS':
            arr = numpy.exp(-((mx -imagesize/2.+x)**2+(my-imagesize/2.+y)**2)/2./psf_scale**2)
            return arr

        if psf is 'Moffat':
            hstep = 3*psf_scale-numpy.sqrt((mx-imagesize/2.+x)**2+(my-imagesize/2.+y)**2)
            idx = hstep < 0.
            hstep[idx] = 0.
            idx = hstep != 0.
            hstep[idx] = 1.
            arr = (1+((mx-imagesize/2.+x)**2+(my-imagesize/2.+y)**2)/psf_scale**2)**(-3.5)*hstep
            return arr

    def gethlr(self, image, disp=1):  # the sum of the image array must be positive!!! or it will fail
        size = image.shape[0]
        ra = size / 2 * 10
        half = 0.
        maxi = numpy.max(image)
        flux = numpy.sum(image)
        y, x = numpy.where(image == maxi)
        y = int(y[0])
        x = int(x[0])
        yy, xx = numpy.mgrid[0:size, 0:size]
        xx = numpy.abs(xx - x)
        yy = numpy.abs(yy - y)
        if flux > 0:
            for r in range(10, ra, 2):
                if half < flux / 2.:
                    cir = xx ** 2 + yy ** 2 - (r / 10.) ** 2
                    idx = cir <= 0.
                    cir[idx] = 1.
                    idx = cir > 1.
                    cir[idx] = 0.
                    half = numpy.sum(cir * image)
                    hlr = r / 10.
                else:
                    break
        else:
            if disp == 1:
                print ("Failed! The sum of image array is non-positive!")
            hlr = 20.
        return hlr

    def get_radius(self, image, scale):  # get the radius of the flux descends to the maximum/scale
        arr = copy.deepcopy(image)
        maxi = numpy.max(arr)
        y, x = numpy.where(arr == maxi)
        idx = arr < maxi / scale
        arr[idx] = 0.
        idx = arr > 0.
        arr[idx] = 1.
        pool = []
        pool.append((int(y[0]), int(x[0])))

        def check(x):  # x is a two components  tuple -like input
            if (x[0] - 1, x[1]) not in pool and arr[x[0] - 1, x[1]] == 1:
                pool.append((x[0] - 1, x[1]))
            if (x[0] + 1, x[1]) not in pool and arr[x[0] + 1, x[1]] == 1:
                pool.append((x[0] + 1, x[1]))
            if (x[0], x[1] + 1) not in pool and arr[x[0], x[1] + 1] == 1:
                pool.append((x[0], x[1] + 1))
            if (x[0], x[1] - 1) not in pool and arr[x[0], x[1] - 1] == 1:
                pool.append((x[0], x[1] - 1))
            return len(pool)

        while True:
            for cor in pool:
                num0 = len(pool)
                num1 = check(cor)
            if num0 == num1:
                break
        return numpy.sqrt(len(pool) / numpy.pi)

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

    def get_centroid(self, image):
        m0 = numpy.sum(image)
        y, x = numpy.mgrid[0:image.shape[0], 0:image.shape[1]]
        x0 = numpy.sum(x * image) / m0
        y0 = numpy.sum(y * image) / m0
        return x0, y0

    def psf_align(self, image):
        imagesize = image.shape[0]
        arr = self.move(image, 0, 0)
        image_f = fft.fft2(arr)
        yy, xx = numpy.mgrid[0:48, 0:48]
        xx = numpy.mod(xx + 24, 48) - 24
        yy = numpy.mod(yy + 24, 48) - 24
        fk = numpy.abs(image_f) ** 2
        line = numpy.sort(numpy.array(fk.flat))
        idx = fk < 4 * line[imagesize ** 2/ 2]
        fk[idx] = 0
        weight = fk / line[imagesize ** 2/ 2]
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

    def mfpoly(self, psfimage):
        p = numpy.where(psfimage == numpy.max(psfimage))
        xr, yr = p[1][0], p[0][0]
        xp, yp = numpy.mgrid[xr - 1:xr + 2, yr - 1:yr + 2]
        patchim = psfimage[xr - 1:xr + 2, yr - 1:yr + 2]
        zz = patchim.reshape(9)
        xx = xp.reshape(9)
        yy = yp.reshape(9)
        xy = xx * yy
        x2 = xx * xx
        y2 = yy * yy
        mA = numpy.array([numpy.ones_like(zz), xx, yy, x2, xy, y2]).T
        cov = numpy.linalg.inv(numpy.dot(mA.T, mA))
        a, b, c, d, e, f = numpy.dot(cov, numpy.dot(mA.T, zz))
        coeffs = numpy.array([[2.0 * d, e], [e, 2.0 * f]])
        mult = numpy.array([-b, -c])
        xc, yc = numpy.dot(numpy.linalg.inv(coeffs), mult)
        return xc, yc

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

    def image_stack(self, image_list, stampsize, columns):  # the inverse operation of divide_stamps
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

    def fit(self, star_stamps, star_noise, star_data, stampsize):
        psf_pool = self.divide_stamps(star_stamps, stampsize)
        noise_pool = self.divide_stamps(star_noise, stampsize)
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
            # arr = psf_fit_pool[i]/numpy.max(psf_fit_pool[i])
            # conv = self.gaussfilter(arr)
            # dx,dy = self.mfpoly(conv)
            # psf = ndimage.shift(arr,(24-dx,24-dy))
            psf = self.pow_spec(psf_pool[i])
            noise = self.pow_spec(noise_pool[i])
            psf_pnoise = numpy.sum(rim*psf)/n
            noise_pnoise = numpy.sum(rim*noise)/n
            psf =psf - noise -psf_pnoise+noise_pnoise
            pmax = (psf[23,24]+psf[24,23]+psf[24,25]+psf[25,24])/4
            psf = psf / pmax
            psf[24,24]=1
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

    def smooth(self, image, edge, ra):
        xl = image.shape[0]
        image = numpy.log10(image)
        arr = copy.deepcopy(image)
        mx, my = numpy.mgrid[0: 2 * ra + 1, 0: 2 * ra + 1]
        nx = mx.flatten()
        x = numpy.delete(nx, (0, 2 * ra, 4 * ra ** 2 + 2 * ra, 4 * ra ** 2 + 4 * ra))
        ny = my.flatten()
        y = numpy.delete(ny, (0, 2 * ra, 4 * ra ** 2 + 2 * ra, 4 * ra ** 2 + 4 * ra))
        x4 = numpy.sum(x ** 4)
        x2y2 = numpy.sum(x ** 2 * (y ** 2))
        x3y = numpy.sum(x ** 3 * y)
        x3 = numpy.sum(x ** 3)
        x2y = numpy.sum(x ** 2 * y)
        x2 = numpy.sum(x ** 2)
        y4 = numpy.sum(y ** 4)
        xy3 = numpy.sum(x * (y ** 3))
        xy2 = numpy.sum(y ** 2 * x)
        y3 = numpy.sum(y ** 3)
        y2 = numpy.sum(y ** 2)
        xy = numpy.sum(x * y)
        sx = numpy.sum(x)
        sy = numpy.sum(y)
        n = len(y)
        for i in range(edge, xl - edge):
            for k in range(edge, xl - edge):
                farr = image[i - ra:i + ra+1, k-ra:k + ra + 1].flatten()
                farr = numpy.delete(farr, (0, 2*ra, 4*ra**2+2*ra,4*ra**2+4*ra))
                szx2 = numpy.sum(farr * (x ** 2))
                szy2 = numpy.sum(farr * (y ** 2))
                szxy = numpy.sum(farr * x * y)
                szx = numpy.sum(farr * x)
                szy = numpy.sum(farr * y)
                sz = numpy.sum(farr)
                co_matr = numpy.array(
                    [[x4, x2y2, x3y, x3, x2y, x2], [x2y2, y4, xy3, xy2, y3, y2],
                     [x3y, xy3, x2y2, x2y, xy2, xy],[x3, xy2, x2y, x2, xy, sx],
                     [x2y, y3, xy2, xy, y2, sy], [x2, y2, xy, sx, sy, n]])
                val = numpy.array([szx2, szy2, szxy, szx, szy, sz])
                re = numpy.linalg.solve(co_matr, val)
                arr[i, k] = ra**2*(re[0]+re[1]+re[2])+(re[3]+re[4])*ra+re[5]
        return numpy.power(10, arr)

    def border(self,edge,size):
        if edge>=size/2. :
            print("Edge must be smaller than half of  the size!")
        else:
            arr = numpy.ones((size,size))
            arr[edge:size-edge,edge:size-edge] = 0.
            return arr

    def set_bins(self,data_array,bin_num,model=1,sym=False):# The input must be one dimensional array.(1,n)
        array=copy.copy(data_array)
        if model ==1:   #set up bins according to the  maximum and minimum of the data
            mi = numpy.min(array)
            ma = numpy.max(array)
            if sym==True: # if the data are symmetric respect to zero
                radius = max(numpy.abs(mi),ma)
                bins0 = numpy.linspace(-radius,radius,bin_num+1)
            else:
                bins0 = numpy.linspace(mi,ma,bin_num+1)
            bin_size = bins0[1]-bins0[0]
            bins = numpy.delete(bins0,-1)
            arr = numpy.digitize(array,bins)
            tag = numpy.linspace(1,bin_num,bin_num)
            num_in_bin = numpy.zeros((bin_num))
            for i in range(bin_num):
                idx = arr ==tag[i]
                num_in_bin[i] = len(arr[idx])
        else:  # set up bins for the middle part  of the data which account for 80% ,then add the data to the two sides of the bins
            datalen = len(array)
            if datalen> 10000:
                idx = numpy.random.randint(0, datalen, 5000)
            else:
                idx = numpy.random.randint(0, datalen, int(len(array) / 10))
            temp = array[idx]
            ma = numpy.max(temp)
            bins = []
            for i in range(int(bin_num / 2)):
                bins.append(ma * i * 2 / bin_num)
                if i != 0:
                    bins.append(-ma * i * 2 / bin_num)
            bins = numpy.sort(bins)
            arr = numpy.digitize(array, bins)
            num_in_bin = numpy.zeros((bin_num))
            tag = numpy.arange(0, bin_num)
            for i in range(bin_num):
                idx = arr == tag[i]
                num_in_bin[i] = len(arr[idx])
            bin_size = bins[1] - bins[0]
            bins = numpy.append(bins[0] - bin_size, bins)
        return bins, num_in_bin, bin_size








