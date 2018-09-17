import platform
if platform.system() == 'Linux':
    import matplotlib
    matplotlib.use('Agg')
import numpy
from numpy import fft
from scipy.optimize import fmin_cg
from scipy import ndimage, signal
import copy
import matplotlib.pyplot as plt
import tool_box
import time
import datetime


class Fourier_Quad:

    def __init__(self, size, seed):
        self.rng = numpy.random.RandomState(seed)
        self.size = size
        self.alpha = (2.*numpy.pi/size)**4
        self.my = numpy.mgrid[0: size, 0: size][0] - size/2.
        self.mx = numpy.mgrid[0: size, 0: size][1] - size/2.
        self.kx2 = self.mx*self.mx
        self.ky2 = self.my*self.my
        self.kxy = self.mx*self.my
        self.k2 = self.kx2+self.ky2
        self.k4 = self.k2*self.k2
        self.mn1 = (-0.5)*(self.kx2 - self.ky2)
        self.mn2 = -self.mx*self.my
        self.mn4 = self.k4 - 8*self.kx2*self.ky2
        self.mn5 = self.kxy*(self.kx2 - self.ky2)
        self.rim = self.border(1)
        self.flux2 = -1.
        self.hlr = 0

    def draw_noise(self, mean, sigma):
        noise_img = self.rng.normal(loc=mean, scale=sigma, size=self.size * self.size).reshape(self.size, self.size)
        return noise_img

    def pow_spec(self, image):
        image_ps = fft.fftshift((numpy.abs(fft.fft2(image)))**2)
        return image_ps

    def shear_est(self, gal_image, psf_image, noise=None, F=False):
        gal_ps = self.pow_spec(gal_image)
        # gal_ps = tool_box.smooth(gal_ps,self.size)
        if noise is not None:
            nbg = self.pow_spec(noise)
            self.flux2 = numpy.sqrt(gal_ps[int(self.size/2), int(self.size/2)]/numpy.sum(self.rim*gal_ps))
            # nbg = tool_box.smooth(nbg,self.size)
            # rim = self.border(2, size)
            # n = numpy.sum(rim)
            # gal_pn = numpy.sum(gal_ps*rim)/n                # the Possion noise of galaxy image
            # nbg_pn = numpy.sum(nbg*rim)/n                   # the  Possion noise of background noise image
            gal_ps = gal_ps - nbg# + nbg_pn - gal_pn

        if F:
            psf_ps = psf_image
        else:
            psf_ps = self.pow_spec(psf_image)
        # hlr = self.get_radius_new(psf_ps, 2)[0]
        wb, beta = self.wbeta(self.hlr)
        maxi = numpy.max(psf_ps)
        idx = psf_ps < maxi / 10000.
        wb[idx] = 0
        psf_ps[idx] = 1.
        tk = wb/psf_ps * gal_ps

        # ky, kx = self.my, self.mx
        # #
        # kx2 = kx*kx
        # ky2 = ky*ky
        # kxy = kx*ky
        # k2 = kx2 + ky2
        # k4 = k2*k2
        # mn1 = (-0.5)*(kx2 - ky2)    # (-0.5)*(kx**2 - ky**2)
        # mn2 = -kxy                  # -kx*ky
        # mn3 = k2 - 0.5*beta**2*k4   # kx**2 + ky**2 - 0.5*beta**2*(kx**2 + ky**2)**2
        # mn4 = k4 - 8*kx2*ky2        # kx**4 - 6*kx**2*ky**2 + ky**4
        # mn5 = kxy*(kx2 - ky2)       # kx**3*ky - kx*ky**3

        # mn1 = self.mn1
        # mn2 = self.mn2
        mn3 = self.k2 - 0.5*beta**2*self.k4
        # mn4 = self.mn4
        # mn5 = self.mn5

        mg1 = numpy.sum(self.mn1 * tk)*self.alpha
        mg2 = numpy.sum(self.mn2 * tk)*self.alpha
        mn  = numpy.sum(mn3 * tk)*self.alpha
        mu  = numpy.sum(self.mn4 * tk)*(-0.5*beta**2)*self.alpha
        mv  = numpy.sum(self.mn5 * tk)*(-2.*beta**2)*self.alpha

        return mg1, mg2, mn, mu, mv

    def wbeta(self, radius):
        w_temp = numpy.exp(-(self.mx**2 + self.my**2)/radius**2)
        return w_temp, 1./radius

    def ran_pos(self, num, radius, g=None):
        xy_coord = numpy.zeros((2, num))
        theta = self.rng.uniform(0., 2 * numpy.pi, num)
        xn = numpy.cos(theta)
        yn = numpy.sin(theta)
        x = 0.
        y = 0.
        for n in range(num):
            x += xn[n]
            y += yn[n]
            if x * x + y * y > radius ** 2:
                x = xn[n]
                y = yn[n]
            xy_coord[0, n] = x
            xy_coord[1, n] = y

        xy_coord[0] = xy_coord[0] - numpy.mean(xy_coord[0])
        xy_coord[1] = xy_coord[1] - numpy.mean(xy_coord[1])

        if g:
            sheared = numpy.dot(numpy.array(([(1 + g[0]), g[1]], [g[1], (1 - g[0])])), xy_coord)
            return xy_coord, sheared
        else:
            return xy_coord

    def rotate(self, pos, theta):
        rot_matrix = numpy.array([[numpy.cos(theta), -numpy.sin(theta)], [numpy.sin(theta), numpy.cos(theta)]])
        return numpy.dot(rot_matrix, pos)

    def shear(self, pos, g1, g2):
        return numpy.dot(numpy.array(([(1+g1), g2], [g2, (1-g1)])), pos)

    def convolve_psf(self, pos, psf_scale, flux=1., psf="GAUSS"):
        x = pos.shape[1]
        arr = numpy.zeros((self.size, self.size))

        if psf == 'GAUSS':
            factor = flux/2/numpy.pi/psf_scale**2
            for i in range(x):
                arr += factor*numpy.exp(-((self.mx-pos[0, i])**2+(self.my-pos[1, i])**2)/2./psf_scale**2)

        elif psf == "Moffat":
            r_scale_sq = 9
            m = 3.5
            factor = flux/(numpy.pi*psf_scale**2*((1. + r_scale_sq)**(1.-m) - 1.)/(1.-m))
            for l in range(x):
                rsq = ((self.mx-pos[0, l])**2+(self.my-pos[1, l])**2)/psf_scale**2
                idx = rsq > r_scale_sq
                pfunction = factor*(1. + rsq)**(-m)
                pfunction[idx] = 0.
                arr += pfunction
        return arr

    def cre_psf(self, psf_scale, flux=1., model="GAUSS"):
        if model is 'GAUSS':
            factor = flux*1./2/numpy.pi/psf_scale**2
            arr = factor*numpy.exp(-(self.mx**2 + self.my**2)/2./psf_scale**2)
            return arr

        elif model is 'Moffat':
            r_scale_sq = 9
            m = 3.5
            factor = flux*1./(numpy.pi*psf_scale**2*((1. + r_scale_sq)**(1.-m) - 1.)/(1.-m))
            rsq = (self.mx**2 + self.my**2) / psf_scale**2
            idx = rsq > r_scale_sq
            rsq[idx] = 0.
            arr = factor*(1. + rsq)**(-m)
            arr[idx] = 0.
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
        self.hlr = numpy.sqrt(len(half_radi_pool) / numpy.pi)
        return numpy.sqrt(len(half_radi_pool) / numpy.pi)

    def get_radius_new(self, image, scale):
        # get the radius of the flux descends to the maximum/scale
        radi_arr = copy.copy(image)
        maxi = numpy.max(radi_arr)
        y, x = numpy.where(radi_arr == maxi)
        idx = radi_arr < maxi / scale
        radi_arr[idx] = 0.
        half_radius_pool = []
        flux = []

        def detect(mask, ini_y, ini_x, signal, signal_val, size):
            if mask[ini_y, ini_x] > 0:
                signal.append((ini_y, ini_x))
                signal_val.append(mask[ini_y, ini_x])
                mask[ini_y, ini_x] = 0
                for cor in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                    if -1 < ini_y + cor[0] < size and -1 < ini_x + cor[1] < size \
                            and mask[ini_y + cor[0], ini_x + cor[1]] > 0:
                        detect(mask, ini_y + cor[0], ini_x + cor[1], signal, signal_val, size)
            return signal, signal_val

        half_radius_pool, flux = detect(radi_arr, y[0], x[0], half_radius_pool, flux, self.size)
        self.hlr = numpy.sqrt(len(half_radius_pool)/numpy.pi)
        return numpy.sqrt(len(half_radius_pool)/numpy.pi), half_radius_pool, numpy.sum(flux), maxi, (y, x)

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

    def get_centroid(self, image, filt=False, radius=2):
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

    def segment(self, image):
        shape = image.shape
        y = int(shape[0] / self.size)
        x = int(shape[1] / self.size)
        star = [image[iy * self.size:(iy + 1) * self.size, ix * self.size:(ix + 1) * self.size]
                for iy in range(y) for ix in range(x)]
        for i in range(x):
            if numpy.sum(star[-1]) == 0:
                star.pop()
        return star  # a list of psfs

    def stack(self, image_array, columns):
        # the inverse operation of divide_stamps
        # the image_array is a three dimensional array of which the length equals the number of the stamps
        num = len(image_array)
        row_num, c = divmod(num, columns)
        if c != 0:
            row_num += 1
        arr = numpy.zeros((row_num*self.size, columns * self.size))
        for j in range(row_num):
            for i in range(columns):
                tag = i + j * columns
                if tag > num - 1:
                    break
                arr[j*self.size:(j+1)*self.size, i*self.size:(i+1)*self.size] = image_array[tag]
        return arr

    def fit(self, star_stamp, noise_stamp, star_data, mode=1):
        psf_pool = self.segment(star_stamp)
        noise_pool = self.segment(noise_stamp)
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
        rim = self.border(d)
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
        a = numpy.zeros((self.size, self.size))
        b = numpy.zeros((self.size, self.size))
        c = numpy.zeros((self.size, self.size))
        co_matr = numpy.array([[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, len(x)]])
        for m in range(self.size):
            for n in range(self.size):
                re = numpy.linalg.solve(co_matr, numpy.array([szx[m, n], szy[m, n], sz[m, n]]))
                a[m, n] = re[0]
                b[m, n] = re[1]
                c[m, n] = re[2]
        return a, b, c

    def border(self, edge):
        if edge >= self.size/2.:
            print("Edge must be smaller than half of  the size!")
        else:
            arr = numpy.ones((self.size, self.size))
            arr[edge: self.size - edge, edge: self.size - edge] = 0.
            return arr

    def set_bin(self, data, bin_num, abs_sort=True):
        if abs_sort:
            temp_data = numpy.sort(numpy.abs(data))
        else:
            temp_data = numpy.sort(data[data>0])
        bin_size = len(temp_data) / bin_num * 2
        bins = numpy.array([temp_data[int(i * bin_size)] for i in range(1, int(bin_num / 2))])
        bins = numpy.sort(numpy.append(numpy.append(-bins, [0.]), bins))
        bound = numpy.max(numpy.abs(data)) * 100.
        bins = numpy.append(-bound, numpy.append(bins, bound))
        return bins

    def G_bin(self, g, nu, g_h, bins, ig_num):  # checked 2017-7-9!!!
        r"""
        to calculate the symmetry the shear estimators
        :param g: estimators from Fourier quad, 1-D numpy array
        :param nu: N + U for g1, N - U for g2, 1-D numpy array
        :param g_h: pseudo shear (guess)
        :param bins: bin of g for calculation of the symmetry, 1-D numpy array
        :param ig_num: the number of inner grid of bin to be neglected
        :return: chi square
        """
        bin_num = len(bins) - 1
        inverse = range(int(bin_num / 2 - 1), -1, -1)
        G_h = g - nu * g_h
        num = numpy.histogram(G_h, bins)[0]
        n1 = num[0:int(bin_num / 2)]
        n2 = num[int(bin_num / 2):][inverse]
        xi = (n1 - n2) ** 2 / (n1 + n2)
        return numpy.sum(xi[:len(xi)-ig_num]) * 0.5

    def G_bin2d(self, mgs, mnus, g_corr, bins, resample=1, ig_nums=0):
        r"""
        to calculate the symmetry of two sets of shear estimators
        :param mgs: a two components list, two 1-D numpy arrays, contains two sets of
                    first shear estimators of Fourier quad ,
        :param mnus: a two components list, two 1-D numpy arrays, contains two sets of
                    shear estimators of Fourier quad (N,U), N + U for g1, N - U for g2
        :param g_corr: float, the correlation between the two sets of shear
        :param bins: a two components list, two 1-D numpy arrays, contains two bins for
                    the shear estimator
        :param resample: repeat the calculation for one shear PDF to reduct the noise
        :param ig_nums: a two components list, [num1, num2], the numbers of the inner grid
                        of each bin to be neglected (developing)
        :return: chi square
        """
        # half of the bin number
        ny, nx = int((len(bins[0]) - 1)/2), int((len(bins[1]) - 1)/2)
        # print("ny, nx:", ny,nx)
        chi_sq = 0
        mu = [0,0]
        cov = [[0.07,g_corr],[g_corr,0.07]]
        data_len = len(mgs[0])
        for i in range(resample):
            t1 = time.time()
            g_distri = numpy.random.multivariate_normal(mu,cov,data_len)
            t2 = time.time()
            mg1 = mgs[0] - mnus[0]*g_distri[:,0]
            mg2 = mgs[1] - mnus[1]*g_distri[:,1]
            num_arr = numpy.histogram2d(mg1, mg2, bins)[0]
            # | arr_1 | arr_2 |
            # | arr_3 | arr_4 |
            # chi square = 0.2*SUM[(arr_2 + arr_3 - arr_1 - arr_4)**2/(arr_1 + arr_2 + arr_3 + arr_4)]
            arr_1 = num_arr[0:ny, 0:nx][:,range(ny-1,-1,-1)]
            arr_2 = num_arr[0:ny, nx:2*nx]
            arr_3 = num_arr[ny:2*ny, 0:nx][range(ny-1,-1,-1)][:,range(nx-1,-1,-1)]
            arr_4 = num_arr[ny:2*ny, nx:2*nx][range(ny-1,-1,-1)]
            chi_sq += 0.5 * numpy.sum(((arr_2 + arr_3 - arr_1 - arr_4) ** 2) / (arr_1 + arr_2 + arr_3 + arr_4))
            t3 = time.time()
            # print("%d resample %d %.2f %.2f" %(self.counter,i, t2-t1, t3-t2))
        return chi_sq/resample, arr_1, arr_2, arr_3, arr_4

    def fmin_g2d(self, mgs, mnus, bin_num, ig_nums=0, left=-0.001, right=0.001,pic_path=False):
        r"""
        to find the true correlation between two sets of galaxies
        :param mgs: a two components list, two 1-D numpy arrays, contains two sets of
                    first shear estimators of Fourier quad ,
        :param mnus: a two components list, two 1-D numpy arrays, contains two sets of
                    shear estimators of Fourier quad (N,U), N + U for g1, N - U for g2
        :param bin_num: [num1, num2] for [mg1, mg2]
        :param ig_nums: a two components list, [num1, num2], the numbers of the inner grid
                        of each bin to be neglected
        :param pic_path: the path for saving the chi square picture
        :param left: initial guess of the correlation
        :param right: initial guess of the correlation
        :return: the correlation between mg1 and mg2
        """
        bins = [self.set_bin(mgs[0], bin_num), self.set_bin(mgs[1], bin_num)]
        same = 0
        iters = 0
        # m1 chi square & left & left chi square & right & right chi square
        records = numpy.zeros((15, 5))
        while True:
            templ = left
            tempr = right
            m1 = (left + right) / 2.
            m2 = (m1 + left) / 2.
            m3 = (m1 + right) / 2.
            fL = self.G_bin2d(mgs, mnus, left, bins, ig_nums=ig_nums)
            fm2 = self.G_bin2d(mgs, mnus, m2, bins, ig_nums=ig_nums)
            fm1 = self.G_bin2d(mgs, mnus, m1, bins, ig_nums=ig_nums)
            fm3 = self.G_bin2d(mgs, mnus, m3, bins, ig_nums=ig_nums)
            fR = self.G_bin2d(mgs, mnus, right, bins, ig_nums=ig_nums)
            values = [fL, fm2, fm1, fm3, fR]
            print("iteration \n",iters,values,"\n", left, m2, m1, m3, right)
            records[iters, ] = fm1, left, fL, right, fR
            if max(values) - min(values) < 20:
                print("BREAK ", iters, left, right, abs(left - right))
                break
            if fL > max(fm1, fm2, fm3) and fR > max(fm1, fm2, fm3):
                if fm1 == fm2:
                    left = m2
                    right = m1
                elif fm1 == fm3:
                    left = m1
                    right = m3
                elif fm2 == fm3:
                    left = m2
                    right = m3
                elif fm1 < fm2 and fm1 < fm3:
                    left = m2
                    right = m3
                elif fm2 < fm1 and fm2 < fm3:
                    right = m1
                elif fm3 < fm1 and fm3 < fm2:
                    left = m1
            elif fR > fm3 >= fL:
                if fL == fm3:
                    right = m3
                elif fL == fm1:
                    right = m1
                elif fL == fm2:
                    right = m2
                elif fm1 == fm2:
                    right = right
                elif fm1 < fm2 and fm1 < fL:
                    left = m2
                    right = m3
                elif fm2 < fL and fm2 < fm1:
                    right = m1
                elif fL < fm1 and fL < fm2:
                    right = m2
            elif fL > fm2 >= fR:
                if fR == fm2:
                    left = m2
                elif fR == fm1:
                    left = m1
                elif fR == fm3:
                    left = m3
                elif fm1 < fR and fm1 < fm3:
                    left = m2
                    right = m3
                elif fm3 < fm1 and fm3 < fR:
                    left = m1
                elif fR < fm1 and fR < fm3:
                    left = m3
                elif fm1 == fm3:
                    left = m1
                    right = m3

            if abs(left-right) < 1.e-5:
                g_h = (left+right)/2.
                print("BREAK ",iters, left, right, abs(left - right))
                break
            iters += 1
            if left == templ and right == tempr:
                same += 1
            if iters > 12 and same > 2 or iters > 14:
                g_h = (left+right)/2.
                print("BREAK ",iters,left, right, abs(left - right))
                break
        # fitting
        # left_x2 = numpy.min(numpy.abs(records[:iters, 2] - fm1 - 30))
        # label_l = numpy.where(left_x2 == numpy.abs(records[:iters, 2] - fm1 - 30))[0]
        # if len(label_l > 1):
        #     label_l = label_l[0]
        #
        # right_x2 = numpy.min(numpy.abs(records[:iters, 4] - fm1 - 30))
        # label_r = numpy.where(right_x2 == numpy.abs(records[:iters, 4] - fm1 - 30))[0]
        # if len(label_r > 1):
        #     label_r = label_r[0]
        #
        # if left_x2 > right_x2:
        #     right = records[label_l, 3]
        #     left = 2 * m1 - right
        # else:
        #     left = records[label_r, 1]
        #     right = 2 * m1 - left

        fit_range = numpy.linspace(left, right, 20)
        chi_sq = [self.G_bin2d(mgs, mnus, fit_range[i], bins, ig_nums=ig_nums) for i in range(len(fit_range))]
        coeff = tool_box.fit_1d(fit_range, chi_sq, 2, "scipy")
        corr_sig = numpy.sqrt(1 / 2. / coeff[2])
        g_corr = -coeff[1] / 2. / coeff[2]
        print("Fitting finish", left, right,g_corr, corr_sig)
        plt.scatter(fit_range,chi_sq)
        plt.plot(fit_range, coeff[0]+coeff[1]*fit_range+coeff[2]*fit_range**2)
        plt.show()
        return g_corr, corr_sig

    def fmin_g(self, g, nu, bin_num, ig_num=0, pic_path=False, left=-0.1, right=0.1):  # checked 2017-7-9!!!
        # nu = N + U for g1
        # nu = N - U for g2
        # temp_data = numpy.sort(numpy.abs(g))[:int(len(g)*0.99)]
        # bin_size = len(temp_data)/bin_num*2
        # bins = numpy.array([temp_data[int(i*bin_size)] for i in range(1, int(bin_num / 2))])
        # bins = numpy.sort(numpy.append(numpy.append(-bins, [0.]), bins))
        # bound = numpy.max(numpy.abs(g)) * 100.
        # bins = numpy.append(-bound, numpy.append(bins, bound))
        bins = self.set_bin(g, bin_num)
        same = 0
        iters = 0
        # m1 chi square & left & left chi square & right & right chi square
        records = numpy.zeros((15, 5))
        while True:
            templ = left
            tempr = right
            m1 = (left + right) / 2.
            m2 = (m1 + left) / 2.
            m3 = (m1 + right) / 2.
            fL = self.G_bin(g, nu, left, bins, ig_num)
            fR = self.G_bin(g, nu, right, bins, ig_num)
            fm1 = self.G_bin(g, nu, m1, bins, ig_num)
            fm2 = self.G_bin(g, nu, m2, bins, ig_num)
            fm3 = self.G_bin(g, nu, m3, bins, ig_num)
            values = [fL, fm2, fm1, fm3, fR]
            points = [left, m2, m1, m3, right]
            records[iters, ] = fm1, left, fL, right, fR
            if max(values) < 30:
                temp_left = left
                temp_right = right
            if fL > max(fm1, fm2, fm3) and fR > max(fm1, fm2, fm3):
                if fm1 == fm2:
                    left = m2
                    right = m1
                elif fm1 == fm3:
                    left = m1
                    right = m3
                elif fm2 == fm3:
                    left = m2
                    right = m3
                elif fm1 < fm2 and fm1 < fm3:
                    left = m2
                    right = m3
                elif fm2 < fm1 and fm2 < fm3:
                    right = m1
                elif fm3 < fm1 and fm3 < fm2:
                    left = m1
            elif fR > fm3 >= fL:
                if fL == fm3:
                    right = m3
                elif fL == fm1:
                    right = m1
                elif fL == fm2:
                    right = m2
                elif fm1 == fm2:
                    right = right
                elif fm1 < fm2 and fm1 < fL:
                    left = m2
                    right = m3
                elif fm2 < fL and fm2 < fm1:
                    right = m1
                elif fL < fm1 and fL < fm2:
                    right = m2
            elif fL > fm2 >= fR:
                if fR == fm2:
                    left = m2
                elif fR == fm1:
                    left = m1
                elif fR == fm3:
                    left = m3
                elif fm1 < fR and fm1 < fm3:
                    left = m2
                    right = m3
                elif fm3 < fm1 and fm3 < fR:
                    left = m1
                elif fR < fm1 and fR < fm3:
                    left = m3
                elif fm1 == fm3:
                    left = m1
                    right = m3

            if abs(left-right) < 1.e-5:
                g_h = (left+right)/2.
                break
            iters += 1
            if left == templ and right == tempr:
                same += 1
            if iters > 12 and same > 2 or iters > 14:
                g_h = (left+right)/2.
                break
                # print(left,right,abs(left-right))
        # fitting
        left_x2 = numpy.min(numpy.abs(records[:iters, 2] - fm1 - 20))
        label_l = numpy.where(left_x2 == numpy.abs(records[:iters, 2] - fm1 - 20))[0]
        if len(label_l > 1):
            label_l = label_l[0]

        right_x2 = numpy.min(numpy.abs(records[:iters, 4] - fm1 - 20))
        label_r = numpy.where(right_x2 == numpy.abs(records[:iters, 4] - fm1 - 20))[0]
        if len(label_r > 1):
            label_r = label_r[0]

        if left_x2 > right_x2:
            right = records[label_l, 3]
            left = 2*m1 - right
        else:
            left = records[label_r, 1]
            right = 2*m1 - left

        g_range = numpy.linspace(left, right, 80)
        xi2 = numpy.array([self.G_bin(g, nu, g_hat, bins, ig_num) for g_hat in g_range])

        gg4 = numpy.sum(g_range ** 4)
        gg3 = numpy.sum(g_range ** 3)
        gg2 = numpy.sum(g_range ** 2)
        gg1 = numpy.sum(g_range)
        xigg2 = numpy.sum(xi2 * (g_range ** 2))
        xigg1 = numpy.sum(xi2 * g_range)
        xigg0 = numpy.sum(xi2)
        cov = numpy.linalg.inv(numpy.array([[gg4, gg3, gg2], [gg3, gg2, gg1], [gg2, gg1, len(g_range)]]))
        paras = numpy.dot(cov, numpy.array([xigg2, xigg1, xigg0]))

        if pic_path:
            plt.scatter(g_range, xi2)
            plt.plot(g_range, paras[0]*g_range**2+paras[1]*g_range+paras[2])
            s = str(round(paras[0],3)) + " " + str(round(paras[1],3)) + " " + str(round(paras[2],3))
            plt.title(s)
            plt.savefig(pic_path)
            plt.close()

        g_sig = numpy.sqrt(1 / 2. / paras[0])
        g_h = -paras[1] / 2 / paras[0]

        return g_h, g_sig

    def fmin_g_new(self, g, nu, bin_num, ig_num=0, pic_path=False, left=-0.1, right=0.1):
        # nu = N + U for g1
        # nu = N - U for g2
        # temp_data = numpy.sort(numpy.abs(g))#[:int(len(g)*0.99)]
        # bin_size = len(temp_data)/bin_num*2
        # bins = numpy.array([temp_data[int(i*bin_size)] for i in range(1, int(bin_num / 2))])
        # bins = numpy.sort(numpy.append(numpy.append(-bins, [0.]), bins))
        # bound = numpy.max(numpy.abs(g)) * 100.
        # bins = numpy.append(-bound, numpy.append(bins, bound))
        bins = self.set_bin(g,bin_num)
        iters = 0
        change = 1
        while change == 1:
            change = 0
            mc = (left + right) / 2.
            mcl = (mc + left) / 2.
            mcr = (mc + right) / 2.
            fmc = self.G_bin(g, nu, mc, bins, ig_num)
            fmcl = self.G_bin(g, nu, mcl, bins, ig_num)
            fmcr = self.G_bin(g, nu, mcr, bins, ig_num)
            temp = fmc + 20
            if fmcl > temp:
                left = (mc + mcl)/2.
                change = 1
            if fmcr > temp:
                right = (mc + mcr)/2.
                change = 1
            iters += 1
            if iters > 16:
                break
        g_range = numpy.linspace(left, right, 200)
        xi2 = numpy.array([self.G_bin(g, nu, g_hat, bins, ig_num) for g_hat in g_range])

        gg4 = numpy.sum(g_range ** 4)
        gg3 = numpy.sum(g_range ** 3)
        gg2 = numpy.sum(g_range ** 2)
        gg1 = numpy.sum(g_range)
        xigg2 = numpy.sum(xi2 * (g_range ** 2))
        xigg1 = numpy.sum(xi2 * g_range)
        xigg0 = numpy.sum(xi2)
        cov = numpy.linalg.inv(numpy.array([[gg4, gg3, gg2], [gg3, gg2, gg1], [gg2, gg1, len(g_range)]]))
        paras = numpy.dot(cov, numpy.array([xigg2, xigg1, xigg0]))

        if pic_path:
            plt.scatter(g_range, xi2)
            plt.plot(g_range, paras[0] * g_range ** 2 + paras[1] * g_range + paras[2])
            s = str(round(paras[0], 3)) + " " + str(round(paras[1], 3)) + " " + str(round(paras[2], 3))
            plt.title(s)
            plt.savefig(pic_path)
            plt.close()

        g_sig = numpy.sqrt(1 / 2. / paras[0])
        g_h = -paras[1] / 2 / paras[0]
        return g_h, g_sig


