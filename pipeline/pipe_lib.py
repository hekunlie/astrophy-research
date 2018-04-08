import numpy
import numpy.ctypeslib as nclib
from scipy.optimize import least_squares


class Fresh:

    def __init__(self, image_x, image_y, size, bad_pixel, seed):
        # the size is the size of each stamp of galaxy subtracted from the exposure image.
        # image_x (y) is the size of each exposure image.
        # bad_pixel is the threshold less than which the pixel will be treated as bad pixel
        # rng is a initialed random number generator
        # clib is the object to call the functions in the library cpp_lib.so
        self.size = size
        self.img_x = image_x
        self.img_y = image_y
        self.rng = numpy.random.RandomState(seed)
        self.bad_pix = bad_pixel
        self.clib = nclib.load_library("/lib/cpp_lib", ".")

    def fit_templare(self, x, para):
        return para[0]*numpy.exp(-(x-para[1])/2/para[2]**2)

    def noise_fit(self, image, bin_num=200):
        # to fit the sigma of the white noise using the least-square method
        pixels = numpy.sort(image[image > self.bad_pix])
        nums, bins = numpy.histogram(pixels[:int(len(pixels)*0.9)], bin_num)

    def detector(self, exposure_list):
        pass

