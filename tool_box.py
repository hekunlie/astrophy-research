import Fourier_Quad
import pandas
from multiprocessing import Pool, Manager,freeze_support
import numpy
import time
from scipy import optimize

def task_distri(target_list, cpu_num):
    # it will divide the target_list into some piece (small lists in a diction) of which the number depends on the specific cpu core number
    # target_list must be a python list
    m, n = divmod(len(target_list), cpu_num)
    distri_pool = {}
    if m == 0 and n != 0:
        for i in range(n):
            distri_pool[i] = map(int, str(target_list[i]))
    elif m != 0:
        for i in range(cpu_num):
            distri_pool[i] = target_list[i * m:(i + 1) * m]
        if n != 0:
            for i in range(n):
                distri_pool[i].append(target_list[-i - 1])
    else:
        print("Caution! Something goes wrong!!!")
    return distri_pool

def list_add(target_list,files_paths):
    # this function is designed for 'function 'classify'
    # target_list is the target list that this function will put data array into
    # files_paths is a list of paths of excel files
    for i in range(len(files_paths)):
        data = pandas.read_excel(files_paths[i]).values
        if i==0:
            temp_data = data
        else:
            temp_data = numpy.row_stack((temp_data,data))
    target_list.append(temp_data)
    target_list.reverse()

def classify(files_path_list, cpu_num):
    # the data will be assembled into same big arrays of which the number equals to cpu_num
    paths_distri = task_distri(files_path_list,cpu_num)
    final_data_list = Manager().list()
    p = Pool()
    cl_ts = time.clock()
    for i in range(cpu_num):
        p.apply_async(list_add, args=(final_data_list,paths_distri[i],))
    p.close()
    p.join()
    cl_te = time.clock()
    return final_data_list,cl_te-cl_ts

def get_threshold(image, bin_num,upper_percent, noise_scale, bad_pix_num):
    # return the sigma of the Gaussian noise
    # bin_num is the number of the bins for estimating the distribution of the noise
    # upper_percent decides the number of chosen pixels
    # noise_scale is number of the corresponding bins interval around zero point
    # bad_pix_num is the number of bad pixels which will be ignored

    image = image.flatten()
    image = numpy.sort(image)
    noise = image[int(bad_pix_num+1):int(len(image)*upper_percent)]

    num, bins = numpy.histogram(noise, bin_num, normed=False)[0:2]
    # let the middle value of the bin interval  present the corresponding noise value
    size = bins[1]-bins[0]
    bins +=size/2.
    abs_bins = numpy.abs(bins)
    bin_min = numpy.min(abs_bins)
    mid = numpy.where(abs_bins==bin_min)[0][0]
    bins = bins[mid- noise_scale:mid+ noise_scale]
    pxs = num[mid - noise_scale:mid + noise_scale]
    pxs = pxs/numpy.sum(pxs)
    def fun(x,*para):
        return para[0]*numpy.exp(-(x-para[1])**2/2./para[2]**2)
    paras = optimize.curve_fit(fun,bins,pxs,p0=[50,20,100])[0]
    # plt.figure()
    # plt.plot(bins, pxs, color='g')
    # plt.plot(bins, paras[0]*numpy.exp(-(bins-paras[1])**2/2/paras[2]**2),color='r')
    # plt.show()
    return paras

def detect( mask, ini_y, ini_x,signal,maxi,max_p):
    # this function will add the source coordinates to the input signal list and return the cooridnates of  maximum value of the source
    # mask is the masked image of the original of which the values of pixels are either 1 (>signal threshold) or 0 (<signal threshold)
    # (ini_x,ini_y) is the initial coordinate of the signal
    # the maxi,max_p are the initial guesses of the maximum and the corresponding coordinates
    y_shape,x_shape = mask.shape
    if len(signal)>250 or ini_y<0 or ini_y>y_shape-1 or ini_x<0 or ini_x>x_shape-1:
        return None
    if mask[ini_y,ini_x]>0:
        signal.append((ini_y,ini_x))
        if mask[ini_y,ini_x]>maxi:
            max_p = (ini_y,ini_x)
            maxi = mask[ini_y,ini_x]
        mask[ini_y,ini_x]=0
        for cor in ((-1,0),(1,0),(0,-1),(0,1)):
            if ini_y+cor[0]<= y_shape-1 and ini_x+cor[1]<= x_shape-1 and mask[ini_y+cor[0],ini_x+cor[1]]>0:
                max_p = detect(mask,ini_y+cor[0],ini_x+cor[1],signal,maxi,max_p)
    return max_p

# if __name__=='__main__':
#     freeze_support()