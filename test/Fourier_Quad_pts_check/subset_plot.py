from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
from plot_tool import Image_Plot
import h5py


data_path = "E:/data/new_pdf/component_separation/stack"

h5f = h5py.File(data_path + "/data_0_noisy.hdf5", "r")
data = h5f["/data"][()]
h5f.close()

print(data.shape)

mg1 = data[:,0]
mg2 = data[:,1]
mg = numpy.sqrt(mg1**2+mg2**2)
mn = data[:,2]
mu = data[:,3]
pk0 = data[:,4]/48/60
pk0_fit = data[:,5]/48/60
max_pk0_fit = data[:,6]/48/60
flux = data[:,7]

print(mg.shape)

img = Image_Plot(xpad=0.25)
img.subplots(1,3)

bin_num =300
num = 500000
sub_num = int(num/5)
hist_temp = numpy.zeros((num, ))
hist1_temp = numpy.zeros((num, ))
hist2_temp = numpy.zeros((num, ))
for i in range(5):
    # hist_temp[i*sub_num:(i+1)*sub_num] = mg[i*1400000:i*1400000+sub_num]
    # img.axs[0][0].hist(hist_temp[i * sub_num:(i + 1) * sub_num], bin_num, histtype="step", label="Component %d" % i)
    # hist1_temp[i*sub_num:(i+1)*sub_num] = mg1[i*1400000:i*1400000+sub_num]
    # hist2_temp[i*sub_num:(i+1)*sub_num] = mg2[i*1400000:i*1400000+sub_num]

    hist_temp[i*sub_num:(i+1)*sub_num] = pk0[i*1400000:i*1400000+sub_num]
    hist1_temp[i*sub_num:(i+1)*sub_num] = pk0_fit[i*1400000:i*1400000+sub_num]
    hist2_temp[i*sub_num:(i+1)*sub_num] = max_pk0_fit[i*1400000:i*1400000+sub_num]

    img.axs[0][0].hist(hist_temp[i*sub_num:(i+1)*sub_num],bin_num,histtype="step", label="Component %d"%i)
    img.axs[0][1].hist(hist1_temp[i*sub_num:(i+1)*sub_num], bin_num,histtype="step", label="Component %d"%i)
    img.axs[0][2].hist(hist2_temp[i*sub_num:(i+1)*sub_num], bin_num,histtype="step", label="Component %d"%i)

img.axs[0][0].hist(hist_temp, bin_num,histtype="step", label="Composite")
img.axs[0][1].hist(hist1_temp, bin_num,histtype="step", label="Composite")
img.axs[0][2].hist(hist2_temp, bin_num,histtype="step", label="Composite")
# fig = img.axs[0][1].hist2d(hist1_temp, hist2_temp,[100,100],cmin=1)[3]
# img.figure.colorbar(fig,ax=img.axs[0][1])

for i in range(3):
    img.axs[0][i].legend()
# img.save_img(data_path + "/pk_hist.png")
img.show_img()