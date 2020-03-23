from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
from plot_tool import Image_Plot
import h5py

data_path = "E:/data/new_pdf/component_separation/stack"

h5f = h5py.File(data_path + "/shear.hdf5", "r")
g1 = h5f["/g1"][()]
g2 = h5f["/g2"][()]
h5f.close()

h5f = h5py.File(data_path + "/new_pdf_noise_free.hdf5", "r")

ave_result = h5f["/average"][()].T
pdf_result = h5f["/PDF"][()].T

new_pdf_result1 = h5f["/new_PDF/g1"][()]
new_pdf_result2 = h5f["/new_PDF/g2"][()]

h5f.close()

shear_num = ave_result.shape[1]

img = Image_Plot(cap_size=4,xpad=0.2, ypad=0.2)
img.subplots(1, 2)

print(ave_result)

mc1 = numpy.array(tool_box.data_fit(g1, ave_result[0],ave_result[1]))
mc2 = numpy.array(tool_box.data_fit(g2, ave_result[2],ave_result[3]))
mc1[0] = mc1[0] - 1
mc2[0] = mc2[0] - 1
print(mc1)
print(mc2)
img.axs[0][0].errorbar(0,mc1[0],mc1[1], label="ave m1", marker="s", capsize=img.cap_size)
img.axs[0][0].errorbar(0,mc2[0],mc2[1], label="ave m2", marker="s", capsize=img.cap_size)
img.axs[0][1].errorbar(0,mc1[2],mc1[3], label="ave c1", marker="s", capsize=img.cap_size)
img.axs[0][1].errorbar(0,mc2[2],mc2[3], label="ave c2", marker="s", capsize=img.cap_size)

print(pdf_result)
mc1 = numpy.array(tool_box.data_fit(g1, pdf_result[0], pdf_result[1]))
mc2 = numpy.array(tool_box.data_fit(g2, pdf_result[2], pdf_result[3]))
mc1[0] = mc1[0] - 1
mc2[0] = mc2[0] - 1
print(mc1)
print(mc2)
img.axs[0][0].errorbar(0,mc1[0],mc1[1], label="PDF m1", marker="s", capsize=img.cap_size)
img.axs[0][0].errorbar(0,mc2[0],mc2[1], label="PDF m2", marker="s", capsize=img.cap_size)
img.axs[0][1].errorbar(0,mc1[2],mc1[3], label="PDF c1", marker="s", capsize=img.cap_size)
img.axs[0][1].errorbar(0,mc2[2],mc2[3], label="PDF c2", marker="s", capsize=img.cap_size)

print(new_pdf_result1.shape)
data_row = 6
mc1s = numpy.zeros((4,4))
mc2s = numpy.zeros((4,4))
for i in range(4):
    mc1 = numpy.array(tool_box.data_fit(g1, new_pdf_result1[:,i*data_row], new_pdf_result1[:,i*data_row+1]))
    mc2 = numpy.array(tool_box.data_fit(g2, new_pdf_result2[:,i*data_row], new_pdf_result2[:,i*data_row+1]))
    mc1[0] = mc1[0] - 1
    mc2[0] = mc2[0] - 1
    mc1s[i] = mc1
    mc2s[i] = mc2
    print(mc1)
    print(mc2)
img.axs[0][0].errorbar(range(4), mc1s[:,0], mc1s[:,1], label="PDF m1-iter", marker="s", capsize=img.cap_size)
img.axs[0][0].errorbar(range(4), mc2s[:,0], mc2s[:,1], label="PDF m2-iter", marker="s", capsize=img.cap_size)
img.axs[0][1].errorbar(range(4), mc1s[:,2], mc1s[:,3], label="PDF c1-iter", marker="s", capsize=img.cap_size)
img.axs[0][1].errorbar(range(4), mc2s[:,2], mc2s[:,3], label="PDF c2-iter", marker="s", capsize=img.cap_size)
img.axs[0][0].legend(ncol=3)
img.axs[0][1].legend(ncol=3)
img.save_img(data_path + "/noise_free.png")
img.show_img()