from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
from plot_tool import Image_Plot
import h5py

iters = 4

data_path = "E:/"

h5f = h5py.File(data_path + "/shear.hdf5", "r")
g1 = h5f["/g1"][()]
g2 = h5f["/g2"][()]
h5f.close()
shear_num = g1.shape[0]
print(shear_num)

file_name = "new_pdf_5.hdf5"
pic_nm1 = data_path + "/%s_1.png"%file_name
pic_nm2 = data_path + "/%s_2.png"%file_name
pic_nm3 = data_path + "/%s_3.png"%file_name

h5f = h5py.File(data_path + "/%s"%file_name, "r")
pdf_result = h5f["/PDF"][()].T
new_pdf_result1 = h5f["/new_PDF/g1"][()].T
new_pdf_result2 = h5f["/new_PDF/g2"][()].T
h5f.close()

print(new_pdf_result1.shape)

pdf_mc1 = numpy.array(tool_box.data_fit(g1, pdf_result[0], pdf_result[1]))
pdf_mc2 = numpy.array(tool_box.data_fit(g2, pdf_result[3], pdf_result[4]))
pdf_mc1[0] = pdf_mc1[0] - 1
pdf_mc2[0] = pdf_mc2[0] - 1

print(pdf_mc1)
print(pdf_mc2)

sub_row = int(new_pdf_result1.shape[0]/iters)
print(sub_row)

mc1s = numpy.zeros((iters, 4))
mc2s = numpy.zeros((iters, 4))

# # g\hat vs g_true
# img = Image_Plot(cap_size=4, xpad=0.2, ypad=0.2)
# img.subplots(1, 2)
#
for i in range(iters):
    mc1 = numpy.array(tool_box.data_fit(g1, new_pdf_result1[i * sub_row], new_pdf_result1[i * sub_row + 1]))
    mc2 = numpy.array(tool_box.data_fit(g2, new_pdf_result2[i * sub_row], new_pdf_result2[i * sub_row + 1]))
    mc1[0] = mc1[0] - 1
    mc2[0] = mc2[0] - 1
    mc1s[i] = mc1
    mc2s[i] = mc2
#     img.axs[0][0].errorbar(g1, new_pdf_result1[i * sub_row], new_pdf_result1[i * sub_row + 1], fmt=" ",
#                            capsize=img.cap_size, label="g1 iter-%d" % i, marker="s")
#     img.axs[0][1].errorbar(g2, new_pdf_result2[i * sub_row], new_pdf_result2[i * sub_row + 1], fmt=" ",
#                            capsize=img.cap_size, label="g2 iter-%d" % i, marker="s")
#     # img.axs[1][0].errorbar(g1, new_pdf_result1[i*sub_row], new_pdf_result1[i*sub_row + 1])
#     # img.axs[1][1].errorbar(g2, new_pdf_result2[i*sub_row], new_pdf_result2[i*sub_row + 1])
#
# img.axs[0][0].errorbar(g1, pdf_result[0], pdf_result[1], fmt=" ", capsize=img.cap_size, label="g1 ori_PDF", marker="v")
# img.axs[0][1].errorbar(g2, pdf_result[3], pdf_result[4], fmt=" ", capsize=img.cap_size, label="g2 ori_PDF", marker="v")
# img.axs[0][0].legend()
# img.axs[0][1].legend()
# img.show_img()



img = Image_Plot(cap_size=4, xpad=0.2, ypad=0.2)
img.subplots(1, 2)
shear_row_idx = [i*sub_row for i in range(iters)]
shear_sig_row_idx = [i*sub_row+1 for i in range(iters)]
print(shear_row_idx)
new_pdf_g1 = numpy.row_stack((pdf_result[0], new_pdf_result1[shear_row_idx]))
new_pdf_g1_sig = numpy.row_stack((pdf_result[1], new_pdf_result1[shear_sig_row_idx]))
new_pdf_g2 = numpy.row_stack((pdf_result[3], new_pdf_result2[shear_row_idx]))
new_pdf_g2_sig = numpy.row_stack((pdf_result[4], new_pdf_result2[shear_sig_row_idx]))
# print(new_pdf_g1)
# print(new_pdf_result1[shear_row_idx])
for i in range(shear_num):
    pass
    # img.axs[0][0].errorbar(range(iters+1), new_pdf_g1[:,i]-g1[i], new_pdf_g1_sig[:,i],
    #                        capsize=img.cap_size, label="$\delta$ g1-%d"%i, marker="s")
    # img.axs[0][1].errorbar(range(iters+1), new_pdf_g2[:,i]-g2[i], new_pdf_g1_sig[:,i],
    #                        capsize=img.cap_size, label="$\delta$g2-%d"%i, marker="s")

for i in range(iters+1):
    img.axs[0][0].errorbar(g1, new_pdf_g1[i]-g1, new_pdf_g1_sig[i],
                           capsize=img.cap_size, label="$\delta$ g1 iter-%d"%i, fmt=" ",marker="s",mfc="none")
    img.axs[0][1].errorbar(g2, new_pdf_g2[i]-g2, new_pdf_g2_sig[i],
                           capsize=img.cap_size, label="$\delta$g2 iter-%d"%i, fmt=" ",marker="s",mfc="none")

# img.axs[0][0].errorbar(g1, pdf_result[0], pdf_result[1], fmt=" ", capsize=img.cap_size, label="g1 ori_PDF", marker="v")
# img.axs[0][1].errorbar(g2, pdf_result[3], pdf_result[4], fmt=" ", capsize=img.cap_size, label="g2 ori_PDF", marker="v")
img.axs[0][0].legend(ncol=3)
img.axs[0][1].legend(ncol=3)
img.save_img(pic_nm1)
img.show_img()


# m&c vs iteration
img = Image_Plot(cap_size=4, xpad=0.2, ypad=0.2)
img.subplots(1, 2)


img.axs[0][0].errorbar(0,pdf_mc1[0],pdf_mc1[1], label="PDF m1", marker="s", capsize=img.cap_size)
img.axs[0][0].errorbar(0,pdf_mc2[0],pdf_mc2[1], label="PDF m2", marker="s", capsize=img.cap_size)
img.axs[0][1].errorbar(0,pdf_mc1[2],pdf_mc1[3], label="PDF c1", marker="s", capsize=img.cap_size)
img.axs[0][1].errorbar(0,pdf_mc2[2],pdf_mc2[3], label="PDF c2", marker="s", capsize=img.cap_size)

img.axs[0][0].errorbar(range(iters), mc1s[:,0], mc1s[:,1], label="PDF m1-iter", marker="s", capsize=img.cap_size)
img.axs[0][0].errorbar(range(iters), mc2s[:,0], mc2s[:,1], label="PDF m2-iter", marker="s", capsize=img.cap_size)
img.axs[0][1].errorbar(range(iters), mc1s[:,2], mc1s[:,3], label="PDF c1-iter", marker="s", capsize=img.cap_size)
img.axs[0][1].errorbar(range(iters), mc2s[:,2], mc2s[:,3], label="PDF c2-iter", marker="s", capsize=img.cap_size)
img.axs[0][0].legend(ncol=2)
img.axs[0][1].legend(ncol=2)
img.save_img(pic_nm2)
img.show_img()