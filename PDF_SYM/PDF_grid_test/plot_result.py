import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from Fourier_Quad import Fourier_Quad
import h5py
from plot_tool import Image_Plot
import tool_box


img = Image_Plot(xpad=0.28, ypad=0.2)
img.subplots(2,2)

data_path = "D:/TEMP/PDF_grid_test"
h5f = h5py.File(data_path + "/shear.hdf5","r")
g1t = h5f["/g1"][()]
g2t = h5f["/g2"][()]
h5f.close()

result = numpy.load(data_path + "/cache.npz")["arr_0"]

print(result.shape)

# noise free g1 ori PDF
mc = tool_box.data_fit(g1t,result[:,0],result[:,1])
label = "ori PDF_SYM noise-free\nm: %.4f(%.4f) c: %.5f(%.5f)"%(mc[0]-1, mc[1], mc[2], mc[3])
img.axs[0][0].errorbar(g1t, result[:,0],result[:,1],c="C0",fmt=" ", label=label,marker="s",ms=5)
# noise free g2 ori PDF
mc = tool_box.data_fit(g2t,result[:,2],result[:,3])
label = "ori PDF_SYM noise-free\nm: %.4f(%.4f) c: %.5f(%.5f)"%(mc[0]-1, mc[1], mc[2], mc[3])
img.axs[0][1].errorbar(g2t, result[:,2],result[:,3],c="C0",fmt=" ", label=label,marker="s",ms=5)

# noise free g1 grid PDF
mc = tool_box.data_fit(g1t,result[:,4],result[:,5])
label = "grid PDF_SYM noise-free\nm: %.4f(%.4f) c: %.5f(%.5f)"%(mc[0]-1, mc[1], mc[2], mc[3])
img.axs[0][0].errorbar(g1t, result[:,4],result[:,5],c="C1",fmt=" ", label=label,marker="s",ms=5)
# noise free g2 grid PDF
mc = tool_box.data_fit(g2t,result[:,6],result[:,7])
label = "grid PDF_SYM noise-free\nm: %.4f(%.4f) c: %.5f(%.5f)"%(mc[0]-1, mc[1], mc[2], mc[3])
img.axs[0][1].errorbar(g2t, result[:,6],result[:,7],c="C1",fmt=" ", label=label,marker="s",ms=5)

# noisy g1 ori PDF
mc = tool_box.data_fit(g1t,result[:,8],result[:,9])
label = "grid PDF_SYM noisy\nm: %.4f(%.4f) c: %.5f(%.5f)"%(mc[0]-1, mc[1], mc[2], mc[3])
img.axs[1][0].errorbar(g1t, result[:,8],result[:,9],c="C0",fmt=" ", label=label,marker="s",ms=5)
# noisy g2 ori PDF
mc = tool_box.data_fit(g2t,result[:,10],result[:,11])
label = "grid PDF_SYM noisy\nm: %.4f(%.4f) c: %.5f(%.5f)"%(mc[0]-1, mc[1], mc[2], mc[3])
img.axs[1][1].errorbar(g2t, result[:,10],result[:,11],c="C0",fmt=" ", label=label,marker="s",ms=5)

# noisy g1 ori PDF
mc = tool_box.data_fit(g1t,result[:,12],result[:,13])
label = "grid PDF_SYM noisy\nm: %.4f(%.4f) c: %.5f(%.5f)"%(mc[0]-1, mc[1], mc[2], mc[3])
img.axs[1][0].errorbar(g1t, result[:,12],result[:,13],c="C1",fmt=" ", label=label,marker="s",ms=5)
# noisy g2 ori PDF
mc = tool_box.data_fit(g2t,result[:,14],result[:,15])
label = "grid PDF_SYM noisy\nm: %.4f(%.4f) c: %.5f(%.5f)"%(mc[0]-1, mc[1], mc[2], mc[3])
img.axs[1][1].errorbar(g2t, result[:,14],result[:,15],c="C1",fmt=" ", label=label,marker="s",ms=5)


for i in range(2):
    img.axs[0][i].legend()
    img.axs[1][i].legend()

    img.axs[0][i].plot([-0.03,0.03],[-0.03,0.03],alpha=0.3,ls="--")
    img.axs[1][i].plot([-0.03,0.03],[-0.03,0.03],alpha=0.3,ls="--")
    for j in range(2):
        img.set_label(i, j, 1, "true g%d" % (i + 1))
        img.set_label(i,j, 0, "measured g%d" % (i + 1))
img.save_img(data_path + "/result.png")
img.show_img()