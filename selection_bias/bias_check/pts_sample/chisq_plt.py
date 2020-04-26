from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
import h5py
from plot_tool import Image_Plot


nm = ["cross_term", "noise_free", "noisy_cpp", "cross_term_noise_residual", "noise_residual"]

for j in range(5):
    total_path = "E:/data/new_pdf/galsim_sample/imgs_%d"%j
    data_path = total_path + "/chisq"

    h5f = h5py.File(data_path + "/shear.hdf5", "r")
    g1 = h5f["/g1"][()]
    g2 = h5f["/g2"][()]
    h5f.close()
    for n in nm:
        for i in range(10):
            h5f = h5py.File(data_path + "/chisq_%d_%s.hdf5"%(i,n),"r")
            chisq1 = h5f["/chisq1"][()]
            chisq2 = h5f["/chisq2"][()]
            shear_h = h5f["/shear"][()]
            h5f.close()
            print(chisq1.min(), chisq2.min())
            idx1 = chisq1 < 10000
            idx2 = chisq2 < 10000
            img = Image_Plot(xpad=0.15, ypad=0.15)
            img.subplots(1,1)
            img.axs[0][0].plot(shear_h[idx1], chisq1[idx1], label="g1")
            img.axs[0][0].plot(shear_h[idx2], chisq2[idx2], label="g2")
            img.axs[0][0].set_xlim(-0.11, 0.11)
            img.axs[0][0].legend()
            img.axs[0][0].set_title("g1=%.4f, g2=%.4f"%(g1[i], g2[i]))
            img.set_label(0,0,0,"$\chi^2$")
            img.set_label(0,0,1,"$g$")
            img.save_img(data_path + "/chisq_%s_%d.png"%(n,i))
            # img.show_img()
            img.close_img()


