import matplotlib
matplotlib.use("Agg")
from sys import path,argv
path.append("/home/hklee/work/mylib")
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import numpy
import h5py


parent_path = argv[1]
h5f = h5py.File(parent_path+"/result.hdf5","r")
mc = h5f["/mc"].value
data = h5f["/data"].value
h5f.close()
print(mc)
h5f = h5py.File(parent_path +"/shear.hdf5","r")
g1 = h5f["/g1"].value
g2 = h5f["/g2"].value
h5f.close()

dg1 = g1 - data[:,0]
dg2 = g2 - data[:,2]
dg1_n = g1 - data[:,4]
dg2_n = g2 - data[:,6]

plot_data = [[data[:,0],data[:,1],data[:,2],data[:,3]],
            [data[:,4],data[:,5],data[:,6],data[:,7]]]

dgs = [[dg1*1000,dg2*1000],[dg1_n*1000,dg2_n*1000]]
pts_s = 40
img = Image_Plot(cap_size=1.5,fig_x=6,fig_y=4,xpad=0.4,ypad=0.2)
img.subplots(1,2)

x = numpy.linspace(-0.1,0.1,10)
for i in range(0,2):

    img.axs[0][i].errorbar(g1, plot_data[i][0], plot_data[i][1], c="C1", fmt=" ", capsize=img.cap_size)
    img.axs[0][i].errorbar(g2, plot_data[i][2], plot_data[i][3], c="C4", fmt=" ", capsize=img.cap_size)
    img.axs[0][i].scatter(g1, plot_data[i][0], c="C1",label="g1", s=pts_s, marker="s")
    img.axs[0][i].scatter(g2, plot_data[i][2], c="C4", label="g2", s=pts_s, marker="s")

    img.axs[0][i].text(0.05, 0.95, "$10^2 m_1 = %.3f \pm %.3f$"%(100*(mc[i,2]-1),100*mc[i,3]), color='green', ha='left',
            va='center', transform=img.axs[0][i].transAxes, fontsize=18)
    img.axs[0][i].text(0.05, 0.85, "$10^4 c_1 = %.3f \pm %.3f$"%(10000*mc[i,0],10000*mc[i,1]), color='green', ha='left',
            va='center', transform=img.axs[0][i].transAxes, fontsize=18)
    img.axs[0][i].text(0.05, 0.75, "$10^2 m_2 = %.3f \pm %.3f$"%(100*(mc[i,6]-1),100*mc[i,7]), color='green', ha='left',
            va='center', transform=img.axs[0][i].transAxes, fontsize=18)
    img.axs[0][i].text(0.05, 0.65, "$10^4 c_2 = %.3f \pm %.3f$"%(10000*mc[i,4],10000*mc[i,5]), color='green', ha='left',
            va='center', transform=img.axs[0][i].transAxes, fontsize=18)

    img.set_label(0,i,1, "True g")
    img.set_label(0,i,0, "Est g")

    yu = max(dgs[i][0].max(),dgs[i][1].max())
    yd = min(dgs[i][0].min(),dgs[i][1].min())
    dy = (yu - yd)*0.3
    ax_share = img.share_axis(0,i,1)
    # img.axs[1][i].scatter(g1, dgs[i][0], label="$g_1 - \hat{g}_1$")
    # img.axs[1][i].scatter(g2, dgs[i][1], label="$g_2 - \hat{g}_2$")
    # img.axs[1][i].legend(fontsize=img.legend_size,ncol=2,loc="upper left")

    ax_share.scatter(g1, dgs[i][0], edgecolor="C1",facecolor="none",label="$\\Delta g_1 = g_1 - \hat{g}_1$",marker="s",s=pts_s)
    ax_share.scatter(g2, dgs[i][1], edgecolor="C4",facecolor="none",label="$\\Delta g_2 = g_2 - \hat{g}_2$",marker="s",s=pts_s)

    ax_share.set_ylabel("$10^3\\Delta g$",fontsize=img.xy_lb_size)
    ax_share.set_ylim(yd-dy,yu+dy)

    img.axs[0][i].plot(x,x,c="grey",ls="dashed")
    img.axs[0][i].set_xlim(-0.06,0.06)
    img.axs[0][i].set_ylim(-0.06,0.06)
    if i== 0:
        img.axs[0][i].legend(fontsize=img.legend_size, ncol=2, frameon=False, loc="upper left",
                             bbox_to_anchor=(0.7, 1.15))
        ax_share.legend(fontsize=img.legend_size, ncol=2, loc="upper left", frameon=False, bbox_to_anchor=(1.13, 1.16))
lb = parent_path.split("_")[-1]
img.save_img(parent_path + "result.png")