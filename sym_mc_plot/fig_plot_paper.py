import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import ScalarFormatter
from sys import path
path.append('D:/GitHub/astrophy-research/my_lib')
import tool_box


fmt='%2.f%%'
fig_x = 8
fig_y = fig_x*4/6
figs = (fig_x, fig_y)
fonts = 22
xy_lb_size = 20
legend_size = fonts - 4
axis_linewidth = 1.5

xticks = mtick.FormatStrFormatter(fmt)

mag = tool_box.mag_generator(5000, 19, 25)
radii = tool_box.radii_from_mags(mag, 0, 5)

fig = plt.figure(figsize=figs)
ax = fig.add_subplot(111)
ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
for axis in ["bottom", "left", "top", "right"]:
    ax.spines[axis].set_linewidth(axis_linewidth)
ax.xaxis.set_tick_params(which="both",direction="in",length=8, width=axis_linewidth)
ax.yaxis.set_tick_params(which="major",direction="in",length=8, width=axis_linewidth)
ax.yaxis.set_tick_params(which="minor",direction="in",length=4, width=axis_linewidth)


ax.scatter(mag,radii, s=3, color="black")
ax.set_ylim(0.05, 3)
x1, x2 = ax.set_xlim(18.9,25.1)
ax.plot([x1,x2],[0.187, 0.187], color="red",linestyle="--")

ax.set_yscale("log")
formattor = ScalarFormatter()
formattor.set_scientific(False)
ax.set_yticks([0.1, 1])
ax.set_yticklabels([r"0.1", r"1"])
ax.yaxis.set_major_formatter(formattor)
ax.set_xlabel("Magnitude",fontsize=xy_lb_size)
ax.set_ylabel("Scale length (arcsec)",fontsize=xy_lb_size)
plt.savefig("F:/works/figs/mag_radius.pdf",bbox_inches='tight')
plt.show()

