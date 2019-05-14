import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy
from matplotlib.colors import ListedColormap

class Image_Plot:

    def __init__(self, fig_x=8, fig_y=6, fontsize=20, xy_lb_size=22, xy_tick_size=17,
                 legend_size=17, axis_linewidth=2, plt_line_width=2, cap_size=5, tick_len=6, pad_size=7):
        self.fig_x = fig_x
        self.fig_y = fig_y
        self.fontsize = fontsize
        self.xy_lb_size = xy_lb_size
        self.xy_tick_size = xy_tick_size
        self.legend_size = legend_size
        self.axis_linewidth = axis_linewidth
        self.plt_line_width = plt_line_width
        self.cap_size = cap_size
        self.tick_len = tick_len
        self.pad_size = pad_size
        self.figure = None
        self.axs = None

    def set_style_default(self):
        matplotlib.style.use('default')

    def set_style(self, font_style="serif"):
        self.set_style_default()
        plt.rcParams['font.family'] = font_style

    def create_subfig(self, ny, nx):
        # fig, sub_fig = plt.subplots(ny, nx, figsize=(int(nx*self.fig_x), int(ny*self.fig_y)))
        fig = plt.figure(figsize=(int(nx * self.fig_x), int(ny * self.fig_y)))
        sub_fig = [[] for i in range(ny)]
        for i in range(ny):
            for j in range(nx):
                ax = fig.add_subplot(ny, nx, i*nx+j+1)
                sub_fig[i].append(ax)
                sub_fig[i][j].tick_params(direction='in', labelsize=self.xy_tick_size, top=True, right=True, pad=self.pad_size)

                for axis in ["bottom", "left", "top", "right"]:
                    # the line width of the frame
                    sub_fig[i][j].spines[axis].set_linewidth(self.axis_linewidth)
                sub_fig[i][j].xaxis.set_tick_params(which="both", direction="in", length=self.tick_len, width=self.axis_linewidth)
                sub_fig[i][j].yaxis.set_tick_params(which="major", direction="in", length=self.tick_len, width=self.axis_linewidth)
                sub_fig[i][j].yaxis.set_tick_params(which="minor", direction="in", length=int(self.tick_len*0.6), width=self.axis_linewidth)
        self.figure = fig
        self.axs = sub_fig

    def axis_major_formatter(self, iy, ix, axis_nm, fmt):
        ticks_form = mtick.FormatStrFormatter(fmt)
        if axis_nm == 1:
            self.axs[iy][ix].xaxis.set_major_formatter(ticks_form)
        else:
            self.axs[iy][ix].yaxis.set_major_formatter(ticks_form)

    def tick_label(self, iy, ix, axis_nm, label, size=None):
        if not size:
            size = self.xy_lb_size
        if axis_nm == 1:
            self.axs[iy][ix].set_xlabel(label, fontsize=size)
        else:
            self.axs[iy][ix].set_ylabel(label, fontsize=size)

    def save_img(self, pic_path, tight=True):
        if tight:
            plt.savefig(pic_path, bbox_inches='tight')
        else:
            plt.savefig(pic_path)

    def show_img(self):
        plt.show()

    def subimg_adjust(self, h, w):
        plt.subplots_adjust(hspace=h, wspace=w)

    def get_colormap(self):
        N = 256
        vals = numpy.ones((N, 4))
        vals[:, 0] = 1
        vals[:, 1] = 0#numpy.linspace(0, 0.01, N)
        vals[:, 2] = 0#numpy.linspace(0, 0.04, N)
        vals[:, 3] = numpy.linspace(0, 1, N)
        return ListedColormap(vals)

    def imgshow(self,iy, ix, arr, colorbar=False):
        ax_ = self.axs[iy][ix].imshow(arr, interpolation='nearest', cmap=self.get_colormap())
        if colorbar:
            norm = matplotlib.colors.Normalize(vmin=arr.min(), vmax=arr.max())
            plt.colorbar(ax_, cmap=self.get_colormap(), norm=norm)
