import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy
from matplotlib.colors import ListedColormap

class Image_Plot:

    def __init__(self, fig_x=5, fig_y=4, xpad=0.1, ypad=0.1, fontsize=18, xy_lb_size=18, xy_tick_size=18,
                 legend_size=18, axis_linewidth=2, plt_line_width=2, cap_size=4, tick_len=6, pad_size=7, pts_size=2):
        self.fig_x = fig_x
        self.fig_y = fig_y
        self.xpad = xpad
        self.ypad = ypad
        self.fontsize = fontsize
        self.xy_lb_size = xy_lb_size
        self.xy_tick_size = xy_tick_size
        self.legend_size = legend_size
        self.axis_linewidth = axis_linewidth
        self.plt_line_width = plt_line_width
        self.cap_size = cap_size
        self.tick_len = tick_len
        self.pad_size = pad_size
        self.pts_size = pts_size
        self.figure = None
        self.axs = None
        self.row = None
        self.col = None

    def set_style_default(self):
        matplotlib.style.use('default')

    def set_style(self, font_style="serif"):
        self.set_style_default()
        plt.rcParams['font.family'] = font_style

    def subplots(self, ny, nx):

        dx = 1 / ( nx * (1 + self.xpad) + self.xpad)
        dy = 1 / (ny * (1 + self.ypad) + self.ypad)

        # fig = plt.figure( figsize=(fx/dx, fy/dy))
        # axs = [ fig.add_axes( [ j * dx* (1+xpad)+ dx*xpad,\
        #                         (ny-1-i) * dy*(1+ypad)+dy*ypad, \
        #                       dx, dy ]) \
        #         for i in range(ny)\
        #             for j in range(nx) ]

        self.row = ny
        self.col = nx
        # fig, sub_fig = plt.subplots(ny, nx, figsize=(int(nx*self.fig_x), int(ny*self.fig_y)))

        fig = plt.figure(figsize=(self.fig_x/dx, self.fig_y/dy))
        sub_fig = [[] for i in range(ny)]
        for i in range(ny):
            for j in range(nx):
                ax = fig.add_axes([j * dx* (1+self.xpad)+ dx*self.xpad,(ny-1-i)*dy*(1+self.ypad)+dy*self.ypad, dx, dy])
                sub_fig[i].append(ax)
                sub_fig[i][j].tick_params(direction='in', labelsize=self.xy_tick_size, top=True, right=True, pad=self.pad_size)

        self.figure = fig
        self.axs = sub_fig

    def scatter_pts(self, iy,ix, x, y, fxy, pts_size=3, color_map='YlOrRd',sci_cb=True):
        norm = plt.Normalize(vmin=numpy.min(fxy), vmax=numpy.max(fxy))
        cmap = plt.get_cmap(color_map)
        cl = cmap(norm(fxy))

        self.axs[iy][ix].scatter(x, y, color=cl, s=pts_size)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        cb = self.figure.colorbar(sm, ax=self.axs[iy][ix])
        if sci_cb:
            cb.formatter.set_powerlimits((0, 0))
            cb.update_ticks()

    def axis_type(self, axis_nm, spine, tick_len=None, tick_width=None, direction="in"):

        if not tick_len:
            tick_len = self.tick_len
        if not tick_width:
            tick_width = self.axis_linewidth

        for axis in ["bottom", "left", "top", "right"]:
            # the line width of the frame
            for i in range(self.row):
                for j in range(self.col):
                    self.axs[i][j].spines[axis].set_linewidth(tick_width)
                    if axis_nm == 0:
                        self.axs[i][j].yaxis.set_tick_params(which=spine, direction=direction,
                                                             length=tick_len, width=tick_width)
                    else:
                        self.axs[i][j].xaxis.set_tick_params(which=spine, direction=direction,
                                                             length=tick_len, width=tick_width)

                    self.axs[i][j].tick_params(right=True, which='minor')
                    self.axs[i][j].tick_params(top=True, which='minor')

    def axis_sci_ticklabel(self, iy, ix, axis_nm, limits=(0,1)):
        if axis_nm == 0:
            self.axs[iy][ix].yaxis.get_major_formatter().set_powerlimits(limits)
        else:
            self.axs[iy][ix].xaxis.get_major_formatter().set_powerlimits(limits)

    def axis_major_formatter(self, iy, ix, axis_nm, fmt):
        ticks_form = mtick.FormatStrFormatter(fmt)
        if axis_nm == 1:
            self.axs[iy][ix].xaxis.set_major_formatter(ticks_form)
        else:
            self.axs[iy][ix].yaxis.set_major_formatter(ticks_form)

    def set_ticklabel(self, iy, ix, axis_nm, n, labels):
        loc = mtick.LinearLocator(n)
        fmt = mtick.FixedFormatter(labels)
        # y-axis
        if axis_nm == 0:
            self.axs[iy][ix].yaxis.set_major_locator(loc)
            self.axs[iy][ix].yaxis.set_major_formatter(fmt)
        # x-axis
        else:
            self.axs[iy][ix].xaxis.set_major_locator(loc)
            self.axs[iy][ix].xaxis.set_major_formatter(fmt)

    def set_ticklabel_str(self, iy, ix, axis_nm, loc, labels):
        # y-axis
        if axis_nm == 0:
            self.axs[iy][ix].set_yticks(loc)
            self.axs[iy][ix].set_yticklabels(labels)
        # x-axis
        else:
            self.axs[iy][ix].set_xticks(loc)
            self.axs[iy][ix].set_xticklabels(labels)

    def share_axis(self,iy,ix,axis_nm):
        if axis_nm == 0:
            share_ax = self.axs[iy][ix].twiny()
            share_ax.xaxis.set_tick_params(which="major", direction="in", length=self.tick_len,
                                            width=self.axis_linewidth)
            share_ax.xaxis.set_tick_params(which="minor", direction="in", length=int(self.tick_len * 0.6),
                                            width=self.axis_linewidth)
        else:
            share_ax = self.axs[iy][ix].twinx()
            share_ax.yaxis.set_tick_params(which="major", direction="in", length=self.tick_len,
                                            width=self.axis_linewidth)
            share_ax.yaxis.set_tick_params(which="minor", direction="in", length=int(self.tick_len * 0.6),
                                            width=self.axis_linewidth)
        share_ax.tick_params(direction='in', labelsize=self.xy_tick_size, top=True, right=True, pad=self.pad_size)
        return share_ax

    def axs_text(self, iy, ix, y, x, text, text_fontsize=None, text_color="green", ha="left",
                 va="center",rotation="horizontal", ax_trans=True):
        if not text_fontsize:
            text_fontsize = self.legend_size
        if ax_trans:
            self.axs[iy][ix].text(x, y, text, color=text_color, ha=ha,  va=va,
                              transform=self.axs[iy][ix].transAxes, rotation=rotation, fontsize=text_fontsize)
        else:
            self.axs[iy][ix].text(x, y, text, color=text_color, ha=ha, va=va, rotation=rotation, fontsize=text_fontsize)


    def del_axis(self, iy, ix, axis_nm, box=None):
        """not to show the axis"""
        # delete the axis, but the box of axis will be preserved
        # axis_nm must be a list, [0] or [0,1]...
        for nm in axis_nm:
            if nm == 0:
                self.axs[iy][ix].get_yaxis().set_visible(False)
            else:
                self.axs[iy][ix].get_xaxis().set_visible(False)
        # delete the side of the box of axis
        # box must be a list of digits, 0 ~ 3
        if box:
            labels = ["left", "bottom", "right", "top"]
            for nm in box:
                self.axs[iy][ix].spines[labels[nm]].set_visible(False)

    def del_ticks(self, iy, ix, axis_nm):
        for nm in axis_nm:
            if nm == 0:
                self.axs[iy][ix].set_yticks([])
            else:
                self.axs[iy][ix].set_xticks([])

    def del_ticklabel(self, iy, ix, axis_nm):
        for nm in axis_nm:
            if nm == 0:
                self.axs[iy][ix].set_yticklabels([])
            else:
                self.axs[iy][ix].set_xticklabels([])

    def set_label(self, iy, ix, axis_nm, label, font="serif", fontsize=None):
        if not fontsize:
            fontsize = self.xy_lb_size
        fontdict = {'family': font, 'size': fontsize}
        if axis_nm == 1:
            self.axs[iy][ix].set_xlabel(label, fontdict=fontdict)
        else:
            self.axs[iy][ix].set_ylabel(label, fontdict=fontdict)

    def set_sci_ax(self,iy, ix, axis_nm,lower=0, upper=0):

        if axis_nm == 0:
            self.axs[iy][ix].yaxis.major.formatter.set_powerlimits((lower,upper))
        else:
            self.axs[iy][ix].xaxis.major.formatter.set_powerlimits((lower,upper))

    def set_legend(self, iy, ix, font="serif", size=None, loc="best", bbox_to_anchor=None):
        if not size:
            size = self.legend_size
        leg_prop = {'family': font, 'size': size}
        if bbox_to_anchor:
            self.axs[iy][ix].legend(loc=loc, bbox_to_anchor=bbox_to_anchor,prop=leg_prop)
        else:
            self.axs[iy][ix].legend(loc=loc, prop=leg_prop)

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

    def close_img(self):
        plt.close()