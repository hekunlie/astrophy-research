img.axs[i][j].set_yticks([int(nx - k * nx / 5) for k in range(6)])
img.axs[i][j].set_yticklabels(["%.1f" % (k * nx / 5. * pixel_scale) for k in range(6)])
img.axs[i][j].set_xticks([int(k * ny / 5) for k in range(6)])
img.axs[i][j].set_xticklabels(["%.1f" % (k * ny / 5. * pixel_scale) for k in range(6)])