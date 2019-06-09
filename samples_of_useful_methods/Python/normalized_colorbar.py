import matplotlib.pyplot as plt
import numpy

n200 = ..
idx = ..

norm = plt.Normalize(vmin=numpy.min(n200[idx]), vmax=numpy.max(n200[idx]))
cmap = plt.get_cmap('YlOrRd')
cl = cmap(norm(n200[idx]))

ax.scatter(ra[idx], dec[idx], color=cl, s=3)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm._A = []

plt.colorbar(sm, ax=ax)


# colorbar of subfigure
img = Image_Plot()
img.create_subfig(1,2)
fig1 = img.axs[0][0].imshow(kappa)
plt.colorbar(fig1, ax=img.axs[0][0])
fig2 = img.axs[0][1].imshow(kappa)
plt.colorbar(fig2, ax=img.axs[0][1])
img.show_img()
