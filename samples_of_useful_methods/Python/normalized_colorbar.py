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
