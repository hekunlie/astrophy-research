import matplotlib
import matplotlib.pyplot as plt
import numpy

fig, ax = plt.subplots(figsize=(15,15))

cmap = plt.cm.get_cmap('YlOrRd')
xy = numpy.arange(0,10,1)
z = numpy.arange(0,10,1)

normalize = matplotlib.colors.Normalize(vmin=numpy.min(z), vmax=numpy.max(z))
# colors = [cmap(normalize(value)) for value in z]
colors = cmap(normalize(z))

ax.scatter(x=xy, y=xy,s=5,color=colors)

cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)

plt.show()