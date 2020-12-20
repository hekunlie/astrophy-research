from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import numpy


def find_cent(x, y, cent_xy):
    pts_num = x.shape[0]
    cent_label = numpy.zeros((pts_num, ),dtype=numpy.intc)
    print(cent_xy)
    for i in range(pts_num):
        rsq = (x[i] - cent_xy[0])**2 + (y[i] - cent_xy[1])**2
        npw = numpy.where(rsq == rsq.min())[0]
        print(npw,rsq)
        cent_label[i] = npw
    return cent_label


def kmeans(x, y ,ncent, max_iter=100, damp=0.4):
    pts_num = x.shape[0]

    xmi, ymi = x.mean(), y.mean()
    xrange = x.max() - x.min()
    yrange = y.max() - y.min()

    cent_xy = numpy.zeros((2, ncent))
    print(xmi,ymi)
    # initialize the cents
    cent_xy[0] = numpy.random.normal(xmi,xrange/8,ncent)
    cent_xy[1] = numpy.random.normal(ymi,yrange/8,ncent)
    cent_label = find_cent(x,y,cent_xy)

    img = Image_Plot()
    img.subplots(1,1)
    for i in range(ncent):
        idx = cent_label == i
        img.axs[0][0].scatter(x[idx],y[idx],label="%d"%i)
    img.axs[0][0].legend()
    img.show_img()

xcents = [-20,4,25]
ycents = [18,3,18]
ncent = len(xcents)
num = 200
x,y = [], []
xy = numpy.zeros((2,int(num*ncent)))
for i in range(ncent):
    r = numpy.abs(numpy.random.normal(0,2,num))
    theta = numpy.random.uniform(0,numpy.pi*2,num)

    xi = xcents[i] + r*numpy.cos(theta)
    yi = ycents[i] + r*numpy.sin(theta)
    x.append(xi)
    y.append(yi)
    xy[0,int(i*num):int((i+1)*num)] = xi
    xy[1,int(i*num):int((i+1)*num)] = yi


img = Image_Plot()
img.subplots(1,1)
for i in range(ncent):
    img.axs[0][0].scatter(x[i],y[i],label="%d"%i,s=5)
img.axs[0][0].legend()
img.show_img()

kmeans(xy[0],xy[1],3)