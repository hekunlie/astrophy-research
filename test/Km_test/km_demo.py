from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import numpy
import os
import shutil


def kmeans(x, y, ncent, max_iter=200, damp=0.7):

    xmi, ymi = x.mean(), y.mean()
    xrange = x.max() - x.min()
    yrange = y.max() - y.min()

    cent_xy = numpy.zeros((2, ncent))

    # initialize the cents
    cent_xy[0] = numpy.random.normal(xmi,xrange/8,ncent)
    cent_xy[1] = numpy.random.normal(ymi,yrange/8,ncent)
    # center label of each point
    cent_label = find_cent(x, y, cent_xy)
    # show the initial classification
    plot_label(x, y, cent_xy, cent_label[0], "Initial",pic_path="./pic/0.png")

    for i in range(max_iter):
        loss_0 = cent_label[1].sum()

        # update the centers
        update_cent(x,y,cent_xy, cent_label[0],damp)
        cent_label = find_cent(x, y, cent_xy)

        # show them
        plot_label(x, y, cent_xy, cent_label[0],"%d'th loop"%i,pic_path="./pic/%d.png"%(i+1))
        loss_1 = cent_label[1].sum()

        diff_ratio = (loss_0 - loss_1)/loss_0
        print("Loss function relative change: ",diff_ratio)
        if diff_ratio <= 0.00005:
            break

    return cent_xy,cent_label


def find_cent(x, y, cent_xy):
    pts_num = x.shape[0]
    ncent = cent_xy.shape[1]

    # ncent rows x pts_num cols
    # if the i'th point belongs to the j'th cent,
    # then, cent_label[j,i] = 1, else =0.
    cent_label = numpy.zeros((ncent, pts_num),dtype=numpy.intc)
    loss_fun = numpy.zeros((ncent, pts_num))

    # find the nearest center to the points
    for i in range(pts_num):
        rsq = (x[i] - cent_xy[0])**2 + (y[i] - cent_xy[1])**2
        rsq_min = rsq.min()
        npw = numpy.where(rsq == rsq_min)[0]
        loss_fun[npw, i] = rsq_min
        cent_label[npw, i] = 1
    return cent_label, loss_fun


def update_cent(x, y, cent_xy, cent_label, damp):
    ncent = cent_label.shape[0]
    for i in range(ncent):
        idx = cent_label[i] == 1
        if idx.sum() == 0:
            # sometimes, there's no points belong to some centers
            # assign these centers to some where near the data center
            cent_xy[0, i] = x.mean() + (x.max() - x.min())*numpy.random.normal(0,0.05,1)
            cent_xy[1, i] = y.mean()+ (y.max() - y.min())*numpy.random.normal(0,0.05,1)
        else:
            # update the centers
            # new center = the mean position of the group
            # the "damp" is used for suppressing the oscillation which may happen sometimes
            cent_xy[0,i] = x[idx].mean()*(1 - damp) + cent_xy[0,i]*damp
            cent_xy[1,i] = y[idx].mean()*(1 - damp) + cent_xy[1,i]*damp


def plot_label(x, y, cent_xy, cent_label, title, show=None, pic_path=None):
    ncent = cent_label.shape[0]
    img = Image_Plot()
    img.subplots(1,1)
    for i in range(ncent):
        img.axs[0][0].scatter(cent_xy[0][i],cent_xy[1][i], color="C%d"%i,marker="*",s=80)
        idx = cent_label[i] == 1
        img.axs[0][0].scatter(x[idx],y[idx],edgecolors="C%d"%i, color="none",label="%d, %d"%(i, idx.sum()),s=5)
    img.axs[0][0].set_title(title)
    img.axs[0][0].legend()
    if pic_path:
        img.save_img(pic_path)
    if show:
        img.show_img()
    img.close_img()


try:
    shutil.rmtree("./pic")
except:
    pass
os.makedirs("./pic")

ncent = 15
xcents = numpy.random.uniform(-20,20,ncent)#[-20, 4, 25]
ycents = numpy.random.uniform(-20,20,ncent)#[18,3,18]

num = 500
x,y = [], []
xy = numpy.zeros((2,int(num*ncent)))
for i in range(ncent):
    r = numpy.abs(numpy.random.normal(0,8,num))
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
img.save_img("./pic/00.png")
img.show_img()

cent_xy,cent_label = kmeans(xy[0],xy[1],ncent)
