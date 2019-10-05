import numpy
from sys import path
path.append("/home/hklee/work/mylib/")
path.append("/home/hkli/work/mylib/")
path.append("/home/hklee/work/multi_shear_detect/")
path.append("D:/Github/astrophy-research/multi_shear_detect")
path.append("D:/Github/astrophy-research/mylib")
import tool_box
from plot_tool import Image_Plot
import gauss_fit_fun
import scipy


def fun_min(para, x,y):
    a, b, c = para[0], para[1], para[2]
    fx = gauss_fit_fun.gauss_coeff(x,a,b,c)
    return numpy.sum((fx - y)**2)

def fx_sub_true(x, w1, s1, w2, s2, sigma):
    f1 = w1*(gauss_fit_fun.df1(x, sigma)*s1 + gauss_fit_fun.df2(x, sigma)*s1**2 + gauss_fit_fun.df3(x,sigma)*s1**3+
                gauss_fit_fun.df4(x,sigma)*s1**4)
    f2 = w2*(gauss_fit_fun.df1(x, sigma)*s2 + gauss_fit_fun.df2(x, sigma)*s2**2 + gauss_fit_fun.df3(x,sigma)*s2**3 +
            gauss_fit_fun.df4(x,sigma)*s2**4)
    return f1+f2, f1,f2

def fx_sub_approx(n, x,w1,mu1,w2,mu2,sig):
    orders = [gauss_fit_fun.df1(x,sig), gauss_fit_fun.df2(x,sig), gauss_fit_fun.df3(x,sig),
              gauss_fit_fun.df4(x,sig), gauss_fit_fun.df5(x,sig), gauss_fit_fun.df6(x,sig)]
    f1 = 0
    f2 = 0
    for i in range(n):
        f1 += w1*(mu1**(i+1))*orders[i]
        f2 += w2*(mu2**(i+1))*orders[i]
    return f1+f2, f1, f2

def min_diff_1(x,a1,a2,a3,sig):
    # n = 2
    # sig = 0.3011
    orders = [gauss_fit_fun.df1(x,sig), gauss_fit_fun.df2(x,sig), gauss_fit_fun.df3(x,sig),
              gauss_fit_fun.df4(x,sig), gauss_fit_fun.df5(x,sig), gauss_fit_fun.df6(x,sig)]
    f = orders[0]*a1 + orders[1]*a2 + orders[2]*a3
    return f

def min_diff_2(x,w1, mu1, mu2):
    n = 3
    sig = 0.3011
    orders = [gauss_fit_fun.df1(x,sig), gauss_fit_fun.df2(x,sig), gauss_fit_fun.df3(x,sig),
              gauss_fit_fun.df4(x,sig), gauss_fit_fun.df5(x,sig), gauss_fit_fun.df6(x,sig)]
    f1 = 0
    f2 = 0
    for i in range(n):
        f1 += w1*(mu1**(i+1))*orders[i]
        f2 += (1-w1)*(mu2**(i+1))*orders[i]
    return f1+f2

w1, mu1, sig1 = 0.5, 0.01, 0.3
w2, mu2, sig2 = 0.5, -0.04, 0.3
print("     weight, signal, sigma")
print("True 1:",w1, mu1, sig1)
print("True 2:",w2, mu2, sig2)


# the cache is got from the servers
npz_f = numpy.load("E:/works/multi_shear/cache.npz")
x_fit, y_fit = npz_f["arr_0"], npz_f["arr_1"]

# fit single gauss
fit_res = scipy.optimize.curve_fit(gauss_fit_fun.gauss_coeff, x_fit, y_fit)[0]
# fx_fit = gauss_fit_fun.gauss_coeff(x_fit, fit_res[0],fit_res[1],fit_res[2])
fx_fit = gauss_fit_fun.gauss_coeff(x_fit, 1, 0, fit_res[2])
print("Single gauss fit: %.4f, %.4f, %.4f"%(fit_res[0],fit_res[1],fit_res[2]))

fx_ratio = y_fit/fx_fit - 1
fx_sub, fx_sub1, fx_sub2 = fx_sub_approx(3,x_fit,w1,mu1,w2,mu2,sig1)

# fit the ratio expression
idx1 = x_fit < 0.1
idx2 = x_fit > -0.1
idx = idx1 & idx2
# # one fitting method
# ratio_fit_result = scipy.optimize.curve_fit(min_diff_1,x_fit[idx],fx_ratio[idx],p0=[0.1,0.2,0.3,0.3])[0]
# a1,a2,a3,sig = ratio_fit_result
# print("%.7f  %.7f  %.7f  %.4f"%(a1,a2,a3,sig))
# fx_sub_fit = min_diff_1(x_fit, a1, a2, a3,sig)
# another method
ratio_fit_result = scipy.optimize.curve_fit(min_diff_2,x_fit[idx],fx_ratio[idx],p0=[0.01,0.01,0.01])[0]
a1,a2,a3 = ratio_fit_result
print("%.7f  %.7f  %.7f "%(a1,a2,a3))
fx_sub_fit = min_diff_2(x_fit, a1, a2, a3)


img = Image_Plot(plt_line_width=3)
img.subplots(1,1)


img.axs[0][0].scatter(x_fit,fx_ratio,c="k",s=10)#,linewidth=img.plt_line_width)
img.set_label(0,0,0,"$\\frac{P(\hat g)}{P(\hat g)_{fit}} - 1$")
img.set_label(0,0,1,"$\hat g$")
img.axs[0][0].plot(x_fit,fx_sub,label="True $\\frac{P(\hat g)}{P(\hat g)_{fit}} - 1$",linewidth=img.plt_line_width)

img.axs[0][0].plot(x_fit,fx_sub_fit,label="Fit $\\frac{P(\hat g)}{P(\hat g)_{fit}} - 1$",linewidth=img.plt_line_width)
img.axs[0][0].plot(x_fit[idx],fx_ratio[idx],label="fit range",linewidth=img.plt_line_width)

# img.axs[0][1].plot(x_fit,fx_sub1_fit,label="Fit component 1,g=%.2f"%mu1,linewidth=img.plt_line_width,linestyle="-.")
# img.axs[0][1].plot(x_fit,fx_sub2_fit,label="Fit component 2,g=%.2f"%mu2,linewidth=img.plt_line_width,linestyle="-.")

for i in range(1):
    ny, nx = divmod(i,2)
    img.axs[ny][nx].legend(fontsize=img.legend_size, ncol=1)
# img.subimg_adjust(h=0,w=0.25)
img.save_img("D:/fit.png")
img.show_img()

