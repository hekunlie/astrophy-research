import numpy
from sys import path
path.append("/home/hklee/work/mylib/")
path.append("D:/GitHub/astrophy-research/mylib/")
from plot_tool import Image_Plot
import scipy
import gauss_fit_fun



# scale = 1
gh = numpy.linspace(-0.2, 0.2, 100)
dg = gh[1] - gh[0]
rng = numpy.random.RandomState(1140)
total_num = 1000000
signal_ini = [0.04, -0.04]
mix_data = numpy.zeros((total_num,))
for i in range(len(signal_ini)):
    num_sub = int(total_num/len(signal_ini))
    mix_data[i*num_sub:(i+1)*num_sub] = rng.normal(signal_ini[i], 0.02, num_sub)

i = 1
st = mix_data.min()

ed_ = mix_data.max()
bins = [st]
while True:
    scale = st + i*dg
    if bins[-1] <= ed_:
        bins.append(scale)
    else:
        break
    i += 1

bins = numpy.array(bins)

# calculate the flow
num_move = gauss_fit_fun.get_flow(mix_data, 1, gh)
scale = num_move.sum()
# N/dg
x_fit, y_fit = gh[1:], num_move/scale/dg


img = Image_Plot(plt_line_width=2.5, fig_x=12, fig_y=9)

img.subplots(1, 1)
img.axs[0][0].scatter(x_fit, y_fit, c="black", s=15, label="data flow (normalized)")

# img.axs[0][1].hist(mix_data, bins, histtype="step", linewidth=img.plt_line_width,label="PDF")

# fitting
text_content = []

fit_res_1 = scipy.optimize.curve_fit(gauss_fit_fun.gauss_fun_2_, x_fit, y_fit)[0]
w1, w2 = fit_res_1

fx1 = gauss_fit_fun.gauss_fun(x_fit, 0.04, 0.02)*w1
fx2 = gauss_fit_fun.gauss_fun(x_fit, -0.04, 0.02)*w2

f_true = gauss_fit_fun.gauss_fun(x_fit, 0.04, 0.02)*0.5 + gauss_fit_fun.gauss_fun(x_fit, -0.04, 0.02)*0.5

text_content.append([w1, 0.04, 0.02])
text_content.append([w2, -0.04, 0.02])



ty = 0.8
gauss_fit_fun.img_text(img.axs[0][0], 0.05, ty, strs_content="weight,  $\mu$,  $\sigma$", size=20)
# print(text_content)
cs = ["C1","C2", "C3"]
for i in range(len(text_content)):
    gauss_fit_fun.img_text(img.axs[0][0], 0.05, ty-0.1-i*0.1, paras=text_content[i], size=20,c=cs[i])

img.axs[0][0].plot(x_fit, fx1, label="fx 1", linewidth=img.plt_line_width,c="C1",linestyle="-.")
img.axs[0][0].plot(x_fit, fx2, label="fx 2", linewidth=img.plt_line_width,c="C2",linestyle="-.")
img.axs[0][0].plot(x_fit, fx2+fx1, label="fx 1 + fx 2", linewidth=img.plt_line_width,c="C3")
img.axs[0][0].plot(x_fit, f_true, label="fx true", linewidth=img.plt_line_width, c="C4", linestyle="dotted")

img.set_label(0,0,0,"Normalized num", size=img.xy_lb_size)
img.set_label(0,0,1,"G", size=img.xy_lb_size)
img.axs[0][0].set_title("Fix $\sigma$ and $\mu$", fontsize=img.xy_lb_size)
img.axs[0][0].legend(fontsize=img.legend_size)
img.show_img()
