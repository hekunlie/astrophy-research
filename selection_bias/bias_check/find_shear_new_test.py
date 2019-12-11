import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
import scipy
import h5py
from plot_tool import Image_Plot
import tool_box

def mc_f(x, c2, k, c1, b, a):
    return a * x ** 2 + (b + k) * x + c1 + c2
def new_fit(x,y):
    popt, pcov = scipy.optimize.curve_fit(mc_f, x, y)
    return popt[2:],popt[:2]


result_path = "D:/shear_result.hdf5"
h5f = h5py.File(result_path, "r")
mean_result_array = h5f["/mean_result"][()]
mc1 = h5f["/mc1"][()]
mc2 = h5f["/mc2"][()]
sym_result_array = h5f["/sym_result"][()]
chisq_g1 = h5f["/chisq_g1"][()]
chisq_g2 = h5f["/chisq_g2"][()]
print(mc1)
print(mc2)
h5f.close()

result_data = sym_result_array.copy()
g1 = mean_result_array[0]
g2 = mean_result_array[3]
shear_num = g1.shape[0]

print(g1.shape)

mg_l = numpy.zeros((4,shear_num))
mg_l[0] = result_data[1]
mg_l[1] = result_data[2]
mg_l[2] = result_data[4]
mg_l[3] = result_data[5]

mc1 = tool_box.data_fit(g1, mg_l[0], mg_l[1])
mc1_n = tool_box.data_fit_scipy(g1, mg_l[0], mg_l[1])
mc2 = tool_box.data_fit(g2, mg_l[2], mg_l[3])
mc2_n = tool_box.data_fit_scipy(g2, mg_l[2], mg_l[3])
print("From CPP")
print("%.5f (%.5f), %.5f (%.5f)"%(mc1[0]-1,mc1[1],mc1[2],mc1[3]))
# print("%.5f (%.5f), %.5f (%.5f)"%(mc1_n[0],mc1_n[1],mc1_n[2],mc1_n[3]))
print("%.5f (%.5f), %.5f (%.5f)"%(mc2[0]-1,mc2[1],mc2[2],mc2[3]))
# print("%.5f (%.5f), %.5f (%.5f)"%(mc2_n[0],mc2_n[1],mc2_n[2],mc2_n[3]))

check_result = numpy.zeros_like(mg_l)
new_result = numpy.zeros_like(mg_l)

for i in range(shear_num):
    fit_x1 = chisq_g1[i,20:]
    fit_y1 = chisq_g1[i,:20]
    fit_x2 = chisq_g2[i,20:]
    fit_y2 = chisq_g2[i,:20]

    coeff = tool_box.fit_1d(fit_x1, fit_y1, 2, "scipy")
    fit_curve_g1 = coeff[0] + coeff[1] * fit_x1 + coeff[2] * fit_x1 ** 2
    print(-coeff[1] / 2. / coeff[2])
    fit_slope_g1 = fit_y1 - fit_curve_g1
    coeff = tool_box.fit_1d(fit_x1, fit_y1-fit_slope_g1, 2, "scipy")
    print(-coeff[1] / 2. / coeff[2])
    fit_curve_g1 = coeff[0] + coeff[1] * fit_x1 + coeff[2] * fit_x1 ** 2
    coeff_new, slope_coeff = new_fit(fit_x1, fit_y1)
    print(coeff_new)#,slope_coeff)
    print(coeff)


    new_fit_curve_g1 = coeff_new[0] + coeff_new[1]*fit_x1 + coeff_new[2]*fit_x1**2



    check_result[0,i] = -coeff[1] / 2. / coeff[2]
    check_result[1,i] = 0.70710678118 / numpy.sqrt(coeff[2])
    new_result[0,i] = -coeff_new[1] / 2. / coeff_new[2]
    new_result[1,i] = 0.70710678118 / numpy.sqrt(coeff_new[2])


    coeff = tool_box.fit_1d(fit_x2, fit_y2, 2, "scipy")
    coeff_new, slope_coeff = new_fit(fit_x2, fit_y2)

    fit_curve_g2 = coeff[0] + coeff[1] * fit_x2 + coeff[2] * fit_x2 ** 2
    new_fit_curve_g2 = coeff_new[0] + coeff_new[1]*fit_x2 + coeff_new[2]*fit_x2**2
    fit_slope_g2 = fit_y2-fit_curve_g2

    check_result[2,i] = -coeff[1] / 2. / coeff[2]
    check_result[3,i] = 0.70710678118 / numpy.sqrt(coeff[2])
    new_result[2,i] = -coeff_new[1] / 2. / coeff_new[2]
    new_result[3,i] = 0.70710678118 / numpy.sqrt(coeff_new[2])

    img = Image_Plot()
    img.subplots(1,2)
    # img.axs[0][0].scatter(fit_x1, fit_y1)
    # img.axs[0][0].plot(fit_x1, fit_curve_g1)
    # img.axs[0][0].plot(fit_x1, new_fit_curve_g1)
    img.axs[0][0].plot(fit_x1, fit_slope_g1)
    img.axs[0][0].set_title("g1=%.4f"%g1[i])

    # img.axs[0][1].scatter(fit_x2, fit_y2)
    # img.axs[0][1].plot(fit_x2, fit_curve_g2)
    # img.axs[0][1].plot(fit_x2, new_fit_curve_g2)
    img.axs[0][1].plot(fit_x2, fit_slope_g2)
    img.axs[0][1].set_title("g2=%.4f"%g2[i])
    img.show_img()
    img.close_img()

print("Python check")
mc1 = tool_box.data_fit(g1, check_result[0], check_result[1])
mc2 = tool_box.data_fit(g2, check_result[2], check_result[3])
print("%.5f (%.5f), %.5f (%.5f)"%(mc1[0]-1,mc1[1],mc1[2],mc1[3]))
print("%.5f (%.5f), %.5f (%.5f)"%(mc2[0]-1,mc2[1],mc2[2],mc2[3]))

print("Python new")
mc1 = tool_box.data_fit(g1, new_result[0], new_result[1])
mc2 = tool_box.data_fit(g2, new_result[2], new_result[3])
print("%.5f (%.5f), %.5f (%.5f)"%(mc1[0]-1,mc1[1],mc1[2],mc1[3]))
print("%.5f (%.5f), %.5f (%.5f)"%(mc2[0]-1,mc2[1],mc2[2],mc2[3]))
