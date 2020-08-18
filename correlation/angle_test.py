from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from plot_tool import Image_Plot
from astropy.coordinates import SkyCoord
import astropy.units as units
import numpy


def separation(RA1, DEC1, RA2, DEC2):
    dec1_rad = DEC1/180*numpy.pi
    dec2_rad = DEC2/180*numpy.pi

    diff_ra_rad = (RA1 - RA2)/180*numpy.pi

    cos_dec1 = numpy.cos(dec1_rad)
    cos_dec2 = numpy.cos(dec2_rad)
    sin_dec1 = numpy.sin(dec1_rad)
    sin_dec2 = numpy.sin(dec2_rad)

    sin_diff_ra = numpy.sin(diff_ra_rad)
    cos_diff_ra = numpy.cos(diff_ra_rad)

    m1 = cos_dec2 * sin_diff_ra
    m2 = cos_dec1 * sin_dec2 - sin_dec1 * cos_dec2*cos_diff_ra
    m = numpy.sqrt(m1*m1 + m2 * m2)
    n = sin_dec1 * sin_dec2 + cos_dec1 * cos_dec2*cos_diff_ra
    return numpy.abs(numpy.arctan2(m, n))/numpy.pi*180

# separation angle
# rng = numpy.random.RandomState(123)
# num = 20
# x1,y1 = rng.uniform(3, 350,num), rng.uniform(20,60,num)
# # x2,y2 = numpy.random.uniform(0, 360,num), numpy.random.uniform(-90,90,num)
# cos_y1 = numpy.cos(y1/180*numpy.pi)
#
#
# delta_theta_1 = numpy.zeros((num,num))
# delta_theta_2 = numpy.zeros((num,num))
#
#
# for i in range(num):
#     c1 = SkyCoord(x1[i]*units.deg, y1[i]*units.deg, frame='fk5')
#
#     dx = rng.uniform(-2, 2, num)
#     dy = rng.uniform(-2, 2, num)
#
#     for j in range(num):
#
#         x2, y2 = x1[i]+dx[j], y1[i]+dy[j]
#         cos_y2 = numpy.cos(y2/180*numpy.pi)
#         c2 = SkyCoord(x2*units.deg, y2*units.deg, frame='fk5')
#
#         sep_1 = numpy.sqrt((x2*cos_y1[i] - x1[i]* cos_y1[i]) ** 2 + (y2- y1[i]) ** 2)
#         # sep_1 = numpy.sqrt((x2* cos_y2 - x1[i]* cos_y1[i])**2 + (y2 - y1[i]) ** 2)
#         # sep_1 = numpy.sqrt((x2 - x1[i]) ** 2 + (y2 - y1[i]) ** 2)
#         sep_2 = c1.separation(c2).deg
#         sep_3 = separation(x1[i],y1[i],x2,y2)
#         # print(sep_1)
#         delta_theta_1[i,j] = sep_1
#         delta_theta_2[i,j] = sep_2
#         print("Point1: %.4f  %.4f"%(x1[i],y1[i]))
#         print("Point2: %.4f  %.4f"%(x2, y2))
#         print("Degree: %.4f  %.4f  %.4f  %.4f  %.4f"%(sep_1, sep_2, sep_3, sep_3-sep_2,sep_1-sep_2))
#         print("Arc minute: %.4f  %.4f  %.4f  %.4f  %.4f"%(sep_1*60, sep_2*60, sep_3*60, (sep_3-sep_2)*60,(sep_1-sep_2)*60))
#
# exit()
# img = Image_Plot()
# img.subplots(1,2)
# for i in range(num):
#     # img.axs[0][0].plot(x1,y1,c="C0")
#     img.axs[0][0].errorbar(y1[i], delta_theta_1[i].mean(), delta_theta_1[i].std()*2,c="C0")
#     img.axs[0][1].errorbar(y1[i], delta_theta_2[i].mean(), delta_theta_2[i].std()*2,c="C0")
# img.show_img()
# # cos_theta_1 = x2-x1


rng = numpy.random.RandomState(123)
num = 20
x1,y1 = rng.uniform(3, 350,num), rng.uniform(20,60,num)
# x1[0], y1[0] = 0, 10
# x2,y2 = numpy.random.uniform(0, 360,num), numpy.random.uniform(-90,90,num)
cos_y1 = numpy.cos(y1/180*numpy.pi)



for i in range(num):
    c1 = SkyCoord(x1[i]*units.deg, y1[i]*units.deg, frame='fk5')

    dx = rng.uniform(-2, 2, num)
    dy = rng.uniform(-2, 2, num)

    for j in range(num):

        x2, y2 = x1[i]+dx[j], y1[i]+dy[j]
        cos_y2 = numpy.cos(y2/180*numpy.pi)
        c2 = SkyCoord(x2*units.deg, y2*units.deg, frame='fk5')

        sep_1 = numpy.sqrt((x2*cos_y1[i] - x1[i]*cos_y1[i]) ** 2 + (y2- y1[i]) ** 2)
        cos_ang_1 = (y2 - y1[i])/sep_1
        sin_ang_1 = (x2*cos_y1[i] - x1[i]*cos_y1[i])/sep_1
        ang_1 = 90-numpy.arctan2(y2 - y1[i],(x2*cos_y1[i] - x1[i]*cos_y1[i]))/numpy.pi*180


        ang = c1.position_angle(c2).deg
        cos_ang = numpy.cos(ang*numpy.pi/180)
        sin_ang = numpy.sin(ang*numpy.pi/180)
        sep = c1.separation(c2).deg

        # print(sep_1)
        print("Point1: %.4f  %.4f"%(x1[i],y1[i]))
        print("Point2: %.4f  %.4f"%(x2, y2))
        print("Separation: %.4f  %.4f  %.4f"%(sep, sep_1, sep-sep_1))
        print("Separation: %.4f  %.4f  %.4f"%(sep*60, sep_1*60,sep*60-sep_1*60))
        print("Position angle: %.4f  sin: %.4f  cos: %.4f"%(ang, sin_ang, cos_ang))
        print("Position angle: %.4f  sin: %.4f  cos: %.4f"%(ang_1, sin_ang_1, cos_ang_1))
        print("diff Position angle: %.4f  sin: %.4f  cos: %.4f"%(ang_1 -ang, sin_ang_1 -sin_ang, cos_ang_1-cos_ang))

