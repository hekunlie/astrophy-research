import matplotlib
matplotlib.use("Agg")
import numpy
import h5py
from sys import path
path.append("/home/hklee/work/mylib/")
import tool_box
from plot_tool import Image_Plot
from astropy.cosmology import FlatLambdaCDM
from astropy import units
from astropy.coordinates import SkyCoord


area_id = 3

coeff_integ = 554.6821355281792
coeff_dist = 1662895.2081868195

h = 0.7
Om0 = 0.30
C_0_hat = 2.99792458
cosmos = FlatLambdaCDM(100*h,Om0)

# foreground
data_path = "/mnt/perc/hklee/CFHT/gg_lensing/data/foreground/PSZ2LenS/"
h5f_src = h5py.File(data_path+"w_3_3.hdf5","r")
keys = list(h5f_src.keys())

ra_f = h5f_src["RA"].value
dec_f = h5f_src["DEC"].value

co_dist_f = h5f_src["COM_DISTANCE"].value[0,0]
co_dist_integ_f = h5f_src["COM_DISTANCE_INTEG"].value[0,0]

phy_dist_f = h5f_src["PHY_DISTANCE"].value[0,0]
phy_dist_integ_f = h5f_src["PHY_DISTANCE_INTEG"].value[0,0]

z_f = h5f_src["Z"].value
h5f_src.close()

print("Position:",ra_f, dec_f,co_dist_f,phy_dist_f,z_f)

print("Max radius:",numpy.pi/180*co_dist_f*1.5)
print("Max radius:",numpy.pi/180*phy_dist_f*1.5)

fore_coord = SkyCoord(ra=ra_f*units.deg,dec=dec_f*units.deg,frame="fk5")


# CFHT source
h5f = h5py.File("/mnt/perc/hklee/CFHT/gg_lensing/data/cfht_cata/cfht_cata_cut.hdf5","r")

ra_s = h5f["/w_%d/RA"%area_id].value
dec_s = h5f["/w_%d/DEC"%area_id].value

co_dist_s = h5f["/w_%d/COM_DISTANCE"%area_id].value
co_dist_integ_s = h5f["/w_%d/COM_DISTANCE_INTEG"%area_id].value

phy_dist_s = h5f["/w_%d/PHY_DISTANCE"%area_id].value
phy_dist_integ_s = h5f["/w_%d/PHY_DISTANCE_INTEG"%area_id].value

z_s = h5f["/w_%d/Z"%area_id].value
zmin_s = h5f["/w_%d/Z_MIN"%area_id].value
zmax_s = h5f["/w_%d/Z_MAX"%area_id].value
odds = h5f["/w_%d/ODDS"%area_id].value
mag = h5f["/w_%d/MAG"%area_id].value

fitclass = h5f["/w_%d/FITCLASS"%area_id].value
mask_s = h5f["/w_%d/MASK"%area_id].value

m_bias = h5f["/w_%d/M"%area_id].value
c_bias = h5f["/w_%d/C"%area_id].value
e1 = h5f["/w_%d/E1"%area_id].value
e2 = h5f["/w_%d/E2"%area_id].value
weight = h5f["/w_%d/WEIGHT"%area_id].value
h5f.close()

z_sig95 = ((zmax_s - zmin_s)/2)

epsilon = 0.0000001

idx_mask = mask_s <= 1 + epsilon
idx_fit = numpy.abs(fitclass) < epsilon
idx_mag = mag <= 24.7
idx_w = weight > 0 - epsilon

idx_measure = idx_mask & idx_fit & idx_mag & idx_w

idx_odd = odds > 0.8
# idx_z = z_s > z_f + 0.05
idx_95 = zmin_s > z_f + 0.05
idx_z1 = z_s >= z_f
idx_z2 = z_s <= 1.2

idx_select = idx_odd & idx_95 & idx_z1 & idx_z2 & idx_measure


z_s = z_s[idx_select]

ra_s = ra_s[idx_select]
dec_s = dec_s[idx_select]
e1,e2 = e1[idx_select], e2[idx_select]
weight = weight[idx_select]
m_bias, c_bias = m_bias[idx_select], c_bias[idx_select]

co_dist_s, co_dist_integ_s = co_dist_s[idx_select], co_dist_integ_s[idx_select]
phy_dist_s, phy_dist_integ_s = phy_dist_s[idx_select], phy_dist_integ_s[idx_select]

source_num = z_s.shape[0]


print("Raw source:",source_num)

# separation
source_coord = SkyCoord(ra=ra_s * units.deg, dec=dec_s * units.deg, frame="fk5")
sep_ang = source_coord.separation(fore_coord).radian
sep_dist = sep_ang * phy_dist_f
# position angle
pos_ang = source_coord.position_angle(fore_coord).radian

cos_2theta = numpy.cos(2 * pos_ang)
sin_2theta = numpy.sin(2 * pos_ang)

et = e1 * cos_2theta - e2 * sin_2theta
ex = e1 * sin_2theta + e2 * cos_2theta

# crit = co_dist_s/co_dist_f/(co_dist_s-co_dist_f)/(1+z_f)
# crit_integ = co_dist_integ_s/co_dist_integ_f/(co_dist_integ_s-co_dist_integ_f)/(1+z_f)

crit = phy_dist_s / phy_dist_f / (phy_dist_s - phy_dist_f) * (1 + z_f)
crit_integ = phy_dist_integ_s / phy_dist_integ_f / (phy_dist_integ_s - phy_dist_integ_f) * (1 + z_f)

# radius bin
radius_num = 24
radius_bin = tool_box.set_bin_log(0.1, 25.12, radius_num + 1)
print("Radius bin:", radius_bin)

theta = numpy.linspace(0, 2 * numpy.pi, 1000)

# mask
mask_py = numpy.zeros((source_num, radius_num), dtype=numpy.intc) - 1

img = Image_Plot()
img.subplots(1, 1)

result = numpy.zeros((3, radius_num))

for i in range(radius_num):
    idx1 = sep_dist >= radius_bin[i]
    idx2 = sep_dist < radius_bin[i + 1]

    idx = idx1 & idx2
    pair_num = idx.sum()

    if pair_num > 0:

        mask_py[:, i][idx] = i

        wei = weight[idx] / (crit_integ[idx] ** 2)
        wei_sum = wei.sum()

        corr_m = 1 + numpy.sum(wei * m_bias[idx]) / wei_sum

        sub = crit_integ[idx] * coeff_integ * wei

        delta_sigma_t = numpy.sum(et[idx] * sub) / wei_sum / corr_m

        sigma = numpy.sum(et[idx].std() * sub) / wei_sum / corr_m / numpy.sqrt(pair_num)

        r_mean = sep_dist[idx].mean()

        result[0, i] = r_mean
        result[1, i] = delta_sigma_t
        result[2, i] = sigma

    print("%d in [%.4f, %.4f], %f, %f, %f" % (
    pair_num, radius_bin[i], radius_bin[i + 1], result[0, i], result[1, i], result[2, i]))
#     if i < 10:
#         print("%d in %d [%.4f, %.4f]"%
#               (pair_num,i,radius_bin[i],radius_bin[i+1]),ra_s[idx],dec_s[idx],sep_ang[idx],sep_dist[idx])


astro_data = numpy.loadtxt("/home/hklee/work/CFHT/gg_lensing/dens_cluster/dens_cluster.txt")
img.axs[0][0].errorbar(astro_data[0], astro_data[1], astro_data[2:4], c="k", marker="s", capsize=4, mfc="none", fmt=" ",
                       label="Dens",alpha=0.6)
img.axs[0][0].errorbar(result[0], result[1], c="C1", marker="o", capsize=4,  fmt=" ", label="T")
img.axs[0][0].legend()
img.axs[0][0].set_ylim(1, 5000)
img.axs[0][0].set_xlim(0.09, 30)

img.axs[0][0].set_xscale("log")
img.axs[0][0].set_yscale("log")
img.save_img("dens.png")