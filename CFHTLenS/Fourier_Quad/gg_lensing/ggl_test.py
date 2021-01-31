from sys import path, argv
path.append("/home/hkli/work/mylib")
path.append("/home/hklee/work/mylib")
import numpy
import h5py
from plot_tool import Image_Plot
import numpy
import tool_box
from astropy.coordinates import SkyCoord
from astropy import units
from Fourier_Quad import Fourier_Quad

result_cata_path = "/mnt/perc/hklee/CFHT/gg_lensing/cata"
h5f = h5py.File(result_cata_path + "/pdf_inform.hdf5", "r")
sep_bins = h5f["/separation_bin"][()]
mg_bins = h5f["/mg_gt_bin"][()]
pdf_guess = h5f["/gt_guess"][()]
h5f.close()

bin_num = mg_bins.shape[0] - 1
print(bin_num)
inverse = range(int(bin_num / 2 - 1), -1, -1)
bin_num2 = int(bin_num / 2)
fq = Fourier_Quad(12, 12312)

h5f = h5py.File(result_cata_path + "/background/stack_data.hdf5", "r")
back_data_ori = h5f["/data"][()]
h5f.close()
idx = back_data_ori[:, 5] < 100
back_data = back_data_ori[idx]

print(idx.sum())
mg1, mg2, mn, mu, mv = back_data[:, 0], back_data[:, 1], back_data[:, 2], back_data[:, 3], back_data[:, 4]

back_ra, back_dec = back_data[:, 5], back_data[:, 6]
back_z, back_z_err = back_data[:, 8], back_data[:, 9]
back_com_dist = back_data[:, 10]

back_skypos = SkyCoord(ra=back_ra * units.deg, dec=back_dec * units.deg, frame="fk5")

print(back_z)
print(back_com_dist)
print(mg_bins)


ir = int(argv[1])
files = [argv[i] for i in range(2,len(argv))]
print(ir, files)

Gts = []
NUts = []
for file_tag in files:
    fore_expo_path = result_cata_path + "/foreground/%s-0.hdf5"%file_tag
    h5f = h5py.File(fore_expo_path, "r")
    fore_data = h5f["/data"][()]
    h5f.close()

    count = 0

    fore_c_num = fore_data.shape[0]

    for ic in range(fore_c_num):
        ic_z = fore_data[ic, 3]
        ic_com_dist = fore_data[ic, 4]
        ic_skypos = SkyCoord(ra=fore_data[ic, 0] * units.deg, dec=fore_data[ic, 1] * units.deg, frame="fk5")

        sep_radian = ic_skypos.separation(back_skypos).radian
        sep_radius = sep_radian * ic_com_dist

        idxz1 = back_z > ic_z + 0.1
        idxz2 = back_z - back_z_err > ic_z
        idxr1 = sep_radius >= sep_bins[ir]
        idxr2 = sep_radius < sep_bins[ir + 1]
        idx = idxz1 & idxz2 & idxr1 & idxr2
        sub_num = idx.sum()
        if sub_num > 0:
            count += sub_num

            pos_angle = ic_skypos.position_angle(back_skypos[idx]).radian
            sub_data = back_data[idx]
            mgt = sub_data[:, 0] * numpy.cos(pos_angle * 2) - sub_data[:, 1] * numpy.sin(pos_angle * 2)
            mnut = sub_data[:, 2] + sub_data[:, 3] * numpy.cos(4 * pos_angle) - sub_data[:, 4] * numpy.sin(4 * pos_angle)
            Gts.extend(mgt.tolist())
            NUts.extend(mnut.tolist())
    print(file_tag, count)
mgt = numpy.array(Gts)
mnut = numpy.array(NUts)

print(mgt.shape)

gh, gh_sig = fq.find_shear(mgt, mnut, bin_num)[:2]
print(ir, gh, gh_sig)

mg_bins_spec = fq.set_bin(mgt, bin_num, 1.1)

chi_sq_spec = numpy.array([fq.get_chisq(mgt, mnut, g_hat, mg_bins_spec, bin_num2, inverse, 0) for g_hat in pdf_guess])
chi_sq = numpy.array([fq.get_chisq(mgt, mnut, g_hat, mg_bins, bin_num2, inverse, 0) for g_hat in pdf_guess])

img = Image_Plot()
img.subplots(1, 1)
img.axs[0][0].plot(pdf_guess, chi_sq_spec)
img.axs[0][0].plot(pdf_guess, chi_sq)
img.save_img("./raidus_%d.png"%ir)
# img.show_img()
# img.close_img()