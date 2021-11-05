import numpy
from sklearn.cluster import KMeans
import healpy
from sklearn.cluster import KMeans


def get_healpy_pix(ra, dec, NSIDE):
    deg2rad = numpy.pi/180

    NPIX = healpy.nside2npix(NSIDE)
    resol = healpy.nside2resol(NSIDE) / deg2rad

    pixel_label = numpy.arange(NPIX)
    theta, phi = healpy.pix2ang(NSIDE, pixel_label)
    pix_dec = (-theta + numpy.pi / 2) / deg2rad
    pix_ra = phi / deg2rad

    ra_radian = ra * deg2rad
    # north pole is 0, pi is the south pole
    dec_radian = -dec * deg2rad + numpy.pi / 2

    pixel_indices = healpy.ang2pix(NSIDE, dec_radian, ra_radian)
    pixel_count = numpy.zeros(NPIX)

    for i in pixel_indices:
        pixel_count[i] += 1

    if pixel_count.sum() != ra.shape[0]:
        print("Number doesn't match!")
        exit()

    count_idx = pixel_count > 0

    # col: ra, dec
    pix_ra_dec = numpy.zeros((count_idx.sum(), 2), dtype=numpy.float32)
    pix_ra_dec[:, 0] = pix_ra[count_idx]
    pix_ra_dec[:, 1] = pix_dec[count_idx]

    src_pix_label = numpy.zeros_like(ra, dtype=numpy.intc) - 1

    for tag, i in enumerate(pixel_label[count_idx]):
        pix_ra_dec[tag, 0] = pix_ra[i]
        pix_ra_dec[tag, 1] = pix_dec[i]

        idx = pixel_indices == i
        src_pix_label[idx] = tag

    return pix_ra_dec, pixel_count[count_idx], src_pix_label, count_idx.sum()


def kmeans_pix(pixel_ra_dec, pixel_count, src_pixel_label, eff_pix_num, cluster_ncent, random_seed):

    weight = pixel_count/pixel_count.sum()
    group_label_pixel = KMeans(n_clusters=cluster_ncent, random_state=random_seed).fit_predict(pixel_ra_dec, sample_weight=weight)

    group_label = numpy.zeros_like(src_pixel_label, dtype=numpy.intc)
    pixel_label = numpy.arange(eff_pix_num)
    for i in range(cluster_ncent):
        idxi = group_label_pixel == i
        for j in pixel_label[idxi]:
            idxj = src_pixel_label == j
            group_label[idxj] = i

    return group_label, group_label_pixel