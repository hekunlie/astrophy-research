from sys import path, argv
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import numpy
import h5py
from Fourier_Quad import Fourier_Quad
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()


psf_tag = argv[1]
gal_num = int(argv[2])
# max_radius = float(argv[3])
# step_size = float(argv[4])

shear_num = 17
g = numpy.linspace(-0.04,0.04,shear_num)
theta = numpy.random.uniform(0,numpy.pi*2,shear_num)
g1 = g*numpy.cos(2*theta)
g2 = g*numpy.sin(2*theta)

stamp_size = 128

fq = Fourier_Quad(stamp_size,432)

# parameters
mr_num = 11
ss_num = 11
max_radius = numpy.linspace(0.4, 10, mr_num)
step_size = numpy.linspace(0.2, 2, ss_num)

pixel_scale = 0.187

psf_e = 0.1
psf_scale = 4


itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes = 4*mr_num*ss_num*itemsize
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
mcs = numpy.ndarray(buffer=buf1, dtype='d', shape=(4, ss_num, mr_num))


# PSF
if psf_tag == "c":
    psf_img = fq.cre_psf(psf_scale)
    if rank == 0:
        print("Circle PSF")
else:
    psf_img = fq.cre_psf(psf_scale, ellip_theta=(psf_e, 0))
    if rank == 0:
        print("Ellipticle PSF")

psf_pow = fq.pow_spec(psf_img)

fq.get_radius_new(psf_pow, 2)

# show PSF
# img = Image_Plot()
# img.subplots(1,2)
# img.axs[0][0].set_title("PSF")
# fig = img.axs[0][0].imshow(psf_img)
# img.figure.colorbar(fig, ax=img.axs[0][0])
# img.axs[0][1].set_title("e-PSF")
# fig = img.axs[0][1].imshow(epsf_img)
# img.figure.colorbar(fig, ax=img.axs[0][1])
# img.show_img()

task_list = [[i,j] for i in range(ss_num) for j in range(mr_num)]
my_task_list = tool_box.alloc(task_list, numprocs)[rank]

# simulate galaxy images

for ij in my_task_list:
    i,j = ij
    # measured shear
    mgs = numpy.zeros((2, shear_num))

    if psf_tag == "c":
        for k in range(shear_num):

            pts = fq.ran_pts(num=30, radius=max_radius[j], step=step_size[i])

            GN = numpy.zeros((3,))

            for ik in range(gal_num):

                for ikk in range(4):
                    #                     print(i,j,k,ik,ikk)
                    pts_r = fq.rotate(pos=pts, theta=ikk * numpy.pi / 4)
                    pts_s = fq.shear(pts_r, g1[k], g2[k])

                    # circle PSF
                    gal_img = fq.convolve_psf(pts_s, psf_scale)

                    gal_pow = fq.pow_spec(gal_img)

                    mg1, mg2, mn = fq.shear_est(gal_pow, psf_pow, F=True)[:3]
                    GN[0] += mg1
                    GN[1] += mg2
                    GN[2] += mn

            mgs[0, k] = GN[0] / GN[2]
            mgs[1, k] = GN[1] / GN[2]
        if rank == 0:
            print(g1 - mgs[0])
            print(g2 - mgs[1])
    else:
        for k in range(shear_num):

            pts = fq.ran_pts(num=30, radius=max_radius[j], step=step_size[i])

            GN = numpy.zeros((3,))

            for ik in range(gal_num):

                for ikk in range(4):
                    #                     print(i,j,k,ik,ikk)
                    pts_r = fq.rotate(pos=pts, theta=ikk * numpy.pi / 4)
                    pts_s = fq.shear(pts_r, g1[k], g2[k])

                    gal_img = fq.convolve_psf(pts_s, psf_scale, ellip_theta=(psf_e, 0))

                    gal_pow = fq.pow_spec(gal_img)

                    mg1, mg2, mn = fq.shear_est(gal_pow, psf_pow, F=True)[:3]
                    GN[0] += mg1
                    GN[1] += mg2
                    GN[2] += mn

            mgs[0, k] = GN[0] / GN[2]
            mgs[1, k] = GN[1] / GN[2]
        if rank == 0:
            print(g1 - mgs[0])
            print(g2 - mgs[1])

    mc1 = tool_box.fit_1d(g1, mgs[0], 1, method="scipy")
    mc2 = tool_box.fit_1d(g2, mgs[1], 1, method="scipy")

    mcs[0, i, j] = mc1[1] -1 # m1
    mcs[1, i, j] = mc1[0] # c1
    mcs[2, i, j] = mc2[1]-1 # m2
    mcs[3, i, j] = mc2[0] # c2

comm.Barrier()

if rank == 0:

    numpy.savez("./result/pts/cache_%s.npz"%psf_tag, mcs)
    pic_nm = "./result/pts/mc_%s.png"%psf_tag
    pic_nm_pdf = "./result/pts/mc_%s.pdf"%psf_tag

    y = numpy.zeros((mr_num*ss_num,))
    x = numpy.zeros((mr_num*ss_num,))
    for i in range(ss_num):
        for j in range(mr_num):
            tag = int(i*mr_num + j)
            y[tag] = step_size[i]
            x[tag] = max_radius[j]

    titles = ["$m_1$","$c_1$","$m_2$", "$c_2$"]
    img = Image_Plot(xpad=0.3,ypad=0.3,fig_x=4,fig_y=3)
    img.subplots(2, 2)
    for i in range(2):
        for j in range(2):
            tag = int(i*2 + j)
            z = mcs[tag].flatten()
            img.scatter_pts(i,j,x,y,z,pts_size=30,color_map="bwr")

            img.set_label(i,j,0,"Step Size")
            img.set_label(i,j,1,"Max Radius")
            img.axs[i][j].set_title(titles[tag],fontsize=img.xy_lb_size)
    img.save_img(pic_nm)
    img.save_img(pic_nm_pdf)
    print(pic_nm)
    # img.show_img()
comm.Barrier()