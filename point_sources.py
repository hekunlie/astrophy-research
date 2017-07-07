from Fourier_Quad import Fourier_Quad
import numpy
import pandas

og1 = numpy.linspace(-0.01,0.01,21)
og2 = numpy.linspace(-0.01,0.01,21)
numpy.random.shuffle(og2)

size = 60
psf_ori = Fourier_Quad().cre_psf(4,size)
p=0
data_col = ["KSB_e1","BJ_e1","RG_e1","FQ_G1","FG_N","fg1", "KSB_e2","BJ_e2","RG_e2","FQ_G2","FG_N","fg2","FQ_U","FQ_V"]
for i in range(len(og1)):
    g1 = og1[i]
    g2 = og2[i]

    for m in range(2000):
        points= Fourier_Quad().ran_pos(50,size)
        for k in range(4):
            angle = numpy.pi*k/4
            points_rotate = Fourier_Quad().rotate(points,angle)
            points_shear = Fourier_Quad().shear(points_rotate,g1,g2)
            final = Fourier_Quad().convolve_psf(points_shear,4,size)
            G1,G2,N,U = Fourier_Quad().shear_est(final,psf_ori,size,F=False)[0:4]
            ith_row = numpy.array([0, 0, 0, G1, N, g1, 0, 0, 0, G2, N, g2, U, 0])
            if p == 0:
                data_matrix = ith_row
            else:
                data_matrix = numpy.row_stack((data_matrix, ith_row))
            p = 1
    print(numpy.sum(data_matrix[:,3])/numpy.sum(data_matrix[:,4]),g1,numpy.sum(data_matrix[:,9])/numpy.sum(data_matrix[:,4]),g2)
df = pandas.DataFrame(data_matrix,columns=data_col )
df.to_excel('/home/hklee/result/points/data.xlsx')