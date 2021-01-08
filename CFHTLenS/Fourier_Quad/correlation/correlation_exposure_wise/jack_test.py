import h5py
from sys import argv
import numpy

def jack_label_check(labels):
    n = int(len(labels)/2)
    print("Labels: ", labels)
    np_labels = numpy.zeros((n,3))
    pl_all = []
    pl = []
    for i in range(n):
        tag = int(2*i)
        lb = [labels[tag], labels[tag+1]]
        np_labels[i,:2] = labels[tag], labels[tag+1]
        pl_all.append(lb)
        if lb not in pl:
            pl.append(lb)
            np_labels[i,2] = 0
        else:
            np_labels[i, 2] = 99
            for j in range(len(pl_all)):
                a,b = pl_all[j]
                c,d = lb
                if (a == c and b == d) or (a == d and b == c):
                    np_labels[j, 2] = 99
            print("duplicate: ", lb)

    print("Labels: ",np_labels)

name = argv[1]

h5f = h5py.File("%s"%name,"r")
buffer_num = h5f["/buffer_num"][()][0]
print(buffer_num)
for buffer_tag in range(buffer_num):
    print("Buffer %d"%buffer_tag)
    jack_label = h5f["/%s/jack_label"%buffer_tag][()]
    jack_label_check(jack_label)
h5f.close()


