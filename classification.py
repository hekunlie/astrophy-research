import  numpy
import os

def classify(g1num,g2num,g1s,g1e,g2s,g2e,path):
    fg1 = numpy.linspace(g1s, g1e, g1num)
    fg2 = numpy.linspace(g2s, g2e, g2num)
    paths = []
    files = os.listdir(path)
    print('starting...')
    for i in files:
        if ".txt" in i:
            paths.append(path + i)

    g1 = {}
    g2 = {}
    fn1 = {}  # for FQ method
    fn2 = {}
    # name: 0:KSB,1:BJ,2:RG,3:F_Q
    for name in range(4):  # dict g1/g2 = {'0':..{fg1/2[i]:[]},{fg1/2[i+1]:[]}..,'1':..,'2"...}
        for i in fg1:
            g1.setdefault(name, {})[i] = []
            fn1[i] = []
        for i in fg2:
            g2.setdefault(name, {})[i] = []
            fn2[i] = []

    for k in paths:  # put the data into the corresponding list
        data = numpy.loadtxt(k, skiprows=1)
        for i in range(len(data)):
            a1 = numpy.abs(fg1 - data[i, 5])
            a1m = numpy.where(a1 == numpy.min(a1))[0][0]
            a2 = numpy.abs(fg2 - data[i, 11])
            a2m = numpy.where(a2 == numpy.min(a2))[0][0]
            for m in range(3):
                if data[i, m] != -10.:
                    g1[m][fg1[a1m]].append(data[i, m])
                if data[i, m + 6] != -10.:
                    g2[m][fg2[a2m]].append(data[i, m + 6])
            g1[3][fg1[a1m]].append(data[i, 3])
            fn1[fg1[a1m]].append(data[i, 4])
            g2[3][fg2[a2m]].append(data[i, 9])
            fn2[fg2[a2m]].append(data[i, 10])

    res_arr1 = numpy.zeros((12, g1num))  # the first 4 rows are the ellipticity,
    res_arr2 = numpy.zeros((12, g2num))  # the second 4 rows are the correspongding error bar,
    # the third 4 rows are the correspongding number of samples.

    for i in range(4):
        for m in range(len(fg1)):
            if i != 3:
                num1 = len(g1[i][fg1[m]])
                arr1 = numpy.array(g1[i][fg1[m]])
                ava1 = numpy.sum(arr1) / num1 / 1.6  # g1
                err1 = numpy.std(arr1) / numpy.sqrt(num1)  # error bar
                res_arr1[i, m] = ava1
                res_arr1[i + 4, m] = err1
                res_arr1[i + 8, m] = num1
            else:
                num1 = len(g1[i][fg1[m]])
                arr1 = numpy.array(g1[i][fg1[m]])
                narr1 = numpy.array(fn1[fg1[m]])
                ava1 = numpy.sum(arr1) / numpy.sum(narr1)  # g1
                err1 = numpy.std(arr1) / numpy.mean(narr1) / numpy.sqrt(num1)
                res_arr1[i, m] = ava1
                res_arr1[i + 4, m] = err1
                res_arr1[i + 8, m] = num1

        for m in range(len(fg2)):
            if i != 3:
                num2 = len(g2[i][fg2[m]])
                arr2 = numpy.array(g2[i][fg2[m]])
                ava2 = numpy.sum(arr2) / num2 / 1.6  # g2
                err2 = numpy.std(arr2) / numpy.sqrt(num2)  # error bar
                res_arr2[i, m] = ava2
                res_arr2[i + 4, m] = err2
                res_arr2[i + 8, m] = num2
            else:
                num2 = len(g2[i][fg2[m]])
                arr2 = numpy.array(g2[i][fg2[m]])
                narr2 = numpy.array(fn2[fg2[m]])
                ava2 = numpy.sum(arr2) / numpy.sum(narr2)  # g2
                err2 = numpy.std(arr2) / numpy.mean(narr2) / numpy.sqrt(num2)
                res_arr2[i, m] = ava2
                res_arr2[i + 4, m] = err2
                res_arr2[i + 8, m] = num2
    print("Classification complete")

    with open(path + 'cache.dat', 'wb') as cache:
        numpy.savetxt(cache, numpy.column_stack((res_arr1, res_arr2)))
