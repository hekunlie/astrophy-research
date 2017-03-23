# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 14:29:31 2017

@author: 33471
"""
import os
import time
import numpy
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy import optimize
g1num = 9
g2num = 9
paths   = []
path  = "/home/lihekun/Downloads/result/elliptical/"
files = os.listdir(path)

for i in files:
    if ".dat" in i:
        paths.append(path+i)
fg1 = numpy.linspace(-0.008,0.008,g1num)
fg2 = numpy.linspace(-0.008,0.008,g2num)
g1={}
g2={}
fn1={} #for FQ method
fn2={}
#name: 0:KSB,1:BJ,2:RG,3:F_Q
for name in range(4): #dict g1/g2 = {'0':..{fg1/2[i]:[]},{fg1/2[i+1]:[]}..,'1':..,'2"...}
    for i in fg1:
        g1.setdefault(name,{})[i]=[]
        fn1[i]=[]
    for i in fg2:
        g2.setdefault(name,{})[i]=[]
        fn2[i]=[]

for k in paths: #put the data into the corresponding list    
    data = numpy.loadtxt(k,skiprows=1)
    for i in range(len(data)):
        a1 = numpy.abs(fg1-data[i,5])
        a1m= numpy.where(a1==numpy.min(a1))[0][0]
        a2 = numpy.abs(fg2-data[i,11])
        a2m= numpy.where(a2==numpy.min(a2))[0][0]
        for m in range(3):
            if data[i,m]!=-10.:
                g1[m][fg1[a1m]].append(data[i,m])
            if data[i,m+6]!=-10.:
                g2[m][fg2[a2m]].append(data[i,m+6])
        g1[3][fg1[a1m]].append(data[i,3])
        fn1[fg1[a1m]].append(data[i,4])
        g2[3][fg2[a2m]].append(data[i,9])
        fn2[fg2[a2m]].append(data[i,10])        

res_arr1 = numpy.zeros((12,g1num))# the first 4 rows are the ellipticity, 
                                 #the second 4 rows are the correspongding error bar,
                                 #the third 4 rows are the correspongding number of samples.
res_arr2 = numpy.zeros((12,g2num))
for i in range(4):
    for m in range(len(fg1)):
        if i!=3:
            num1 = len(g1[i][fg1[m]])            
            arr1 = numpy.array(g1[i][fg1[m]])
            ava1 = numpy.sum(arr1)/num1/1.6  #g1
            err1 = numpy.sqrt(numpy.sum((arr1/1.6-ava1)**2)/num1/(num1-1)) #error bar
            res_arr1[i,m] = ava1
            res_arr1[i+4,m] = err1
            res_arr1[i+8,m] = num1
        else:
            num1 = len(g1[i][fg1[m]])
            arr1 = numpy.array(g1[i][fg1[m]])
            ava1 = numpy.sum(arr1)/numpy.sum(numpy.array(fn1[fg1[m]])) #g1
            err1 = numpy.sqrt((numpy.sum((arr1/numpy.array(fn1[fg1[m]])-ava1)**2))/num1/(num1-1))
            res_arr1[i,m] = ava1
            res_arr1[i+4,m] = err1           
            res_arr1[i+8,m] = num1

    for m in range(len(fg2)):
        if i!=3:
            num2 = len(g2[i][fg2[m]])
            arr2 = numpy.array(g2[i][fg2[m]])
            ava2 = numpy.sum(arr2)/num2/1.6  #g2
            err2 = numpy.sqrt(numpy.sum((arr2/1.6-ava2)**2)/num2/(num2-1)) #error bar
            res_arr2[i,m] = ava2
            res_arr2[i+4,m] = err2
            res_arr2[i+8,m] = num2
        else:
            num2 = len(g2[i][fg2[m]])            
            arr2 = numpy.array(g2[i][fg2[m]])            
            ava2 = numpy.sum(arr2)/numpy.sum(numpy.array(fn2[fg2[m]])) #g2
            err2 = numpy.sqrt((numpy.sum((arr2/numpy.array(fn2[fg2[m]])-ava2)**2))/num2/(num2-1))
            res_arr2[i,m] = ava2
            res_arr2[i+4,m] = err2            
            res_arr2[i+8,m] = num2   
print "Classification complete"
# fitting the line
def fun(a,x,y):
    return a[0]*x+a[1]-y
a = numpy.array([1,0])
name = ['KSB','BJ','REGAUSS','Fourier_Quad']
for i in range(4):
    e1mc=optimize.least_squares(fun,a,args=(fg1,res_arr1[i,:])).x
    e2mc=optimize.least_squares(fun,a,args=(fg2,res_arr2[i,:])).x

    fig=plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111)
    xmajorLocator   = MultipleLocator(0.01) 
    xminorLocator   = MultipleLocator(0.001)
    ymajorLocator   = MultipleLocator(0.01) 
    yminorLocator   = MultipleLocator(0.001)
    ax.xaxis.set_major_locator(xmajorLocator)
#    plt.xaxis.set_major_formatter(xmajorFormatter)
    ax.yaxis.set_major_locator(ymajorLocator)
#    plt.yaxis.set_major_formatter(ymajorFormatter)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.errorbar(fg1,res_arr1[i,:],res_arr1[i+4,:],ecolor='black',elinewidth='1',fmt='none')
    ax.plot(fg1,e1mc[0]*fg1+e1mc[1],label=name[i],color='red')
    ax.plot(fg1,fg1,label='y=x',color='blue')
    ax.scatter(fg1,res_arr1[i,:],c='black')
#    for j in range(g1num):
#        ax.text(fg1[j],res_arr1[i,j],str(round(res_arr1[i+8,j]/1000,1))+"K",color = "red")
#    plt.xlim(-0.01,0.01)
#    plt.ylim(-0.01,0.01)
    plt.xlabel('True  g1',fontsize = 18)
    plt.ylabel('Estimated  g1',fontsize = 18)
    plt.title(name[i],fontsize = 18)
    plt.legend(fontsize = 14) 
    nm1 = name[i]+"_g1.png"
    plt.savefig(nm1)

    
    fig=plt.figure(figsize=(12,12))
    ax =fig.add_subplot(111)
    xmajorLocator   = MultipleLocator(0.01) 
    xminorLocator   = MultipleLocator(0.001)
    ymajorLocator   = MultipleLocator(0.01) 
    yminorLocator   = MultipleLocator(0.001)
    ax.xaxis.set_major_locator(xmajorLocator)
#    plt.xaxis.set_major_formatter(xmajorFormatter)
    ax.yaxis.set_major_locator(ymajorLocator)
#    plt.yaxis.set_major_formatter(ymajorFormatter)   
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)   
    ax.errorbar(fg2,res_arr2[i,:],res_arr2[i+4,:],ecolor='black',elinewidth='1',fmt='none')
    ax.plot(fg2,e2mc[0]*fg2+e2mc[1],label=name[i],color='red')
    ax.plot(fg2,fg2,label='y=x',color='blue')
    ax.scatter(fg2,res_arr2[i,:],c='black')  
#    for j in range(g2num):
#        ax.text(fg2[j],res_arr2[i,j],str(round(res_arr2[i+8,j]/1000,1))+"K",color = "red")
#    plt.xlim(-0.01,0.01)
#    plt.ylim(-0.01,0.01)
    plt.xlabel('True  g2',fontsize = 18)
    plt.ylabel('Estimated  g2',fontsize = 18)
    plt.title(name[i],fontsize = 18)
    plt.legend(fontsize = 14)    
    nm2 = name[i]+"_g2.png"
    plt.savefig(nm2)
print "Complete"    

# add something
                   
        
        
