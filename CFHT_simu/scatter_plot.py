import matplotlib.pyplot as plt
import numpy
import h5py


h5f = h5py.File("F:/result.hdf5","r")
data = h5f["/data"].value
h5f.close()
print(data.shape)
print(data)
idx = data[:,19] > -1
# plt.scatter(data[:,0], data[:,8])
# plt.show()

flux2s = [data[:,0]]
snrs = []
fig = plt.figure(figsize=(18, 9))
for i in range(4, 9):
    flux2s.append(data[:,i])
# original SNR_F
for i in range(len(flux2s)):
    ax = fig.add_subplot(3,6,1+i)
    ax.hist(flux2s[i], 50)
    # ax.xaxis.get_major_formatter().set_powerlimits((0,1))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
# flux_ext
for i in range(5):
    ax = fig.add_subplot(3,6,8+i)
    ax.hist(data[:,9+i], 50)
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
# SNR_ext
for i in range(5):
    size = data[:,14+i][idx]
    flux = data[:,9+i][idx]
    snr = flux/numpy.sqrt(size)/60
    ax = fig.add_subplot(3,6,14+i)
    ax.hist(snr, 50)
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    snrs.append(snr)
plt.savefig("F:/mag21_r0.8_btr_0.4_e(0.6-0.7).pdf",bbox_inches='tight')
plt.show()
plt.close()

fig = plt.figure(figsize=(20,12))
for i in range(1,len(flux2s)):
    ax = fig.add_subplot(3,5,i)
    ax.scatter(flux2s[0][:5000], flux2s[i][:5000], s=0.1)
    # if i == 1:
    #     ys = ax.set_ylim()
    #     xs = ax.set_xlim()
    # else:
    #     ax.set_xlim(xs[0], xs[1])
    #     ax.set_ylim(ys[0], ys[1])
    # ax.xaxis.get_major_formatter().set_powerlimits((0,1))
    # ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
for i in range(1,5):
    ax = fig.add_subplot(3,5,6+i)
    ax.scatter(data[:5000,9], data[:5000,9+i], s=0.1)
    # if i == 1:
    #     ys = ax.set_ylim()
    #     xs = ax.set_xlim()
    # else:
    #     ax.set_xlim(xs[0], xs[1])
    #     ax.set_ylim(ys[0], ys[1])
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
for i in range(1,5):
    ax = fig.add_subplot(3,5,11+i)
    ax.scatter(snrs[0][:5000], snrs[i][:5000], s=0.1)
    # if i == 1:
    #     ys = ax.set_ylim()
    #     xs = ax.set_xlim()
    # else:
    #     ax.set_xlim(xs[0], xs[1])
    #     ax.set_ylim(ys[0], ys[1])
plt.savefig("F:/mag21_r0.8_btr_0.4_e(0.6-0.7)_scatter.pdf",bbox_inches='tight')
plt.show()
plt.close()
