import numpy as np
import matplotlib.pyplot as plt

#####################################################
# There are some demos to show how the matpotlib works
# delete the '#' ahead of the codes of each demo to see how it works
#####################################################

# Demo 1
# step 1
# let's show the 'x+y' array
# my and mx are the corresponding coordinates array of y-axis and x-axis
# you can put any 2-d array into the 'imshow('any 2-d array')
#####################################################

# my,mx =np.mgrid[0:15,0:15]
# plt.imshow(my+mx)
# plt.colorbar()
# plt.show()

#####################################################
# step 2
# show a Gaussian function e^(-r^2/2), r^2 = (x-7)^2+(y-7)^2,
# of which the center is at (7,7)
#####################################################

# my,mx =np.mgrid[0:15,0:15]
# gaussian = np.exp(-((mx-7)**2+(my-7)**2)/2)
# plt.imshow(gaussian)
# plt.colorbar()
# plt.show()

#####################################################

# Demo 2
# let's show the my and mx separately
# google the method of subplot(), i.e. assembling sub-figures into one figure
# take notice of the number in the bracket of subplot().
#####################################################

# my,mx =np.mgrid[0:15,0:15]
# plt.subplot(121)
# plt.imshow(mx)
# plt.colorbar()
# plt.subplot(122)
# plt.imshow(my)
# plt.colorbar()
# plt.show()

#####################################################

# Demo 3
# step 1
# let's plot some line, i.e. y = 2*x+10
# x is a 1-d array
# y is a 1-d array
# take notice of the sequence of x and y in plot()
# you can specify the color ,style of the line style and other parameters
# google 'matplotlib' and search 'plot' for details
#####################################################

# x = np.linspace(1,10,10)
# y = 2*x+10
# plt.plot(x,y,color='g')
# plt.show()

#####################################################
# step 2
# plot the dots not the line
#####################################################

# plt.scatter(x,y)
# plt.show()

#####################################################