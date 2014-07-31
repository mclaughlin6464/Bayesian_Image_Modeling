#Generates a 3D contour plot of the passed in image

import pyfits

from matplotlib import pyplot as pyplt
from matplotlib.colors import LogNorm as lm
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

from sys import argv, exit

from math import log10

def threeDplot(image): #plots a 3D contour plot of the fits image

	fig = pyplt.figure()
	ax = fig.gca(projection = '3d')
	
	ylen, xlen = image.shape	

	threshold = image.mean() + image.std()

	X = np.arange(0, xlen,1)
	Y = np.arange(0, ylen, 1)
	X,Y = np.meshgrid(X,Y)

	Z = [ [ 0 for y in xrange(ylen)] for x in xrange(xlen)] 

	for i in xrange(xlen):
		for j in xrange(ylen):
			if image[j,i] <threshold:
				continue
			Z[i][j] = log10(image[j,i]+10)

	Z = np.array(Z)

	surf = ax.plot_surface(X,Y,Z, rstride =1, cstride = 1, cmap=cm.coolwarm, 
			linewidth = 0, antialiased = False, vmin = log10(image.mean()+image.std()+10))

	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

	#ax.set_zscale('log')

	fig.colorbar(surf, shrink = .5, aspect = 5)

	pyplt.show()


usage = '''
USAGE: python %s image

WHERE 

image is the name of a fits file of the image of interest

PURPOSE:

This code creates a 3-D contour plot in matplotlib from the passed in image.
'''

if __name__ == '__main__':
	if len(argv) !=2:
		print usage&(argv[0])

	else:

		try:

			fitsImage = pyfits.open(argv[1])
			image = fitsImage[0].data

			threeDplot(image)

			fitsImage.close()

		except IOError:
			print "ERROR: file", argv[1], "does not exist.\n"
			exit(-1)
