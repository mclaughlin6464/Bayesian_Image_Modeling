'''
This module models a galaxy in a fits image with a MCMC approach.

This current version is a test of a simple case.
@author Sean McLaughlin
@date 6/5/2014
'''

desc = '''
This module takes a fits file as an input and creates a model using MCMC.
'''
#process imputs cleanly
import argparse
parser = argparse.ArgumentParser(description = desc)
parser.add_argument('filename', metavar = 'fname', type = str, help = 'The filename')
parser.add_argument('nGaussians', metavar = 'n', type = int, help = 'The number of Gaussians to use in this fit')
parser.add_argument('digitsOfPrecision', metavar = 'DoP', type = int, help = 'The digits of precision to use for the Gaussians')

args = parser.parse_args()

nGaussians = args.nGaussians
digitsOfPrecision = args.digitsOfPrecision

import pyfits
import numpy as np
import pymc as pm 
from pymc.Matplot import plot as mcplot

#attempt to open the file
try:
	fitsImage = pyfits.open(args.filename)
except IOError:
	print 'ERROR: Invalid filename'
	from sys import exit
	exit(-1)

#below I take the image, and normalize it so that the total sum of all pixels is 1
image = fitsImage[0].data
img_y, img_x = image.shape

norm_image = np.zeros(image.shape)
total_flux = float(image.sum())
for i in xrange(img_x):
	for j in xrange(img_y):
		norm_image[j,i] = image[j,i]/total_flux
		
image = norm_image 

#this part is a bit hackish
#this inserts the coordinates into img_coords with a weight according to their brightness
#this allows the variable to read the data
N = 10**digitsOfPrecision
img_coords = np.array([[0,0] for i in xrange(N)]) 
idx = 0
for y in xrange(img_y):
	for x in xrange(img_x):
		for i in xrange(int(image[y,x]*N)):
			img_coords[idx] = [x,y]
			idx+=1

dd = pm.Dirichlet('dd', theta = (1,)*nGaussians)
category = pm.Categorical('category', p = dd, size = nGaussians)

mean_x = pm.Uniform("mean_x", 0, img_x, size = nGaussians)
mean_y = pm.Uniform('mean_y',0, img_y, size = nGaussians)

@pm.deterministic
def means(mean_x = mean_x, mean_y = mean_y):
	return np.array([mean_x, mean_y])

@pm.deterministic
def mean_i(means = means, category = category):
	return means[category]

taus_ = []
for i in xrange(nGaussians):
	taus_.append(pm.Wishart('tau%d'%i, n=4, Tau = np.eye(2))) 

taus_= np.asarray(taus_)
@pm.deterministic
def taus(taus = taus_):
	return taus

@pm.deterministic
def tau_i(taus = taus, category = category):
	return taus[category]

data = pm.MvNormal("data", mean_i, tau_i, value =img_coords , observed = True)

model= pm.Model([means, mean_i,mean_x, mean_y, taus, tau_i, data])

mcmc = pm.MCMC(model)
mcmc.sample(20000,burn = 2000, thin = 2)

mean_samples = mcmc.trace('mean')[:]
tau_samples = mcmc.trace('tau')[:]

mcplot(mcmc, common_scale = False)
