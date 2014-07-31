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

args = parser.parse_args()

import pyfits
import numpy as np
import pymc as pm 
from pymc.Matplot import plot as mcplot

try:
	fitsImage = pyfits.open(args.filename)
except IOError:
	print 'ERROR: Invalid filename'
	from sys import exit
	exit(-1)

image = fitsImage[0].data
img_y, img_x = image.shape

norm_image = np.zeros(image.shape)
total_flux = float(image.sum())
for i in xrange(img_x):
	for j in xrange(img_y):
		norm_image[j,i] = image[j,i]/total_flux
		
image = norm_image 

digitsOfPrecision = 1
N = 10**digitsOfPrecision
img_coords = np.array([[0,0] for i in xrange(N)]) 
idx = 0
for y in xrange(img_y):
	for x in xrange(img_x):
		for i in xrange(int(image[y,x]*N)):
			img_coords[idx] = [x,y]
			idx+=1

mean_x = pm.Uniform("mean_x", 0, img_x)
mean_y = pm.Uniform('mean_y',0, img_y)

@pm.deterministic
def mean(mean_x = mean_x, mean_y = mean_y):
	return np.array([mean_x, mean_y])

tau = pm.Wishart('tau', n=4, Tau = np.eye(2))

data = pm.MvNormal("data", mean, tau, value =img_coords , observed = True)

model= pm.Model([mean,mean_x, mean_y, tau, data])

map_ = pm.MAP(model)
map_.fit()

mcmc = pm.MCMC(model)
mcmc.sample(20000,burn = 2000, thin = 2)

mean_samples = mcmc.trace('mean')[:]
tau_samples = mcmc.trace('tau')[:]

mcplot(mcmc, common_scale = False)
