'''
I'm simplifying to a 1-D fit
Since I successfully fit to one Gaussian, I'm going to fit to a sum now
'''

import numpy as np
import pymc as pm
from scipy.stats import norm
from time import time
t0 = time()

SIZE = 100
N = 10 
DEC = 1

np.random.seed(1)
true_mean = np.random.randint(20,81) 
true_sigs = [np.random.rand()*10+5 for n in xrange(N)]
true_weights = np.asarray([np.random.randint(1,10) for n in xrange(N)])

x = np.arange(SIZE)
image = np.zeros(SIZE/DEC)
for n in xrange(N):
	image += true_weights[n]*(norm(loc = true_mean,scale = true_sigs[n]).pdf(np.arange(0,SIZE, DEC))) 

#image = image/image.sum() #normalize

#add some noise
NOISE_LEVEL = .005*true_weights.mean()
noise = np.random.normal(size = SIZE/DEC)*NOISE_LEVEL

image = image+noise

sigs = pm.Chi2('sigs',7, size = N)

mu = pm.Uniform('mu', 0, SIZE, value = SIZE/2.0 )

weights = pm.Uniform('weights', 1, 10, value = [5 for i in xrange(N)], size = N)

#model
x = np.arange(0,SIZE,DEC)
@pm.deterministic
def normal(x = x, mu = mu, sigs = sigs, weights = weights):
	return  sum(weights[i]*norm(loc = mu, scale = sigs[i]).pdf(x) for i in xrange(N))

tau = pm.Chi2('tau', 4)

y = pm.Normal('y', normal, tau, value = image, observed = True) 

model = pm.Model([y, tau, normal, mu, sigs, weights])

map_ = pm.MAP(model)
map_.fit()

R = pm.MCMC(model)
	
R.sample(10000, 2000, 2)

from matplotlib import pyplot as plot

sig_means = []
weight_means = []

mu_trace = R.trace('mu')[:]
mu_mean = mu_trace.mean()
sig_trace = R.trace('sigs')[:]
weight_trace = R.trace('weights')[:]
for i in xrange(N):
	sig_means.append((sig_trace[:,i]).mean())
	weight_means.append((weight_trace[:,i]).mean())

for i in xrange(N):
	print mu_mean, sig_means[i] , weight_means[i]
	print true_mean, true_sigs[i], true_weights[i] 
	print 

plot.scatter(x,image, color = 'r')

plot.plot(x, sum(weight_means[i]*norm(loc= mu_mean, scale=sig_means[i]).pdf(x) for i in xrange(N)), color = 'b') 

print 'Time: %f'%(time()-t0)
plot.show()

#pm.Matplot.plot(R, common_scale = False)

