'''
I'm simplifying to a 1-D fit
'''

import numpy as np
import pymc as pm
from scipy.stats import norm

SIZE = 100
CENTER = SIZE/2
RAD = .25

def circle(x):
	a = (SIZE*RAD)**2-(x-CENTER)**2 
	return a if a>0 else 0

image = np.asarray([ circle(x) for x in np.arange(0,SIZE, .5)])
image = 2*image/image.sum() #normalize
'''
import time
np.random.seed(int(time.time()))
true_mean = np.random.randint(20, 81) 
true_sig = np.random.rand()*8+2 
image = norm(loc = true_mean,scale = true_sig).pdf(np.arange(0,SIZE, .5)) 

#add some noise
NOISE_LEVEL = .005
noise = np.random.normal(size = SIZE*2)*NOISE_LEVEL

image = image+noise
'''

sig = pm.Chi2('sig',7)

mu = pm.Uniform('mu', 0, SIZE,value = SIZE/2.0)

#model
x = np.arange(0,SIZE,.5)
@pm.deterministic
def normal(x = x, mu = mu, sig = sig):
	return  norm(loc = mu, scale = sig).pdf(x)	

tau = pm.Chi2('tau', 4)

print normal.value.shape
print tau.value.shape
print image.shape
y = pm.Normal('y', normal, tau, value = image, observed = True) 

model = pm.Model([y, tau, normal, mu, sig])

map_ = pm.MAP(model)
map_.fit()

R = pm.MCMC(model)
	
R.sample(50000, 2000)

from matplotlib import pyplot as plot

mu_trace = R.trace('mu')[:]
sig_trace = R.trace('sig')[:]
tau_trace = R.trace('tau')[:]

mu_mean = mu_trace.mean(axis = 0)
sig_mean= sig_trace.mean(axis = 0)
tau_mean = tau_trace.mean(axis =0)
'''
print mu_mean, sig_mean, 1/tau_mean
print true_mean, true_sig, NOISE_LEVEL 
'''

plot.scatter(x,image, color = 'r')

total = 0
b = norm(loc = mu_mean, scale=sig_mean)
for i in np.arange(0, 100, .5):
	total+=b.pdf(i)
print total

plot.scatter(x, norm(loc= mu_mean, scale=sig_mean ).pdf(x), color = 'b', alpha = .6) 

plot.show()

pm.Matplot.plot(R, common_scale = False)
