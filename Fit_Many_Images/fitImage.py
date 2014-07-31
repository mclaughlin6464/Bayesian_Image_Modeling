desc ='''
This module contains the function fit_image. It can be run as main or imported.
The function fit_image fits to a fits image using emcee.
'''
import numpy as np
import emcee as mc
from time import time
from matplotlib import pyplot as plt
from scipy.stats import mode,chisquare
from multiprocessing import cpu_count

def gaussian(x,y, cx, cy, a, r):
    return a*np.exp(-1./(r**2)*((x-cx)**2+(y-cy)**2))
    
#theta contain the variables for the amplitude and width
#theta = [A1,A2...An,R1...Rn]
def lnprior(theta):
    #log uniform priors
    #simplified down to uniform priors, but these represent the exponents not the values themselves
    N = len(theta)/2 #save us from having to put N in the global scope
    amps = theta[:N]
    rads = theta[N:]
    if not all(-1<a<3 for a in amps):
        return -np.inf
    if not all(-1<r<2 for r in rads):
        return -np.inf
    return 0

def lnlike(theta, image, xx,yy,c_x, c_y,inv_sigma2):
    N = len(theta)/2
    amps = theta[:N]
    rads = theta[N:]

    model = np.zeros(image.shape) 
    for a,r in zip(amps, rads):
        if a<0: #if a<0 makes the amplitdue 0.
            continue 
        model+=gaussian(xx,yy,c_x, c_y, 10**a, 10**r)

    diff = image-model
    #largerThanImage = diff<0
    #diff[largerThanImage]*=10 #punish going over more than under

    return -.5*(np.sum(((diff)**2)*inv_sigma2-np.log(inv_sigma2)))

def lnprob(theta, image, xx, yy, c_x, c_y, inv_sigma2):
    lp = lnprior(theta)
    if np.isfinite(lp):
        return lp+lnlike(theta, image, xx, yy, c_x, c_y, inv_sigma2)
    return -np.inf

def cropImage(image, c_x,c_y, boxRatio = .03):
    #takes a passed in image and crops it around the passed in center
    #returns the image along with the new center points
    img_y, img_x = image.shape

    xLims = [c_x+img_x*boxRatio*(-1)**(1+i) for i in xrange(2)]
    yLims = [c_y+img_y*boxRatio*(-1)**(1+i) for i in xrange(2)]

    for i in xrange(2):
        if xLims[i] <0:
            xLims[i] = 0
        if xLims[i]>img_x:
            xLims[i] = img_x
        if yLims[i] <0:
            yLims[i] = 0
        if yLims[i]>img_y:
            yLims[i] = img_y

    image = image[yLims[0]:yLims[1]+1, xLims[0]:xLims[1]+1]
    #subtracting one seems to keep the centers on center.
    c_y,c_x = c_y-yLims[0], c_x-xLims[0]
    #find the max around the center
    dy, dx = np.unravel_index(image[c_y-1:c_y+2, c_x-1:c_x+2].argmax(), (3,3)) 
    dy,dx = dy-1, dx-1
    c_y, c_x = c_y+dy, c_x+dx
    return image, c_x, c_y

def calcsCounter(samples, N, decimals = 2):
#determine the calculated parameters using Counters
#esentially accepts the most common of the parameters
    from collections import Counter
    sampled_as = samples[:, :N].reshape((-1))
    sampled_rs = samples[:, N:].reshape((-1))

    Amp_counter = Counter(np.round(sampled_as, decimals = decimals).tolist())
    Rad_counter = Counter(np.round(sampled_rs, decimals = decimals).tolist())
    calc_as, a_counts= zip(*Amp_counter.most_common(N))
    calc_rs, r_counts = zip(*Rad_counter.most_common(N)) 

    return calc_as, calc_rs

def calcsSliceModes(samples, N, decimals = 2):
    #takes the modes of slices
    from scipy.stats import mode
    NGauss = 2
    sampled_as = np.round(samples[:, :N].reshape((-1)), decimals = decimals)
    sampled_rs = np.round(samples[:, N:].reshape((-1)), decimals = decimals)

    #uncomment to do corner plot of distributions
    #from triangle import corner
    #x = np.c_[sampled_as, sampled_rs]
    #fig = corner(x, labels = ['As', 'Rs'])
    calc_as, calc_rs = [],[]
    #i'm making a very large structure here. unecessary.
    sortSampledAs = np.sort(sampled_as)
    sortSampledRs = np.sort(sampled_rs)
    lengthOverN = sortSampledAs.shape[0]/float(NGauss)
    calc_as.extend([mode(sortSampledAs[n*lengthOverN:(n+1)*lengthOverN])[0][0] for n in xrange(NGauss)]) 
    calc_rs.extend([mode(sortSampledRs[n*lengthOverN:(n+1)*lengthOverN])[0][0] for n in xrange(NGauss)]) 

    return calc_as, calc_rs     

def calcsMixture(mc_samples, N, decimals = 2):
#determine the calculated parameters using the mixture package and expectation maximation
#attempst ot fit a mixture of Gaussians to the distributions.
    import mixture
    sampled_as = np.round(mc_samples[:, :N].reshape((-1)), decimals = decimals)
    sampled_rs = np.round(mc_samples[:, N:].reshape((-1)), decimals = decimals)

    #start by making 2 gaussians. Not worth making N yet.
    NGauss = 2
    calc_as, calc_rs = [], []
    for samples, calc in ((sampled_as,calc_as),(sampled_rs, calc_rs)):
        data = mixture.DataSet()
        data.fromArray(samples)
        sortSamples = np.sort(samples)
        lengthOverN = sortSamples.shape[0]/float(NGauss)
        mus = [sortSamples[n*lengthOverN:(n+1)*lengthOverN].mean() for n in xrange(NGauss)] 
        sigs =[sortSamples[n*lengthOverN:(n+1)*lengthOverN].std() for n in xrange(NGauss)] 

        gaussians = [mixture.NormalDistribution(mu, sig) for mu,sig in zip(mus, sigs)]

        model = mixture.MixtureModel(NGauss, [1./NGauss for n in xrange(NGauss)], gaussians)
        for g in model.components:
            g.distList[0].min_sigma = 0.001 #screw their safeguards!

        #randMaxEM was throwing assertion errors with no explanation
        #NOTE can get pick uncertainty on parameters with this.
        model.EM(data, 500, .001, silent = True) #perform the maximation
        sigs = []
        for g in model.components:
            calc.append(g.distList[0].mu)
            sigs.append(g.distList[0].sigma)
            
    return calc_as, calc_rs

def calcsCluster(samples, N, decimals = 2, n_clusters = 3):
#select the parameters with a clustering approach
    from scipy.stats import mode
    from sklearn.cluster import KMeans
    from triangle import corner

    sampled_as = samples[:, :N].reshape((-1))
    sampled_rs = samples[:, N:].reshape((-1))

    data = np.c_[sampled_as, sampled_rs]
    #reshape the data such that it's 2-D, with amps on one axis and radii on the other.

    n_clusters = N+1
    k_means = KMeans(init = 'k-means++',n_clusters = n_clusters, n_init = 50)
    labels = k_means.fit_predict(data)

    #plotting/debugging stuff
    '''
    colors = ['r', 'b', 'g', 'c','m','k']
    for i,color in enumerate(colors):
        g = data[labels==i]
        plt.scatter(g[:,0][:5000], g[:,1][:5000], color = color)
    '''
    roundData = np.round_(data, decimals = decimals)
    clusters = [roundData[labels == i] for i in xrange(n_clusters)]

    #get a list of the highest points and their average count
    allModes = []
    for cluster in clusters:
           point, counts = mode(cluster)
           allModes.append((point[0][0], point[0][1], counts.mean()))

    allModes.sort(key = lambda x:x[-1], reverse = True) #sort by highest count
    allModes = np.array(allModes)
    modes = allModes[:N,:2]

    calc_as, calc_rs = modes[:,0], modes[:,1]
    '''
    plt.scatter(calc_as, calc_rs, s = 40, color = 'y', label = 'Mode')
    plt.legend(loc = 3)
    plt.title('Cluster Fit')

    figs= [corner(data, truths = modes[i]) for i in xrange(N)]

    plt.show()
    '''
    return calc_as, calc_rs 

def fitImage(image, N, c_x=-1, c_y=-1, n_walkers = 600, crop = True, ddof = 0, filename = None):

    np.random.seed(int(time()))

    if -1 in (c_x, c_y) : #one's the default
        c_y, c_x = np.unravel_index(image.argmax(), image.shape)

    if crop: #recrop if necessary
        image, c_x, c_y = cropImage(image, c_x, c_y)
    img_y, img_x = image.shape

    #y,x = [np.arange(i) for i in image.shape]
    #xx, yy = np.meshgrid(x,y)
    yy, xx = np.indices(image.shape)

    err = .1
    inv_sigma2 = 1./(err**2)

    ndim = N*2
    nburn = n_walkers*.1
    nsteps = 200

    pos = []
    for walk in xrange(n_walkers):
        row = np.zeros(ndim)
        for n in xrange(N):
            row[n] = 4*np.random.rand()-1
            row[n+N] = 2*np.random.rand()-1
        pos.append(row)

    args = (image, xx, yy, c_x, c_y, inv_sigma2) 
    sampler = mc.EnsembleSampler(n_walkers, ndim, lnprob, args = args, \
                                                    threads = cpu_count()) 
    sampler.run_mcmc(pos, nsteps)
    samples = sampler.chain[:,nburn:,:].reshape((-1, ndim))
    sampler.pool.terminate()
    del(sampler)
    #should I save all to file? not a bad idea.
    #a lot of memory though
    #if filename is not None:
    #    np.savetxt('Samples/'+filename, samples, delimiter = ',')

    calc_as, calc_rs = calcsCluster(samples, N)
    del(samples)

    calc_img = sum(gaussian(xx,yy,c_x,c_y,10**a,10**r) for a,r in zip(calc_as, calc_rs))
    #Matias told me that the DoF is the number of pixels PLUS the number of parameters. Not sure if that's right but I'll give it a go.
    ddof = -2*N + ddof 
    chi2stat, b = chisquare(image, f_exp = calc_img, ddof = ddof,axis = None)

    #TODO the crop function makes returing the original unecessary
    return image, calc_img, chi2stat

if __name__ == '__main__':
    import argparse
    import pyfits

    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('filename', metavar = 'fname', type = str, help = 'The name of the fits file to be read in.')
    parser.add_argument('nGaussians', metavar = 'N', type = int, help = 'Number of Gaussians to use in the fit.', default = 2)
    parser.add_argument('center_x', metavar = 'cx', type = int, help = 'The center in x', default = -1)
    parser.add_argument('center_y', metavar = 'cy', type = int, help = 'The center in y', default = -1)

    args = parser.parse_args()

    filename = args.filename

    try:
        fitsImage = pyfits.open(filename)
    except IOError:
        print 'ERROR: Invalid filename.'
        from sys import exit
        exit(-1)

    image = fitsImage[0].data
    c_y, c_x = args.center_y, args.center_x

    N = args.nGaussians
    image, calc_img, chi2stat = fitImage(image, N, c_x, c_y, 1000) 
    plt.subplot(131)
    plt.imshow(image)
    plt.colorbar()
    plt.subplot(132)
    plt.imshow(image-calc_img)
    plt.colorbar()
    plt.subplot(133)
    plt.imshow(calc_img)
    plt.colorbar()
    plt.show()
