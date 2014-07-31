#! /usr/bin/env python

# This program produces the stellar orbit plot for Homework #2, 
# Problem 10, Astronomy 596, Spring 2011 semester.

# Robert J. Brunner

import pyfits

def center(image):

    threshold = image.mean() + image.std()
    
    xmean = ymean = totval = 0
    
    ylen, xlen = image.shape
    
    x_pixels = []
    y_pixels = []

    for i in range(xlen):
        for j in range(ylen):
            if (image[j, i] > threshold):
                
                xmean += i * image[j, i]
                ymean += j * image[j, i]
                totval += image[j, i]      
                
                x_pixels.append(i)
                y_pixels.append(j)
    
    # Note could handle div by zero here, but we have bigger 
    # problems if no pixels are above threshold.

    return(xmean/float(totval), ymean/float(totval), x_pixels, y_pixels)

# This is an extra credit function, to calculate ellipticity and
# orientation angle

from math import sqrt, atan, degrees, cos, sin, pi

def moments(image):

    threshold = image.mean() + image.std()
    
    # We cast to long here to avoid integer overflow during the summation process.
    Ix = Iy = Ixx = Iyy = Ixy = I = long(0)
    
    ylen, xlen = image.shape

    for i in range(xlen):
        for j in range(ylen):
            if (image[j, i] > threshold):
                I += image[j, i]      
                Ix += i * image[j, i]
                Iy += j * image[j, i]
                Ixx += i * i * image[j, i]
                Iyy += j * j * image[j, i]
                Ixy += i * j * image[j, i]


    # Note could handle div by zero here, but we have bigger 
    # problems if no pixels are above threshold.

    xmean = Ix / I
    ymean = Iy / I
    
    Uxx = Ixx/I - xmean **2
    Uyy = Iyy/I - ymean **2
    Uxy = Ixy/I - xmean * ymean
    
    theta = 0.5 * atan(2.0 * Uxy / (Uxx - Uyy))

    lambda1 = 0.5 * (Uxx + Uyy) + 0.5 * sqrt(4 * Uxy**2 + (Uxx - Uyy)**2)
    lambda2 = 0.5 * (Uxx + Uyy) - 0.5 * sqrt(4 * Uxy**2 + (Uxx - Uyy)**2)

    print 'a = %5.2f' % sqrt(lambda1)
    print 'b = %5.2f' % sqrt(lambda2)

    print degrees(theta), sqrt(1 - lambda2 / lambda1)
    
    return (xmean, ymean, sqrt(lambda1), sqrt(lambda2), theta)
    
    # this values can be verified in a number of places, but are
    # essentially solutions
    # www.stsci.edu/jwst/externaldocs/technicalreports/JWST-STScI-000396.pdf

    # term1 = 0.5 * (xx + yy)
    # term2 = 0.5 * sqrt((xx - yy)**2 + 4.0 * xy**2)
	
    # a = term1 - term2
    # b = term1 + term2

    # ellipticity = sqrt((xx - yy)**2 + 4.0 * xy**2)/(xx + yy)
    # angle = 0.5 * atan(2.0 * xy/(xx - yy))

    # print a, b, (1. - b/a), ellipticity
    # return(ellipticity, angle)

from matplotlib import pyplot as pyplt
from matplotlib.colors import LogNorm as lm
from matplotlib import cm as cm

def my_plot(image): 

    (xc, yc, x, y) = center(image)

    (xc, yc, a, b, theta) = moments(image)

    fig = pyplt.figure()
    plt = fig.add_subplot(111)

    # Note I am using Latex in the plot labels
    
    pyplt.scatter(x, y)
        
    plt.scatter(xc, yc, c='r')

    plt.set_xlabel("X")
    plt.set_ylabel("Y")
        
    plt.set_title("ASTR 406: Informatics #2, Problem #3EC (Image Pixel Plot)")

    my,mx = image.shape
    plt.set_xlim(0, mx)
    plt.set_ylim(0, my)

    pyplt.plot((xc - a * sin(theta), xc + a * sin(theta)), \
        (yc - a * cos(theta), yc + a * cos(theta)), color='g')

    pyplt.plot((xc - b * cos(theta), xc + b * cos(theta)), \
        (yc + b * sin(theta), yc - b * sin(theta)), color='y')
        
    # Uncomment next line to show plot on screen
    #pyplt.show()

    #pyplt.savefig("info2-3ec-Pixels.pdf")
    
def my_image_plot(image):

    (xc, yc, a, b, theta) = moments(image)

    fig = pyplt.figure()
    plt = fig.add_subplot(111)

    # Here we display the image. Note that I am using a Log
    # Normalization, where we map the image mean to 0 and the image max to
    # 1, and clip any points outside this range (any points less than the
    # mean should be set to zero. This makes the image dynamic range more
    # realistic. If I knew how, I would use a histogram equalization
    # technique.

    pyplt.imshow(image, \
        norm=lm(image.mean() + 0.5 * image.std(), image.max(), clip='True'), \
        cmap=cm.gray, origin="lower")
        
    plt.scatter(xc, yc, c='r')

    plt.set_xlabel("X")
    plt.set_ylabel("Y")


    pyplt.plot((xc - a * sin(theta), xc + a * sin(theta)), \
        (yc - a * cos(theta), yc + a * cos(theta)), color='g')

    pyplt.plot((xc - b * cos(theta), xc + b * cos(theta)), \
        (yc + b * sin(theta), yc - b * sin(theta)), color='y')
        
    plt.set_title("ASTR 406: Informatics #2, Problem #3EC (Image Center Plot)")

    my,mx = image.shape
    plt.set_xlim(0, mx)
    plt.set_ylim(0, my)
    
    # Uncomment next line to show plot on screen
    pyplt.show()

    #pyplt.savefig("info2-3ec-Image.pdf")

# General Main method.

if __name__ == '__main__':

    from sys import argv, exit
    
    usage = """
    USAGE: python %s image

    WHERE

    image is the name of the fits file containing the image of interest
    (e.g., M86.fits)

    PURPOSE:
	
	This code thresholds an image, and calulates the image center by
	using a weighted mean algorithm. In addition, all pixels that are
	above the threshold value are collected and plotted as a scatter
	plot (to approximate the image shape).
    
    """

    # Now if we do not call the function correctly, we should exit and print
    # useful message

    if (len(argv) != 2):
        print usage % (argv[0])
 
    else:
    
        try:
            fitsImage = pyfits.open(argv[1])
            image = fitsImage[0].data
               
            my_plot(image)
            my_image_plot(image)

            fitsImage.close()

        except IOError:
            print "ERROR: file ", argv[1], "does not exist.\n"
            exit(-1) # So exit with an error code to the system.    	

#    	my_plot(x, y)


            
