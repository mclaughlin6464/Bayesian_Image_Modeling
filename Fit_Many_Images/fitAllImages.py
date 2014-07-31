desc='''
This is the main module for this package. It performs an MCMC fit on many simulated lens images.

Input: The name of the directory where the images are located.
       The name of the catalog file.
       The name of the directory to save the output plots. 

Output: For each collection of 5 ugriz images, produces 1 png image.
        The image will hold gri crops of the galaxy in question.
        It will also show 2 fits, and the difference between them.
'''
import matplotlib
#matplotlib.use('Agg') #allows to be run on server.
from fitImage import fitImage, cropImage
from galaxyDict import getGalaxyDict
import os
import argparse
import pyfits
from matplotlib import pyplot as plt

#get inputs cleanly
parser = argparse.ArgumentParser(description = desc)
parser.add_argument('imageDirectory', metavar = 'dName', type = str, help = 'The name of the directory holding the fits files.')
parser.add_argument('catalogFilename', metavar = 'cfname', type = str, help = 'The name of the catalog file for the images.')
parser.add_argument('outputDirectory', metavar = 'outputName', type = str, \
                         help = 'The location to store the output images.')

args = parser.parse_args()

#get the locations of the centers of the galaxies.
galaxyDict = getGalaxyDict(args.catalogFilename)
#TODO check all are valid directories, esp output.
#output isn't run until after 2 fits.
imageDirectory = args.imageDirectory
outputDirectory = args.outputDirectory
if imageDirectory[-1] != '/':
    imageDirectory+='/'
if outputDirectory[-1] != '/':
    outputDirectory+='/'
   
fileList = os.listdir(imageDirectory)
baseNames = set()#set prevents duplicates
trackedObjs = []
for fName in fileList:
    #NOTE may have to change when moving to other files/types
    if fName[:6] == 'CFHTLS': #one of the ones from the simulation.
        baseNames.add(fName[:-7])#take off the band and suffix

#sort here?
baseNames = list(baseNames)
for baseName in baseNames:
    try:
        fullCx, fullCy = galaxyDict[baseName[7:]] #get the centers
        cx, cy = 0,0
        images = {}
        bands = ['g', 'r', 'i']
        #get the images and crop them
        for band in bands:
            fitsImage = pyfits.open(imageDirectory + baseName+ '_'+band+'.fits')
            fullImage = fitsImage[0].data
            image, cx, cy = cropImage(fullImage, fullCx, fullCy)
            images[band] = image

        #pefrom the fits on the g and i bands
        nGaussians = 4
        #measure by the statistic, over multiple N's?
        g_image,g_fit,g_stat = fitImage(images['g'],nGaussians,cx,cy,1000, crop = False, filename = baseName+'_g')
        i_image,i_fit,i_stat = fitImage(images['i'],nGaussians,cx,cy,1000, crop = False, ddof = -1, filename = baseName+'_i')
        #this conversion fails if the center of the g image is not the max.
        #second version will use the centers as reference.
        #i_fit = i_fit*images['g'].max()/images['i'].max()#rescale the i fit.
        c = (int(cy), int(cx))
        i_fit = i_fit*images['g'][c]/images['i'][c]#rescale the i fit.

        #below, plot and save.
        fig = plt.figure()
        fig.suptitle(baseName)
        for i, band in enumerate(bands):
            plt.subplot(2,3,i+1)
            plt.title(band+'-band original')
            im = plt.imshow(images[band])    
            plt.scatter(cx, cy, c = 'k')
            plt.colorbar(im)
        
        plt.subplot(2,3,4)
        plt.title('g-band fit')
        im =plt.imshow(images['g']-g_fit)
        plt.scatter(cx, cy, c= 'k')
        plt.colorbar(im) 
        plt.subplot(2,3,5)
        plt.title('i-band fit')
        im =plt.imshow(images['g']-i_fit)
        plt.scatter(cx, cy, c= 'k')
        plt.colorbar(im) 
        plt.subplot(2,3,6)
        plt.title('g-fit - i-fit')
        im =plt.imshow(g_fit-i_fit)
        plt.scatter(cx, cy, c= 'k')
        plt.colorbar(im)
        plt.show()
        #plt.savefig(outputDirectory+baseName+'_fit.png')
        #plt.clf()
    except KeyboardInterrupt:
        raise
    
    except Exception as e: #catch any error
        print "ERROR %s with %s"%(e,baseName)
        raise
        #continue #keep going.
