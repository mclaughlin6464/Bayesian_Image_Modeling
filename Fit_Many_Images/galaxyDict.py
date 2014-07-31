'''
This module is part of a larger package that fits to a larget catalog of images. 
This module contains a function that takes the filename of a catalog that marks the center of the lensed galaxy in each image. It reads in the file and createss a dictionary of the x-y points and returns it. 
'''

def getGalaxyDict(filename):
    
    galaxyDict = {}
    try:
        with open(filename) as f: 
            for line in f:
                splitLine = line.strip().split(' ')
                coords = [float(x) for x in splitLine[2:]]
                galaxyDict[splitLine[0]] = coords

    except IOError:
        print 'Invalid filename passed into getGalaxyDict'
        raise

    return galaxyDict
                    
