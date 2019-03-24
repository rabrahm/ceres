# -*- coding: utf-8 -*-
import pyfits
from pylab import *
import Marsh
import numpy
import scipy
from scipy.signal import medfilt

def getSpectrum(filename,b,Aperture,minimum_column,maximum_column):

    hdulist = pyfits.open(filename)                  # Here we obtain the image...
    data=hdulist[0].data                             # ... and we obtain the image matrix.
    D=data
    Result=Marsh.SimpleExtraction((data.flatten()).astype(double),scipy.polyval(b,numpy.arange(D.shape[1])).astype(double),data.shape[0],data.shape[1],data.shape[1],Aperture,minimum_column,maximum_column)
    FinalMatrix=asarray(Result)                      # After the function, we convert our list to a Numpy array.
    return FinalMatrix


FileName='Image.fits'
FN = 'Background.fits'
O=pyfits.getdata('Orders.fits')
Aperture=5.0
D=pyfits.getdata(FileName)
Spectrum=getSpectrum(FileName,O[30],Aperture,0,len(D[0,:])) # Simple spectrum of the object.
BSpectrum=getSpectrum(FN,O[30],Aperture,0,len(D[0,:]))      # Simple spectrum of the background.
x=arange(len(Spectrum))
#plot(x,Spectrum,'r-')
#plot(x,BSpectrum,'b-')
#M = medfilt(BSpectrum,int(sqrt(len(BSpectrum)))+1)
#plot(x,M,color='black')
#plot(S[0],S[1],'g-')                                             # Plot the results.
#show()
