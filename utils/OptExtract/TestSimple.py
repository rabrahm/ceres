# -*- coding: utf-8 -*-
import pyfits
from pylab import *
import Marsh
import numpy
import scipy

def getSpectrum(filename,b,Aperture,minimum_column,maximum_column):

    hdulist = pyfits.open(filename)                  # Here we obtain the image...
    data=hdulist[0].data                             # ... and we obtain the image matrix.
    Result=Marsh.SimpleExtraction((data.flatten()).astype(double),scipy.polyval(b,numpy.arange(data.shape[1])).astype(double),data.shape[0],data.shape[1],data.shape[1],Aperture,minimum_column,maximum_column)
    FinalMatrix=asarray(Result)                      # After the function, we convert our list to a Numpy array.
    return FinalMatrix

# Main function:

FileName='../../../transmission_spectroscopy/WASP6/data.fits'                                          # Filename of the image of the spectrum.
b=pyfits.getdata('../../../transmission_spectroscopy/WASP6/trace_coeffs.fits')                                        # We write our trace...
Aperture=15
Spectrum=getSpectrum(FileName,b,Aperture,0,0)                # In this example we take our spectrum from x=50 to 1500,
                                                                 # where x is in our matrix coordinates.
x=arange(len(Spectrum))
plot(x,Spectrum,'-')                                             # Plot the results.
show()
