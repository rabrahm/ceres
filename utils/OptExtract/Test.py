# -*- coding: utf-8 -*-
import pyfits
from pylab import *
import Marsh
import numpy
import scipy

"""
----------------------------------------------------------------
The P Function.

Returns an array with the P[i][j] coefficients, the modeled
light fractions.
----------------------------------------------------------------

Here we pass the matrix, coefficients and lengths to the C programmed module. To call this
function of the module we have to do:
 
  Marsh.ObtainP(M,T,MRows,MCols,TLength,Aperture,RON,GAIN,NSigma,S,N,MODE)
  
Where the inputs are:
  
  M:       Matrix of the image containing the pixel values with the spectrum's wavelength direction on the column direction.
           It is assumed that the image passed is sky-substracted. 
  
  T:       Trace vector's coefficients (D_N,D_(N-1),D_(N-2)...,D_0) where the trace vector is in the form sum(D_n*x^n), n=0...N.
           The trace vector is assumed to be on matrix coordinates, with the horizontal 0 (i.e. x/column 0) on the first horizontal
           pixel (from left to right) to have spectrum data.
	   
  MRows:   Number of rows of the M matrix.
  
  MCols:   Number of columns of the M matrix.
  
  TLength: Number of elements (length) of the trace vector that contains the coefficients.
  
  Aperture:Number of pixels, from the center pixel, of the aperture radius.
  
  RON:     Read Out Noise of the measurements (in electrons).
  
  GAIN:    Gain of the CCD in electrons/ADU
  
  NSigma:  Number of times we multiply sigma to get a confidence interval.
  
  S:       The spacing between the polynomials to be fitted (see Marsh, 1989). Set this to one
           if Horne's algorithm is going to be applied.
  
  N:       The number of coefficients of the polynomials to be fitted (see Marsh, 1989, Horne,
           1986).
           
  MODE:    Set this to 0 if Marsh's algorithm is desired (curved trace). Any other number will
           produce Horne's algorithm to be applied to the data.
	   
  min_col: Pixel coordinates of the input data to start operating from. If zero, the MRows input is used instead.
  
  max_col: Pixel coordinates of the input data to end operating from.  If zero, the MCols input is used instead.

The output is a flattened python list that we convert to a Numpy Array (with the asarray() function) and then 
convert in matrix form doing a resize.
"""

def PCoeff(filename,b,Aperture,RON,Gain,NSigma,S,N,marsh_alg,min_col,max_col):
 hdulist = pyfits.open(filename)                  # Here we obtain the image...
 data=hdulist[0].data                             # ... and we obtain the image matrix.
 Result=Marsh.ObtainP((data.flatten()).astype(double),scipy.polyval(b,numpy.arange(data.shape[1])).astype(double),data.shape[0],data.shape[1],data.shape[1],Aperture,RON,Gain,NSigma,S,N,marsh_alg,min_col,max_col)
 FinalMatrix=asarray(Result)                      # After the function, we convert our list to a Numpy array.  
 FinalMatrix.resize(data.shape[0],data.shape[1])  # And return the array in matrix-form.
 return FinalMatrix

"""
----------------------------------------------------------------
The getSpectrum function.

Returns a matrix with three rows:

	First row:        Pixel coordinates (taken from the first horizontal/column/x point with data) of the fluxes and variances.
	Second row:       Flux values for the coordinates on the first row, same column.
	Third row:        Variances of the respective flux values.
----------------------------------------------------------------
Here we pass the image matrix, the P Coefficients, coefficients of trace vector and lengths to C. To call this
function we do:
 
  Marsh.ObtainSpectrum(M,T,P,MRows,MCols,TLength,Length,RON,GAIN,S,CSigma)
  
Where the inputs of the ObtainSpectrum function are:
  
  M:       Matrix of the image containing the pixel values with the spectrum's wavelength direction on the column direction.
  
  T:       Trace vector's coefficients (D_N,D_(N-1),D_(N-2)...,D_0) where the trace vector is in the form sum(D_n*x^n), n=0...N.
           The trace vector is assumed to be on matrix coordinates, with the horizontal 0 (i.e. x/column 0) on the first horizontal
           pixel (from left to right) to have spectrum data.

  P:       Coefficients of each image value obtained on the ObtainP function. 
	   
  MRows:   Number of rows of the M/P matrix.
  
  MCols:   Number of columns of the M/P matrix.
  
  TLength: Number of elements (length) of the trace vector that contains the coefficients.
  
  Length:  Number of pixels, from the center pixel, of the aperture radius.
  
  RON:     Read Out Noise of the measurements (in electrons).
  
  GAIN:    Gain of the CCD in electrons/ADU
  
  S:       The spacing between the fitted polynomials. Set this to one
           if Horne's algorithm was applied.
  
  CSIgma:  Number, in sigmas, for the cosmic ray rejection.
  
  min_col: Pixel coordinates of the input data to start operating from. If zero, the MRows input is used instead.
  
  max_col: Pixel coordinates of the input data to end operating from.  If zero, the MCols input is used instead.

The output is a flattened python list that we convert to a Numpy Array (with the asarray() function) and then 
convert in matrix form doing a resize.
"""
 
def getSpectrum(P,filename,b,Aperture,RON,Gain,S,NCosmic,min_col,max_col):
  
 hdulist = pyfits.open(filename)                  # Here we obtain the image...
 data=hdulist[0].data                             # ... and we obtain the image matrix.
 Result,size =Marsh.ObtainSpectrum((data.flatten()).astype(double),scipy.polyval(b,numpy.arange(data.shape[1])).astype(double),P.flatten(),data.shape[0],data.shape[1],data.shape[1],Aperture,RON,Gain,S,NCosmic,min_col,max_col)
 FinalMatrix=asarray(Result)                      # After the function, we convert our list to a Numpy array.  
 FinalMatrix.resize(3,size)                       # And return the array in matrix-form.
 return FinalMatrix

# Main function:

FileName='../../../transmission_spectroscopy/WASP6/data.fits'                                          # Filename of the image of the spectrum.
b=pyfits.getdata('../../../transmission_spectroscopy/WASP6/trace_coeffs.fits')                                        # We write our trace...
Aperture=15                                                             # Aperture, in pixels, from the center.
RON=4.9                                                                 # Read out Noise, in electrons.
Gain=0.83                                                                 # Gain, in e/ADU
NSigma=6                                                                 # Sigma-times to reject outliers.
S=0.4                                                                    # Spacing between the polynomials.
N=3                                                                      # Number of coefficients in our polynomials.
P=PCoeff(FileName,b,Aperture,RON,Gain,NSigma,S,N,0,0,0)                  # We get our P-coefficient matrix...that easy!
                                                                         # (see the PCoeff function for more info). 
								         # In this example we take our spectrum from x=50 to
                                                                         # 1500, where x is in our matrix coordinates.
NCosmic=5.0                                                              # Sigma-times to reject cosmic rays.
Spectrum=getSpectrum(P,FileName,b,Aperture,RON,Gain,S,NCosmic,0,0)   # We obtain the spectrum and variances!
plot(Spectrum[0],Spectrum[1],'-')                                        # Plot the spectrum to check if everything went ok.
show()

