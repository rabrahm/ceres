from __future__ import print_function
#!/usr/bin/python

import numpy
import math
import correlation
import os

path = '/home/rabrahm/Desktop/spec2/'
spectra = os.listdir(path)

cant = len(spectra)

m=0
Tv = numpy.zeros(cant,float)
Gv = numpy.zeros(cant,float)
Zv = numpy.zeros(cant,float)
Vv = numpy.zeros(cant,float)
"""print spectra"""
m=5
print(spectra)
while m<cant:

    Tv[m],Gv[m],Zv[m],Vv[m] = correlation.CCF(path,spectra[m])

    m=m+1
m=0
while m<cant:

    print(spectra[m],Tv[m],Gv[m],Zv[m],Vv[m])

    m=m+1
