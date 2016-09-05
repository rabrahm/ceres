import matplotlib.pylab as plt
import os
import sys
print sys.path

import numpy
import scipy
import os
from scipy import interpolate
from scipy import integrate

def intercec(A,B):
	
	if A[0] < B[0]:
		ini = B[0]
	else:
		ini = A[0]

	if A[-1] < B[-1]:
		fin = A[-1]
	else:
		fin = B[-1]
	return ini,fin

			

def CCF(L1,F1,L2,F2,vi,vf):
	lux = 299792.458
	vel = vi
	delta = L1[1]-L1[0]
	CF = []
	vels = []
	while vel <=vf:
		
		L2p = L2*(1-vel/lux)
		ini,fin = intercec(L1,L2p) 
		
		I = numpy.where((L1 >= ini) & (L1 <= fin))[0]
		II = numpy.where((L2p >= ini) & (L2p <= fin))[0]
		if len(I)==0 or len(II)==0:	
			print 'Problem: no wavelenght intersection'
		
		wav = numpy.arange(ini,fin,delta)
		tck1 = interpolate.splrep(L1,F1,k=3,s=0)
		tck2 = interpolate.splrep(L2p,F2,k=3,s=0)
		F1s = interpolate.splev(wav,tck1,der=0)
		F2s = interpolate.splev(wav,tck2,der=0)

		CF.append(numpy.add.reduce(F1s*F2s)/numpy.sqrt(numpy.add.reduce(F1s*F1s)*numpy.add.reduce(F2s*F2s)))
		vels.append(vel)
		vel = vel + 1
	
	return vels,CF
