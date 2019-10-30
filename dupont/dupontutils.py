import sys
import matplotlib
matplotlib.use("Agg") 

base = '../'
sys.path.append(base+"utils/GLOBALutils")
import GLOBALutils

import numpy as np
import scipy
from astropy.io import fits as pyfits
import os
import glob
import scipy.signal
from scipy.signal import medfilt
from scipy import interpolate
import copy


from pylab import *



def milk_comb(ImgList, darks, zero='Bias.fits'):
	n = len(ImgList)
	if n==0:
		raise ValueError("empty list provided!")
	h = pyfits.open(ImgList[0])[0]
	d = h.data

	head = pyfits.getheader(ImgList[0])
	expt = head['EXPTIME']
	
	d = OverscanTrim(d,h.header['BIASSEC'])

	Master = pyfits.getdata(zero)

	if len(darks) > 0:
		Dark = get_dark(darks,expt)
	else:
		Dark = np.zeros((d.shape[0],d.shape[1]),float)

	d -= Master
	d -= Dark

	factor = 1.25
	if (n < 3):
		factor = 1
	
	ron1 = h.header['ENOISE']
	gain = h.header['EGAIN']
	

	ronoise = factor * h.header['ENOISE'] / np.sqrt(n)

	if (n == 1):
		return d, ron1, gain
	
	else:
		for i in range(n-1):
			h = pyfits.open(ImgList[i+1])[0]
			head = pyfits.getheader(ImgList[i+1])
			expt = head['EXPTIME']
			if len(darks) > 0:
				Dark = get_dark(darks,expt)
			else:
				Dark = np.zeros((d.shape[0],d.shape[1]),float)
			rd = OverscanTrim(h.data,h.header['BIASSEC']) - Master - Dark 
			d = np.dstack((d,rd/np.median(rd)))
		out = np.median(d,axis=2)
		return out, ronoise, gain
	

def FileClassify(path,log):
	biases   = []
	milkflat = []
	objects  = []
	darks    = []
	thars    = []
	skflat   = []
	qflat    = []
	f = open(log,'w')
	archs = glob.glob(path+'ccd*.fits')
	if os.access(path+'bad_files.txt',os.F_OK):
		ff = open(path + 'bad_files.txt', 'r')
		bfiles = ff.readlines()
	else:
		bfiles = []
	for arch in archs:
		use = True
		for bf in bfiles:
			if arch == path + bf[:-1]:
				print 'Dumped file', arch
				use = False
				break
		if use:	
			h = pyfits.open(arch)
			header = pyfits.getheader(arch)
			name = header['OBJECT']

			if True:
				if header['EXPTYPE'] == 'Object':
			
					if name.count('sky')>0:
						skflat.append(arch)
					elif (name.lower()).count('milky')>0:
						milkflat.append(arch)
					elif name.count('qFlat')>0:
						qflat.append(arch)
					else:
						expt = header['EXPTIME']
						ra   = header['RA']
						dec  = header['DEC']
						ams  = header['AIRMASS']
						date = header['DATE-OBS']
						UT   = header['UT-TIME']
						line = "%-15s %10s %10s %8.2f %4.2f %8s %8s %s\n" % (name, ra, dec, expt, ams, date, UT, arch) 
						f.write(line)
						cos = name.split('ThAr')
						if 'thar' in name.lower() or len(cos)> 1 or 'comp' in name.lower():
							thars.append(arch)
						elif name.count('milky')>0 or name.count('Milky')>0 or name.count('flat')>0:
							milkflat.append(arch)
						else:
							objects.append(arch)	
				elif  header['EXPTYPE'] == 'Bias':
					biases.append(arch)
				elif  header['EXPTYPE'] == 'Flat':
					nam = header['OBJECT']
					if nam.count('sky')>0 or nam.count('Sky')>0 :
						skflat.append(arch)
					elif nam.count('milky')>0 or nam.count('Milky')>0 or nam.count('flat')>0:
						milkflat.append(arch)
				elif  header['EXPTYPE'] == 'Dark':
					darks.append(arch)
	
			h.close()
	f.close()
	return biases, milkflat, skflat, objects, thars, darks

def MedianCombine(ImgList, zero_bo=False, zero='Bias.fits', dark_bo=False, darks=[], flat_bo=False, flat='Flat.fits'):
	"""
	Median combine a list of images
	"""
	n = len(ImgList)
	if n==0:
		raise ValueError("empty list provided!")

	h = pyfits.open(ImgList[0])[0]
	d = h.data
	d = OverscanTrim(d,h.header['BIASSEC'])

	if zero_bo:
		Master = pyfits.getdata(zero)
	else:
		Master = np.zeros((d.shape[0],d.shape[1]),float)

	if dark_bo and len(darks)!=0:
		hd = pyfits.getheader(ImgList[0])
		time = hd['EXPTIME']
		Dark = get_dark(darks, time)
	else:
		Dark = np.zeros((d.shape[0],d.shape[1]),float)

	if flat_bo:
		Flat = pyfits.getdata(flat)
	else:
		Flat = np.zeros((d.shape[0],d.shape[1]),float) + 1.0

	if flat_bo:
		d = (d - Master - Dark)/Flat
	else:
		d = (d - Master - Dark)

	factor = 1.25

	if (n < 3):
		factor = 1

	ronoise = factor * h.header['ENOISE'] / np.sqrt(n)
	gain = h.header['EGAIN']

	if (n == 1):
		return d, ronoise, gain

	else:
		for i in range(n-1):
			h = pyfits.open(ImgList[i+1])[0]
			if flat_bo:
				d = np.dstack((d,(OverscanTrim(h.data,h.header['BIASSEC'])-Master-Dark)/Flat))
			else:
				
				d = np.dstack((d,OverscanTrim(h.data,h.header['BIASSEC'])-Master-Dark))
	
	return np.median(d,axis=2), ronoise, gain


def OverscanTrim(d,bsec):
	"""
	Overscan correct and Trim a refurbished DuPont image
	"""
	bsec = bsec[1:-1]
	bsec1 = bsec.split(',')[0]
	bsec2 = bsec.split(',')[1]
	b11 = int(bsec1.split(':')[0])
	b12 = int(bsec1.split(':')[1])
	b21 = int(bsec2.split(':')[0])
	b22 = int(bsec2.split(':')[1])

	t1 = d[:1500,:b11-1]
	t2 = d[b21-1:,:b11-1]
	
	nd = np.zeros((t1.shape[0],t1.shape[1]),float)
	nd[:1500,:] = t1

	overscan1 = np.median(t2,axis=0)
	newdata = nd - overscan1

	return newdata

def get_dark(darks,t):
	exact = 0
	dts = []
	for dark in darks:
		hd = pyfits.getheader(dark)
		dt = hd['EXPTIME']
		dts.append(dt)
		if dt == t:
			DARK = pyfits.getdata(dark)
			exact = 1
	
	dts = np.array(dts)
	if exact == 0:
		if t < dts.min():
			I = np.where( dts == dts.min() )[0]
			DARK = pyfits.getdata(darks[I[0]])*t/dts[I[0]]
		elif t > dts.max():
			I = np.where( dts == dts.max() )[0]
			DARK = pyfits.getdata(darks[I[0]])*t/dts[I[0]]
		else:
			tmin = dts.min()
			tmax = dts.max()
			I = np.where( dts == dts.min() )[0]
			Dmin = pyfits.getdata(darks[I[0]])
			Dminname=darks[I[0]]
			I = np.where( dts == dts.max() )[0]
			Dmax = pyfits.getdata(darks[I[0]])
			Dmaxname = darks[I[0]]
			
			i = 0
			while i < len(dts):
				if dts[i] < t and dts[i] > tmin:
					tmin = dts[i]
					Dminname = darks[i]
					Dmin = pyfits.getdata(darks[i])
				elif dts[i] > t and dts[i] < tmax:
					tmax = dts[i]
					Dmaxname = darks[i]
					Dmax = pyfits.getdata(darks[i])
				i+=1
			
			num = Dmax - Dmin
			den = tmax-tmin
			m = num/den
			n = Dmax - m*tmax
			DARK = m*t+n
					
	return DARK

def get_blaze(LL,FF, low=1.0, hi=3.0, n = 6):
	NF = FF.copy()
	for j in range(LL.shape[0]):
		L = LL[j]
		F = FF[j]
		ejex = np.arange(len(F))
		F[:150]  = 0.0
		F[-150:] = 0.0
		Z = np.where(F!=0)[0]
		F = scipy.signal.medfilt(F[Z],31)
		ejexx = ejex.copy()
		ejex = ejex[Z]
		L = L[Z]
		I = np.where((L>5870) & (L<5890))[0]
		if len(I)>0:
			W = np.where(L<5870)[0]
			R = np.where(L>5890)[0]
			ejetemp = np.hstack((ejex[W],ejex[R]))
			Ftemp = np.hstack((F[W],F[R]))	
			coefs = np.polyfit(ejetemp,Ftemp,n)
			fit = np.polyval(coefs,ejetemp)
			
		else:
			ejetemp=ejex
			Ftemp=F
			coefs = np.polyfit(ejex,F,n)
			fit = np.polyval(coefs,ejex)
		i = 0
		while i < 30:
			res = Ftemp - fit
			IP = np.where((res>=0) & (Ftemp!=0.0))[0]
			IN = np.where((res<0) & (Ftemp!=0.0))[0]
			devp = np.mean(res[IP])
			devn = np.mean(res[IN])
			I = np.where((res > -low*abs(devn)) & (res < hi*abs(devp)) & (Ftemp!=0))[0]
			coefs = np.polyfit(ejetemp[I],Ftemp[I],n)
			fit = np.polyval(coefs,ejetemp)
			i+=1
		fit = np.polyval(coefs,ejexx)

		NF[j]=fit
	NNF = NF.copy()
	for j in range(LL.shape[0]):

		L = LL[j]
		I = np.where((L>6520) & (L<6600))[0]
		if len(I)>0:
			
			if j+2 < LL.shape[0]:
				for i in range(len(L)):
					vec = np.array([NF[j-2,i],NF[j-1,i],NF[j+1,i],NF[j+2,i]])
					tck = scipy.interpolate.splrep(np.array([0.0,1.0,3.0,4.0]),vec,k=2)
					NNF[j,i] = scipy.interpolate.splev(2.0,tck,der=0)
			elif j+1 < LL.shape[0]:
				for i in range(len(L)):
					vec = np.array([NF[j-2,i],NF[j-1,i],NF[j+1,i]])
					tck = scipy.interpolate.splrep(np.array([0.0,1.0,3.0]),vec,k=1)
					NNF[j,i] = scipy.interpolate.splev(2.0,tck,der=0)
			elif j < LL.shape[0]:
				for i in range(len(L)):
					vec = np.array([NF[j-3,i],NF[j-2,i],NF[j-1,i]])
					tck = scipy.interpolate.splrep(np.array([0.0,1.0,2.0]),vec,k=1)
					NNF[j,i] = scipy.interpolate.splev(3.0,tck,der=0)
		
			
		I = np.where((L>4870) & (L<4880))[0]
		if len(I)>0:
			
			if j+2 < LL.shape[0]:
				for i in range(len(L)):
					vec = np.array([NF[j-2,i],NF[j-1,i],NF[j+1,i],NF[j+2,i]])
					tck = scipy.interpolate.splrep(np.array([0.0,1.0,3.0,4.0]),vec,k=2)
					NNF[j,i] = scipy.interpolate.splev(2.0,tck,der=0)
			elif j+1 < LL.shape[0]:
				for i in range(len(L)):
					vec = np.array([NF[j-2,i],NF[j-1,i],NF[j+1,i]])
					tck = scipy.interpolate.splrep(np.array([0.0,1.0,3.0]),vec,k=1)
					NNF[j,i] = scipy.interpolate.splev(2.0,tck,der=0)
			else:
				for i in range(len(L)):
					vec = np.array([NF[j-3,i],NF[j-2,i],NF[j-1,i]])
					tck = scipy.interpolate.splrep(np.array([0.0,1.0,2.0]),vec,k=1)
					NNF[j,i] = scipy.interpolate.splev(3.0,tck,der=0)

		I = np.where((L>4320) & (L<4325))[0]
		if len(I)>0:
			if j+2 < LL.shape[0]:
				for i in range(len(L)):
					vec = np.array([NF[j-2,i],NF[j-1,i],NF[j+1,i],NF[j+2,i]])
					tck = scipy.interpolate.splrep(np.array([0.0,1.0,3.0,4.0]),vec,k=2)
					NNF[j,i] = scipy.interpolate.splev(2.0,tck,der=0)
	return NNF

def get_close(tht,rat,dect,fits):
	t0 = 1000000.
	close = fits[0]
	for fit in fits:
		#print close
		hd = pyfits.getheader(fit)
		sct,mjd0 = mjd_fromheader(hd)
		expt = hd['EXPTIME']/(3600.*24.)
		dec = hd['DEC-D']
		ra = hd['RA-D']
		if abs(dec - dect)<0.05 and abs(ra - rat)<0.05:
			#print sct+expt,tht
			if abs(sct+expt-tht) < t0:
				t0 = abs(sct+expt-tht)
				close = fit
	return close

			
def b_col(d):
	d[:,746] = 0.5*(d[:,745]+d[:,748])
	d[:,747] = 0.5*(d[:,745]+d[:,748])		
	return d

def mjd_fromheader(h):
    """
    return modified Julian date from header
    """  
    datetu = h['UT-DATE']
    timetu = h['UT-TIME']
    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[:4]),int(datetu[5:7]),int(datetu[8:]))
    ho = int(timetu[:2])
    mi = int(timetu[3:5])
    se = float(timetu[7:])
    ut = float(ho) + float(mi)/60.0 + float(se)/3600.0
    mjd_start = mjd + ut/24.0

    secinday = 24*3600.0
    fraction = 0.5
    texp     = h['EXPTIME'] #sec

    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0
