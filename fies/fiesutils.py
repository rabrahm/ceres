import pyfits
import numpy as np
import scipy
import copy
import glob
import os
import matplotlib.pyplot as plt
import sys
from pylab import *

base = "../"
sys.path.append(base+"utils/GLOBALutils")
import GLOBALutils

from rpy2 import robjects
import rpy2.robjects.numpy2ri
r = robjects.r

def ra_from_sec(ra,time=True):
	ra = float(ra)
	sign = ' '
	if ra < 0:
		sign = '-' 
		ra *= -1

	hh = ra/3600.
	mm = (hh - int(hh))*60.
	ss = (mm - int(mm))*60.
	shh = str(int(hh))
	smm = str(int(mm))
	sss = str(np.around(ss,2))
	if hh<10:
		shh = '0' + shh
	if mm<10:
		smm = '0' + smm
	if ss<10:
		sss = '0' + sss
	return sign + shh + ':' + smm + ':' + sss

def FileClassify(diri, log,binning=1):
    """
    
    Classifies all files in a directory and writes a night log of science images

    """

    # define output lists
    sim_sci        = []
    biases         = []
    flats          = []
    ThAr_ref       = []
    ThAr_ref_dates = []
    flat_ref_dates = []
    bias_ref_dates = []
    obnames        = []
    exptimes       = []

    f = open(log,'w')

    #Do not consider the images specified in dir+badfiles.txt
    bad_files = []
    if os.access(diri+'bad_files.txt',os.F_OK):
    	bf = open(diri+'bad_files.txt')
	linesbf = bf.readlines()
	for line in linesbf:
		bad_files.append(diri+line[:-1])
	bf.close()
    
    all_files = glob.glob(diri+"/*fits")
    for archivo in all_files:
	#print archivo
	dump = False
	for bf in bad_files:
		if archivo == bf:
			dump = True
			break
	if dump == False:
		h = pyfits.open(archivo)
		hd = pyfits.getheader(archivo)
		if int(h[0].header['DETXBIN']) == binning and int(h[0].header['DETYBIN']) == binning:
			print archivo, h[0].header['IMAGETYP'], h[0].header['SHSTAT'], h[0].header['EXPTIME'], h[0].header['OBJECT'], h[0].header['TCSTGT'] 
			if h[0].header['IMAGETYP'].replace(' ','') == '' or h[0].header['IMAGETYP'] == 'science':
			    sim_sci.append(archivo)
			    obname = h[0].header['OBJECT']
			    obnames.append( obname )
			    ra     = ra_from_sec(h[0].header['RA']*3600.*24./360.)
			    delta  = ra_from_sec(h[0].header['DEC']*3600.)
			    airmass= float(h[0].header['AIRMASS'])
			    texp   = float(h[0].header['EXPTIME'])

			    date   =  h[0].header['DATE-OBS']		    
			    hour   = date[11:]
			    date    = date[:10]
			    exptimes.append( texp ) 
			    line = "%-15s %10s %10s %8.2f %4.2f %8s %11s %s\n" % (obname, ra, delta, texp, airmass, date, hour, archivo)
			    f.write(line)
			elif h[0].header['IMAGETYP'] == 'BIAS':
			    biases.append(archivo)
			    mjd, mjd0 = mjd_fromheader2(h)
			    bias_ref_dates.append( mjd )
			elif h[0].header['IMAGETYP'] == 'FLAT':
			    flats.append(archivo)
			    mjd, mjd0 = mjd_fromheader2(h)
			    flat_ref_dates.append( mjd )
			elif h[0].header['IMAGETYP'] == 'WAVE':
			    ThAr_ref.append(archivo)
			    mjd, mjd0 = mjd_fromheader2(h)
			    ThAr_ref_dates.append( mjd )
		
    
    flat_ref_dates = np.array(flat_ref_dates)
    flats = np.array(flats)
    IS = np.argsort(flat_ref_dates)
    flat_ref_dates = flat_ref_dates[IS]
    flats = flats[IS]
    #for i in range(len(flats)):
    #	print 'flat',flats[i], flat_ref_dates[i]

    bias_ref_dates = np.array(bias_ref_dates)
    biases = np.array(biases)
    IS = np.argsort(bias_ref_dates)
    bias_ref_dates = bias_ref_dates[IS]
    biases = biases[IS]
    #for i in range(len(biases)):
    #	print 'bias',biases[i], bias_ref_dates[i]
    f.close()

    return biases, np.array(flats), np.array(ThAr_ref), sim_sci, np.array(ThAr_ref_dates), obnames, exptimes

def mjd_fromheader2(h):
	"""
	return modified Julian date from header
	"""

	datetu = h[0].header['DATE-OBS']
	mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[:4]),int(datetu[5:7]),int(datetu[8:10]))
	ho = int(datetu[11:13])
	mi = int(datetu[14:16])
	se = float(datetu[17:])
	ut = float(ho) + float(mi)/60.0 + float(se)/3600.0
	mjd_start = mjd + ut/24.0

	secinday = 24*3600.0
	fraction = 0.5
	texp     = h[0].header['EXPTIME'] #sec

	mjd = mjd_start + (fraction * texp) / secinday

	return mjd, mjd0

def get_RONGAIN(hd):
	return 2.8, 0.737

def MedianCombine(ImgList, zero='none', binning=1):
	"""
	Median combine a list of images
	"""

	n = len(ImgList)
	if n==0:
		raise ValueError("empty list provided!")

	h = pyfits.open(ImgList[0])

	d1 = h[1].data
	h1 = h[1].header
	d1 = OverscanTrim(d1, binning=binning)
	if zero != 'none':
		z = pyfits.open(zero)[0]
		d1 -= z.data
    
	factor = 1.25
	if (n < 3):
		factor = 1

	ron1,gain1 = get_RONGAIN(h[0].header)

	ron1 = factor * ron1 / np.sqrt(n)
	if n>1:
		for i in range(n-1):
			td = pyfits.open(ImgList[i+1])
			if zero == 'none':
				d1 = np.dstack((d1,OverscanTrim(td[1].data, binning=binning)))
		
			else:
				d1 = np.dstack((d1,OverscanTrim(td[1].data, binning=binning)-z.data))
		d1 = np.median(d1,axis=2)
	return d1, ron1, gain1

def OverscanTrim(dat,binning=1):
	"""
	Overscan correct and Trim a refurbished FEROS image
	"""
	ff = 2098
	ii = 50
	ff = int(np.around(ff/binning))
	ii = int(np.around(ii/binning))
	os = dat[:,ff:]
	s = np.median(os)
	newdata = dat[:,ii:ff].copy() - s

	return newdata