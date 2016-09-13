import pyfits
import scipy
import glob
import os
import sys

sys.path.append("../utils/GLOBALutils")
import GLOBALutils

from pylab import *

from rpy2 import robjects
import rpy2.robjects.numpy2ri
r = robjects.r

def FileClassify(dir, log,binning):
    """
    
    Classifies all files in a directory and writes a night log of science images

    """

    # define output lists
    thars      = []
    biases     = []
    quartzB    = []
    quartzR    = []
    science    = []
    obnames    = []
    exptimes   = []
    thar_dates = []

    f = open(log,'w')

    bad_files = []
    if os.access(dir+'bad_files.txt',os.F_OK):
    	bf = open(dir+'bad_files.txt')
	linesbf = bf.readlines()
	for line in linesbf:
		bad_files.append(dir+line[:-1])
	bf.close()
    
    all_files = glob.glob(dir+"/*fits")
    for archivo in all_files:
	dump = False
	for bf in bad_files:
		if archivo == bf:
			dump = True
			break
	h = pyfits.open(archivo)
	if dump == False and h[0].header['INSTRUME'] == 'echelle':
		if h[0].header['IMAGETYP'] == 'object':
		    science.append(archivo)
		    obname  = h[0].header['OBJNAME']
		    ra      = h[0].header['RA']
		    delta   = h[0].header['DEC']
		    airmass = h[0].header['AIRMASS']
		    texp    = h[0].header['EXPTIME']
		    date    = h[0].header['DATE-OBS'][:10]
		    hour    = h[0].header['DATE-OBS'][11:]
		    obnames.append( obname )
		    exptimes.append( texp ) 
		    line = "%-15s %10s %10s %8.2f %4.2f %8s %8s %s\n" % (obname, ra, delta, texp, airmass, date, hour, archivo)
		    f.write(line)
		elif h[0].header['IMAGETYP'] == 'zero':
		    biases.append(archivo)
		elif h[0].header['IMAGETYP'] == 'flat':
		    if h[0].header['FILTER'] == 'Blue':
		    	quartzB.append(archivo)
		    elif h[0].header['FILTER'] == 'Open':
			quartzR.append(archivo)
		elif h[0].header['IMAGETYP'] == 'comp':
		    thars.append(archivo)
		    mjd, mjd0 = mjd_fromheader(h)
		    thar_dates.append( mjd )
      
    f.close()

    return biases, quartzB, quartzR, science, thars, thar_dates, obnames, exptimes

def mjd_fromheader(h):
    """
    return modified Julian date from header
    """
    
    datetu = h[0].header['DATE-OBS'][:10]
    ut     = h[0].header['DATE-OBS'][11:]
    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[5:7]),int(datetu[8:]))
    ut = float(ut[:2]) + float(ut[3:5])/60. + float(ut[6:])/3600.
    mjd_start = mjd + ut/24.0
    secinday = 24*3600.0
    fraction = .5
    texp     = h[0].header['EXPTIME'] #sec
    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0

def bad_col_corr(dat):
	mat = np.array([[80  ,81  ,199 ,2048],\
			[240 ,241 ,530 ,2048],\
			[469 ,470 ,1555,1666],\
			[547 ,548 ,1690,1703],\
			[632 ,634 ,1985,1962],\
			[766 ,767 ,782 ,1980],\
			[832 ,833 ,1850,1906],\
			[960 ,961, 1590,1871],\
			[1263,1264,1772,1816],\
			[1334,1335,1454,2048],\
			[1361,1362,1885,1924],\
			[1395,1396,1906,1955],\
			[1460,1461,488 ,2048],\
			[1477,1478,758 ,2048],\
			[1510,1511,771 ,2048],\
			[1581,1582,1398,1763],\
			[1629,1630,420 ,2048],\
			[1661,1662,643 ,1980]])
	for bd in mat:
		new = .5*( dat[bd[2]:bd[3], bd[0]-1] + dat[bd[2]:bd[3],bd[1]] )
		bci = bd[0]
		while bci < bd[1]:
			dat[ bd[2]:bd[3], bci ] = new
			bci += 1 
	return dat

def MedianCombine(ImgList, bias = 0.):
    """
    Median combine a list of images
    """
    n = len(ImgList)
    if n == 0:
	print "\t\tWarning: 0 biases"
	return 0, 0, 0
    h = pyfits.open(ImgList[0])[0]
    d = h.data
    d = OverscanTrim(d)
    d = bad_col_corr(d)
    d -= bias
    
    factor = 1.25
    if (n < 3):
        factor = 1

    ronoise = factor * h.header['RDNOISE'] / np.sqrt(n)
    gain    = h.header['GAIN']

    
    if (n == 1):
        return d, ronoise, gain
    else:
        for i in range(n-1):
            h = pyfits.open(ImgList[i+1])[0]
	    #print ImgList[i+1], h.data.shape
	    d2 = OverscanTrim(h.data)
	    d2 = bad_col_corr(d2)
            d = np.dstack( (d ,d2 - bias ))
        return np.median(d,axis=2), ronoise, gain

def OverscanTrim(d):
    """
    Overscan correct and Trim a refurbished APO image
    """
    ov = d[:20,21:-59]
    ov = np.median(ov,axis = 0)
    newdata = d[20:,21:-59]
    overscan = np.zeros(newdata.shape)
    for i in range(overscan.shape[1]):
	overscan[:,i] = ov

    newdata = newdata - overscan
    return newdata

