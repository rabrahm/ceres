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


def getObName(h):
    """
    Get name of object under consideration
    """
    obname = h[0].header['OBJECT'].upper().replace(' ','')
    return obname

def ra_from_sec(ra):
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
	
def yesterday(date):
	y,m,d = int(date[:4]),int(date[5:7]),int(date[8:])
	pd = d - 1
	pm = m
	py = y
	if pd == 0:
		pm = m - 1
		py = y
		if pm in [1,3,5,7,8,10,12]:
			pd = 31
		elif pm in [4,6,9,11]:
			pd = 30
		else:
			pd = 28
			p,q,r = False,False,False
			p = y%4
			q = y%100
			r = y%400
			if p == 0:
				p = True
			if q == 0:
				q = True
			if r == 0:
				r = True
			if p and (not q or r):
 				pd = 29
		
		if pm == 0:
			py = y - 1
			pm = 12
			pd = 31
		 	
	sy,sm,sd = str(py),str(pm),str(pd)
	if pm < 10:
		sm = '0'+sm
	if pd < 10:
		sd = '0'+sd
	return sy + '-' + sm + '-' + sd

def FileClassify(diri, log):
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
		if h[0].header['IMAGETYP'] == 'Science':
		    sim_sci.append(archivo)
		    obname = getObName(h)
		    obnames.append( obname )
		    ra     = ra_from_sec(h[0].header['RA'])
		    delta  = ra_from_sec(h[0].header['DEC'])
		    airmass= float(h[0].header['AIRMASS'])
		    texp   = float(h[0].header['EXPTIME'])

		    date   =  h[0].header['DATE']		    
		    stime   = h[0].header['UT_START']
		    etime   = h[0].header['UT_END']
		    date    = date[:10]
		    if stime > etime:
			date = yesterday(date)
		    hour   = ra_from_sec(stime)
		    exptimes.append( texp ) 
		    line = "%-15s %10s %10s %8.2f %4.2f %8s %11s %s\n" % (obname, ra, delta, texp, airmass, date, hour, archivo)
		    f.write(line)
		elif h[0].header['IMAGETYP'] == 'Bias':
		    biases.append(archivo)
		    mjd, mjd0 = mjd_fromheader2(h)
		    bias_ref_dates.append( mjd )
		elif h[0].header['IMAGETYP'] == 'Flat':
		    flats.append(archivo)
		    mjd, mjd0 = mjd_fromheader2(h)
		    flat_ref_dates.append( mjd )
		elif h[0].header['IMAGETYP'] == 'Calibration':
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

def mjd_from_data(datetu,ut,texp,fraction):
    secinday = 24*3600.0
    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[5:7]),int(datetu[8:10]))

    ut        = (float(ut[:2])*3600. + float(ut[3:5])*60. + float(ut[6:]))
    mjd_start = mjd + ut / secinday

    mjd = mjd_start + (fraction * texp) / secinday
    return mjd, mjd0

def mjd_fromheader(h):
    """
    return modified Julian date from header
    """
    secinday = 24*3600.0

    date   = h[0].header['DATE']		    
    stime   = h[0].header['UT_START']
    etime   = h[0].header['UT_END']
    date    = date[:10]

    if stime > etime:
	date = yesterday(date)
    hour   = ra_from_sec(stime)

    datetu = date.replace('-',':')
    ut = hour.replace('-',':')
    datetu = date.replace(' ','')
    ut = hour.replace(' ','')
    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[5:7]),int(datetu[8:10]))
    ut        = (float(ut[:2])*3600. + float(ut[3:5])*60. + float(ut[6:]))
    mjd_start = mjd + ut / secinday

def mjd_fromheader2(h):
    """
    return modified Julian date from header
    """
    rtime = 90.

    secinday = 24*3600.0

    date    = h[0].header['DATE']
    hour   = date[11:]		    
    date    = date[:10]
    

    datetu = date.replace('-',':')
    ut = hour.replace('-',':')
    datetu = date.replace(' ','')
    ut = hour.replace(' ','')
    #print datetu
    #print ut
    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[5:7]),int(datetu[8:10]))
    ut        = (float(ut[:2])*3600. + float(ut[3:5])*60. + float(ut[6:]))
    mjd_start = mjd + ut / secinday

    fraction = 0.5
    texp     = float(h[0].header['EXPTIME']) #sec

    mjd = mjd_start - (fraction * texp + rtime) / secinday

    return mjd, mjd0

def MedianCombine(ImgList, zero='none'):
    """

    Median combine a list of images

    """

    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    h = pyfits.open(ImgList[0])

    d1 = h[0].data
    h1 = h[0].header
    d1 = b_col(d1)
    d1 = OverscanTrim(d1).astype('float')

    if zero != 'none':
	z = pyfits.open(zero)[0]
	d1 -= z.data.astype('float')
    
    factor = 1.25
    if (n < 3):
        factor = 1

    ron1,gain1 = float(h1['CCDRON']),float(h1['CCDGAIN'])

    ron1 = factor * ron1 / np.sqrt(n)
    if n>1:
	    for i in range(n-1):
		td = pyfits.open(ImgList[i+1])
		if zero == 'none':
		    d1 = np.dstack((d1,OverscanTrim(b_col(td[0].data))))
		
		else:
		    d1 = np.dstack((d1,OverscanTrim(b_col(td[0].data))-z.data))
	    d1 = np.median(d1,axis=2)

    return d1, ron1, gain1

def OverscanTrim(d):
    """
    Overscan correct and Trim a refurbished CORALIE image
    """
    """
    # bias has no significant structure, so a single median suffices, I think
    # overscan = [0:49] [2097:2145]
    data = d[:,51:2099]
    ov1  = d[:,:50]
    ov2  = d[:,2098:]
    ov = np.vstack((ov1.transpose(),ov2.transpose()))
    overscan = median(ov,axis=0)
    overscan = np.array([overscan]*2048)
    overscan = overscan.transpose()
    newdata = data - overscan
    """
    newdata = d.copy()
    return newdata

def b_col(d):
	d[:,233] = 0.5*(d[:,232]+d[:,235])
	d[:,234] = 0.5*(d[:,232]+d[:,235])	
	return d
