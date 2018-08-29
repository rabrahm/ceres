import matplotlib
matplotlib.use("Agg")
import numpy as np
import scipy
from scipy import signal,interpolate
import scipy.special as sp
import copy
import glob
import os
import sys
from pylab import *
sys.path.append("../utils/GLOBALutils")
import GLOBALutils
from astropy.io import fits as pyfits

import statsmodels.api as sm
lowess = sm.nonparametric.lowess

def FileClassify(diri, log):
    """
    
    Classifies all files in a directory and writes a night log of science images

    """

    # define output lists
    science        = []
    flats          = []
    ThAr_ref       = []
    ThAr_ref_dates = []
    flat_ref_dates = []
    obnames        = []
    exptimes       = []
    orderref       = []
    biases         = []

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
		print archivo,  h[0].header['HIERARCH ESO DPR CATG'],h[0].header['HIERARCH ESO DPR TYPE']
		if 'SCIENCE' in h[0].header['HIERARCH ESO DPR CATG']:
		    science.append(archivo)
		    obname = h[0].header['OBJECT']
		    obnames.append( obname )
		    ra     = h[0].header['RA']
		    delta  = h[0].header['DEC']
		    airmass= h[0].header['HIERARCH ESO TEL AIRM START']
		    texp   = h[0].header['EXPTIME']
		    date   = h[0].header['DATE-OBS'].split('T')[0]	    
		    hour   = h[0].header['DATE-OBS'].split('T')[1]
		    exptimes.append( texp )
		    #rint obname, ra, delta, texp, airmass, date, hour, archivo
		    line = "%-15s %10s %10s %8.2f %4.2f %8s %11s %s\n" % (obname, ra, delta, float(texp), float(airmass), date, hour, archivo)
		    f.write(line)
		elif 'CALIB' in h[0].header['HIERARCH ESO DPR CATG']:
			if h[0].header['HIERARCH ESO DPR TYPE'] == 'LAMP,FLAT':
				flats.append(archivo)
				mjd, mjd0 = mjd_fromheader2(h[0].header)
				flat_ref_dates.append( mjd )
			elif h[0].header['HIERARCH ESO DPR TYPE'] == 'LAMP,WAVE':
				ThAr_ref.append(archivo)
				mjd, mjd0 = mjd_fromheader2(h[0].header)
				ThAr_ref_dates.append( mjd )
			elif h[0].header['HIERARCH ESO DPR TYPE'] == 'LAMP,ORDERDEF':
				orderref.append(archivo)
			elif h[0].header['HIERARCH ESO DPR TYPE'] == 'BIAS':
				biases.append(archivo)

    flat_ref_dates,ThAr_ref,ThAr_ref_dates = np.array(flat_ref_dates),np.array(ThAr_ref),np.array(ThAr_ref_dates)
    flats = np.array(flats)
    IS = np.argsort(flat_ref_dates)
    flat_ref_dates = flat_ref_dates[IS]
    flats = flats[IS]
    IS = np.argsort(ThAr_ref_dates)
    ThAr_ref_dates = ThAr_ref_dates[IS]
    ThAr_ref = ThAr_ref[IS]


    return np.array(biases),np.array(flats), np.array(orderref), np.array(ThAr_ref), np.array(science), np.array(ThAr_ref_dates), obnames, exptimes

def mjd_fromheader2(hd):
    """
    return modified Julian date from header
    """

    secinday = 24*3600.0

    date   = hd['DATE-OBS'].split('T')[0]	    
    ut     = hd['DATE-OBS'].split('T')[1]	    
    
    datetu = date.replace('-',':')

    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[5:7]),int(datetu[8:10]))
    ut        = (float(ut[:2])*3600. + float(ut[3:5])*60. + float(ut[6:]))
    mjd_start = mjd + ut / secinday

    fraction = 0.5
    texp     = float(hd['EXPTIME']) #sec

    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0

def MedianCombine(ImgList, bias=0.):
    """
    Median combine a list of images
    """

    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    print 0, ImgList[0]
    h = pyfits.open(ImgList[0])

    d1 = h[1].data
    d2 = h[2].data
    
    factor = 1.25
    if (n < 3):
        factor = 1

    ron1,gain1 = float(h[1].header['HIERARCH ESO DET OUT1 RON']),float(h[1].header['HIERARCH ESO DET OUT1 GAIN'])
    ron2,gain2 = float(h[2].header['HIERARCH ESO DET OUT1 RON']),float(h[2].header['HIERARCH ESO DET OUT1 GAIN'])

    ron1 = factor * ron1 / np.sqrt(n)
    ron2 = factor * ron2 / np.sqrt(n)
    if n>1:
    	for i in range(n-1):
		print i,ImgList[i+1]
    		td = pyfits.open(ImgList[i+1])
    		td = np.dstack((td[1].data,td[2].data)) - bias
    		d1 = np.dstack((d1,td[:,:,0]))
    		d2 = np.dstack((d2,td[:,:,1]))
    	d1 = np.median(d1,axis=2)
    	d2 = np.median(d2,axis=2)

	combined = np.dstack((d1,d2))
    return combined, [ron1,ron2], [gain1,gain2]

def OverscanTrim(d,ds):
    """
    Overscan correct and Trim a refurbished HIRES image
    """
    ds = ds.split(',')
    dsy = ds[0][1:]
    dsx = ds[1][:-1]
    dsy = dsy.split(':')
    dsx = dsx.split(':')
    x0 = int(dsx[0])-1
    x1 = int(dsx[1])
    y0 = int(dsy[0])-1
    y1 = int(dsy[1])
    newdata = d[x0:x1]
    newdata = newdata[:,y0:y1]
    newdata = newdata[:4000,:]
    return newdata

def bac_flat(flat,traces,ext):
	ext = int(np.around(.5*ext))
	ejex = np.arange(flat.shape[0])
	peaks = np.zeros((len(traces),flat.shape[0]))
	bac = np.zeros(flat.shape)
	for i in range(len(traces)):
		peaks[i] = np.around(np.polyval(traces[i],ejex))
	peaks = peaks.astype('int')
	for px in np.arange(flat.shape[0]):
		cut = flat[px]
		cut = scipy.signal.medfilt(cut,11)
		refx,refy = [],[]
		for i in np.arange(len(traces)-1)+1:
			"""
			if i == 0:
				ncut = cut[peaks[i,px]-40:peaks[i,px]]
				im = peaks[i,px]-40 + np.argmin(ncut)
				refx.append(im)
				refy.append(cut[im])
			"""
			ncut = cut[peaks[i-1,px]:peaks[i,px]]
			im = peaks[i-1,px] + np.argmin(ncut)
			refx.append(im)
			refy.append(cut[im])
		refx,refy = np.array(refx),np.array(refy)
		tck = scipy.interpolate.splrep(refx,refy,k=1)
		bac[px] = scipy.interpolate.splev(np.arange(flat.shape[1]),tck)
		bac[px,:refx[0]] = refy[0]
		bac[px,refx[-1]:] = refy[-1]
	bac = scipy.signal.medfilt(bac,[11,1])
	return bac
			
def ccd_flat(flat,traces,ext):
	print ext
	nflat = np.ones(flat.shape)
	ext = int(np.around(ext))
	ejex = np.arange(flat.shape[0])
	peaks = np.zeros((len(traces),flat.shape[0]))
	for i in range(len(traces)):
		peaks[i] = np.around(np.polyval(traces[i],ejex))
	peaks = peaks.astype('int')
	for px in np.arange(flat.shape[0]):
		for py in range(len(traces)):
			tf = flat[px,peaks[py,px]-ext:peaks[py,px]+ext+1]
			nflat[px,peaks[py,px]-ext:peaks[py,px]+ext+1] = tf

	return nflat
			

def normalize_flat(flat,traces,ext):
	nflat = np.ones(flat.shape)
	ext = int(np.around(.5*ext))+2
	ejex = np.arange(flat.shape[0])
	peaks = np.zeros((len(traces),flat.shape[0]))
	for i in range(len(traces)):
		peaks[i] = np.around(np.polyval(traces[i],ejex))
	peaks = peaks.astype('int')
	spec_flat = np.zeros((len(traces),flat.shape[0]))

	for px in np.arange(flat.shape[0]):
		for py in range(len(traces)):
			tf = flat[px,peaks[py,px]-ext:peaks[py,px]+ext+1]
			tf2 = flat[px,peaks[py,px]-int(np.around(ext/2.)):peaks[py,px]+int(np.around(ext/2.))+1]
			nflat[px,peaks[py,px]-ext:peaks[py,px]+ext+1] = tf / np.median(tf)
			spec_flat[py,px] = np.sum(tf2)
	#for i in range(len(traces)):
	#	plot(spec_flat[i])
	#	show()
	for i in range(len(traces)):
		I = np.where(np.isnan(spec_flat[i]))[0]
		spec_flat[i][I] = 0.
		ss =  scipy.signal.medfilt(spec_flat[i],15)
		spec_flat[i] /= ss.max()

	#spec_flat_o = spec_flat.copy()
	#for i in range(len(traces)):
	#	spec_flat[i] = spec_flat[i]/scipy.signal.medfilt(spec_flat[i],301)
	#	spec_flat[i] = spec_flat[i]/sm.max()
	return nflat,spec_flat


def get_sky(sc,lim,span1,span2, typ='median'):
	sky = np.zeros(sc.shape)
	ejeX = np.arange(sc.shape[0])

	for y in range(sc.shape[1]):
		lims = np.around(lim[:,y]).astype('int')
		for j in range(len(lims)):
			med1 = np.median(sc[lims[j]-span2:lims[j]-span1+1,y])
			med2 = np.median(sc[lims[j]+span1:lims[j]+span2+1,y])
			med = np.mean([med1,med2])
			sky[lims[j]-span2:lims[j]+span2+1,y] = med

	return sky