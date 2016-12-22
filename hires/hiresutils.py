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
import pyfits

import statsmodels.api as sm
lowess = sm.nonparametric.lowess

def FileClassify(diri, log):
    """
    
    Classifies all files in a directory and writes a night log of science images

    """

    # define output lists
    sim_sci        = []
    flats          = []
    ThAr_ref       = []
    ThAr_ref_dates = []
    flat_ref_dates = []
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
		if h[0].header['IMAGETYP'] == 'object':
		    sim_sci.append(archivo)
		    obname = hd['TARGNAME']
		    obnames.append( obname )
		    ra     = hd['RA']
		    delta  = hd['DEC']
		    airmass= hd['AIRMASS']
		    texp   = hd['EXPTIME']
		    date   = hd['DATE-OBS']		    
		    hour   = hd['UTC']
		    exptimes.append( texp )
		    #rint obname, ra, delta, texp, airmass, date, hour, archivo
		    line = "%-15s %10s %10s %8.2f %4.2f %8s %11s %s\n" % (obname, ra, delta, float(texp), float(airmass), date, hour, archivo)
		    f.write(line)
		elif h[0].header['IMAGETYP'] == 'flatlamp':
		    flats.append(archivo)
		    mjd, mjd0 = mjd_fromheader2(hd)
		    flat_ref_dates.append( mjd )
		elif h[0].header['IMAGETYP'] == 'arclamp':
		    ThAr_ref.append(archivo)
		    mjd, mjd0 = mjd_fromheader2(hd)
		    ThAr_ref_dates.append( mjd )
    
    flat_ref_dates,ThAr_ref,ThAr_ref_dates = np.array(flat_ref_dates),np.array(ThAr_ref),np.array(ThAr_ref_dates)
    flats = np.array(flats)
    IS = np.argsort(flat_ref_dates)
    flat_ref_dates = flat_ref_dates[IS]
    flats = flats[IS]
    #for i in range(len(flats)):
    #	print 'flat',flats[i], flat_ref_dates[i]
    IS = np.argsort(ThAr_ref_dates)
    ThAr_ref_dates = ThAr_ref_dates[IS]
    ThAr_ref = ThAr_ref[IS]


    return np.array(flats), np.array(ThAr_ref), sim_sci, np.array(ThAr_ref_dates), obnames, exptimes

def mjd_fromheader2(hd):
    """
    return modified Julian date from header
    """

    secinday = 24*3600.0

    date    = hd['DATE-OBS']
    ut      = hd['UTC']		    
    
    datetu = date.replace('-',':')

    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[5:7]),int(datetu[8:10]))
    ut        = (float(ut[:2])*3600. + float(ut[3:5])*60. + float(ut[6:]))
    mjd_start = mjd + ut / secinday

    fraction = 0.5
    texp     = float(hd['EXPTIME']) #sec

    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0

def MedianCombine(ImgList, chip=1):
    """
    Median combine a list of images
    """

    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    h = pyfits.open(ImgList[0])

    d1 = h[chip].data
    h1 = h[0].header
    #d1 = b_col(d1)
    d1 = OverscanTrim(d1,h[chip].header['DATASEC'])
    
    factor = 1.25
    if (n < 3):
        factor = 1

    ron1,gain1 = float(h1['CCDRN0'+str(int(chip))]),float(h1['CCDGN0'+str(int(chip))])

    ron1 = factor * ron1 / np.sqrt(n)
    if n>1:
	    for i in range(n-1):
		td = pyfits.open(ImgList[i+1])
		#print ImgList[i+1],td[2].data.shape
		td = OverscanTrim(td[chip].data,td[chip].header['DATASEC'])
		d1 = np.dstack((d1,td))
		
	    d1 = np.median(d1,axis=2)

    return d1, ron1, gain1

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
			spec_flat[py,px] = np.median(tf2)
	#for i in range(len(traces)):
	#	plot(spec_flat[i])
	#	show()
	for i in range(len(traces)):
		ss =  scipy.signal.medfilt(spec_flat[i],15)
		spec_flat[i] /= ss.max()
	#spec_flat_o = spec_flat.copy()
	#for i in range(len(traces)):
	#	spec_flat[i] = spec_flat[i]/scipy.signal.medfilt(spec_flat[i],301)
	#	spec_flat[i] = spec_flat[i]/sm.max()
	return nflat,spec_flat


def get_cont(x,y,n=1,sl=1.,sh=5.):
	orilen = len(x)
	coef = np.polyfit(x,y,n)
	res = y - np.polyval(coef,x)
	IH = np.where(res>0)[0]
	IL = np.where(res<0)[0]
	dev = np.mean(res[IH])
	I = np.where((res>-sl*dev) & (res<sh*dev))[0]
	J1 = np.where(res<=-sl*dev)[0]
	J2 = np.where(res>=sh*dev)[0]
	J = np.unique(np.hstack((J1,J2)))
	cond = True
	if len(J)==0 or len(x)< .3*orilen:
		cond=False
	while cond:
		x = np.delete(x,J)
		y = np.delete(y,J)
		coef = np.polyfit(x,y,n)
		res = y - np.polyval(coef,x)
		IH = np.where(res>0)[0]
		IL = np.where(res<0)[0]
		dev = np.mean(res[IH])
		I = np.where((res>-sl*dev) & (res<sh*dev))[0]
		J1 = np.where(res<=-sl*dev)[0]
		J2 = np.where(res>=sh*dev)[0]
		J = np.unique(np.hstack((J1,J2)))
		cond = True
		if len(J)==0 or len(x)< .1*orilen:
			cond=False
	return coef,x,y

def get_cont2(x,y,sl=2.,sh=5.):
	orilen = len(x)
	xo = x.copy()
	yo = y.copy()
	#J1 = np.where((x>5870) & (x<5910))[0]
	#x = np.delete(x,J1)
	#y = np.delete(y,J1)

	#J2 = np.where((x>5160) & (x<5210))[0]
	#x = np.delete(x,J2)
	#y = np.delete(y,J2)
	
	yy     = scipy.signal.medfilt(y,11)
	#lowess = robjects.r("lowess")
	#approx = robjects.r("approx")
	#Temp = lowess(x,yy,f=0.1,iter=10)
	#pred = np.array( approx(Temp[0],Temp[1],xout=x, method="linear", rule=2) )[1]

	pred = lowess(yy, x,frac=0.2,it=10,return_sorted=False)

	res = y - pred
	IH = np.where(res>0)[0]
	IL = np.where(res<0)[0]
	dev = np.mean(res[IH])
	I = np.where((res>-sl*dev) & (res<sh*dev))[0]
	J1 = np.where(res<=-sl*dev)[0]
	J2 = np.where(res>=sh*dev)[0]
	J = np.unique(np.hstack((J1,J2)))
	cond = True
	if len(J)==0 or len(x)< .3*orilen:
		cond=False
	#cond = False

	while cond:
		x = np.delete(x,J)
		y = np.delete(y,J)

		yy     = scipy.signal.medfilt(y,11)
		#lowess = robjects.r("lowess")
		#approx = robjects.r("approx")
		#Temp = lowess(x,yy,f=0.2,iter=10)
		#pred = np.array( approx(Temp[0],Temp[1],xout=x, method="linear", rule=2) )[1]
		pred = lowess(yy, x,frac=0.2,it=10,return_sorted=False)
		
		res = y - pred
		IH = np.where(res>0)[0]
		IL = np.where(res<0)[0]
		dev = np.mean(res[IH])
		I = np.where((res>-sl*dev) & (res<sh*dev))[0]
		J1 = np.where(res<=-sl*dev)[0]
		J2 = np.where(res>=sh*dev)[0]
		J = np.sort(np.unique(np.hstack((J1,J2))))
		cond = True
		if len(J)==0 or len(x)< .1*orilen:
			cond=False
	#pred = np.array( approx(Temp[0],Temp[1],xout=xo, method="linear", rule=2) )[1]
	tck1 = scipy.interpolate.splrep(x,pred,k=1)
	pred = scipy.interpolate.splev(xo,tck1)
	#plot(xo,yo)
	#plot(xo,pred)
	#show()
	#pred = lowess(yy, x,frac=0.4,it=10,return_sorted=False)
	return pred,x,y

def bspline(points, n=100, degree=3, periodic=False):

	x = points[:,0]
	y = points[:,1]

	t = range(len(points))
	ipl_t = np.linspace(0.0, len(points) - 1, n)

	x_tup = scipy.interpolate.splrep(t, x, k=3)
	y_tup = scipy.interpolate.splrep(t, y, k=3)

	x_list = list(x_tup)
	xl = x.tolist()
	x_list[1] = xl + [0.0, 0.0, 0.0, 0.0]

	y_list = list(y_tup)
	yl = y.tolist()
	y_list[1] = yl + [0.0, 0.0, 0.0, 0.0]

	x_i = scipy.interpolate.splev(ipl_t, x_list)
	y_i = scipy.interpolate.splev(ipl_t, y_list)
	return x_i,y_i
