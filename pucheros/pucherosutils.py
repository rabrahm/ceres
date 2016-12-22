import sys
base = '../'
sys.path.append(base+"utils/GLOBALutils")
import GLOBALutils

import numpy as np
import scipy
import pyfits
import os
import glob
import tempfile
import StringIO
import pycurl
from pylab import *


def is_there(string, word):
	l=len(word)
	i=0
	ist = False
	while i < len(string)-l:
		if string[i:i+l] == word:
			ist = True
		i+=1
	return ist

def search_name(obj):
	name = obj.split('/')[-1]
        try:
	  name = name.split('_')[1]
        except:
          name = name.split('_')[1]
	#print 'NAME:', name
	return name

def FileClassify(path,log):
	biases   = []
	flats = []
        img_flats = []
        fib_flats = []
	objects  = []
	darks    = []
	thars    = []
	
	
	lines = []
	dates = []
	archs = glob.glob(path+'*.fit')

	bad_files = []
	if os.access(path+'bad_files.txt',os.F_OK):
		bf = open(path+'bad_files.txt')
		linesbf = bf.readlines()
		for line in linesbf:
			bad_files.append(path+line[:-1])

		bf.close()
	for arch in archs:
		dump = False
		for bf in bad_files:
			if arch == bf:
				dump = True
				break
		if not dump:
			
			h = pyfits.open(arch)
			
			#print h[0].header['XBINNING'], h[0].header['YBINNING'], arch
			if h[0].header['XBINNING'] == 1 and h[0].header['YBINNING'] == 1:
				if h[0].header['IMAGETYP'] == 'Light Frame' or h[0].header['IMAGETYP'] == 'LIGHT':
					if 'flat' in arch:
						flats.append(arch)
					else:
						name = h[0].header['OBJECT']
						expt = h[0].header['EXPTIME']
						date = h[0].header['DATE-OBS']
						line = "%-15s %8.2f %8s %s\n" % (name, expt, date, arch)

						ye = float(date[:4])
						mo = float(date[5:7])
						da = float(date[8:10])
						ho = float(date[11:13])-4.0
						mi = float(date[14:15])
						se = float(date[17:])

						lines.append(line)
						dates.append( jd( ye,mo,da,ho,mi,se ) )
						#f.write(line)
					
						if is_there(arch.lower(),'thar')  or is_there(arch.lower(),'th_ar'):
							thars.append(arch)				
						else:
							objects.append(arch)
	
				elif  h[0].header['IMAGETYP'] == 'Bias Frame' or h[0].header['IMAGETYP'] == 'BIAS':
					biases.append(arch)
				elif  (h[0].header['IMAGETYP'] == 'Flat Frame' or h[0].header['IMAGETYP'] == 'FLAT') and arch != 'MasterFlat.fits':
		                        # Now check which kind of flat it is.
		                        # Maybe a "surface" flat...
		                        if(is_there(arch.lower(),'imgflat')):
		                          img_flats.append(arch)
		                        # ...a fibre flat...
		                        elif(is_there(arch.lower(),'fibre')):
		                          fib_flats.append(arch)
		                          # (use them for traces, blaze and col-to-col)
		                          flats.append(arch)
		                        # ...else, it is a screen flat (w/difussor):
		                        else:
					  flats.append(arch)
				elif  h[0].header['IMAGETYP'] == 'Dark Frame' or h[0].header['IMAGETYP'] == 'DARK':
					if h[0].header['EXPTIME']!=0.0:
						darks.append(arch)
	
		h.close()
	lines = np.array(lines)
	dates = np.array(dates)
	I = np.argsort(dates)
	lines = lines[I]
	f = open(log,'w')
	for line in lines:
		f.write(line)
	f.close()
	return biases,flats,img_flats,fib_flats,objects,thars,darks

def get_rg():
	return 9.6,1.6

def MedianCombine(ImgList,zero_bo,zero,dark_bo=False, dlist = []):
	"""
	Median combine a list of images
	"""

	hf = pyfits.getheader(ImgList[0])

	if zero_bo:
		Master = pyfits.getdata(zero)
	if dark_bo:
		Dark = get_dark(dlist,hf['EXPTIME'])

	n = len(ImgList)

	if n==0:
		raise ValueError("empty list provided!")

	d = pyfits.getdata(ImgList[0])
	
	if zero_bo:
		d = d - Master
	if dark_bo:
		d = d - Dark
	
	factor = 1.25

	if (n < 3):
		factor = 1

	#ronoise = factor * h.header['ENOISE'] / np.sqrt(n)
	#gain = h.header['EGAIN']
	ronoise,gain=get_rg()
	if (n == 1):
		return d, ronoise, gain

	else:
		for i in range(n-1):
			h = pyfits.getdata(ImgList[i+1])
			if zero_bo:
				h = h-Master
			if dark_bo:
				h = h-Dark
			d = np.dstack((d,h))

		return np.median(d,axis=2), ronoise/np.sqrt(n), gain

def get_dark(darks,t):
	exact = 0
	dts = []
	for dark in darks:
		hd = pyfits.getheader(dark)
		dt = hd['EXPTIME']
		dts.append(dt)
		if dt == t:
			#print 'dark:',dark
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
			#print darks[I[0]]
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


def jd(y,m,d,h,mins,s):
       "Julian day is calculated here if it's needed"
       MY = (m-14)/12
       y = MY+y
       return ( 1461 * ( y + 4800 ) ) / 4  + (  367 * ( m - 2 - 12*MY ) ) / 12 - ( 3 * ( ( y + 4900 ) / 100 ) ) / 4 + d -32077.5 + htosec(h,mins,s)/86400.0

def htosec(h,m,s):
       "transform from hour,minute and seconds, to seconds"
       return s+60.0*(m+60.0*h)

def fit_blaze(w,f,n=5):
	warnings.simplefilter('ignore', np.RankWarning)
	li = len(w)
	co = np.polyfit(w,f,n)
	res = f - np.polyval(co,w)
	dev = np.sqrt(np.var(res))
	J1 = np.where(res < -1.5*dev)[0]
	J2 = np.where(res > 3*dev)[0]
	J = np.hstack((J1,J2))
	J = np.sort(J)
	I = np.where( (res >= -1.5*dev) & (res <= 3*dev) )[0]
	cond = True
	if len(J)==0 or len(I) < .3*li:
		cond = False
	while cond:
		w,f = w[I],f[I]
		co = np.polyfit(w,f,n)
		res = f - np.polyval(co,w)
		dev = np.sqrt(np.var(res))
		J1 = np.where(res < -1.5*dev)[0]
		J2 = np.where(res > 3*dev)[0]
		J = np.hstack((J1,J2))
		J = np.sort(J)
		I = np.where( (res >= -1.5*dev) & (res <= 3*dev) )[0]
		cond = True
		if len(J)==0 or len(I) < .3*li:
			cond = False
	return co

def mjd_fromheader(h):
    """
    return modified Julian date from header
    """
    
    datetu = h[0].header['DATE-OBS'] 
    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[:4]),int(datetu[5:7]),int(datetu[8:10]))
    ut = float(datetu[11:13]) + float(datetu[14:16])/60. + float(datetu[17:])/3600.
    mjd_start = mjd + ut/24.0

    secinday = 24*3600.0
    fraction = 0.5
    texp     = h[0].header['EXPTIME'] #sec

    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0

def get_coords(obname,mjd):
	if obname.lower() == 'alphacent':
		obname = 'alpha cent'
	elif obname.lower() == 'alphaboo':
		obname = 'alpha boo'
	elif obname.lower() == 'hadar' or obname.lower() == 'betacen':
		obname = 'beta cen'
	elif obname.lower() == 'diphda':
		obname = 'bet cet'
	elif obname.lower() == 'betacar':
		obname = 'beta car'
	elif obname.lower() == 'betscl':
		obname = 'bet scl'
	elif obname.lower() == 'bvel':
		obname = 'b vel'
	elif obname.lower() == 'deltasco':
		obname = 'del sco'
	elif obname.lower() == 'delcen':
		obname = 'del cen'
	elif obname.lower() == 'epsilonaqr':
		obname = 'eps aqr'
	elif obname.lower() == 'epspsa':
		obname = 'eps psa'
	elif obname.lower() == 'etahya' or obname.lower() == 'ethahydra':
		obname = 'eta Hya'
	elif obname.lower() == 'etapsa':
		obname = 'eta psa'
	elif obname.lower() == 'etacen':
		obname = 'eta cen'
	elif obname.lower() == 'opup':
		obname = 'o Pup'
	elif obname.lower() == 'etacar':
		obname = 'eta Car'
	elif obname.lower() == 'agcar':
		obname = 'ag Car'
	elif obname.lower() == 'hrcar':
		obname = 'hr Car'
	elif obname.lower() == 'sslep':
		obname = 'ss lep'
	elif obname.lower() == 'thetavir':
		obname = 'theta vir'
	elif obname.lower() == 'mucen':
		obname = 'mu cen'
	elif obname.lower() == 'lesath':
		obname = 'ups sco'
	elif obname.lower() == 'mulup':
		obname = 'mu lup'
	elif obname.lower() == 'chioph':
		obname = 'chi oph'
	elif obname.lower() == 'dlup':
		obname = 'd lup'
	elif obname.lower() == '48lib':
		obname = '48 lib'
	elif obname.lower() == 'iotara':
		obname = 'iot ara'
	elif obname.lower() == 'qvtel':
		obname = 'qv tel'			
	elif obname.lower() == 'taucet':
		obname = 'tau cet'
	elif obname.lower() == 'pi2ori':
		obname = 'pi2 ori'
	elif obname.lower() == 'zetapeg':
		obname = 'zet peg'
	elif obname.lower() == 'tpyx':
		obname = 't pyx'
	elif obname.lower() == 'omicronpup':
		obname = 'omi pup'

	sp,ra,dec = 0,0,0
	(th,tfile) = tempfile.mkstemp(prefix='CP', text=True)
        tf = open(tfile,'w')
	tf.write("output console=off\n")
	tf.write("output script=off\n")
	tf.write("output error=merge\n")
	tf.write("set limit 1\n")
	tf.write("format object fmt1 \"%IDLIST(1) | %OTYPELIST(S) | %SP(S) | %COO(A) | %COO(D) | %PM(A) | %PM(D)\"\n")
	tf.write("result full\n")
	tf.write("query id %s\n" % ( obname ) )
	tf.close()
	values = [("scriptFIle", (pycurl.FORM_FILE, tfile))]
	output = StringIO.StringIO() 
	c = pycurl.Curl()
	c.setopt(pycurl.URL, "http://simbad.harvard.edu/simbad/sim-script")
	c.setopt(c.HTTPPOST, values)
	c.setopt(pycurl.WRITEFUNCTION, output.write)
	cond = True
	while cond:
		try:
			c.perform()
		except:
			print 'Trying again to perform query to SIMBAD'
		else:
			cond = False
	c.close()
	result = output.getvalue() 
	lines = result.split('\n')
	info = lines[6].split('|')

	if 'Unrecogniezd' in info[0] or 'not' in info[0]:
		know = False
	else:
		know = True
		sp,ra,dec,pmra,pmdec = info[2],info[3],info[4],info[5],info[6]

		if '~' in pmra:
			pmra = '0.'
		if '~' in pmdec:
			pmdec = '0.'

		rad = ra.split()
		decd = dec.split()
		ra = float(rad[0])*360./24. + float(rad[1])*6./24. + float(rad[2])/240. +  (float(pmra)/(3600*1000.))*((mjd-51544.5)/365.)
		if float(decd[0])<0:
			dec = -(np.absolute(float(decd[0])) + float(decd[1])/60. + float(decd[2])/3600.) +  (float(pmdec)/(3600*1000.))*((mjd-51544.5)/365.)
		else:
			dec = float(decd[0]) + float(decd[1])/60. + float(decd[2])/3600. +  (float(pmdec)/(3600*1000.))*((mjd-51544.5)/365.)
	return sp,ra,dec,know
