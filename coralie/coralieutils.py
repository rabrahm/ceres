import pyfits
import numpy as np
import scipy
from scipy import optimize
import copy
import glob
import os
import matplotlib.pyplot as plt
import sys
sys.path.append("../utils/GLOBALutils")
import GLOBALutils

from rpy2 import robjects
import rpy2.robjects.numpy2ri
r = robjects.r

import pycurl


def MedianCombine(ImgList,ZF=0.):
    """

    Median combine a list of images

    """

    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    h = pyfits.open(ImgList[0])[0]
    d = h.data
    d = OverscanTrim(d) - ZF
    
    factor = 1.25
    if (n < 3):
        factor = 1

    ronoise = factor * h.header['HIERARCH ESO CORA CCD RON'] / np.sqrt(n)
    gain    = h.header['HIERARCH ESO CORA CCD GAIN']

    if (n == 1):
        return d, ronoise, gain
    else:
        for i in range(n-1):
            h = pyfits.open(ImgList[i+1])[0]
            d = np.dstack((d,OverscanTrim(h.data)-ZF))
        return np.median(d,axis=2), ronoise, gain

def OverscanTrim(d):
    """
    Overscan correct and Trim a refurbished CORALIE image
    """
    # bias has no significant structure, so a single median suffices, I think
    # overscan = [0:49] [2097:2145]
    overscan = 0.5*(np.median(d[:,0:49]) + np.median(d[:,2097:2146]))
    newdata = d[:,50:2097] - overscan
    return newdata

def getObName(h):
    """
    Get name of object under consideration
    """
    obname = h[0].header['HIERARCH ESO OBS REFNOCOD'].upper().replace(' ','')
    observer = h[0].header['OBSERVER']
    if (obname == 'H726273'):
        obname = 'HD72673'
    if (obname == 'HAT'):
        obname = 'HATS563-036'
    if (obname == '9999.999'):
        # try the header that James usually uses:
        obname =  h[0].header['HIERARCH ESO OBS TARG CODE'].upper().replace(' ','')
    if (obname == '9999.999'):
        # seems like observer was too lazy to type name, so use SIMBAD
        (th,tfile) = tempfile.mkstemp(prefix='CP', text=True)
        tf = open(tfile,'w')
        tf.write("output console=off\n")
        tf.write("output script=off\n")
        tf.write("output error=merge\n")
        tf.write("set limit 1\n")
        tf.write("format object fmt1 \"%IDLIST(HD|HIP|1) | %OTYPELIST(S) | %SP(S)\"\n")
        tf.write("result full\n")
        tf.write("query sample region(circle, %s %s,5s) & otype='Star'\n" % (h[0].header['HIERARCH ESO TEL TARG ALPHA'],h[0].header['HIERARCH ESO TEL TARG DELTA'] ))
        #tf.write("set radius 5s\n")
        tf.close()
        values = [
            ("scriptFIle", (pycurl.FORM_FILE, tfile))
            ]
        output = StringIO.StringIO() 
        c = pycurl.Curl()
        
        c.setopt(pycurl.URL, "http://simbad.harvard.edu/simbad/sim-script")
        c.setopt(c.HTTPPOST, values)
        c.setopt(pycurl.WRITEFUNCTION, output.write) 
        c.perform()
        c.close()
        
        result = output.getvalue() 
        lines = result.split('\n') 
        result = lines[len(lines)-3]
        if (result.count('No') > 0):
            # build alternate obname based on ra-dec
            ra_s  = h[0].header['HIERARCH ESO OBS ALPHACAT'].replace('h','').replace('m','')
            dec_s = h[0].header['HIERARCH ESO OBS DELTACAT'].replace(':','')
            hour_l = len( h[0].header['HIERARCH ESO OBS ALPHACAT'].split('h')[0] )
            if (hour_l == 1):
                obname = 'J_0'+ra_s+dec_s
            elif (hour_l == 2):
                obname = 'J_'+ra_s+dec_s
            else:
                raise ValueError("Unexpected length for RA string from header")
        else:
            obname = lines[len(lines)-3].split('|')[0].replace(" ","")
        os.remove(tfile)

    return obname

def FileClassify(dir, log):
    """
    
    Classifies all files in a directory and writes a night log of science images

    """

    # define output lists
    simThAr_sci = []
    simFP_sci = []
    biases = []
    ob_flats = []
    co_flats = []
    ob_loc = []
    co_loc = []
    ThAr_ref = []
    FP_ref = []
    ThAr_ref_dates = []
    ThFP_ref_dates = []
    obnames = []
    obnames_FP = []
    exptimes = []
    exptimes_FP = []
    f = open(log,'w')

    bad_files = []
    if os.access(dir+'bad_files.txt',os.F_OK):
    	bf = open(dir+'bad_files.txt')
	linesbf = bf.readlines()
	for line in linesbf:
		bad_files.append(dir+line[:-1])
	bf.close()
    
    all_files = glob.glob(dir+"/CORALIE*fits")
    for archivo in all_files:

	dump = False
	for bf in bad_files:
		if archivo == bf:
			dump = True
			break
	if dump == False:
		h = pyfits.open(archivo)
		hd = pyfits.getheader(archivo)
		if h[0].header['HIERARCH ESO TPL TYPE'] == 'OBTH' or h[0].header['HIERARCH ESO TPL TYPE'] == 'OBFP':
		    obname = getObName(h)
		    ra     = h[0].header['HIERARCH ESO OBS ALPHACAT']
		    delta  = h[0].header['HIERARCH ESO OBS DELTACAT']
		    airmass= h[0].header['HIERARCH ESO OBS TARG AIRMASS']
		    texp   = h[0].header['HIERARCH ESO OBS TEXP']
		    date   = h[0].header['HIERARCH ESO CORA SHUTTER START DATE']
		    hour   = h[0].header['HIERARCH ESO CORA SHUTTER START HOUR'] 
		    line = "%-15s %10s %10s %8.2f %4.2f %8s %7.4f %s\n" % (obname, ra, delta, texp, airmass, date, hour, archivo)
		    f.write(line)
		    simThAr_sci.append(archivo)
		    obnames.append( obname )
		    exptimes.append( texp )
	 
		elif h[0].header['HIERARCH ESO TPL TYPE'] == 'BIAS':
		    biases.append(archivo)
		elif h[0].header['HIERARCH ESO TPL TYPE'] == 'FFO':
		    ob_flats.append(archivo)
		elif h[0].header['HIERARCH ESO TPL TYPE'] == 'FFC':
		    co_flats.append(archivo)
		elif h[0].header['HIERARCH ESO TPL TYPE'] == 'LOCO':
		    ob_loc.append(archivo)
		elif h[0].header['HIERARCH ESO TPL TYPE'] == 'LOCC':
		    co_loc.append(archivo)
		elif h[0].header['HIERARCH ESO TPL TYPE'] == 'THA2':
		    ThAr_ref.append(archivo)
		    mjd, mjd0 = mjd_fromheader(h)
		    ThAr_ref_dates.append( mjd )
		elif h[0].header['HIERARCH ESO TPL TYPE'] == 'THFP':
		    FP_ref.append(archivo)
		    mjd, mjd0 = mjd_fromheader(h)
		    ThFP_ref_dates.append( mjd )
    f.close()

    return biases, ob_flats, co_flats, ob_loc, co_loc, ThAr_ref, FP_ref, simThAr_sci, simFP_sci, ThAr_ref_dates, ThFP_ref_dates, obnames, obnames_FP, exptimes, exptimes_FP

def mjd_fromheader(h):
    """
    return modified Julian date from header
    """
    
    datetu = h[0].header['HIERARCH ESO CORA SHUTTER START DATE'] 
    ut     = h[0].header['HIERARCH ESO CORA SHUTTER START HOUR']
    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[4:6]),int(datetu[6:8]))

    mjd_start = mjd + ut/24.0

    secinday = 24*3600.0
    fraction = h[0].header['HIERARCH ESO CORA PM FLUX TMMEAN']
    texp     = h[0].header['HIERARCH ESO OBS TEXP'] #sec

    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0

def XC_Final_Fit_Rot( X, Y, ldc = 0.8, vsini = 10.0 ):
    """
    Fits a Gaussian and Gauss-Hermite series
    Higher order in the expansion is horder

    """
    f0 = 0.1
    vel0 = X[len(X)/2]
	 
    def conv(x,g,x1,r1):
	xi = np.argmin(x1**2)
	r1 = np.hstack((r1[xi:],r1[:xi]))
	tg = np.zeros(len(g))
	for i in range(len(x)):
		tg[i] = np.add.reduce(g*r1)*(x[1]-x[0])
		r1 = np.hstack((r1[-1:],r1[:-1]))
	return tg

    def fitfunc(p,x,e):
	A    =  p[0]
	rv   =  p[1]
	vrot =  p[2]
	s    =  p[3]
	
	g1 = np.exp(-0.5*((x-rv)/s)**2)/np.sqrt(2*np.pi*s*s)
	d = x[1]-x[0]
	x1 = np.arange(len(x))*d
	x1 -= int(np.round(x1.mean()))
	I = np.where(x1**2 < vrot**2)[0]
	c1 = 2.*(1. - e) / (np.pi * vrot * (1. - e/3.))
	c2 = .5 * e / (vrot * (1.-e/3.))
	r1 = np.zeros(len(x1))
	r1[I] = (c1*np.sqrt(1.-(x1[I]/vrot)**2) + c2*(1. - (x1[I]/vrot)**2))

	prof = conv(x,g1,x1,r1)

	ret = 1. - A*prof
	return ret
        
    def errfunc(p, x, y, ldc):
        clutch = 0.0
        mean = p[1]
        if (mean < np.min(x)):
            clutch = 1e10*(1.0 - np.exp(-np.abs(mean-np.min(x)) / 3) )
        if (mean > np.max(x)): 
            clutch = 1e10*(1.0 - np.exp(-np.abs(mean-np.max(x)) / 3) )
        return np.ravel( (fitfunc(p,x,ldc) - y) )  + clutch


    p0 = np.array( [1.,vel0,vsini,1.] )

    
    p1, success = scipy.optimize.leastsq(errfunc,p0, args=(X,Y,ldc))
    #plt.plot(X,Y)
    #plt.plot(X, 1. - p1[0]*np.exp(-0.5*((X-p1[1])/p1[3])**2)/np.sqrt(2*np.pi*p1[3]*p1[3]))
    c1 = 2.*(1. - ldc) / (np.pi * p1[2] * (1. - ldc/3.))
    c2 = .5 * ldc / (p1[2] * (1.-ldc/3.))
    I = np.where(X**2 < p1[2]**2)[0]
    r1 = np.zeros(len(X))
    r1[I] = (c1*np.sqrt(1.-(X[I]/p1[2])**2) + c2*(1. - (X[I]/p1[2])**2))

    #plt.plot(X, 1. - p1[0]*r1)
    #plt.show()
        
    return p1, fitfunc(p1,X,ldc)

def get_ldc(T,G,Z,M,ldfile = 'lin_coe_sloan2.dat'):
	f = np.loadtxt(ldfile)
	
	I = np.argmin( np.sqrt((T - f[:,2])**2)/T + np.sqrt((G - f[:,1])**2)/G + np.sqrt((Z - f[:,3])**2)/np.sqrt(Z**2) )
	return f[I][5]
