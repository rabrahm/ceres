import matplotlib
matplotlib.use("Agg")
from astropy.io import fits as pyfits
import numpy as np
import scipy
import glob
import os
import matplotlib.pyplot as plt
import sys
sys.path.append("../utils/GLOBALutils")
import GLOBALutils

from pylab import *

import pycurl

def FileClassify(dir, log,binning):
    """
    Classifies all files in a directory and writes a night log of science images
    """

    # define output lists
    thars      = []
    biases     = []
    quartz     = []
    iodines    = []
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

    all_files = glob.glob(dir+"/pfs*fits")

    for archivo in all_files:
        dump = False
        for bf in bad_files:
            if archivo == bf:
                dump = True
                break

        h = pyfits.open(archivo)
        if dump == False and h[0].header['BINNING']==binning:
            if h[0].header['EXPTYPE'] == 'Object' and h[0].header['OBJECT'] != 'ThAr' and h[0].header['OBJECT'] != 'WideFlat':
                science.append(archivo)
                obname  = h[0].header['OBJECT']
                ra      = h[0].header['RA']
                delta   = h[0].header['DEC']
                airmass = h[0].header['AIRMASS']
                texp    = h[0].header['EXPTIME']
                date    = h[0].header['UT-DATE']
                hour    = h[0].header['UT-TIME']
                obnames.append( obname )
                exptimes.append( texp )
                line = "%-15s %10s %10s %8.2f %4.2f %8s %8s %s\n" % (obname, ra, delta, texp, airmass, date, hour, archivo)
                f.write(line)
            elif h[0].header['EXPTYPE'] == 'Bias':
                print h[0].data.shape
                biases.append(archivo)
            elif h[0].header['EXPTYPE'] == 'Quartz' or h[0].header['EXPTYPE'] == 'Flat':
                if h[0].header['IODINE'] == 'Window':
                    quartz.append(archivo)
                elif h[0].header['IODINE'] == 'I-Cell':
                    iodines.append(archivo)
            elif h[0].header['EXPTYPE'] == 'ThAr' or h[0].header['OBJECT'] == 'ThAr':
                thars.append(archivo)
                mjd, mjd0 = mjd_fromheader(h)
                thar_dates.append( mjd )

    f.close()
    return biases, quartz, iodines, science, thars, thar_dates, obnames, exptimes

def mjd_fromheader(h):
    """
    return modified Julian date from header
    """

    datetu = h[0].header['UT-DATE']
    ut     = h[0].header['UT-TIME']
    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[5:7]),int(datetu[8:]))
    ut = float(ut[:2]) + float(ut[3:5])/60. + float(ut[6:])/3600.
    mjd_start = mjd + ut/24.0
    secinday = 24*3600.0
    fraction = .5
    texp     = h[0].header['EXPTIME'] #sec
    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0

def MedianCombine(ImgList, bs, os, bias = 0.):
    """

    Median combine a list of images

    """

    n = len(ImgList)
    #if n==0:
    #    raise ValueError("empty list provided!")
    if n == 0:
        print "\t\t\tWarning: 0 biases"
        return 0, 0, 0
    h = pyfits.open(ImgList[0])[0]
    d = h.data
    d = OverscanTrim(d,bs,os)
    d -= bias

    factor = 1.25
    if (n < 3):
        factor = 1

    ronoise = factor * h.header['ENOISE'] / np.sqrt(n)
    gain    = h.header['EGAIN']


    if (n == 1):
        return d, ronoise, gain
    else:
        for i in range(n-1):
            h = pyfits.open(ImgList[i+1])[0]
            d2 = OverscanTrim(h.data,bs,os)
            d = np.dstack( (d ,OverscanTrim(h.data,bs,os) -bias ))
        return np.median(d,axis=2), ronoise, gain

def OverscanTrim(d,bs,ts):
    """
    Overscan correct and Trim a refurbished PFS image
    """
    # bias has no significant structure, so a single median suffices, I think
    # overscan = [0:49] [2097:2145]
    #bs = [0,2032,2160]
    #ts = [1,2032,2160]
    if bs[0] == 0:
        bss = d[bs[1]:bs[2],:]
        bss = scipy.median(bss,axis=0)[:ts[1]]
        bss = np.tile(bss,bs[1])
        bss = np.reshape(bss,[bs[1],ts[1]])
        newdata = d[:bs[1],:ts[1]] - bss
    else:
        bss = d[:,bs[1]:bs[2]]
        bss = scipy.median(bss,axis=1)[:ts[1]]
        bss = np.tile(bss,ts[1])
        bss = np.reshape(bss,[ts[1],bs[1]]).T
        newdata = d[:ts[1],:bs[1]] - bss

    return newdata
