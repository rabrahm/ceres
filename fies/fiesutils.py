from __future__ import print_function
import matplotlib
matplotlib.use("Agg")
from astropy.io import fits as pyfits
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


def get_thar_offsets(lines_thar, order_dir='wavcals/', pref='order_', suf='.iwdat', delt_or=10, del_width=200.,binning=1):
    start_or = int(.5*delt_or)
    xcs = []
    for ii in range(delt_or,len(lines_thar)-delt_or):
        thar_order = lines_thar[ii]
        xct = []
        for order in range(ii-start_or,ii+start_or):
            order_s = str(order)
            if (order < 10):
                order_s = '0' + order_s
            if os.access(order_dir+pref+order_s+suf,os.F_OK):
                f = open(order_dir+pref+order_s+suf,'r')
                llins = f.readlines()
                if True:
                    pixel_centers_0 = []
                    for line in llins:
                        w = line.split()
                        nlines = int(w[0])
                        for j in range(nlines):
                            pixel_centers_0.append(float(w[2*j+1])*1./float(binning))
                    pixel_centers_0 = np.array(pixel_centers_0).astype('int')
                    #plot(thar_order)
                    #plot(pixel_centers_0,thar_order[pixel_centers_0],'ro')
                    #print order, order_s
                    #show()
                    ml = np.array(pixel_centers_0) - 2
                    mh = np.array(pixel_centers_0) + 2
                    if len(ml)>0:
                        xc,offs = GLOBALutils.XCorPix( thar_order, ml, mh, del_width=del_width)
                    else:
                        xc = np.zeros(len(offs))

            if len(xct) == 0:
                xct = xc.copy()
            else:
                xct = np.vstack((xct,xc))

        if len(xcs) == 0:
            xcs = xct.copy()
        else:
            xcs += xct
    maxes, maxvels = [],[]
    for i in range(xcs.shape[0]):
        maxes.append(xcs[i].max())
        maxvels.append(offs[np.argmax(xcs[i])])
        #plot(offs,xcs[i])
    #show()

    maxes,maxvels = np.array(maxes),np.array(maxvels)
    orders_offset = -start_or + np.argmax(maxes)
    rough_shift = maxvels[np.argmax(maxes)]

    return orders_offset, rough_shift

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

def FileClassify(diri, log,binning=1,mode='F1', dark_corr=False):
    """

    Classifies all files in a directory and writes a night log of science images

    """

    # define output lists
    sim_sci        = []
    biases         = []
    flats          = []
    ThAr_ref       = []
    ThAr_ref_dates = []
    ThAr_co       = []
    ThAr_co_dates = []
    ThAr_sim       = []
    ThAr_sim_dates = []
    flat_ref_dates = []
    bias_ref_dates = []
    obnames        = []
    exptimes       = []
    darks          = []
    flats_co               = []
    flats_co_dates  = []

    sdarks = []
    if dark_corr and os.access(diri+'/darks.txt',os.F_OK):
        fd = open(diri+'/darks.txt','r')
        ds = fd.readlines()
        for dk in ds:
            sdarks.append(diri+dk[:-1])
    sdarks = np.array(sdarks)

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
        isdark=False
        for df in sdarks:
            if archivo == df:
                darks.append(archivo)
                isdark=True

        if dump == False and isdark == False:
            print(archivo)
            h = pyfits.open(archivo)
            hd = pyfits.getheader(archivo)
            if int(h[0].header['DETXBIN']) == binning and int(h[0].header['DETYBIN']) == binning and (mode in h[0].header['FIFMSKNM']) and h[0].header['IMAGETYP'] != 'COUNTTEST':
                print(archivo, h[0].header['IMAGETYP'], h[0].header['SHSTAT'], h[0].header['EXPTIME'], h[0].header['OBJECT'], h[0].header['TCSTGT'], int(h[0].header['DETYBIN']))

                if h[0].header['IMAGETYP'] == 'BIAS':
                    biases.append(archivo)
                    mjd, mjd0 = mjd_fromheader2(h)
                    bias_ref_dates.append( mjd )
                elif h[0].header['IMAGETYP'] == 'FLAT':
                    flats.append(archivo)
                    mjd, mjd0 = mjd_fromheader2(h)
                    flat_ref_dates.append( mjd )
                    if h[0].header['FICARMID'] == 6 and h[0].header['FILMP1'] == 1 and h[0].header['FILMP6']==0:
                        flats_co.append(archivo)
                        mjd, mjd0 = mjd_fromheader2(h)
                        flats_co_dates.append( mjd )
                    else:
                        flats.append(archivo)
                        mjd, mjd0 = mjd_fromheader2(h)
                        flat_ref_dates.append( mjd )
                        sc = pyfits.getdata(archivo)
                        #plot(sc[1000])
                elif h[0].header['IMAGETYP'] == 'WAVE':
                    ThAr_ref.append(archivo)
                    mjd, mjd0 = mjd_fromheader2(h)
                    ThAr_ref_dates.append( mjd )


                elif ((mode=='F3' or mode=='F4') and h[0].header['FICARMID'] == 6 and h[0].header['FILMP4'] == 0 and h[0].header['FILMP7']==1)\
                   or (mode=='F1' and h[0].header['FICARMID'] == 2 and h[0].header['FILMP4'] == 0 and h[0].header['FILMP7']==1):
                    ThAr_ref.append(archivo)
                    mjd, mjd0 = mjd_fromheader2(h)
                    ThAr_ref_dates.append( mjd )

                elif h[0].header['FICARMID'] == 6 and h[0].header['FILMP4'] == 1 and h[0].header['FILMP7']==0:
                    ThAr_co.append(archivo)
                    mjd, mjd0 = mjd_fromheader2(h)
                    ThAr_co_dates.append( mjd )

                elif h[0].header['FICARMID'] == 6 and h[0].header['FILMP4'] == 1 and h[0].header['FILMP7']==1:
                    ThAr_sim.append(archivo)
                    mjd, mjd0 = mjd_fromheader2(h)
                    ThAr_sim_dates.append( mjd )

                elif (mode=='F3' and h[0].header['FICARMID'] == 2) or (mode == 'F1' and h[0].header['FICARMID'] == 5)\
                   or (mode=='F4' and (h[0].header['FICARMID'] == 5 or h[0].header['FICARMID'] == 4)):
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
                    if  h[0].header['FILMP4'] == 1:
                        simult = 'SIMULT'
                    else:
                        simult = 'NO_SIMULT'
                    line = "%-15s %10s %10s %8.2f %4.2f %8s %11s %s %s\n" % (obname, ra, delta, texp, airmass, date, hour, archivo, simult)
                    f.write(line)

    #show()
    flat_ref_dates = np.array(flat_ref_dates)
    flats = np.array(flats)
    IS = np.argsort(flat_ref_dates)
    flat_ref_dates = flat_ref_dates[IS]
    flats = flats[IS]
    #for i in range(len(flats)):
    #       print 'flat',flats[i], flat_ref_dates[i]

    bias_ref_dates = np.array(bias_ref_dates)
    biases = np.array(biases)
    IS = np.argsort(bias_ref_dates)
    bias_ref_dates = bias_ref_dates[IS]
    biases = biases[IS]
    #for i in range(len(biases)):
    #       print 'bias',biases[i], bias_ref_dates[i]
    f.close()

    return biases, np.array(flats), np.array(ThAr_ref), sim_sci, np.array(ThAr_ref_dates), obnames, exptimes, np.array(darks), np.array(flats_co), np.array(flats_co_dates),np.array(ThAr_sim), np.array(ThAr_sim_dates),np.array(ThAr_co), np.array(ThAr_co_dates)

def get_darktimes(darks):
    times = []
    for dark in darks:
        hd = pyfits.getheader(dark)
        times.append(hd['EXPTIME'])
    return np.unique(np.sort(np.array(times))), np.array(times)

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
    return hd['RDNOISE'], hd['GAIN']

def MedianCombine(ImgList, zero='none', binning=1, oii=100, off=2148):
    """
    Median combine a list of images
    """

    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    h = pyfits.open(ImgList[0])

    d1 = h[1].data
    h1 = h[1].header
    d1 = OverscanTrim(d1, binning=binning,ii=oii, ff=off)
    if zero != 'none':
        z = pyfits.open(zero)[0]
        d1 -= z.data

    factor = 1.25
    if (n < 3):
        factor = 1

    ron1,gain1 = get_RONGAIN(h[1].header)

    ron1 = factor * ron1 / np.sqrt(n)
    if n>1:
        for i in range(n-1):
            td = pyfits.open(ImgList[i+1])
            if zero == 'none':
                d1 = np.dstack((d1,OverscanTrim(td[1].data, binning=binning, ii=oii, ff=off)))

            else:
                d1 = np.dstack((d1,OverscanTrim(td[1].data, binning=binning, ii=oii, ff=off)-z.data))
        d1 = np.median(d1,axis=2)
    return d1, ron1, gain1

def OverscanTrim(dat,binning=1,ii=100,ff=2148):
    """
    Overscan correct and Trim a refurbished FEROS image
    """
    #ff = 2098
    #ii = 50
    ff = int(np.around(ff/binning))
    ii = int(np.around(ii/binning))
    os = dat[:,ff:]
    s = np.median(os)
    newdata = dat[:,ii:ff].copy() - s

    return newdata
