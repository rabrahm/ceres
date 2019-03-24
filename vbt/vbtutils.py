from __future__ import print_function
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

def get_thar_offsets(lines_thar, order_dir='wavcals/', pref='order_', suf='.iwdat', ior=20,fior=45, delt_or=3, del_width=200.,binning=1):
    xcs = []
    for ii in range(ior,fior):
        thar_order = lines_thar[ii]
        xct = []
        for order in range(ii-delt_or,ii+delt_or):
            order_s = str(order)
            if (order < 10):
                order_s = '0' + order_s
            if os.access(order_dir+pref+order_s+suf,os.F_OK):
                f = open(order_dir+pref+order_s+suf,'r')
                llins = f.readlines()
                if len(llins)>5:
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
                else:
                    xc = np.array([])
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
    orders_offset = -delt_or + np.argmax(maxes)
    rough_shift = maxvels[np.argmax(maxes)]

    return orders_offset, rough_shift

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
    objects  = []
    darks    = []
    thars    = []
    flat     = []
    f = open(log,'w')
    archs = glob.glob(path+'*.fits')
    if os.access(path+'bad_files.txt',os.F_OK):
        ff = open(path + 'bad_files.txt', 'r')
        bfiles = ff.readlines()
    else:
        bfiles = []
    for arch in archs:
        use = True
        for bf in bfiles:
            if arch == path + bf[:-1]:
                print('Dumped file', arch)
                use = False
                break
        if use:
            h = pyfits.open(arch)
            header = pyfits.getheader(arch)
            name = header['OBJECT']

            if True:
                if header['IMAGETYP'] == 'object' or header['IMAGETYP'] == 'obj' or header['IMAGETYP'] == 'solar':
                    expt = header['EXPTIME']
                    try:
                        ra   = header['TELRA']
                    except:
                        ra   = header['RA']
                    try:
                        dec  = header['TELDEC']
                    except:
                        dec  = header['DEC']
                    ams  = header['AIRMASS']
                    date = header['DATE-OBS'].split('T')[0]
                    UT   = header['DATE-OBS'].split('T')[1]
                    line = "%-15s %10s %10s %8.2f %4.2f %8s %8s %s\n" % (name, ra, dec, expt, ams, date, UT, arch)
                    f.write(line)
                    objects.append(arch)
                elif  header['IMAGETYP'] == 'comp':
                    thars.append(arch)
                elif  header['IMAGETYP'] == 'zero':
                    biases.append(arch)
                elif  header['IMAGETYP'] == 'flat':
                    flat.append(arch)
                elif  header['IMAGETYP'] == 'dark':
                    darks.append(arch)

            h.close()
    f.close()
    return biases, flat, objects, thars, darks

def MedianCombine(ImgList, zero_bo=False, zero='Bias.fits', dark_bo=False, darks=[], flat_bo=False, flat='Flat.fits',bsec=[0,50,4146,4196]):
    """
    Median combine a list of images
    """
    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    h = pyfits.open(ImgList[0])[0]
    d = h.data[0]
    d = OverscanTrim(d,bsec)

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

    ronoise = factor * h.header['RDNOISE'] / np.sqrt(n)
    gain = h.header['GAIN']

    if (n == 1):
        return d, ronoise, gain

    else:
        for i in range(n-1):
            h = pyfits.open(ImgList[i+1])[0]
            if flat_bo:
                d = np.dstack((d,(OverscanTrim(h.data[0],bsec)-Master-Dark)/Flat))
            else:

                d = np.dstack((d,OverscanTrim(h.data[0],bsec)-Master-Dark))

    return np.median(d,axis=2), ronoise, gain


def OverscanTrim(d,bsec,bin=1.):
    """
    Overscan correct and Trim a refurbished DuPont image
    """


    t1 = d[:,bsec[0]:bsec[1]]
    t2 = d[:,bsec[2]:bsec[3]]

    nd = d[:,bsec[1]:bsec[2]].copy()

    overscan1 = np.median(np.dstack((t1,t2)))
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
    datetu = h['DATE-OBS'].split('T')[0]
    timetu = h['DATE-OBS'].split('T')[1]
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
