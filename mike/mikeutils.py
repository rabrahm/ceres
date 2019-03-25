from __future__ import print_function
import matplotlib
matplotlib.use("Agg")
from astropy.io import fits as pyfits
import numpy as np
import scipy
from scipy import signal,optimize
import copy
import glob
import os
import matplotlib.pyplot as plt
import sys
sys.path.append("../utils/GLOBALutils")
import GLOBALutils

from pylab import *

def get_hour(hour):
    hh = int(hour/3600.)
    mm = int((hour/3600. - hh) * 60.)
    ss = int(((hour/3600. - hh) * 60. - mm)*60)
    shh,smm,sss = str(hh),str(mm),str(ss)
    if hh<10:
        shh = '0'+shh
    if mm<10:
        smm = '0'+smm
    if ss<10:
        sss = '0'+sss
    hour = shh + ':' + smm + ':' + sss
    return hour

def FileClassify(dir, log,binning):
    """

    Classifies all files in a directory and writes a night log of science images

    """

    # define output lists
    thars      = []
    biases     = []
    milkys     = []
    flatsB     = []
    flatsR     = []
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

    fix_files, ffEXP, ffOBJ = [],[],[]
    if os.access(dir+'fix_files.txt',os.F_OK):
        ff = open(dir+'fix_files.txt')
        linesff = ff.readlines()
        for line in linesff:
            cos = line.split()
            fix_files.append(dir+cos[0])
            ffEXP.append(cos[1])
            ffOBJ.append(cos[2])
        ff.close()
    fix_files, ffEXP, ffOBJ = np.array(fix_files), np.array(ffEXP), np.array(ffOBJ)
    all_files = glob.glob(dir+"/r*fits*")

    for archivo in all_files:
        dump = False
        for bf in bad_files:
            if archivo == bf:
                dump = True
                break
        h = pyfits.open(archivo)
        if archivo[-2:] == 'fz':
            h[0].header = h[1].header

        if dump == False and h[0].header['BINNING']==binning:

            I = np.where(archivo == fix_files)[0]
            if len(I)>0:
                h[0].header['EXPTYPE'] = ffEXP[I[0]]
                h[0].header['OBJECT']  = ffOBJ[I[0]]

            if h[0].header['EXPTYPE'] == 'Object':
                if 'milky' in str(h[0].header['OBJECT']).lower():
                    milkys.append(archivo)
                elif 'flat' in h[0].header['OBJECT'].lower():
                    if 'blue' in h[0].header['OBJECT'].lower():
                        flatsB.append(archivo)
                    else:
                        flatsR.append(archivo)
                else:
                    science.append(archivo)
                    obname  = h[0].header['OBJECT']
                    ra      = h[0].header['RA']
                    delta   = h[0].header['DEC']
                    airmass = h[0].header['AIRMASS']
                    texp    = h[0].header['EXPTIME']
                    date    = h[0].header['UT-DATE']
                    hour    = float(h[0].header['UT-TIME'])
                    hh = int(hour/3600.)
                    mm = int((hour/3600. - hh) * 60.)
                    ss = int(((hour/3600. - hh) * 60. - mm)*60)
                    shh = str(hh)
                    smm = str(mm)
                    sss = str(ss)
                    if hh < 10:
                        shh = '0' + shh
                    if mm < 10:
                        smm = '0' + smm
                    if ss < 10:
                        sss = '0' + sss
                    hour = shh + ':' + smm + ':' + sss
                    obnames.append( obname )
                    exptimes.append( texp )
                    line = "%-15s %10s %10s %8.2f %4.2f %8s %8s %s\n" % (obname, ra, delta, texp, airmass, date, hour, archivo)
                    f.write(line)
            elif h[0].header['EXPTYPE'] == 'Comp':
                if 'milky' in h[0].header['OBJECT'].lower():
                    milkys.append(archivo)
                elif 'flat' in h[0].header['OBJECT'].lower():
                    if 'blue' in h[0].header['OBJECT'].lower():
                        flatsB.append(archivo)
                    else:
                        flatsR.append(archivo)
                else:
                    thars.append(archivo)
                    mjd, mjd0 = mjd_fromheader(h)
                    thar_dates.append( mjd )
            elif h[0].header['EXPTYPE'] == 'Bias':
                biases.append(archivo)
            elif h[0].header['EXPTYPE'] == 'Flat':
                if 'milky' in h[0].header['OBJECT'].lower():
                    milkys.append(archivo)
                else:
                    if 'blue' in h[0].header['OBJECT'].lower():
                        flatsB.append(archivo)
                    else:
                        flatsR.append(archivo)


    f.close()

    return biases, milkys, flatsR, flatsB, science, thars, thar_dates, obnames, exptimes

def mjd_fromheader(h):
    """
    return modified Julian date from header
    """
    ih = 0
    if len(h)>1:
        ih = 1
    datetu = h[ih].header['UT-DATE']
    ut     = h[ih].header['UT-TIME']
    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[5:7]),int(datetu[8:]))
    ut = ut/3600.
    mjd_start = mjd + ut/24.0
    secinday = 24*3600.0
    fraction = .5
    texp     = h[ih].header['EXPTIME'] #sec
    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0

def MedianCombine(ImgList,bs,os, bias = 0.):
    """
    Median combine a list of images
    """
    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

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
            d = np.dstack( (d ,OverscanTrim(h.data,bs,os) -bias ))
        return np.median(d,axis=2), ronoise, gain

def MilkyCombine(ImgList,bs,os, bias = 0.):
    """
    Median combine a list of images
    """
    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    h = pyfits.open(ImgList[0])
    if len(h)>1:
        ih = 1
    else:
        ih = 0

    d = h[ih].data
    d = OverscanTrim(d,bs,os)
    d -= bias
    md = signal.medfilt(d,[21,5])
    d /= md.max()
    factor = 1.25
    if (n < 3):
        factor = 1

    ronoise = factor * h[ih].header['ENOISE'] / np.sqrt(n)
    gain    = h[ih].header['EGAIN']


    if (n == 1):
        return d, ronoise, gain
    else:
        for i in range(n-1):
            h = pyfits.open(ImgList[i+1])[ih]
            d2 = OverscanTrim(h.data,bs,os) - bias
            md2 = signal.medfilt(d2,[21,5])
            d = np.dstack( (d, d2/md2.max() ))
            d = np.median(d,axis=2)
    mm = signal.medfilt(d,[21,5])
    d /= mm
    return d, ronoise, gain

def FlatCombine(ImgList,bs,os, bias = 0.):
    """
    Median combine a list of images
    """
    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    h = pyfits.open(ImgList[0])
    if len(h)>1:
        ih = 1
    else:
        ih = 0

    d = h[ih].data
    d = OverscanTrim(d,bs,os)
    d -= bias

    factor = 1.25
    if (n < 3):
        factor = 1

    ronoise = factor * h[ih].header['ENOISE'] / np.sqrt(n)
    gain    = h[ih].header['EGAIN']


    if (n == 1):
        return d, ronoise, gain
    else:
        for i in range(n-1):
            h = pyfits.open(ImgList[i+1])
            d2 = OverscanTrim(h[ih].data,bs,os) - bias
            md2 = signal.medfilt(d,[5,5])
            d = np.dstack( (d ,d2/(md2.max()) - bias ))
            d = np.median(d,axis=2)
    return d, ronoise, gain

def OverscanTrim(d,bs,os):
    """
    Overscan correct and Trim a refurbished image
    """

    bs[2] = d.shape[0]
    if bs[0] == 1:
        ov = d[:,os[1]:os[2]]
        ov = np.median(ov,axis = 1)
        ov = np.tile(ov,os[1])
        ov = np.reshape(ov,[os[1],bs[2]]).T
        newdata  = d[:,:os[1]] - ov
        newdata  = newdata[:bs[1],:]
    return newdata

def make_blue(lista):
    nlista = []
    for img in lista:
        indx = img.rfind('/')
        nlista.append(img[:indx]+'/b'+img[indx+2:])
    return nlista

def get_sky_shift(spec,lbary_ltopo,path):
    lux = 2.99792458E8
    def gauss_sky(params,x):
        med = params[0]
        sig = params[1]
        amp = params[3]
        bas = params[2]
        g = bas + amp * np.exp(-0.5*(x-med)*(x-med)/(sig*sig))
        return g

    def res_gauss_sky(params,g,x):
        return g-gauss_sky(params,x)

    sky_lines = np.loadtxt(path)
    sky_waves = sky_lines[:,0]
    sky_ints = sky_lines[:,1]
    I = np.where((sky_waves<10000) & (sky_ints>.5))[0]
    sky_waves=sky_waves[I]
    sky_ints=sky_ints[I]
    del_vels = []
    spec[0] /= lbary_ltopo
    for order in range(spec.shape[1]):
        I = np.where((sky_waves>spec[0,order,10]) & (sky_waves<spec[0,order,-10]))[0]
        #plot(spec[0,order,:],spec[5,order,:])
        liml,limh,wavr,intr = [],[],[],[]
        iwr = 0
        for il in sky_waves[I]:
            J = np.where((spec[0,order]>il - 1.) & (spec[0,order]<il + 1.))[0]
            liml.append(J[0])
            limh.append(J[-1])
            wavr.append(il)
            intr.append(sky_ints[I][iwr])
            #plot(spec[0,order,J],spec[5,order,J])
            iwr+=1
        iwr = 0

        while iwr < len(wavr):
            med = [wavr[iwr]]
            amp = [intr[iwr]]
            LIML = liml[iwr]
            LIMH = limh[iwr]
            cond = True
            iwr += 1
            if iwr < len(wavr):
                while cond:
                    if liml[iwr] <= LIMH:
                        med.append(wavr[iwr])
                        amp.append(intr[iwr])
                        LIMH = limh[iwr]
                        iwr += 1
                        if iwr >= len(wavr):
                            break
                    else:
                        cond = False
                x = spec[0,order,LIML:LIMH+1]
                y = spec[5,order,LIML:LIMH+1]
            else:
                x = spec[0,order,LIML:LIMH+1]
                y = spec[5,order,LIML:LIMH+1]
            if len(med) == 1:
                mm0 = med[0]
                guess = [mm0,0.1,y.min(),(y-y.min()).max()]
                Gfit = optimize.leastsq(res_gauss_sky,guess,args=(y,x))
                del_vel = lux * (Gfit[0][0] - mm0) / mm0
                del_vels.append(del_vel)
                #print del_vel, med
            elif len(med) == 2 and amp[0]==amp[1]:
                mm0 = .5*(med[0]+med[1])
                guess = [mm0,0.1,y.min(),(y-y.min()).max()]
                Gfit = optimize.leastsq(res_gauss_sky,guess,args=(y,x))
                del_vel = lux * (Gfit[0][0] - mm0) / mm0
                del_vels.append(del_vel)
    del_vels = np.array(del_vels)
    del_sky_ms, rms_sky_ms, ndel_vels = sig_cl_sky(del_vels)
    print('sky_shift =', del_sky_ms, rms_sky_ms, rms_sky_ms / np.sqrt(float(len(ndel_vels))))
    #hist(ndel_vels,20)
    #show()
    return del_sky_ms, rms_sky_ms, rms_sky_ms / np.sqrt(float(len(ndel_vels)))

def sig_cl_sky(vec):
    m = np.median(vec)
    res = vec - m
    sig = np.sqrt( np.sum(res**2) / float(len(vec)-1) )
    cond = True
    I = np.where(np.abs(res)>3*sig)[0]
    if len(I) == 0 and sig < 2000:
        cond = False
    while cond:
        iw = np.argmax(res**2)
        vec = np.delete(vec,iw)
        m = np.median(vec)
        res = vec - m
        sig = np.sqrt( np.sum(res**2) / float(len(vec)-1) )
        I = np.where(np.abs(res)>3*sig)[0]
        if len(I) == 0 and sig < 2000:
            cond = False
    return m,sig,vec
