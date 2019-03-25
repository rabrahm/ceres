from __future__ import print_function
import matplotlib
matplotlib.use("Agg")
from astropy.io import fits as pyfits
import numpy as np
from numpy import median,sqrt,array,exp
import scipy
from scipy import signal,special,optimize,interpolate
import scipy.special as sp
import copy
import glob
import os
from pylab import *
import sys
base = '../'
sys.path.append(base+"utils/GLOBALutils")
import GLOBALutils



def FileClassify(diri, log,mode='so',amps='a'):
    """
    Classifies all files in a directory and writes a night log of science images
    """
    if mode == 'so':
        md = 'star only'
    elif mode == 'ss':
        md = 'star+sky'
    else:
        md = 'Polarimetry'

    if amps == 'a':
        samps = 'a'
    elif amps == 'ab':
        samps = 'a,b'

    # define output lists
    simThAr_sci    = []
    simSky_sci     = []
    biases         = []
    bias_dates     = []
    flat_dates     = []
    flats          = []
    ThAr_ref       = []
    ThAr_ref_dates = []
    mode           = []

    f = open(log,'w')
    bad_files = []
    if os.access(diri+'bad_files.txt',os.F_OK):
        bf = open(diri+'bad_files.txt')
        linesbf = bf.readlines()
        for line in linesbf:
            bad_files.append(diri+line[:-1])
        bf.close()

    all_files = glob.glob(diri+"*fits")

    jj = 0

    for archivo in all_files:
        jj+=1
        dump = False
        for bf in bad_files:
            if archivo == bf:
                dump = True
                break

        if not dump:
            h = pyfits.open(archivo)
            #print archivo, h[0].header['EXPTYPE'], h[0].header['INSTMODE']
            if md in h[0].header['INSTMODE'] and samps == h[0].header['AMPLIST']:
                print(archivo, h[0].header['EXPTYPE'], h[0].header['INSTMODE'], h[0].data.shape, h[0].header['AMPLIST'])
                if h[0].header['EXPTYPE'] == 'OBJECT':
                    simThAr_sci.append(archivo)
                    obname = h[0].header['OBJNAME']
                    ra     = h[0].header['OBJRA']
                    delta  = h[0].header['OBJDEC']
                    airmass= h[0].header['AIRMASS']
                    texp   = h[0].header['EXPTIME']
                    date   = h[0].header['DATE-OBS'] + h[0].header['UTC-OBS']
                    line = "%-15s %10s %10s %8.2f %4.2f %s %8s %s\n" % (obname, ra, delta, texp, airmass, h[0].header['EXPTYPE'], date, archivo)
                    f.write(line)

                elif h[0].header['EXPTYPE'] == 'BIAS':
                    biases.append(archivo)
                    mjd,mjd0 = mjd_fromheader(h)
                    bias_dates.append(mjd)

                elif h[0].header['EXPTYPE'] == 'FLAT':
                    flats.append(archivo)
                    mjd,mjd0 = mjd_fromheader(h)
                    flat_dates.append(mjd)

                elif h[0].header['EXPTYPE'] == 'COMPARISON':
                    mjd, mjd0 = mjd_fromheader(h)
                    ThAr_ref.append(archivo)
                    ThAr_ref_dates.append( mjd )

    f.close()
    biases, bias_dates = np.array(biases), np.array(bias_dates)
    flats, flat_dates  = np.array(flats), np.array(flat_dates)
    IS = np.argsort(bias_dates)
    biases, bias_dates = biases[IS], bias_dates[IS]
    IS = np.argsort(flat_dates)
    flats, flat_dates = flats[IS], flat_dates[IS]

    return biases, flats, ThAr_ref, simThAr_sci, ThAr_ref_dates

def mjd_fromheader(h):
    """
    return modified Julian date from header
    """

    datetu = h[0].header['DATE-OBS']
    timetu = h[0].header['UTC-OBS'].split(':')

    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[:4]),int(datetu[5:7]),int(datetu[8:10]))

    ho = int(timetu[0])
    mi = int(timetu[1])
    se = float(timetu[2])
    ut = float(ho) + float(mi)/60.0 + float(se)/3600.0
    mjd_start = mjd + ut/24.0

    secinday = 24*3600.0
    fraction = 0.5
    texp     = h[0].header['EXPTIME'] #sec

    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0

def MedianCombine(ImgList, zero_bo=False, zero='MasterBias.fits'):
    """
    Median combine a list of images
    """
    if zero_bo:
        BIAS = pyfits.getdata(zero)

    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    h = pyfits.open(ImgList[0])[0]
    d = h.data
    d = OverscanTrim(d)
    if zero_bo:
        d -= BIAS
        d = np.round(d).astype('int')
    factor = 1.25
    if (n < 3):
        factor = 1
    ronoise = factor * h.header['RDNOISEA'] / np.sqrt(n)
    gain    = h.header['GAINA']
    if (n == 1):
        return d, ronoise, gain
    else:
        for i in range(n-1):
            h = pyfits.open(ImgList[i+1])[0]
            ot = OverscanTrim(h.data)
            if zero_bo:
                d = np.dstack((d,np.round((ot - BIAS)).astype('int')))
            else:
                d = np.dstack((d,np.round(ot).astype('int')))
    return np.median(d,axis=2), ronoise, gain

def OverscanTrim(d):
    """
    Overscan correct and Trim a refurbished FEROS image
    """
    os = d[:,2060:]
    overscan = np.median(os,axis=1)[:4608]
    newdata = np.zeros(d[:4608,:2049].shape)
    for i in range(len(overscan)):
        newdata[i,:] = d[i,:2049] - overscan[i]
    return newdata

def get_them(sc,exap,ncoef,mode=1,shift=0.):
    def gauss(params,x):
        amp1 = params[0]
        med1 = params[1]
        sig1 = params[2]
        g1 = amp1 * np.exp(-0.5*((x-med1)/sig1)**2)
        return g1
    def res_gauss(params,g,x):
        return g-gauss(params,x)

    def gauss2(params,x):
        amp1 = params[0]
        amp2 = params[1]
        med1 = params[2]
        med2 = params[3]
        sig1 = params[4]
        sig2 = params[5]
        g1 = amp1 * np.exp(-0.5*((x-med1)/sig1)**2)
        g2 = amp2 * np.exp(-0.5*((x-med2)/sig2)**2)
        return g1 + g2

    def res_gauss2(params,g,x):
        return g-gauss2(params,x)

    def fit_drift(ref,dif,oref):
        coef  = np.polyfit(ref,dif,2)

        residuals = dif - np.polyval(coef,ref)
        rms = np.sqrt(np.var(residuals))
        I = np.where(np.absolute(residuals)>3*rms)[0]
        cond = True
        if len(I)==0:
            cond = False
        while cond:
            im  = np.argmax(np.absolute(residuals))
            dif = np.delete(dif,im)
            ref = np.delete(ref,im)
            coef = np.polyfit(ref,dif,2)
            residuals = dif - np.polyval(coef,ref)
            rms = np.sqrt(np.var(residuals))
            I = np.where(np.absolute(residuals)>3*rms)[0]
            if len(I)==0:
                cond = False

        cdif = np.polyval(coef,oref)
        ref = oref + cdif
        return ref

    sc_or = sc.copy()
    medc = int(.5*sc.shape[1])
    d = np.median(sc[:,medc-exap:medc+exap+1],axis=1)
    #plot(d)
    ejx = np.arange(len(d))

    ref = []
    cents = np.loadtxt('ref_ords.txt').astype('int') - shift
    #plot(d)
    #plot(ejx[cents[:,0]],d[cents[:,0]],'ro')
    #plot(ejx[cents[:,1]],d[cents[:,1]],'go')
    #show()

    posO, posC = cents[:,0],cents[:,1]
    if mode == 0:
        posO = cents[:,2]
        posC = cents[:,2]
    refO,refC = [],[]
    exap2 = exap + .2*exap
    dev   = exap2/4.
    for i in range(len(posO)):
        if posO[i]-exap2 < 0:
            x = ejx[:posC[i]+exap2+1]
            y = d[:posC[i]+exap2+1]
        elif posC[i]+exap2+1 > len(d):
            x = ejx[posO[i]-exap2:]
            y = d[posO[i]-exap2:]
        else:
            x = ejx[posO[i]-exap2:posC[i]+exap2+1]
            y = d[posO[i]-exap2:posC[i]+exap2+1]
        y -= y.min()
        tx1 = np.arange(x[0]-dev,x[0],1)
        tx2 = np.arange(x[-1]+1,x[-1]+dev+1,1)
        ty1 = np.zeros(len(tx1))
        ty2 = np.zeros(len(tx2))
        x = np.hstack((tx1,x,tx2))
        y = np.hstack((ty1,y,ty2))

        midi = int(0.5*len(x))

        if mode == 1:
            guess = [np.max(y[:midi]),np.max(y[midi:]),x[0]+np.argmax(y[:midi]),x[0]+midi+np.argmax(y[midi:]),dev,dev]
            p, success =  scipy.optimize.leastsq(res_gauss2, guess, args=(y,x))
            refO.append(p[2])
            refC.append(p[3])
        else:
            #plot(x,y)
            #show()
            #print fcds
            guess = [np.max(y),x[0]+np.argmax(y),dev]
            p, success =  scipy.optimize.leastsq(res_gauss, guess, args=(y,x))
            refO.append(p[1])
            refC.append(p[1])

    refO = np.array(refO)
    refOt = np.around(refO).astype('int')
    matO = np.zeros((len(refO),sc.shape[1]))
    matO[:,medc] = refO
    refC = np.array(refC)
    refCt = np.around(refC).astype('int')
    matC = np.zeros((len(refC),sc.shape[1]))
    matC[:,medc] = refC

    #plot(d)
    #plot(cents[:,0],d[cents[:,2]],'ro')
    #plot(cents[:,1],d[cents[:,1]],'go')
    #plot(refOt,d[refOt],'ko')
    #plot(refCt,d[refCt],'yo')
    #show()
    #print gfd5
    i = medc -1
    while i >=0:
        d = sc[:,i]
        #if i<1000:
        #       plot(d)
        j = 0
        posO = np.around(refO).astype('int')
        posC = np.around(refC).astype('int')
        trefO,trefC = [],[]
        exap2 = exap + .2*exap
        dev   = exap2/4.

        while j < len(posO):
            #print j
            if posO[j]-exap2 < 0:
                x = ejx[:posC[j]+exap2+1]
                y = d[:posC[j]+exap2+1]
            elif posC[j]+exap2+1 > len(d):
                x = ejx[posO[j]-exap2:]
                y = d[posO[j]-exap2:]
            else:
                x = ejx[posO[j]-exap2:posC[j]+exap2+1]
                y = d[posO[j]-exap2:posC[j]+exap2+1]
            if len(x) < 7:
                trefO.append(refO[j])
                trefC.append(refC[j])
            else:
                y -= y.min()
                tx1 = np.arange(x[0]-dev,x[0],1)
                tx2 = np.arange(x[-1]+1,x[-1]+dev+1,1)
                ty1 = np.zeros(len(tx1))
                ty2 = np.zeros(len(tx2))
                x = np.hstack((tx1,x,tx2))
                y = np.hstack((ty1,y,ty2))
                midi = int(0.5*len(x))
                if mode == 1:
                    guess = [np.max(y[:midi]),np.max(y[midi:]),x[0]+np.argmax(y[:midi]),x[0]+midi+np.argmax(y[midi:]),dev,dev]
                    p, success =  scipy.optimize.leastsq(res_gauss2, guess, args=(y,x))
                    trefO.append(p[2])
                    trefC.append(p[3])
                else:
                    guess = [np.max(y),x[0]+np.argmax(y),dev]
                    p, success =  scipy.optimize.leastsq(res_gauss, guess, args=(y,x))
                    trefO.append(p[1])
                    trefC.append(p[1])

            j+=1

        orefO = refO.copy()
        trefO = np.array(trefO)
        difO = trefO-refO
        refO = fit_drift(refO,difO,orefO)
        matO[:,i] = refO
        orefC = refC.copy()
        trefC = np.array(trefC)
        difC = trefC-refC
        refC = fit_drift(refC,difC,orefC)
        matC[:,i] = refC

        #if i<500:
        #       plot(d)
        #       I = np.where(refO<2049)[0]
        #       plot(refO[I].astype('int'),d[refO[I].astype('int')],'ko')
        #       show()
        #       print gfd

        i-=4

    refO = matO[:,medc]
    refC = matC[:,medc]

    i = medc +1
    while i < sc.shape[1]:
        d = sc[:,i]
        j = 0
        posO = np.around(refO).astype('int')
        posC = np.around(refC).astype('int')
        trefO,trefC = [],[]
        exap2 = exap + .2*exap
        dev   = exap2/4.

        while j < len(posO):
            #print j
            if posO[j]-exap2 < 0:
                x = ejx[:posC[j]+exap2+1]
                y = d[:posC[j]+exap2+1]
            elif posC[j]+exap2+1 > len(d):
                x = ejx[posO[j]-exap2:]
                y = d[posO[j]-exap2:]
            else:
                x = ejx[posO[j]-exap2:posC[j]+exap2+1]
                y = d[posO[j]-exap2:posC[j]+exap2+1]
            if len(x) < 7:
                trefO.append(refO[j])
                trefC.append(refC[j])
            else:
                y -= y.min()
                tx1 = np.arange(x[0]-dev,x[0],1)
                tx2 = np.arange(x[-1]+1,x[-1]+dev+1,1)
                ty1 = np.zeros(len(tx1))
                ty2 = np.zeros(len(tx2))
                x = np.hstack((tx1,x,tx2))
                y = np.hstack((ty1,y,ty2))
                midi = int(0.5*len(x))
                if mode == 1:
                    guess = [np.max(y[:midi]),np.max(y[midi:]),x[0]+np.argmax(y[:midi]),x[0]+midi+np.argmax(y[midi:]),dev,dev]
                    p, success =  scipy.optimize.leastsq(res_gauss2, guess, args=(y,x))
                    trefO.append(p[2])
                    trefC.append(p[3])
                else:
                    guess = [np.max(y),x[0]+np.argmax(y),dev]
                    p, success =  scipy.optimize.leastsq(res_gauss, guess, args=(y,x))
                    trefO.append(p[1])
                    trefC.append(p[1])
            j+=1

        orefO = refO.copy()
        trefO = np.array(trefO)
        difO = trefO-refO
        refO = fit_drift(refO,difO,orefO)
        matO[:,i] = refO

        orefC = refC.copy()
        trefC = np.array(trefC)
        difC = trefC-refC
        refC = fit_drift(refC,difC,orefC)
        matC[:,i] = refC

        i+=4

    #imshow(sc_or,vmax=5000)
    for i in range(matO.shape[0]):
        y = matO[i]
        x = np.arange(len(y))
        I = np.where(y!=0)[0]
        x,y = x[I],y[I]
        coef = np.polyfit(x,y,ncoef)
        #plot(x,y,'r')
        residuals = y - np.polyval(coef,x)
        rms = np.sqrt(np.var(residuals))
        I = np.where(np.absolute(residuals)>3*rms)[0]
        cond = True
        if len(I)==0:
            cond = False
        while cond:
            im = np.argmax(np.absolute(residuals))
            x = np.delete(x,im)
            y = np.delete(y,im)
            coef = np.polyfit(x,y,ncoef)
            residuals = y - np.polyval(coef,x)
            rms = np.sqrt(np.var(residuals))
            I = np.where(np.absolute(residuals)>3*rms)[0]
            if len(I)==0:
                cond = False
        if i == 0:
            acoefsO = coef
        else:
            acoefsO = np.vstack((acoefsO,coef))
        #plot(np.polyval(coef,np.arange(len(matO[i]))),'g')
    #show()

    if mode == 1:
        for i in range(matC.shape[0]):
            y = matC[i]
            x = np.arange(len(y))
            I = np.where(y!=0)[0]
            x,y = x[I],y[I]
            coef = np.polyfit(x,y,ncoef)
            #plot(x,y,'r')
            residuals = y - np.polyval(coef,x)
            rms = np.sqrt(np.var(residuals))
            I = np.where(np.absolute(residuals)>3*rms)[0]
            cond = True
            if len(I)==0:
                cond = False
            while cond:
                im = np.argmax(np.absolute(residuals))
                x = np.delete(x,im)
                y = np.delete(y,im)
                coef = np.polyfit(x,y,ncoef)
                residuals = y - np.polyval(coef,x)
                rms = np.sqrt(np.var(residuals))
                I = np.where(np.absolute(residuals)>3*rms)[0]
                if len(I)==0:
                    cond = False
            if i == 0:
                acoefsC = coef
            else:
                acoefsC = np.vstack((acoefsC,coef))
            #plot(np.polyval(coef,np.arange(len(matC[i]))),'b')
        #show()
        #print gfds
        #print acoefs.shape
        return acoefsO, acoefsC
    else:
        return acoefsO
