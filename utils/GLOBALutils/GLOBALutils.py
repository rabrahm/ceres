from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

from astropy.io import fits as pyfits
import numpy as np
from numpy import median,sqrt,array,exp
import scipy
from scipy import signal,special,optimize,interpolate,integrate
import scipy.special as sp
import copy
import glob
import os
import matplotlib.pyplot as plt
import sys
import tempfile
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import string
import pycurl
from multiprocessing import Pool
sys.path.append("../utils/CCF")
sys.path.append("../utils/OptExtract")
import Marsh
import CCF
from pylab import *
from PyAstronomy import pyasl
import time

from matplotlib.backends.backend_pdf import PdfPages

global GDATA,P

import statsmodels.api as sm
lowess = sm.nonparametric.lowess

# Some functions to be used by the pipeline
class Constants:
    "Here I declare the constants I will use in the different functions"
    "G,c: the gravity and speed of light constants in SI units; mearth and mmoon: the earth and moon masses in kg"
    G,c = 6.673E-11,2.99792458E8
    mearth, mmoon = 5.9736E24,7.349E22
    "the mass of the planets,moon and sun in earth masses, the ennumeration is sun, moon, mercury,venus,mars,jupiter,saturn, uranus,neptune,pluto"
    mplanets = np.array([332981.787,0.0123000371,0.05528,0.815,0.10745,317.83,95.19,14.536,17.147,0.0021])
    "conversion from degrees to radians, and from hours to degrees, and from degrees to hours"
    degtorad = np.pi/180.0
    radtodeg = 180.0/np.pi
    HtoDeg = 360.0/24.0
    DegtoH = 1.0/HtoDeg
    "Req: equatorial radius f: polar flatenning , w angular speedof earth in rad/s"
    Req = 6378136.6 #in m
    f = 1.0/298.256420
    w = 7.2921158554E-5
    DaystoYear = 1.0/365.256363004

def update_header(hdu,k,v,c=''):
    try:
        hdu.header.update(k,v)
    except:
        hdu.header[k] = v
    return hdu


def get_dark_times(darks,key='EXPTIME'):
    """
    Given a list containing the names of all the dark frames,
    this function returns an array with the different exposure times
    """
    times = []
    for dark in darks:
        hd = pyfits.getheader(dark)
        times.append(hd[key])
    times = np.sort(np.unique(np.array(times)))
    return times

def get_tdarks(darks,dtime,key='EXPTIME'):
    tdarks = []
    for dark in darks:
        hd = pyfits.getheader(dark)
        if hd[key] == dtime:
            tdarks.append(dark)
    return tdarks

#wavelength functions
def n_Edlen(l):
    sigma = 1e4 / l
    sigma2 = sigma*sigma
    n = 1 + 1e-8 * (8342.13 + 2406030 / (130-sigma2) + 15997/(38.9-sigma2))
    return n

def n_Morton(l):
    sigma = 1e4 / l
    sigma2 = sigma*sigma
    n = 1 + 6.4328e-5 + 2.94981e-2 / (146.-sigma2) + 2.5540e-4/(41.-sigma2)
    return n

def ToAir(l):
    return (l / n_Edlen(l))

def ToVacuum(l):
    cond = 1
    l_prev = l.copy()
    while(cond):
        l_new = n_Edlen(l_prev) * l
        if (max(np.absolute(l_new - l_prev)) < 1e-10): cond = 0
        l_prev = l_new
    return l_prev

def getcoords(obname,mjd,filen='/data/echelle/feros/coords.txt'):
    RA,DEC = 0.,0.
    mjdJ2000 = 2451545.0 - 2400000.5

    try:
        f = open(filen,'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            cos = line.split(',')
            if cos[0]==obname and int(cos[5])==1:
                ra,dec = cos[1],cos[2]
                cos1 = ra.split(':')
                cos2 = dec.split(':')
                RA0 = (float(cos1[0]) + float(cos1[1])/60. + float(cos1[2])/3600.) * 360. / 24.
                if float(cos2[0]) > 0:
                    DEC0 = np.absolute(float(cos2[0])) + float(cos2[1])/60. + float(cos2[2])/3600.
                else:
                    DEC0 = -(np.absolute(float(cos2[0])) + float(cos2[1])/60. + float(cos2[2])/3600.)
                PMRA = float(cos[3])/1000.
                PMDEC = float(cos[4])/1000.

                RA  = RA0 + (PMRA/3600.)*(mjd-mjdJ2000)/365.
                DEC = DEC0 + (PMDEC/3600.)*(mjd-mjdJ2000)/365.
                break
    except:
        print('\t\tWarning! Problem with reference coordinates files.')
    return RA,DEC

def get_them(sc,exap,ncoef,maxords=-1,startfrom=0,nsigmas=10.,mode=1,endat=-1,nc2=2):
    exap = int(exap)
    def fitfunc(p,x):
        ret = p[0] + p[1] * np.exp(-.5*((x-p[2])/p[3])**2)
        return ret
    errfunc = lambda p,y,x: np.ravel( (fitfunc(p,x)-y) )

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
    sc_or = sc.copy()
    if endat == -1:
        sc = sc[startfrom:,:]
    else:
        sc = sc[startfrom:endat,:]

    medc = int(.5*sc.shape[1])
    d = np.median(sc[:,medc-exap:medc+exap+1],axis=1)
    #print sc[:,2000]
    #plot(sc[:,1000])
    #show()
    #print sdf
    ejx = np.arange(len(d))
    ccf=[]
    refw = 3*exap
    sigc = .5*exap
    i = 0
    while i < sc.shape[0]:
        if i-refw < 0:
            x = ejx[:i+refw+1]
            y = d[:i+refw+1]
        elif i + refw +1 > sc.shape[0]:
            x = ejx[i-refw:]
            y = d[i-refw:]
        else:
            x = ejx[i-refw:i+refw+1]
            y = d[i-refw:i+refw+1]

        g = np.exp(-0.5*((x-i)/sigc)**2)
        ccf.append(np.add.reduce(y*g))
        i+=1

    i = 1
    maxs = []
    while i < len(ccf)-2:
        if ccf[i]>ccf[i-1] and ccf[i]>ccf[i+1]:
            maxs.append(i)
        i+=1

    maxs = np.array(maxs)
    ccf = np.array(ccf)

    pos = np.arange(len(ccf))[maxs] #+ refw
    dt = d.copy()
    sp = pos[1] - pos[0] - 2*exap
    vx,vy = [],[]
    #plot(d)
    #plot(ccf)
    #plot(pos,d[pos],'ro')
    #show()
    #print gfd
    exap2 = exap + int(exap*2./3.)
    tbase = np.array([])
    for i in range(len(pos)):
        exs,vec = [],[]
        if i == 0:
            if pos[i] - exap2 - sp < 0:
                exs = np.arange(0, pos[i] - exap2+1,1)
                vec = d[: pos[i] - exap2+1]
            else:
                exs = np.arange(pos[i] - exap2 - sp, pos[i] - exap2+1,1)
                vec = d[pos[i] - exap2 - sp: pos[i] - exap2+1]
        else:
            if pos[i-1] + exap2 < pos[i] - exap2+1:
                #print pos[i-1] + exap2 , pos[i] - exap2+1
                exs = np.arange(pos[i-1] + exap2 , pos[i] - exap2+1,1)
                vec = d[pos[i-1] + exap2 : pos[i] - exap2 + 1]

        if len(exs)>0:
            tbase = np.hstack((tbase,exs))
            vx.append(np.median(exs))
            vy.append(np.median(vec))

    tbase = tbase.astype('int')
    vx,vy = np.array(vx),np.array(vy)
    tck = interpolate.splrep(vx,vy,k=1)
    #print interpolate.splev(np.arange(len(d)),tck)
    #plot(d)
    #plot(vx,vy,'ko')
    #plot(np.arange(len(d)),interpolate.splev(np.arange(len(d)),tck))
    #show()
    #print gfds
    dtemp = d[tbase] - interpolate.splev(tbase,tck)
    ddev = np.sqrt(np.var(dtemp[5:-5]))
    dt = d-interpolate.splev(np.arange(len(d)),tck)
    #plot(dt)
    #axhline(3*ddev)
    #axhline(dtemp.mean() + nsigmas*ddev)
    #show()
    #print gfds
    I = np.where(dt[pos]>dtemp.mean() + nsigmas*ddev)[0]
    #print vy

    #plot(d)
    #plot(tbase,d[tbase])
    #plot(np.arange(len(d))[pos], d[pos],'bo')
    #show()
    #print gfds
    #plot(np.arange(len(d))[pos], d[pos],'ro')
    pos = pos[I]
    #print len(pos)
    #print gfds
    #plot(np.arange(len(d))[pos], d[pos],'go')
    #show()
    #print jfhedslja
    #axhline(dtemp.mean() + nsigmas*ddev)
    #xlabel('pixel', fontsize=18)
    #ylabel('Flux [ADUs]', fontsize=18)
    #show()
    #print gfd
    #print len(pos)
    if maxords >0:
        pos = pos[::-1]
        pos = pos[:maxords]
        pos = pos[::-1]

    I = np.where(pos < exap)[0]
    pos = np.delete(pos,I)
    I = np.where(pos > sc.shape[0]-exap)[0]
    pos = np.delete(pos,I)
    #print len(pos)
    #plot(d)
    #plot(np.arange(len(d))[pos], d[pos],'ro')
    #show()
    #plot(d)
    ref = []
    if mode == 1 or mode == 2:
        if mode == 1:
            exap2 = exap + .5*exap
            dev = exap2/3.
        else:
            exap2 = exap + .2*exap
            dev = exap2/4.
        exap2 = int(exap2)
        for i in range(len(pos)):
            if pos[i]-exap2 < 0:
                x = ejx[:int(pos[i]+exap2+1)]
                y = d[:int(pos[i]+exap2+1)]
            elif pos[i]+exap2+1 > len(d):
                x = ejx[int(pos[i]-exap2):]
                y = d[int(pos[i]-exap2):]
            else:
                x = ejx[int(pos[i]-exap2):int(pos[i]+exap2+1)]
                y = d[int(pos[i]-exap2):int(pos[i]+exap2+1)]
            tx1 = np.arange(x[0]-dev,x[0],1)
            tx2 = np.arange(x[-1]+1,x[-1]+dev+1,1)
            ty1 = np.zeros(len(tx1))
            ty2 = np.zeros(len(tx2))
            x = np.hstack((tx1,x,tx2))
            y = np.hstack((ty1,y,ty2))
            y -= y.min()

            if mode == 1:
                if len(x) < 4:
                    tref.append(ref[j])
                else:
                    p, success =  scipy.optimize.leastsq(errfunc, [y.min(),y.max()-y.min(),x.mean(),dev], args=(y,x))
                ref.append(p[2])
            else:
                midi = int(0.5*len(x))
                if len(x) < 7:
                    tref.append(ref[j])
                else:
                    guess = [np.max(y[:midi]),np.max(y[midi:]),x[0]+np.argmax(y[:midi]),x[0]+midi+np.argmax(y[midi:]),dev,dev]
                    p, success =  scipy.optimize.leastsq(res_gauss2, guess, args=(y,x))
                    #print guess
                    #print p
                    #plot(x,y)
                    #axvline(p[2])
                    #axvline(p[3])
                    #axvline(0.5*(p[2]+p[3]))
                    #show()
                ref.append(0.5*(p[2]+p[3]))

    ref = np.array(ref)
    #plot(d)
    #plot(np.arange(len(d))[np.around(ref).astype('int')], d[np.around(ref).astype('int')],'ro')
    #show()
    mat = np.zeros((len(ref),sc.shape[1]))
    mat[:,medc] = ref
    i = medc -1

    while i >=0:
        #print i
        d = sc[:,i]
        j = 0
        pos = np.around(ref).astype('int')
        tref = []
        if mode == 1 or mode == 2:
            if mode == 1:
                exap2 = exap + .5*exap
                dev = exap2/3.
            else:
                exap2 = exap + .2*exap
                dev = exap2/4.
            exap2 = int(exap2)
            while j < len(pos):
                if pos[j]-exap2 < 0:
                    x = ejx[:int(pos[j]+exap2+1)]
                    y = d[:int(pos[j]+exap2+1)]
                elif pos[j]+exap2+1 > len(d):
                    x = ejx[int(pos[j]-exap2):]
                    y = d[int(pos[j]-exap2):]
                else:
                    x = ejx[int(pos[j]-exap2):int(pos[j]+exap2+1)]
                    y = d[int(pos[j]-exap2):int(pos[j]+exap2+1)]

                if mode==1:
                    if len(x) < 4:
                        tref.append(ref[j])
                    else:
                        tx1 = np.arange(x[0]-dev,x[0],1)
                        tx2 = np.arange(x[-1]+1,x[-1]+dev+1,1)
                        ty1 = np.zeros(len(tx1)) + y.min()
                        ty2 = np.zeros(len(tx2)) + y.min()
                        x = np.hstack((tx1,x,tx2))
                        y = np.hstack((ty1,y,ty2))
                        p, success =  scipy.optimize.leastsq(errfunc, [y.min(),y.max()-y.min(),x.mean(),dev], args=(y,x))
                        #plot(x,y)
                        #plot(x,fitfunc(p,x))
                        tref.append(p[2])
                else:
                    if len(x) < 7:
                        tref.append(ref[j])
                    else:
                        tx1 = np.arange(x[0]-dev,x[0],1)
                        tx2 = np.arange(x[-1]+1,x[-1]+dev+1,1)
                        ty1 = np.zeros(len(tx1)) + y.min()
                        ty2 = np.zeros(len(tx2)) + y.min()
                        x = np.hstack((tx1,x,tx2))
                        y = np.hstack((ty1,y,ty2))
                        y -= y.min()
                        midi = int(0.5*len(x))
                        guess = [np.max(y[:midi]),np.max(y[midi:]),x[0]+np.argmax(y[:midi]),x[0]+midi+np.argmax(y[midi:]),dev,dev]
                        p, success =  scipy.optimize.leastsq(res_gauss2, guess, args=(y,x))
                        tref.append(0.5*(p[2]+p[3]))
                j+=1
        #show()
        oref = ref.copy()
        tref = np.array(tref)
        dif = tref-ref
        coef = np.polyfit(ref,dif,nc2)
        #plot(ref,dif,'ro')
        #plot(oref,np.polyval(coef,oref))
        coef2 = np.polyfit(np.arange(len(dif)),dif,1)
        #if i < 500:
        #       plot(np.arange(len(dif)),dif,'.')
        #       plot(np.arange(len(dif)),np.polyval(coef2,np.arange(len(dif))))
        #       xlabel('echelle order number',fontsize=18)
        #       ylabel('shift [px]',fontsize=18)
        #       show()
        #       print gfds
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
            coef = np.polyfit(ref,dif,nc2)
            residuals = dif - np.polyval(coef,ref)
            rms = np.sqrt(np.var(residuals))
            I = np.where(np.absolute(residuals)>3*rms)[0]
            if len(I)==0:
                cond = False
        #plot(oref,np.polyval(coef,oref))
        cdif = np.polyval(coef,oref)
        ref = oref + cdif
        #plot(ref,np.polyval(coef,ref))
        #show()

        mat[:,i] = ref
        i-=4

    i = medc+1
    ref = mat[:,medc]
    while i < sc.shape[1]:
        #print i
        d = sc[:,i]
        j = 0
        pos = np.around(ref).astype('int')
        tref = []
        if mode == 1 or mode == 2:
            if mode == 1:
                exap2 = exap + .5*exap
                dev = exap2/3.
            else:
                exap2 = exap + .2*exap
                dev = exap2/4.
            while j < len(pos):
                if pos[j]-exap2 < 0:
                    x = ejx[:int(pos[j]+exap2+1)]
                    y = d[:int(pos[j]+exap2+1)]
                elif pos[j]+exap2+1 > len(d):
                    x = ejx[int(pos[j]-exap2):]
                    y = d[int(pos[j]-exap2):]
                else:
                    x = ejx[int(pos[j]-exap2):int(pos[j]+exap2+1)]
                    y = d[int(pos[j]-exap2):int(pos[j]+exap2+1)]

                if mode == 1:
                    if len(x) < 4:
                        tref.append(ref[j])
                    else:
                        tx1 = np.arange(x[0]-dev,x[0],1)
                        tx2 = np.arange(x[-1]+1,x[-1]+dev+1,1)
                        ty1 = np.zeros(len(tx1)) + y.min()
                        ty2 = np.zeros(len(tx2)) + y.min()
                        x = np.hstack((tx1,x,tx2))
                        y = np.hstack((ty1,y,ty2))
                        p, success =  scipy.optimize.leastsq(errfunc, [y.min(),y.max()-y.min(),x.mean(),dev], args=(y,x))
                        tref.append(p[2])
                else:
                    if len(x) < 7:
                        tref.append(ref[j])
                    else:
                        tx1 = np.arange(x[0]-dev,x[0],1)
                        tx2 = np.arange(x[-1]+1,x[-1]+dev+1,1)
                        ty1 = np.zeros(len(tx1)) + y.min()
                        ty2 = np.zeros(len(tx2)) + y.min()
                        x = np.hstack((tx1,x,tx2))
                        y = np.hstack((ty1,y,ty2))
                        y -= y.min()
                        midi = int(0.5*len(x))
                        guess = [np.max(y[:midi]),np.max(y[midi:]),x[0]+np.argmax(y[:midi]),x[0]+midi+np.argmax(y[midi:]),dev,dev]
                        p, success =  scipy.optimize.leastsq(res_gauss2, guess, args=(y,x))
                        tref.append(0.5*(p[2]+p[3]))

                j+=1

        oref = ref.copy()
        tref = np.array(tref)
        dif = tref-ref

        coef = np.polyfit(ref,dif,nc2)

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
            coef = np.polyfit(ref,dif,nc2)
            residuals = dif - np.polyval(coef,ref)
            rms = np.sqrt(np.var(residuals))
            I = np.where(np.absolute(residuals)>3*rms)[0]
            if len(I)==0:
                cond = False

        cdif = np.polyval(coef,oref)
        ref = oref + cdif
        #plot(ref,np.polyval(coef,ref))
        #show()
        mat[:,i] = ref
        i+=4

    #imshow(sc_or,vmax=5000)
    for i in range(mat.shape[0]):
        y = mat[i]
        x = np.arange(len(y))
        I = np.where(y!=0)[0]
        x,y = x[I],y[I]+startfrom
        coef = np.polyfit(x,y,ncoef)
        #plot(x,y,'bo')
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
        #coef[-1] += startfrom
        if i == 0:
            acoefs = coef
        else:
            acoefs = np.vstack((acoefs,coef))
        #plot(np.polyval(coef,np.arange(len(mat[i]))),'k')

    #show()
    #print gfds
    #print acoefs.shape
    return acoefs, len(acoefs)

def good_orders(coef,nord,ny,nx,ext_aperture):
    Centers = np.zeros((nord,nx))
    ejx = np.arange(nx)
    bad_inx = []
    for i in range(nord):
        Centers[i,:]=scipy.polyval(coef[i],ejx)
        I = np.where(Centers[i,:]+ext_aperture>ny)[0]
        if len(I)>0:
            bad_inx.append(i)
    bad_inx = np.array(bad_inx)
    if len(bad_inx) > 0:
        im = np.min(bad_inx)
        return coef[:im], nord-len(bad_inx)
    else:
        return coef, nord

def get_zero_order_number(ords,wavs):
    ords,wvas = np.array(ords),np.array(wavs)
    i = 0
    o0 = 0
    pends = []
    o0s = np.arange(200)
    figure(1)
    subplot(211)
    for i in o0s:
        orders = ords + i
        val = orders*wavs
        val2 = val/np.add.reduce(val)
        coef = np.polyfit(ords,val2,1)
        print(coef)
        pends.append(coef[0])
        plot(ords,val2)
    pends = np.array(pends)
    I = np.argmin(pends**2)
    print(o0s[I])
    xlabel('raw order number')
    ylabel(' $\lambda_c$*(order number + zero order)')
    subplot(212)
    plot(o0s,np.absolute(pends))
    xlabel('zero order')
    ylabel('slope')
    show()
    print(gfdsa)

def Mesh(ycpoly_ob,ycpoly_co):
    """
    Mesh trace products into one single array.
    """
    ycpoly = np.vstack((ycpoly_ob,ycpoly_co))
    lastcol = ycpoly[:,-1]
    I = np.argsort(lastcol)
    for i in range(len(I)):
        if i == 0:
            npoly = ycpoly[I[i]]
        else:
            npoly = np.vstack((npoly,ycpoly[I[i]]))
    return npoly

def get_mask(sc,coefs,spa):
    for i in range(len(coefs)):
        vec = np.polyval(coefs[i],np.arange(sc.shape[1])).astype('int')
        if i == 0:
            traces = vec.copy()
        else:
            traces = np.vstack((traces,vec))

    mask = np.zeros(sc.shape)
    for i in range(traces.shape[1]):
        for j in range(traces.shape[0]):
            mask[traces[j,i]-spa:traces[j,i]+spa+1,i] = 1
    return mask

def get_scat(sc,lim,span = 7, typ='median', allow_neg=False,option=0):
    scat = np.zeros(sc.shape)
    ejeX = np.arange(sc.shape[0])

    for y in range(sc.shape[1]):
        lims = np.around(lim[:,y]).astype('int')
        nejX,nejY = np.array([]),np.array([])
        #plot(sc[:,y])
        for j in range(len(lims)):
            if j == 0:
                #print lims[j] - span
                if lims[j] - span < 0:
                    ejx, ejy = [],[]
                elif lims[j] - 2 * span < 0:
                    ejx=ejeX[:lims[j]-span]
                    ejy=sc[:lims[j]-span,y]
                else:
                    ejx=ejeX[lims[j]- 2 * span:lims[j]-span+1]
                    ejy=sc[lims[j]- 2 * span:lims[j]-span+1,y]


            else:
                if lims[j-1] + span >= sc.shape[0] or lims[j-1] + span < 0:
                    ejx,ejy = [],[]

                elif lims[j] - span + 1 > sc.shape[0]:
                    ejx=ejeX[lims[j-1] + span: ]
                    ejy=sc[lims[j-1] + span:, y]
                elif lims[j-1] + span >= lims[j]- span + 1:
                    ejx,ejy = [],[]
                else:
                    ejx=ejeX[lims[j-1] + span:lims[j]- span + 1 ]
                    ejy=sc[lims[j-1] + span:lims[j]- span + 1, y]

                if option == 1 and len(ejx) == 0:
                    tpos = int(np.around(0.5*(lims[j-1] + lims[j])))
                    if tpos >= 0 and tpos < sc.shape[0]:
                        ejx = np.array([ejeX[tpos]])
                        ejy = np.array([sc[tpos,y]])

            #plot(sc[:,y])
            #plot(lims[j],sc[lims[j],y],'ro')
            #plot(ejx,ejy)
            #show()
            #print fd

            if len(ejy)>0:
                if typ== 'median':
                    value = np.median(ejy)
                elif typ == 'min':
                    value = np.min(ejy)

                if np.isnan(value)==True or np.isnan(-value)==True:
                    value = 0.
                if value < 0 and (not allow_neg):
                    value = 0.

            if len(ejx) > 0:
                if len(nejX) == 0:
                    nejX = np.hstack((nejX,np.median(ejx)))
                    nejY = np.hstack((nejY,value))
                elif np.median(ejx) > nejX[-1]:
                    nejX = np.hstack((nejX,np.median(ejx)))
                    nejY = np.hstack((nejY,value))
                if j == 1 and len(nejY)>1:
                    nejY[0] = nejY[1]

        if lims[-1]+span >= sc.shape[0]:
            ejx,ejy = [],[]
        elif lims[-1]+2*span > sc.shape[0]:
            ejx,ejy = ejeX[lims[-1]+span:],sc[lims[-1]+span:,y]
        else:
            ejx=ejeX[lims[-1]+span:lims[-1]+2*span]
            ejy=sc[lims[-1]+span:lims[-1]+2*span,y]

        if len(ejx)>0:
            value = np.median(ejy)
            if value < 0 or np.isnan(value)==True or np.isnan(-value)==True:
                value = 0.
            nejX = np.hstack((nejX,np.median(ejx)))
            nejY = np.hstack((nejY,value))
        tck = interpolate.splrep(nejX,nejY,k=1)

        scat[:lims[-1]+2*span,y] = interpolate.splev(ejeX,tck)[:lims[-1]+2*span]

    #plot(np.arange(sc.shape[0])[lim[:,1000].astype('int')],sc[lim[:,1000].astype('int'),1000],'ro')
    #plot(sc[:,1000])
    #plot(scat[:,1000])
    #show()
    #plot(scat[1000])
    scat = scipy.signal.medfilt(scat,[5,15])
    #plot(scat[1000])
    #show()
    return scat

def get_scat2(sc,ps):
    ejex = np.arange(sc.shape[0])
    bac = sc.copy()
    #print sc.shape
    #print ps.shape
    for col in range(sc.shape[1]):
        #print col
        #col = 2000
        J = np.where(ps[:,col] == 0)[0]
        Z = np.zeros((sc.shape[0]))
        Z[J] = 1.
        #plot(sc[:,col])
        #plot(sc[:,col]*Z,linewidth=2.0)
        #show()
        #print gfd
        Z2 = np.append(Z[-1],Z[:-1])
        I = np.where((Z!=0) & (Z2==0))[0]
        J = np.where((Z2!=0) & (Z==0))[0]
        J = J[1:]
        I = I[:-1]
        points = []
        vals = []
        for i in range(len(I)):
            #plt.plot(np.mean(ejex[I[i]:J[i]]),np.median(sc[I[i]:J[i],col]),'ro')
            points.append(np.mean(ejex[I[i]:J[i]]))
            vals.append(np.median(sc[I[i]:J[i],col]))
        #plt.show()

        if len(points)>2:
            points,vals = np.array(points),np.array(vals)
            F = np.where(vals < 10000)[0]
            vals = vals[F]
            points = points[F]
            tck = interpolate.splrep(points,vals,k=1)
            scat = interpolate.splev(ejex,tck)
            scat[:I[0]] = scat[I[0]]
            scat[J[-1]:] = scat[J[-1]]
            bac[:,col] = scat
            #plot(ejex,scat,'r')
            #show()
        else:
            bac[:,col] = 0.

    bacm = signal.medfilt2d(bac,[1,21])
    #plot(bac[1000])
    #plot(bacm[1000])
    #show()
    return bacm

def clean(x,y):
    lI = len(x)
    coef = np.polyfit(x,y,1)
    res = y - np.polyval(coef,x)
    dev = np.sqrt( np.sum(res**2) / float(len(res) - 2) )
    I = np.where(np.abs(res) < 2*dev)[0]
    cond = True
    if len(I) == lI:
        cond = False

    while cond:
        iw = np.argmax(res**2)
        x = np.delete(x,iw)
        y = np.delete(y,iw)
        coef = np.polyfit(x,y,1)
        res = y - np.polyval(coef,x)
        dev = np.sqrt( np.sum(res**2) / float(len(res) - 2) )
        #plot(x,y,'.')
        #plot(x,np.polyval(coef,x))
        #plot(x,np.polyval(coef,x)+2*dev)
        #plot(x,np.polyval(coef,x)-2*dev)
        #show()
        I = np.where(np.abs(res) < 2*dev)[0]
        J = np.where(np.abs(res) >= 2*dev)[0]
        if len(J) == 0:
            cond = False
    return x,y
def sig_cli(v,ns=2.):
    vm = np.median(v)
    res = v - vm
    dev = np.sqrt(np.sum(res**2) / float(len(res)-1))
    I = np.where(np.abs(res)>ns*dev)[0]
    cond = True
    if len(I)==0:
        cond = False
    while cond:
        iw = np.argmax(res**2)
        v = np.delete(v,iw)
        vm = np.median(v)
        res = v - vm
        dev = np.sqrt(np.sum(res**2) / float(len(res)-1))
        I = np.where(np.absolute(res)>ns*dev)[0]
        if len(I)==0:
            cond = False
    return vm,dev

def sig_cli2(v,ns=3.):
    vm = np.median(v)
    res = v - vm
    dev = np.sqrt(np.sum(res**2) / float(len(res)-1))
    I = np.where(np.abs(res)>ns*dev)[0]
    cond = True
    if len(I)==0:
        cond = False
    while cond:
        iw = np.argmax(res**2)
        v = np.delete(v,iw)
        vm = np.median(v)
        res = v - vm
        dev = np.sqrt(np.sum(res**2) / float(len(res)-1))
        I = np.where(np.absolute(res)>ns*dev)[0]
        if len(I)==0:
            cond = False
    return v

def fit(x,y,n):
    xo = x.copy()
    coef = np.polyfit(x,y,n)
    res = y - np.polyval(coef,x)
    dev = np.sqrt(np.sum(res**2)/float(len(x)-n-1))
    I = np.where(np.absolute(res)>3*dev)[0]
    cond = True
    if len(I) == 0:
        cond = False
    while cond:
        iw = np.argmax(res**2)
        x = np.delete(x,iw)
        y = np.delete(y,iw)
        coef = np.polyfit(x,y,n)
        res = y - np.polyval(coef,x)
        dev = np.sqrt(np.sum(res**2)/float(len(x)-n-1))
        I = np.where(np.absolute(res)>3*dev)[0]
        if len(I) == 0:
            cond = False
    return np.around(np.polyval(coef,xo)).astype('int')

def bkg_flat(data,lim,span = 7):

    bac = data.copy()
    ejx = np.arange(data.shape[0])
    for i in range(data.shape[1]):
        #print i
        #i = 4098
        line = data[:,i]
        #plot(line)
        #show()
        tline = line.copy()
        lims = lim[:,i].astype('int')
        dys = []
        L1s,L2s=[],[]
        for j in range(len(lims)):
            #print j
            x = ejx[lims[j]-span:lims[j]+span+1]
            y = line[lims[j]-span:lims[j]+span+1]
            my,dy = sig_cli(y)
            #dy = np.sqrt(np.var(y))
            I = np.where(line[:lims[j]] < my-3*dy)[0]
            if len(I)==0:
                L1 = 0
            else:
                L1 = I[-1]
            I = np.where(line[lims[j]:] < my-3*dy)[0]
            if len(I)==0:
                L2 = 0
            else:
                L2 = lims[j]+I[0]
            L1s.append(L1)
            L2s.append(L2)

            #if j == 0:
            #       tline[:L1 - span] = 0
            #if j == len(lims)-1:
            #       tline[L2+span:] = 0
            #print L1,L2
            #plot(ejx[L1:L2+1],line[L1:L2+1])
            #dys.append(my)
            #plot(x,y,'r')
            #plot(x,np.zeros(len(x))+ my + 3*dy,'k')
            #plot(x,np.zeros(len(x))+ my - 3*dy,'k')
        #show()
        L1s,L2s = np.array(L1s),np.array(L2s)
        ref = np.arange(len(L1s))
        #print L1s,L2s
        L1s = fit(ref,L1s,3)
        L2s = fit(ref,L2s,3)
        #print L1s,L2s
        for j in range(len(lims)):
            tline[L1s[j]:L2s[j]+1] = 0
        tline[:L1s[0] - span] = 0
        tline[L2s[-1] + span:] = 0
        I = np.where(tline!=0)[0]
        tx = ejx[I]
        ty = tline[I]
        tx2 = np.append(tx[-1],tx[:-1])
        I = np.where(tx2 +1 != tx )[0]
        vx,vy = [],[]
        for j in range(len(I)):
            if j != len(I)-1:
                vy.append(np.median(ty[I[j]:I[j+1]]))
                vx.append(np.median(tx[I[j]:I[j+1]]))
            else:
                vy.append(np.median(ty[I[j]:]))
                vx.append(np.median(tx[I[j]:]))
        tck = interpolate.splrep(vx,vy,k=1)
        bkg = interpolate.splev(ejx,tck)
        bkg[:L1s[0] - span] = bkg[L1s[0] - span]
        bkg[L2s[-1] + span:] = bkg[L2s[-1] + span]
        bac[:,i] = bkg
        #plot(data[:,i])
        #plot(bac[:,i])
        #show()
    span1 = (span/2)*2 + 1
    bac = signal.medfilt(bac,[span1,span*2+1])
    return bac

def MedianCombine(ImgList,ZF=0.):
    """

    Median combine a list of images

    """

    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    h = pyfits.open(ImgList[0])[0]
    d = h.data
    d = OverscanTrim(d) - ZP

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
            d = np.dstack((d,OverscanTrim(h.data)-ZP))
        return np.median(d,axis=2), ronoise, gain

def MedianCombine_simple(ImgList,ZF=0.):
    """

    Median combine a list of images

    """

    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    h = pyfits.open(ImgList[0])[0]
    d = h.data - ZF

    factor = 1.25
    if (n < 3):
        factor = 1

    #ronoise = factor * h.header['HIERARCH ESO CORA CCD RON'] / np.sqrt(n)
    #gain    = h.header['HIERARCH ESO CORA CCD GAIN']

    if (n == 1):
        return d
    else:
        for i in range(n-1):
            h = pyfits.open(ImgList[i+1])[0]
            dtemp = h.data-ZF
            #plot(dtemp[:,1000])
            d = np.dstack((d,h.data-ZF))
        #show()
        return np.median(d,axis=2)


def FlatNormalize(S_flat_ob, S_flat_ob_simple, mid=1023):
    """
    Given flat-field spectra return normalized version and normalized+smooth version
    """
    norders = S_flat_ob_simple.shape[0]
    span = 500
    S_flat_ob_n  = np.zeros( np.shape(S_flat_ob) )
    S_flat_ob_simple_n  = np.zeros( np.shape(S_flat_ob_simple) )

    for j in np.arange(norders).astype('int'):
        max_val        = (scipy.signal.medfilt(S_flat_ob[j,1,int(np.around(mid-span)):int(np.around(mid+span+1))],21)).max()
        S_flat_ob_n[j,1,:]     = (S_flat_ob[j,1,:] / max_val)
        S_flat_ob_n[j,2,:]  = S_flat_ob[j,1,:] * max_val**2

        max_val        = (scipy.signal.medfilt(S_flat_ob_simple[j,int(np.around(mid-span)):int(np.around(mid+span+1))],21)).max()
        S_flat_ob_simple_n[j,:]     = (S_flat_ob_simple[j,:] / max_val)

        # fill the trivial parts now
        S_flat_ob_n[j,0,:]         = S_flat_ob[j,0,:]

    return S_flat_ob_n,S_flat_ob_simple_n

def FlatNormalize_single(S_flat,mid=1023,span=500):
    """
    Given flat-field spectra return normalized version and normalized+smooth version
    """
    S_flat_n  = np.zeros( np.shape(S_flat) )
    norders = S_flat.shape[0]
    max_vals = []
    if len(S_flat.shape)==2:
        for j in np.arange(norders):
            max_val = (scipy.signal.medfilt(S_flat[j,int(mid-span):int(mid+span)+1],21)).max()
            S_flat_n[j,:]     = S_flat[j,:] / max_val
            max_vals.append(max_val)
    else:
        for j in np.arange(norders):
            max_val = (scipy.signal.medfilt(S_flat[j,1,int(mid-span):int(mid+span+1)],21)).max()
            S_flat_n[j,1,:] = S_flat[j,1,:] / max_val
            S_flat_n[j,2,:] = S_flat[j,2,:] * max_val**2
            S_flat_n[j,0,:] = S_flat[j,0,:]
            max_vals.append(max_val)
    return S_flat_n,np.array(max_vals)

def retrace(dat, c_all,span=9):
    def gauss(params,x):
        med = params[0]
        sig = params[1]
        g = np.exp(-0.5*(x-med)*(x-med)/(sig*sig))
        return g

    def res_gauss(params,g,x):
        return g-gauss(params,x)

    #for i in range(len(c_all)):
    #       px = int(np.around(np.polyval(c_all[i],2000)))
    #       plot(dat[px-7:px+8,2000])
    #show()

    med = int(.5*dat.shape[1])
    #print dat.shape
    Centers = np.zeros((len(c_all),dat.shape[1]))
    Cen = np.zeros((len(c_all),dat.shape[1]))
    for i in range(c_all.shape[0]):
        Centers[i,:] = np.around(scipy.polyval(c_all[i,:],np.arange(dat.shape[1])))
        Cen [i,:]    = scipy.polyval(c_all[i,:],np.arange(dat.shape[1]))
    i = -span
    CCF = []
    pix = []
    #plot(dat[:,1000])
    #show()
    while i < span+1:

        mat = np.zeros(dat.shape)
        for j in range(dat.shape[1]):
            vec = Centers[:,j]
            vec = vec.astype('int')+i
            IN = np.where((vec>0) & (vec<dat.shape[0]))[0]
            mat[vec[IN],j] = 1.0
        #print np.add.reduce((mat*dat))[med - 100 : med + 101]
        vvv = np.add.reduce((mat*dat))[med - 100 : med + 101]
        II = np.where(np.isnan(vvv)==False)[0]
        vvv = vvv[II]
        CCF.append(np.add.reduce(vvv))
        pix.append(i)


        i+=1
    pix = np.array(pix)
    CCF = np.array(CCF)

    CCF -= np.min(CCF)
    CCF /= np.max(CCF)
    #print CCF
    #plot(pix,CCF)
    #show()
    im = np.where(CCF == CCF.max())[0][0]

    jj = im
    mi=0
    while jj > 0:
        if CCF[jj-1] >= CCF[jj] and  CCF[jj+1] >=CCF[jj]:
            mi = jj
            break
        jj-=1
    jj = im
    ma = len(CCF)-1
    while jj < len(CCF)-1:

        if CCF[jj-1] >= CCF[jj] and  CCF[jj+1]>= CCF[jj]:
            ma = jj
            break
        jj+=1

    if pix[im]<0:
        guess = [pix[im],3.0]
    else:
        guess = [pix[im],3.0]
    #plot(pix,CCF)
    #show()
    Gfit = optimize.leastsq(res_gauss,guess,args=(CCF[mi:ma+1],pix[mi:ma+1]))
    shift = Gfit[0][0]

    #print 'yshift:', shift

    c_new = c_all.copy()
    for j in range(c_all.shape[0]):
        c_new[j] = np.polyfit( np.arange(dat.shape[1]), Cen[j]+shift, c_all.shape[1]-1 )

    return c_new,int(round(shift))

def get_drift(data,P,c_all,pii,win=10):
    d = np.median(data[:,pii-win:pii+win+1],axis=1)
    p = np.median(P[:,pii-win:pii+win+1],axis=1)
    d -= d.min()

    for i in range(len(c_all)):
        j = int(np.around(np.polyval(c_all[i],pii)))
        p[j-win:j+win+1] /= p[j-win:j+win+1].max()
        p[j-win:j+win+1] *= d[j-win:j+win+1].max()
    ejxo = np.arange(len(p))
    ps = -4
    dels,ccf = [],[]
    while ps <=4:
        ejxt = ejxo + ps
        tck = scipy.interpolate.splrep(ejxt,p,k=3)
        tp = scipy.interpolate.splev(ejxo,tck)
        ccf.append(np.add.reduce(d[50:-100]*tp[50:-100]))
        dels.append(ps)
        ps+=0.1
    dels,ccf = np.array(dels),np.array(ccf)
    tck = scipy.interpolate.splrep(dels,ccf,k=3)
    ndels = np.arange(dels[0],dels[-1],0.01)
    nccf = scipy.interpolate.splev(ndels,tck)
    shift = ndels[np.argmax(nccf)]

    c_new = c_all.copy()
    c_new[:,-1] += shift
    return shift, c_new

def invert(spec):
    if len(spec.shape)==2:
        for i in range(spec.shape[0]):
            spec[i] = spec[i][::-1]
    else:
        for i in range(spec.shape[0]):
            spec[i,0] = spec[i,0][::-1]
            spec[i,1] = spec[i,1][::-1]
            spec[i,2] = spec[i,2][::-1]
    return spec

def shift_P(P,shift,c_new,ap):
    nP = np.zeros(P.shape)
    for i in range(P.shape[1]):
        d = P [:,i]
        ej = np.arange(len(d))+shift
        ej2 = np.arange(len(d))
        ej, ej2 = ej.astype('float'),ej2.astype('float')
        tck = scipy.interpolate.splrep(ej,d,k=3)

        nP[:,i] = scipy.interpolate.splev(ej2,tck)
        I = np.where(nP[:,i]< 0)[0]
        nP[I,i]=0
        inv = np.array([])
        for j in range(len(c_new)):
            pi = np.polyval(c_new[j],i)
            pi = int(np.around(pi))
            #print pi
            tv = nP[pi - ap:pi+ap+1,i]
            tv /= np.sum(tv)
            nP[pi - ap:pi+ap+1,i] = tv
            inv = np.hstack( (inv,np.arange(pi-ap,pi+ap+1,1)) )

        inv = inv.astype('int')
        I = np.where(inv>P.shape[0]-1)[0]
        if len(I)>0:
            inv = np.delete(inv,I)
        I = np.where(inv<0)[0]
        if len(I)>0:
            inv = np.delete(inv,I)
        niv = np.zeros(P.shape[0]).astype('int')
        niv[inv] = 1.0
        nP[:,i] *= niv
    return nP

#geographical functions
def obspos(longitude,obsradius,R0):
    """
    Set the observatory position respect to geocenter in the coordinates(x,y,z)required by jplepem,
    x to equator/greenwich intersection,
    y 90 degrees east,
    z positive to north
    """
    obpos = []
    x = obsradius*np.cos( (np.pi / 180.0) * longitude )
    obpos.append(x)
    y = obsradius*np.sin( (np.pi / 180.0) * longitude )
    obpos.append(y)
    z = R0
    obpos.append(z)
    return obpos

def JPLR0(lat, altitude):
    "the corrections due to earth rotation"
    "this function returns the velocity in m/s, the projected distance of the observatory in the equator axis and in the earth spin axis and \
    also returns the distance from the rotation axis, and the distance from the equator plane to the observatory position"
    "The arguments are latitude,altitude and hour angle at observatory , dec: the declination of the star\
    the variables are: \
    Raxis: the distance from the observatory to the earth axis, and Psid: the period of sidereal day"
    lat = Constants.degtorad*lat
    e2 = Constants.f*(2.0-Constants.f)
    c1 = 1.0 - e2*(2.0 - e2)*np.sin(lat)**2
    c2 = 1.0 - e2*np.sin(lat)**2

    #radius at 0 elevation
    R0 = Constants.Req*np.sqrt(c1/c2)

    #the geocentric latitude
    c1 = e2*np.sin(2.0*lat)
    c2 = 2.0*c2
    geolat = lat - np.arctan(c1/c2)
    #Calculate geocentric radius at altitude of the observatory
    GeoR = R0*np.cos(geolat) + altitude*np.cos(lat)

    # the R0 vector is now the distance from the observatory to the declination 0 deg plane
    R0 = R0*np.sin(abs(geolat))+altitude*np.sin(lat)
    return GeoR,R0

def JPLiers(path, mjdini, mjdend):
    output    = open(path+'iers.tab','w')
    filename  = path+'finals2000A.data'
    finaldata = open(filename,'r')

    for line in finaldata:
        mj = line[7:15]
        if float(mj) >= float(mjdini) and float(mj) <= float(mjdend) and len(line.split()) > 5:
            c1 = line[18:27]
            c2 = line[37:46]
            c3 = line[58:68]
            l  = ' '+mj+' '+c1+' '+c2+' '+c3+' '+'\n'
            output.write(l)
            if mj == mjdini+999: print("estoy en el dia D")
        if float(mj) > float(mjdend):
            break
    finaldata.close()
    output.close()
    pass

# Extraction Functions
def PCoeff(data, trace_coeffs, Aperture, RON, Gain, NSigma, S, N, Marsh_alg,min_col,max_col):
    Result      = Marsh.ObtainP((data.flatten()).astype('double'), \
                                   scipy.polyval(trace_coeffs,np.arange(data.shape[1])).astype('double'), \
                                   data.shape[0], data.shape[1], data.shape[1], Aperture, RON, Gain, \
                                   NSigma, S, N, Marsh_alg,min_col,max_col)
    FinalMatrix = np.asarray(Result)                      # After the function, we convert our list to a Numpy array.
    FinalMatrix.resize(data.shape[0],data.shape[1])   # And return the array in matrix-form.
    return FinalMatrix

def PCoeff2(pars):
    trace_coeffs = pars[0]
    Aperture = pars[1]
    RON=pars[2]
    Gain = pars[3]
    NSigma = pars[4]
    S = pars[5]
    N = pars[6]
    Marsh_alg = pars[7]
    min_col = pars[8]
    max_col = pars[9]
    Result      = Marsh.ObtainP((GDATA.flatten()).astype('double'), \
                                   scipy.polyval(trace_coeffs,np.arange(GDATA.shape[1])).astype('double'), \
                                   GDATA.shape[0], GDATA.shape[1], GDATA.shape[1], Aperture, RON, Gain, \
                                   NSigma, S, N, Marsh_alg,min_col,max_col)
    FinalMatrix = np.asarray(Result)                      # After the function, we convert our list to a Numpy array.
    FinalMatrix.resize(GDATA.shape[0],GDATA.shape[1])   # And return the array in matrix-form.
    return FinalMatrix

def obtain_P(data, trace_coeffs, Aperture, RON, Gain, NSigma, S, N, Marsh_alg,min_col,max_col,npools):
    global GDATA
    GDATA = data
    npars_paralel = []

    if 'int' in str(type(min_col)) or 'float' in str(type(min_col)):
        min_col = np.zeros(len(trace_coeffs)) + int(min_col)
    if 'int' in str(type(max_col)) or 'float' in str(type(max_col)):
        max_col = np.zeros(len(trace_coeffs)) + int(max_col)

    for i in range(len(trace_coeffs)):
        npars_paralel.append([trace_coeffs[i,:],Aperture,RON,Gain,NSigma,S,N,Marsh_alg,int(min_col[i]),int(max_col[i])])
    p = Pool(npools)
    spec = np.array((p.map(PCoeff2, npars_paralel)))
    p.terminate()
    return np.sum(spec,axis=0)

def getSpectrum(P,data,trace_coeffs,Aperture,RON,Gain,S,NCosmic, min_col,max_col):
    Result,size = Marsh.ObtainSpectrum( (data.flatten()).astype('double'), \
                                            scipy.polyval(trace_coeffs,np.arange(data.shape[1])).astype('double'), \
                                            P.flatten().astype('double'), data.shape[0],\
                                            data.shape[1],data.shape[1],Aperture,RON,\
                                            Gain,S,NCosmic,min_col,max_col)
    FinalMatrix = np.asarray(Result)                      # After the function, we convert our list to a Numpy array.
    FinalMatrix.resize(3,size)                            # And return the array in matrix-form.
    return FinalMatrix

def getSpectrum2(pars):
    trace_coeffs = pars[0]
    Aperture = pars[1]
    RON=pars[2]
    Gain = pars[3]
    S = pars[4]
    NCosmic = pars[5]
    min_col = pars[6]
    max_col = pars[7]
    Result,size = Marsh.ObtainSpectrum( (GDATA.flatten()).astype('double'), \
                                            scipy.polyval(trace_coeffs,np.arange(GDATA.shape[1])).astype('double'), \
                                            P.flatten().astype('double'), GDATA.shape[0],\
                                            GDATA.shape[1],GDATA.shape[1],Aperture,RON,\
                                            Gain,S,NCosmic,min_col,max_col)
    FinalMatrix = np.asarray(Result)                      # After the function, we convert our list to a Numpy array.
    FinalMatrix.resize(3,size)                            # And return the array in matrix-form.
    return FinalMatrix

def getSimpleSpectrum(data,trace_coeffs,Aperture,min_col,max_col):
    Result   =   Marsh.SimpleExtraction((data.flatten()).astype('double'),\
                                            scipy.polyval(trace_coeffs,np.arange(data.shape[1])).astype('double'), \
                                            data.shape[0],data.shape[1],\
                                            data.shape[1],Aperture,min_col,max_col)
    FinalMatrix = np.asarray(Result) # After the function, we convert our list to a Numpy array.
    return FinalMatrix

def getSimpleSpectrum2(pars):
    trace_coeffs = pars[0]
    Aperture = pars[1]
    min_col = pars[2]
    max_col = pars[3]
    Result   =   Marsh.SimpleExtraction((GDATA.flatten()).astype('double'),\
                                            scipy.polyval(trace_coeffs,np.arange(GDATA.shape[1])).astype('double'), \
                                            GDATA.shape[0],GDATA.shape[1],\
                                            GDATA.shape[1],Aperture,min_col,max_col)
    FinalMatrix = np.asarray(Result) # After the function, we convert our list to a Numpy array.
    return FinalMatrix

def simple_extraction(data,coefs,ext_aperture,min_extract_col,max_extract_col,npools):
    global GDATA
    GDATA = data
    npars_paralel = []
    if 'int' in str(type(min_extract_col)) or 'float' in str(type(min_extract_col)):
        min_extract_col = np.zeros(len(coefs)) + int(min_extract_col)
    if 'int' in str(type(max_extract_col)) or 'float' in str(type(max_extract_col)):
        max_extract_col = np.zeros(len(coefs)) + int(max_extract_col)
    for i in range(len(coefs)):
        npars_paralel.append([coefs[i,:],ext_aperture,int(min_extract_col[i]),int(max_extract_col[i])])
    p = Pool(npools)
    spec = np.array((p.map(getSimpleSpectrum2, npars_paralel)))
    p.terminate()
    return spec

def optimal_extraction(data,Pin,coefs,ext_aperture,RON,GAIN,MARSH,COSMIC,min_extract_col,max_extract_col,npools):
    global GDATA,P
    P = Pin
    GDATA = data
    npars_paralel = []
    if 'int' in str(type(min_extract_col)) or 'float' in str(type(min_extract_col)):
        min_extract_col = np.zeros(len(coefs)) + int(min_extract_col)
    if 'int' in str(type(max_extract_col)) or 'float' in str(type(max_extract_col)):
        max_extract_col = np.zeros(len(coefs)) + int(max_extract_col)

    for i in range(len(coefs)):
        npars_paralel.append([coefs[i,:],ext_aperture,RON,GAIN,MARSH,COSMIC,int(min_extract_col[i]),int(max_extract_col[i])])
    p = Pool(npools)
    spec = np.array((p.map(getSpectrum2, npars_paralel)))
    p.terminate()
    return spec

# CCF functions
def get_herms(horder):
    norms  = []
    herms = []
    for l in range(3,horder+1):
        norms.append( 1.0 / np.sqrt(scipy.misc.factorial(l) * (2 ** l) ) )
        herms.append( scipy.special.hermite( l ) )
    return norms, herms

def XCor(spectra, mask_l, mask_h, mask_w, vel, lbary_ltopo, vel_width=30,\
             vel_step=0.3,  start_order=0, spec_order=9,iv_order=10,sn_order=8,max_vel_rough=300.):
    """
    Calculates the cross-correlation function for a Coralie Spectra
    """
    # speed of light, km/s
    c = 2.99792458E5

    # loop over orders
    norders = spectra.shape[1]

    # determine minimum velocities
    vel_min = vel - vel_width
    vel_max = vel + vel_width
    N = int(np.ceil( (2*vel_width) / vel_step ))

    Xcor_full   = np.zeros( (N, norders+1) )
    sn          = np.zeros( (norders) )
    nlines_used = np.zeros( (norders) )

    velocities = vel_min + np.arange( N ) * vel_step

    Xcor_full[:,0] = velocities

    weight=0.0
    mask_middle = 0.5*(mask_l + mask_h)
    W = np.zeros( norders )

    vwt = 300
    for j in range(start_order,norders):
        t1 = time.time()
        LL = np.where( spectra[spec_order,j,:] != 0 )
        if len(LL[0]) > 0:
            x1 = np.min( LL )
            x2 = np.max( LL )
            w1 = np.argmin( np.absolute( spectra[0,j,:] - spectra[0,j,x1] ) )
            w2 = np.argmin( np.absolute( spectra[0,j,:] - spectra[0,j,x2] ) )
            l1_0 = spectra[0,j,w1] / lbary_ltopo
            l2_0 = spectra[0,j,w2] / lbary_ltopo
            ww1 = np.argmin( np.abs( spectra[0,j,:] - l1_0*(1+(31+max_vel_rough)/c) ) )
            ww2 = np.argmin( np.abs( spectra[0,j,:] - l2_0*(1-(31+max_vel_rough)/c) ) )
            # should not happen, but hey, just in case...
            if (ww1 < w1):
                ww1 = w1
            if (ww2 > w2):
                ww2 = w2
            l1 = spectra[0,j,ww1]
            l2 = spectra[0,j,ww2]
            II = np.where( (mask_l > l1) & (mask_h < l2) )
            #if len(II[0])>0:
            #print j,II[0][0],II[0][-1]
            nlu = len(II[0])
            nlines_used[j] = nlu
            snw1 = int(0.25*spectra.shape[2])
            snw2 = int(0.75*spectra.shape[2])
            if (nlu > 0):
                # calculate median S/N
                #median_sn = np.median( spectra[5,j,w1:w2] * np.sqrt( spectra[6,j,w1:w2] ) )
                median_sn = np.median( spectra[sn_order,j,snw1:snw2] )
                sn[j]     = median_sn
                S = spectra[spec_order,j,w1:w2]
                #iv = spectra[iv_order,j,w1:w2]
                signal2noise = spectra[sn_order,j,w1:w2]
                snwa = np.zeros(N)
                for k in range(N):
                    #print k
                    Xcor_full[k,j+1], snw = CCF.ccfcos(mask_l[II], mask_h[II], spectra[0,j,w1:w2], S,\
                                           mask_w[II], signal2noise, vel_min + k*vel_step)
                    snwa[k] = snw

                    if np.isnan(Xcor_full[k,j+1]):
                        Xcor_full[k,j+1] = Xcor_full[k-1,j+1]
                        snwa[k] = snwa[k-1]
                    #if k ==182 and j==35:
                    #       #print mask_l[II], mask_h[II], spectra[0,j,w1:w2], S,mask_w[II], signal2noise, vel_min + k*vel_step
                    #       #for z in range(len(mask_l[II])):
                    #       #       III =  np.where((spectra[0,j,w1:w2]>=mask_l[II][z])&(spectra[0,j,w1:w2]<=mask_h[II][z]))[0]
                    #       #       print spectra[0,j,w1:w2][III],S[III]
                    #       #print Xcor_full[k,j+1]
                    #       #print snw
                    #       #print gfd

                xc_weight  = np.median( snwa )
                Xcor_full[:,j+1] /= snwa #xc_weight
                W[j] = xc_weight

    return velocities, Xcor_full, sn, nlines_used, W

def XC_Herm_Fit(X,Y,back_lag=5, usemin=True, horder=20, sigma_res = 2, horder_res_herm=10, sigma_res_herm = 2.5):
    """
    Fits a Gauss-Hermite series to the XC function
    Higher order in the expansion is horder

    """
    vel0 = X[len(X)/2]
    def fitfunc(p,x,horder, norms, herms):
        gau = p[0]*CorGaussian(x,p[1],p[2])
        her = 0.0
        xnorm = (x-p[1])/p[2]
        for l in range(3,horder+1):
#            norm  = 1.0 / np.sqrt(scipy.misc.factorial(l) * (2 ** l) )
#            her  += p[4 + (l-3)] * norm*scipy.polyval( scipy.special.hermite( l ), x )
            her  += p[4 + (l-3) ] * norms[l-3] * scipy.polyval( herms[l-3], xnorm )
        return gau*(1.0 + her) + 1.0 + p[3]

    def errfunc(p, x, y, h, norms, herms):
        clutch = 0.0
        mean = p[1]
        if (mean < np.min(x)):
            clutch = 1e10*(1.0 - exp(-np.abs(mean-np.min(x)) / 3) )
        if (mean > np.max(x)):
            clutch = 1e10*(1.0 - exp(-np.abs(mean-np.max(x)) / 3) )
        return np.ravel( (fitfunc(p,x,h,norms, herms) - y) )  + clutch

    if (horder >= 3):
        n = 4 + (horder-2)
    else:
        n = 4


    norms, herms = get_herms(horder)

    p0 = np.zeros( n )
    p0[3] = 1.0 #0.5*( np.median(Y[:back_lag]) + np.median(Y[-back_lag:]) )
    p0[2] = 5
    if (usemin):
        p0[1] = vel0
        p0[0] = np.min(Y) - p0[3]
    else:
        p0[1] = vel0
        p0[0] = np.max(Y) - p0[3]

    p1, success = scipy.optimize.leastsq(errfunc,p0, args=(X,Y,horder,norms, herms))
    predicted = fitfunc(p1,X,horder,norms,herms)

    if (horder < horder_res_herm):
        horder_res_herm = horder

    # now run h=horder_res_herm and Gaussian fit to +- sigma_res* sigma of mean
    mean = p1[1]
    sigma = p1[2]

    L1 = np.where( np.abs(X - mean) <= sigma_res_herm* sigma)
    horder=horder_res_herm
    if (horder >= 3):
        n = 4 + (horder-2)
    else:
        n = 4

    if (len(L1[0]) > 0):
        norms, herms = get_herms(horder)
        p0 = p0[0:n]
        p1_s, success_s = scipy.optimize.leastsq(errfunc,p0, args=(X[L1],Y[L1],horder,norms,herms))
        predicted_s   = fitfunc(p1_s,X[L1],horder,norms,herms)
    else:
        p1_s = np.zeros(n)
        predicted_s = np.zeros( len(X[L1]) )

    L2 = np.where( np.abs(X - mean) <= sigma_res * sigma)
    # the following is just a Gaussian fit
    horder = 2
    if (len(L2[0]) > 0):
        norms, herms = get_herms(horder)
        p0 = p0[0:4]
        p1_gau, success_gau = scipy.optimize.leastsq(errfunc,p0, args=(X[L2],Y[L2],horder,norms,herms))
        predicted_gau = fitfunc(p1_gau,X[L2],horder,norms,herms)
    else:
        p1_gau = np.zeros(n)
        predicted_gau = np.zeros( len(X[L2]) )

    return p1, predicted, p1_s, predicted_s, p1_gau, predicted_gau, L1, L2
    #return p1, predicted, p1_gau, predicted_gau, L2

def XC_Final_Fit( X, Y, usemin=True, sigma_res = 1.5, horder=20, moonv = 0., moons = 1., moon = False):
    """
    Fits a Gaussian and Gauss-Hermite series
    Higher order in the expansion is horder

    """
    f0 = 0.1
    vel0 = X[len(X)/2]

    def fitfunc(p,x,horder, norms, herms):
        gau = p[0]*CorGaussian(x,p[1],p[2])
        her = 0.0
        xnorm = (x-p[1])/p[2]
        for l in range(3,horder+1):
            her  += p[4 + (l-3) - 1] * norms[l-3] * scipy.polyval( herms[l-3], xnorm )
        return gau*(1.0 + her) + 1.0 # + p[3]

    def errfunc(p, x, y, h, norms, herms):
        clutch = 0.0
        mean = p[1]
        if (mean < np.min(x)):
            clutch = 1e10*(1.0 - exp(-np.abs(mean-np.min(x)) / 3) )
        if (mean > np.max(x)):
            clutch = 1e10*(1.0 - exp(-np.abs(mean-np.max(x)) / 3) )
        return np.ravel( (fitfunc(p,x,h,norms, herms) - y) )  + clutch

    def fitfunc2(p,x,horder, norms, herms, moonv,moons):
        gau = p[0]*CorGaussian2(x, p[1], p[2], moonv, moons,p[3])
        her = 0.0
        xnorm = (x-p[1])/p[2]
        for l in range(3,horder+1):
            her  += p[4 + (l-3) - 1] * norms[l-3] * scipy.polyval( herms[l-3], xnorm )
        return gau*(1.0 + her) + 1.0 # + p[3]

    def errfunc2(p, x, y, h, norms, herms, moonv,moons):
        clutch = 0.0
        mean = p[1]
        if (mean < np.min(x)):
            clutch = 1e10*(1.0 - exp(-np.abs(mean-np.min(x)) / 3) )
        if (mean > np.max(x)):
            clutch = 1e10*(1.0 - exp(-np.abs(mean-np.max(x)) / 3) )
        return np.ravel( (fitfunc2(p,x,h,norms, herms, moonv,moons) - y) )  + clutch

    # Gauss-Hermit fit
    if (horder >= 3):
        n = 4 + (horder-2) - 1
    else:
        n = 4 - 1

    norms, herms = get_herms(horder)
    p0 = np.zeros( n )
    p0[3] = 0.0
    p0[2] = 3.5
    p0[1] = vel0

    if (usemin):
        p0[0] = np.min(Y) - (1.0 + p0[3])
    else:
        p0[0] = np.max(Y) - (1.0 + p0[3])

    p1, success = scipy.optimize.leastsq(errfunc,p0, args=(X,Y,horder,norms, herms))
    predicted = fitfunc(p1,X,horder,norms,herms)
    mean = p1[1]


    # the following are pars for a Gaussian fit
    horder = 2
    if (horder >= 3):
        n = 4 + (horder-2) - 1
    else:
        n = 4 - 1
    p0 = p0[0:n]

    # first pass, assume sigma = np.ceil(vel_width/6)
    sigma0 = np.ceil( np.abs(vel0-X[0])/6.0 )

    L2 = np.where( np.abs(X - mean) <= sigma_res * sigma0)

    if (len(L2[0]) > 0):
        norms, herms = get_herms(horder)
        p1_gau0, success_gau = scipy.optimize.leastsq(errfunc,p0, args=(X[L2],Y[L2],horder,norms,herms))
    else:
        p1_gau0 = np.zeros(n)

    # second pass, use actual measured sigma and mean
    mean  = p1_gau0[1]
    sigma = p1_gau0[2]
    L2 = np.where( np.abs(X - mean) <= sigma_res * sigma)

    if moon:
        p1_gau0 = np.append(p1_gau0, f0)
        if (len(L2[0]) > 0):
            norms, herms = get_herms(horder)
            p1_gau, success_gau = scipy.optimize.leastsq(errfunc2,p1_gau0, args=(X[L2],Y[L2],horder,norms,herms,moonv,moons))
            predicted_gau = fitfunc2(p1_gau,X[L2],horder,norms,herms,moonv,moons)
        else:
            p1_gau = np.zeros(n)
            predicted_gau = np.zeros( len(X[L2]) )
    else:
        if (len(L2[0]) > 0):
            norms, herms = get_herms(horder)
            p1_gau, success_gau = scipy.optimize.leastsq(errfunc,p1_gau0, args=(X[L2],Y[L2],horder,norms,herms))
            predicted_gau = fitfunc(p1_gau,X[L2],horder,norms,herms)
        else:
            p1_gau = np.zeros(n)
            predicted_gau = np.zeros( len(X[L2]) )

    return p1, predicted, p1_gau, predicted_gau, L2


def Average_CCF(xc_full, sn, start_order=0,sn_min=0.15, Simple=False, W=None, boot_ind=[]):
    """
    Given the CCF of all orders, construct an average CCF
    """
    norders = xc_full.shape[1] - 1
    N = xc_full.shape[0]

    xc_av = np.zeros( N )

    ws = 0.0

    combine_indices = range(start_order, norders)
    #print combine_indices
    if len(boot_ind)>0:
        combine_indices = boot_ind

    for order in combine_indices:
        #I = np.where(np.isnan(xc_full[:,order+1]))
        #print I
        sum_xc = np.sum( xc_full[:,order+1] )
        #print order, sum_xc, xc_full[:,order+1]
        if (sum_xc > 0):
            if (Simple):
                norm = np.median( xc_full[:,order+1] )
                if (norm > 0) and (sn[order] > sn_min):
                    if (W is None):
                        ws += 1 #sn[order]
                        xc_av += xc_full[:,order+1]#/norm
                    else:
                        ws    += W[order]
                        xc_av += W[order] * xc_full[:,order+1]

            else:
                p1,XCmodel,p1s,XCmodels,p1gau,XCmodelgau, Ls = XC_Herm_Fit(xc_full[:,0], xc_full[:,order+1],horder=8, sigma_res=3)
                CCF_norm = xc_full[:,order+1] / p1[3]
                weight_sn    = sn[order]
                if (sn[order] > sn_min):
                    ws       += weight_sn
                    xc_av    += weight_sn * CCF_norm

    IN = np.where(xc_av!=0)[0]
    if len(IN) == 0 or ws == 0:
        xc_av = np.ones(len(xc_av))
    else:
        xc_av /= ws

    return xc_av

def IntGaussian(x,mu,sigma):
    """

    Returns Gaussian integrated over a pixel

    """
    s2 = sqrt(2)
    arg1 = (x+0.5-mu)/(s2*sigma)
    arg2 = (x-0.5-mu)/(s2*sigma)
    ret = 0.5*(special.erf(arg1) - special.erf(arg2))
    return ret

def CorGaussian(x,mu,sigma):
    """ Returns Gaussian """
    norm = 1.0 / (np.sqrt(2 * np.pi) * sigma)
    clutch = 0
    if (sigma < 0):
        clutch = 1e32
    return norm * np.exp(-0.5*((x-mu)**2) / sigma**2) + clutch

def CorGaussian2(x,mu,sigma,mu2,sigma2,f):
    """ Returns Gaussian """
    norm = 1.0 / (np.sqrt(2 * np.pi) * sigma)
    norm2 = 1.0 / (np.sqrt(2 * np.pi) * sigma2)
    clutch = 0
    if (sigma < 0):
        clutch = 1e32
    return norm * np.exp(-0.5*((x-mu)**2) / sigma**2) + clutch + f * np.exp(-0.5*((x-mu2)**2) / sigma2**2)

def get_rough_offset(sc,files,window=100):
    i = 0
    xct = []
    while i < len(files):
        spec = sc[i]
        f = open(files[i]).readlines()
        pixel_centers_0 = []
        for line in f:
            w = line.split()
            nlines = int(w[0])
            for j in range(nlines):
                pixel_centers_0.append(float(w[2*j+1]))

        ml = array(pixel_centers_0) - 2
        mh = array(pixel_centers_0) + 2
        xc,offs = XCorPix( spec, ml, mh,del_width=window)
        if len(xct) == 0:
            xct = xc
        else:
            xct += xc
        i+=1

    ind_max = np.argmax( xct )
    delta   = offs[ind_max]
    #plot(offs,xct)
    #show()
    return delta

def fit_these_lines(waves_ob,filename,spec,order,wei, rough_shift = 0.0, del_width=5.0, binning=1,\
                                        line_width=4, fact=1,do_xc=True,sigmai=2.2):

    f = open(filename).readlines()
    pixel_centers = array([])
    wavelengths   = array([])
    sigmas        = array([])
    centroids     = array([])
    intensities   = array([])

    if do_xc:
        pixel_centers_0 = []
        for line in f:
            w = line.split()
            nlines = int(w[0])
            for j in range(nlines):
                pixel_centers_0.append(float(w[2*j+1])*fact/float(binning) + rough_shift)
        ml = array(pixel_centers_0) - 2
        mh = array(pixel_centers_0) + 2
        xc,offs = XCorPix( spec, ml, mh, del_width=del_width)
        ind_max = np.argmax( xc )
        delta   = offs[ind_max]
    else:
        delta=0.

    N_l = 0
    bad_indices = []
    bad_indices_ct = 0
    #plot(spec)
    for line in f:
        #print line
        w = line.split()
        # extract info, and fit line(s)
        nlines = int(w[0])
        pix = []
        wav = []
        for j in range(nlines):
            if float(w[2*j+1])*fact/float(binning) + rough_shift+delta > 20 and float(w[2*j+1])*fact/float(binning) + rough_shift+delta < len(spec)-20:
                pix.append(float(w[2*j+1])*fact/float(binning) + delta + rough_shift)
                wav.append(float(w[2*j+2]))
        if len(pix) > 0:
            N_l += len(pix)
            pix = np.array(pix)
            #pix2=np.around(pix).astype('int')
            #plot(pix2,spec[pix2],'ro')
            wav = np.array(wav)
            xmin = int(round(min(pix)))
            xmax = int(round(max(pix)))
            X = array(range(xmin-line_width,xmax+line_width+1))
            Y = spec[xmin-line_width:xmax+line_width+1]
            if (nlines == 1):
                num       = np.sum(X*Y)
                den       = np.sum(Y)
                if (den > 0):
                    Cent = num/den
                else:
                    Cent = -1
            weight = wei[xmin-line_width:xmax+line_width+1]
            kk = np.where( weight == 0)
            # Input Spectrum is background subtracted ---> B=0
            B = np.zeros(len(X))
            mu = pix
            sigma = np.zeros(nlines) + sigmai * fact / float(binning)
            #plot(X,Y)
            #show()
            #print X, Y, B, mu, sigma, weight
            p1, suc = LineFit_SingleSigma( X, Y, B, mu, sigma, weight)
            #print p1
            #if (suc<1) or (suc > 4):
            #    print "Problem", order, X, delta
            # collect fit information
            #plot(X,Y)
            reto =  p1[0]*np.exp((X-p1[1])**2/(-0.5*p1[2]**2))
            #print reto
            #plot(X,reto,'r')

            wavelenghts = np.append(wavelengths,wav)
            for j in range(nlines):
                pixel_centers = np.append(pixel_centers,p1[3*j + 1])
                sigmas        = np.append(sigmas,p1[3*j + 2])
                wavelengths   = np.append(wavelengths,wav[j])
                intensities   = np.append(intensities,p1[3*j])
                if (nlines == 1):
                    centroids = np.append(centroids, Cent)
                else:
                    centroids = np.append(centroids, -1)
    wavelengths_co,pixel_centers_co,intensities_co,sigmas_co,centroids_co = [],[],[],[],[]
    for i in range(len(wavelengths)):
        if np.around(wavelengths[i],4) in np.around(waves_ob,4):
            wavelengths_co.append(wavelengths[i])
            pixel_centers_co.append(pixel_centers[i])
            intensities_co.append(intensities[i])
            sigmas_co.append(sigmas[i])
            centroids_co.append(centroids[i])
    wavelengths_co,pixel_centers_co,intensities_co,sigmas_co,centroids_co = \
       np.array(wavelengths_co),np.array(pixel_centers_co),np.array(intensities_co), \
       np.array(sigmas_co),np.array(centroids_co)
    #plot(np.around(pixel_centers_co).astype('int'),spec[np.around(pixel_centers_co).astype('int')],'ro')
    #show()

    return wavelengths_co,pixel_centers_co,intensities_co,sigmas_co,centroids_co



def Initial_Wav_Calibration(filename,spec,order,wei, porder=3, rmsmax=75, minlines=10, FixEnds=True,\
            Dump_Argon=False, Dump_AllLines=False, Cheby=False, rough_shift = 0.0,\
            del_width=5.0, binning=1,line_width=4, fact=1,do_xc=True,sigmai=2.2,pixelization=False):

    f = open(filename).readlines()
    pixel_centers = array([])
    wavelengths   = array([])
    sigmas        = array([])
    centroids     = array([])
    intensities   = array([])

    if do_xc:
        pixel_centers_0 = []
        for line in f:
            w = line.split()
            nlines = int(w[0])
            for j in range(nlines):
                pixel_centers_0.append(float(w[2*j+1])*fact/float(binning) + rough_shift)

        ml = array(pixel_centers_0) - 2
        mh = array(pixel_centers_0) + 2
        xc,offs = XCorPix( spec, ml, mh, del_width=del_width)
        ind_max = np.argmax( xc )
        delta   = offs[ind_max]
    else:
        delta=0.

    #print "Computed offset for order ", order, " is ", delta
    N_l = 0
    bad_indices = []
    bad_indices_ct = 0
    #print order
    #plot(spec)
    for line in f:
        if line[0]!='#':
            w = line.split()
            # extract info, and fit line(s)
            nlines = int(w[0])
            pix = []
            wav = []
            for j in range(nlines):
                if float(w[2*j+1])*fact/float(binning) + rough_shift+delta > 20 and \
                        float(w[2*j+1])*fact/float(binning) + rough_shift+delta < len(spec)-20:
                    pix.append(float(w[2*j+1])*fact/float(binning) + delta + rough_shift)
                    wav.append(float(w[2*j+2]))

            if len(pix) > 0:
                pix = np.array(pix)
                wav = np.array(wav)
                xmin = int(round(min(pix)))
                xmax = int(round(max(pix)))
                X = array(range(xmin-line_width,xmax+line_width+1))
                Y = spec[xmin-line_width:xmax+line_width+1]
                if len(np.where(Y!=0)[0])>0:
                    N_l += len(pix)
                    if (nlines == 1):
                        num       = np.sum(X*Y)
                        den       = np.sum(Y)
                        if (den > 0):
                            Cent = num/den
                        else:
                            Cent = -1

                    weight = wei[xmin-line_width:xmax+line_width+1]
                    kk = np.where( weight == 0)
                    # Input Spectrum is background subtracted ---> B=0
                    B = np.zeros(len(X))
                    mu = pix
                    sigma = np.zeros(nlines) + sigmai * fact / float(binning)

                    p1, suc = LineFit_SingleSigma( X, Y, B, mu, sigma, weight,pixelization=pixelization)
                    #print p1
                    #if (suc<1) or (suc > 4):
                    #    print "Problem", order, X, delta
                    # collect fit information
                    #reto =  fitfunc_temp(p1, X, len(pix),pixelization= pixelization)
                    #plot(X,reto,'r')
                    wavelenghts = np.append(wavelengths,wav)
                    for j in range(len(pix)):
                        pixel_centers = np.append(pixel_centers,p1[3*j + 1])
                        sigmas        = np.append(sigmas,p1[3*j + 2])
                        wavelengths   = np.append(wavelengths,wav[j])
                        intensities   = np.append(intensities,p1[3*j])
                        if (nlines == 1):
                            centroids = np.append(centroids, Cent)
                        else:
                            centroids = np.append(centroids, -1)

    #print len(pixel_centers)
    pixel_centers2 = np.around(pixel_centers).astype('int')
    I = np.where((pixel_centers2>=0) & (pixel_centers2<len(spec)))
    pixel_centers2 = pixel_centers2[I]
    #plot(pixel_centers2,spec[pixel_centers2],'go')
    #show()
    #print gfd

    #I = np.where((pixel_centers>0) & (pixel_centers<2048))[0]
    #plot(np.around(pixel_centers[I]).astype('int'),spec[np.around(pixel_centers[I]).astype('int')],'ro')
    #show()
    I1 = np.where(pixel_centers<50)[0]
    I2 = np.where(pixel_centers>len(spec)-50)[0]
    II = np.hstack((I1,I2))
    bad_indices = np.hstack((np.array(bad_indices),II))
    bad_indices = list(np.unique(bad_indices))
    # now, do the polynomial fit, rejecting some lines until RMS is below rmsmax
    I = range( N_l )

    for bi in bad_indices:
        I.remove( bi )
        N_l -= 1

    if (Cheby):
        coeffs_pix2wav   = Cheby_Fit(pixel_centers[I], wavelengths[I], porder,len(spec))
        coeffs_pix2sigma = Cheby_Fit(pixel_centers[I], sigmas[I], porder,len(spec))
    else:
        coeffs_pix2wav   = scipy.polyfit(pixel_centers[I], wavelengths[I], porder)
        coeffs_pix2sigma = scipy.polyfit(pixel_centers[I], sigmas[I], porder)

    rmsms, residuals = rms_ms(coeffs_pix2wav, pixel_centers[I], wavelengths[I], len(spec), Cheby=Cheby)

    if (FixEnds):
        minI = np.min( I ) + 1
        maxI = np.max( I ) - 1
    else:
        minI = np.min( I )
        maxI = np.max( I )
    #if order==26:
    #           plot(pixel_centers[I],residuals,'ro')
    #           plot([0,4096],[0,0])
    #plot(np.arange(4096),Cheby_eval(coeffs_pix2wav,np.arange(4096),len(spec)))
    #show()
    #print dfgh
    count = 0
    while ((N_l > minlines) and (rmsms > rmsmax)):
        rmsms, residuals = rms_ms(coeffs_pix2wav, pixel_centers[I], wavelengths[I], len(spec), Cheby=Cheby)
        index_worst = np.argmax( np.absolute(residuals) )
        I.pop( index_worst)
        N_l -= 1
        if (Cheby):
            coeffs_pix2wav   = Cheby_Fit(pixel_centers[I], wavelengths[I], porder,len(spec))
            coeffs_pix2sigma = Cheby_Fit(pixel_centers[I], sigmas[I], porder,len(spec))
        else:
            coeffs_pix2wav   = scipy.polyfit(pixel_centers[I], wavelengths[I], porder)
            coeffs_pix2sigma = scipy.polyfit(pixel_centers[I], sigmas[I], porder)
        count +=1

    rmsms, residuals = rms_ms(coeffs_pix2wav, pixel_centers[I], wavelengths[I], len(spec), Cheby=Cheby)

    pci = np.around(pixel_centers).astype('int')
    #plot(spec)
    #plot(pci,spec[pci],'ro')
    #plot(pci[I],spec[pci[I]],'bo')
    #show()
    #plot(wavelengths[I],residuals,'ro')
    #show()
    #print "RMS is ", rmsms, "using ", N_l, " lines at indices ", I
    #plot(pixel_centers[I],wavelengths[I],'ro')
    #if order == 26:
    #           plot(pixel_centers[I],residuals-0.1,'bo')
    #           plot([0,4096],[-0.1,-0.1])
    #           #plot(np.arange(4096),Cheby_eval(coeffs_pix2wav,np.arange(4096),len(spec)))
    #           show()
    #print order, len(pixel_centers), len(I)
    return coeffs_pix2wav, coeffs_pix2sigma, pixel_centers[I], wavelengths[I], \
    rmsms, residuals, centroids[I], sigmas[I], intensities[I]

def XCorPix(spectra, mask_l, mask_h, del0=0, del_width=5, del_step=0.1):
    """
    Calculates the cross-correlation function for a Spectra, pixel space, linear shifts
    Single order

    Depends on CCF module
    """

    # determine minimum velocities
    del_min = del0 - del_width
    del_max = del0 + del_width
    N = int(np.ceil( (2*del_width) / del_step ).astype('int'))

    Xcor = np.zeros( N )

    deltas = del_min + np.arange( N ) * del_step

    LL = np.where( spectra != 0 )
    X = np.arange( len(spectra) )
    x1 = np.min( LL )
    x2 = np.max( LL )
    #print x1, x2
    for k in range(N):
        Xcor[k] += CCF.ccfpix(mask_l, mask_h, X[x1:x2], spectra[x1:x2], del_min + k*del_step)

    return Xcor, deltas

def rms_ms(coeffs_pix2wav, pixel_centers, wavelengths, npix, Cheby=False):
    " Returns rms deviation of best fit in m/s"

    if (Cheby):
        residuals = Cheby_eval(coeffs_pix2wav,pixel_centers,npix) - wavelengths
        central_wav = 0.5 * (Cheby_eval(coeffs_pix2wav,50.,npix) + Cheby_eval(coeffs_pix2wav,npix-50,npix))
    else:
        residuals = scipy.polyval(coeffs_pix2wav,pixel_centers) - wavelengths
        central_wav = 0.5 * (scipy.polyval(coeffs_pix2wav,50.) + scipy.polyval(coeffs_pix2wav,npix-50))

    rms_ms = np.sqrt( np.var( residuals ) ) * 299792458.0 / central_wav

    return rms_ms, residuals

def Cheby_Fit(x,y,order,npix):
    """
    Fits Chebyshev polynomials to y as a function of x

    """
    med = .5*npix
    #normalize x
    x_norm = (x-med) / med

    def fitfunc(p,chebs,order):
        ret_val = 0.0
        for i in range(0,order+1):
            ret_val += p[order-i]*chebs[i]
        return ret_val
    errfunc = lambda p,chebs,y,order: np.ravel( (fitfunc(p,chebs,order)-y) )

    def get_chebs(x,order):
        chebs = []
        for i in range(0,order+1):
            chebs.append( scipy.special.chebyt(i)(x) )
        return chebs

    p0 = np.zeros( order + 1 )
    p0[order] = np.mean( y )
    chebs = get_chebs( x_norm, order)
    p1, success = scipy.optimize.leastsq(errfunc, p0, args=(chebs, y, order))
    return p1

def Cheby_eval(p,x,npix):
    """
    evaluates Chebyshev polynomial fit at x given best-fit parameters p
    """
    med = .5*npix
    x_norm = (x-med) / med
    order = len(p) - 1
    ret_val = 0.0
    for i in range(order + 1):
        ret_val += p[order - i]*scipy.special.chebyt(i)(x_norm)

    return ret_val

def fitfunc_temp(p, x, n, pixelization=True):

    if pixelization:
        lxo = len(x)
        xo = x.copy()
        x = np.arange(x[0]-0.5,x[-1]+0.5,0.01)

    ret = np.zeros(len(x))
    for i in range(n):
        ret += ( p[i*2+1] * IntGaussian(x,p[i*2+2],p[0]) )

    if pixelization:
        ret = ret.reshape((lxo,100))
        ret = np.mean(ret,axis=1)

    return ret

def LineFit_SingleSigma(X, Y, B, mu, sigma, weight,pixelization=False):
    """
    This function fits a series of Gaussians simultaneously, given a set
    of input pixels, sigmas, and intensities

    Sigma is the same for all lines

    Returns (mu_i, sigm_i, int_i), i=1,.,n where n is the number of components
    """

    # get number of components to be fit
    n = len(mu)

    def fitfunc(p, x, n):

        if pixelization:
            lxo = len(x)
            xo = x.copy()
            x = np.arange(x[0]-0.5,x[-1]+0.5,0.01)

        ret = np.zeros(len(x))
        for i in range(n):
            ret += ( p[i*2+1] * IntGaussian(x,p[i*2+2],p[0]) )

        if pixelization:
            ret = ret.reshape((lxo,100))
            ret = np.mean(ret,axis=1)


        return ret
    errfunc = lambda p,x,n,y,weight: np.ravel( (fitfunc(p,x,n)-y)*weight )

    # Build error array
    #err = sqrt((B+Y)/gain + readnoise**2)

    # Build input parameters list
    # length of inout parameter list: 3 * n
    p0 = np.zeros(2*n+1)
    p0[0] = sigma[0]
    for i in range(n):
        p0[i*2+1] = Y[int(round(mu[i])-X[0])]
        p0[i*2+2] = mu[i]

    # perform fit
    #plot(X,Y-B,'b')
    p1, success = scipy.optimize.leastsq(errfunc,p0, args=(X,n,Y-B, weight))
    #print p1
    #plot(X,fitfunc(p1,X,n),'r')
    # build output consistent with LineFit
    p_output = np.zeros(3*n)
    for i in range(n):
        p_output[i*3]   = p1[i*2 + 1]
        p_output[i*3+1] = p1[i*2 + 2]
        p_output[i*3+2] = p1[0]

    return p_output, success

def Fit_Global_Wav_Solution(pix_centers, wavelengths, orders, Wgt, p0, minlines=1000, maxrms=150, order0=89, ntotal= 70, npix=2048, Cheby=False, Inv = False,nx=5,nm=6):
    """
    Given x_i, lamda_i and m_i fit for a solution of the form

    lambda_i = (1/m_i) * Joint_Polynomial(x_i,m_i)

    where m_i = raw_order_numer + 89
    [89 comes from running Find_m on an image solved with order-by-order wav cal]

    """
    def fitfunc(p, x, m):
        ret = (1.0/m) * Joint_Polynomial(p,x,m)
        return ret
    errfunc = lambda p,x,m,y: np.ravel( (fitfunc(p,x,m)-y) )

    def fitfunc_cheb(p,chebs,m):
        ret = (1.0/m) * Joint_Polynomial_Cheby(p,chebs,nx=nx,nm=nm)
        return ret
    errfunc_cheb = lambda p,chebs,y,m, w: np.ravel( w*(fitfunc_cheb(p,chebs,m)-y) )
    errfunc_cheb_nw = lambda p,chebs,y,m: np.ravel( (fitfunc_cheb(p,chebs,m)-y) )

    #print order0, ntotal, npix, nx, nm, minlines, maxrms, Cheby, Inv
    #print gfd

    if (Cheby):
        chebs = Calculate_chebs(pix_centers, orders+order0, Inverse=Inv,order0=order0,ntotal=ntotal,npix=npix,nx=nx,nm=nm)
        p1, success =  scipy.optimize.leastsq(errfunc_cheb, p0, args=(chebs, wavelengths, orders + order0, Wgt))
        residuals    = errfunc_cheb_nw(p1, chebs, wavelengths, orders + order0)
    else:
        p1, success = scipy.optimize.leastsq(errfunc, p0, args=(pix_centers, orders + order0, wavelengths))
        residuals    = errfunc(p1, pix_centers, orders + order0, wavelengths)

    #for oo in np.unique(orders):
    #   III = np.where(orders==oo)[0]
    #   plot(wavelengths[III],residuals[III],'.')
    #show()

    residuals_ms = 299792458.0 * residuals / wavelengths
    rms_ms       = np.sqrt( np.var( residuals_ms ) )

    N_l = len( pix_centers )
    I = range( N_l )

    cond = 1
    L = np.where( np.absolute(residuals_ms) > 4.0*rms_ms )
    if ( (len(L[0]) == 0) and (rms_ms < maxrms) ) or (N_l < minlines):
        cond=0

    print("\t\t\tStart Global culling with ", N_l, " number of lines")
    bad_wavs, bad_ords = [],[]
    while (cond):
        index_worst = np.argmax( np.absolute(residuals) )
        bad_wavs.append(wavelengths[I[index_worst]])
        bad_ords.append(orders[I[index_worst]])
        #print orders[I[index_worst]],wavelengths[I[index_worst]],rms_ms
        I.pop( index_worst )
        N_l -= 1
        if (Cheby):
            chebs = Calculate_chebs(pix_centers[I], orders[I]+order0, Inverse=Inv,nx=nx,nm=nm,ntotal=ntotal,npix=npix,order0=order0)
            p1, success =  scipy.optimize.leastsq(errfunc_cheb, p0, args=(chebs, wavelengths[I], orders[I] + order0, Wgt[I]))
            residuals    = errfunc_cheb_nw(p1, chebs, wavelengths[I], orders[I] + order0 )
        else:
            p1, success = scipy.optimize.leastsq(errfunc, p0, args=(pix_centers[I], orders[I] + order0, wavelengths[I]))
            residuals    = errfunc(p1, pix_centers[I], orders[I] + order0, wavelengths[I])

        residuals_ms = 299792458.0 * residuals / wavelengths[I]
        rms_ms       = np.sqrt( np.var( residuals_ms ) )

        L = np.where( np.absolute(residuals_ms) > 4.*rms_ms )
        #print rms_ms, maxrms, len(L[0])
        if ( (len(L[0]) == 0) and (rms_ms < maxrms) ) or (N_l < minlines):
            cond=0
        #print "Eliminated line ", index_worst, " at order ", orders[index_worst]
        #print "RMS is ", rms_ms, " after elimination", len(L[0]), rms_ms, maxrms, N_l
    #for oo in np.unique(orders[I]):
    #   III = np.where(orders[I]==oo)[0]
    #   plot(wavelengths[I][III],residuals[III],'.')
    #   plot(np.median(wavelengths[I][III]),np.median(residuals[III]),'ko')
    #show()
    """
    bad_wavs,bad_ords = np.array(bad_wavs), np.array(bad_ords)
    for i in np.unique(bad_ords):
        J = np.where(bad_ords==i)[0]
        tmpwvs,tmpords = bad_wavs[J],bad_ords[J]
        ist = np.argsort(tmpwvs)
        for j in ist:
                print tmpords[j],tmpwvs[j]
    """
    print("\t\t\tFinal RMS is ", rms_ms)
    print("\t\t\tNumber of lines is ", N_l)
    print("\t\t\t--> Achievable RV precision is ", rms_ms/np.sqrt(N_l))

    return p1, pix_centers[I], orders[I], wavelengths[I], I, rms_ms, residuals

def Global_Wav_Solution_vel_shift(pix_centers, wavelengths, orders, Wgt, p_ref, minlines=1000, maxrms=150, order0=89, ntotal= 70, npix=2048, Cheby=False, Inv = False,nx=5,nm=6):
    """
    Given x_i, lamda_i and m_i fit for a solution of the form

    lambda_i = (1/m_i) * Joint_Polynomial(x_i,m_i)

    where m_i = raw_order_numer + 89
    [89 comes from running Find_m on an image solved with order-by-order wav cal]

    """
    def fitfunc(p, p_ref,x, m):
        ret = ((1+1e-6*p)/m) * Joint_Polynomial(p_ref,x,m)
        return ret
    errfunc = lambda p,p_ref,x,m,y: np.ravel( (fitfunc(p,p_ref,x,m)-y) )

    def fitfunc_cheb(p,p_ref,chebs,m):
        ret = ((1+1e-6*p)/m) * Joint_Polynomial_Cheby(p_ref,chebs,nx=nx,nm=nm)
        return ret
    errfunc_cheb    = lambda p,p_ref,chebs,y,m, w: np.ravel( w*(fitfunc_cheb(p,p_ref,chebs,m)-y) )
    errfunc_cheb_nw = lambda p,p_ref,chebs,y,m: np.ravel( (fitfunc_cheb(p,p_ref,chebs,m)-y) )

    #print order0, ntotal, npix, nx, nm, minlines, maxrms, Cheby, Inv
    p0 = np.array( [0] )

    if (Cheby):
        chebs = Calculate_chebs(pix_centers, orders+order0, order0=order0, ntotal=ntotal, npix=npix, Inverse=Inv,nx=nx,nm=nm)
        p1, success =  scipy.optimize.leastsq(errfunc_cheb, p0, args=(p_ref,chebs, wavelengths, orders + order0, Wgt))
        residuals    = errfunc_cheb_nw(p1, p_ref, chebs, wavelengths, orders + order0)
    else:
        p1, success = scipy.optimize.leastsq(errfunc, p0, args=(p_ref, pix_centers, orders + order0, wavelengths))
        residuals    = errfunc(p1, p_ref, pix_centers, orders + order0, wavelengths)

    residuals_ms = 299792458.0 * residuals / wavelengths
    rms_ms       = np.sqrt( np.var( residuals_ms ) )

    N_l = len( pix_centers )
    I = range( N_l )

    cond = 1
    L = np.where( np.absolute(residuals_ms) > 4.0*rms_ms )
    if ( (len(L[0]) == 0) and (rms_ms < maxrms) ) or (N_l < minlines):
        cond=0

    #for oo in np.unique(orders[I]):
    #   III = np.where(orders[I]==oo)[0]
    #   plot(wavelengths[I][III],residuals[III],'.')
    #show()

    print("\t\t\tStart Global culling with ", N_l, " number of lines")
    while (cond):
        index_worst = np.argmax( np.absolute(residuals) )
        I.pop( index_worst )
        N_l -= 1
        if (Cheby):
            chebs = Calculate_chebs(pix_centers[I], orders[I]+order0, order0=order0, ntotal=ntotal, npix=npix, Inverse=Inv,nx=nx,nm=nm)
            p1, success =  scipy.optimize.leastsq(errfunc_cheb, p0, args=(p_ref,chebs, wavelengths[I], orders[I] + order0, Wgt[I]))
            residuals    = errfunc_cheb_nw(p1, p_ref, chebs, wavelengths[I], orders[I] + order0 )
        else:
            p1, success = scipy.optimize.leastsq(errfunc, p0, args=(p_ref, pix_centers[I], orders[I] + order0, wavelengths[I]))
            residuals    = errfunc(p1, p_ref,  pix_centers[I], orders[I] + order0, wavelengths[I])

        residuals_ms = 299792458.0 * residuals / wavelengths[I]
        rms_ms       = np.sqrt( np.var( residuals_ms ) )
        #print 'p1',(1e-6*p1)*299792458.0
        L = np.where( np.absolute(residuals_ms) > 4.*rms_ms )
        if ( (len(L[0]) == 0) and (rms_ms < maxrms) ) or (N_l < minlines):
            cond=0
        #print "Eliminated line ", index_worst, " at order ", orders[index_worst]
        #print "RMS is ", rms_ms, " after elimination", len(L[0]), rms_ms, maxrms, N_l

    print("\t\t\tFinal RMS is ", rms_ms)
    print("\t\t\tNumber of lines is ", N_l)
    print("\t\t\t--> Achievable RV precision is ", rms_ms/np.sqrt(N_l))
    print("\t\t\tVelocity of ThAr w/r to solution provided is ", (1e-6*p1)*299792458.0)

    #for oo in np.unique(orders[I]):
    #   III = np.where(orders[I]==oo)[0]
    #   plot(wavelengths[I][III],residuals[III],'.')
    #show()

    return p1, pix_centers[I], orders[I], wavelengths[I], I, rms_ms, residuals

def Calculate_chebs(x,m, order0=89, ntotal=70,npix=2048.,Inverse=False,nx=5,nm=6):
    #normalize x
    medp = float(int(0.5*npix))
    x_norm = (x-medp) / medp

    if (Inverse):
        # invert and normalize m
        im = 1.0 / m
        im_max = 1.0 / float(order0)
        im_min = 1.0 / float(order0+ntotal)
    else:
        im = m
        im_max = float(order0+ntotal)
        im_min = float(order0)

    delta = 0.5*(im_max - im_min)
    m_norm = (im-im_min-delta)/delta

    coefs = []
    for i in range(nx):
        coefs.append(sp.chebyt(i+1)(x_norm))
    for i in range(nm):
        coefs.append(sp.chebyt(i+1)(m_norm))
    """
    u = sp.chebyt(1)(x_norm)
    u2 = sp.chebyt(2)(x_norm)
    u3 = sp.chebyt(3)(x_norm)
    u4 = sp.chebyt(4)(x_norm)
    v  = sp.chebyt(1)(m_norm)
    v2 = sp.chebyt(2)(m_norm)
    v3 = sp.chebyt(3)(m_norm)
    v4 = sp.chebyt(4)(m_norm)
    v5 = sp.chebyt(5)(m_norm)
    v6 = sp.chebyt(6)(m_norm)
    return u,u2,u3,u4,v,v2,v3,v4,v5,v6
    """
    return np.array(coefs)

def Joint_Polynomial_Cheby(p,chebs,nx,nm):
    """
    Evaluates a *product* of Chebyshev polynomials in x and m
    Polynomial is tailored.
    """
    xvec = chebs[:nx]
    mvec = chebs[nx:]
    ret_val = p[0]
    k=1
    for i in range(nx):
        ret_val += p[k]*xvec[i]
        k+=1
    for i in range(nm):
        ret_val += p[k]*mvec[i]
        k+=1

    if nx >= nm:
        for i in range(nx):
            for j in range(min(nx-i,nm)):
                ret_val += p[k]*xvec[i]*mvec[j]
                k+=1
    else:
        for j in range(nm):
            for i in range(min(nm-j-1,nx)):
                ret_val += p[k]*xvec[i]*mvec[j]
                k+=1
    return ret_val

def fp_base(f,n=3):
    of = f.copy()
    x = np.arange(len(f))
    I = np.where(f==0)[0]
    x = np.delete(x,I)
    f = np.delete(f,I)
    f1 = np.hstack((f[1:],f[0]))
    f2 = np.hstack((f[-1],f[:-1]))
    I = np.where((f<f1)&(f<f2))[0]
    xmin = x[I]
    fmin = f[I]

    coef = np.polyfit(xmin,fmin,n)
    res = fmin - np.polyval(coef,xmin)
    rms = np.sqrt(np.mean(res**2))
    I = np.where(np.absolute(res)>3*rms)[0]
    cond = True
    if len(I)==0:
        cond = False
    while cond:
        #print I
        iw = np.argmax(res**2)
        xmin = np.delete(xmin,iw)
        fmin = np.delete(fmin,iw)
        coef = np.polyfit(xmin,fmin,n)
        res = fmin - np.polyval(coef,xmin)
        rms = np.sqrt(np.mean(res**2))
        I = np.where(res>3*rms)[0]
        if len(I)==0:
            cond = False

    ox = np.arange(len(of))
    tbase = np.polyval(coef,ox)
    I = np.where(tbase<0)[0]
    tbase[I] = 0.
    ret = of - tbase
    I = np.where(of==0)[0]
    ret[I] = 0.
    return ret

def ccf_fp(fp,fpr,p1,order,order0=89,ntotal=70,npix=2048,Inv=True,nx=5,nm=6):
    def gaufp(p,v):
        retval = np.exp(-((v-p[0])**2)/(2*p[1]**2))
        return retval
    errgaufp = lambda p,c,v: np.ravel( (gaufp(p,v)-c) )
    pix_centers = np.arange(len(fp))
    chebs = Calculate_chebs(pix_centers, np.zeros(len(fp))+order+order0,order0=order0,ntotal=ntotal,npix=npix,Inverse=Inv,nx=nx,nm=nm)
    wav = (1.0/float(order+order0)) * Joint_Polynomial_Cheby(p1,chebs,nx,nm)
    #plt.plot(wav,fpr)
    #plt.show()
    vels = np.arange(-10000.,10000.,30.)
    ccf = []
    for v in vels:
        twav = wav*(1+v/299792458.)
        tck = interpolate.splrep(twav,fp,k=3)
        tfp = interpolate.splev(wav,tck)[500:-500]
        tfp /= np.sum(tfp)
        tfpr = fpr[500:-500]
        tfpr /= np.sum(tfpr)
        ccf.append(np.sum(tfpr*tfp))
    ccf = np.array(ccf)
    try:
        am = np.argmax(ccf)
        imin = np.argmin(ccf[:am])
        rmin = am + np.argmin(ccf[am:])
        vels = vels[imin:rmin+1]
        ccf  = ccf[imin:rmin+1]
        ccf -= .5*(ccf[0]+ccf[-1])
        ccf /= ccf.max()

        pg, success = scipy.optimize.leastsq(errgaufp, [0.,10.], args=(ccf, vels))

        #plt.plot(vels,ccf)
        #plt.plot(vels,gaufp(pg,vels))
        #plt.show()
        return pg[0]

    except:
        return -999.

def simbad_query_obname(obname):
    (th,tfile) = tempfile.mkstemp(prefix='CP', text=True)
    tf = open(tfile,'w')
    tf.write("output console=off\n")
    tf.write("output script=off\n")
    tf.write("output error=merge\n")
    tf.write("set limit 1\n")
    tf.write("format object fmt1 \"%IDLIST(1) | %OTYPELIST(S) | %SP(S)\"\n")
    tf.write("result full\n")
    tf.write("query id %s\n" % ( obname ) )
    tf.close()
    values = [
        ("scriptFIle", (pycurl.FORM_FILE, tfile))
    ]
    output = StringIO()
    c = pycurl.Curl()

    c.setopt(pycurl.URL, "http://simbad.u-strasbg.fr/simbad/sim-script")
    c.setopt(c.HTTPPOST, values)
    c.setopt(pycurl.WRITEFUNCTION, output.write)
    c.perform()
    c.close()

    result = output.getvalue()
    lines = result.split('\n')
    result = lines[len(lines)-3]
    query_success = False
    if (len(result) == 0):
        query_success = False
    elif (result.count('No known') > 0) or (result.count('not')>0) or (result.count('Unrecogniezd')>0):
        query_success = False
    else:
        query_success = True
        if len(result.split('|'))>=3:
            sp_type_query = result.split('|')[2]
        else:
            query_success = False
            sp_type_query = 'None'
    if not query_success:
        sp_type_query = None
    os.remove(tfile)
    return query_success, sp_type_query

def simbad_query_coords(ra,dec):
    (th,tfile) = tempfile.mkstemp(prefix='CP', text=True)
    tf = open(tfile,'w')
    tf.write("output console=off\n")
    tf.write("output script=off\n")
    tf.write("output error=merge\n")
    tf.write("set limit 1\n")
    tf.write("format object fmt1 \"%IDLIST(1) | %OTYPELIST(S) | %SP(S)\"\n")
    tf.write("result full\n")
    tf.write("set radius 5s\n")
    tf.write("query coo " + str(ra)+ " " + str(dec) + "\n")
    tf.close()
    values = [
        ("scriptFIle", (pycurl.FORM_FILE, tfile))
    ]
    output = StringIO()
    c = pycurl.Curl()

    c.setopt(pycurl.URL, "http://simbad.u-strasbg.fr/simbad/sim-script")
    c.setopt(c.HTTPPOST, values)
    c.setopt(pycurl.WRITEFUNCTION, output.write)
    c.perform()
    c.close()

    result = output.getvalue()
    lines = result.split('\n')
    result = lines[len(lines)-3]
    query_success = False

    if (result.count('No') > 0):
        query_success = False
    else:
        query_success = True
        try:
            sp_type_query = result.split('|')[2]
        except:
            sp_type_query = 'None'
            query_success = False
    if not query_success:
        sp_type_query = 'None'
    os.remove(tfile)
    return query_success, sp_type_query

def Lines_mBack(thar, sd, thres_rel=3, line_w=10):
    """
    Given an extracted ThAr order, return a version where the background has been removed
    """
    L = np.where(sd > 0)

    # initial background estimate, simple 25% percentile on a 100-pixel filter
    bkg = np.zeros( len(sd) )

    bkg[L] = scipy.signal.order_filter(thar[L],np.ones(101),25)

    d = thar - bkg
    lines = FindLines_simple_sigma(d,sd, thres=thres_rel)
    # Now, mask these lines
    mask = np.ones( len(sd) )
    for kk in lines:
        mask[kk-line_w:kk+line_w+1] = 0

    # New, final background estimnate
    X = np.array( range( len( d ) ) )
    K = np.where((sd > 0) & (mask > 0))
    if len(K[0])>0:
        bkg = np.zeros( len(sd) )
        bkg_T = lowess(thar[K].astype('double'), X[K],frac=0.2,it=3,return_sorted=False)
        tck1 = scipy.interpolate.splrep(X[K],bkg_T,k=1)
        bkg[L] = scipy.interpolate.splev(X[L],tck1)
        return bkg
    else:
        return np.zeros( len(sd) )

def FindLines_simple_sigma(d,sd,thres=3):
    """
    Given an array, find lines above a sigma-threshold
    """
    L = np.where( (d > (thres*sd)) & (d > 0) )
    lines = []
    for i in range(np.shape(L)[1]):
        j = L[0][i]
        if ((d[j] > d[j-1]) & (d[j-1] > d[j-2]) & (d[j] > d[j+1]) & (d[j+1] > d[j+2])):
            lines.append(j)

    return lines

def XC_ThAr(thar_ref, thar_comp, pixel_span):
    """
    Calculates the cross-correlation of two ThAr spectra
    """

    I1 = np.where( thar_ref != 0 )
    I2 = np.where( thar_comp != 0 )
    x1 = np.maximum( np.min( I1 ), np.min( I2) ) + pixel_span
    x2 = np.minimum( np.max( I1 ), np.max( I2) ) - pixel_span + 1

    N = int(2*pixel_span +1)
    shifts = -pixel_span + np.arange( N )
    XCor = np.zeros( N )

    mu1 = np.average( thar_ref[I1] )
    sd1 = np.sqrt( np.var( thar_ref[I1] ) )
    mu2 = np.average( thar_comp[I2] )
    sd2 = np.sqrt( np.var( thar_comp[I2] ) )

    norm = np.shape( thar_ref[I1] )[0]

    #print x1, x2
    for k in range(N):
        l1 = x1 + shifts[k]
        l2 = x2 + shifts[k]
        XCor[k] = np.sum( (thar_ref[x1:x2] - mu1)* (thar_comp[l1:l2] - mu2) / (sd1 * sd2)  ) / (norm-1)

    return XCor, shifts

def XC_Gau_Fit(X,Y,back_lag=5, usemin=1):
    """
    Fits a Gaussian to a XC
    """

    fitfunc = lambda p, x: p[0]*CorGaussian(x,p[1],p[2]) + p[3]
    errfunc = lambda p, x, y: np.ravel( (fitfunc(p,x) - y) )

    n = 4
    p0 = np.zeros( n )
    p0[3] = 0.5*( np.median(Y[:back_lag]) + np.median(Y[-back_lag:]) )
    p0[2] = 5
    if (usemin == 1):
        p0[1] = X[ np.argmin( Y ) ]
        p0[0] = np.min(Y) - p0[3]
    else:
        p0[1] = X[ np.argmax( Y ) ]
        p0[0] = np.max(Y) - p0[3]

    p1, success = scipy.optimize.leastsq(errfunc,p0, args=(X,Y))

    return p1

def get_mask(sp_type_query,T_eff, query_success):
    if (query_success):
        if (sp_type_query.count('G') > 0):
            print("Using G mask")
            sp_type = 'G2'
        elif (sp_type_query.count('K') > 0):
            print("Using K mask")
            sp_type = 'K5'
        elif (sp_type_query.count('M') > 0):
            print("Using M mask")
            sp_type = 'K5'
        elif (sp_type_query.count('F') > 0) or (sp_type_query.count('A') > 0):
            print("Rather hot star, using G mask but this target is probably earlier than that..")
            sp_type = 'G2'
        else:
            print("Probable problem, spectral type not in the books, using G-mask")
            sp_type = 'G2'
    else:
        if (T_eff >= 4300) and (T_eff < 5350):
            print("Using K mask")
            sp_type = 'K5'
        elif (T_eff >= 5350) and (T_eff < 6100):
            print("Using G mask")
            sp_type = 'G2'
        elif (T_eff >= 6100):
            print("Rather hot star, using G mask but this target is probably earlier than that..")
            sp_type = 'G2'
        else:
            print("Using M mask")
            sp_type = 'K5'
    return sp_type

def get_cont_single(W,F,E,nc=3,ll=3,lu=3,span=10,fact=3.,frac=0.3):
    I = np.where(F==0)[0]
    wav,flx,err = np.delete(W,I),np.delete(F,I),1./np.sqrt(np.delete(E,I))
    I = np.where(np.isnan(flx)==True)[0]
    wav,flx,err = np.delete(wav,I),np.delete(flx,I),np.delete(err,I)

    i = 0
    rw,re,rd = [],[],[]
    good_w,good_f = [],[]
    cond = True
    while cond:
        while i < len(wav):
            try:
                w,f,e = wav[i:int(i+span)],flx[i:int(i+span)],err[i:int(i+span)]
                dev = np.sqrt(np.var(f))
                rw.append(np.mean(w))
                re.append(np.median(e))
                rd.append(dev)
                if dev<fact*np.median(err) and np.median(f)>0.5*np.median(flx):
                    good_w.append(w[int(0.5*span)])
                    good_f.append(f[int(0.5*span)])
            except:
                None
            i+=1

        #print len(good_w),0.1*len(wav),len(wav)
        if len(good_w)>0.2*len(wav):
            cond = False
        else:
            i = 0
            rw,re,rd = [],[],[]
            good_w,good_f = [],[]
            fact +=1
        if fact>20:
            return np.array([0,np.max(scipy.signal.medfilt(F,21))])

    gw,gf = np.array(good_w),np.array(good_f)

    rw,re,rd = np.array(rw),np.array(re),np.array(rd)
    #plot(rw,rd,'ro')
    #plot(rw,fact*re,'bo')
    #show()
    lori = len(gw)
    cond = True
    counter = 0
    while cond:
        coef = np.polyfit(gw,gf,nc)
        res = gf - np.polyval(coef,gw)
        I = np.where(res>0)[0]
        dev = np.mean(res[I])
        IO1 = np.where(res> lu*dev)[0]
        IO2 = np.where(res<-ll*dev)[0]
        IO = np.sort(np.hstack((IO1,IO2)))
        if len(IO)==0 or len(gw)<frac*lori:
            cond = False
        else:
            gw = np.delete(gw,IO)
            gf = np.delete(gf,IO)
        counter +=1
    return coef



def get_cont(W,F,nc=3,ll=1.,lu = 5.,frac=0.1,window=21):

    blns = [[6755,6769],[6530,6600],[4840,4880],[4320,4360],[4085,4120],[3950,3990],[3880,3910],[3825,3850]]
    flx = F[0]
    wav = W[0]
    I = np.where(flx!=0)[0]
    wav,flx = wav[I],flx[I]

    for i in range(F.shape[0]-1):
        f = F[i+1]
        w = W[i+1]
        I = np.where(f!=0)[0]
        w,f = w[I],f[I]
        I = np.where(w>=wav[0])[0]

        if len(I)>5:
            J = np.where(wav<w[-1])[0]
            #print w[I]
            #print f[I]
            #plot(w[I],f[I],'r',linewidth=2.0)
            #plot(wav,flx,'k')
            #plot(w,f,'g')
            #show()
            tck = scipy.interpolate.splrep(w[I],f[I],k=3)
            nf = scipy.interpolate.splev(wav[J],tck)
            flx[J] = .5*(flx[J]+nf)
        I = np.where(w<wav[0])[0]
        wav = np.hstack((w[I],wav))
        flx = np.hstack((f[I],flx))

    #plot(wav,flx)
    #show()
    #medf = scipy.signal.medfilt(flx,31)
    #maxF = medf.max()
    #flx /= medf.max()
    wavo,flxo = wav.copy(),flx.copy()
    for lns in blns:
        I = np.where((wav>lns[0])&(wav<lns[1]))[0]
        wav,flx = np.delete(wav,I),np.delete(flx,I)
    for i in range(F.shape[0]):
        if i == 0:
            wi = W[i+2,0]
            wf = W[i,-1]
        elif i == 1:
            wi = W[i+2,0]
            wf = W[i-1,-2]
        elif i == W.shape[0] - 1:
            wi = W[i,0]
            wf = W[i-2,-1]
        elif i == W.shape[0] - 2:
            wi = W[i+1,0]
            wf = W[i-2,-1]
        else:
            wi = W[i+2,0]
            wf = W[i-2,-1]

        IM = np.where((wav>wi) & (wav<wf))[0]
        tw,tf = wav[IM],scipy.signal.medfilt(flx[IM],21)
        #plot(tw,flx[IM])
        JJJ = np.where((wav>W[i,0])&(wav<W[i,-1]))[0]
        #plot(wav[JJJ],flx[JJJ])

        ori = len(tw)
        coef = np.polyfit(tw,tf,nc)
        #show()
        #print gfds
        res = tf - np.polyval(coef,tw)
        IU = np.where(res>0)[0]
        dev = np.mean(res[IU])
        #plot(tw,tf)
        #plot(tw, np.polyval(coef,tw))
        #plot(tw, np.polyval(coef,tw) +lu*dev)
        #plot(tw, np.polyval(coef,tw) -ll*dev)
        #show()
        #print gfds
        I = np.where((res<lu*dev)&(res>-ll*dev))[0]
        cond = True

        while cond:
            tw,tf = tw[I],tf[I]
            coef = np.polyfit(tw,tf,nc)
            res = tf - np.polyval(coef,tw)
            IU = np.where(res>0)[0]
            dev = np.mean(res[IU])
            I = np.where((res<lu*dev)&(res>-ll*dev))[0]
            J1 = np.where(res>=lu*dev)[0]
            J2 = np.where(res<=-ll*dev)[0]

            if (len(J1)==0 and len(J2)==0) or len(tw)<frac*ori:
                cond = False

        #plot(tw, np.polyval(coef,tw))
        #plot(tw, np.polyval(coef,tw) +lu*dev)
        #plot(tw, np.polyval(coef,tw) -ll*dev)
        #show()
        #print gfds
        if i == 0:
            coefs = coef
        else:
            coefs = np.vstack((coefs,coef))

    #coefs /= maxF
    return coefs


def get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,RA,DEC):
    mephem = ephem.Moon()
    mephem.compute(gobs)
    pfm = ephem.previous_full_moon(gobs.date)
    nfm = ephem.next_full_moon(gobs.date)
    pnm = ephem.previous_new_moon(gobs.date)
    nnm = ephem.next_new_moon(gobs.date)
    if gobs.date - pnm < gobs.date - pfm:
        moon_state = 'crescent'
        lunation   = (gobs.date-pnm)/(nfm-pnm)
    else:
        moon_state = 'waning'
        lunation   = 1. - (gobs.date-pfm)/(nnm-pfm)

    # now compute the radial velocity sun-moon
    ram  = Mcoo['ra'][0]
    decm = Mcoo['dec'][0]
    Svx  = Sp['x_rate'][0] - Mp['x_rate'][0]
    Svy  = Sp['y_rate'][0] - Mp['y_rate'][0]
    Svz  = Sp['z_rate'][0] - Mp['z_rate'][0]
    Sx   = Sp['x'][0] - Mp['x'][0]
    Sy   = Sp['y'][0] - Mp['y'][0]
    Sz   = Sp['z'][0] - Mp['z'][0]
    sun_moon = (Sx*Svx + Sy*Svy +Sz*Svz)/np.sqrt(Sx**2 + Sy**2 +Sz**2)  # This is the radial velocity between the sun and the moon
    # now compute the radial velocity observer-moon
    moon_obs = res['frac'][0]* 2.99792458E5
    rasep    = np.sqrt((RA - ram*360./24.)**2)
    decsep   = np.sqrt((DEC - decm)**2)
    if rasep > 180:
        rasep = 360 - rasep
    if decsep > 180:
        decsep = 360 -decsep
    moonsep2 = np.sqrt( (rasep)**2 + (decsep)**2 )
    moonvel =  sun_moon + moon_obs
    return lunation,moon_state,moonsep2,moonvel

def calc_bss(vels, xc_av, bot_i=0.1, bot_f=0.4, top_i=0.6, top_f=0.85):
        #plot(vels, xc_av)
        #show()
        #dic = pickle.load(open('/data/echelle/coralie/red/20090908/CORALIE.2009-09-09T07:09:49.WASP-18_XC_G2.pkl','r'))
        #dic = pickle.load(open('/data/echelle/coralie/red/20121110/CORALIE.2012-11-11T08:38:58.HD72673_XC_G2.pkl','r'))
        #vels = dic['vels']
        #xc_av = dic['xc_av']

    xc_av = xc_av - xc_av.min()

    min_index = np.argmin(xc_av) #index with near minimum of CCF

    lim_inf, lim_sup = 0, len(vels)-1

    i = min_index
    while i >= 1:
        if xc_av[i] >= xc_av[i-1]:
            lim_inf = i
            break
        i-=1

    i = min_index
    while i <= len(vels)-2:
        if xc_av[i] >= xc_av[i+1]:
            lim_sup = i
            break
        i+=1

    if xc_av[lim_inf] < xc_av[lim_sup]:
        xc_av = xc_av / xc_av[lim_inf]
    else:
        xc_av = xc_av / xc_av[lim_sup]
    #print xc_av[lim_inf:min_index+1][::-1]
    #print vels[lim_inf:min_index+1][::-1]
    #plot(xc_av[lim_inf:min_index+1][::-1], vels[lim_inf:min_index+1][::-1])
    #show()
    """
    finf = open('/home/rabrahm/testinf.txt','w')
    for i in range(len(xc_av[lim_inf:min_index+1][::-1])):
            finf.write(str(vels[lim_inf:min_index+1][::-1][i])+'\t'+str(xc_av[lim_inf:min_index+1][::-1][i])+'\n')
    finf.close()

    fsup = open('/home/rabrahm/testsup.txt','w')
    for i in range(len(xc_av[min_index:lim_sup+1])):
            fsup.write(str(vels[min_index:lim_sup+1][i])+'\t'+str(xc_av[min_index:lim_sup+1][i])+'\n')
    fsup.close()
    """
    if len(xc_av[lim_inf:min_index+1][::-1]) <= 3 or len(xc_av[min_index:lim_sup+1])<=3:
        stat = 0
    else:
        stat = 1

    if stat == 1:
        Ii = np.argsort(xc_av[lim_inf:min_index+1])
        Is = np.argsort(xc_av[min_index:lim_sup+1])
        tck_inf = interpolate.splrep( xc_av[lim_inf:min_index+1][Ii], vels[lim_inf:min_index+1][Ii], k=3, s=0 )
        tck_sup = interpolate.splrep( xc_av[min_index:lim_sup+1][Is], vels[min_index:lim_sup+1][Is], k=3, s=0 )

        if xc_av[lim_inf] < xc_av[lim_sup]:
            xc_vec = np.arange(xc_av[min_index],xc_av[lim_inf],0.01)
        else:
            xc_vec = np.arange(xc_av[min_index],xc_av[lim_sup],0.01)

        interp_vels_inf = interpolate.splev(xc_vec, tck_inf, der=0)
        interp_vels_sup = interpolate.splev(xc_vec, tck_sup, der=0)
        der_interp_vels_inf = interpolate.splev(xc_vec, tck_inf, der=1)
        der_interp_vels_sup = interpolate.splev(xc_vec, tck_sup, der=1)
        #print interp_vels_inf
        #print interp_vels_sup
        bisect = 0.5 * (interp_vels_inf + interp_vels_sup)

        I_bottom = np.where((xc_vec > bot_i) & (xc_vec < bot_f))[0]
        I_top = np.where((xc_vec > top_i) & (xc_vec < top_f))[0]
        span = bisect[I_top].mean() - bisect[I_bottom].mean()
        der_bottom = 0.5*(np.absolute(der_interp_vels_inf[I_bottom].mean()) + np.absolute(der_interp_vels_sup[I_bottom].mean()))
        der_top = 0.5*(np.absolute(der_interp_vels_inf[I_top].mean()) + np.absolute(der_interp_vels_sup[I_top].mean()))
        #print der_bottom, der_top
        #print span
        #plot(vels, xc_av)
        #plot(bisect, xc_vec,'r.')
        #show()
        I2 = np.where( (xc_vec > 0.4) & (xc_vec < 0.81) )[0]
        aju = np.polyfit(xc_vec[I2], bisect[I2],1)
        slope = aju[0]
    else:
        span = -999.0
        der_bottom = 0
        der_top = 0
        slope = 0
    return span,der_bottom,der_top, slope, stat

def calc_bss2(vels,xc,coef, bot_i=0.15, bot_f=0.4, top_i=0.6, top_f=0.9, dt=0.01):
    try:

        I1 = np.where((vels>coef[1]-3*coef[2]) & (vels<coef[1]) )[0]
        I2 = np.where((vels<coef[1]+3*coef[2]) & (vels>coef[1]) )[0]
        I3 = np.where(vels<coef[1]-4*coef[2])[0]
        I4 = np.where(vels>coef[1]+4*coef[2])[0]
        I = np.hstack((I3,I4))
        base = np.median(xc[I])

        xc = base - xc
        xc /= xc.max()


        v1,x1 = vels[I1],xc[I1]
        v2,x2 = vels[I2],xc[I2]
        #plot(v1,x1)
        #plot(v2,x2)
        #show()
        dp = top_f
        vect = []
        while dp >= top_i:
            lb = np.where(x1>dp)[0][0]
            m = (v1[lb] - v1[lb-1])/(x1[lb]-x1[lb-1])
            n = v1[lb] - m*x1[lb]
            bs1 = m*dp+n

            lb = np.where(x2>dp)[0][-1]
            m = (v2[lb] - v2[lb+1])/(x2[lb]-x2[lb+1])
            n = v2[lb] - m*x2[lb]
            bs2 = m*dp+n
            vect.append(0.5*(bs2+bs1))
            dp-=dt
        vect = np.array(vect)

        dp = bot_f
        vecb = []
        while dp >= bot_i:

            lb = np.where(x1>dp)[0][0]
            m = (v1[lb] - v1[lb-1])/(x1[lb]-x1[lb-1])
            n = v1[lb] - m*x1[lb]
            bs1 = m*dp+n

            lb = np.where(x2>dp)[0][-1]
            m = (v2[lb] - v2[lb+1])/(x2[lb]-x2[lb+1])
            n = v2[lb] - m*x2[lb]
            bs2 = m*dp+n
            vecb.append(0.5*(bs2+bs1))
            dp-=dt
        vecb = np.array(vecb)

        return np.median(vecb) - np.median(vect)
    except:
        return -999.0

def gauss_samp(p,l,mu):
    A,B,sig = p[0],p[1],p[2]
    px = np.arange(len(l)).astype('float')
    tck = interpolate.splrep(px,l,k=3)
    pxs = np.arange(-0.5,len(l)-0.5,0.01)
    I = np.where(pxs<0)[0]
    pxs[I] *= -1.
    I = np.where(pxs>=len(l)-1)[0]
    temp = -np.arange(0.,0.5,0.01)
    pxs[I] = temp + len(l) - 1
    ls = interpolate.splev(pxs,tck)
    G = B + A * np.exp(-.5*((ls - mu)/sig)**2)
    pxs = pxs.reshape((len(l),100))
    ls = ls.reshape((len(l),100))
    G = G.reshape((len(l),100))
    ret = np.mean(G,axis=1)
    return ret

def err_samp(p,y,x,mu):
    return np.ravel( (gauss_samp(p,x,mu) - y) )

def convolve(wav,flx,R):
    devs = wav/(2.35482*R)
    lims1 = wav - 3*devs
    lims2 = wav + 3*devs
    NF = []
    for i in range(len(wav)):
        I = np.where((wav>lims1[i]) & (wav<lims2[i]))[0]
        W = wav[I]
        F = flx[I]
        G = np.exp(-(wav[I] - wav[i])**2/(0.5*devs[i]**2))
        G = G / np.sum(G)
        NF.append(np.sum(F*G))
    NF = np.array(NF)
    return NF

def calc_bss(vels, xc_av, bot_i=0.1, bot_f=0.4, top_i=0.6, top_f=0.85):

    xc_av = xc_av - xc_av.min()

    min_index = np.argmin(xc_av) #index with near minimum of CCF

    lim_inf, lim_sup = 0, len(vels)-1

    i = min_index
    while i >= 1:
        if xc_av[i] >= xc_av[i-1]:
            lim_inf = i
            break
        i-=1

    i = min_index
    while i <= len(vels)-2:
        if xc_av[i] >= xc_av[i+1]:
            lim_sup = i
            break
        i+=1

    if xc_av[lim_inf] < xc_av[lim_sup]:
        xc_av = xc_av / xc_av[lim_inf]
    else:
        xc_av = xc_av / xc_av[lim_sup]
    #print xc_av[lim_inf:min_index+1][::-1]
    #print vels[lim_inf:min_index+1][::-1]
    #plot(xc_av[lim_inf:min_index+1][::-1], vels[lim_inf:min_index+1][::-1])
    #show()
    """
    finf = open('/home/rabrahm/testinf.txt','w')
    for i in range(len(xc_av[lim_inf:min_index+1][::-1])):
            finf.write(str(vels[lim_inf:min_index+1][::-1][i])+'\t'+str(xc_av[lim_inf:min_index+1][::-1][i])+'\n')
    finf.close()

    fsup = open('/home/rabrahm/testsup.txt','w')
    for i in range(len(xc_av[min_index:lim_sup+1])):
            fsup.write(str(vels[min_index:lim_sup+1][i])+'\t'+str(xc_av[min_index:lim_sup+1][i])+'\n')
    fsup.close()
    """
    if len(xc_av[lim_inf:min_index+1][::-1]) <= 3 or len(xc_av[min_index:lim_sup+1])<=3:
        stat = 0
    else:
        stat = 1

    if stat == 1:
        Ii = np.argsort(xc_av[lim_inf:min_index+1])
        Is = np.argsort(xc_av[min_index:lim_sup+1])
        tck_inf = interpolate.splrep( xc_av[lim_inf:min_index+1][Ii], vels[lim_inf:min_index+1][Ii], k=3, s=0 )
        tck_sup = interpolate.splrep( xc_av[min_index:lim_sup+1][Is], vels[min_index:lim_sup+1][Is], k=3, s=0 )

        if xc_av[lim_inf] < xc_av[lim_sup]:
            xc_vec = np.arange(xc_av[min_index],xc_av[lim_inf],0.01)
        else:
            xc_vec = np.arange(xc_av[min_index],xc_av[lim_sup],0.01)

        interp_vels_inf = interpolate.splev(xc_vec, tck_inf, der=0)
        interp_vels_sup = interpolate.splev(xc_vec, tck_sup, der=0)
        der_interp_vels_inf = interpolate.splev(xc_vec, tck_inf, der=1)
        der_interp_vels_sup = interpolate.splev(xc_vec, tck_sup, der=1)
        #print interp_vels_inf
        #print interp_vels_sup
        bisect = 0.5 * (interp_vels_inf + interp_vels_sup)

        I_bottom = np.where((xc_vec > bot_i) & (xc_vec < bot_f))[0]
        I_top = np.where((xc_vec > top_i) & (xc_vec < top_f))[0]
        span = bisect[I_top].mean() - bisect[I_bottom].mean()
        der_bottom = 0.5*(np.absolute(der_interp_vels_inf[I_bottom].mean()) + np.absolute(der_interp_vels_sup[I_bottom].mean()))
        der_top = 0.5*(np.absolute(der_interp_vels_inf[I_top].mean()) + np.absolute(der_interp_vels_sup[I_top].mean()))
        #print der_bottom, der_top
        #print span
        #plot(vels, xc_av)
        #plot(bisect, xc_vec,'r.')
        #show()
        I2 = np.where( (xc_vec > 0.4) & (xc_vec < 0.81) )[0]
        aju = np.polyfit(xc_vec[I2], bisect[I2],1)
        slope = aju[0]
    else:
        span = -999.0
        der_bottom = 0
        der_top = 0
        slope = 0
    return span,der_bottom,der_top, slope, stat

def cor_thar(spec, span=10, filename='/data/echelle/ecpipe/DuPont/wavcals/',binning=1, di=0.1):

    f = open(filename).readlines()

    pixel_centers_0 = []

    for line in f:
        w = line.split()
        nlines = int(w[0])
        for j in range(nlines):
            pixel_centers_0.append(float(w[2*j+1])/binning)

    ml = array(pixel_centers_0) - 1
    mh = array(pixel_centers_0) + 1

    vi = -span
    vf = span

    NM = np.sqrt(np.add.reduce((mh - ml)**2))
    NS = np.sqrt(integrate.simps(spec*spec) )

    MP = np.arange(0,len(spec)-1+di,di)
    MV = np.zeros(len(MP))

    for lin in range(len(ml)):
        I = np.where((MP >= ml[lin]) & (MP <= mh[lin]))[0]
        MV[I] = 1.0

    tck = interpolate.splrep(np.arange(len(spec)),spec,k=3,s=0)
    SV = interpolate.splev(MP,tck,der=0)

    CCF = []
    DEL = []

    while vi <= vf:
        mp = MP.copy() + vi

        mv = MV.copy()

        I = np.where(mp < 0.0)[0]

        if len(I) > 0:

            comp = mv[I]
            mp = np.delete(mp,I)
            mv = np.delete(mv,I)
            mv = np.hstack([mv,comp])
        I = np.where(mp > float(len(spec)-1))[0]

        if len(I) > 0:

            comp = mv[I]
            mp = np.delete(mp,I)
            mv = np.delete(mv,I)
            mv = np.hstack([comp,mv])
        CCF.append(np.add.reduce(SV*mv)*di/(NM*NS))
        DEL.append(vi)

        vi = vi + di

    CCF = np.array(CCF)
    DEL = np.array(DEL)
    #plot(DEL,CCF,'ro')
    #show()
    I = np.argmax(CCF)
    #plot(DEL,CCF)
    return CCF[I],DEL[I]

def plot_CCF(xc_dict,moon_dict,path='XC.pdf'):

    rvels        = xc_dict['rvels']
    rxc_av       = xc_dict['rxc_av']
    rpred        = xc_dict['rpred']
    rxc_av_orig  = xc_dict['rxc_av_orig']
    rvel0_xc     = xc_dict['rvel0_xc']
    refvel       = xc_dict['refvel']
    vels         = xc_dict['vels']
    xc_av        = xc_dict['xc_av']
    XCmodelgau   = xc_dict['XCmodelgau']
    XCmodelgau_m = xc_dict['XCmodelgau_m']
    Ls2          = xc_dict['Ls2']
    Ls2_m        = xc_dict['Ls2_m']
    p1gau        = xc_dict['p1gau']
    p1gau_m      = xc_dict['p1gau_m']

    moon_state  = moon_dict['moon_state']
    moonsep     = moon_dict['moonsep']
    moonmatters = moon_dict['moonmatters']
    lunation    = moon_dict['lunation']
    texp        = moon_dict['texp']
    mephem      = moon_dict['mephem']

    pp = PdfPages(path)
    plot(rvels, rxc_av,'b.')
    plot(rvels, rpred,'g-')
    plot(rvels, rxc_av_orig,'r.')
    if refvel > rvels[0] and refvel < rvels[-1]:
        plot([refvel,refvel],[rxc_av.min(), rxc_av.max()],'r--')
    axvline(rvel0_xc,linestyle=':')
    axhline(1.0,linestyle='-')
    pp.savefig()

    # plot of XC function
    f1 = figure()
    ax1 = f1.add_subplot(111)
    ax1.plot(vels, xc_av,'b.', label='CCF')
    ax1.plot(vels[Ls2], XCmodelgau,'r-',label='Gaussian fit')
    if moonmatters:
        ax1.plot(vels[Ls2_m], XCmodelgau_m,'g-',label='Double Gaussian fit')

    if refvel > vels[0] and refvel < vels[-1]:
        ax1.axvline(refvel,linestyle='--',color='k', label = 'RV Moon')
    ax1.text(vels[0], xc_av.min(), moon_state + ' moon:\nseparation: '+str(round(moonsep,1))+' deg\nillumination: '+str(int(round(lunation,2)*100))+' %\naltitude: '+ str(mephem.alt)+'\nTexp = '+str(int(round(texp))) +' s', style='italic',fontsize=7)

    xlabel('Velocity (km/s)')
    ylabel('XC')
    ax1.axvline(p1gau[1],linestyle=':',color='r')
    ax1.axvline(p1gau_m[1],linestyle=':',color='g')
    ax1.axhline(1.0,linestyle='-')
    title('Average Cross-Correlation Function + Fit')
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles[::-1], labels[::-1],prop={'size':6})
    pp.savefig()
    pp.close()
    clf()
    close('all')

    pass

def get_disp(obname,reffile='reffile.txt'):
    disp = 0.
    try:
        f = open(reffile,'r')
        lines = f.readlines()
        f.close()
        found = False
        for line in lines:
            cos = line.split(',')
            if cos[0] == obname:
                disp = float(cos[7])
                found = True
                break
        if not found:
            print('\t\tWarning! There is no predefined dispersion of the CCF.')
    except:
        print('\t\tWarning! There is no predefined dispersion of the CCF.')

    return disp

def get_mask_reffile(obname,reffile='reffile.txt',base='../../xc_masks/'):
    xc_masks = [base+'G2.mas',\
                base+'K5.mas',\
                base+'M2.mas' ]
    sp_type = 'G2'

    try:
        f = open(reffile,'r')
        lines = f.readlines()
        f.close()

        found = False
        for line in lines:
            cos = line.split(',')
            msk = cos[6].strip()
            if cos[0] == obname and (msk=='G2' or msk=='K5' or msk=='M2'):
                sp_type = cos[6][:2]
                found = True
                break
        if not found:
            print('\t\tWarning! Target not found in reference mask file. Using default mask (G2)')
    except:
        print('\t\tWarning! Problem with reference mask file. Forcing to G2 mask')

    if sp_type == 'G2':
        xc_mask = xc_masks[0]
    elif sp_type == 'K5':
        xc_mask = xc_masks[1]
    else:
        xc_mask = xc_masks[2]

    return sp_type, xc_mask

def get_mask_query(sp_type_query,base='../../xc_masks/'):
    xc_masks = [base+'G2.mas',\
                base+'K5.mas',\
                base+'M2.mas' ]

    if (sp_type_query.count('G') > 0):
        print("Using G mask according to simbad")
        mask = xc_masks[0]
        sp_type = 'G2'
    elif (sp_type_query.count('K') > 0):
        print("Using K mask according to simbad")
        mask = xc_masks[1]
        sp_type = 'K5'
    elif (sp_type_query.count('M') > 0):
        print("Using M mask according to simbad")
        mask = xc_masks[2]
        sp_type = 'M5'
    elif (sp_type_query.count('F') > 0) or (sp_type_query.count('A') > 0):
        print("Rather hot star according to simbad, using G mask but this target is probably earlier than that..")
        mask = xc_masks[0]
        sp_type = 'G2'
    else:
        print("Probable problem, spectral type not in the books according to simbad, using G-mask")
        mask = xc_masks[0]
        sp_type = 'G2'
    return sp_type, mask

def get_mask_teff(T_eff,base='../../xc_masks/'):
    xc_masks = [base+'G2.mas',\
                base+'K5.mas',\
                base+'M2.mas' ]
    if (T_eff >= 4300) and (T_eff < 5350):
        print("Using K mask according to T_eff")
        mask = xc_masks[1]
        sp_type = 'K5'
    elif (T_eff >= 5350) and (T_eff < 6100):
        print("Using G mask according to T_eff")
        mask = xc_masks[0]
        sp_type = 'G2'
    elif (T_eff >= 6100):
        print("Rather hot star according to T_eff, using G mask but this target is probably earlier than that..")
        mask = xc_masks[0]
        sp_type = 'G2'
    else:
        print("Using M mask according to T_eff")
        mask = xc_masks[2]
        sp_type = 'M5'
    return sp_type, mask


def iau_cal2jd(IY,IM,ID):
    IYMIN = -4799.
    MTAB = np.array([ 31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
    J = 0
    if IY < IYMIN:
        J = -1
    else:
        if IM>=1 and IM <= 12:
            if IY%4 == 0:
                MTAB[1] = 29.
            else:
                MTAB[1] = 28.

            if IY%100 == 0 and IY%400!=0:
                MTAB[1] = 28.
            if ID < 1 or ID > MTAB[IM-1]:
                J = -3
            a = ( 14 - IM ) / 12
            y = IY + 4800 - a
            m = IM + 12*a -3
            DJM0 = 2400000.5
            DJM = ID + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045 - 2400001.
        else:
            J = -2
    return DJM0, DJM, J

def identify_order(thar, path, window=100,di=0.5):
    ccfs,shifts,ors = [],[],[]
    for o in range(len(thar)):
        tharo = thar[o] - scipy.signal.medfilt(thar[o],101)
        ccf_max, shift = cor_thar(tharo,filename=path,span=window,di=di)
        ccfs.append(ccf_max)
        shifts.append(shift)
        ors.append(o)
    ccfs,shifts,ors = np.array(ccfs),np.array(shifts),np.array(ors)
    im = np.argmax(ccfs)
    return ors[im],shifts[im]


def simbad_coords(obname,mjd):
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

    paths = ['http://simbad.u-strasbg.fr/simbad/sim-script','http://simbad.harvard.edu/simbad/sim-script']

    cond = True
    i = 0
    while cond:
        try:
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
            output = StringIO()
            c = pycurl.Curl()
            c.setopt(pycurl.URL, paths[i])
            c.setopt(c.HTTPPOST, values)
            c.setopt(pycurl.WRITEFUNCTION, output.write)
            c.perform()
        except:
            if i == 0:
                i = 1
            else:
                i=0
            print('Trying again to perform query to SIMBAD, changing to', paths[i])

        else:
            cond = False
    c.close()
    result = output.getvalue()
    lines = result.split('\n')
    info = lines[6].split('|')
    print(info)
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

def spec_ccf(sw,sf,mw,mf,vi,vf,dv):
    mf = mf -1
    mf = -mf
    tck = interpolate.splrep(mw,mf,k=1)
    v = vi
    retccf = []
    vels = []
    while v<=vf:
        swt = sw * (1 + v/299792.458)
        mft = interpolate.splev(swt,tck)
        mft -= np.mean(mft)
        sft = sf - np.mean(sf)
        #sft = sf.copy()
        retccf.append(np.sum(mft*sft)/np.sqrt(np.sum(mft**2)*np.sum(sft**2)))
        vels.append(v)
        v+=dv
    return np.array(vels),np.array(retccf)


def RVforFR(data,teff=6750,logg=4.5,feh=-1.0,vsini=40.,model_path='../data/COELHO_MODELS/R_40000b/'):
    def fitfunc(p,x):
        ret = p[3] + p[0] * np.exp(-.5*((x-p[1])/p[2])**2)
        return ret
    errfunc = lambda p,x,y: np.ravel( (fitfunc(p,x)-y) )

    if feh == -1:
        sfeh = 'm10'
    elif feh == -0.5:
        sfeh = 'm05'
    elif feh == 0:
        sfeh = 'p00'
    else:
        sfeh = 'p05'


    model = model_path+'/vsini_0.0/R_0.0_'+str(int(teff))+'_'+str(int(logg*10))+'_'+sfeh+'p00.ms.fits'
    sc = pyfits.getdata(model)
    hd = pyfits.getheader(model)
    mw = np.arange(len(sc))*hd['CDELT1']+hd['CRVAL1']
    #""""
    I = np.where((mw>6520)&(mw<6600))[0]
    sc[I] = 1.
    I = np.where((mw>5888)&(mw<5897))[0]
    sc[I] = 1.
    I = np.where((mw>4310)&(mw<4360))[0]
    sc[I] = 1.
    I = np.where((mw>4070)&(mw<4130))[0]
    sc[I] = 1.
    I = np.where((mw>3875)&(mw<3900))[0]
    sc[I] = 1.
    I = np.where((mw>3920)&(mw<3945))[0]
    sc[I] = 1.
    I = np.where((mw>3955)&(mw<3980))[0]
    sc[I] = 1.
    I = np.where(mw<3850)[0]
    sc[I] = 1.
    #"""

    mw = ToVacuum(mw)
    ccftot = []

    for i in range(data.shape[1]):
        scf = data[5,i]
        scw = data[0,i]
        J = np.where(scf!=0)[0]
        scw,scf = scw[J],scf[J]
        I = np.where((mw>scw[0]-100) & (mw<scw[-1]+100))
        tmf = pyasl.fastRotBroad(mw[I], sc[I], 0.5, vsini)
        ccv,ccf = spec_ccf(scw,scf,mw[I],tmf,-1000,1000,10.)
        ccf = np.array(ccf)
        if len(ccftot)==0:
            ccftot = ccf.copy()
        else:
            ccftot = np.vstack((ccftot,ccf))
        #plot(ccv,ccf/np.sum(ccf))
        #show()+
    ccftot = np.mean(ccftot,axis=0)
    ccftot += 1.


    p0 = [-(1-ccftot.min()),ccv[np.argmin(ccftot)],vsini,ccftot[0]]
    p1, success = scipy.optimize.leastsq(errfunc,p0, args=(ccv,ccftot))

    #plot(ccv,ccftot)
    #plot(ccv,fitfunc(p0,ccv))
    #plot(ccv,fitfunc(p1,ccv))
    #print p0
    #print p1
    #show()

    return p1,ccv,ccftot,fitfunc(p1,ccv)


def new_ccf(data,model_path):
    teffs = [6000,6250,6500,6750,7000]
    loggs = [3.0,3.5,4.0,4.5]
    fehs  = [-1.0,-0.5,0.0,0.5]
    rots  = [1.,10.,50.,100.,200.,300.]

    mods,pars = [],[]
    for tef in teffs:
        stef = str(int(tef))
        for log in loggs:
            slog = str(int(log*10))
            for feh in fehs:
                if feh == -1:
                    sfeh = 'm10'
                elif feh == -0.5:
                    sfeh = 'm05'
                elif feh == 0.0:
                    sfeh = 'p00'
                else:
                    sfeh = 'p05'
                mod = '/vsini_0.0/R_0.0_'+stef+'_'+slog+'_'+sfeh+'p00.ms.fits'
                mods.append(model_path+mod)
                if len(pars) == 0:
                    pars = np.array([tef,log,feh])
                else:
                    pars = np.vstack((pars,np.array([tef,log,feh])))
    mods = np.array(mods)
    mods = [model_path+'/vsini_0.0/R_0.0_6750_45_m10p00.ms.fits']
    modmin = -1
    mini   = 1000
    rotmin = -1
    for i in range(len(mods)):
        model = mods[i]
        par = pars[i]

        sc = pyfits.getdata(model)
        hd = pyfits.getheader(model)
        mw = np.arange(len(sc))*hd['CDELT1']+hd['CRVAL1']
        I = np.where((mw>6520)&(mw<6600))[0]
        sc[I] = 1.
        I = np.where((mw>5888)&(mw<5897))[0]
        sc[I] = 1.
        I = np.where((mw>4310)&(mw<4360))[0]
        sc[I] = 1.
        I = np.where((mw>4070)&(mw<4130))[0]
        sc[I] = 1.
        I = np.where((mw>3875)&(mw<3900))[0]
        sc[I] = 1.
        I = np.where((mw>3920)&(mw<3945))[0]
        sc[I] = 1.
        I = np.where((mw>3955)&(mw<3980))[0]
        sc[I] = 1.
        I = np.where(mw<3850)[0]
        sc[I] = 1.

        mw = ToVacuum(mw)

        for rot in rots:
            ccftot = []
            for i in range(data.shape[1]):
                scf = data[5,i]
                scw = data[0,i]
                J = np.where(scf!=0)[0]
                scw,scf = scw[J],scf[J]
                I = np.where((mw>scw[0]-100) & (mw<scw[-1]+100))
                tmf = pyasl.fastRotBroad(mw[I], sc[I], 0.5, rot)
                #plot(scw,scf)
                #plot(mw[I],tmf)
                #show()
                ccv,ccf = spec_ccf(scw,scf,mw[I],tmf,-1000,1000,10.)
                ccf = np.array(ccf)
                if len(ccftot)==0:
                    ccftot = ccf.copy()
                else:
                    ccftot += ccf
                #plot(ccv,ccf/np.sum(ccf))
                #show()
            print(model, rot, ccftot.min())
            if ccftot.min() < mini:
                mini = ccftot.min()
                modmin = par
                rotmin = rot
            #plot(ccv,ccftot)
    #show()
    print(modmin,rotmin,mini)
