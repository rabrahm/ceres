from __future__ import print_function
import matplotlib
matplotlib.use("Agg")
from astropy.io import fits as pyfits
from pylab import *
import scipy
import glob
import os
import sys
base = '../'
sys.path.append(base+"utils/GLOBALutils")
import GLOBALutils

def getObName(h):
    """
    Get name of object under consideration
    """
    obname = h[0].header['HIERARCH ESO OBS TARG NAME'].upper().replace(' ','')
    if 'ARCFILE' in h[0].header.keys():
        if h[0].header['ARCFILE'] == 'HARPS.2011-03-29T01:43:58.502.fits':
            obname = 'HD77338'
    return obname

def FileClassify(diri, log, mode='HARPS'):
    """

    Classifies all files in a directory and writes a night log of science images

    """

    # define output lists
    sim_sci        = []
    biases         = []
    ob_flats       = []
    co_flats       = []
    flats          = []
    ThAr_ref       = []
    ThAr_ref_dates = []
    obnames        = []
    exptimes       = []
    co_types       = []
    f = open(log,'w')

    #Do not consider the images specified in dir+badfiles.txt
    bad_files = []
    if os.access(diri+'bad_files.txt',os.F_OK):
        bf = open(diri+'bad_files.txt')
        linesbf = bf.readlines()
        for line in linesbf:
            bad_files.append(diri+line[:-1])
        bf.close()

    all_files = glob.glob(diri+"/HARPS*fits")
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
            if h[0].header['HIERARCH ESO INS MODE'] == mode:
                if h[0].header['HIERARCH ESO DPR TYPE'][:9] == 'STAR,WAVE' or h[0].header['HIERARCH ESO DPR TYPE'][:8] == 'STAR,SKY' or h[0].header['HIERARCH ESO DPR TYPE'][:9] == 'STAR,DARK':
                    co_type = h[0].header['HIERARCH ESO DPR TYPE'].split(',')[1]
                    if co_type == 'DARK':
                        co_type = 'SKY'
                    co_types.append(co_type)
                    sim_sci.append(archivo)
                    obname = getObName(h)
                    obnames.append( obname )
                    ra     = h[0].header['RA']
                    delta  = h[0].header['DEC']
                    airmass= h[0].header['HIERARCH ESO TEL AIRM START']
                    texp   = h[0].header['EXPTIME']
                    date   = h[0].header['DATE-OBS'][:10]
                    hour   = h[0].header['DATE-OBS'][11:]
                    exptimes.append( texp )
                    line = "%-15s %8s %10s %10s %8.2f %4.2f %8s %11s %s\n" % (obname, co_type, ra, delta, texp, airmass, date, hour, archivo)
                    f.write(line)
                elif h[0].header['HIERARCH ESO DPR TYPE'][:9] == 'BIAS,BIAS':
                    biases.append(archivo)
                elif h[0].header['HIERARCH ESO DPR TYPE'][:9] == 'LAMP,DARK':
                    ob_flats.append(archivo)
                elif h[0].header['HIERARCH ESO DPR TYPE'][:9] == 'DARK,LAMP':
                    co_flats.append(archivo)
                elif h[0].header['HIERARCH ESO DPR TYPE'][:9] == 'LAMP,LAMP':
                    flats.append(archivo)
                elif h[0].header['HIERARCH ESO DPR TYPE'] == 'WAVE,WAVE,THAR2':
                    ThAr_ref.append(archivo)
                    mjd, mjd0 = mjd_fromheader(h)
                    ThAr_ref_dates.append( mjd )

    f.close()

    return biases, flats, ob_flats, co_flats, ThAr_ref, sim_sci, ThAr_ref_dates, obnames, exptimes, co_types

def mjd_fromheader(h):
    """
    return modified Julian date from header
    """
    secinday = 24*3600.0

    datetu   = h[0].header['DATE-OBS'][:10]
    ut       = h[0].header['DATE-OBS'][11:]

    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[5:7]),int(datetu[8:10]))

    ut        = (float(ut[:2])*3600. + float(ut[3:5])*60. + float(ut[6:]))
    mjd_start = mjd + ut / secinday

    fraction = h[0].header['HIERARCH ESO INS DET1 TMMEAN']
    texp     = h[0].header['EXPTIME'] #sec

    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0

def MedianCombine(ImgList, zero='none'):
    """

    Median combine a list of images

    """

    n = len(ImgList)
    if n==0:
        raise ValueError("empty list provided!")

    h = pyfits.open(ImgList[0])
    d1 = h[1].data
    h1 = h[1].header
    d2 = h[2].data
    h2 = h[2].header

    d1 = OverscanTrim(d1)
    d2 = OverscanTrim(d2)

    if zero != 'none':
        z = pyfits.open(zero)[0]
        d1 -= z.data[:,:,0]
        d2 -= z.data[:,:,1]

    factor = 1.25
    if (n < 3):
        factor = 1

    ron1,gain1 = h1['HIERARCH ESO DET OUT1 RON'],h1['HIERARCH ESO DET OUT1 GAIN']
    ron2,gain2 = h2['HIERARCH ESO DET OUT1 RON'],h2['HIERARCH ESO DET OUT1 GAIN']
    ron1 = factor * ron1 / np.sqrt(n)
    ron2 = factor * ron2 / np.sqrt(n)
    if n>1:
        for i in range(n-1):
            td = pyfits.open(ImgList[i+1])
            if zero == 'none':
                d1 = np.dstack((d1,OverscanTrim(td[1].data)))
                d2 = np.dstack((d2,OverscanTrim(td[2].data)))

            else:
                d1 = np.dstack((d1,OverscanTrim(td[1].data)-z.data[:,:,0]))
                d2 = np.dstack((d2,OverscanTrim(td[2].data)-z.data[:,:,1]))
        d1 = np.median(d1,axis=2)
        d2 = np.median(d2,axis=2)

    return np.dstack((d1,d2)), np.array([ron1,ron2]), np.array([gain1,gain2])

def OverscanTrim(d):
    """
    Overscan correct and Trim a refurbished CORALIE image
    """
    # bias has no significant structure, so a single median suffices, I think
    # overscan = [0:49] [2097:2145]
    data = d[:,51:2099]
    ov1  = d[:,:50]
    ov2  = d[:,2098:]
    ov = np.vstack((ov1.transpose(),ov2.transpose()))
    overscan = np.median(ov,axis=0)
    overscan = np.array([overscan]*2048)
    overscan = overscan.transpose()
    newdata = data - overscan
    return newdata.transpose()
