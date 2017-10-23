import sys
from pylab import *

base = '../'
sys.path.append(base+"utils/Correlation")
sys.path.append(base+"utils/GLOBALutils")
sys.path.append(base+"utils/OptExtract")


baryc_dir= base+'utils/SSEphem/'
sys.path.append(baryc_dir)
ephemeris='DEc403'

import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ceres modules
import uvesutils
import correlation
import GLOBALutils
import Marsh

# other useful modules
import argparse
import ephem
import jplephem
from math import radians as rad
import pyfits
import pickle
import os
import scipy
import scipy.interpolate
from scipy import interpolate

import statsmodels.api as sm
lowess = sm.nonparametric.lowess


# Recive input parameters
parser = argparse.ArgumentParser()
parser.add_argument('directorio')
parser.add_argument('-avoid_plot', action="store_true", default=False)
parser.add_argument('-dirout',default='default')
parser.add_argument('-do_class', action="store_true", default=False)
parser.add_argument('-just_extract', action="store_true", default=False)
parser.add_argument('-npools', default=1)
parser.add_argument('-o2do',default='all')
parser.add_argument('-reffile',default='default')
parser.add_argument('-ofind', default='last')

args = parser.parse_args()
dirin            = args.directorio
avoid_plot       = args.avoid_plot
dirout           = args.dirout
DoClass          = args.do_class
JustExtract      = args.just_extract
npools           = int(args.npools)
object2do        = args.o2do
reffile          = args.reffile
stst            = args.ofind

if dirin[-1] != '/':
    dirin = dirin + '/'

if dirout == 'default':
    dirout = dirin[:-1]+'_red/'

if not os.access(dirout,os.F_OK):
    os.system('mkdir '+dirout)
if os.access(dirout+'proc',os.F_OK):
    os.system('rm -r '+dirout+'proc')
os.system('mkdir '+dirout+'proc')

f_res = open(dirout+'proc/'+'results.txt','w')

if reffile == 'default':
    reffile = dirin+'reffile.txt'

####### GLOBAL VARIABLES #####
## perhaps put into options ##
force_pre_process  = False
force_flat_extract = False
force_flat_nor     = False
force_thar_extract = False
force_thar_wavcal  = False
force_tharxc       = False	
force_sci_extract  = False
force_spectral_file_build = True
force_stellar_pars = False
force_bac          = False
dumpargon          = False
minlines_glob      = 400

Inverse_m          = True
use_cheby          = True
MRMS               = 100   # max rms in m/s, global wav solution

trace_degree       = 6
Marsh_alg          = 0
ext_aperture       = 21
sky_aperture       = 31
NSigma_Marsh       = 10
NCosmic_Marsh      = 10
S_Marsh            = 0.4
N_Marsh            = 3      # grado polinomio 
min_extract_col      = 0
min_extract_cols1    = [920,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
min_extract_cols2    = [1710,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
max_extract_col      = 4000
max_extract_cols1    = [3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895]
max_extract_cols2    = [3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3895,3390]

xbinning           = 1
ybinning           = 1

ext_aperture = int(float(ext_aperture)/float(ybinning)) 

ncoef_x   = 5
ncoef_m   = 6
npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2
ro0       = 90
order_gap = 2
order_dir  = base+"/uves/wavcals/"

#############################

print "\n\n\tUVES @VLT  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'

# file containing the log
log = dirout+'night.log'
biases, flats, orderref, ThAr_ref, sim_sci, ThAr_ref_dates, obnames, exptimes = uvesutils.FileClassify(dirin,log)

print '\tThis in the log of the night:\n'
f = open(log)
flines = f.readlines()
for line in flines:
	print '\t\t'+line[:-1]
print '\n'

if ((os.access(dirout+'Flat.fits',os.F_OK) == False) or \
    (os.access(dirout+'MasterBias.fits',os.F_OK) == False) or \
    (os.access(dirout+'trace.pkl',os.F_OK) == False)  or \
    (force_pre_process) ):
    print "\tNo previous pre-processing files or found"
    pre_process = 1
else:
    print "\tPre-processing files found, going straight to extraction"
    pre_process = 0

if (pre_process == 1):
    # median combine list of flats
    print "\t\tGenerating Master calibration frames..."

    MasterBias, RO_bias, GA_bias = uvesutils.MedianCombine(biases)
    hdu = pyfits.PrimaryHDU(MasterBias)
    if (os.access(dirout+'MasterBias.fits',os.F_OK)):
        os.remove(dirout+'MasterBias.fits')
    hdu.writeto(dirout+'MasterBias.fits')
    print "\t\t-> MasterBias: done!"

    Flat, RO_fl, GA_fl = uvesutils.MedianCombine(flats,bias=MasterBias)
    hdu = pyfits.PrimaryHDU(Flat)
    if (os.access(dirout+'Flat.fits',os.F_OK)):
        os.remove(dirout+'Flat.fits')
    hdu.writeto(dirout+'Flat.fits')
    print "\t\t-> MasterFlat: done!"

    hdu = pyfits.open(orderref[0])
    d1 = hdu[1].data - MasterBias[:,:,0]
    d2 = hdu[2].data - MasterBias[:,:,1]
    d1 = d1.T
    d2 = d2.T
    print "\tTracing echelle orders..."
    c_all1, nord1 = GLOBALutils.get_them(d1,int(0.2*ext_aperture),trace_degree,mode=1,nsigmas=10)
    print '\t\t'+str(nord1)+' orders for chip1 found...'
    c_all2, nord2 = GLOBALutils.get_them(d2,int(0.2*ext_aperture),trace_degree,mode=1,nsigmas=10, endat=2030)
    print '\t\t'+str(nord2)+' orders for chip2 found...'

    trace_dict = {'c_all1':c_all1,'nord1':nord1,\
                  'c_all2':c_all2,'nord2':nord2,\
                  'GA_bias': GA_bias, 'RO_bias': RO_bias,\
                  'GA_fl': GA_fl, 'RO_fl': RO_fl}
    pickle.dump( trace_dict, open( dirout+'trace.pkl', 'w' ) )

else:
    trace_dict = pickle.load( open( dirout+'trace.pkl' ) )
    c_all1 = trace_dict['c_all1']
    nord1  = trace_dict['nord1']
    c_all2 = trace_dict['c_all2']
    nord2  = trace_dict['nord2']

    GA_fl = trace_dict['GA_fl']
    RO_fl = trace_dict['RO_fl']
    GA_bias = trace_dict['GA_bias']
    RO_bias = trace_dict['RO_bias']

    Flat = pyfits.getdata(dirout+'Flat.fits')
    MasterBias = pyfits.getdata(dirout+'MasterBias.fits')

print '\n\tNormalization of Flat calibration frame:'
# Normalization of flat frame
if (os.access(dirout+'NFlat.fits',os.F_OK)==False) or (os.access(dirout+'Flat_spec.pkl',os.F_OK)==False) or force_flat_nor:

    ccdflat1 = uvesutils.ccd_flat(Flat[:,:,0],c_all1,sky_aperture)
    ccdflat2 = uvesutils.ccd_flat(Flat[:,:,1],c_all2,sky_aperture)
    ccdflat  = np.dstack((ccdflat1,ccdflat2))

    nflat1,sflat1 = uvesutils.normalize_flat(Flat[:,:,0],c_all1,sky_aperture)
    nflat2,sflat2 = uvesutils.normalize_flat(Flat[:,:,1],c_all2,sky_aperture)

    nflat = np.dstack((nflat1,nflat2))
    sflat = {'chip1':sflat1, 'chip2':sflat2}

    hdu = pyfits.PrimaryHDU(nflat)
    if (os.access(dirout+'NFlat.fits',os.F_OK)):
        os.remove(dirout+'NFlat.fits')
    hdu.writeto(dirout+'NFlat.fits')

    pickle.dump( sflat, open( dirout+'Flat_spec.pkl', 'w' ) )
	
else:
    nflat = pyfits.getdata(dirout+'NFlat.fits')
    sflat = pickle.load( open( dirout+'Flat_spec.pkl' ) )

print '\n\tExtraction of ThAr calibration frames:'
for fsim in ThAr_ref:
    thar_fits = dirout + fsim.split('/')[-1][:-4]+'spec.pkl'
    if ( os.access(thar_fits,os.F_OK) == False ) or ( force_thar_extract ):
        hthar = pyfits.open( fsim )
        dthar = (np.dstack((hthar[1].data,hthar[2].data)) -  MasterBias) / nflat
        print "\t\tNo previous extraction or extraction forced for ThAr file", fsim, "extracting..."
        dthar1 = dthar[:,:,0].T
        dthar2 = dthar[:,:,1].T 

        thar1_S = GLOBALutils.simple_extraction(dthar1,c_all1,ext_aperture,min_extract_cols1,max_extract_cols1,npools)
        thar2_S = GLOBALutils.simple_extraction(dthar2,c_all2,ext_aperture,min_extract_cols2,max_extract_cols2,npools)
        thar1_S = GLOBALutils.invert(thar1_S)
        thar2_S = GLOBALutils.invert(thar2_S)
        thar_S = {'chip1':thar1_S, 'chip2':thar2_S}
        pickle.dump( thar_S, open( thar_fits, 'w' ) )

    else:
        print "\t\tThAr file", fsim, "all ready extracted, loading..."

print "\n\tWavelength solution of ThAr calibration spectra:"
for fsim in ThAr_ref:

    hthar = pyfits.open( fsim )
    thar_fits  = dirout + fsim.split('/')[-1][:-4]+'spec.pkl'
    wavsol_pkl = dirout + fsim.split('/')[-1][:-4]+'spec.wavsolpars.pkl'

    if ( os.access(wavsol_pkl,os.F_OK) == False ) or ( force_thar_wavcal ):

        print "\t\tWorking on initial ThAr file", fsim
        mjd, mjd0 = uvesutils.mjd_fromheader2( hthar[0].header )
        thar_hd   = hthar[0].header
        thar_S    = pickle.load(open(thar_fits,'r'))

        lines_thar2       = thar_S['chip1']
        All_Pixel_Centers = np.array([])
        All_Wavelengths   = np.array([])
        All_Orders        = np.array([])
        All_Centroids     = np.array([])
        All_Sigmas        = np.array([])
        All_Intensities   = np.array([])

        for order in range(len(c_all1)):
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)
            thar_order_orig = lines_thar2[order]
            I = np.where(np.isnan(thar_order_orig))
            J = np.where(np.isnan(thar_order_orig)==False)
            thar_order_orig[I] = np.median(thar_order_orig[J])
            thar_order         = thar_order_orig - scipy.signal.medfilt(thar_order_orig,101)
            wei = np.ones(len(thar_order))
            if order == 0:
                do_xc = False
            else:
                do_xc = True
            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths,\
            rms_ms, residuals, centroids, sigmas, intensities =\
                     GLOBALutils.Initial_Wav_Calibration(order_dir+'nborder_'+\
                     order_s+'.iwdat', thar_order, order, wei, rmsmax=100, \
                     minlines=20,FixEnds=False,Dump_Argon=dumpargon,sigmai=2.,\
                     Dump_AllLines=True, Cheby=use_cheby,porder=ncoef_x,do_xc=False, line_width=6,pixelization=True)

            if order == 7: 
                Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 0.5*(len(thar_order)-1), len(thar_order))

            All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
            All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
            All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order )
            All_Centroids     = np.append( All_Centroids, centroids)
            All_Sigmas        = np.append( All_Sigmas, sigmas)
            All_Intensities   = np.append( All_Intensities, intensities )

        p0 = np.zeros( npar_wsol )
        p0[0] =  (7+ro0) * Global_ZP 
        p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
            GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders,\
                                                np.ones(All_Intensities.shape), p0, Cheby=use_cheby,\
                                                maxrms=MRMS, Inv=Inverse_m,minlines=200,order0=ro0, \
                                                ntotal=nord1,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
        """
        equis = np.arange(len(thar_order))
        order = 0
        while order < nord1:
            m   = order + ro0
            chebs = GLOBALutils.Calculate_chebs(equis, m, order0=ro0,ntotal=nord1, \
                npix=len(equis), Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)

            WavSol = (1.0/float(m)) * GLOBALutils.Joint_Polynomial_Cheby(p1,chebs,ncoef_x,ncoef_m)
            plot(WavSol,lines_thar2[order])
            order+=1
        show()
        """
        lines_thar2       = thar_S['chip2']
        All_Pixel_Centers2 = np.array([])
        All_Wavelengths2   = np.array([])
        All_Orders2        = np.array([])
        All_Centroids2     = np.array([])
        All_Sigmas2        = np.array([])
        All_Intensities2   = np.array([])

        for order in range(len(c_all2)):
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)

            thar_order_orig = lines_thar2[order]
            I = np.where(np.isnan(thar_order_orig))
            J = np.where(np.isnan(thar_order_orig)==False)
            thar_order_orig[I] = np.median(thar_order_orig[J])
            thar_order         = thar_order_orig - scipy.signal.medfilt(thar_order_orig,101)
            wei = np.ones(len(thar_order))

            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths,\
            rms_ms, residuals, centroids, sigmas, intensities =\
					 GLOBALutils.Initial_Wav_Calibration(order_dir+'nrorder_'+\
					 order_s+'.iwdat', thar_order, order, wei, rmsmax=100, \
					 minlines=20,FixEnds=False,Dump_Argon=dumpargon,sigmai=2.,\
					 Dump_AllLines=True, Cheby=use_cheby,porder=ncoef_x,do_xc=False, line_width=6,pixelization=True)

            if order == 11: 
                Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 0.5*(len(thar_order)-1), len(thar_order))

            All_Pixel_Centers2 = np.append( All_Pixel_Centers2, pixel_centers )
            All_Wavelengths2   = np.append( All_Wavelengths2, wavelengths )
            All_Orders2        = np.append( All_Orders2, np.zeros( len(pixel_centers) ) + order )
            All_Centroids2     = np.append( All_Centroids2, centroids)
            All_Sigmas2        = np.append( All_Sigmas2, sigmas)
            All_Intensities2   = np.append( All_Intensities2, intensities )
	
			#isz  = pixel_centers - sigmas
			#der  = pixel_centers + sigmas
			#isz  = GLOBALutils.Cheby_eval( coeffs_pix2wav, isz, len(thar_order))
			#der  = GLOBALutils.Cheby_eval( coeffs_pix2wav, der, len(thar_order))
			#sig  = 0.5*(der-isz)
			#fwhm = 2.35 * sig
			#resol = wavelengths / fwhm
			#print order, np.median(resol)
			#plot(pixel_centers,resol,'o')
		#show()
        #print dfgh

        p0 = np.zeros( npar_wsol )
        p0[0] =  (11 + ro0 + len(c_all1) + order_gap) * Global_ZP 
        p2, G_pix2, G_ord2, G_wav2, II2, rms_ms2, G_res2 = \
			GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers2, All_Wavelengths2, All_Orders2,\
                                                np.ones(All_Intensities2.shape), p0, Cheby=use_cheby,\
                                                maxrms=MRMS, Inv=Inverse_m,minlines=400,order0=ro0 + len(c_all1) + order_gap, \
                                                ntotal=nord2,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
        """
        equis = np.arange(len(thar_order))
        order = 0
        while order < nord2:
            m   = order + nord1 + order_gap + ro0
            chebs = GLOBALutils.Calculate_chebs(equis, m, order0=nord1 + order_gap + ro0,ntotal=nord2, \
                npix=len(equis), Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
            WavSol = (1.0/float(m)) * GLOBALutils.Joint_Polynomial_Cheby(p2,chebs,ncoef_x,ncoef_m)
            plot(WavSol,lines_thar2[order])
            order+=1
        show()
        """
        #for i in np.unique(G_ord2):
        #    I = np.where(G_ord2 == i)[0]
        #    plot(G_wav2[I],G_res2[I],'.')
        #show()
        pdict = {'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                 'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Orders':All_Orders, 'All_Sigmas':All_Sigmas,\
                 'p2':p2, 'G_pix2':G_pix2, 'G_ord2':G_ord2, 'G_wav2':G_wav2, 'II2':II2, 'rms_ms2':rms_ms2,\
                 'G_res2':G_res2, 'All_Centroids2':All_Centroids2, 'All_Orders2':All_Orders2, 'All_Sigmas2':All_Sigmas2}
        pickle.dump( pdict, open( wavsol_pkl, 'w' ) )

    else:
        print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl			

print '\n\tThe following targets will be processed:'
new_list = []
if object2do == 'all':
	new_list = sim_sci
else:
	for fsim in sim_sci:
		hd = pyfits.getheader(fsim)
		if object2do in hd['OBJECT']:
			new_list.append(fsim)

for fsim in new_list:
	hd = pyfits.getheader(fsim)
	print '\t\t'+str(hd['OBJECT'])

# Does any image have a special requirement for dealing with the moonlight?
if os.access(dirin + 'moon_corr.txt', os.F_OK):
    fmoon = open(dirin + 'moon_corr.txt','r')
    moon_lns = fmoon.readlines()
    spec_moon = []
    use_moon = []
    for line in moon_lns:
        spec_moon.append(line.split()[0])
        if line.split()[1] == '0':
            use_moon.append(False)
        else:
            use_moon.append(True)
else:
    spec_moon = []
    use_moon = []

spec_moon = np.array(spec_moon)
use_moon  = np.array(use_moon)

# Extract all Science files
for fsim in new_list:
    print '\n'
    print "\t--> Working on image: ", fsim
    know_moon = False

    if fsim.split('/')[-1] in spec_moon:
        I = np.where(fsim.split('/')[-1] == spec_moon)[0]
        know_moon = True
        here_moon = use_moon[I]

    h = pyfits.open(fsim)
    mjd,mjd0      = uvesutils.mjd_fromheader2(h[0].header)    
    ron1,gain1 = float(h[1].header['HIERARCH ESO DET OUT1 RON']),float(h[1].header['HIERARCH ESO DET OUT1 GAIN'])
    ron2,gain2 = float(h[2].header['HIERARCH ESO DET OUT1 RON']),float(h[2].header['HIERARCH ESO DET OUT1 GAIN'])

    # Object name
    obname    = h[0].header['OBJECT'].replace(' ','')

    print "\t\tObject name:",obname
    # Open file, trim, overscan subtract and MasterBias subtract
    bac_fits = dirout + 'BAC_' + fsim.split('/')[-1][:-4]+'.fits'
    sky_fits = dirout + 'SKY_' + fsim.split('/')[-1][:-4]+'.fits'
    #data = hiresutils.OverscanTrim( h[chip].data, h[chip].header['DATASEC'] ).T
    #c_new,pshift = GLOBALutils.retrace( data, c_all, span=15 )
    #print '\t\ty shift = ', pshift, 'pxs'

    data = np.dstack((h[1].data,h[2].data)) -  MasterBias
    centers1 = np.zeros((len(c_all1),data.shape[0]))
    centers2 = np.zeros((len(c_all2),data.shape[0]))

    for i in range(len(c_all1)):
        centers1[i] = np.polyval(c_all1[i],np.arange(data.shape[0]))
    for i in range(len(c_all2)):
        centers2[i] = np.polyval(c_all2[i],np.arange(data.shape[0]))

    print '\t\tEstimating background of scattered light...'
    force_bac=False
    if (os.access(bac_fits,os.F_OK) == False) or force_bac:
        bac1 = GLOBALutils.get_scat(data[:,:,0].T,centers1,span=sky_aperture)
        bac2 = GLOBALutils.get_scat(data[:,:,1].T,centers2,span=sky_aperture)
        bac = np.dstack((bac1.T,bac2.T))
        bach = pyfits.PrimaryHDU(bac)
        if os.access(bac_fits,os.F_OK):
            os.system('rm ' + bac_fits)
        bach.writeto(bac_fits)
    else:
        bac = pyfits.getdata(bac_fits)

    data = data - bac

    if (os.access(sky_fits,os.F_OK) == False) or force_bac:
        sky1 = uvesutils.get_sky(data[:,:,0].T,centers1,ext_aperture,sky_aperture)
        sky2 = uvesutils.get_sky(data[:,:,1].T,centers2,ext_aperture,sky_aperture)
        sky = np.dstack((sky1.T,sky2.T))

        bach = pyfits.PrimaryHDU(sky)
        if os.access(sky_fits,os.F_OK):
            os.system('rm ' + sky_fits)
        bach.writeto(sky_fits)
    else:
        sky = pyfits.getdata(sky_fits)

    data = data - sky 
    data = data / nflat

    ra  = float(h[0].header['RA'])
    dec = float(h[0].header['DEC'])

    ra2,dec2 = GLOBALutils.getcoords(obname,mjd,filen=reffile)
    if ra2 !=0 and dec2 != 0:
        ra = ra2
        dec = dec2
    else:
        print '\t\tUsing the coordinates found in the image header.'

    # set observatory parameters
    altitude    = float(h[0].header['ESO TEL GEOELEV'])
    latitude    = float(h[0].header['ESO TEL GEOLAT'])
    longitude   = float(h[0].header['ESO TEL GEOLON'])
    epoch       = 2000.

    iers                    = GLOBALutils.JPLiers( baryc_dir, mjd-999.0, mjd+999.0 )
    obsradius, R0           = GLOBALutils.JPLR0( latitude, altitude)
    obpos                   = GLOBALutils.obspos( longitude, obsradius, R0 )

    jplephem.set_ephemeris_dir( baryc_dir , ephemeris )
    jplephem.set_observer_coordinates( float(obpos[0]), float(obpos[1]), float(obpos[2]) )

    res         = jplephem.doppler_fraction(float(ra/15.0), float(dec), long(mjd), float(mjd%1), 1, 0.0)
    lbary_ltopo = 1.0 + res['frac'][0]
    bcvel_baryc = ( lbary_ltopo - 1.0 ) * 2.99792458E5

    print "\t\tBarycentric velocity:", bcvel_baryc

    res  = jplephem.pulse_delay(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
    mbjd = mjd + res['delay'][0] / (3600.0 * 24.0)

    # Moon Phase Calculations
    gobs      = ephem.Observer()  
    gobs.name = 'VLT'
    gobs.lat  = rad(latitude)  # lat/long in decimal degrees  
    gobs.long = rad(longitude) 
    gobs.date = h[0].header['DATE-OBS'].replace('T',' ')
    mephem    = ephem.Moon()
    mephem.compute(gobs)
    Mcoo        = jplephem.object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Mp   = jplephem.barycentric_object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Sp   = jplephem.barycentric_object_track("Sun", int(mjd), float(mjd%1), 1, 0.0)
    res  = jplephem.object_doppler("Moon", int(mjd), mjd%1, 1, 0.0)
    lunation,moon_state,moonsep,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,ra,dec)
    refvel = bcvel_baryc + moonvel
    print '\t\tRadial Velocity of sacttered moonlight:',refvel

    sci_fits        = dirout + fsim.split('/')[-1][:-4]+'spec.pkl'
    sci_fits_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.pkl'
    P_fits          = dirout + 'P_' + fsim.split('/')[-1][:-4]+'spec.fits'

    if ( os.access(sci_fits,os.F_OK) == False ) or ( os.access(sci_fits_simple,os.F_OK) == False ) or \
       ( force_sci_extract ):
        print "\t\tNo previous extraction or extraction forced for science file", fsim, "extracting..."
        #"""
        print "\t\t\tweights for chip1..."
        P1 = GLOBALutils.obtain_P(data[:,:,0].T,c_all1,ext_aperture,ron1,gain1,NSigma_Marsh, S_Marsh, \
                                    N_Marsh, Marsh_alg, min_extract_cols1,max_extract_cols1, npools)
        print "\t\t\tweights for chip2..."
        P2 = GLOBALutils.obtain_P(data[:,:,1].T,c_all2,ext_aperture,ron2,gain2,NSigma_Marsh, S_Marsh, \
                                    N_Marsh, Marsh_alg, min_extract_cols2,max_extract_cols2, npools)
        P = np.dstack((P1,P2))
        #"""
        #P = pyfits.getdata(P_fits)
        print "\t\t\tsimple extraction for chip1..."
        sci_Ss1 = GLOBALutils.simple_extraction(data[:,:,0].T,c_all1,ext_aperture,min_extract_cols1,max_extract_cols1,npools)
        print "\t\t\tsimple extraction for chip2..."
        sci_Ss2 = GLOBALutils.simple_extraction(data[:,:,1].T,c_all2,ext_aperture,min_extract_cols2,max_extract_cols2,npools)
        
        sci_Ss1 = GLOBALutils.invert(sci_Ss1)
        sci_Ss2 = GLOBALutils.invert(sci_Ss2)
        #"""
        print "\t\t\t\optimal extraction for chip1..."
        sci_S1  = GLOBALutils.optimal_extraction(data[:,:,0].T,P1,c_all1,ext_aperture,ron1,gain1,S_Marsh,NCosmic_Marsh,\
		 min_extract_cols1,max_extract_cols1,npools)
        print "\t\t\t\optimal extraction for chip2..."
        sci_S2  = GLOBALutils.optimal_extraction(data[:,:,1].T,P2,c_all2,ext_aperture,ron2,gain2,S_Marsh,NCosmic_Marsh,\
         min_extract_cols2,max_extract_cols2,npools)
        #"""
        #sci_S = pickle.load(open(sci_fits,'r'))
        #sci_S1 = sci_S['chip1']
        #sci_S2 = sci_S['chip2']

        sci_S1 = GLOBALutils.invert(sci_S1)
        sci_S2 = GLOBALutils.invert(sci_S2)
        
        if (os.access(P_fits,os.F_OK)):
            os.remove( P_fits )
        hdu = pyfits.PrimaryHDU( P )
        hdu.writeto( P_fits )
        sci_S  = {'chip1':sci_S1, 'chip2':sci_S2}
        pickle.dump(sci_S,open(sci_fits,'w'))
        
        sci_Ss = {'chip1':sci_Ss1, 'chip2':sci_Ss2}
        pickle.dump(sci_Ss,open(sci_fits_simple,'w'))
        
    else:
        print '\t\t'+fsim+" has already been extracted, reading in product fits files..."
        sci_S  = pickle.load(open(sci_fits,'r'))
        sci_Ss = pickle.load(open(sci_fits_simple,'r'))

    fout = 'proc/'+ obname + '_' + h[0].header['DATE-OBS'] + '_' + 'sp.fits'
    if ( os.access(dirout+fout ,os.F_OK) == False ) or (force_spectral_file_build):

        spec = np.zeros((11, nord1 + nord2, data.shape[0]))
        hdu = pyfits.PrimaryHDU( spec )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MJD', mjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MBJD', mbjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START DATE', h[0].header['DATE-OBS'].split('T')[0] )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START UT',  h[0].header['DATE-OBS'].split('T')[1])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TEXP (S)',h[0].header['EXPTIME'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (KM/S)', bcvel_baryc)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH (LAMBDA_BARY / LAMBDA_TOPO)', lbary_ltopo)    
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARGET NAME', obname)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA',h[0].header['RA'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC',h[0].header['DEC'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA BARY',ra)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC BARY',dec)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH EQUINOX',2000.)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LATITUDE',latitude)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LONGITUDE',longitude)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS ALTITUDE',altitude)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARG AIRMASS',h[0].header['ESO OBS AIRM'])

        models_path = base+"data/COELHO_MODELS/R_40000b/"
        ref_RES = 60000.
        # get ThAr closest in time
        im = np.argmin(np.absolute(np.array(ThAr_ref_dates) - mjd))
        pkl_wsol = dirout + ThAr_ref[im].split('/')[-1][:-4]+'spec.wavsolpars.pkl'
        print "\t\tUnpickling wavelength solution from", pkl_wsol, " ..."
        wsol_dict = pickle.load(open(pkl_wsol,'r'))

        # Apply new wavelength solution including barycentric correction
        equis = np.arange( data.shape[0] )        
        for order in range(nord1):
            #order = 16
            #print order
            m = order + ro0
            chebs = GLOBALutils.Calculate_chebs(equis, m, npix=data.shape[0],\
		    order0=ro0, ntotal=nord1, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
            WavSol = lbary_ltopo * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wsol_dict['p1'],chebs,ncoef_x,ncoef_m)   
            spec[0,order,:] = GLOBALutils.ToVacuum(WavSol)
            spec[1,order,:] = sci_S['chip1'][order,1,:]
            #spec[1,order,:300] = 0.
            spec[2,order,:] = sci_S['chip1'][order,2,:]
            spec[3,order,:] = spec[1,order,:] / sflat['chip1'][order][::-1]
            spec[4,order,:] = spec[2,order,:] * sflat['chip1'][order][::-1] ** 2
            ccoef = GLOBALutils.get_cont_single(spec[0,order],spec[3,order],spec[4,order],ll=1.5,lu=5,nc=3)
            ratio = np.polyval(ccoef, spec[0,order])
            L     = np.where( spec[1,order,:] != 0 )
            spec[5,order,:][L] = spec[3,order,:][L] / ratio[L]
            spec[3,order,:][L] = spec[3,order,:][L] / ratio[L]
            spec[6,order,:][L] = spec[4,order,:][L] * (ratio[L] ** 2 )
            spec[7,order,:][L] = ratio[L]
            spec[8,order,:][L] = ratio[L] * sflat['chip1'][order][::-1][L] / np.sqrt( ratio[L] * sflat['chip1'][order][::-1][L] / gain1 + (ron1/gain1)**2 )
            rI = np.where(spec[5,order] > 1. + 5./spec[8,order])
            spec[5,order,rI] = 1.
            lI = np.where(spec[5,order] < -5./spec[8,order])
            spec[5,order,lI] = 0.

            spl           = scipy.interpolate.splrep(np.arange(len(WavSol)), WavSol,k=3)
            dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
            NN            = np.average(dlambda_dx)
            dlambda_dx    /= NN
            # clean-up of CRs in continuum-normalized spectrum. Set troublesome pixels to 1
            medflx = np.zeros(len(spec[5,order]))
            nI = np.where(np.isnan(spec[5,order])==False)[0]
            nJ = np.where(np.isnan(spec[5,order])==True)[0]
            spec[6,order,:][L] = spec[8,order,:][L] ** 2
            spec[9,order,:][L] = spec[5,order,:][L] * (dlambda_dx[L] ** 1) 
            spec[10,order,:][L] = spec[6,order,:][L] / (dlambda_dx[L] ** 2)


        for order in range(nord2):
            order2 = order + nord1
            m = order + ro0 + nord1 + order_gap
            chebs = GLOBALutils.Calculate_chebs(equis, m, npix=data.shape[0],\
            order0=ro0 + nord1 + order_gap, ntotal=nord2, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
            WavSol = lbary_ltopo * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wsol_dict['p2'],chebs,ncoef_x,ncoef_m)   
            spec[0,order2,:] = GLOBALutils.ToVacuum(WavSol)
            spec[1,order2,:] = sci_S['chip2'][order,1,:]
            #spec[1,order2,:300] = 0.
            spec[2,order2,:] = sci_S['chip2'][order,2,:]
            spec[3,order2,:] = spec[1,order2,:] / sflat['chip2'][order][::-1]
            spec[4,order2,:] = spec[2,order2,:] * sflat['chip2'][order][::-1] ** 2
            ccoef = GLOBALutils.get_cont_single(spec[0,order2],spec[3,order2],spec[4,order2],ll=1.5,lu=5,nc=3)
            ratio = np.polyval(ccoef, spec[0,order2])
            L     = np.where( spec[1,order2,:] != 0 )
            spec[5,order2,:][L] = spec[3,order2,:][L] / ratio[L]
            spec[3,order2,:][L] = spec[3,order2,:][L] / ratio[L]
            spec[6,order2,:][L] = spec[4,order2,:][L] * (ratio[L] ** 2 )
            spec[7,order2,:][L] = ratio[L]
            spec[8,order2,:][L] = ratio[L] * sflat['chip2'][order][::-1][L] / np.sqrt( ratio[L] * sflat['chip2'][order][::-1][L] / gain2 + (ron2/gain2)**2 )
            rI = np.where(spec[5,order2] > 1. + 5./spec[8,order2])
            #lI = np.where(spec[5,order2] < -8./spec[8,order2])
            lI = np.where(spec[5,order2] < -5./spec[8,order2])
            spec[5,order2,rI] = 1.
            spec[5,order2,lI] = 0.

            spl           = scipy.interpolate.splrep(np.arange(len(WavSol)), WavSol,k=3)
            dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
            NN            = np.average(dlambda_dx)
            dlambda_dx    /= NN
            # clean-up of CRs in continuum-normalized spectrum. Set troublesome pixels to 1
            medflx = np.zeros(len(spec[5,order2]))
            nI = np.where(np.isnan(spec[5,order2])==False)[0]
            nJ = np.where(np.isnan(spec[5,order2])==True)[0]
            spec[6,order2,:][L] = spec[8,order2,:][L] ** 2
            spec[9,order2,:][L] = spec[5,order2,:][L] * (dlambda_dx[L] ** 1) 
            spec[10,order2,:][L] = spec[6,order2,:][L] / (dlambda_dx[L] ** 2)

    if (not JustExtract):
        if DoClass:
            # spectral analysis
            print "\t\tSpectral Analysis..."
            query_success = False
            # First, query SIMBAD with the object name
            query_success,sp_type_query = GLOBALutils.simbad_query_obname(obname)
            # Now, query SIMBAD by coordinates if above not successful
            if (not query_success):
                query_success,sp_type_query = GLOBALutils.simbad_query_coords('12:00:00','00:00:00')
            print "\t\t\tSpectral type returned by SIMBAD query:",sp_type_query

            hdu = GLOBALutils.update_header(hdu,'HIERARCH SIMBAD SPTYP', sp_type_query)

            pars_file = dirout + fsim.split('/')[-1][:-4]+'_stellar_pars.txt'

            if os.access(pars_file,os.F_OK) == False or force_stellar_pars:
                print "\t\t\tEstimating atmospheric parameters:"
                Rx = np.around(1./np.sqrt(1./40000.**2 - 1./ref_RES**2))
                spec2 = spec.copy()
                for i in range(spec.shape[1]):
                    IJ = np.where(spec[5,i]!=0.)[0]
                    spec2[5,i,IJ] = GLOBALutils.convolve(spec[0,i,IJ],spec[5,i,IJ],Rx)
                T_eff, logg, Z, vsini, vel0, ccf = correlation.CCF(spec2,model_path=models_path,npools=npools)
                line = "%6d %4.1f %4.1f %8.1f %8.1f\n" % (T_eff,logg, Z, vsini, vel0)
                f = open(pars_file,'w')
                f.write(line)
                f.close()
		       
            else:
                print "\t\t\tAtmospheric parameters loaded from file:"
                T_eff, logg, Z, vsini, vel0 = np.loadtxt(pars_file,unpack=True)

            print "\t\t\t\tT_eff=",T_eff,"log(g)=",logg,"Z=",Z,"vsin(i)=",vsini,"vel0",vel0

        else:
            T_eff, logg, Z, vsini, vel0 = -999,-999,-999,-999,-999

        T_eff_epoch = T_eff
        logg_epoch  = logg
        Z_epoch     = Z
        vsini_epoch = vsini
        vel0_epoch  = vel0
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TEFF', float(T_eff))
        hdu = GLOBALutils.update_header(hdu,'HIERARCH LOGG', float(logg))
        hdu = GLOBALutils.update_header(hdu,'HIERARCH Z', Z)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH VSINI', vsini)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH VEL0', vel0)

        print "\t\tRadial Velocity analysis:"
        # assign mask
        sp_type, mask = GLOBALutils.get_mask_reffile(obname,reffile=reffile,base='../data/xc_masks/')
        print "\t\t\tWill use",sp_type,"mask for CCF."
		
        velw  = 300
        velsh = 3.
            
        # Read in mask
        ml, mh, weight = np.loadtxt(mask,unpack=True)
        ml_v = GLOBALutils.ToVacuum( ml )
        mh_v = GLOBALutils.ToVacuum( mh )
       
        # make mask larger accounting for factor ~2.5 lower res in FEROS w/r to HARPS
        av_m = 0.5*( ml_v + mh_v )
        ml_v -= 1.5*(av_m - ml_v)
        mh_v += 1.5*(mh_v - av_m)
        mask_hw_kms = (GLOBALutils.Constants.c/1e3) * 0.5*(mh_v - ml_v) / av_m

        disp = GLOBALutils.get_disp(obname, reffile=reffile)
        if disp == 0:
            known_sigma = False
            if vsini != -999 and vsini != 0.:
                disp = vsini
            else:
                disp = 3.
        else:
            known_sigma = True

        mask_hw_wide = av_m * disp / (GLOBALutils.Constants.c/1.0e3)
        ml_v = av_m - mask_hw_wide
        mh_v = av_m + mask_hw_wide 

        print '\t\t\tComputing the CCF...'
        cond = True
        while (cond):
            # first rough correlation to find the minimum
            vels, xc_full, sn, nlines_ccf, W_ccf = \
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight, 0, lbary_ltopo, vel_width=velw,vel_step=velsh,\
                                          spec_order=9,iv_order=10,sn_order=8,max_vel_rough=velw)

            xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3., Simple=True, W=W_ccf)
            yy     = scipy.signal.medfilt(xc_av,11)
            pred = lowess(yy, vels,frac=0.4,it=10,return_sorted=False)
            tck1 = scipy.interpolate.splrep(vels,pred,k=1)
            xc_av_orig = xc_av.copy()
            xc_av /= pred
            vel0_xc = vels[ np.argmin( xc_av ) ]
            rvels, rxc_av, rpred, rxc_av_orig, rvel0_xc = vels.copy(), xc_av.copy(), pred.copy(), xc_av_orig.copy(), vel0_xc
            xc_av_rough = xc_av
            vels_rough  = vels
                
            vel_width = np.maximum( 20.0, 6*disp )

            vels, xc_full, sn, nlines_ccf, W_ccf =\
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight, vel0_xc, lbary_ltopo, vel_width=vel_width,vel_step=0.2,\
                                          spec_order=9,iv_order=10,sn_order=8,max_vel_rough=velw)
		
            xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3., Simple=True, W=W_ccf)
            pred = scipy.interpolate.splev(vels,tck1)
            xc_av /= pred
		
            if sp_type == 'M5':
                moon_sig = 2.5
            elif sp_type == 'K5':
                moon_sig = 3.3
            else:
                moon_sig = 4.5

            p1,XCmodel,p1gau,XCmodelgau,Ls2 = GLOBALutils.XC_Final_Fit( vels, xc_av , sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = False)

            moonmatters = False
            if (know_moon and here_moon):
                moonmatters = True
                ismoon = True
                confused = False
                p1_m,XCmodel_m,p1gau_m,XCmodelgau_m,Ls2_m = GLOBALutils.XC_Final_Fit( vels, xc_av , sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = True)
                moon_flag = 1
            else:
                confused = False
                ismoon = False
                p1_m,XCmodel_m,p1gau_m,XCmodelgau_m,Ls2_m = p1,XCmodel,p1gau,XCmodelgau,Ls2
                moon_flag = 0

            SP = GLOBALutils.calc_bss2(vels,xc_av,p1gau)
            #SP = bspan[0]

            #print p1gau[1]
            if (not known_sigma):
                disp = np.floor(p1gau[2])
                if (disp < 3.0): 
                    disp = 3.0
                mask_hw_wide = av_m * disp / (GLOBALutils.Constants.c/1.0e3)
                ml_v = av_m - mask_hw_wide
                mh_v = av_m + mask_hw_wide            
                known_sigma = True
            else:
                cond = False

            BSerr = -999.00

        xc_dict = {'vels':vels,'xc_av':xc_av,'XCmodelgau':XCmodelgau,'Ls2':Ls2,'refvel':refvel,\
		       'rvels':rvels,'rxc_av':rxc_av,'rpred':rpred,'rxc_av_orig':rxc_av_orig,\
		       'rvel0_xc':rvel0_xc,'xc_full':xc_full, 'p1':p1, 'sn':sn, 'p1gau':p1gau,\
		       'p1_m':p1_m,'XCmodel_m':XCmodel_m,'p1gau_m':p1gau_m,'Ls2_m':Ls2_m,\
		       'XCmodelgau_m':XCmodelgau_m}

        moon_dict = {'moonmatters':moonmatters,'moon_state':moon_state,'moonsep':moonsep,\
		         'lunation':lunation,'mephem':mephem,'texp':h[0].header['EXPTIME']}


        pkl_xc = dirout + fsim.split('/')[-1][:-4]+obname+'_XC_'+sp_type+'.pkl'
        pickle.dump( xc_dict, open( pkl_xc, 'w' ) )

        ccf_pdf = dirout + 'proc/' + fsim.split('/')[-1][:-4] + obname + '_XCs_' + sp_type + '.pdf'

        if not avoid_plot:
            GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)


        airmass  = h[0].header['ESO OBS AIRM']
        seeing   = -999
        TEXP = h[0].header['EXPTIME']    
        SNR_5130 = np.median(spec[8,13,1900:2101] )

        if sp_type=='G2':
            if T_eff < 6000:
                D = 0.32815
                C = 0.00453
            else:
                D = 0.32815
                C = 0.00453
        elif  sp_type == 'K5':
            D = 0.27404
            C = 0.00433
        else:
            D = 0.27404
            C = 0.00433

        BSerr = D / float(np.round(SNR_5130)) + C
        RVerr = 0.5

        RV     = np.around(p1gau_m[1],4)  
        BS     = np.around(SP,4)   
        RVerr2 = np.around(RVerr,4)
        BSerr  = np.around(BSerr,4)
        print '\t\t\tRV = '+str(RV)+' +- '+str(RVerr2)
        print '\t\t\tBS = '+str(BS)+' +- '+str(BSerr)

        bjd_out = 2400000.5 + mbjd
        T_eff_err = 100
        logg_err = 0.5
        Z_err = 0.5
        vsini_err = 2
        XC_min = np.abs(np.around(np.min(XCmodel),2))

        SNR_5130 = np.around(SNR_5130)
        SNR_5130_R = np.around(SNR_5130*np.sqrt(2.3))
        # write to output
        disp_epoch = np.around(p1gau_m[2],1)
        hdu = GLOBALutils.update_header(hdu,'RV', RV)
        hdu = GLOBALutils.update_header(hdu,'RV_E', RVerr2)
        hdu = GLOBALutils.update_header(hdu,'BS', BS)
        hdu = GLOBALutils.update_header(hdu,'BS_E', BSerr)
        hdu = GLOBALutils.update_header(hdu,'DISP', disp_epoch)
        hdu = GLOBALutils.update_header(hdu,'SNR', SNR_5130)
        hdu = GLOBALutils.update_header(hdu,'SNR_R', SNR_5130_R)
        hdu = GLOBALutils.update_header(hdu,'INST', 'UVES')
        hdu = GLOBALutils.update_header(hdu,'RESOL', ref_RES)
        hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')
        hdu = GLOBALutils.update_header(hdu,'XC_MIN', XC_min)
        hdu = GLOBALutils.update_header(hdu,'BJD_OUT', bjd_out)
        line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f   uves   ceres   %8d %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
                      (obname, bjd_out, RV, RVerr2, BS, BSerr, ref_RES, T_eff_epoch, logg_epoch, Z_epoch, vsini_epoch, XC_min, disp_epoch,\
		       TEXP, SNR_5130_R, ccf_pdf)
        f_res.write(line_out)
        if (os.access( dirout + fout,os.F_OK)):
            os.remove( dirout + fout)
        hdu.writeto( dirout + fout )

f_res.close()
