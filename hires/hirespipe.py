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
import hiresutils
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

# interface to R
from rpy2 import robjects
import rpy2.robjects.numpy2ri
try:
	rpy2.robjects.numpy2ri.activate()
except:
	None
#import rpy2.robjects.numpy2ri
r = robjects.r
r.library("MASS")

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
minlines_glob      = 700

Inverse_m          = True
use_cheby          = True
MRMS               = 50   # max rms in m/s, global wav solution

trace_degree       = 4
Marsh_alg          = 0
ext_aperture       = 30
NSigma_Marsh       = 5
NCosmic_Marsh      = 5
S_Marsh            = 0.4
N_Marsh            = 3      # grado polinomio 
min_extract_col    = 0
max_extract_col    = 4000
xbinning           = 1
ybinning           = 3
chip               = 2

if chip == 2:
	ro0 = 120

ext_aperture = int(float(ext_aperture)/float(ybinning)) 

ncoef_x   = 4
ncoef_m   = 5
npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2

order_dir  = base+"/hires/wavcals/"

#############################

print "\n\n\tHIRES KECK10m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'

# file containing the log
log = dirout+'night.log'
flats, ThAr_ref, sim_sci, ThAr_ref_dates, obnames, exptimes = hiresutils.FileClassify(dirin,log)

print '\tThis in the log of the night:\n'
f = open(log)
flines = f.readlines()
for line in flines:
	print '\t\t'+line[:-1]
print '\n'

if ((os.access(dirout+'Flat_' + str(int(chip)) + '.fits',os.F_OK) == False) or \
    (os.access(dirout+'trace_' + str(int(chip)) + '.pkl',os.F_OK) == False)  or \
    (force_pre_process) ):
    print "\tNo previous pre-processing files or found"
    pre_process = 1
else:
    print "\tPre-processing files found, going straight to extraction"
    pre_process = 0

if (pre_process == 1):
    # median combine list of flats
    print "\t\tGenerating Master calibration frames..."
    Flat, RO_fl, GA_fl = hiresutils.MedianCombine(flats,chip=chip)
    hdu = pyfits.PrimaryHDU(Flat)
    if (os.access(dirout+'Flat_' + str(int(chip)) + '.fits',os.F_OK)):
        os.remove(dirout+'Flat_' + str(int(chip)) + '.fits')
    hdu.writeto(dirout+'Flat_' + str(int(chip)) + '.fits')
    print "\t\t-> Masterflat: done!"

    hdu = pyfits.open(dirin+stst)
    d = hdu[chip].data
    d = hiresutils.OverscanTrim(d,hdu[chip].header['DATASEC'])
    d = d.T
    print "\tTracing echelle orders..."
    c_all, nord = GLOBALutils.get_them(d,ext_aperture,trace_degree,mode=1)
    print '\t\t'+str(nord)+' orders found...'

    # pickle traces
    trace_dict = {'c_all':c_all,
                  'nord':nord,
                  'GA_fl': GA_fl, 'RO_fl': RO_fl}
    pickle.dump( trace_dict, open( dirout+'trace_' + str(int(chip)) + '.pkl', 'w' ) )

else:
    trace_dict = pickle.load( open( dirout+'trace_' + str(int(chip)) + '.pkl' ) )
    c_all = trace_dict['c_all']
    nord  = trace_dict['nord']
    # recover GA*, RO*
    GA_fl = trace_dict['GA_fl']
    RO_fl = trace_dict['RO_fl']
    # recover flats & master bias
    h = pyfits.open(dirout+'Flat_' + str(int(chip)) + '.fits')
    Flat = h[0].data

print '\n\tNormalization of Flat calibration frame:'
# Normalization of flat frame
if (os.access(dirout+'NFlat_' + str(int(chip)) + '.fits',os.F_OK)==False) or (os.access(dirout+'Flat_spec_' + str(int(chip)) + '.fits',os.F_OK)==False) or force_flat_nor:
	if (os.access(dirout+'BAC_Flat_' + str(int(chip)) + '.fits',os.F_OK)==False) or force_bac:
		bac = hiresutils.bac_flat(Flat,c_all,ext_aperture)
		if os.access(dirout+'BAC_Flat_' + str(int(chip)) + '.fits',os.F_OK):
			os.system('rm '+dirout+'BAC_Flat_' + str(int(chip)) + '.fits')
		hdu = pyfits.PrimaryHDU(bac)
		hdu.writeto(dirout+'BAC_Flat_' + str(int(chip)) + '.fits')
	else:
		bac = pyfits.getdata(dirout+'BAC_Flat_' + str(int(chip)) + '.fits')

	ccdflat = hiresutils.ccd_flat(Flat,c_all,ext_aperture)
	nflat,sflat = hiresutils.normalize_flat(Flat,c_all,ext_aperture)
	hdu = pyfits.PrimaryHDU(nflat)
	if (os.access(dirout+'NFlat_' + str(int(chip)) + '.fits',os.F_OK)):
        	os.remove(dirout+'NFlat_' + str(int(chip)) + '.fits')
	hdu.writeto(dirout+'NFlat_' + str(int(chip)) + '.fits')
	hdu = pyfits.PrimaryHDU(sflat)
	if (os.access(dirout+'Flat_spec_' + str(int(chip)) + '.fits',os.F_OK)):
        	os.remove(dirout+'Flat_spec_' + str(int(chip)) + '.fits')
	hdu.writeto(dirout+'Flat_spec_' + str(int(chip)) + '.fits')
	
else:
	nflat = pyfits.getdata(dirout+'NFlat_' + str(int(chip)) + '.fits')
	sflat = pyfits.getdata(dirout+'Flat_spec_' + str(int(chip)) + '.fits')

print '\n\tExtraction of ThAr calibration frames:'
# Extract all ThAr files
for fsim in ThAr_ref:
    thar_fits        = dirout + fsim.split('/')[-1][:-4]+'spec_'+str(int(chip))+'.fits'
    if ( os.access(thar_fits,os.F_OK) == False ) or ( force_thar_extract ):
	hthar = pyfits.open( fsim )
	dthar = hiresutils.OverscanTrim( hthar[chip].data, hthar[chip].header['DATASEC'] ) / nflat
        print "\t\tNo previous extraction or extraction forced for ThAr file", fsim, "extracting..."
	dthar = dthar.T 
	thar_S = np.zeros( (nord,dthar.shape[1]) )
	tR,tG = hthar[0].header['CCDRN0'+str(int(chip))],hthar[0].header['CCDGN0'+str(int(chip))]

	thar_S = GLOBALutils.simple_extraction(dthar,c_all,ext_aperture,min_extract_col,max_extract_col,npools)

        # save as fits file
        if (os.access(thar_fits,os.F_OK)):
            os.remove( thar_fits )            
        hdu = pyfits.PrimaryHDU( thar_S )
        hdu.writeto( thar_fits )
    else:
        print "\t\tThAr file", fsim, "all ready extracted, loading..."

print "\n\tWavelength solution of ThAr calibration spectra:"
for fsim in ThAr_ref:
	hthar = pyfits.open( fsim )
	thar_fits = dirout + fsim.split('/')[-1][:-4]+'spec_'+str(int(chip))+'.fits'
	wavsol_pkl = dirout + fsim.split('/')[-1][:-4]+'spec_'+str(int(chip))+ '.wavsolpars.pkl'
	if ( os.access(wavsol_pkl,os.F_OK) == False ) or ( force_thar_wavcal ):
		print "\t\tWorking on initial ThAr file", fsim
		mjd, mjd0 = hiresutils.mjd_fromheader2( hthar[0].header )
		thar_hd  = pyfits.getheader(fsim)
		print "\t\tDecher Name:",thar_hd['DECKNAME']
		thar_S = pyfits.getdata( thar_fits )
		lines_thar  = thar_S.copy()
		All_Pixel_Centers = np.array([])
		All_Wavelengths   = np.array([])
		All_Orders        = np.array([])
		All_Centroids     = np.array([])
		All_Sigmas        = np.array([])
		All_Intensities   = np.array([])
		
		for order in range(len(c_all)):
			order_s = str(order)
			if (order < 10):
				order_s = '0'+str(order)
			thar_order_orig = lines_thar[order]
			I = np.where(np.isnan(thar_order_orig))
			J = np.where(np.isnan(thar_order_orig)==False)
			thar_order_orig[I] = np.median(thar_order_orig[J])
			thar_order      = thar_order_orig - scipy.signal.medfilt(thar_order_orig,101)
			wei = np.ones(len(thar_order))

			coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths,\
			rms_ms, residuals, centroids, sigmas, intensities =\
					 GLOBALutils.Initial_Wav_Calibration(order_dir+'echorder_'+\
					 order_s+'_'+str(int(chip))+'.iwdat', thar_order, order, wei, rmsmax=100, \
					 minlines=20,FixEnds=False,Dump_Argon=dumpargon,\
					 Dump_AllLines=True, Cheby=use_cheby,porder=ncoef_x)
	
			if (chip == 2)  and (order == 8): 
				Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 0.5*(len(thar_order)-1), len(thar_order))

			All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
			All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
			All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order )
			All_Centroids     = np.append( All_Centroids, centroids)
			All_Sigmas        = np.append( All_Sigmas, sigmas)
			All_Intensities   = np.append( All_Intensities, intensities )

	
			isz  = pixel_centers - sigmas
			der  = pixel_centers + sigmas
			isz  = GLOBALutils.Cheby_eval( coeffs_pix2wav, isz, len(thar_order))
			der  = GLOBALutils.Cheby_eval( coeffs_pix2wav, der, len(thar_order))
			sig  = 0.5*(der-isz)
			fwhm = 2.35 * sig
			resol = wavelengths / fwhm
			#print order, np.median(resol)
			#plot(pixel_centers,resol,'o')
		#show()
		
		p0 = np.zeros( npar_wsol )
		p0[0] =  (8+ro0) * Global_ZP 
		p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
			GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders,\
                                                np.ones(All_Intensities.shape), p0, Cheby=use_cheby,\
                                                maxrms=MRMS, Inv=Inverse_m,minlines=500,order0=ro0, \
                                                ntotal=nord,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

		pdict = {'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                     'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Orders':All_Orders, 'All_Sigmas':All_Sigmas}
		pickle.dump( pdict, open( wavsol_pkl, 'w' ) )

	else:
		print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl			

print '\n\tThe following targets will be processed:'
new_list = []
if object2do == 'all':
	new_list = sim_sci
elif object2do == 'no_I2':
	for fsim in sim_sci:
		h = pyfits.open(fsim)
		if h[0].header['IODIN'] == False:
			new_list.append(fsim)
else:
	for fsim in sim_sci:
		hd = pyfits.getheader(fsim)
		if object2do in hd['TARGNAME']:
			new_list.append(fsim)

for fsim in new_list:
	hd = pyfits.getheader(fsim)
	print '\t\t'+str(hd['TARGNAME'])


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
    mjd,mjd0      = hiresutils.mjd_fromheader2(h[0].header)
    ronoise, gain = float(h[0].header['CCDRN0'+str(int(chip))]),float(h[0].header['CCDGN0'+str(int(chip))])
    
    # Object name
    obname    = h[0].header['TARGNAME'].replace(' ','')

    print "\t\tObject name:",obname
    print "\t\tDeckr Name = ", h[0].header['DECKNAME']
    # Open file, trim, overscan subtract and MasterBias subtract
    bac_fits = dirout + 'BAC_' + fsim.split('/')[-1][:-4]+'spec_'+str(int(chip))+'.fits'
    data = hiresutils.OverscanTrim( h[chip].data, h[chip].header['DATASEC'] ).T
    c_new,pshift = GLOBALutils.retrace( data, c_all, span=15 )
    print '\t\ty shift = ', pshift, 'pxs'

    centers = np.zeros((len(c_new),nflat.shape[0]))
    for i in range(len(c_new)):
    	centers[i] = np.polyval(c_new[i],np.arange(nflat.shape[0]))
    
    force_bac=False
    if (os.access(bac_fits,os.F_OK) == False) or force_bac:
    	bac = GLOBALutils.get_scat(data,centers,span=10)
	bach = pyfits.PrimaryHDU(bac)
	if os.access(bac_fits,os.F_OK):
		os.system('rm ' + bac_fits)
	bach.writeto(bac_fits)
    else:
	bac = pyfits.getdata(bac_fits)

    data = data - bac
    data = data / nflat.T

    ra = h[0].header['RA'].split(':')
    ra = float(ra[0]) + float(ra[1])/60. + float(ra[2])/3600.
    ra = 360.*ra/24.
    dec = h[0].header['DEC'].split(':')
    dec = float(dec[0]) + float(dec[1])/60. + float(dec[2])/3600.

    ra2,dec2 = GLOBALutils.getcoords(obname,mjd,filen=reffile)
    if ra2 !=0 and dec2 != 0:
	ra = ra2
	dec = dec2
    else:
	print '\t\tUsing the coordinates found in the image header.'

    # set observatory parameters
    altitude    = 4145.
    latitude    = 19.82636
    longitude   = -155.47501
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
    gobs.name = 'Keck'
    gobs.lat  = rad(latitude)  # lat/long in decimal degrees  
    gobs.long = rad(longitude) 
    gobs.date = h[0].header['DATE-OBS'] + ' ' + h[0].header['UTC'].replace(':','-')
    mephem    = ephem.Moon()
    mephem.compute(gobs)
    Mcoo        = jplephem.object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Mp   = jplephem.barycentric_object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Sp   = jplephem.barycentric_object_track("Sun", int(mjd), float(mjd%1), 1, 0.0)
    res  = jplephem.object_doppler("Moon", int(mjd), mjd%1, 1, 0.0)
    lunation,moon_state,moonsep,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,ra,dec)
    refvel = bcvel_baryc + moonvel
    print '\t\tRadial Velocity of sacttered moonlight:',refvel
 
    sci_fits        = dirout + fsim.split('/')[-1][:-4]+'spec_'+str(int(chip))+'.fits'
    sci_fits_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple_'+str(int(chip))+'.fits'
    P_fits          = dirout + 'P_' + fsim.split('/')[-1][:-4]+'spec_'+str(int(chip))+'.fits'

    if ( os.access(sci_fits,os.F_OK) == False ) or ( os.access(sci_fits_simple,os.F_OK) == False ) or \
       ( force_sci_extract ):
	
        print "\t\tNo previous extraction or extraction forced for science file", fsim, "extracting..."

	P = GLOBALutils.obtain_P(data,c_new,ext_aperture,ronoise,\
                                    gain,NSigma_Marsh, S_Marsh, \
				    N_Marsh, Marsh_alg, min_extract_col,\
				    max_extract_col, npools)

	sci_Ss = GLOBALutils.simple_extraction(data,c_new,ext_aperture,min_extract_col,max_extract_col,npools)
	sci_S  = GLOBALutils.optimal_extraction(data,P,c_new,ext_aperture,ronoise,gain,S_Marsh,NCosmic_Marsh,\
		 min_extract_col,max_extract_col,npools)
            
        if (os.access(sci_fits,os.F_OK)):
            os.remove( sci_fits )
        if (os.access(sci_fits_simple,os.F_OK)):
            os.remove( sci_fits_simple )
	if (os.access(P_fits,os.F_OK)):
            os.remove( P_fits )

        hdu = pyfits.PrimaryHDU( sci_S )
        hdu.writeto( sci_fits )
        hdu = pyfits.PrimaryHDU( sci_Ss )
        hdu.writeto( sci_fits_simple )
        hdu = pyfits.PrimaryHDU( P )
        hdu.writeto( P_fits )
	

    else:
        print '\t\t'+fsim+" has already been extracted, reading in product fits files..."
        sci_S  = pyfits.getdata( sci_fits )
        sci_Ss = pyfits.getdata( sci_fits_simple )

    fout = 'proc/'+ obname + '_' + h[0].header['DATE-OBS'] + '_' + h[0].header['UTC'] + '_' + 'sp.fits'

    if ( os.access(dirout+fout ,os.F_OK) == False ) or (force_spectral_file_build):

        spec = np.zeros((11, nord, data.shape[1]))
        hdu = pyfits.PrimaryHDU( spec )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MJD', mjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MBJD', mbjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START DATE', h[0].header['DATE-OBS'] )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START UT',  h[0].header['UTC'])
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
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARG AIRMASS',h[0].header['AIRMASS'])
	hdu = GLOBALutils.update_header(hdu,'DECKNAME',h[0].header['DECKNAME'])

	models_path = base+"../COELHO_MODELS/R_40000b/"
	deckname = h[0].header['DECKNAME']
	if deckname in ['B1','B2','B3','B4']:
		ref_RES = 72000.		
	elif deckname in ['B5','C1','C2','C3']:
		ref_RES = 48000.
	elif deckname in ['C4','C5','D1','D2']:
		ref_RES = 36000.
	else:
		ref_RES = -999

        # get ThAr closest in time
	im = np.argmin(np.absolute(np.array(ThAr_ref_dates) - mjd))
	pkl_wsol = dirout + ThAr_ref[im].split('/')[-1][:-4]+'spec_'+str(int(chip))+ '.wavsolpars.pkl'
        print "\t\tUnpickling wavelength solution from", pkl_wsol, " ..."
        wsol_dict = pickle.load(open(pkl_wsol,'r'))

        # Apply new wavelength solution including barycentric correction
        equis = np.arange( data.shape[1] )        
        for order in range(nord):
            m = order + ro0
            chebs = GLOBALutils.Calculate_chebs(equis, m, npix=data.shape[1],\
		    order0=ro0, ntotal=nord, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
            WavSol = lbary_ltopo * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wsol_dict['p1'],chebs,ncoef_x,ncoef_m)   
            spec[0,order,:] = GLOBALutils.ToVacuum(WavSol)
            spec[1,order,:] = sci_S[order,1, :]
	    #spec[1,order,-100:] = 0.
            spec[2,order,:] = sci_S[order,2, :]
	    spec[3,order,:] = spec[1,order,:] / sflat[order]
	    spec[4,order,:] = spec[2,order,:] * sflat[order] ** 2
	    #cont_coef = continuum.NORM_single(spec[0,order,:],spec[3,order,:],orden=3)
	    I = np.where(spec[1,order]!=0)[0]
	    cont_coef,x,y = hiresutils.get_cont(spec[0,order,I],scipy.signal.medfilt(spec[3,order,I],11),n=3)
	    #plot(spec[0,order],spec[3,order])
	    #plot(spec[0,order],np.polyval(cont_coef,spec[0,order]))
	    #cv = np.vstack((x,y)).T
	    #x_i,y_i = HIRESutils.bspline(cv, n=1000, degree=3, periodic=False)
	    #plot(x_i,y_i)
	    #show()
	    #print gfd
	    ratio = np.polyval(cont_coef, spec[0,order,:])

	    IN   = np.where(np.isnan(spec[3,order])==False)[0]
	    pred = np.ones(len(spec[3,order]))
	    if len(IN)>100:
		    predt,xu,yu = hiresutils.get_cont2(spec[0,order,IN],spec[3,order,IN])
		    #plot(spec[0,order,IN],spec[3,order,IN])
		    #print spec[0,order,IN].shape ,pred.shape 
		    #plot(spec[0,order,IN],predt,'k',linewidth=2.0)
		    #plot(spec[0,order,IN],ratio[IN],'r',linewidth=2.0)
		    pred[IN] = predt
	    else:
		'Nothing happens'
	    ratio = pred
	    #plot(spec[0,order,I],spec[3,order,I])
	    #plot(spec[0,order],ratio)
            L  = np.where( spec[1,order,:] != 0 )
            spec[5,order,:][L] = spec[3,order,:][L] / ratio[L]
	    spec[3,order,:][L] = spec[3,order,:][L] / ratio[L]
	    #spec[3,order,:][L] = spec[5,order,:][L].copy()
            spec[6,order,:][L] = spec[4,order,:][L] * (ratio[L] ** 2 )
            spec[7,order,:][L] = ratio[L]
            spec[8,order,:][L] = ratio[L] * sflat[order][L] / np.sqrt( ratio[L] * sflat[order][L] / gain + (ronoise/gain)**2 )
	    rI = np.where(spec[5,order] > 1. + 8./spec[8,order])
	    spec[5,order,rI] = 1.
	    #plot(spec[0,order],spec[5,order])
            spl           = scipy.interpolate.splrep(np.arange(len(WavSol)), WavSol,k=3)
            dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
            NN            = np.average(dlambda_dx)
            dlambda_dx    /= NN
            # clean-up of CRs in continuum-normalized spectrum. Set troublesome pixels to 1
	    medflx = np.zeros(len(spec[5,order]))
	    nI = np.where(np.isnan(spec[5,order])==False)[0]
	    nJ = np.where(np.isnan(spec[5,order])==True)[0]
            spec[6,order,:][I] = spec[8,order,:][I] ** 2

            spec[9,order,:][L] = spec[5,order,:][L] * (dlambda_dx[L] ** 1) 
            spec[10,order,:][L] = spec[6,order,:][L] / (dlambda_dx[L] ** 2)
	#show()
	#for order in range(nord):
	#	plot(spec[0,order],spec[5,order])
	#show()
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
            spec2 = spec.copy()
            if ref_RES>45000:
                Rx = np.around(1./np.sqrt(1./40000.**2 - 1./ref_RES**2))
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
                if vsini != -999:
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
                lowess = robjects.r("lowess")
                approx = robjects.r("approx")
                Temp = lowess(vels,yy,f=0.4,iter=10)
                pred = np.array( approx(Temp[0],Temp[1],xout=vels, method="linear", rule=2) )[1]
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
                pred = np.array( approx(Temp[0],Temp[1],xout=vels, method="linear", rule=2) )[1]
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

		bspan = GLOBALutils.calc_bss(vels,xc_av)
		SP = bspan[0]

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


	    airmass  = h[0].header['AIRMASS']
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
	    hdu = GLOBALutils.update_header(hdu,'INST', 'HIRES')
	    hdu = GLOBALutils.update_header(hdu,'RESOL', ref_RES)
	    hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')
	    hdu = GLOBALutils.update_header(hdu,'XC_MIN', XC_min)
	    hdu = GLOBALutils.update_header(hdu,'BJD_OUT', bjd_out)
            line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f   hires   ceres   %8d %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
                      (obname, bjd_out, RV, RVerr2, BS, BSerr, ref_RES, T_eff_epoch, logg_epoch, Z_epoch, vsini_epoch, XC_min, disp_epoch,\
		       TEXP, SNR_5130_R, ccf_pdf)
	    f_res.write(line_out)
	    if (os.access( dirout + fout,os.F_OK)):
            	os.remove( dirout + fout)
            hdu.writeto( dirout + fout )

f_res.close()

