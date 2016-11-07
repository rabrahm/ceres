import sys
from pylab import *
base = '../'
sys.path.append(base+"utils/Continuum")
sys.path.append(base+"utils/Correlation")
sys.path.append(base+"utils/OptExtract")
sys.path.append(base+"utils/GLOBALutils")

baryc_dir= base+'utils/SSEphem/'
sys.path.append(baryc_dir)
ephemeris='DEc403'

import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt

# ceres modules
import continuum
import correlation
import Marsh
import arcesutils
import GLOBALutils

# other useful modules
import pyfits
import pickle
import os
import numpy as np
import scipy
import scipy.interpolate
import string
import argparse
from math import radians as rad
from matplotlib.backends.backend_pdf import PdfPages
import ephem
import jplephem

# interface to R, not clear we'll use it
from rpy2 import robjects
import rpy2.robjects.numpy2ri
try:
	rpy2.robjects.numpy2ri.activate()
except:
	None
import rpy2.robjects.numpy2ri
r = robjects.r
r.library("MASS")

parser = argparse.ArgumentParser()
parser.add_argument('directorio')
parser.add_argument('-ofind', default='last')
parser.add_argument('-o2do',default='all')
parser.add_argument('-just_extract', action="store_true", default=False)
parser.add_argument('-do_class', action="store_true", default=False)
parser.add_argument('-avoid_plot', action="store_true", default=False)
parser.add_argument('-npools', default=1)
parser.add_argument('-reffile',default='default')
parser.add_argument('-dirout',default='default')

args = parser.parse_args()

dirin       = args.directorio
stst        = args.ofind
object2do   = args.o2do
JustExtract = args.just_extract
avoid_plot  = args.avoid_plot
npools      = int(args.npools)
reffile     = args.reffile
dirout      = args.dirout
DoClass     = args.do_class

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
force_pre_process  = False
force_flat	   = False
force_P            = False
force_flat_extract = False
force_thar_extract = False
force_thar_wavcal  = False
force_sci_extract  = False
force_bkg          = False
force_spectral_file_build = True
force_stellar_pars = True
dumpargon          = False
minlines_glob      = 500
Inverse_m          = True
use_cheby          = True
MRMS               = 100   # max rms in m/s, global wav solution

trace_degree       = 5
Marsh_alg          = 0
ext_aperture       = 5
NSigma_Marsh       = 30
NCosmic_Marsh      = 10
S_Marsh            = 0.4
N_Marsh            = 3      # grado polinomio 
min_extract_col    = 200
max_extract_col    = 1847

ncoef_x            = 4
ncoef_m            = 6
npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2
models_path = base+"data/COELHO_MODELS/R_40000b/"
order_dir   = base+"arces/wavcals/"

n_useful = 90    # up to which order do we care?
binning = 1

#############################
log = dirout+'night.log'

print "\n\n\tARCES APO3.5m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'
biases, quartzB, quartzR, science, thars, thar_dates, obnames, exptimes = arcesutils.FileClassify(dirin,log,binning)
nightlog = open(log,'r')
loglines = nightlog.readlines()
print "\tThese are all the images to proccess:"
for line in loglines:
	print '\t'+line[:-1]
print '\n'
if stst == 'last':
	if os.access(dirout+'findstar.txt',os.F_OK):
		fst = open(dirout+'findstar.txt','r')
		stst = fst.readline()
		fst.close()
	else:
		raise ValueError("There is not a previously defined standard star file!!! You have to enter one (i.e. -ofind pfs0001.fits).\n")
else:
	fst = open(dirout+'findstar.txt','w')
	fst.write(stst+'\n')
	fst.close()

if (     (os.access(dirout+'FlatB.fits',       os.F_OK) == False)	or \
	 (os.access(dirout+'FlatR.fits',       os.F_OK) == False)	or \
         (os.access(dirout+'MasterBias.fits', os.F_OK) == False)	or \
         (os.access(dirout+'trace.pkl',       os.F_OK) == False)	or \
         (force_pre_process) ):
    print "\tNo previous pre-processing files or found"
    pre_process = 1
else:
    print "\tPre-processing files found, going straight to extraction"
    pre_process = 0

if (pre_process == 1):
    # median combine Biases
    print "\tGenerating Master calibration frames..."

    if len(biases)!=0:
	MasterBias, RO_bias, GA_bias = arcesutils.MedianCombine(biases)
	print "\t\t-> Masterbias: done!"
    else:
	MasterBias, RO_bias, GA_bias = np.zeros((2048,2048)),0,1
	print "\t\t\tWarning: 0 biases"
    hdu = pyfits.PrimaryHDU( MasterBias )
    if (os.access(dirout+'MasterBias.fits',os.F_OK)):
	    os.remove(dirout+'MasterBias.fits')
    hdu.writeto(dirout+'MasterBias.fits')

    # median combine list of ob flats
    FlatR, RO_flatR, GA_flatR = arcesutils.MedianCombine(quartzR, bias = MasterBias)
    FlatB, RO_flatB, GA_flatB = arcesutils.MedianCombine(quartzB, bias = MasterBias)
    print "\t\t-> Masterflats: done!"

    # save this file for later reference
    hdu = pyfits.PrimaryHDU( FlatR )
    if (os.access(dirout+'FlatR.fits',os.F_OK)):
        os.remove(dirout+'FlatR.fits')
    hdu.writeto(dirout+'FlatR.fits')
    hdu = pyfits.PrimaryHDU( FlatB )
    if (os.access(dirout+'FlatB.fits',os.F_OK)):
        os.remove(dirout+'FlatB.fits')
    hdu.writeto(dirout+'FlatB.fits')
    
    # Find orders & traces
    print "\tTracing echelle orders..."
    c_all, nord = GLOBALutils.get_them( FlatB + FlatR, ext_aperture, trace_degree, maxords=-1,mode=1 )
    print "\t\t\t", nord, 'orders found ...'

    trace_dict = {'c_all':c_all, 'nord':nord, 'GA_bias': GA_bias, 'RO_bias': RO_bias, \
                  'GA_flatB': GA_flatB, 'RO_flatB': RO_flatB, 'GA_flatR': GA_flatR, 'RO_flatR': RO_flatR}
    pickle.dump( trace_dict, open( dirout+"trace.pkl", 'w' ) )


else:
    print '\tLoading Masterbias, Masterflat and traces'
    trace_dict = pickle.load( open( dirout+"trace.pkl", 'r' ) )
    c_all = trace_dict['c_all']
    nord = trace_dict['nord']
    GA_bias = trace_dict['GA_bias']
    RO_bias = trace_dict['RO_bias']
    GA_flatR = trace_dict['GA_flatR']
    RO_flatR = trace_dict['RO_flatR']
    GA_flatB = trace_dict['GA_flatB']
    RO_flatB = trace_dict['RO_flatB']
    # recover flats & master bias
    h = pyfits.open(dirout+'FlatR.fits')
    FlatR = h[0].data
    h = pyfits.open(dirout+'FlatB.fits')
    FlatB = h[0].data
    h = pyfits.open(dirout+'MasterBias.fits')
    MasterBias = h[0].data

Pref_fits = dirout + 'P_ref.fits'
force_Pref = False
if ( os.access(Pref_fits,os.F_OK) == False ) or (force_Pref):
	print "\n\tDetermining reference weights for optimal extraction..."
	h = pyfits.open(dirin+stst)[0]
	ron  = h.header['RDNOISE']
	gain = h.header['GAIN']
	hth = pyfits.getheader(dirin+stst)
	d = h.data
	d = arcesutils.OverscanTrim(d)
	d = arcesutils.bad_col_corr(d)
	d -= MasterBias
	c_alls = c_all.copy()
	Centers = np.zeros((len(c_alls),d.shape[1]))
	for i in range(nord):
		Centers[i,:]=scipy.polyval(c_alls[i,:],np.arange(len(Centers[i,:])))
	bkg_obj_fits = dirout + 'BKG_' + 'Pref.fits'
	if ( os.access(bkg_obj_fits,os.F_OK) == False or force_bkg):
		bkg = GLOBALutils.get_scat(d,Centers,span=4)
		if (os.access(bkg_obj_fits,os.F_OK)):
			os.remove( bkg_obj_fits )
		hdu = pyfits.PrimaryHDU( bkg )
		hdu.writeto( bkg_obj_fits )
	else:
		bkg = pyfits.getdata(bkg_obj_fits)
	d -= bkg
	if os.access(Pref_fits,os.F_OK) == False or force_P:
		P_ref = np.zeros( d.shape )
		for i in range(nord):	
			P_marsh = GLOBALutils.PCoeff( d, c_alls[i,:], ext_aperture, ron, gain, NSigma_Marsh, S_Marsh, N_Marsh, Marsh_alg , min_extract_col, max_extract_col )
			P_ref += P_marsh
		
		if (os.access(Pref_fits,os.F_OK)):
			os.remove( Pref_fits )
		hdu = pyfits.PrimaryHDU( P_ref )
		hdu.writeto( Pref_fits )
else:
	print "\tWeights for optimal extraction loaded..."
	P_ref = pyfits.getdata(Pref_fits)

# Extract Flat
print '\n\tExtraction of Flat calibration frames:'
Flat_spec_fits = dirout + 'Flat_spec.fits'
P_fits = dirout + 'P_flat.fits'
if ( os.access(Flat_spec_fits,os.F_OK) == False ) or (force_flat_extract):
	print "\t\tNo previous Flat extracted or extraction forced, extracting and saving..."
	c_alls = c_all.copy()
	Centers = np.zeros((len(c_alls),FlatB.shape[1]))

	if os.access(P_fits,os.F_OK) == False or force_P:
		P = np.zeros( FlatB.shape )
		for i in range(nord):
			P_marsh = GLOBALutils.PCoeff( FlatR+FlatB, c_alls[i,:], ext_aperture,\
                                  RO_flatR, GA_flatR, NSigma_Marsh*10, S_Marsh, 2, Marsh_alg,\
                                  min_extract_col, max_extract_col )
			P += P_marsh
		if (os.access(P_fits,os.F_OK)):
			os.remove( P_fits )
		hdu = pyfits.PrimaryHDU( P )
		hdu.writeto( P_fits )
	else:
		P = pyfits.getdata(P_fits)

	flat_S  = GLOBALutils.optimal_extraction(FlatB+FlatR,P,c_all,ext_aperture,RO_flatB,GA_flatB,\
                                       S_Marsh,10*NCosmic_Marsh,min_extract_col,max_extract_col,npools) 

	for i in range(nord):
		flat_S[i,1,:] = flat_S[i,1,:][::-1]
		flat_S[i,2,:] = flat_S[i,2,:][::-1]
			
        if (os.access(Flat_spec_fits,os.F_OK)):
            os.remove( Flat_spec_fits )
	hdu = pyfits.PrimaryHDU( flat_S )
        hdu.writeto( Flat_spec_fits )
else: 
	print "\tExtracted flat found, loading..."
	flat_S = pyfits.getdata( Flat_spec_fits )

thar_az,thar_al  = [],[]
thar_ra,thar_dec = [],[]
print '\n\tExtraction of ThAr calibration frames:'
# Extract all ThAr files
for fsim in thars:
    hthar = pyfits.open( fsim )[0]
    dthar = arcesutils.OverscanTrim( hthar.data )
    dthar = arcesutils.bad_col_corr(dthar) - MasterBias
    thar_az.append(hthar.header['TELAZ'])
    thar_al.append(hthar.header['TELALT'])
    ra  = hthar.header['RA']
    ra  = float(ra.split(':')[0]) + float(ra.split(':')[1])/60. + float(ra.split(':')[2])/3600.
    ra  = ra * 360. / 24.
    dec = hthar.header['DEC']
    dec = float(dec.split(':')[0]) + float(dec.split(':')[1])/60. + float(dec.split(':')[2])/3600.
    thar_ra.append(ra)
    thar_dec.append(dec)
    hd = pyfits.getheader(fsim)
    thar_fits = dirout + 'ARCES_' + hthar.header['DATE-OBS'] + '.ThAr.spec.fits'

    if ( os.access(thar_fits,os.F_OK) == False ) or (force_thar_extract):
        print "\t\tNo previous extraction or extraction forced for ThAr file", fsim, "extracting..."
        thar_S = np.zeros( (nord,dthar.shape[1]) )
	thar_S  = GLOBALutils.simple_extraction(dthar,c_all,ext_aperture,min_extract_col,max_extract_col,npools)
        for i in range(nord):
            thar_S[i]  = thar_S[i][::-1]
            
        # save as fits file
        if (os.access(thar_fits,os.F_OK)):
            os.remove( thar_fits )
        hdu = pyfits.PrimaryHDU( thar_S )
        hdu.writeto( thar_fits )
    else:
        print "\t\tThAr file", fsim, "all ready extracted, loading..."

thar_az, thar_al = np.array(thar_az), np.array(thar_al)
# Compute wavelength calibration of ThAr
print '\n\tWavelength solution of ThAr calibration spectra:'
sorted_thar_dates = np.argsort( thar_dates )
badind = []
for i in range(len(thars)):
    index = sorted_thar_dates[i]
    hthar = pyfits.open( thars[index] )
    wavsol_pkl = dirout + 'ARCES_' + '_' + hthar[0].header['DATE-OBS']+'.wavsolpars.pkl'
    if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
        print "\t\tWorking on initial ThAr file", thars[index] 
        
        mjd, mjd0 = arcesutils.mjd_fromheader( hthar )
        thar_fits = dirout + 'ARCES_' + hthar[0].header['DATE-OBS'] + '.ThAr.spec.fits'
        thar_S = pyfits.getdata( thar_fits )
	hd = pyfits.getheader(thar_fits)
	thar_out = dirout + 'ARCES_' + hthar[0].header['DATE-OBS'] + '.ThAr.wav.fits'

	thar_data = np.zeros((2,n_useful,thar_S.shape[1]))

        lines_thar  = thar_S.copy()
        
        All_Pixel_Centers = np.array([])
        All_Wavelengths   = np.array([])
        All_Orders        = np.array([])
        All_Centroids     = np.array([])
        All_Sigmas        = np.array([])
        All_Intensities   = np.array([])

	force_corr = True
	if os.access(dirout+'id_orders.pkl',os.F_OK) == False or force_corr:
		maxes = 0
		or41 = 0
		for order in range(len(lines_thar)):
			ccf_max, shift = GLOBALutils.cor_thar(lines_thar[order],span=10,filename=order_dir+'arces_order41.dat')
			if ccf_max > maxes:
				maxes       = ccf_max
				rough_shift = shift
				or41        =  order
		print '\t\t\tThe real echelle order 41 is order',or41
		print '\t\t\tShift in pixels:',rough_shift
		or0 = or41 - 41
		or10 = 10 + or0
			
		if or0 >= 0:
			orwa = 0
		else:
			orwa = - or0
			or0  = 0
		if os.access(dirout+'id_orders.pkl',os.F_OK):
			os.remove(dirout+'id_orders.pkl')
			
		pdict = {'or0':or0, 'orwa':orwa, 'or10':or10, 'rough_shift':rough_shift}
		pickle.dump( pdict, open(dirout+'id_orders.pkl', 'w' ) )
	else:
		pdict = pickle.load(open(dirout+'id_orders.pkl','r'))
		or10 = pdict['or10']
	
		
	order = or10
	worder = 10
	medl = []
        while order <= or10 + n_useful:
            order_s = str(worder)
            if (worder < 10):
                order_s = '0'+str(worder)
            
            thar_order_orig = lines_thar[order,:]
            thar_order      = thar_order_orig - scipy.signal.medfilt(thar_order_orig, 21)
	    
            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
                = GLOBALutils.Initial_Wav_Calibration(order_dir+'arces_order'+order_s+'.dat', thar_order, order, np.ones(len(thar_order)),rmsmax=1000, minlines=4,FixEnds=False,Dump_Argon=dumpargon,Dump_AllLines=True, Cheby=use_cheby,porder=3,rough_shift=rough_shift)
	    medl.append(GLOBALutils.Cheby_eval(coeffs_pix2wav,.5*len(thar_order),len(thar_order)))

            if (worder == 55): 
                if (use_cheby):
                    Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 0.5*len(thar_order), len(thar_order))
                else:
                    Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

            All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
            All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
            All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + worder )
            All_Centroids     = np.append( All_Centroids, centroids)
            All_Sigmas        = np.append( All_Sigmas, sigmas)
            All_Intensities   = np.append( All_Intensities, intensities )
	   
	    order +=1
	    worder +=1

        p0 = np.zeros( npar_wsol )
        p0[0] =  (55+52) * Global_ZP 
        p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
            GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers,All_Wavelengths,All_Orders,\
                                                     np.ones(All_Intensities.shape),p0,Cheby=use_cheby,\
                                                     maxrms=100, Inv=Inverse_m,minlines=minlines_glob,order0=52,\
						     ntotal=n_useful,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

	equis = np.arange( thar_S.shape[1] )
	order = 0
	while order < n_useful:
            m = order + 52 + 10
            chebs = GLOBALutils.Calculate_chebs(equis, m, order0=52, ntotal=n_useful,\
						npix=thar_S.shape[1],Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
            WavSol = (1./m) * GLOBALutils.Joint_Polynomial_Cheby(p1,chebs,ncoef_x,ncoef_m)   
            thar_data[0,order] = WavSol
	    thar_data[1,order] = lines_thar[order+or10,:]
	    order  += 1

	if (os.access(thar_out,os.F_OK)):
            os.remove( thar_out )
            
        hdu = pyfits.PrimaryHDU( thar_data )
        hdu.writeto( thar_out )

        pdict = {'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'Gobname_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                     'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Orders':All_Orders, 'All_Sigmas':All_Sigmas}
        pickle.dump( pdict, open( wavsol_pkl, 'w' ) )

    else:
        print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl
        pdict           = pickle.load(open(wavsol_pkl,'r'))
    
    if pdict['rms_ms']/np.sqrt(float(len(pdict['II']))) > 10:
	badind.append(index)

pditct2 = pickle.load(open(dirout+'id_orders.pkl','r'))
or10 = pditct2['or10']
thars = np.array(thars)
thar_dates = np.array(thar_dates)
if len(badind)>0:
	thars = np.delete(thars,badind)
	thar_dates = np.delete(thar_dates,badind)
	thar_ra = np.delete(thar_ra,badind)
	thar_dec = np.delete(thar_dec,badind)

print '\n\tStarting science frame reductions'

### start of science frame reductions ###
new_list = []
new_list_obnames = []
new_list_texp = []
for i in range(len(science)):
    fsim   = science[i]
    obname = obnames[i]
    texp   = exptimes[i]
    if (object2do == 'all'):
        new_list.append(fsim)
        new_list_obnames.append( obname )
        new_list_texp.append( texp )
    else:
        if (obname == object2do):
            new_list.append(fsim)
            new_list_obnames.append( obname )
            new_list_texp.append( texp )

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


for nlisti in range(len(new_list)):
    fsim   = new_list[ nlisti ]
    obname = new_list_obnames[ nlisti ]
    TEXP   =  new_list_texp[ nlisti ]

    know_moon = False
    if fsim.split('/')[-1] in spec_moon:
        I = np.where(fsim.split('/')[-1] == spec_moon)[0]
        know_moon = True
        here_moon = use_moon[I]

    h = pyfits.open(fsim)

    print "\n\t-->\tWorking on image: ", fsim

    # get mjd and mjd0
    mjd,mjd0 = arcesutils.mjd_fromheader(h)

    #  get gain and readnoise of object 
    ronoise = h[0].header['RDNOISE']
    gain    = h[0].header['GAIN']

    print "\t\tObject name:",obname

    # Find lambda_bary/lambda_topo using baryc
    altitude    = 2788.
    latitude    = h[0].header['LATITUDE']
    longitude   = h[0].header['LONGITUD']
    ra          = h[0].header['RA']
    ra  = float(ra.split(':')[0]) + float(ra.split(':')[1])/60. + float(ra.split(':')[2])/3600.
    ra  = ra * 360. / 24.
    dec         = h[0].header['DEC']
    dec = float(dec.split(':')[0]) + float(dec.split(':')[1])/60. + float(dec.split(':')[2])/3600.	
    epoch       = h[0].header['EQUINOX']

    ra2,dec2 = GLOBALutils.getcoords(obname,mjd,filen=reffile)
    if ra2 !=0 and dec2 != 0:
        ra = ra2
        dec = dec2
    else:
        print '\t\tUsing the coordinates found in the image header.'

    iers          = GLOBALutils.JPLiers( baryc_dir, mjd-999.0, mjd+999.0 )
    obsradius, R0 = GLOBALutils.JPLR0( latitude, altitude)
    obpos         = GLOBALutils.obspos( longitude, obsradius, R0 )
    jplephem.set_ephemeris_dir( baryc_dir , ephemeris )
    jplephem.set_observer_coordinates( obpos[0], obpos[1], obpos[2] )

    res = jplephem.doppler_fraction(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
    lbary_ltopo = 1.0 + res['frac'][0]
    bcvel_baryc = ( lbary_ltopo - 1.0 ) * 2.99792458E5
    print '\t\tBarycentric velocity:', bcvel_baryc
    res = jplephem.pulse_delay(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
    mbjd = mjd + res['delay'][0] / (3600.0 * 24.0)

    # Moon Phase Calculations
    gobs = ephem.Observer()  
    gobs.name='APO3.5'  
    gobs.lat=rad(latitude)  # lat/long in decimal degrees  
    gobs.long=rad(longitude)

    DDATE = h[0].header['DATE-OBS'].split('T')[0]
    HHOUR = h[0].header['DATE-OBS'].split('T')[1]
    Mho = HHOUR[:2]
    Mmi = HHOUR[3:5]
    Mse = HHOUR[6:]
    gobs.date = str(DDATE[:4]) + '-' +  str(DDATE[5:7]) + '-' + str(DDATE[8:]) + ' ' +  Mho + ':' + Mmi +':' + Mse
    mephem = ephem.Moon()
    mephem.compute(gobs)
    #print "Barycentric moon velocity:", bcvel_baryc_moon, mjd
    Mcoo = jplephem.object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Mp = jplephem.barycentric_object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Sp = jplephem.barycentric_object_track("Sun", int(mjd), float(mjd%1), 1, 0.0)
    res  = jplephem.object_doppler("Moon", int(mjd), mjd%1, 1, 0.0)
    lunation,moon_state,moonsep,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,ra,dec)
    refvel = bcvel_baryc + moonvel
    print '\t\tRadial Velocity of sacttered moonlight:',refvel

    sorted_indices = np.argsort( np.abs( np.array(thar_dates) - mjd ) )

    # optimally and simply extract spectra
    sci_fits = dirout + 'ARCES_' + h[0].header['DATE-OBS'] + '.' + obname +'.spec.fits'
    sci_fits_simple = dirout + 'ARCES_' + h[0].header['DATE-OBS'] +'.'+ obname +'.spec.simple.fits'
    P_fits = dirout + 'P_' +  h[0].header['DATE-OBS'] +'.'+ obname +'.fits'

    # Open file, trim, overscan subtract and MasterBias subtract
    data = h[0].data
    data = arcesutils.OverscanTrim(data)
    data = arcesutils.bad_col_corr(data)

    if ( os.access(sci_fits,os.F_OK) == False ) or ( os.access(sci_fits_simple,os.F_OK) == False ) or (force_sci_extract):
        data -= MasterBias
        #print '\t\t\tRecentering traces...'
        c_alls = c_all.copy()
        c_alls, pshift = GLOBALutils.retrace( data, c_all, span=5 )

        Centers = np.zeros((len(c_alls),data.shape[1]))
        for i in range(nord):
            Centers[i,:]=scipy.polyval(c_alls[i,:],np.arange(len(Centers[i,:])))
        bkg_obj_fits = dirout + 'BKG_' + h[0].header['DATE-OBS'] +'.'+ obname +'.fits'
        if ( os.access(bkg_obj_fits,os.F_OK) == False or force_bkg):
            bkg = GLOBALutils.get_scat(data,Centers,span=4)

            if (os.access(bkg_obj_fits,os.F_OK)):
                os.remove( bkg_obj_fits )
            hdu = pyfits.PrimaryHDU( bkg )
            hdu.writeto( bkg_obj_fits )
        else:
            bkg = pyfits.getdata(bkg_obj_fits)
        data -= bkg

        print '\t\tExtraction:'
	    #print '\t\tComputing weights...'

        if os.access(P_fits,os.F_OK) == False or force_P:
            P = GLOBALutils.obtain_P(data,c_alls,ext_aperture,ronoise,\
                        gain,NSigma_Marsh, S_Marsh, N_Marsh, Marsh_alg, min_extract_col, max_extract_col, npools)

            if (os.access(P_fits,os.F_OK)):
                os.remove( P_fits )
            hdu = pyfits.PrimaryHDU( P )
            hdu.writeto( P_fits )
        else:
            P = pyfits.getdata(P_fits)

        print "\t\t\tNo previous extraction or extraction forced for science file", fsim, "extracting..."

        sci_Ss = GLOBALutils.simple_extraction(data,c_alls,ext_aperture,\
                                                  min_extract_col,max_extract_col,npools)
        sci_S  = GLOBALutils.optimal_extraction(data,P,c_alls,ext_aperture,\
                                                   ronoise,gain,S_Marsh,NCosmic_Marsh,\
                                                   min_extract_col,max_extract_col,npools)

        for i in range(nord):
            sci_S[i,1,:] = sci_S[i,1,:][::-1]
            sci_S[i,2,:] = sci_S[i,2,:][::-1]
            sci_Ss[i,:]  = sci_Ss[i][::-1]

        # save as fits file
        if (os.access(sci_fits,os.F_OK)):
            os.remove( sci_fits )
        if (os.access(sci_fits_simple,os.F_OK)):
            os.remove( sci_fits_simple )

        hdu = pyfits.PrimaryHDU( sci_S )
        hdu.writeto( sci_fits )
        hdu = pyfits.PrimaryHDU( sci_Ss )
        hdu.writeto( sci_fits_simple )

    else:
        print '\t\t'+fsim, "has already been extracted, reading in product fits files..."
        sci_S = pyfits.getdata( sci_fits )
        sci_Ss = pyfits.getdata( sci_fits_simple )

    fout = 'proc/'+ obname + '_' + h[0].header['DATE-OBS'] + '_' + 'sp.fits'

    #Build spectra
    if ( os.access(dirout+fout ,os.F_OK) == False ) or (force_spectral_file_build):
        # initialize file that will have the spectra
        spec = np.zeros((11, n_useful, data.shape[1]))
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
        hdu = GLOBALutils.update_header(hdu,'HIERARCH EQUINOX',h[0].header['EQUINOX'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LATITUDE',h[0].header['LATITUDE'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LONGITUDE',h[0].header['LONGITUD'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS ALTITUDE',2788.)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARG AIRMASS',h[0].header['AIRMASS'])

        print '\t\tWavelength calibration:'
        # get ThAr closest in time and position
        if len(sorted_indices)>1:
            indice1 = sorted_indices[0]
            indice2 = sorted_indices[1]
            dist1ra = np.absolute( thar_ra[indice1] - ra )
            dist2ra = np.absolute( thar_ra[indice2] - ra )
            if dist1ra > 180:
                dist1ra = 360 - dist1ra
            if dist2ra > 180:
                dist2ra = 360 - dist2ra
            dist1 = (dist1ra)**2 + (thar_dec[indice1] - dec)**2
            dist2 = (dist2ra)**2 + (thar_dec[indice2] - dec)**2
            indice = indice1
            if dist2 < dist1:
                indice = indice2
        else:
            indice = sorted_indices[0]

        hthar = pyfits.open(thars[indice])
        thar_fits_ob = dirout + 'ARCES_' +  hthar[0].header['DATE-OBS'] +'.ThAr.spec.fits'
        pkl_wsol = dirout + 'ARCES__' + hthar[0].header['DATE-OBS'] +'.wavsolpars.pkl'
        print "\t\t\tUnpickling wavelength solution from", pkl_wsol, " ..."
        wsol_dict = pickle.load(open(pkl_wsol,'r'))

        # Apply new wavelength solution including barycentric correction
        equis = np.arange( data.shape[1] )        
        #print '\t\tMaking final output...' 
        for order in range(n_useful):
            m = order + 10 + 52
            chebs = GLOBALutils.Calculate_chebs(equis, m, npix=data.shape[1], order0=52, ntotal=n_useful, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
            WavSol = lbary_ltopo * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wsol_dict['p1'],chebs,ncoef_x,ncoef_m)   
            spec[0,order,:] = GLOBALutils.ToVacuum(WavSol)
            spec[1,order,:] = sci_S[order+or10,1, :]
            spec[2,order,:] = sci_S[order+or10,2, :]
            I = np.where(flat_S[order+or10,1]!=0)[0]
            spec[3,order,I] = spec[1,order,I] / flat_S[order+or10,1,I]
            spec[4,order,I] = spec[2,order,I] * flat_S[order+or10,1,I] ** 2
            nJ = np.where(np.isnan(spec[3,order])==True)[0]
            nJ2 = np.where(np.isinf(spec[3,order])==True)[0]
            spec[3,order,nJ] = 1.
            spec[3,order,nJ2] = 1.
            IJJ = np.where(spec[3,order]!=0)[0]
            if len(IJJ)>0:
                cont_coef = GLOBALutils.get_cont_single(spec[0,order],spec[3,order],spec[4,order],ll=1.5,lu=5,nc=3)   
                ratio = np.polyval(cont_coef, spec[0,order,:])
            else:
                ratio = np.ones(len(spec[3,order]))
            L  = np.where( spec[1,order,:] != 0 )
            spec[5,order,:][L] = spec[3,order,:][L] / ratio[L]
            nJ = np.where(np.isnan(spec[5,order])==True)[0]
            nJ2 = np.where(np.isinf(spec[5,order])==True)[0]
            spec[5,order,nJ] = 1.0
            spec[5,order,nJ2] = 1.0
            rI = np.where(spec[5,order] > 1. + 8./spec[8,order])
            spec[5,order,rI] = 1.
            spec[6,order,:][L] = spec[4,order,:][L] * (ratio[L] ** 2 )
            spec[7,order,:][L] = ratio[L]
            spec[8,order,:][L] = ratio[L] * flat_S[order,1][L] / np.sqrt( ratio[L] * flat_S[order,1][L] / gain + (ronoise/gain)**2 )
            spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
            dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
            NN            = np.average(dlambda_dx)
            dlambda_dx    /= NN

            spec[9,order,:][L] = spec[5,order,:][L] * (dlambda_dx[L] ** 1) 
            spec[10,order,:][L] = spec[6,order,:][L] / (dlambda_dx[L] ** 2)

    else:
        spec = pyfits.getdata(dirout+fout)

    imax,oMg = 0,0
    for i in range(spec.shape[1]):
        IM = np.where((spec[0,i]>5160)&(spec[0,i]<5200))[0]
        if len(IM)>imax:
            imax = len(IM)
            oMg  = i
    SNR_5130 = np.median(spec[8,oMg,1000:1101] )

    JustExtract = False
    if (not JustExtract):
        if DoClass:
            print '\t\tSpectral Analysis:'
            # spectral analysis
            query_success = False
            query_success,sp_type_query = GLOBALutils.simbad_query_obname(obname)
            if (not query_success):
                query_success,sp_type_query = GLOBALutils.simbad_query_coords('12:00:00','00:00:00')
            print "\t\t\tSpectral type returned by SIMBAD query:",sp_type_query

            hdu = GLOBALutils.update_header(hdu,'HIERARCH SIMBAD SPTYP', sp_type_query)
            pars_file = dirout + 'ARCES_' + h[0].header['DATE-OBS'] + '.' + obname +'_stellar_pars.txt'

            if os.access(pars_file,os.F_OK) == False or force_stellar_pars:
                print "\t\t\tEstimating atmospheric parameters:"
                T_eff, logg, Z, vsini, vel0, ccf = correlation.CCF(spec,model_path=models_path,npools=npools)
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

        # store the parameters measured for this epoch
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
        # Read in mask
        # make mask larger accounting for factor ~3 lower res in Arces w/r to HARPS
        ml, mh, weight = np.loadtxt(mask,unpack=True)
        ml_v = GLOBALutils.ToVacuum( ml )
        mh_v = GLOBALutils.ToVacuum( mh )
        av_m = 0.5*( ml_v + mh_v )
        ml_v -= 2.*(av_m - ml_v)
        mh_v += 2.*(mh_v - av_m)

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
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight, 0, lbary_ltopo, vel_width=300,vel_step=3,\
                                          spec_order=9,iv_order=10,sn_order=8,max_vel_rough=300)
            xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3.0, Simple=True, W=W_ccf)
            # Normalize the continuum of the CCF robustly with R     
            yy = scipy.signal.medfilt(xc_av,11)
            lowess = robjects.r("lowess")
            approx = robjects.r("approx")
            I = np.where(np.isnan(yy))[0]
            if len(I)==0:
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
		            GLOBALutils.XCor(spec, ml_v, mh_v, weight, vel0_xc, lbary_ltopo, vel_width=vel_width,vel_step=0.3,\
		                                  spec_order=9,iv_order=10,sn_order=8,max_vel_rough=300)

                xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3.0, Simple=True, W=W_ccf)
                pred = np.array( approx(Temp[0],Temp[1],xout=vels, method="linear", rule=2) )[1]
                xc_av /= pred

                if sp_type == 'M5':
                    moon_sig = 2.5
                elif sp_type == 'K5':
                    moon_sig = 3.3
                else:
                    moon_sig = 4.5

                p1,XCmodel,p1gau,XCmodelgau,Ls2 = GLOBALutils.XC_Final_Fit( vels, xc_av , sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = False)
                p1_m,XCmodel_m,p1gau_m,XCmodelgau_m,Ls2_m = p1,XCmodel,p1gau,XCmodelgau,Ls2

                confused = False
                ismoon = False
                moon_flag = 0

                bspan = GLOBALutils.calc_bss(vels,xc_av)
                SP = bspan[0]

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
                problem = False
            else:
                p1,p1gau = 0.,[0,0,0]
                problem = True
                cond = False

        if not problem:

            xc_dict = {'vels':vels,'xc_av':xc_av,'XCmodelgau':XCmodelgau,'Ls2':Ls2,'refvel':refvel,\
			       'rvels':rvels,'rxc_av':rxc_av,'rpred':rpred,'rxc_av_orig':rxc_av_orig,\
			       'rvel0_xc':rvel0_xc,'xc_full':xc_full, 'p1':p1, 'sn':sn, 'p1gau':p1gau,\
			       'p1_m':p1_m,'XCmodel_m':XCmodel_m,'p1gau_m':p1gau_m,'Ls2_m':Ls2_m,\
			       'XCmodelgau_m':XCmodelgau_m}

            moon_dict = {'moonmatters':moonmatters,'moon_state':moon_state,'moonsep':moonsep,\
				 'lunation':lunation,'mephem':mephem,'texp':h[0].header['EXPTIME']}

            pkl_xc = dirout + fsim.split('/')[-1][:-8]+obname+'_XC_'+sp_type+'.pkl'
            pickle.dump( xc_dict, open( pkl_xc, 'w' ) )

            ccf_pdf = dirout + 'proc/' + fsim.split('/')[-1][:-4] + obname + '_XCs_' + sp_type + '.pdf'
            if not avoid_plot:
                GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)

            airmass  = h[0].header['AIRMASS']
            seeing   = -999

            if sp_type == 'G2':
                A = 0.06544
                B = 0.00146
                D = 0.24416
                C = 0.00181
            elif  sp_type == 'K5':
                A = 0.05348
                B = 0.00147	
                D = 0.20695
                C = 0.00321
            else:
                A = 0.05348
                B = 0.00147	
                D = 0.20695
                C = 0.00321


            BSerr = D / float(np.round(SNR_5130)) + C
            RVerr2 = 0.5
            RV     = np.around(p1gau_m[1],3)  
            BS     = np.around(SP,3) 
            BSerr = np.around(BSerr,4)

            print '\t\t\tRV = '+str(RV)+' +- '+str(RVerr2)
            print '\t\t\tBS = '+str(BS)+' +- '+str(BSerr)

            bjd_out = 2400000.5 + mbjd
            T_eff_err = 100
            logg_err = 0.5
            Z_err = 0.5
            vsini_err = 2
            XC_min = np.abs(np.around(np.min(XCmodel),2))

            SNR_5130 = np.around(SNR_5130)
            SNR_5130_R = np.around(SNR_5130*np.sqrt(2.5))

            disp_epoch = np.around(p1gau_m[2],1)
            hdu = GLOBALutils.update_header(hdu,'RV', RV)
            hdu = GLOBALutils.update_header(hdu,'RV_E', RVerr2)
            hdu = GLOBALutils.update_header(hdu,'BS', BS)
            hdu = GLOBALutils.update_header(hdu,'BS_E', BSerr)
            hdu = GLOBALutils.update_header(hdu,'DISP', disp_epoch)
            hdu = GLOBALutils.update_header(hdu,'SNR', SNR_5130)
            hdu = GLOBALutils.update_header(hdu,'SNR_R', SNR_5130_R)
            hdu = GLOBALutils.update_header(hdu,'INST', 'ARCES')
            hdu = GLOBALutils.update_header(hdu,'RESOL', '40000')
            hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')
            hdu = GLOBALutils.update_header(hdu,'XC_MIN', XC_min)
            hdu = GLOBALutils.update_header(hdu,'BJD_OUT', bjd_out)

            line_out = "%-15s %18.8f %8.3f %5.3f %5.3f %5.3f arces ceres  40000 %6d %4.1f %4.1f %5.1f %3.1f %3.1f %6.1f %4d %s\n"%\
                      (obname, bjd_out, RV, RVerr2, BS, BSerr, T_eff_epoch, logg_epoch, Z_epoch, vsini_epoch, XC_min, disp_epoch,\
		       TEXP, SNR_5130_R, ccf_pdf)
            f_res.write(line_out)

    if (os.access( dirout + fout,os.F_OK)):
        os.remove( dirout + fout)
    hdu.writeto( dirout + fout )

f_res.close()
