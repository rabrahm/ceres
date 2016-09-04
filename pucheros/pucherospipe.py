import sys
from pylab import *

base = '../'

sys.path.append(base+"utils/Correlation")
sys.path.append(base+"utils/OptExtract/")
sys.path.append(base+"utils/GLOBALutils")

baryc_dir= base+'utils/JPLEphemx/'
sys.path.append(baryc_dir)
ephemeris='DEc403'

import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ceres modules
import pucherosutils
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
rpy2.robjects.numpy2ri.activate()

import rpy2.robjects.numpy2ri
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


args   = parser.parse_args()
dirin            = args.directorio
avoid_plot       = args.avoid_plot
dirout           = args.dirout
DoClass          = args.do_class
JustExtract      = args.just_extract
npools           = int(args.npools)
object2do        = args.o2do
reffile          = args.reffile

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
force_pre_process = False
force_flat_extract = False
force_thar_extract = False
force_tharxc = False
force_thar_wavcal = False
force_bkg = False
force_sci_extract = False
force_stellar_pars = False
compute_rv = True
dark_corr = True
back_corr = True

Inverse_m = True
use_cheby = True

MRMS = 800 # max rms in m/s, global wav solution

trace_degree  = 4
Marsh_alg     = 0
ext_aperture  = 3
NSigma_Marsh  = 5
NCosmic_Marsh = 5
S_Marsh       = 0.4
N_Marsh       = 3 # grado polinomio
min_extract_col = 20
max_extract_col = 1000

ncoef_x            = 4
ncoef_m            = 6
npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2

models_path = base+"../COELHO_MODELS/R_40000b/"    # path to the synthetic models 
order_dir = "wavcals/"
n_useful  = 41 # up to which order do we care?

def gettime(hd):
	date = hd['DATE-OBS']
	ye = float(date[:4])
	mo = float(date[5:7])
	da = float(date[8:10])
	ho = float(date[11:13])
	mi = float(date[14:16])
	se = float(date[17:])
	jd = pucherosutils.jd(ye,mo,da,ho,mi,se)

	return jd - 4./24.

# file containing the log
log = dirout+'night.log'

print "\n\n\tPUCHEROS ESO0.5m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'
# classify all the images according to its imagetype
biases, flats, img_flats, fib_flats, objects, ThAr_ref, darks = pucherosutils.FileClassify(dirin,log)
#print biases
if len(flats) > 2:
	have_flat = True
	JustExtract = False
else:
	print 'Warning: No sufficient fiber flats found ...'
	have_flat = False
	JustExtract = True

have_darks = False
if len(darks)>0:
   have_darks = True

if len(fib_flats) > 0:
        have_fib_flats = True
else:
        have_fib_flats = False

if ( (os.access(dirout+'trace.pkl',os.F_OK) == False) or \
(os.access(dirout+'MasterBias.fits',os.F_OK) == False) or \
(force_pre_process) ):
        print "\tNo previous pre-processing files or found"
	pre_process = 1
else:
        print "\tPre-processing files found, going straight to extraction"
	pre_process = 0

RON,GAIN = pucherosutils.get_rg()

if (pre_process == 1):
        print "\tGenerating Master calibration frames..."
	# median combine Biases
	MasterBias, roB, gaB = pucherosutils.MedianCombine(biases, False, dirout+'MasterBias.fits')
	hdu = pyfits.PrimaryHDU( MasterBias )
	if (os.access(dirout+'MasterBias.fits',os.F_OK)):
		os.remove(dirout+'MasterBias.fits')
	hdu.writeto(dirout+'MasterBias.fits')
        print "\t\t-> Masterbias: done!"
	# median combine Darks
	
	if len(darks)!=0:
		dark_times = []
		each_time = []

		for dark in darks:	
			dh = pyfits.getheader(dark)
			dtime = dh['EXPTIME']
			each_time.append(dtime)
			esta = False
			for t in dark_times:
				if dtime == t:
					esta = True
					break
				else:
					esta = False
			if esta == False:
				dark_times.append(dtime)

		MasDarl = []
	
		for t in dark_times:
			sirven = []
			i = 0
			while i < len(darks):
				if each_time[i] == t:
					sirven.append(darks[i])
				i+=1
			Mdark, roD, gaF = pucherosutils.MedianCombine(sirven, True, dirout+'MasterBias.fits')
			nMdark = 'MasterDark_'+str(t)+'s.fits'
			if (os.access(dirout+nMdark,os.F_OK)):
				os.remove(dirout+nMdark)
			hdu = pyfits.PrimaryHDU( Mdark )
			hdu.header.update('EXPTIME',t)
			hdu.writeto(dirout+nMdark)
			MasDarl.append(dirout+nMdark)

		darks_dict = {'d_names': MasDarl, 'd_times':dark_times}
		pickle.dump( darks_dict, open( dirout+"darks.pkl", 'w' ) )
		print "\t\t-> Masterdarks: done!"

	else:
		dark_times = []
		MasDarl = []
		darks_dict = {'d_names': MasDarl, 'd_times':dark_times}
		pickle.dump( darks_dict, open( dirout+"darks.pkl", 'w' ) )
		print "\t\t-> 0 Masterdarks found!"

	# median combine list of flats
	if len(flats)>0:
		hd = pyfits.getdata(flats[0])
		Flat, roF, gaF = pucherosutils.MedianCombine(flats, True, dirout+'MasterBias.fits', dark_bo=have_darks, dlist=MasDarl)
		hdu = pyfits.PrimaryHDU( Flat )
		if (os.access(dirout+'MasterFlat.fits',os.F_OK)):
			os.remove(dirout+'MasterFlat.fits')
		hdu.writeto(dirout+'MasterFlat.fits')

	print "\t\t-> Masterflats: done!"
        print "\tTracing echelle orders..."
	h = pyfits.open(dirout+'MasterFlat.fits')[0]
	d = h.data
	d = d.T
	c_all, nord = GLOBALutils.get_them(d,ext_aperture+1,trace_degree,maxords=45)
	print '\t\t'+str(nord)+' orders found ...'
	trace_dict = {'c_all':c_all, 'nord':nord,'roF':roF,'gaF':gaF}
	pickle.dump( trace_dict, open( dirout+"trace.pkl", 'w' ) )		

else:
	h = pyfits.open(dirout+'MasterBias.fits')
	MasterBias = h[0].data
	# load orders & tracers
	trace_dict = pickle.load( open( dirout+"trace.pkl", 'r' ) )
	#print trace_dict['c_all'].shape
	c_all = trace_dict['c_all']
	nord = trace_dict['nord']
	roF  = trace_dict['roF']
	gaF  = trace_dict['gaF']
	darks_dict = pickle.load( open( dirout+"darks.pkl", 'r' ) )
	dark_times = darks_dict['d_times']
	MasDarl = darks_dict['d_names']
	h = pyfits.open(dirout+'MasterFlat.fits')[0]
	d = h.data
	d = d.T
	Flat = d.copy()

################################ ThAr spectra extraction & calibration ##################################################
thtimes   = []
nThAr_ref = []
nthtimes  = []
print '\n\tExtraction of ThAr calibration frames:'
#force_thar_extract = True
for thar in ThAr_ref:
	#print thar
	h = pyfits.open(thar)
	mjd,mjd0 = pucherosutils.mjd_fromheader(h)
	hdth     = pyfits.getheader(thar)
	thtime   = mjd
	thtimes.append( thtime )
	
	dthar = pyfits.getdata( thar ) - MasterBias
	dthar = dthar.T

	#force_thar_extract=False
	thar_fits_simple = dirout+'ThAr'+hdth['DATE-OBS'][:10]+'_'+hdth['DATE-OBS'][11:13]+'-'+hdth['DATE-OBS'][14:16]+'-'+hdth['DATE-OBS'][17:]+'.spec.simple.fits.S'

	if ( os.access(thar_fits_simple,os.F_OK) == False ) or (force_thar_extract):

		print "\t\tNo previous extraction or extraction forced for ThAr file", thar, "extracting..."
		thar_Ss = GLOBALutils.simple_extraction(dthar,c_all,ext_aperture,\
                                                  0,1023,npools)		
		thar_Ss = thar_Ss[::-1]
		
		# save as fits file
		if (os.access(thar_fits_simple,os.F_OK)):
			os.remove( thar_fits_simple )

		hdu = pyfits.PrimaryHDU( thar_Ss )
		hdu.writeto( thar_fits_simple )
				
	else:
	        print "\t\tThAr file", thar, "all ready extracted, loading..."
		thar_Ss = pyfits.getdata(thar_fits_simple)

numt = 0
print "\n\tWavelength solution of ThAr calibration spectra:"
for thar in ThAr_ref:
	h = pyfits.open(thar)
	mjd,mjd0 = pucherosutils.mjd_fromheader(h)
	hdth     = pyfits.getheader(thar)
	thtime   = mjd
	thar_fits_simple = dirout+'ThAr'+hdth['DATE-OBS'][:10]+'_'+hdth['DATE-OBS'][11:13]+'-'+hdth['DATE-OBS'][14:16]+'-'+hdth['DATE-OBS'][17:]+'.spec.simple.fits.S'
	thar_fits_simple_wav = dirout+'ThAr'+hdth['DATE-OBS'][:10]+'_'+hdth['DATE-OBS'][11:13]+'-'+hdth['DATE-OBS'][14:16]+'-'+hdth['DATE-OBS'][17:]+'.spec.wav.fits.S'
	wavsol_pkl = dirout+'ThAr'+hdth['DATE-OBS'][:10]+'_'+hdth['DATE-OBS'][11:13]+'-'+hdth['DATE-OBS'][14:16]+'-'+hdth['DATE-OBS'][17:]+'.wavsolpars.pkl'

	#force_thar_wavcal = True
	if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
		thar_Ss = pyfits.getdata(thar_fits_simple)
		print "\t\tWorking on ThAr file", thar
	
		lines_thar = thar_Ss.copy()
		delta = 0.0
		force_ofind=True
		if numt == 0:
			if os.access(dirout+'order_find.pkl',os.F_OK)==False or force_ofind:
				maxes = 0
				or14 = 0
				for order in range(len(lines_thar)):
					ccf_max, deltar = GLOBALutils.cor_thar(lines_thar[order],span=15,filename=order_dir+'gorder14.dat')
					#print order, ccf_max
					if ccf_max > maxes:
						maxes = ccf_max
						or14  =  order
						delta = deltar

				or0 = or14 - 14
				or40 = or14 + 26
				#print or14
				if or0 >= 0:
					orwa = 0
				else:
					orwa = - or0
					or0  = 0

				if or40 <= nord - 1:
					orwb = 40
				else:
					orwb = 40 - (or40 - nord - 1)
					or40 = nord - 1
				#print or0,or40, orwa, orwb
				
				pdict = {'orwa':orwa, 'or0':or0, 'orwb':orwb, 'or40':or40, 'delta':delta}
				pickle.dump( pdict, open( dirout+'order_find.pkl', 'w' ) )
			else:
				pdict = pickle.load(open(dirout+'order_find.pkl','r'))
				orwa = pdict['orwa']
				or0 = pdict['or0']
				orwb = pdict['orwb']
				or40 = pdict['or40']
				delta = pdict['delta']

		iv_thar = 1/((lines_thar/GAIN) + (RON**2/GAIN**2))

		All_Pixel_Centers = np.array([])
		All_Wavelengths = np.array([])
		All_Orders = np.array([])
		All_Centroids = np.array([])
		All_Sigmas = np.array([])
		All_Intensities = np.array([])
		All_Residuals = np.array([])
		All_Sigmas = np.array([])
		orre = or0
		order = orwa

		OK = []
		OW = []
		nup = or40 - or0 + 1
		trans = np.zeros([nup,4])
		while orre <= or40:
			order_s = str(order)
			if (order < 10):
				order_s = '0'+str(order)

			thar_order_orig = lines_thar[orre,:]
			IV = iv_thar[orre,:]
			wei = np.sqrt( IV )
			#bkg = utils.Lines_mBack(thar_order_orig, IV, thres_rel=3)
			#thar_order = thar_order_orig - bkg
			thar_order = thar_order_orig
            		coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths,\
            		rms_ms, residuals, centroids, sigmas, intensities \
                		= GLOBALutils.Initial_Wav_Calibration(order_dir+'gorder'+order_s+'.dat',\
                                                      thar_order,order,wei,rmsmax=5000000,\
                                                      minlines=6,FixEnds=True,Dump_Argon=False,\
                                                      Dump_AllLines=True, Cheby=use_cheby)
			
			if (order == 20): 
				if (use_cheby):
					Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 511, len(thar_order) )
				else:
					Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )
			#print residuals
			All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
			All_Wavelengths = np.append( All_Wavelengths, wavelengths )
			All_Orders = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order )
			All_Centroids = np.append( All_Centroids, centroids)
			All_Sigmas = np.append( All_Sigmas, sigmas)
			All_Intensities = np.append( All_Intensities, intensities )
			All_Residuals = np.append( All_Residuals, residuals)
			All_Sigmas = np.append( All_Sigmas,sigmas)
			trans[orre,:] =  coeffs_pix2wav
			
			order += 1
			orre += 1
		"""
		JJ = np.unique(All_Orders)
		for order in JJ:
			I = np.where(All_Orders == order)[0]
			plot(All_Wavelengths[I],All_Residuals[I],'o')
		show()
		"""
		p0 = np.zeros( npar_wsol )
		oro0 = 58
		p0[0] =  (20+oro0) * Global_ZP
        	p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
            		GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders,\
                                                np.ones(All_Intensities.shape), p0, Cheby=use_cheby,\
                                                maxrms=150, Inv=Inverse_m,minlines=250,order0=oro0, \
                                                ntotal=nup,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

		#print p1
		#for odd in np.unique(G_ord):
		#	Io = np.where(G_ord == odd)[0]
		#	plot(G_wav[Io], G_res[Io],'.')
		#show()
		#print gfd
		thar_wav_Ss = np.zeros( (2,nup,dthar.shape[1]) )
		equis = np.arange( np.shape(thar_wav_Ss)[2] ) 
		order = orwa
		orre = 0
		
		while orre < nup:
			m = order + oro0
			chebs = GLOBALutils.Calculate_chebs(equis, m, order0=oro0,ntotal=nup, npix=len(thar_order), Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
			thar_wav_Ss[0,orre,:] = GLOBALutils.ToVacuum( (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(p1,chebs,ncoef_x,ncoef_m) )[::-1]
			thar_wav_Ss[1,orre,:] = thar_Ss[orre][::-1]
			orre += 1
			order+=1

		if (os.access(thar_fits_simple_wav,os.F_OK)):
			os.remove( thar_fits_simple_wav )

		hdu = pyfits.PrimaryHDU( thar_wav_Ss )
		hdu.writeto( thar_fits_simple_wav )

		#if rms_ms/np.sqrt(NL)<50:
		nThAr_ref.append(thar)
		nthtimes.append(thtime)
		pdict = {'p1':p1, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms, 'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Sigmas':All_Sigmas, 'trans':trans, 'or0':or0, 'orwa':orwa, 'oro0':oro0}
		pickle.dump( pdict, open( wavsol_pkl, 'w' ) )
		numt+=1

pdict = pickle.load(open(dirout+'order_find.pkl','r'))
orwa = pdict['orwa']
or0 = pdict['or0']
orwb = pdict['orwb']
or40 = pdict['or40']
delta = pdict['delta']
nup = or40 - or0 + 1
#ThAr_ref= nThAr_ref
#print ThAr_ref
ThAr_ref = np.array(ThAr_ref)
thtimes = np.array(thtimes)
I = np.argsort(thtimes)

thtimes = thtimes[I]
ThAr_ref= ThAr_ref[I]
ThAr_ref = list(ThAr_ref)

min_rms = 100000
thar_min = ThAr_ref[0]
for thar in ThAr_ref:
	hdth = pyfits.getheader(thar)
	wavsol_pkl = dirout+'ThAr'+hdth['DATE-OBS'][:10]+'_'+hdth['DATE-OBS'][11:13]+'-'+hdth['DATE-OBS'][14:16]+'-'+hdth['DATE-OBS'][17:]+'.wavsolpars.pkl'
	wavsol = pickle.load( open( wavsol_pkl, 'r' ) )
	if wavsol['rms_ms'] < min_rms:
		min_rms = wavsol['rms_ms']
		thar_min = thar
thar_min = ThAr_ref[0]
hdth = pyfits.getheader(thar_min)
thtime = gettime(hdth)
wavsol_pkl = dirout+'ThAr'+hdth['DATE-OBS'][:10]+'_'+hdth['DATE-OBS'][11:13]+'-'+hdth['DATE-OBS'][14:16]+'-'+hdth['DATE-OBS'][17:]+'.wavsolpars.pkl'
wavsol = pickle.load( open( wavsol_pkl, 'r' ) )
p1_ref = wavsol['p1']
best_p1 =  wavsol['p1']
f = open(dirout + 'ThAr_shifts.txt', 'w')
thshifts = []
thtimes = []
for thar in ThAr_ref:
	h = pyfits.open(thar)
	mjd,mjd0 = pucherosutils.mjd_fromheader(h)
	hdth = pyfits.getheader(thar)
	thtimes.append(mjd)
	wavsol_pkl = dirout+'ThAr'+hdth['DATE-OBS'][:10]+'_'+hdth['DATE-OBS'][11:13]+'-'+hdth['DATE-OBS'][14:16]+'-'+hdth['DATE-OBS'][17:]+'.wavsolpars.pkl'
	wavsol = pickle.load( open( wavsol_pkl, 'r' ) )
	G_pix = wavsol['G_pix']
	G_wav = wavsol['G_wav']
	G_ord = wavsol['G_ord']
	oro0 = wavsol['oro0']
	p_shift, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(G_pix, G_wav, G_ord,\
		        np.ones(G_wav.shape), p1_ref, Cheby=True, Inv=True, \
			maxrms=150,minlines=250,order0=oro0,ntotal=nup,\
			npix=h[0].data.shape[0],nx=ncoef_x,nm=ncoef_m)
	thshifts.append(p_shift)
	f.write(thar+'\t'+str(thtime)+'\t'+str((1e-6*p_shift)*299792.458)+'\n')
	#print p_shift
thshifts,thtimes = np.array(thshifts),np.array(thtimes)

Is = np.argsort(thtimes)
thtimes,thshifts = thtimes[Is],thshifts[Is]
if len(thtimes) > 3:
	thshifts = (1e-6*thshifts)*299792.458
	thtck = scipy.interpolate.splrep(thtimes,thshifts,k=3,s=0)
	ejeje = np.arange(thtimes[0],thtimes[-1],0.0001)
	ejeyy = scipy.interpolate.splev(ejeje,thtck,der=0)  
	#if not avoid_plot:
	#	plot(thtimes-56413,thshifts,'ro')
	#	plot(ejeje-56413,ejeyy)
	#	show()
else:
	thtck = 0
f.close()

#print ThAr_ref
c_all = trace_dict['c_all'][::-1]
c_all = c_all[or0:or40+1]
c_all = c_all[::-1]
print '\n\tExtraction of Flat calibration frames:'
if have_flat:
	# names of extracted Masterflat (optimal & simple) & of the P flat matrix
	S_flat_fits = dirout+'Masterflat.spec.fits'
	S_flat_fits_simple = dirout+'Masterflat.spec.fits.S'
	flat_P = dirout+'P_Flat.fits'
	sm_flat_fits = dirout+'SmoothMasterflat.spec.fits'
	# name of the extracted background image
	bacfile = dirout+'BkgFlatImage_Flat.fits'
	if ( os.access(S_flat_fits,os.F_OK) == False )  or ( os.access(S_flat_fits_simple,os.F_OK) == False ) or (force_flat_extract) or ( os.access(flat_P,os.F_OK) == False ):

		print "\t\tNo previous extraction or extraction forced for flat file",  "extracting..."
		#Flat = utils.invert(Flat)
		
		Flat  = pyfits.getdata(dirout+'MasterFlat.fits')
		Flat  = Flat.T
		if (os.access(bacfile,os.F_OK))== False or True:
		    Centers = np.zeros((len(c_all),Flat.shape[0]))
		    for i in range(c_all.shape[0]):
		        Centers[i,:]=scipy.polyval(c_all[i,:],np.arange(len(Centers[i,:])))
		    bac = GLOBALutils.get_scat(Flat,Centers,span=5)
		    hdbac = pyfits.PrimaryHDU( bac )
		    if os.access(bacfile,os.F_OK):
			os.system('rm -r '+bacfile)
		    hdbac.writeto(bacfile)
		Flat -= bac

		# Determination the P matrix
		if os.access(flat_P,os.F_OK):
			P = pyfits.getdata(flat_P)
		else:
			P = GLOBALutils.obtain_P(Flat,c_all,ext_aperture,roF,\
                                    gaF,NSigma_Marsh, S_Marsh, \
				    N_Marsh, Marsh_alg, 0,1023, npools)
			hdu = pyfits.PrimaryHDU(P)
			hdu.writeto(flat_P)

		flat_S  = GLOBALutils.optimal_extraction(Flat,P,c_all,ext_aperture,roF,\
			  gaF,S_Marsh,NCosmic_Marsh,min_extract_col,max_extract_col,npools)
		sm_flat = flat_S.copy()
		flat_Ss = GLOBALutils.simple_extraction(Flat,c_,ext_aperture,\
                                                  min_extract_col,max_extract_col,npools)

		flat_S  = flat_S[::-1]
		flat_Ss = flat_Ss[::-1]
		sm_flat = sm_flat[::-1]
		# save as fits file
		if (os.access(S_flat_fits,os.F_OK)):
			os.remove( S_flat_fits )
		if (os.access(S_flat_fits_simple,os.F_OK)):
			os.remove( S_flat_fits_simple )
		if (os.access(flat_P,os.F_OK)):
			os.remove( flat_P )
		if (os.access(sm_flat_fits,os.F_OK)):
			os.remove( sm_flat_fits )
	
		hdu = pyfits.PrimaryHDU( flat_S )
		hdu.writeto( S_flat_fits )
		hdu = pyfits.PrimaryHDU( flat_Ss )
		hdu.writeto( S_flat_fits_simple )
		hdu = pyfits.PrimaryHDU( P )
		hdu.writeto( flat_P )
		hdu = pyfits.PrimaryHDU( sm_flat )
		hdu.writeto( sm_flat_fits )

	# recover Flat spectra and P matrix	
	else:
		flat_S = pyfits.getdata(S_flat_fits)
		flat_Ss = pyfits.getdata(S_flat_fits_simple)
		P = pyfits.getdata(flat_P)
		sm_flat = pyfits.getdata( sm_flat_fits )

##################################################### Science images extraction ######################################################################
sm_flat, norms = GLOBALutils.FlatNormalize_single( sm_flat[:,1,:], mid=int(0.5*sm_flat.shape[2]),span=200)

new_list = []
new_list_obnames = []

for i in range(len(objects)):
    fsim = objects[i]
    obname = pucherosutils.search_name(fsim)
    if (object2do == 'all'):
        new_list.append(fsim)
        new_list_obnames.append( obname )
    else:
        if (obname == object2do):
            new_list.append(fsim)
            new_list_obnames.append( obname )

print '\n\tThe following targets will be processed:'
for nlisti in range(len(new_list)):
	print '\t\t'+new_list_obnames[nlisti]

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


fin_spec = []
have_specph = False
for obj in new_list:
	
	know_moon = False
	if obj.split('/')[-1] in spec_moon:
	    I = np.where(obj.split('/')[-1] == spec_moon)[0]
	    know_moon = True
	    here_moon = use_moon[I]

	print "\t--> Working on image: ", obj
	index1,index2 = -1,-1
	h = pyfits.open(obj)
	mjd,mjd0 = pucherosutils.mjd_fromheader(h)
	hd = pyfits.getheader(obj)
	nombre = pucherosutils.search_name(obj)
	exptime = hd['EXPTIME']
	
	print "\t\tObject name:",nombre
	
	altitude    = 1450.
	latitude    = -(33. + 16./60. + 9./3600.)
    	longitude   = -(70. + 32./60. + 4./3600.)
	known_coords = False	
	sp,ra,dec,known_coords = pucherosutils.get_coords(nombre,mjd)
	epoch = 2000.

	ra2,dec2 = GLOBALutils.getcoords(obname,mjd,filen=reffile)
	if ra2 !=0 and dec2 != 0:
		ra = ra2
		dec = dec2
		known_coords = True
	elif known_coords:
		print '\t\tUsing the coordinates from symbad query.'
	else:
		print '\t\tUnknown coordinate for this object.'
	
	if known_coords:
		iers                    = GLOBALutils.JPLiers( baryc_dir, mjd-999.0, mjd+999.0 )
		obsradius, R0           = GLOBALutils.JPLR0( latitude, altitude)
		obpos                   = GLOBALutils.obspos( longitude, obsradius, R0 )
	    	jplephem.set_ephemeris_dir( baryc_dir , ephemeris )
	   	jplephem.set_observer_coordinates( obpos[0], obpos[1], obpos[2] )
		res = jplephem.doppler_fraction(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
		lbary_ltopo = 1.0 + res['frac'][0]
		bcvel_baryc = ( lbary_ltopo - 1.0 ) * 2.99792458E5
		print "\t\tBarycentric velocity:", bcvel_baryc
		res = jplephem.pulse_delay(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
		mbjd = mjd + res['delay'][0] / (3600.0 * 24.0)
		
	else:
		bcvel_baryc = 0.
		mbjd = mjd

	# Moon Phase Calculations
	gobs      = ephem.Observer()  
	gobs.name = 'ODUC'  
	gobs.lat  = rad(latitude)  # lat/long in decimal degrees  
	gobs.long = rad(longitude)

	gobs.date = hd['DATE-OBS'].replace('T',' ')

	mephem = ephem.Moon()
	mephem.compute(gobs)
	Mcoo = jplephem.object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
	Mp   = jplephem.barycentric_object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
	Sp   = jplephem.barycentric_object_track("Sun", int(mjd), float(mjd%1), 1, 0.0)
	res  = jplephem.object_doppler("Moon", int(mjd), mjd%1, 1, 0.0)
	lunation,moon_state,moonsep,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,ra,dec)
	refvel = bcvel_baryc + moonvel
	print '\t\tRadial Velocity of sacttered moonlight:',refvel


	#print utils.search_name(obj)
	nama = pucherosutils.search_name(obj)+'_'+hd['DATE-OBS'][:10]+'_'+hd['DATE-OBS'][11:13]+'-'+hd['DATE-OBS'][14:16]+'-'+hd['DATE-OBS'][17:]		
	obj_fits = dirout+nama+'.spec.fits.S'
	bkg_obj_fits = dirout+'Bkg_'+nama+'.fits'
	obj_fits_simple = dirout+nama+'.spec.simple.fits.S'
	

	if ( os.access(obj_fits,os.F_OK) == False )  or \
	   ( os.access(obj_fits_simple,os.F_OK) == False ) or \
	   (force_sci_extract):
		print "\t\tNo previous extraction or extraction forced for science file", obj, "extracting..."

		dat = pyfits.getdata(obj)
		dat -= MasterBias
		if len(MasDarl)>0 and dark_corr:
			dat -= pucherosutils.get_dark(MasDarl, hd['EXPTIME'])
		dat = dat.T

		drift, c_new = GLOBALutils.get_drift(dat,P,c_all,pii=512,win=5)
		print '\t\ty drift:', drift

		P_new = GLOBALutils.shift_P(P,drift,c_new,ext_aperture)

		if back_corr:
			Centers = np.zeros((len(c_new),dat.shape[1]))
			for i in range(c_new.shape[0]):
				Centers[i,:]=scipy.polyval(c_new[i,:],np.arange(len(Centers[i,:])))
			bac   = GLOBALutils.get_scat(dat, Centers,span=5,allow_neg=True)
			dat -= bac

		if not have_flat:
			P_new = GLOBALutils.obtain_P(dat,c_new,ext_aperture,RON,\
                                    GAIN,NSigma_Marsh, S_Marsh, \
				    N_Marsh, Marsh_alg, 0,1023, npools)
		obj_Ss = GLOBALutils.simple_extraction(dat,c_new,ext_aperture,0,1023,npools)
		obj_S  = GLOBALutils.optimal_extraction(dat,P_new,c_new,ext_aperture,\
                                                    RON,GAIN,S_Marsh,NCosmic_Marsh,0,1023,npools)

		obj_Ss = obj_Ss[::-1]
		obj_S  = obj_S[::-1]

		# save as fits file
		if (os.access(obj_fits,os.F_OK)):
			os.remove( obj_fits )
		if (os.access(obj_fits_simple,os.F_OK)):
			os.remove( obj_fits_simple )
			
		hdu = pyfits.PrimaryHDU( obj_S )
		hdu.writeto( obj_fits )
		hdu = pyfits.PrimaryHDU( obj_Ss )
		hdu.writeto( obj_fits_simple )
	
	else:
		obj_S = pyfits.getdata(obj_fits)
		obj_Ss = pyfits.getdata(obj_fits_simple)

	deltat=1000
	i = 0
	while i < len(thtimes):
		if abs(thtimes[i]-mjd) < deltat:
			index1 = i
			deltat = abs(thtimes[i]-mjd)
		i+=1


	if mjd < thtimes[0]:
		print "Problem with ThAr and science times"
		index1 = 0
		index2 = 0
		indexx = 0
	elif mjd > thtimes[-1]:
		print "Problem with ThAr and science times"
		index1 = -1
		index2 = -1
		indexx = -1
	else:
		for i in range(len(thtimes)-1):
			if mjd >= thtimes[i] and mjd < thtimes[i+1]:
				index1 = i
				index2 = i+1
				if abs(mjd - thtimes[i]) < abs(mjd - thtimes[i+1]):
					indexx = i
				else:
					indexx = i+1
				break
	
	#print ThAr_ref[index1], obj, ThAr_ref[index2]
	hdth = pyfits.getheader(ThAr_ref[index1])
	wavsol_pkl = dirout+'ThAr'+hdth['DATE-OBS'][:10]+'_'+hdth['DATE-OBS'][11:13]+'-'+hdth['DATE-OBS'][14:16]+'-'+hdth['DATE-OBS'][17:]+'.wavsolpars.pkl'
	pdict = pickle.load(open(wavsol_pkl,'r'))
	global1 = pdict['p1']
	All_Pixel_Centers = pdict['G_pix']
	All_Orders = pdict['G_ord']
	All_Wavelengths = pdict['G_wav']
	rms_ms = pdict['rms_ms']
	All_Residuals = pdict['G_res']
	All_Centroids = pdict['All_Centroids']
	All_Sigmas = pdict['All_Sigmas']
	trans = pdict['trans']
	or0 = pdict['or0']
	orwa = pdict['orwa']
	oro0 = pdict['oro0']

	hdth = pyfits.getheader(ThAr_ref[index2])
	wavsol_pkl2 = dirout+'ThAr'+hdth['DATE-OBS'][:10]+'_'+hdth['DATE-OBS'][11:13]+'-'+hdth['DATE-OBS'][14:16]+'-'+hdth['DATE-OBS'][17:]+'.wavsolpars.pkl'
	pdict2 = pickle.load(open(wavsol_pkl2,'r'))
	global2 = pdict2['p1']

	hdth = pyfits.getheader(ThAr_ref[indexx])
	wavsol_pklx = dirout+'ThAr'+hdth['DATE-OBS'][:10]+'_'+hdth['DATE-OBS'][11:13]+'-'+hdth['DATE-OBS'][14:16]+'-'+hdth['DATE-OBS'][17:]+'.wavsolpars.pkl'
	pdictx = pickle.load(open(wavsol_pklx,'r'))
	globalx = pdictx['p1']

	if JustExtract:
		final = np.zeros( [6, nup,np.shape(obj_S)[2]] )
	else:
		final = np.zeros( [11, nup,np.shape(obj_S)[2]] )
	hdu = pyfits.PrimaryHDU( final )
        hdu.header.update('HIERARCH MJD', mjd)
        hdu.header.update('HIERARCH MBJD', mbjd)
        hdu.header.update('HIERARCH SHUTTER START DATE', hdth['DATE-OBS'][:10] )
        hdu.header.update('HIERARCH SHUTTER START UT',  hdth['DATE-OBS'][:11])
        hdu.header.update('HIERARCH TEXP (s)',hdth['EXPTIME'])
        hdu.header.update('HIERARCH BARYCENTRIC CORRECTION (km/s)', bcvel_baryc)
        hdu.header.update('HIERARCH (lambda_bary / lambda_topo)', lbary_ltopo)    
        hdu.header.update('HIERARCH TARGET NAME', nombre)
	try:
		hdu.header.update('HIERARCH RA',ra)
		hdu.header.update('HIERARCH DEC',dec)
		hdu.header.update('HIERARCH RA BARY',ra)
		hdu.header.update('HIERARCH DEC BARY',dec)
	except:
		None
        hdu.header.update('HIERARCH EQUINOX',2000)
        hdu.header.update('HIERARCH OBS LATITUDE',latitude)
        hdu.header.update('HIERARCH OBS LONGITUDE',longitude)
        hdu.header.update('HIERARCH OBS ALTITUDE',altitude)

	equis = np.arange( np.shape(obj_S)[2] ) 
	order = orwa
	orre = 0
	
	while orre < nup:
		
		m = order + oro0
		chebs = GLOBALutils.Calculate_chebs(equis, m, order0=oro0, npix=len(equis), ntotal=nup, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
		WavSol1 = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(global1,chebs,ncoef_x,ncoef_m)
		WavSol2 = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(global2,chebs,ncoef_x,ncoef_m)
		WavSolx = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(globalx,chebs,ncoef_x,ncoef_m)
		bestWavSol = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(best_p1,chebs,ncoef_x,ncoef_m)

		if thtimes[index2] == thtimes[index1]:
			WavSol = WavSol1
		else:
			pen = (WavSol2 - WavSol1) / (thtimes[index2] - thtimes[index1])
			coe = WavSol2 - pen* thtimes[index2]
			WavSol = pen*mjd + coe
			if thtck == 0.:
				p_shift = 0.
			else:
				p_shift = scipy.interpolate.splev(mjd, thtck, der=0)
				p_shift = 0
			WavSol = bestWavSol * (1.0 + 1.0e-6*p_shift)
		final[0,orre,:] = GLOBALutils.ToVacuum(WavSol[::-1])

		#final[0,orre,:] = WavSolx[::-1]
		if len(np.where(np.isnan(obj_S[orre,2,:][::-1]) == True) [0]) < 100:	
			final[1,orre,:] = obj_S[orre,1,:][::-1]   # ...flux in ADU (optimal extraction)...
			final[2,orre,:] = obj_S[orre,2,:][::-1]   # ...and 1/variance...
			if JustExtract==False:
				final[3,orre,:] = final[1,orre,:] / sm_flat[orre,:][::-1]   # ...flux (optimal)/flat...
				final[4,orre,:] = final[2,orre,:] * (sm_flat[orre,:][::-1] ** 2)   # ...and 1/variance * flat**2...
				cont = GLOBALutils.get_cont_single(final[0,orre,:],final[3,orre,:],final[4,orre,:], nc = 3)
				ratio = np.polyval(cont,final[0,orre])
				final[6,orre,:] = final[4,orre,:] * (ratio**2)
				final[7,orre,:] = ratio
				final[8,orre,:] = ratio * sm_flat[orre,:][::-1] / np.sqrt( ratio * sm_flat[orre,:][::-1] / GAIN + (RON/GAIN)**2 )
				final[5,orre,:] = final[3,orre,:] / ratio
				nJ = np.where(np.isnan(final[5,orre])==True)[0]
				nJ2 = np.where(np.isinf(final[5,orre])==True)[0]
				final[5,orre,nJ] = 1.0
				final[5,orre,nJ2] = 1.0
				#plot(final[8,orre])
				rI = np.where(final[5,orre] > 1. + 8./final[8,orre])
				final[5,orre,rI] = 1.
				spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), final[0,orre,:],k=3)
				dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
				NN            = np.average(dlambda_dx)
				dlambda_dx    /= NN
				final[9,orre] = final[5,orre] * (dlambda_dx ** 1) 
				final[10,orre] = final[6,orre] / (dlambda_dx ** 2)

			else:
				rw,rf = final[0,orre],final[1,orre]
				cbl = pucherosutils.fit_blaze(rw,rf)
				ratio = np.polyval(cbl,rw)
				final[3,orre] = rf/ratio
				final[4,orre] = final[2,orre]*(ratio**2)
				final[5,orre] = ratio / np.sqrt( ratio / GAIN + (RON/GAIN)**2 )
				medflx = np.zeros(len(final[3,orre,:]))
				nI = np.where(np.isnan(final[3,orre])==False)[0]
				nJ = np.where(np.isnan(final[3,orre])==True)[0]
				medflx = scipy.signal.medfilt(final[3,orre,:][nI],3)
				res = final[3,orre,:][nI] - medflx
				dev = np.sqrt(np.var(res))
				I = np.where(final[3,orre,:] > 1. + 5*dev)[0]
				final[3,orre][I]=1.
				final[3,orre][nJ]=1.

		
		order += 1
		orre += 1

        if (not JustExtract):
	    spec = final.copy()
	    #DoClass = False
	    if DoClass:
		    print '\t\tSpectral Analysis:'
		    # spectral analysis
		    # First, query SIMBAD with the object name
		    query_success = False
		    query_success,sp_type_query = GLOBALutils.simbad_query_obname(nombre)
		    # Now, query SIMBAD by coordinates if above not successful
		    if (not query_success):
		        query_success,sp_type_query = GLOBALutils.simbad_query_coords('12:00:00','00:00:00')
		    print "\t\t\tSpectral type returned by SIMBAD query:",sp_type_query

		    hdu.header.update('HIERARCH SIMBAD SPTYP', sp_type_query)

                    pars_file = dirout + obj.split('/')[-1][:-4]+'_stellar_pars.txt'

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

	    T_eff_epoch = T_eff
	    logg_epoch  = logg
	    Z_epoch     = Z
	    vsini_epoch = vsini
	    vel0_epoch  = vel0
	    hdu.header.update('HIERARCH TEFF', float(T_eff))
	    hdu.header.update('HIERARCH LOGG', float(logg))
	    hdu.header.update('HIERARCH Z', Z)
	    hdu.header.update('HIERARCH VSINI', vsini)
	    hdu.header.update('HIERARCH VEL0', vel0)

            print "\t\tRadial Velocity analysis:"
            # assign mask
	    sp_type, mask = GLOBALutils.get_mask_reffile(obname,reffile=reffile,base='../data/xc_masks/')
	    print "\t\t\tWill use",sp_type,"mask for CCF."

            # Read in mask
            ml, mh, weight = np.loadtxt(mask,unpack=True)
            ml_v = GLOBALutils.ToVacuum( ml )
            mh_v = GLOBALutils.ToVacuum( mh )
       
            # make mask larger accounting for factor ~2 lower res in CORALIE w/r to HARPS
            av_m = 0.5*( ml_v + mh_v )
            ml_v -= 5*(av_m - ml_v)
            mh_v += 5*(mh_v - av_m)
            mask_hw_kms = (GLOBALutils.Constants.c/1e3) * 0.5*(mh_v - ml_v) / av_m

            #sigma_fout = stellar_pars_dir + obname + '_' +'sigma.txt'

	    disp = GLOBALutils.get_disp(obname, reffile=reffile)
	    
	    if disp == 0:
                known_sigma = False
                if vsini != -999:
	            disp = vsini
	        else:
	            disp = 3.
            else:
                known_sigma = True

	    if disp < 10:
		disp = 5
            mask_hw_wide = av_m * disp / (GLOBALutils.Constants.c/1.0e3)
            ml_v = av_m - mask_hw_wide
            mh_v = av_m + mask_hw_wide 

            print '\t\t\tComputing the CCF...'
            cond = True
            while (cond):
                # first rough correlation to find the minimum
                vels, xc_full, sn, nlines_ccf, W_ccf = \
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight, 0, 1., vel_width=300,vel_step=1.,\
                                          spec_order=9,iv_order=10,sn_order=8,max_vel_rough=300)
                xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=0.0, Simple=True, W=W_ccf)
		#print W_ccf
                # Normalize the continuum of the CCF robustly with R     
                yy = scipy.signal.medfilt(xc_av,11)
                lowess = robjects.r("lowess")
                approx = robjects.r("approx")
                Temp = lowess(vels,yy,f=0.4,iter=10)
                pred = np.array( approx(Temp[0],Temp[1],xout=vels, method="linear", rule=2) )[1]
                xc_av_orig = xc_av.copy()
                xc_av /= pred
                vel0_xc = vels[ np.argmin( xc_av ) ]
		rvels, rxc_av, rpred, rxc_av_orig, rvel0_xc = vels.copy(), \
			xc_av.copy(), pred.copy(), xc_av_orig.copy(), vel0_xc

                xc_av_rough = xc_av
                vels_rough  = vels
                if disp > 30:
			disp = 30.
                vel_width = np.maximum( 40.0, 6*disp )

                vels, xc_full, sn, nlines_ccf, W_ccf =\
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight, vel0_xc, 1., vel_width=vel_width,vel_step=1.,\
                                          spec_order=9,iv_order=10,sn_order=8,max_vel_rough=300)
		#print W_ccf

                xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=0.0, Simple=True, W=W_ccf)
                pred = np.array( approx(Temp[0],Temp[1],xout=vels, method="linear", rule=2) )[1]
                xc_av /= pred

		if sp_type == 'M5':
			moon_sig = 5.0
		elif sp_type == 'K5':
			moon_sig = 6.6
		else:
			moon_sig = 9.0

                p1,XCmodel,p1gau,XCmodelgau,Ls2 = GLOBALutils.XC_Final_Fit( vels, xc_av ,\
                                  sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = False)

		#ldc = CoralieUtils.get_ldc(T_eff, logg, Z, 1.0, ldfile = 'lin_coe_sloan2.dat')
		#p1R, ROTmodel = CoralieUtils.XC_Final_Fit_Rot( vels, xc_av, ldc = ldc, vsini = vsini )
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

		#print 'Bisector span:', SP
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

            xc_dict = {'vels':vels,'xc_av':xc_av,'XCmodelgau':XCmodelgau,'Ls2':Ls2,'refvel':refvel,\
		       'rvels':rvels,'rxc_av':rxc_av,'rpred':rpred,'rxc_av_orig':rxc_av_orig,\
		       'rvel0_xc':rvel0_xc,'xc_full':xc_full, 'p1':p1, 'sn':sn, 'p1gau':p1gau,\
		       'p1_m':p1_m,'XCmodel_m':XCmodel_m,'p1gau_m':p1gau_m,'Ls2_m':Ls2_m,\
		       'XCmodelgau_m':XCmodelgau_m}

	    moon_dict = {'moonmatters':moonmatters,'moon_state':moon_state,'moonsep':moonsep,\
		         'lunation':lunation,'mephem':mephem,'texp':exptime}

            pkl_xc = dirout + obj.split('/')[-1][:-4]+nombre+'_XC_'+sp_type+'.pkl'
            pickle.dump( xc_dict, open( pkl_xc, 'w' ) )

	    ccf_pdf = dirout + 'proc/' + obj.split('/')[-1][:-4] + nombre + '_XCs_' + sp_type + '.pdf'

	    if not avoid_plot:
	        GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)

            SNR_5130 = np.median(spec[8,30,400:601] )
            airmass  = -999
            seeing   = -999
        
            if sp_type == 'G2':
		if T_eff < 6000:
			A = 0.06544
			B = 0.00146
			D = 0.24416
			C = 0.00181
		else:
			A = 0.09821
			B = 0.00014
			D = 0.33491
			C = 0.00113
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

            RVerr =  B + ( 1.6 + 0.2 * p1gau[2] ) * A / np.round(SNR_5130)
	    depth_fact = 1. + p1gau[0]/(p1gau[2]*np.sqrt(2*np.pi))

	    if depth_fact >= 1.:
		RVerr2 = -999.000
	    else:
                if sp_type == 'G2':
                    depth_fact = (1 - 0.62) / (1 - depth_fact)
                else:
                    depth_fact = (1 - 0.59) / (1 - depth_fact)
                RVerr2 = RVerr * depth_fact
                if (RVerr2 <= 0.009):
                    RVerr2 = 0.009

	    if RVerr2 < 0.1:
		RVerr2 = 0.1
	    BSerr = D / float(np.round(SNR_5130)) + C

	    RV     = np.around(p1gau_m[1],4)  
	    BS     = np.around(SP,4)   
            RVerr2 = np.around(RVerr2,4)
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
	    SNR_5130_R = np.around(SNR_5130*np.sqrt(2.9))

	    disp_epoch = np.around(p1gau_m[2],1)
            hdu.header.update('RV', RV)
            hdu.header.update('RV_E', RVerr2)
            hdu.header.update('BS', BS)
            hdu.header.update('BS_E', BSerr)
            hdu.header.update('DISP', disp_epoch)
            hdu.header.update('SNR', SNR_5130)
            hdu.header.update('SNR_R', SNR_5130_R)
	    hdu.header.update('INST', 'PUCHEROS')
	    hdu.header.update('RESOL', '20000')
	    hdu.header.update('PIPELINE', 'CERES')
	    hdu.header.update('XC_MIN', XC_min)
	    hdu.header.update('BJD_OUT', bjd_out)
	    
            line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f  pucheros   ceres   20000 %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
                      (obname, bjd_out, RV, RVerr2, BS, BSerr, T_eff_epoch, logg_epoch, Z_epoch, vsini_epoch, XC_min, disp_epoch,\
		       exptime, SNR_5130_R, ccf_pdf)
	    f_res.write(line_out)

	fout = 'proc/'+nama+'_final.fits'
	if (os.access( dirout + fout,os.F_OK)):
            os.remove( dirout + fout)
        hdu.writeto( dirout + fout )
	
f_res.close()
