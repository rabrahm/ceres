import sys
from pylab import *

base = '../'
sys.path.append(base+"utils/Continuum")
sys.path.append(base+"utils/Correlation")
sys.path.append(base+"utils/GLOBALutils")

baryc_dir= base+'utils/SSEphem/'
sys.path.append(baryc_dir)
ephemeris='DEc403'

import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt

# ecpipe modules
import continuum
import correlation
import fiesutils
import GLOBALutils

# other useful modules
import pyfits
import pickle
import os
import numpy as np
import scipy
import scipy.interpolate
from math import radians as rad
import argparse
import warnings
warnings.filterwarnings("ignore")

import ephem
import jplephem
from matplotlib.backends.backend_pdf import PdfPages

import statsmodels.api as sm
lowess = sm.nonparametric.lowess


parser = argparse.ArgumentParser()
parser.add_argument('directorio')
parser.add_argument('-o2do',default='all')
parser.add_argument('-just_extract', action="store_true", default=False)
parser.add_argument('-do_class', action="store_true", default=False)
parser.add_argument('-avoid_plot', action="store_true", default=False)
parser.add_argument('-npools', default=1)
parser.add_argument('-reffile',default='default')
parser.add_argument('-dirout',default='default')
parser.add_argument('-binning', default=1)
parser.add_argument('-fibre',default='default')
#parser.add_argument('-simult', action="store_true", default=False)

args = parser.parse_args()
DoClass     = args.do_class
avoid_plot  = args.avoid_plot
dirin       = args.directorio
object2do   = args.o2do
JustExtract = args.just_extract
#simult      = args.simult
npools      = int(args.npools)
reffile     = args.reffile
dirout      = args.dirout
binning     = int(args.binning)
mode        = args.fibre

if dirin[-1] != '/':
    dirin = dirin + '/'

if dirout == 'default':
    dirout = dirin[:-1]+'_red_'+mode+'_B'+str(int(binning))+'/'

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
force_flat_extract = False
force_sci_extract  = False
force_thar_extract = False	
force_tharxc       = False
force_thar_wavcal  = True
force_spectral_file_build = True
force_stellar_pars = True
dumpargon          = False
dark_corr		   = True
minlines_glob      = 700

Inverse_m          = True
use_cheby          = True
MRMS               = 90   # max rms in m/s, global wav solution

trace_degree       = 4
Marsh_alg          = 0
ext_aperture       = 4
NSigma_Marsh       = 5
NCosmic_Marsh      = 5
S_Marsh            = 0.4
N_Marsh            = 3      # grado polinomio 
min_extract_col    = 5
max_extract_col    = 2048

ncoef_x            = 4
ncoef_m            = 6
npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2

models_path = base+"data/COELHO_MODELS/R_40000b/"
order_dir   = base+"fies/wavcals/"

n_useful = 80
ro0      = 79

bacap=8
startfrom = 0
if mode == 'F1':
	resol = 25000
	MRMS  = 200.
	ext_aperture = 6
elif mode == 'F3':
	resol = 50000
elif mode == 'F4':
	resol = 67000
	startfrom = 10

if binning == 2:
	ext_aperture /= 2
	min_extract_col = int(np.around(min_extract_col/2.))
	max_extract_col = int(np.around(max_extract_col/2.))
	bacap /= bacap
#############################

# file containing the log
log = dirout+'night.log'

print "\n\n\tFIES NOT2.5m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'

biases, flats, ThAr_ref, sim_sci, ThAr_ref_dates, obnames, exptimes, darks, flats_co, flats_co_dates, ThAr_sim, ThAr_sim_dates, ThAr_co, ThAr_co_dates = fiesutils.FileClassify(dirin,log,binning=binning, mode=mode, dark_corr=dark_corr)
IS = np.argsort(ThAr_ref_dates)
ThAr_ref_dates = ThAr_ref_dates[IS]
ThAr_ref = ThAr_ref[IS]

print '\tThese are all the images to proccess:'
f = open(log)
flines = f.readlines()
for line in flines:
	print '\t'+line[:-1]
print '\n'

if (     (os.access(dirout+'Flat.fits',os.F_OK) == False) or \
         (os.access(dirout+'trace.pkl',os.F_OK) == False)  or \
         (os.access(dirout+'MasterBias.fits',os.F_OK) == False)  or \
         (force_pre_process) ):
    print "\tNo previous pre-processing files or found"
    pre_process = 1
else:
    print "\tPre-processing files found, going straight to extraction"
    pre_process = 0

if (pre_process == 1):
	# median combine Biases
	print "\tGenerating Master calibration frames..."
	if len(biases)>0:
		MasterBias, RO_bias, GA_bias = fiesutils.MedianCombine(biases, binning=binning)
	else:
		MasterBias = np.zeros(fiesutils.OverscanTrim(pyfits.getdata(flats[0]),binning=binning).shape)
		RO_bias = 0.
		GA_bias = 1.
	print "\t\t-> Masterbias: done!"
	hdu = pyfits.PrimaryHDU( MasterBias )
	if (os.access(dirout+'MasterBias.fits',os.F_OK)):
		os.remove(dirout+'MasterBias.fits')
	hdu.writeto(dirout+'MasterBias.fits')

	dark_names = []
	dark_utimes = []
	if dark_corr and len(darks)>0:
		dark_utimes, dark_times = fiesutils.get_darktimes(darks)
		for tim in dark_utimes:
			I = np.where(dark_times == tim)[0]
			dark,ron_d,gain_d = fiesutils.MedianCombine(darks[I], zero=dirout+'MasterBias.fits',binning=binning)
			hdu = pyfits.PrimaryHDU(dark)
			dark_names.append(dirout+'Dark_'+str(int(tim))+'.fits')
			if (os.access(dark_names[-1],os.F_OK)):
				os.remove(dark_names[-1])
			hdu.writeto(dark_names[-1])
		print "\t\t-> MasterDarks: done!"
	dark_names, dark_utimes = np.array(dark_names), np.array(dark_utimes)

	# median combine list of flats
	Flat, RO_fl, GA_fl = fiesutils.MedianCombine(flats, zero=dirout+'MasterBias.fits',binning=binning)
	hdu = pyfits.PrimaryHDU(Flat)
	if (os.access(dirout+'Flat.fits',os.F_OK)):
		os.remove(dirout+'Flat.fits')
	hdu.writeto(dirout+'Flat.fits')


	if len(flats_co)>0:
		Flat_co, RO_fl_co, GA_fl_co = fiesutils.MedianCombine(flats_co, zero=dirout+'MasterBias.fits',binning=binning)
		hdu = pyfits.PrimaryHDU(Flat_co)
		if (os.access(dirout+'Flat_co.fits',os.F_OK)):
			os.remove(dirout+'Flat_co.fits')
		hdu.writeto(dirout+'Flat_co.fits')
	else:
		c_co = []
		nord_co = 0
	print "\t\t-> Masterflats: done!"

	print "\tTracing echelle orders..."
	c_all,nord = GLOBALutils.get_them(Flat.T, ext_aperture, trace_degree, maxords=-1, mode=1,nsigmas=5.)
	c_all,nord = GLOBALutils.good_orders(c_all,nord,Flat.shape[1],Flat.shape[0],ext_aperture)
	print '\t\t'+ str(nord)+' orders found.'
	if len(flats_co)>0:
		c_co, nord_co = GLOBALutils.get_them(Flat_co.T, ext_aperture, trace_degree, maxords=-1, mode=1,nsigmas=5.)
		c_co, nord_co = GLOBALutils.good_orders(c_co,nord_co,Flat_co.shape[1],Flat_co.shape[0],ext_aperture)
		print '\t\t'+ str(nord_co)+' comparison orders found.'
	# pickle traces
	trace_dict = {'c_all':c_all, 'c_co':c_co, 
                  'nord':nord, 'nord_co':nord_co,
                  'GA_bias': GA_bias, 'RO_bias' : RO_bias,
                  'GA_fl': GA_fl, 'RO_fl': RO_fl,
                  'dark_names':dark_names, 'dark_utimes':dark_utimes}
	pickle.dump( trace_dict, open( dirout+"trace.pkl", 'w' ) )

else:
	trace_dict = pickle.load( open( dirout+"trace.pkl", 'r' ) )
	c_all = trace_dict['c_all']
	nord  = trace_dict['nord']
	c_co = trace_dict['c_co']
	nord_co  = trace_dict['nord_co']
	# recover GA*, RO*
	GA_bias = trace_dict['GA_bias']
	RO_bias = trace_dict['RO_bias']
	GA_fl = trace_dict['GA_fl']
	RO_fl = trace_dict['RO_fl']
	if dark_corr and len(darks)>0:
		dark_utimes = trace_dict['dark_utimes']
		dark_names = trace_dict['dark_names']
	# recover flats & master bias
	h = pyfits.open(dirout+'Flat.fits')
	Flat = h[0].data
	if len(c_co)>0:
		h = pyfits.open(dirout+'Flat_co.fits')
		Flat_co = h[0].data
	h = pyfits.open(dirout+'MasterBias.fits')
	MasterBias = h[0].data

print '\n\tExtraction of Flat calibration frames:'
if len(c_co)>0:
	c_tot = GLOBALutils.Mesh(c_all,c_co)
P_fits      = dirout + 'P.fits'
S_flat_fits = dirout +'flat.fits'
S_flat      = np.zeros((nord, 3, Flat.shape[1]) )

if ( os.access(P_fits,os.F_OK) == False ) or \
   ( os.access(S_flat_fits,os.F_OK) == False ) or \
   (force_flat_extract):
    print "\t\tNo extracted flat object spectra found or extraction forced, extracting and saving..."

    Centers = np.zeros((len(c_all),Flat.shape[0]))
    for i in range(nord):
        Centers[i,:]=scipy.polyval(c_all[i,:],np.arange(len(Centers[i,:])))

    bac = GLOBALutils.get_scat(Flat.T,Centers,span=bacap)
    fl = Flat.T - bac
    #plot(fl[:,1000])
    #plot(np.around(Centers[:,1000]).astype('int'),fl[np.around(Centers[:,1000].astype('int')),1000],'ro')
    #show()
    #print gfd
    bacfile = dirout + 'BAC_FLAT.fits'
    if (os.access(bacfile,os.F_OK)):
            os.remove( bacfile )
    hdbac = pyfits.PrimaryHDU( bac )
    hdbac.writeto(bacfile)

    print "\t\tWill extract",nord,"orders for object fibre..."
    P = GLOBALutils.obtain_P(fl,c_all,ext_aperture,RO_fl,\
                                    GA_fl,NSigma_Marsh, S_Marsh, \
                    N_Marsh, Marsh_alg, min_extract_col,\
                    max_extract_col, npools)

    if (os.access(P_fits,os.F_OK)):
        os.remove( P_fits )
    hdu = pyfits.PrimaryHDU( P )
    hdu.writeto( P_fits )

    S_flat  = GLOBALutils.optimal_extraction(fl,P,c_all,ext_aperture,RO_fl,GA_fl,\
              S_Marsh,10*NCosmic_Marsh,min_extract_col,max_extract_col,npools)
    S_flat = S_flat[::-1]
    S_flat = GLOBALutils.invert(S_flat)
    if (os.access(S_flat_fits,os.F_OK)):
        os.remove( S_flat_fits )
    hdu = pyfits.PrimaryHDU( S_flat )
    hdu.writeto( S_flat_fits )

else:
    print "\t\tExtracted flat object spectra found, loading..."
    P    = pyfits.getdata( P_fits )
    S_flat = pyfits.getdata( S_flat_fits )

if nord_co>0:
	P_co_fits      = dirout + 'P_co.fits'
	S_flat_co_fits = dirout +'flat_co.fits'
	S_flat_co      = np.zeros((nord_co, 3, Flat_co.shape[1]) )

	if ( os.access(P_co_fits,os.F_OK) == False ) or \
	   ( os.access(S_flat_co_fits,os.F_OK) == False ) or \
	   (force_flat_extract):
	    print "\t\tNo extracted flat comparison spectra found or extraction forced, extracting and saving..."

	    Centers = np.zeros((len(c_co),Flat_co.shape[0]))
	    for i in range(nord_co):
		Centers[i,:]=scipy.polyval(c_co[i,:],np.arange(len(Centers[i,:])))

	    bac = GLOBALutils.get_scat(Flat_co.T,Centers,span=bacap)
	    fl = Flat_co.T - bac
	    #plot(fl[:,1000])
	    #plot(np.around(Centers[:,1000]).astype('int'),fl[np.around(Centers[:,1000].astype('int')),1000],'ro')
	    #show()
	    #print gfd
	    bacfile = dirout + 'BAC_FLAT_CO.fits'
	    if (os.access(bacfile,os.F_OK)):
		    os.remove( bacfile )
	    hdbac = pyfits.PrimaryHDU( bac )
	    hdbac.writeto(bacfile)

	    print "\t\tWill extract",nord_co,"orders for comparison fibre..."
	    P_co = GLOBALutils.obtain_P(fl,c_co,ext_aperture,RO_fl,\
		                            GA_fl,NSigma_Marsh, S_Marsh, \
		            N_Marsh, Marsh_alg, min_extract_col,\
		            max_extract_col, npools)

	    if (os.access(P_co_fits,os.F_OK)):
		os.remove( P_co_fits )
	    hdu = pyfits.PrimaryHDU( P_co )
	    hdu.writeto( P_co_fits )

	    S_flat_co  = GLOBALutils.optimal_extraction(fl,P_co,c_co,ext_aperture,RO_fl,GA_fl,\
		      S_Marsh,10*NCosmic_Marsh,min_extract_col,max_extract_col,npools)
	    S_flat_co = S_flat_co[::-1]
	    S_flat_co = GLOBALutils.invert(S_flat_co)
	    if (os.access(S_flat_co_fits,os.F_OK)):
		os.remove( S_flat_co_fits )
	    hdu = pyfits.PrimaryHDU( S_flat_co )
	    hdu.writeto( S_flat_co_fits )

	else:
	    print "\t\tExtracted flat object spectra found, loading..."
	    P_co    = pyfits.getdata( P_co_fits )
	    S_flat_co = pyfits.getdata( S_flat_co_fits )

# Normalize flat field spectra.
S_flat_n, Snorms = GLOBALutils.FlatNormalize_single( S_flat, mid = int(.5*S_flat.shape[2]) )
if nord_co>0:
	S_flat_co_n, Snorms_co = GLOBALutils.FlatNormalize_single( S_flat_co, mid = int(.5*S_flat_co.shape[2]) )

print '\n\tExtraction of ThAr calibration frames:'

# Extract all ThAr files
for fsim in ThAr_ref:
	thar_fits        = dirout + fsim.split('/')[-1][:-4]+'spec.fits.S'
	thar_simple_fits = dirout + fsim.split('/')[-1][:-4]+'spec.simple.fits.S'
	if ( os.access(thar_simple_fits,os.F_OK) == False ) or ( force_thar_extract ):
		hthar = pyfits.open( fsim )
		dthar = fiesutils.OverscanTrim( hthar[1].data, binning=binning ) - MasterBias
		Centers = np.zeros((len(c_all),dthar.shape[0]))
		for i in range(nord):
			Centers[i,:]=scipy.polyval(c_all[i,:],np.arange(len(Centers[i,:])))
		bac = GLOBALutils.get_scat(dthar.T,Centers,span=ext_aperture,option=1,allow_neg=True)
		sdthar = dthar.T - bac
		#plot(dthar.T[:,1000])
		#plot(bac[:,1000])
		#show()
		print "\t\tNo previous extraction or extraction forced for ThAr file", fsim, "extracting..."
		thar_S  = np.zeros( (nord,3,dthar.shape[0]) )
		thar_Ss = np.zeros( (nord,dthar.shape[0]) )

		tR,tG = fiesutils.get_RONGAIN(hthar[1].header)
		#print tR,tG
		thar_Ss = GLOBALutils.simple_extraction(sdthar,c_all,ext_aperture,\
                                                  min_extract_col,max_extract_col,npools)
		thar_Ss = thar_Ss[::-1]
		thar_Ss = GLOBALutils.invert(thar_Ss)
		if (os.access(thar_fits,os.F_OK)):
			os.remove( thar_fits )
		if (os.access(thar_simple_fits,os.F_OK)):
			os.remove( thar_simple_fits )

		hdu = pyfits.PrimaryHDU( thar_S )
		hdu.writeto( thar_fits )
		hdu = pyfits.PrimaryHDU( thar_Ss )
		hdu.writeto( thar_simple_fits )
	else:
		print "\t\tThAr file", fsim, "all ready extracted, loading..."

# Extract all ThAr files
for fsim in ThAr_co:
	thar_co_fits        = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
	thar_co_simple_fits = dirout + fsim.split('/')[-1][:-4]+'spec.co.simple.fits.S'
	if ( os.access(thar_co_simple_fits,os.F_OK) == False ) or ( force_thar_extract ):
		hthar = pyfits.open( fsim )
		dthar = fiesutils.OverscanTrim( hthar[1].data, binning=binning ) - MasterBias
		Centers = np.zeros((len(c_co),dthar.shape[0]))
		for i in range(nord_co):
			Centers[i,:]=scipy.polyval(c_co[i,:],np.arange(len(Centers[i,:])))
		bac = GLOBALutils.get_scat(dthar.T,Centers,span=ext_aperture+2)
		sdthar = dthar.T - bac
		print "\t\tNo previous extraction or extraction forced for ThAr file", fsim, "extracting..."
		thar_S  = np.zeros( (nord_co,3,dthar.shape[0]) )
		thar_Ss = np.zeros( (nord_co,dthar.shape[0]) )

		tR,tG = fiesutils.get_RONGAIN(hthar[1].header)
		#print tR,tG
		thar_Ss = GLOBALutils.simple_extraction(sdthar,c_co,ext_aperture,\
                                                  min_extract_col,max_extract_col,npools)
		thar_Ss = thar_Ss[::-1]
		thar_Ss = GLOBALutils.invert(thar_Ss)
		if (os.access(thar_co_fits,os.F_OK)):
			os.remove( thar_co_fits )
		if (os.access(thar_co_simple_fits,os.F_OK)):
			os.remove( thar_co_simple_fits )

		hdu = pyfits.PrimaryHDU( thar_S )
		hdu.writeto( thar_co_fits )
		hdu = pyfits.PrimaryHDU( thar_Ss )
		hdu.writeto( thar_co_simple_fits )
	else:
		print "\t\tThAr file", fsim, "all ready extracted, loading..."

for fsim in ThAr_sim:
		thar_fits        = dirout + fsim.split('/')[-1][:-4]+'spec.fits.S'
		thar_simple_fits = dirout + fsim.split('/')[-1][:-4]+'spec.simple.fits.S'
		thar_co_fits        = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
		thar_co_simple_fits = dirout + fsim.split('/')[-1][:-4]+'spec.co.simple.fits.S'
		if ( os.access(thar_simple_fits,os.F_OK) == False ) or ( force_thar_extract ):
			hthar = pyfits.open( fsim )
			dthar = fiesutils.OverscanTrim( hthar[1].data, binning=binning ) - MasterBias
			Centers = np.zeros((len(c_tot),dthar.shape[0]))
			ccc = c_tot.copy()

			for i in range(Centers.shape[0]):
				Centers[i,:]=scipy.polyval(ccc[i,:],np.arange(len(Centers[i,:])))

			bac = GLOBALutils.get_scat(dthar.T,Centers,span=ext_aperture,option=1)
			sdthar = dthar.T - bac

			print "\t\tNo previous extraction or extraction forced for simultaneous ThAr file", fsim, "extracting..."
			thar_S  = np.zeros( (nord,3,dthar.shape[0]) )
			thar_Ss = np.zeros( (nord,dthar.shape[0]) )
			thar_S_co  = np.zeros( (nord_co,3,dthar.shape[0]) )
			thar_S_co  = np.zeros( (nord_co,3,dthar.shape[0]) )

			tR,tG = fiesutils.get_RONGAIN(hthar[1].header)
			#print tR,tG
			thar_Ss = GLOBALutils.simple_extraction(sdthar,c_all,ext_aperture,\
		                                          min_extract_col,max_extract_col,npools)
			thar_Ss_co = GLOBALutils.simple_extraction(sdthar,c_co,ext_aperture,\
		                                          min_extract_col,max_extract_col,npools)
			thar_Ss = thar_Ss[::-1]
			thar_Ss = GLOBALutils.invert(thar_Ss)
			if (os.access(thar_fits,os.F_OK)):
				os.remove( thar_fits )
			if (os.access(thar_simple_fits,os.F_OK)):
				os.remove( thar_simple_fits )

			thar_Ss_co = thar_Ss_co[::-1]
			thar_Ss_co = GLOBALutils.invert(thar_Ss_co)
			if (os.access(thar_co_fits,os.F_OK)):
				os.remove( thar_co_fits )
			if (os.access(thar_co_simple_fits,os.F_OK)):
				os.remove( thar_co_simple_fits )

			hdu = pyfits.PrimaryHDU( thar_S )
			hdu.writeto( thar_fits )
			hdu = pyfits.PrimaryHDU( thar_Ss )
			hdu.writeto( thar_simple_fits )

			hdu = pyfits.PrimaryHDU( thar_S_co )
			hdu.writeto( thar_co_fits )
			hdu = pyfits.PrimaryHDU( thar_Ss_co )
			hdu.writeto( thar_co_simple_fits )
		else:
			print "\t\tThAr file", fsim, "all ready extracted, loading..."

"""
p0    = np.zeros( npar_wsol )
#p0[0] =  (int(.5*n_useful)+ro0) * Global_ZP 
dat = np.loadtxt('wavcals/initial.txt')
p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = GLOBALutils.Fit_Global_Wav_Solution(dat[:,2], dat[:,1],\
						     dat[:,0], np.ones(len(dat[:,0])), p0, Cheby=use_cheby,       \
						     order0=79, ntotal=nord, maxrms=100000, Inv=Inverse_m, minlines=10,  \
						     npix=S_flat_n.shape[2],nx=3,nm=4)
#f = open('thar_list.txt','r')
f = open('lovis.txt','r')
lines = f.readlines()
wv,tp = [],[]
for line in lines:
	cos = line.split()
	wv.append(float(cos[0]))
	if len(cos)==4:
		tp.append(cos[3])
	else:
		tp.append(' ')
wv,tp = np.array(wv),np.array(tp)

sc =pyfits.getdata(dirout + 'FIyh230033.spec.simple.fits.S')
for i in range(sc.shape[0]):
	print i
	sorder = str(i)
	if i <10:
		sorder = '0'+sorder
	fout = open('wavcals/order_'+sorder+'.iwdat','w')
	m = i + 79
	ejx = np.arange(sc.shape[1])
	ejy = sc[i]
	chebs = GLOBALutils.Calculate_chebs(ejx, np.zeros(len(ejx))+i+79,Inverse=True,order0=79,ntotal=nord,npix=len(ejx),nx=ncoef_x,nm=ncoef_m)
	wavsol = GLOBALutils.ToVacuum(1.0/float(m) * GLOBALutils.Joint_Polynomial_Cheby(p1,chebs,ncoef_x,ncoef_m))
	print wavsol
	tck = scipy.interpolate.splrep(wavsol,ejx,k=3)
	I = np.where((wv>wavsol[50])&(wv<wavsol[-50]))[0]

	wlist,plist,tlist = [],[],[]
	write= False
	if len(I)>0:
		px = np.around(scipy.interpolate.splev(wv[I],tck)).astype('int')
		tpx,twv = [],[]
		for o in range(len(px)):
			if ejy[px[o]]> 100:
				if len(tlist) == 0:
					wlist.append(wv[I][o])
					plist.append(px[o])
					tlist.append(tp[I][o])
				else:
					if px[o] - plist[-1] < 8:
						wlist.append(wv[I][o])
						plist.append(px[o])
						tlist.append(tp[I][o])
					else:
						write = True

			if write:
				lout = str(len(tlist))
				for ix in range(len(wlist)):
					lout += '\t'+str(plist[ix])+'\t'+str(wlist[ix])
				for ix in range(len(wlist)):
					lout += '\t'+str(tlist[ix])
				lout+='\n'
				fout.write(lout)
				wlist,plist,tlist = [],[],[]
				write=False

		#plot(ejx,ejy)
		#plot(tpx,ejy[tpx],'ro')
		#show()
		#print ghj
	fout.close()

plot(G_wav,G_res,'r.')
show()
print p1
"""

print "\n\tWavelength solution of ThAr calibration spectra of object fibre:"
for fsim in ThAr_ref:
	wavsol_pkl = dirout + fsim.split('/')[-1][:-4]+'wavsolpars.pkl'
	h = pyfits.open(fsim)
	hd = pyfits.getheader(fsim)
	ron,gain = fiesutils.get_RONGAIN(h[0].header)
	if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
		print "\t\tWorking on simple ThAr file", fsim
		hthar = pyfits.open( fsim )
		mjd, mjd0 = fiesutils.mjd_fromheader2( hthar )
		thar_fits = dirout + fsim.split('/')[-1][:-4]+'spec.simple.fits.S'
		thar_S = pyfits.getdata( thar_fits )
		lines_thar  = thar_S.copy()

		All_Pixel_Centers = np.array([])
		All_Wavelengths   = np.array([])
		All_Orders        = np.array([])
		All_Centroids     = np.array([])
		All_Sigmas        = np.array([])
		All_Intensities   = np.array([])
		All_residuals     = np.array([])
		
		orders_offset, rough_shift = fiesutils.get_thar_offsets(lines_thar,binning=binning, delt_or=20)
		print 'orders_ofset:',orders_offset
		print 'rough_shift:',rough_shift

		orderi = 0
		if orders_offset < 0:
			orderi = - orders_offset
		orderf = nord - 1
		if orderf + orders_offset >= n_useful:
			orderf = n_useful - orders_offset - 1

		for order in range(orderi, orderf+1):
			#print order
			order_s = str(order+orders_offset)
			if (order + orders_offset < 10):
				order_s = '0' + str(order+orders_offset)
			f = open(order_dir+'order_'+order_s+'.iwdat','r')
			llins = f.readlines()
			if len(llins)>5:
				thar_order_orig = lines_thar[order]
				#IV              = iv_thar_ob_R[order,:]
				L = np.where(thar_order_orig != 0)[0]
				IV = 1. / (thar_order_orig / gain + (ron/gain)**2 )
				IV[L] = 0.
				wei             = np.ones(len(thar_order_orig))  #np.sqrt( IV )
				bkg = scipy.signal.medfilt(thar_order_orig,101)
				thar_order      = thar_order_orig - bkg

				coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids,sigmas, intensities \
	                = GLOBALutils.Initial_Wav_Calibration( order_dir+'order_'+order_s+'.iwdat', thar_order, order, wei, \
	                	rmsmax=500, minlines=10, FixEnds=False,Dump_Argon=dumpargon, Dump_AllLines=True, line_width=6, Cheby=use_cheby,porder=3,rough_shift=rough_shift,binning=binning,del_width=5.,do_xc=False)
				#plot(wavelengths, residuals,'.')	
				#show()		            
				if (order == int(.5*n_useful)): 	
					if (use_cheby):
						Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav,  0.5*len(thar_order), len(thar_order) )
					else:
						Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

				All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
				All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
				All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order)
				All_Centroids     = np.append( All_Centroids, centroids)
				All_Sigmas        = np.append( All_Sigmas, sigmas)
				All_Intensities   = np.append( All_Intensities, intensities )
				All_residuals     = np.append( All_residuals, residuals )
		#show()
		p0    = np.zeros( npar_wsol )
		p0[0] =  int(.5*n_useful) * Global_ZP 
		p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths,\
						     All_Orders, np.ones(All_Intensities.shape), p0, Cheby=use_cheby,       \
						     order0=ro0+orders_offset, ntotal=n_useful, maxrms=MRMS, Inv=Inverse_m, minlines=300,  \
						     npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
		pdict = {'orders_offset':orders_offset,'rough_shift':rough_shift,'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                         'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Wavelengths':All_Wavelengths, \
			 'All_Orders':All_Orders, 'All_Pixel_Centers':All_Pixel_Centers, 'All_Sigmas':All_Sigmas}

		pickle.dump( pdict, open(wavsol_pkl, 'w' ) )
	else:
		print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl
		pdict = pickle.load(open(wavsol_pkl,'r'))

print "\n\tWavelength solution of ThAr calibration spectra of comparison fibre:"
for fsim in ThAr_co:
	wavsol_pkl = dirout + fsim.split('/')[-1][:-4]+'wavsolpars.pkl'
	h = pyfits.open(fsim)
	mjd, mjd0 = fiesutils.mjd_fromheader2( h )
	hd = pyfits.getheader(fsim)
	ron,gain = fiesutils.get_RONGAIN(h[0].header)
	if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
		print "\t\tWorking on simple ThAr file", fsim
		hthar = pyfits.open( fsim )

		thar_fits = dirout + fsim.split('/')[-1][:-4]+'spec.co.simple.fits.S'
		thar_S = pyfits.getdata( thar_fits )
		lines_thar  = thar_S.copy()

		All_Pixel_Centers = np.array([])
		All_Wavelengths   = np.array([])
		All_Orders        = np.array([])
		All_Centroids     = np.array([])
		All_Sigmas        = np.array([])
		All_Intensities   = np.array([])
		All_residuals     = np.array([])
		
		orders_offset, rough_shift = fiesutils.get_thar_offsets(lines_thar,binning=binning)
		print 'orders_ofset:',orders_offset
		print 'rough_shift:',rough_shift

		orderi = 0
		if orders_offset < 0:
			orderi = - orders_offset
		orderf = nord_co - 1
		if orderf + orders_offset >= n_useful:
			orderf = n_useful - orders_offset - 1

		for order in range(orderi,orderf+1):
			#print order
			order_s = str(order+orders_offset)
			if (order + orders_offset < 10):
				order_s = '0' + str(order+orders_offset)
			f = open(order_dir+'order_'+order_s+'.iwdat','r')
			llins = f.readlines()
			if len(llins)>5:
				thar_order_orig = lines_thar[order]
				#IV              = iv_thar_ob_R[order,:]
				L = np.where(thar_order_orig != 0)[0]
				IV = 1. / (thar_order_orig / gain + (ron/gain)**2 )
				IV[L] = 0.
				wei             = np.ones(len(thar_order_orig))  #np.sqrt( IV )
				bkg = scipy.signal.medfilt(thar_order_orig,101)
				thar_order      = thar_order_orig - bkg

				coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids,sigmas, intensities \
	                = GLOBALutils.Initial_Wav_Calibration( order_dir+'order_'+order_s+'.iwdat', thar_order, order, wei, \
	                	rmsmax=500, minlines=10, FixEnds=False,Dump_Argon=dumpargon, Dump_AllLines=True, line_width=6, Cheby=use_cheby,porder=3,rough_shift=rough_shift,binning=binning,del_width=5.,do_xc=False)
				#plot(wavelengths, residuals,'.')	
				#show()		            
				if (order == int(.5*n_useful)): 	
					if (use_cheby):
						Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav,  0.5*len(thar_order), len(thar_order) )
					else:
						Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

				All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
				All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
				All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order)
				All_Centroids     = np.append( All_Centroids, centroids)
				All_Sigmas        = np.append( All_Sigmas, sigmas)
				All_Intensities   = np.append( All_Intensities, intensities )
				All_residuals     = np.append( All_residuals, residuals )
		#show()
		p0    = np.zeros( npar_wsol )
		p0[0] =  int(.5*n_useful) * Global_ZP 
		p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths,\
						     All_Orders, np.ones(All_Intensities.shape), p0, Cheby=use_cheby,       \
						     order0=ro0+orders_offset, ntotal=n_useful, maxrms=MRMS, Inv=Inverse_m, minlines=300,  \
						     npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

		I = np.argmin((ThAr_ref_dates - mjd)**2)
		pdict = {'orders_offset_co':orders_offset,'rough_shift_co':rough_shift,'p1_co':p1,'mjd_co':mjd,\
			 'G_pix_co':G_pix, 'G_ord_co':G_ord, 'G_wav_co':G_wav, 'II_co':II, 'rms_ms_co':rms_ms,\
                         'G_res_co':G_res, 'All_Centroids_co':All_Centroids, 'All_Wavelengths_co':All_Wavelengths,\
			 'All_Orders_co':All_Orders, 'All_Pixel_Centers_co':All_Pixel_Centers, 'All_Sigmas_co':All_Sigmas,\
			 'ref_thar_ob':ThAr_ref[I]}

		pickle.dump( pdict, open(wavsol_pkl, 'w' ) )
	else:
		print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl
		pdict = pickle.load(open(wavsol_pkl,'r'))

for fsim in ThAr_sim:
	wavsol_pkl = dirout + fsim.split('/')[-1][:-4]+'wavsolpars.pkl'
	h = pyfits.open(fsim)
	hd = pyfits.getheader(fsim)
	ron,gain = fiesutils.get_RONGAIN(h[0].header)
	if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
		print "\t\tWorking on sim ThAr file", fsim
		hthar = pyfits.open( fsim )
		mjd, mjd0 = fiesutils.mjd_fromheader2( hthar )
		thar_fits = dirout + fsim.split('/')[-1][:-4]+'spec.simple.fits.S'
		thar_S = pyfits.getdata( thar_fits )
		lines_thar  = thar_S.copy()

		All_Pixel_Centers = np.array([])
		All_Wavelengths   = np.array([])
		All_Orders        = np.array([])
		All_Centroids     = np.array([])
		All_Sigmas        = np.array([])
		All_Intensities   = np.array([])
		All_residuals     = np.array([])
		xcs = []
		
		orders_offset, rough_shift = fiesutils.get_thar_offsets(lines_thar,binning=binning)
		print 'orders_ofset:',orders_offset
		print 'rough_shift:',rough_shift

		orderi = 0
		if orders_offset < 0:
			orderi = - orders_offset
		orderf = nord - 1
		if orderf + orders_offset >= n_useful:
			orderf = n_useful - orders_offset - 1

		for order in range(orderi,orderf+1):
			#print order
			order_s = str(order+orders_offset)
			if (order + orders_offset < 10):
				order_s = '0' + str(order+orders_offset)
			f = open(order_dir+'order_'+order_s+'.iwdat','r')
			llins = f.readlines()
			if len(llins)>5:
				thar_order_orig = lines_thar[order]
				#IV              = iv_thar_ob_R[order,:]
				L = np.where(thar_order_orig != 0)[0]
				IV = 1. / (thar_order_orig / gain + (ron/gain)**2 )
				IV[L] = 0.
				wei             = np.ones(len(thar_order_orig))  #np.sqrt( IV )
				bkg = scipy.signal.medfilt(thar_order_orig,101)
				thar_order      = thar_order_orig - bkg

				coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids,sigmas, intensities \
	                = GLOBALutils.Initial_Wav_Calibration( order_dir+'order_'+order_s+'.iwdat', thar_order, order, wei, \
	                	rmsmax=500, minlines=10, FixEnds=False,Dump_Argon=dumpargon, Dump_AllLines=True, line_width=6, Cheby=use_cheby,porder=3,rough_shift=rough_shift,binning=binning,del_width=5.,do_xc=False)
				#plot(wavelengths, residuals,'.')	
				#show()		            
				if (order == int(.5*n_useful)): 	
					if (use_cheby):
						Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav,  0.5*len(thar_order), len(thar_order) )
					else:
						Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

				All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
				All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
				All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order)
				All_Centroids     = np.append( All_Centroids, centroids)
				All_Sigmas        = np.append( All_Sigmas, sigmas)
				All_Intensities   = np.append( All_Intensities, intensities )
				All_residuals     = np.append( All_residuals, residuals )
		#show()
		p0    = np.zeros( npar_wsol )
		p0[0] =  int(.5*n_useful) * Global_ZP 
		p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths,\
						     All_Orders, np.ones(All_Intensities.shape), p0, Cheby=use_cheby,       \
						     order0=ro0+orders_offset, ntotal=n_useful, maxrms=MRMS, Inv=Inverse_m, minlines=300,  \
						     npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
		"""
		equis = np.arange(len(thar_order))
		for order in range(lines_thar.shape[0]):
			m   = order + ro0 + orders_offset
			chebs  = GLOBALutils.Calculate_chebs(equis, m, npix=len(equis),\
				order0=ro0 + orders_offset, ntotal=nord, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
			WavSol = (1.0/float(m)) * GLOBALutils.Joint_Polynomial_Cheby(p1,chebs,ncoef_x,ncoef_m) 
			plot(WavSol, lines_thar[order])
		show()
		"""

			
		#	I = np.where(G_ord == i)[0]
		#	plot(G_pix[I],G_wav[I],'.')
		#show()

		thar_co_fits = dirout + fsim.split('/')[-1][:-4]+'spec.co.simple.fits.S'
		thar_S = pyfits.getdata( thar_co_fits )
		lines_thar  = thar_S.copy()

		All_Pixel_Centers_co = np.array([])
		All_Wavelengths_co   = np.array([])
		All_Orders_co        = np.array([])
		All_Centroids_co     = np.array([])
		All_Sigmas_co        = np.array([])
		All_Intensities_co   = np.array([])
		All_residuals_co     = np.array([])

		orders_offset_co, rough_shift_co = fiesutils.get_thar_offsets(lines_thar,binning=binning)
		print 'orders_ofset_co:',orders_offset_co
		print 'rough_shift_co:',rough_shift_co

		orderi = 0
		if orders_offset_co < 0:
			orderi = - orders_offset_co
		orderf = nord_co - 1
		if orderf + orders_offset_co >= n_useful:
			orderf = n_useful - orders_offset_co - 1

		for order in range(orderi,orderf+1):
			#print order
			order_s = str(order+orders_offset_co)
			if (order + orders_offset_co < 10):
				order_s = '0' + str(order+orders_offset_co)
			f = open(order_dir+'order_'+order_s+'.iwdat','r')
			llins = f.readlines()
			if len(llins)>5:
					thar_order_orig = lines_thar[order]
					#IV              = iv_thar_ob_R[order,:]
					L = np.where(thar_order_orig != 0)[0]
					IV = 1. / (thar_order_orig / gain + (ron/gain)**2 )
					IV[L] = 0.
					wei             = np.ones(len(thar_order_orig))  #np.sqrt( IV )
					bkg = scipy.signal.medfilt(thar_order_orig,101)
					thar_order      = thar_order_orig - bkg

					coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids,sigmas, intensities \
			        = GLOBALutils.Initial_Wav_Calibration( order_dir+'order_'+order_s+'.iwdat', thar_order, order, wei, \
			        	rmsmax=500, minlines=10, FixEnds=False,Dump_Argon=dumpargon, Dump_AllLines=True, line_width=6, Cheby=use_cheby,porder=3,rough_shift=rough_shift,binning=binning,del_width=5.,do_xc=False)
					#plot(wavelengths, residuals,'.')	
					#show()		            
					if (order == int(.5*n_useful)): 	
						if (use_cheby):
							Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav,  0.5*len(thar_order), len(thar_order) )
						else:
							Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

					All_Pixel_Centers_co = np.append( All_Pixel_Centers_co, pixel_centers )
					All_Wavelengths_co   = np.append( All_Wavelengths_co, wavelengths )
					All_Orders_co        = np.append( All_Orders_co, np.zeros( len(pixel_centers) ) + order)
					All_Centroids_co     = np.append( All_Centroids_co, centroids)
					All_Sigmas_co        = np.append( All_Sigmas_co, sigmas)
					All_Intensities_co   = np.append( All_Intensities_co, intensities )
					All_residuals_co     = np.append( All_residuals_co, residuals )

		p0    = np.zeros( npar_wsol )
		p0[0] =  int(.5*n_useful) * Global_ZP 
		p1_co, G_pix_co, G_ord_co, G_wav_co, II_co, rms_ms_co, G_res_co = \
			GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co, All_Wavelengths_co,\
				All_Orders_co, np.ones(All_Intensities_co.shape), p0, Cheby=use_cheby,\
				order0=ro0+orders_offset_co, ntotal=n_useful, maxrms=MRMS, Inv=Inverse_m,\
				minlines=300, npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

		pdict = {'orders_offset':orders_offset,'orders_offset_co':orders_offset_co, 'rough_shift':rough_shift,\
			 'rough_shift_co':rough_shift_co,'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                         'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Wavelengths':All_Wavelengths, \
			 'All_Orders':All_Orders, 'All_Pixel_Centers':All_Pixel_Centers, 'All_Sigmas':All_Sigmas,\
			 'p1_co':p1_co, 'G_pix_co':G_pix_co, 'G_ord_co':G_ord_co, 'G_wav_co':G_wav_co, 'II_co':II_co,\
			 'rms_ms_co':rms_ms_co, 'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Wavelengths':All_Wavelengths,\
			 'All_Orders_co':All_Orders_co, 'All_Pixel_Centers_co':All_Pixel_Centers_co, 'All_Sigmas_co':All_Sigmas_co}
		pickle.dump( pdict, open(wavsol_pkl, 'w' ) )

	else:
		print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl
		pdict = pickle.load(open(wavsol_pkl,'r'))

if len(ThAr_sim)>0:
	wavsol_pkl = dirout + ThAr_sim[0].split('/')[-1][:-4]+'wavsolpars.pkl'
	dct        = pickle.load(open(wavsol_pkl,'r'))
	p_ref = dct['p1']
	p_ref_co = dct['p1_co']
	orders_offset = dct['orders_offset']
	orders_offset_co = dct['orders_offset_co']
	rough_shift = dct['rough_shift']
	rough_shift_co = dct['rough_shift_co']
else:
	if  len(ThAr_ref)>0 and  len(ThAr_co)>0:
		wavsol_pkl = dirout + ThAr_co[0].split('/')[-1][:-4]+'wavsolpars.pkl'
		dct        = pickle.load(open(wavsol_pkl,'r'))
		wavsol_ob_pkl = dirout + dct['ref_thar_ob'].split('/')[-1][:-4]+'wavsolpars.pkl'
		dct_ob     = pickle.load(open(wavsol_ob_pkl,'r'))	
		p_ref_co   = dct['p1_co']	
		p_ref      = dct_ob['p1']
		orders_offset = dct_ob['orders_offset']
		orders_offset_co = dct['orders_offset_co']
		rough_shift = dct_ob['rough_shift']
		rough_shift_co = dct['rough_shift_co']

	elif len(ThAr_ref)>0:
		wavsol_pkl = dirout + ThAr_ref[0].split('/')[-1][:-4]+'wavsolpars.pkl'
		dct        = pickle.load(open(wavsol_pkl,'r'))
		p_ref = dct['p1']
		orders_offset = dct['orders_offset']
		rough_shift = dct['rough_shift']
	
mjds_thar,shifts,shifts_co = [],[],[]
print '\n\tDetermination of instrumental drift through the night...'
for i in range(len(ThAr_sim)):
	hthar = pyfits.open( ThAr_sim[i] )
	mjd, mjd0 = fiesutils.mjd_fromheader2( hthar )
	wavsol_pkl = dirout + ThAr_sim[i].split('/')[-1][:-4]+'wavsolpars.pkl'
	dthar = pyfits.getdata( ThAr_sim[i] )
	npix = dthar.shape[0]
	dct = pickle.load(open(wavsol_pkl,'r'))
	p_shift, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(dct['G_pix'], dct['G_wav'], dct['G_ord'],\
			np.ones(len(dct['G_ord'])), p_ref, order0=ro0 + orders_offset, npix=npix,\
			Cheby=use_cheby, ntotal=n_useful, maxrms=MRMS, Inv=Inverse_m, minlines=300, nx=ncoef_x,nm=ncoef_m)
	p_shift_co, pix_centers_co, orders_co, wavelengths_co, I_co, rms_ms_co, residuals_co  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(dct['G_pix_co'], dct['G_wav_co'], dct['G_ord_co'],\
			np.ones(len(dct['G_ord_co'])), p_ref_co, order0=ro0 + orders_offset_co, npix=npix,\
			Cheby=use_cheby, ntotal=n_useful, maxrms=MRMS, Inv=Inverse_m, minlines=300, nx=ncoef_x,nm=ncoef_m)

	if rms_ms / np.sqrt(float(len(orders))) < 10. and rms_ms_co / np.sqrt(float(len(orders_co))) < 10.:
		shifts.append(p_shift[0])
		shifts_co.append(p_shift_co[0])
		mjds_thar.append(mjd)

used_thars = []
for i in range(len(ThAr_co)):
	hthar = pyfits.open( ThAr_co[i] )
	mjd, mjd0 = fiesutils.mjd_fromheader2( hthar )
	wavsol_pkl = dirout + ThAr_co[i].split('/')[-1][:-4]+'wavsolpars.pkl'
	dthar = pyfits.getdata( ThAr_co[i] )
	npix = dthar.shape[0]
	dct = pickle.load(open(wavsol_pkl,'r'))

	wavsol_ob_pkl = dirout + dct['ref_thar_ob'].split('/')[-1][:-4]+'wavsolpars.pkl'
	dct_ob     = pickle.load(open(wavsol_ob_pkl,'r'))	
	used_thars.append(dct['ref_thar_ob'])

	p_shift, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(dct_ob['G_pix'], dct_ob['G_wav'], dct_ob['G_ord'],\
			np.ones(len(dct_ob['G_ord'])), p_ref, order0=ro0 + orders_offset, npix=npix,\
			Cheby=use_cheby, ntotal=n_useful, maxrms=MRMS, Inv=Inverse_m, minlines=300, nx=ncoef_x,nm=ncoef_m)
	p_shift_co, pix_centers_co, orders_co, wavelengths_co, I_co, rms_ms_co, residuals_co  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(dct['G_pix_co'], dct['G_wav_co'], dct['G_ord_co'],\
			np.ones(len(dct['G_ord_co'])), p_ref_co, order0=ro0 + orders_offset_co, npix=npix,\
			Cheby=use_cheby, ntotal=n_useful, maxrms=MRMS, Inv=Inverse_m, minlines=300, nx=ncoef_x,nm=ncoef_m)

	if rms_ms / np.sqrt(float(len(orders))) < 10. and rms_ms_co / np.sqrt(float(len(orders_co))) < 10.:
		shifts.append(p_shift[0])
		shifts_co.append(p_shift_co[0])
		mjds_thar.append(mjd)
used_thars = np.array(used_thars)

for i in range(len(ThAr_ref_dates)):
	if not ThAr_ref[i] in used_thars:
		hthar = pyfits.open( ThAr_ref[i] )
		mjd, mjd0 = fiesutils.mjd_fromheader2( hthar )
		wavsol_pkl = dirout + ThAr_ref[i].split('/')[-1][:-4]+'wavsolpars.pkl'
		dthar = pyfits.getdata( ThAr_ref[i] )
		npix = dthar.shape[0]
		dct = pickle.load(open(wavsol_pkl,'r'))
		p_shift, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
			    GLOBALutils.Global_Wav_Solution_vel_shift(dct['G_pix'], dct['G_wav'], dct['G_ord'],\
				                                np.ones(len(dct['G_ord'])), p_ref, order0=ro0 + orders_offset, npix=npix,\
				                                Cheby=use_cheby, ntotal=n_useful, maxrms=MRMS, Inv=Inverse_m, minlines=300, nx=ncoef_x,nm=ncoef_m)

		if rms_ms / np.sqrt(float(len(orders))) < 10.:
			shifts.append(p_shift[0])
			shifts_co.append(p_shift[0])
			mjds_thar.append(mjd)

mjds_thar,shifts = np.array(mjds_thar),np.array(shifts)
I = np.argsort(mjds_thar)
mjds_thar  = mjds_thar[I]
shifts     = shifts[I]
shv = (1e-6*shifts)*299792458.0

if len(shifts)>1:
	tck_v = scipy.interpolate.splrep(mjds_thar,shv,k=1)
	tck_shift = scipy.interpolate.splrep(mjds_thar,shifts,k=1)

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

new_list = []
new_times = []
for fsim in sim_sci:
	h = pyfits.open(fsim)
	if object2do in h[0].header['OBJECT'] or object2do == 'all':
		new_list.append(fsim)
		mjd, mjd0 = fiesutils.mjd_fromheader2( h )
		new_times.append(mjd)

new_list  = np.array(new_list)
new_times = np.array(new_times)
IS        = np.argsort(new_times)
new_list  = new_list[IS]
new_times = new_times[IS]

for fsim in new_list:

	know_moon = False
	if fsim.split('/')[-1] in spec_moon:
		I = np.where(fsim.split('/')[-1] == spec_moon)[0]
		know_moon = True
		here_moon = use_moon[I]

	print '\n'
	print "\t--> Working on image: ", fsim
	h             = pyfits.open(fsim)
	mjd,mjd0      = fiesutils.mjd_fromheader2(h)
	ronoise, gain = fiesutils.get_RONGAIN(h[0].header)

	# Object name
	obname    = h[0].header['OBJECT']
	print "\t\tObject name:",obname

	# Open file, trim, overscan subtract and MasterBias subtract
	data        = h[1].data
	data        = fiesutils.OverscanTrim( data, binning=binning ) - MasterBias
	if dark_corr and len(darks)>0 and int(h[0].header['EXPTIME']) in dark_utimes.astype('int'):
		I = np.where(dark_utimes.astype('int') == int(h[0].header['EXPTIME']))[0]
		data = data - pyfits.getdata(dark_names[I][0])

	#drift,c_new = GLOBALutils.get_drift(data,P,c_all,pii=1024,win=10)
	#P_new       = GLOBALutils.shift_P(P,drift,c_new,ext_aperture)
	#print 'ydrift:',drift

	simult = False
	if h[0].header['FILMP4'] == 1:
		simult = True

	bacfile = dirout + 'BAC_' + fsim.split('/')[-1][:-4]+'fits'
	if os.access(bacfile,os.F_OK)==False:
		if simult:
			Centers = np.zeros((len(c_tot),dthar.shape[0]))
			ccc = c_tot.copy()
		else:
			Centers = np.zeros((len(c_all),dthar.shape[0]))
			ccc = c_all.copy()

		for i in range(Centers.shape[0]):
			Centers[i,:]=scipy.polyval(ccc[i,:],np.arange(len(Centers[i,:])))

		if simult:
			bac = GLOBALutils.get_scat(data.T,Centers,span=ext_aperture,option=1)
		else:
			bac = GLOBALutils.get_scat(data.T,Centers,span=bacap)

		if (os.access(bacfile,os.F_OK)):
			os.remove( bacfile )
		hdbac = pyfits.PrimaryHDU( bac )
		hdbac.writeto(bacfile)
	else:
		bac = pyfits.getdata(bacfile)
	data = data.T - bac

	ra          = float(h[0].header['RA'])
	dec         = float(h[0].header['DEC'])
	altitude    = 2382.
	latitude    = 28.75722
	longitude   = -17.885
	epoch       = 2000.

	ra2,dec2 = GLOBALutils.getcoords(obname,mjd,filen=reffile)
	if ra2 !=0 and dec2 != 0:
		ra = ra2
		dec = dec2
	else:
		print '\t\tUsing the coordinates found in the image header.'

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

	gobs      = ephem.Observer()  
	gobs.name = 'La Palma'
	gobs.lat  = rad(latitude)  # lat/long in decimal degrees  
	gobs.long = rad(longitude) 
	gobs.date = h[0].header['DATE-OBS'][:10] + ' ' + h[0].header['DATE-OBS'][11:]
	mephem    = ephem.Moon()
	mephem.compute(gobs)

	Mcoo = jplephem.object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
	Mp = jplephem.barycentric_object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
	Sp = jplephem.barycentric_object_track("Sun", int(mjd), float(mjd%1), 1, 0.0)
	res  = jplephem.object_doppler("Moon", int(mjd), mjd%1, 1, 0.0)
	lunation,moon_state,moonsep,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,ra,dec)
	refvel = bcvel_baryc + moonvel
	print '\t\tRadial Velocity of sacttered moonlight:',refvel
	#moon_alts.update({fsim:mephem.alt})
	#moon_ills.update({fsim:lunation})

	print '\t\tExtraction:'
	sci_fits        = dirout + fsim.split('/')[-1][:-4]+'spec.fits.S'
	sci_fits_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.fits.S'
	if ( os.access(sci_fits,os.F_OK) == False ) or ( os.access(sci_fits_simple,os.F_OK) == False ) or \
		( force_sci_extract ):

		print "\t\t\tNo previous extraction or extraction forced for science file", fsim, "extracting..."
		sci_Ss = GLOBALutils.simple_extraction(data,c_all,ext_aperture,\
                                                  min_extract_col,max_extract_col,npools)
		sci_S  = GLOBALutils.optimal_extraction(data,P,c_all,ext_aperture,\
                                                   ronoise,gain,S_Marsh,NCosmic_Marsh,\
                                                   min_extract_col,max_extract_col,npools)

		sci_Ss = sci_Ss[::-1]
		sci_Ss = GLOBALutils.invert(sci_Ss)

		for iii in range(3):
			sci_St =  sci_S[:,iii].copy()
			sci_St = sci_St[::-1]
			sci_S[:,iii] = sci_St
		sci_S = GLOBALutils.invert(sci_S)
            
		if (os.access(sci_fits,os.F_OK)):
			os.remove( sci_fits )
		if (os.access(sci_fits_simple,os.F_OK)):
			os.remove( sci_fits_simple )

		hdu = pyfits.PrimaryHDU( sci_S )
		hdu.writeto( sci_fits )
		hdu = pyfits.PrimaryHDU( sci_Ss )
		hdu.writeto( sci_fits_simple )
	

	else:
		print '\t\t\t '+fsim+" has already been extracted, reading in product fits files..."
		sci_S  = pyfits.getdata( sci_fits )
		sci_Ss = pyfits.getdata( sci_fits_simple )

	if simult:
		sci_co_fits        = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
		sci_co_simple_fits = dirout + fsim.split('/')[-1][:-4]+'spec.co.simple.fits.S'
		if ( os.access(sci_co_fits,os.F_OK) == False ) or ( os.access(sci_co_simple_fits,os.F_OK) == False ) or \
			( force_sci_extract ):

			print "\t\t\tNo previous extraction or extraction forced for comparison orders of science file", fsim, "extracting..."
			sci_co_Ss = GLOBALutils.simple_extraction(data,c_co,ext_aperture,\
		                                          min_extract_col,max_extract_col,npools)
			sci_co_S  = GLOBALutils.optimal_extraction(data,P_co,c_co,ext_aperture,\
		                                           ronoise,gain,S_Marsh,NCosmic_Marsh,\
		                                           min_extract_col,max_extract_col,npools)

			sci_co_Ss = sci_co_Ss[::-1]
			sci_co_Ss = GLOBALutils.invert(sci_co_Ss)

			for iii in range(3):
				sci_co_St =  sci_co_S[:,iii].copy()
				sci_co_St = sci_co_St[::-1]
				sci_co_S[:,iii] = sci_co_St
			sci_co_S = GLOBALutils.invert(sci_co_S)
		    
			if (os.access(sci_co_fits,os.F_OK)):
				os.remove( sci_co_fits )
			if (os.access(sci_co_simple_fits,os.F_OK)):
				os.remove( sci_co_simple_fits )

			hdu = pyfits.PrimaryHDU( sci_co_S )
			hdu.writeto( sci_co_fits )
			hdu = pyfits.PrimaryHDU( sci_co_Ss )
			hdu.writeto( sci_co_simple_fits )
	

		else:
			print '\t\t\t '+fsim+" has already been extracted, reading in product fits files..."
			sci_co_S  = pyfits.getdata( sci_co_fits )
			sci_co_Ss = pyfits.getdata( sci_co_simple_fits )


	fout = 'proc/' + obname + '_' + h[0].header['DATE-OBS'] + '_sp.fits'
	if ( os.access(dirout+fout ,os.F_OK) == False ) or (force_spectral_file_build):
		orderi = 0
		if orders_offset < 0:
			orderi = - orders_offset
		orderf = nord - 1
		if orderf + orders_offset >= n_useful:
			orderf = n_useful - orders_offset - 1

		spec = np.zeros((11, orderf - orderi + 1, data.shape[1]))
		hdu  = pyfits.PrimaryHDU( spec )
		hdu = GLOBALutils.update_header(hdu,'HIERARCH MJD', mjd)
		hdu = GLOBALutils.update_header(hdu,'HIERARCH MBJD', mbjd)
		hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START DATE', h[0].header['DATE-OBS'][:10] )
		hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START UT',  h[0].header['DATE-OBS'][11:])
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
		hdu = GLOBALutils.update_header(hdu,'HIERARCH TARG AIRMASS START',h[0].header['AIRMASS'])
		hdu = GLOBALutils.update_header(hdu,'HIERARCH MOON VEL',refvel)
		hdu = GLOBALutils.update_header(hdu,'HIERARCH SIMULT',simult)

	print '\t\tWavelength calibration:'

	if simult and (len(ThAr_co)>0 or len(ThAr_sim)>0):
		lines_thar  = sci_co_S[:,1,:].copy()

		All_Pixel_Centers_co = np.array([])
		All_Wavelengths_co   = np.array([])
		All_Orders_co        = np.array([])
		All_Centroids_co     = np.array([])
		All_Sigmas_co        = np.array([])
		All_Intensities_co   = np.array([])
		All_residuals_co     = np.array([])

		orderi = 0
		if orders_offset_co < 0:
			orderi = - orders_offset_co
		orderf = nord_co - 1
		if orderf + orders_offset_co >= n_useful:
			orderf = n_useful - orders_offset_co - 1

		for order in range(orderi,orderf+1):
			#print order
			order_s = str(order+orders_offset_co)
			if (order + orders_offset_co < 10):
				order_s = '0' + str(order+orders_offset_co)
			f = open(order_dir+'order_'+order_s+'.iwdat','r')
			llins = f.readlines()
			if len(llins)>5:
					thar_order_orig = lines_thar[order]
					#IV              = iv_thar_ob_R[order,:]
					L = np.where(thar_order_orig != 0)[0]
					IV = 1. / (thar_order_orig / gain + (ron/gain)**2 )
					IV[L] = 0.
					wei             = np.ones(len(thar_order_orig))  #np.sqrt( IV )
					bkg = scipy.signal.medfilt(thar_order_orig,101)
					thar_order      = thar_order_orig - bkg

					coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids,sigmas, intensities \
			        = GLOBALutils.Initial_Wav_Calibration( order_dir+'order_'+order_s+'.iwdat', thar_order, order, wei, \
			        	rmsmax=500, minlines=10, FixEnds=False,Dump_Argon=dumpargon, Dump_AllLines=True, line_width=6, Cheby=use_cheby,porder=3,rough_shift=rough_shift,binning=binning,del_width=5.,do_xc=False)
					#plot(wavelengths, residuals,'.')	
					#show()		            
					if (order == int(.5*n_useful)): 	
						if (use_cheby):
							Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav,  0.5*len(thar_order), len(thar_order) )
						else:
							Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

					All_Pixel_Centers_co = np.append( All_Pixel_Centers_co, pixel_centers )
					All_Wavelengths_co   = np.append( All_Wavelengths_co, wavelengths )
					All_Orders_co        = np.append( All_Orders_co, np.zeros( len(pixel_centers) ) + order)
					All_Centroids_co     = np.append( All_Centroids_co, centroids)
					All_Sigmas_co        = np.append( All_Sigmas_co, sigmas)
					All_Intensities_co   = np.append( All_Intensities_co, intensities )
					All_residuals_co     = np.append( All_residuals_co, residuals )

		p0    = np.zeros( npar_wsol )
		p0[0] =  int(.5*n_useful) * Global_ZP 
		p1_co, G_pix_co, G_ord_co, G_wav_co, II_co, rms_ms_co, G_res_co = \
		GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co, All_Wavelengths_co,\
				All_Orders_co, np.ones(All_Intensities_co.shape), p0, Cheby=use_cheby,\
				order0=ro0+orders_offset_co, ntotal=n_useful, maxrms=MRMS, Inv=Inverse_m,\
				minlines=300, npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

		p_shift, pix_centers_co, orders_co, wavelengths_co, I_co, rms_ms_co, residuals_co  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(G_pix_co, G_wav_co, G_ord_co,\
			np.ones(len(G_ord_co)), p_ref_co, order0=ro0 + orders_offset_co, npix=len(thar_order),\
			Cheby=use_cheby, ntotal=n_useful, maxrms=MRMS, Inv=Inverse_m, minlines=300, nx=ncoef_x,nm=ncoef_m)

	else:
		if len(shifts)>1:
			p_shift = scipy.interpolate.splev(mjd,tck_shift)
		else:
			p_shift = 0.

	orderi = 0
	if orders_offset < 0:
		orderi = - orders_offset
	orderf = nord - 1
	if orderf + orders_offset >= n_useful:
		orderf = n_useful - orders_offset - 1

	print "\t\t\tInstrumental drift:",(1e-6*p_shift)*299792458.0, 'm/s'
	# Apply new wavelength solution including barycentric correction
	equis = np.arange( data.shape[1] )
	order = orderi
	while order < orderf+1:
		m   = order + ro0 + orders_offset
		chebs  = GLOBALutils.Calculate_chebs(equis, m, npix=sci_S.shape[2], order0=ro0 + orders_offset, ntotal=n_useful, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
		WavSol = GLOBALutils.ToVacuum( lbary_ltopo * (1.0 + 1.0e-6*p_shift) * (1.0/float(m)) * \
		     GLOBALutils.Joint_Polynomial_Cheby(p_ref,chebs,ncoef_x,ncoef_m) )

		spec[0,order,:] = WavSol
		spec[1,order,:] = sci_S[order,1]
		spec[2,order,:] = sci_S[order,2]
		fn  = S_flat[order,1,:]
		L  = np.where( fn == 0 )[0]
		spec[3,order,:] = spec[1,order,:] / S_flat[order,1,:]
		spec[4,order,:] = spec[2,order] * ( S_flat_n[order,1,:] ** 2 )
		spec[3,order,L] = 0.
		spec[4,order,L] = 0.
		ccoef = GLOBALutils.get_cont_single(spec[0,order],spec[3,order],spec[4,order],ll=1.5,lu=5,nc=3)
		L  = np.where( spec[1,order] != 0 )
		spec[5,order,:][L] = spec[3,order][L] / np.polyval(ccoef,spec[0,order][L])    
		ratio            = np.polyval(ccoef,spec[0,order][L]) * Snorms[order]
		spec[6,order,:][L] = spec[4,order][L] * (ratio ** 2 )
		spec[7,order,:][L] = ratio
		spec[8,order,:][L] = ratio * S_flat_n[order,1][L] / np.sqrt( ratio * S_flat_n[order,1][L] / gain + (ronoise/gain)**2 )
		spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
		dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
		NN            = np.average(dlambda_dx)
		dlambda_dx    /= NN
		LL = np.where(spec[5,order] > 1 + 10. / scipy.signal.medfilt(spec[8,order],21))[0]
		spec[5,order,LL] = 1.
		spec[9,order][L] = spec[5,order][L] * (dlambda_dx[L] ** 1) 
		spec[10,order][L] = spec[6,order][L] / (dlambda_dx[L] ** 2)
		#plot(spec[0,order],spec[5,order])
		order +=1
	#show()
	if os.access(dirout + fout, os.F_OK):
		os.remove(dirout + fout)
	hdu.writeto(dirout + fout)
	
	if (not JustExtract):

		if DoClass:
			print '\t\tSpectral Analysis:'
			# spectral analysis
			# First, query SIMBAD with the object name
			query_success = False
			sp_type_query = 'None'
			#query_success,sp_type_query = GLOBALutils.simbad_query_obname(obname)
			# Now, query SIMBAD by coordinates if above not successful
			#if (not query_success):
			#	query_success,sp_type_query = GLOBALutils.simbad_query_coords('12:00:00','00:00:00')
			print "\t\t\tSpectral type returned by SIMBAD query:",sp_type_query

			hdu = GLOBALutils.update_header(hdu,'HIERARCH SIMBAD SPTYP', sp_type_query)

			pars_file = dirout + fsim.split('/')[-1][:-8]+'_stellar_pars.txt'

			if os.access(pars_file,os.F_OK) == False or force_stellar_pars:
				print "\t\t\tEstimating atmospheric parameters:"
				spec2 = spec.copy()
				#print resol
				if resol > 44000:
					Rx = np.around(1./np.sqrt(1./40000.**2 - 1./resol**2))
					for i in range(spec.shape[1]):
						IJ = np.where(spec[5,i]!=0.)[0]
						spec2[5,i,IJ] = GLOBALutils.convolve(spec[0,i,IJ],spec[5,i,IJ],Rx)
						#plot(spec[0,i],spec[5,i])
						#plot(spec2[0,i],spec2[5,i])
					#show()
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

		# Read in mask
		ml, mh, weight = np.loadtxt(mask,unpack=True)
		ml_v = GLOBALutils.ToVacuum( ml )
		mh_v = GLOBALutils.ToVacuum( mh )
       
		# make mask larger accounting for factor ~2 lower res in CORALIE w/r to HARPS
		av_m = 0.5*( ml_v + mh_v )
		ml_v -= (av_m - ml_v)
		mh_v += (mh_v - av_m)
		mask_hw_kms = (GLOBALutils.Constants.c/1e3) * 0.5*(mh_v - ml_v) / av_m

		#sigma_fout = stellar_pars_dir + obname + '_' +'sigma.txt'

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

			# Normalize the continuum of the CCF robustly with lowess     
			yy = scipy.signal.medfilt(xc_av,11)
			pred = lowess(yy, vels,frac=0.4,it=10,return_sorted=False)
			tck1 = scipy.interpolate.splrep(vels,pred,k=1)
			xc_av_orig = xc_av.copy()
			xc_av /= pred
			vel0_xc = vels[ np.argmin( xc_av ) ]
			rvels, rxc_av, rpred, rxc_av_orig, rvel0_xc = vels.copy(), \
			xc_av.copy(), pred.copy(), xc_av_orig.copy(), vel0_xc

			xc_av_rough = xc_av
			vels_rough  = vels
			if disp > 30:
				disp = 30.
			vel_width = np.maximum( 20.0, 6*disp )

			vels, xc_full, sn, nlines_ccf, W_ccf =\
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight, vel0_xc, lbary_ltopo, vel_width=vel_width,vel_step=0.1,\
                                          spec_order=9,iv_order=10,sn_order=8,max_vel_rough=300)

			xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3.0, Simple=True, W=W_ccf)
			pred = scipy.interpolate.splev(vels,tck1)
			xc_av /= pred
		
			if sp_type == 'M5':
				moon_sig = 2.5
			elif sp_type == 'K5':
				moon_sig = 3.3
			else:
				moon_sig = 4.5

			p1,XCmodel,p1gau,XCmodelgau,Ls2 = GLOBALutils.XC_Final_Fit( vels, xc_av ,\
                                  sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = False)

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
		         'lunation':lunation,'mephem':mephem,'texp':float(h[0].header['EXPTIME'])}

		pkl_xc = dirout + fsim.split('/')[-1][:-8]+obname+'_XC_'+sp_type+'.pkl'
		pickle.dump( xc_dict, open( pkl_xc, 'w' ) )

		ccf_pdf = dirout + 'proc/' + fsim.split('/')[-1][:-4] + obname + '_XCs_' + sp_type + '.pdf'

		if not avoid_plot:
			GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)

		SNR_5130 = np.median(spec[8,32,900:1101] )
		airmass  = float(h[0].header['AIRMASS'])
		seeing   = -999
		TEXP = float(h[0].header['EXPTIME'])

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
		hdu = GLOBALutils.update_header(hdu,'RV', RV)
		hdu = GLOBALutils.update_header(hdu,'RV_E', RVerr2)
		hdu = GLOBALutils.update_header(hdu,'BS', BS)
		hdu = GLOBALutils.update_header(hdu,'BS_E', BSerr)
		hdu = GLOBALutils.update_header(hdu,'DISP', disp_epoch)
		hdu = GLOBALutils.update_header(hdu,'SNR', SNR_5130)
		hdu = GLOBALutils.update_header(hdu,'SNR_R', SNR_5130_R)
		hdu = GLOBALutils.update_header(hdu,'INST', 'FIES')
		hdu = GLOBALutils.update_header(hdu,'RESOL', str(resol))
		hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')
		hdu = GLOBALutils.update_header(hdu,'XC_MIN', XC_min)
		hdu = GLOBALutils.update_header(hdu,'BJD_OUT', bjd_out)

		# write to output
		line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f   fies   ceres   %8d %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
                      (obname, bjd_out, RV, RVerr2, BS, BSerr, resol, T_eff_epoch, logg_epoch, Z_epoch, vsini_epoch, XC_min, disp_epoch,\
			TEXP, SNR_5130_R, ccf_pdf)
		f_res.write(line_out)
    
		if (os.access( dirout + fout,os.F_OK)):
			os.remove( dirout + fout)

		hdu.writeto( dirout + fout )
	else:
		print "Reading spectral file from", fout
		spec = pyfits.getdata( fout )
    
f_res.close()
