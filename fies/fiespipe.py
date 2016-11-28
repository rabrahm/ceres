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

import rpy2.robjects.numpy2ri
try:
	rpy2.robjects.numpy2ri.activate()
except:
	None
from rpy2 import robjects
robjects.r.library("MASS")
"""
import rpy2
from rpy2 import robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
rpy2.robjects.activate()
r = robjects.r
r.library("MASS")
"""
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

args = parser.parse_args()
DoClass     = args.do_class
avoid_plot  = args.avoid_plot
dirin       = args.directorio
object2do   = args.o2do
JustExtract = args.just_extract
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
force_thar_wavcal  = False
force_spectral_file_build = True
force_stellar_pars = False
dumpargon          = False
dark_corr		   = True
minlines_glob      = 700

Inverse_m          = True
use_cheby          = True
MRMS               = 100   # max rms in m/s, global wav solution

trace_degree       = 4
Marsh_alg          = 0
ext_aperture       = 5
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

n_useful = 79
ro0      = 79

bacap=7
orders_offset = 0
if binning == 2:
	ext_aperture = 4
	min_extract_col = int(np.around(min_extract_col/2.))
	max_extract_col = int(np.around(max_extract_col/2.))
	bacap = 4
	orders_offset = 2

if mode == 'F1':
	resol = 25000
	MRMS  = 200.
elif mode == 'F3':
	resol = 46000
elif mode == 'F4':
	resol = 67000
#############################

# file containing the log
log = dirout+'night.log'

print "\n\n\tFIES NOT2.5m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'

biases, flats, ThAr_ref, sim_sci, ThAr_ref_dates, obnames, exptimes, darks = fiesutils.FileClassify(dirin,log,binning=binning, mode=mode, dark_corr=dark_corr)
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
	MasterBias, RO_bias, GA_bias = fiesutils.MedianCombine(biases, binning=binning)
	hdu = pyfits.PrimaryHDU( MasterBias )
	if (os.access(dirout+'MasterBias.fits',os.F_OK)):
		os.remove(dirout+'MasterBias.fits')
	hdu.writeto(dirout+'MasterBias.fits')
	print "\t\t-> Masterbias: done!"

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
	print "\t\t-> Masterflats: done!"

	print "\tTracing echelle orders..."
	c_all,nord = GLOBALutils.get_them(Flat.T, ext_aperture, trace_degree, maxords=-1, mode=1,nsigmas=5.)
	c_all,nord = GLOBALutils.good_orders(c_all,nord,Flat.shape[1],Flat.shape[0],ext_aperture)
	print '\t\t'+ str(nord)+' orders found.'

	# pickle traces
	trace_dict = {'c_all':c_all,
                  'nord':nord,
                  'GA_bias': GA_bias, 'RO_bias' : RO_bias,
                  'GA_fl': GA_fl, 'RO_fl': RO_fl,
                  'dark_names':dark_names, 'dark_utimes':dark_utimes}
	pickle.dump( trace_dict, open( dirout+"trace.pkl", 'w' ) )

else:
	trace_dict = pickle.load( open( dirout+"trace.pkl", 'r' ) )
	c_all = trace_dict['c_all']
	nord  = trace_dict['nord']
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
	h = pyfits.open(dirout+'MasterBias.fits')
	MasterBias = h[0].data

print '\n\tExtraction of Flat calibration frames:'

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

# Normalize flat field spectra.
S_flat_n, Snorms = GLOBALutils.FlatNormalize_single( S_flat, mid = int(.5*S_flat.shape[2]) )

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
		bac = GLOBALutils.get_scat(dthar.T,Centers,span=bacap)
		sdthar = dthar.T - bac
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

print "\n\tWavelength solution of ThAr calibration spectra:"

for i in range(len(ThAr_ref_dates)):
	rough_shift = 0.
	wavsol_pkl = dirout + ThAr_ref[i].split('/')[-1][:-4]+'wavsolpars.pkl'
	h = pyfits.open(ThAr_ref[i])
	hd = pyfits.getheader(ThAr_ref[i])
	ron,gain = fiesutils.get_RONGAIN(h[0].header)
	#force_thar_wavcal = True
	if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
		print "\t\tWorking on initial ThAr file", ThAr_ref[i] 
		hthar = pyfits.open( ThAr_ref[i] )
		mjd, mjd0 = fiesutils.mjd_fromheader2( hthar )
		thar_fits = dirout + ThAr_ref[i].split('/')[-1][:-4]+'spec.simple.fits.S'
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
		for order in range(nord):
			order_s = str(order+orders_offset)
			if (order + orders_offset < 10):
				order_s = '0' + str(order+orders_offset)
			f = open(order_dir+'order_'+order_s+'.iwdat','r')
			llins = f.readlines()
			if len(llins)>5:
				thar_order = lines_thar[order]
				pixel_centers_0 = []
				for line in llins:
					w = line.split()
					nlines = int(w[0])
					for j in range(nlines):
						pixel_centers_0.append(float(w[2*j+1])*1./float(binning))
				pixel_centers_0 = np.array(pixel_centers_0).astype('int')
				ml = np.array(pixel_centers_0) - 2
				mh = np.array(pixel_centers_0) + 2

				xc,offs = GLOBALutils.XCorPix( thar_order, ml, mh, del_width=10.)
				if len(xcs)==0:
					xcs = xc.copy()
				else:
					xcs = xcs + xc
		#plot(offs,xcs)
		#show()
		rough_shift = offs[np.argmax(xcs)]

		for order in range(lines_thar.shape[0]):
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
				thar_order      = thar_order_orig #- bkg

				coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids,sigmas, intensities \
	                = GLOBALutils.Initial_Wav_Calibration( order_dir+'order_'+order_s+'.iwdat', thar_order, order, wei, \
	                	rmsmax=1000, minlines=10, FixEnds=False,Dump_Argon=dumpargon, Dump_AllLines=True, line_width=6, Cheby=use_cheby,porder=3,rough_shift=rough_shift,binning=binning,del_width=5.,do_xc=False)
				#plot(wavelengths, residuals,'.')	
				#show()		            
				if (order == int(.5*nord)): 	
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
		p0[0] =  int(.5*nord) * Global_ZP 
		p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths,\
						     All_Orders, np.ones(All_Intensities.shape), p0, Cheby=use_cheby,       \
						     order0=ro0+orders_offset, ntotal=nord, maxrms=MRMS, Inv=Inverse_m, minlines=300,  \
						     npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

		#for i in np.unique(G_ord):
		#	I = np.where(G_ord == i)[0]
		#	plot(G_wav[I],G_res[I],'.')
		#show()

		pdict = {'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                     'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Wavelengths':All_Wavelengths, 'All_Orders':All_Orders, 'All_Pixel_Centers':All_Pixel_Centers, 'All_Sigmas':All_Sigmas}
		pickle.dump( pdict, open( wavsol_pkl, 'w' ) )

	else:
		print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl
		pdict = pickle.load(open(wavsol_pkl,'r'))

wavsol_pkl = dirout + ThAr_ref[0].split('/')[-1][:-4]+'wavsolpars.pkl'
dct        = pickle.load(open(wavsol_pkl,'r'))
p_ref = dct['p1']
mjds_thar,shifts = [],[]
print '\n\tDetermination of instrumental drift through the night...'
for i in range(len(ThAr_ref_dates)):
	hthar = pyfits.open( ThAr_ref[i] )
	mjd, mjd0 = fiesutils.mjd_fromheader2( hthar )
	wavsol_pkl = dirout + ThAr_ref[i].split('/')[-1][:-4]+'wavsolpars.pkl'
	dthar = pyfits.getdata( ThAr_ref[i] )
	npix = dthar.shape[0]
	mjds_thar.append(mjd)
	dct = pickle.load(open(wavsol_pkl,'r'))
	p_shift, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(dct['G_pix'], dct['G_wav'], dct['G_ord'],\
		                                        np.ones(len(dct['G_ord'])), p_ref, order0=ro0 + orders_offset, npix=npix,\
		                                        Cheby=use_cheby, ntotal=nord, maxrms=MRMS, Inv=Inverse_m, minlines=300, nx=ncoef_x,nm=ncoef_m)
	shifts.append(p_shift[0])

mjds_thar,shifts = np.array(mjds_thar),np.array(shifts)
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

	bacfile = dirout + 'BAC_' + fsim.split('/')[-1][:-4]+'fits'
	if os.access(bacfile,os.F_OK)==False:
		Centers = np.zeros((len(c_all),data.shape[0]))
		for i in range(nord):
			Centers[i,:]=scipy.polyval(c_all[i,:],np.arange(len(Centers[i,:])))
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

	fout = 'proc/' + obname + '_' + h[0].header['DATE-OBS'] + '_sp.fits'
	if ( os.access(dirout+fout ,os.F_OK) == False ) or (force_spectral_file_build):
		spec = np.zeros((11, n_useful, data.shape[1]))
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

	print '\t\tWavelength calibration:'
	print "\t\t\tInstrumental drift:",(1e-6*p_shift)*299792458.0
	# Apply new wavelength solution including barycentric correction
	equis = np.arange( data.shape[1] )
	if len(shifts)>1:
		p_shift = scipy.interpolate.splev(mjd,tck_shift)
	else:
		p_shift = 0.
	order = 0
	while order < nord:
		m   = order + ro0
		chebs  = GLOBALutils.Calculate_chebs(equis, m, npix=data.shape[1], order0=ro0, ntotal=nord, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
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
				if resol > 44000:
					Rx = np.around(1./np.sqrt(1./40000.**2 - 1./resol**2))
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
			vel_width = np.maximum( 20.0, 6*disp )

			vels, xc_full, sn, nlines_ccf, W_ccf =\
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight, vel0_xc, lbary_ltopo, vel_width=vel_width,vel_step=0.1,\
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

