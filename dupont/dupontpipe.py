import sys
from pylab import *
ioff()

base = '../'

sys.path.append(base+"utils/Continuum/")
sys.path.append(base+"utils/Correlation/")
sys.path.append(base+"utils/GLOBALutils/")
sys.path.append(base+"utils/OptExtract/")

baryc_dir= base+'utils/SSEphem/'
sys.path.append(baryc_dir)
ephemeris='DEc403'

import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ceres modules
import dupontutils
import continuum
import correlation
import GLOBALutils
import Marsh

# other useful modules
import argparse
import ephem
import glob
import jplephem
import os
import pickle
import pyfits
import scipy
import scipy.interpolate
from scipy import optimize
from scipy import interpolate
from numpy import radians as rad

import statsmodels.api as sm
lowess = sm.nonparametric.lowess

# Recive input parameters
parser = argparse.ArgumentParser()
parser.add_argument('directorio')
parser.add_argument('-avoid_plot', action="store_true", default=False)
parser.add_argument('-dirout',default='default')
parser.add_argument('-do_class', action="store_true", default=False)
parser.add_argument('-just_extract', action="store_true", default=False)
parser.add_argument('-no_dark_sub',     action = "store_true", default = False)
parser.add_argument('-npools', default=1)
parser.add_argument('-oblaze',   default = 'same')
parser.add_argument('-ofind',    default = 'last')
parser.add_argument('-o2do',default='all')
parser.add_argument('-reffile',default='default')
parser.add_argument('-resolution',default='40000')
args = parser.parse_args()

dirin                = args.directorio
avoid_plot           = args.avoid_plot
dirout               = args.dirout
DoClass              = args.do_class
JustExtract          = args.just_extract
no_dark_substraction = args.no_dark_sub
npools               = int(args.npools)
object2do            = args.o2do
reffile              = args.reffile
stst                 = args.ofind
stblaze              = args.oblaze
resolution           = float(args.resolution)

dark_substraction = not no_dark_substraction

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

if resolution != 40000 and resolution != 50000 and resolution != 60000:
	raise ValueError("Input Resolution not permited. try 40000, 50000 or 60000\n")

####### GLOBAL VARIABLES #####
force_pre_process  = False
force_bl           = False	#
force_bkg          = False
force_P            = False
force_thar_extract = False
force_thar_wavcal  = False 
force_sci_extract  = False
force_sci_proc     = True	#
force_RV           = False	#
force_stellar_pars = False
force_corr         = False

use_ref            = False
ref_traces         = '/data/echelle/duPont/red/20120630/trace.pkl'

bad_colummn        = True
have_skyF          = False

Inverse_m          = True
use_cheby          = True

MRMS              = 200
trace_degree      = 5
Marsh_alg         = 0
ext_aperture      = 6
NSigma_Marsh      = 5
NCosmic_Marsh     = 5
S_Marsh           = 0.4
N_Marsh           = 3
min_extract_col   = 50
max_extract_col   = 2000

#npar_wsol = 27	#number of parameters of wavelength solution
ncoef_x   = 5
ncoef_m   = 6
npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2

n_last    = 55 # last echelle order as defined in the reference wavelength files
ntotal	  = 64
oro0      = 36

order_dir   = base+'dupont/wavcals/'

models_path = base+'data/COELHO_MODELS/R_40000b/'

print "\n\n\tEchelle Du Pont 2.5m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'

# file containing the log
log = dirout+'night.log'
biases, milkyflats, skyflats, objects, ThAr_ref, darks = dupontutils.FileClassify(dirin,log)
ThAr_ref = ThAr_ref[:2]
if dark_substraction == True and len(darks)<3:
    dark_substraction = False

f = open(log,'r')
lines = f.readlines()

print '\tThese are all the images to proccess:'
for bias in biases:
	hd = pyfits.getheader(bias)
	print '\tbias', hd['EXPTYPE'], hd['EXPTIME'], hd['UT-DATE'],hd['UT-TIME'],bias
print '\n'
for dark in darks:
	hd = pyfits.getheader(dark)
	print '\tdark', hd['EXPTYPE'], hd['EXPTIME'], hd['UT-DATE'],hd['UT-TIME'],dark
print '\n'
for milky in milkyflats:
	hd = pyfits.getheader(milky)
	print '\tmilky', hd['EXPTYPE'], hd['EXPTIME'], hd['UT-DATE'],hd['UT-TIME'],milky
print '\n'
for line in lines:
	print '\t'+line[:-1]

if stst == 'last':
    if os.access(dirout+'findstar.txt',os.F_OK):
        fst = open(dirout+'findstar.txt','r')
        stst = fst.readline()
        fst.close()
    else:
        raise ValueError("There is not a previously defined standard star file!!! You have to enter one (i.e. -ofind ccd0001.fits).\n")
else:
    fst = open(dirout+'findstar.txt','w')
    fst.write(stst+'\n')
    fst.close()

if stblaze == 'same':
    print '\n\tThe pipeline will use the image that traces the orders to derive the blaze function ...'
    stblaze = stst

if ( (os.access(dirout+'Flat.fits',os.F_OK) == False) or\
   (os.access(dirout+'trace.pkl',os.F_OK) == False) or \
   (os.access(dirout+'MasterBias.fits',os.F_OK) == False) or \
   (force_pre_process) ):
    print "\tNo previous pre-processing files or found"
    pre_process = 1
else:
    print "\tPre-processing files found, going straight to extraction"
    pre_process = 0

if (pre_process == 1):
    print "\n\tGenerating Master calibration frames..."
    MasterBias, RO_bias, GA_bias = dupontutils.MedianCombine(biases, zero_bo=False, dark_bo=False, flat_bo=False)
    hdu = pyfits.PrimaryHDU( MasterBias )

    if (os.access(dirout+'MasterBias.fits',os.F_OK)):
        os.remove(dirout+'MasterBias.fits')
    hdu.writeto(dirout+'MasterBias.fits')
    print "\t\t-> Masterbias: done!"

    MDARKS = []
    dark_times = []

    if dark_substraction:
		
        for dark in darks:
			hd = pyfits.getheader(dark)
			if len(dark_times) == 0:
				dark_times.append(hd['EXPTIME'])
			else:
				if dark_times.count(hd['EXPTIME']) == 0:
					dark_times.append(hd['EXPTIME'])
        dark_groups = []
        ndark_times = []	
        for time in dark_times:
            group = []
            for dark in darks:
                hd = pyfits.getheader(dark)
                if hd['EXPTIME'] == time:
                    group.append(dark)
            if len(group)>2:
                dark_groups.append(group)
                ndark_times.append(hd['EXPTIME'])

        dark_times = ndark_times
        i = 0
        while i < len(dark_times):
            DARK, RON, GAIN = dupontutils.MedianCombine(dark_groups[i], zero_bo=True, zero=dirout+'MasterBias.fits',dark_bo=False)
            hdu = pyfits.PrimaryHDU( DARK )
            hdu = GLOBALutils.update_header(hdu,'EXPTIME',dark_times[i])
            if os.access(dirout+'DARK_'+str(int(dark_times[i]))+'s.fits',os.F_OK):
                os.remove(dirout+'DARK_'+str(int(dark_times[i]))+'s.fits')
            hdu.writeto( dirout+'DARK_'+str(int(dark_times[i]))+'s.fits' )
            MDARKS.append(dirout+'DARK_'+str(int(dark_times[i]))+'s.fits')
            i+=1
    print "\t\t-> Masterdarks: done!"

    force_flat_comb = True
    if force_flat_comb or os.access(dirout+'Flat.fits',os.F_OK) == False:
        Flat, RON, GAIN = dupontutils.milk_comb(milkyflats, MDARKS, zero=dirout+'MasterBias.fits')
        hdu = pyfits.PrimaryHDU( Flat )
        hdu = GLOBALutils.update_header(hdu,'RON',RON)
        hdu = GLOBALutils.update_header(hdu,'GAIN',GAIN)
        if (os.access(dirout+'Flat_or.fits',os.F_OK)):
            os.remove(dirout+'Flat_or.fits')
        hdu.writeto(dirout+'Flat_or.fits')

        Flat = Flat/scipy.signal.medfilt(Flat,[15,15])
        hdu = pyfits.PrimaryHDU( Flat )
        hdu = GLOBALutils.update_header(hdu,'RON',RON)
        hdu = GLOBALutils.update_header(hdu,'GAIN',GAIN)

        if (os.access(dirout+'Flat.fits',os.F_OK)):
            os.remove(dirout+'Flat.fits')
        hdu.writeto(dirout+'Flat.fits')
    else:
        h = pyfits.open(dirout+'Flat.fits')
        Flat = h[0].data
        RON  = h[0].header['RON']
        GAIN = h[0].header['GAIN']
    print "\t\t-> Masterflat: done!"
	
    ftra = True
    median_filter = False

    # Find orders & traces
    print "\tTracing echelle orders..."
    if ftra or os.access(dirout+'trace.pkl',os.F_OK) == False:
        h = pyfits.open(dirin+stst)[0]
        hth = pyfits.getheader(dirin+stst)
        d = h.data
        d = dupontutils.OverscanTrim(d,hth['BIASSEC'])
        d -= MasterBias
        d /= Flat

        if bad_colummn:
            d = dupontutils.b_col(d)	

        d = d[:1420,:]	
        if use_ref:
            c_comp = pickle.load( open( ref_traces, 'r' ) )['c_all']
            nord = pickle.load( open( ref_traces, 'r' ) )['nord']
            lim = pickle.load( open( ref_traces, 'r' ) )['lims']
            c_all,pshift = GLOBALutils.retrace( d, c_comp )
        else:
            c_all, nord = GLOBALutils.get_them(d,5,trace_degree,mode=1)
            print '\t\t'+str(nord)+' orders found...'
    else:
        trace_dict = pickle.load( open( dirout+"trace.pkl", 'r' ) )
        c_all = trace_dict['c_all']
        nord = trace_dict['nord']
        GAIN = trace_dict['GA_ob']
        RON = trace_dict['RO_ob']
	
    trace_dict = {'c_all':c_all, 'nord':nord, 'GA_ob': GAIN, 'RO_ob': RON, 'DARKS':MDARKS, 'dtimes':dark_times}
    pickle.dump( trace_dict, open( dirout+"trace.pkl", 'w' ) )


else:
    trace_dict = pickle.load( open( dirout+"trace.pkl", 'r' ) )

    c_all = trace_dict['c_all']
    nord = trace_dict['nord']

    GAIN = trace_dict['GA_ob']
    RON = trace_dict['RO_ob']

    h = pyfits.open(dirout+'MasterBias.fits')
    MasterBias = h[0].data

    MDARKS = trace_dict['DARKS']
    dark_times = trace_dict['dtimes']

    h = pyfits.open(dirout+'Flat.fits')
    Flat = h[0].data

print '\n\tExtraction of ThAr calibration frames:'

for fsim in ThAr_ref:
    hth = pyfits.getheader(fsim)
    thmjd,mjd0 = dupontutils.mjd_fromheader(hth)
    dth = pyfits.getdata(fsim)
    dth = dupontutils.OverscanTrim(dth,hth['BIASSEC'])
    dth -= MasterBias
    dth /= Flat

    thar_fits_simple = dirout+'ThAr_'+hth['DATE-OBS']+'_'+hth['UT-TIME'][:2]+ \
                       '-'+hth['UT-TIME'][3:5]+'-'+hth['UT-TIME'][6:]+'.spec.simple.fits'

    if ( os.access(thar_fits_simple,os.F_OK) == False ) or (force_thar_extract): 
        print "\t\tNo previous extraction or extraction forced for ThAr file", fsim, "extracting..."
        close_spec = dupontutils.get_close(thmjd,hth['RA-D'],hth['DEC-D'],objects)
        temp_spec = pyfits.getdata(close_spec)
        temp_hd   = pyfits.getheader(close_spec)
        temp_spec = dupontutils.OverscanTrim(temp_spec,temp_hd['BIASSEC'])
        temp_spec -= MasterBias
        if bad_colummn:
            temp_spec = dupontutils.b_col(temp_spec)
        temp_spec = temp_spec/Flat
        c_alls, pshift = GLOBALutils.retrace( temp_spec, c_all )

        thar_Ss = GLOBALutils.simple_extraction(dth,c_alls,ext_aperture,min_extract_col,max_extract_col,npools)
        thar_Ss = thar_Ss[::-1]
        if (os.access(thar_fits_simple,os.F_OK)):
            os.remove( thar_fits_simple )
        hdu = pyfits.PrimaryHDU( thar_Ss )
        hdu.writeto( thar_fits_simple )
    else:
        print "\t\tThAr file", fsim, "all ready extracted, loading..."

print "\n\tWavelength solution of ThAr calibration spectra:"
# compute wavelength calibration files

thtimes = []
thnames = []
bad_thars = []
thRA = []
thDEC = []
ntr = 0
counter = 0
for thar in ThAr_ref:
    if ntr > 0:
        force_corr = False
    hth = pyfits.getheader(thar)
    thmjd,mjd0 = dupontutils.mjd_fromheader(hth)

    thar_fits_simple = dirout+'ThAr_'+hth['DATE-OBS']+'_'+hth['UT-TIME'][:2]+'-'+\
                       hth['UT-TIME'][3:5]+'-'+hth['UT-TIME'][6:]+'.spec.simple.fits'
    wavsol_pkl       = dirout+'ThAr' +hth['DATE-OBS']+'_'+hth['UT-TIME'][:2]+'-'+\
                       hth['UT-TIME'][3:5]+'-'+hth['UT-TIME'][6:]+'.wavsolpars.pkl'

    thar_Ss = pyfits.getdata(thar_fits_simple)

    if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
        print "	\t\tWorking on ThAr file", thar
	
        RON  = hth['ENOISE']
        GAIN = hth['EGAIN']
	
        lines_thar = thar_Ss[:,:]
        if os.access(dirout+'id_orders.pkl',os.F_OK) == False or force_corr:	
            maxes = 0
            or32 = 0
            for order in range(len(lines_thar)):
                ccf_max, shift = GLOBALutils.cor_thar(lines_thar[order],filename=order_dir+'order_32o.iwdat',span=50)
                if ccf_max > maxes:
                    maxes = ccf_max
                    rough_shift = shift
                    or32  =  order
            or0 = or32 - 32
			
            if or0 >= 0:
                orwa = 0
            else:
                orwa = - or0
                or0  = 0
            #print 'n_lasts',n_last,or0 + nord
            if n_last > or0 + nord:
                n_last = or0 + nord
            if or32-32 >0:
                nn = n_last + 1
            else:
                nn = n_last + 1 + or32 -32
            #print 'NN',nn

            if os.access(dirout+'id_orders.pkl',os.F_OK):
                os.remove(dirout+'id_orders.pkl')
			
            pdict = {'or0':or0, 'orwa':orwa, 'n_last':n_last, 'rough_shift':rough_shift,'nn':nn}
            pickle.dump( pdict, open(dirout+'id_orders.pkl', 'w' ) )
            ntr += 1
		
        else:
            pdict       = pickle.load(open(dirout+'id_orders.pkl', 'r'))
            or0         = pdict['or0']
            orwa        = pdict['orwa']
            n_last      = pdict['n_last']
            nn          = pdict['nn']
            rough_shift = pdict['rough_shift']
		
        iv_thar = 1/((lines_thar/GAIN) + (RON**2/GAIN**2))

        All_Pixel_Centers = np.array([])
        All_Wavelengths   = np.array([])
        All_Orders        = np.array([])
        All_Centroids     = np.array([])
        All_Sigmas        = np.array([])
        All_Intensities   = np.array([])
        All_Residuals     = np.array([])
        All_Sigmas        = np.array([])

        orre  = or0
        order = orwa

        OK = []
        OW = []
	
        while order < n_last:
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)
            thar_order_orig = lines_thar[orre,:]
            IV              = iv_thar[orre,:]
            wei             = np.sqrt( IV )
            bkg             = GLOBALutils.Lines_mBack(thar_order_orig, IV, thres_rel=3, line_w=10)
            thar_order      = thar_order_orig - bkg

            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms,\
                              residuals, centroids, sigmas, intensities = GLOBALutils.Initial_Wav_Calibration( \
                              order_dir+'order_'+order_s+'o.iwdat', thar_order, order, wei, rmsmax=1000, \
			      minlines=10, FixEnds=False, Dump_Argon=False, Dump_AllLines=True, Cheby=use_cheby, \
                              rough_shift = rough_shift)

            fwhms_lns = sigmas*2.355
            inis_lns  = pixel_centers - fwhms_lns*0.5
            fins_lns  = pixel_centers + fwhms_lns*0.5			
            inis_wvs  = GLOBALutils.Cheby_eval(coeffs_pix2wav,inis_lns,float(len(thar_order)))
            fins_wvs  = GLOBALutils.Cheby_eval(coeffs_pix2wav,fins_lns,float(len(thar_order)))
            fwhms_wvs = inis_wvs - fins_wvs
            resolution2 = wavelengths / fwhms_wvs

            if wavelengths.max() > 5500 and wavelengths.min()<5500:
                print "\t\t\tmedian Resolution of order", order, '=', np.around(np.median(resolution2))
			
            if (order == 32): 
                if (use_cheby):
                    Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 1023, len(thar_order) )
                else:
                    Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

            All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
            All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
            All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order )
            All_Centroids     = np.append( All_Centroids, centroids)
            All_Sigmas        = np.append( All_Sigmas, sigmas)
            All_Intensities   = np.append( All_Intensities, intensities )
            All_residuals     = np.append( All_Residuals, residuals)
            All_Sigmas        = np.append( All_Sigmas,sigmas)
			
            order += 1
            orre  += 1

        p0    = np.zeros( npar_wsol )
        p0[0] = (32+oro0) * Global_ZP

        p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
			GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders, np.ones(All_Intensities.shape), p0, npix=len(lines_thar[orre,:]), Cheby=use_cheby, maxrms=150, Inv=Inverse_m, minlines=800,order0=oro0,nx=ncoef_x,nm=ncoef_m,ntotal=ntotal)

        if rms_ms/np.sqrt(float(len(G_pix))) < 20:
            pdict = {'p1':p1, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II,\
                     'rms_ms':rms_ms, 'G_res':G_res, 'All_Centroids':All_Centroids,\
                     'All_Sigmas':All_Sigmas, 'or0':or0, 'orwa':orwa, 'oro0':oro0,\
                     'rough_shift':rough_shift}
            pickle.dump( pdict, open( wavsol_pkl, 'w' ) )
            thtimes.append(thmjd)
            thRA.append(hth['RA-D'])
            thDEC.append(hth['DEC-D'])
            thnames.append(wavsol_pkl)
        else:
            bad_thars.append(counter)
    else:
        print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl
        pdict_temp = pickle.load(open(wavsol_pkl, 'r'))
        if pdict_temp['rms_ms']/np.sqrt(float(len(pdict_temp['G_pix']))) < 20:
            thtimes.append(thmjd)
            thRA.append(hth['RA-D'])
            thDEC.append(hth['DEC-D'])
            thnames.append(wavsol_pkl)
        else:
            bad_thars.append(counter)
    counter += 1

bad_thars = np.array(bad_thars)
if len(bad_thars)>0:
    ThAr_ref = np.delete(ThAr_ref,bad_thars)

pdict = pickle.load(open(dirout+'id_orders.pkl', 'r'))
or0         = pdict['or0']
orwa        = pdict['orwa']
n_last      = pdict['n_last']
nn          = pdict['nn']
rough_shift = pdict['rough_shift']

thtimes,thnames,thRA,thDEC = np.array(thtimes),np.array(thnames),np.array(thRA),np.array(thDEC)
SI = np.argsort(thtimes)
thtimes,thnames,thRA,thDEC = thtimes[SI],thnames[SI],thRA[SI],thDEC[SI]

"""
# get the reference P matrix
bkg_obj_fits = dirout + 'Bkg_ref.fits'
P_fits       = dirout + 'P_ref.fits'

dat = pyfits.getdata(dirin+stst)
hd  = pyfits.getheader(dirin+stst)

RON, GAIN = hd['ENOISE'], hd['EGAIN']
dat = utils.OverscanTrim(dat,hd['BIASSEC'])
dat -= MasterBias

if dark_substraction and hd['EXPTIME'] > 600:
	dat -= utils.get_dark(MDARKS, hd['EXPTIME'])
if bad_colummn:
	dat = utils.b_col(dat)
dat = dat/Flat
		
print 'Recentering traces...'
c_alls,pshift = GLOBALutils.retrace( dat, c_all )

Centers = np.zeros((len(c_alls),dat.shape[1]))
for i in range(nord):
	Centers[i,:]=scipy.polyval(c_alls[i,:],np.arange(len(Centers[i,:])))
print 'Scatter Light Determination...'
if np.median(dat[Centers[0,1024]:Centers[-1,1024],1024])> 0:
	if ( os.access(bkg_obj_fits,os.F_OK) == False or force_bkg):
		#bkg = bkg_simple.get_scat(dat,Centers)
		bkg = GLOBALutils.get_scat(dat,Centers,span=6)
		if (os.access(bkg_obj_fits,os.F_OK)):
			os.remove( bkg_obj_fits )
		hdu = pyfits.PrimaryHDU( bkg )
		hdu.writeto( bkg_obj_fits )
	else:
		bkg = pyfits.getdata(bkg_obj_fits)
	dat -= bkg

print 'P matrix determination...'
force_P = False
if os.access(P_fits,os.F_OK) == False or force_P:
	P  = GLOBALutils.obtain_P(dat,c_alls,ext_aperture,RON,GAIN,NSigma_Marsh, S_Marsh,N_Marsh, Marsh_alg,0,2047,npools)
	if (os.access(P_fits,os.F_OK)):
		os.remove( P_fits )
	hdu = pyfits.PrimaryHDU( P )
	hdu.writeto( P_fits )
else:
	P = pyfits.getdata(P_fits)
"""

print '\n\tProcessing of science images:'

sthd   = pyfits.getheader(dirin+stblaze)
stname = sthd['OBJECT']
new_list         = []
new_list_obnames = []
new_list_texp    = []
for i in range(len(objects)):
    fsim = objects[i]
    hd = pyfits.getheader(objects[i])
    obname = hd['OBJECT']
    texp = hd['EXPTIME']
    if object2do == 'all':
        new_list.append(fsim)
        new_list_obnames.append( obname )
        new_list_texp.append( texp )
    else:
        if (obname == object2do) or obname == stname:
            new_list.append(fsim)
            new_list_obnames.append( obname )
            new_list_texp.append( texp )
objects=new_list
objects = objects[:5]

for obj in objects:

	print '\n'
	print "\t--> Working on image: ", obj

	hd = pyfits.getheader(obj)
	nombre    = hd['OBJECT']
	RON, GAIN = hd['ENOISE'], hd['EGAIN']

	print "\t\tObject name:",nombre
	
	nama = nombre + '_' + hd['DATE-OBS'] + '_' + hd['UT-TIME'][:2] +\
               '-' + hd['UT-TIME'][3:5]+'-' + hd['UT-TIME'][6:]

	obj_fits        = dirout + nama + '.spec.fits.S'
	obj_fits_simple = dirout + nama + '.spec.simple.fits.S'
	bkg_obj_fits    = dirout + 'Bkg_' + nama + '.fits'	
	P_fits          = dirout + 'P_' + nama + '.fits'

	only_simple   = False
	if stname == nombre:
	    stpath = obj_fits_simple

	cond2 = False
	if only_simple:
		if ( os.access(obj_fits_simple,os.F_OK) == False ) or (force_sci_extract) :
			cond2 = True
	else:
		if ( os.access(obj_fits,os.F_OK) == False )  or ( os.access(obj_fits_simple,os.F_OK) == False ) or (force_sci_extract) or ( os.access(P_fits,os.F_OK) == False ):
			cond2 = True

	if cond2:

		print "\t\tNo previous extraction or extraction forced for science file", obj, "extracting..."

		dat = pyfits.getdata(obj)
		hdt = pyfits.getheader(obj)
		dat = dupontutils.OverscanTrim(dat,hdt['BIASSEC'])
		dat -= MasterBias
		if dark_substraction and hd['EXPTIME'] >= 899:
		    dat -= dupontutils.get_dark(MDARKS, hd['EXPTIME'])

		if bad_colummn:
		    dat = dupontutils.b_col(dat)
		dat = dat/Flat
		
		c_alls,pshift = GLOBALutils.retrace( dat, c_all )

		Centers = np.zeros((len(c_alls),dat.shape[1]))
		for i in range(nord):
	  		Centers[i,:]=scipy.polyval(c_alls[i,:],np.arange(len(Centers[i,:])))

		#print 'Scatter Light Determination...'

		#print 'SN',np.median(dat[Centers[0,1024]:Centers[-1,1024],1024])
		if np.median(dat[Centers[0,1024]:Centers[-1,1024],1024])> 0:
			if ( os.access(bkg_obj_fits,os.F_OK) == False or force_bkg):
				bkg = GLOBALutils.get_scat(dat,Centers,span=6)
				if (os.access(bkg_obj_fits,os.F_OK)):
					os.remove( bkg_obj_fits )
				hdu = pyfits.PrimaryHDU( bkg )
				hdu.writeto( bkg_obj_fits )
			else:
				bkg = pyfits.getdata(bkg_obj_fits)
			dat -= bkg
		else:
			NCosmic_Marsh = 6

		#print 'P matrix determination...'
		use_ref_P = False
		if use_ref_P:
			P = pyfits.getdata(dirout + 'P_ref.fits')
			Pn = P.copy()
			for i in range(P.shape[0]):
				j = i - pshift
				if j < 0 or j>=P.shape[0]:
					Pn[i] = 0.
				else:
					Pn[i] = P[j]
			P = Pn.copy()

		else:
			if os.access(P_fits,os.F_OK) == False or force_P:
				P = GLOBALutils.obtain_P(dat,c_alls,ext_aperture,RON,\
                                    GAIN,NSigma_Marsh, S_Marsh,N_Marsh, Marsh_alg,0,2047,npools)
				if (os.access(P_fits,os.F_OK)):
					os.remove( P_fits )
				hdu = pyfits.PrimaryHDU( P )
				hdu.writeto( P_fits )
			else:
				P = pyfits.getdata(P_fits)

		if only_simple==False:
			obj_S  = GLOBALutils.optimal_extraction(dat,P,c_alls,ext_aperture,RON,\
                                 GAIN,S_Marsh,NCosmic_Marsh,50,2000,npools)
			obj_S  = obj_S[::-1]
			if (os.access(obj_fits,os.F_OK)):
				os.remove( obj_fits )
			hdu = pyfits.PrimaryHDU( obj_S )
			hdu.writeto( obj_fits )

		obj_Ss = GLOBALutils.simple_extraction(dat,c_alls,ext_aperture,min_extract_col,max_extract_col,npools)
		obj_Ss = obj_Ss[::-1]
		
		if (os.access(obj_fits_simple,os.F_OK)):
			os.remove( obj_fits_simple )
		
		hdu = pyfits.PrimaryHDU( obj_Ss )
		hdu.writeto( obj_fits_simple )

###################################################################################################################
############################################### get blaze function ################################################
###################################################################################################################
print '\n\tComputing the blaze function using the standard star...'

if os.access(dirout+'blaze_f.fits',os.F_OK) == False or force_bl:		
    index = -1
    hd = pyfits.getheader(dirin+stblaze)
    nombre  = hd['OBJECT']
    exptime = hd['EXPTIME']	
    RA      = hd['RA-D']
    DEC     = hd['DEC-D']
    RON     = hd['ENOISE']
    GAIN    = hd['EGAIN']
    scmjd,scmjd0 = dupontutils.mjd_fromheader(hd)
    altitude     = hd['SITEALT']
    latitude     = hd['SITELAT']
    longitude    = hd['SITELONG']
    epoch        = hd['EPOCH']

    iers          = GLOBALutils.JPLiers( baryc_dir, scmjd-999.0, scmjd+999.0 )
    obsradius, R0 = GLOBALutils.JPLR0( latitude, altitude)
    obpos         = GLOBALutils.obspos( longitude, obsradius, R0 )

    jplephem.set_ephemeris_dir( baryc_dir , ephemeris )
    jplephem.set_observer_coordinates( obpos[0], obpos[1], obpos[2] )

    res           = jplephem.doppler_fraction(RA/15.0, DEC, int(scmjd), scmjd%1, 1, 0.0)
    lbary_ltopo   = 1.0 + res['frac'][0]
    bcvel_baryc   = ( lbary_ltopo - 1.0 ) * 2.99792458E5

    res  = jplephem.pulse_delay(RA/15.0, DEC, int(scmjd), scmjd%1, 1, 0.0)   
    mbjd = scmjd + res['delay'][0] / (3600.0 * 24.0)

    gobs      = ephem.Observer()
    gobs.name = 'DUPONT'
    gobs.lat  = rad(latitude)  # lat/long in decimal degrees  
    gobs.long = rad(longitude)
    gobs.date = hd['UT-DATE'].replace('-','/') + ' ' + hd['UT-TIME']

    mephem = ephem.Moon()
    mephem.compute(gobs)
    Mcoo = jplephem.object_track("Moon", int(scmjd), float(scmjd%1), 1, 0.0)
    Mp   = jplephem.barycentric_object_track("Moon", int(scmjd), float(scmjd%1), 1, 0.0)
    Sp   = jplephem.barycentric_object_track("Sun", int(scmjd), float(scmjd%1), 1, 0.0)
    res  = jplephem.object_doppler("Moon", int(scmjd), scmjd%1, 1, 0.0)
    lunation,moon_state,moonsep2,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,RA,DEC)
    refvel = bcvel_baryc + moonvel
    print '\t\tBarycentric velocity:',refvel

    i=0
    mind = 1000.0
    for t in thtimes:
        if abs(t-scmjd)<mind:
            mind  = abs(t-scmjd)
            index = i
        i+=1

    obj_S = pyfits.getdata(stpath)

    wavsol_pkl1 = thnames[index]
    pdict1      = pickle.load(open(wavsol_pkl1,'r'))
    global1     = pdict1['p1']
    or01        = pdict1['or0']
    orwa1       = pdict1['orwa']
    oro0        = pdict1['oro0']

    final  = np.zeros( [3, n_last-orwa1, np.shape(obj_S)[1]] )
    equis  = np.arange( np.shape(obj_S)[1] ) 
    order  = orwa1
    orre   = or01
    final2 = np.zeros( [nn, np.shape(obj_S)[1]] )
    while order < n_last:
        m = order + oro0
        chebs = GLOBALutils.Calculate_chebs(equis, m, order0=oro0, ntotal=ntotal, npix=obj_S.shape[1], Inverse=Inverse_m, nx=ncoef_x,nm=ncoef_m )
        WavSol = lbary_ltopo*(1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(global1,chebs,ncoef_x,ncoef_m)
        final[0,orre,:] = GLOBALutils.ToVacuum(WavSol[::-1])
        final[1,orre,:] = obj_S[orre,:][::-1]
        order += 1
        orre += 1

    final[2,:,:] = dupontutils.get_blaze(final[0,:,:],final[1,:,:])
		
    blaze = np.zeros((2,final.shape[1],final.shape[2]))
    blaze[0] = final[0]
    blaze[1] = dupontutils.get_blaze(final[0,:,:],final[1,:,:])
	
    if (os.access(dirout+'blaze_f.fits',os.F_OK)):
        os.remove( dirout+'blaze_f.fits' )
    hdu = pyfits.PrimaryHDU( blaze )
    hdu.writeto(dirout+'blaze_f.fits'  )
else:
    blaze = pyfits.getdata(dirout+'blaze_f.fits' )


for i in range(blaze.shape[1]):
    blaze[1,i]/=blaze[1,i].max()

#################################################################################################################
############################################  Final output ######################################################
#################################################################################################################

print "\n\tBuilding the final output spectra..."

for obj in objects:
    hd     = pyfits.getheader(obj)
    nombre = hd['OBJECT']
    nama   = nombre+'_'+hd['DATE-OBS']+'_'+hd['UT-TIME'][:2]+'-'+hd['UT-TIME'][3:5]+'-'+hd['UT-TIME'][6:]
    nf     = nama+'_final.fits'
    nfS    = nama+'_finalS.fits'
    print "\n\t\t-->Building", nf
    if (os.access(dirout+'proc/'+nf,os.F_OK)) == False or force_sci_proc:
			index1,index2 = -1,-1
			
			# Get observing info from header
			hd = pyfits.getheader(obj)
			nombre    = hd['OBJECT']
			exptime   = hd['EXPTIME']
			RA        = hd['RA-D']
			DEC       = hd['DEC-D']
			RON       = hd['ENOISE']
			GAIN      = hd['EGAIN']
			altitude  = hd['SITEALT']
			latitude  = hd['SITELAT']
			longitude = hd['SITELONG']
			epoch     = hd['EPOCH']

			scmjd,scmjd0 = dupontutils.mjd_fromheader(hd)
			ra2,dec2 = GLOBALutils.getcoords(obname,scmjd,filen=reffile)
			if ra2 !=0 and dec2 != 0:
				RA = ra2
				DEC = dec2
			else:
				print '\t\tUsing the coordinates found in the image header.'

			# set info for compute the baricentric correction
			iers          = GLOBALutils.JPLiers( baryc_dir, scmjd-999.0, scmjd+999.0 )
			obsradius, R0 = GLOBALutils.JPLR0( latitude, altitude)
			obpos         = GLOBALutils.obspos( longitude, obsradius, R0 )
			jplephem.set_ephemeris_dir( baryc_dir , ephemeris )
		    	jplephem.set_observer_coordinates( obpos[0], obpos[1], obpos[2] )
			res         = jplephem.doppler_fraction(RA/15.0, DEC, int(scmjd), scmjd%1, 1, 0.0)
			lbary_ltopo = 1.0 + res['frac'][0]
		   	bcvel_baryc = ( lbary_ltopo - 1.0 ) * 2.99792458E5	#This in the barycentric velocity
		   	res         = jplephem.pulse_delay(RA/15.0, DEC, int(scmjd), scmjd%1, 1, 0.0)
		   	scmbjd      = scmjd + res['delay'][0] / (3600.0 * 24.0)	#This is the modified barycentric julian day of the observation

			# set observatory info to retrive info about the moon
			gobs = ephem.Observer()
			gobs.name = 'DUPONT'
			gobs.lat  = rad(latitude)
			gobs.long = rad(longitude)
			#gobs.date = hd['UT-DATE'] + ' ' + hd['UT-TIME'].replace(':','_')
			gobs.date = hd['UT-DATE'].replace('-','/') + ' ' + hd['UT-TIME']
		
			mephem = ephem.Moon()
			mephem.compute(gobs)
			Mcoo = jplephem.object_track("Moon", int(scmjd), float(scmjd%1), 1, 0.0)
			Mp   = jplephem.barycentric_object_track("Moon", int(scmjd), float(scmjd%1), 1, 0.0)
			Sp   = jplephem.barycentric_object_track("Sun", int(scmjd), float(scmjd%1), 1, 0.0)
		   	res      = jplephem.object_doppler("Moon", int(scmjd), scmjd%1, 1, 0.0)
			lunation,moon_state,moonsep2,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,RA,DEC)
			refvel = bcvel_baryc + moonvel	#This is the velocity of the spectrum of the moon with the applied barycentric correction in the direction of the target. 
			
			print '\t\t\tBarycentric velocity:',refvel

			# Set the ThAr lamps for applying the wavelength solution
			if scmjd < thtimes[0]:
				"\t\t\tProblem with ThAr and science times"
				index1 = 0
				index2 = 0
			elif scmjd > thtimes[-1]:
				"\t\t\tProblem with ThAr and science times"
				index1 = -1
				index2 = -1
			else:
				"\t\t\tThAr images taken before and after science image found..."
				for i in range(len(thtimes)-1):
					if scmjd >= thtimes[i] and scmjd < thtimes[i+1]:
						index1 = i
						index2 = i+1
						break
			
			if abs(RA-thRA[index1]) > 0.005 or abs(DEC-thDEC[index1]) > 0.005:
				index1 = index2
			if abs(RA-thRA[index2]) > 0.005 or abs(DEC-thDEC[index2]) > 0.005:
				index2 = index1
		
			obj_fits        = dirout+nama+'.spec.fits.S'
			obj_fits_simple = dirout+nama+'.spec.simple.fits.S'
			obj_S           = pyfits.getdata(obj_fits)
			obj_Ss          = pyfits.getdata(obj_fits_simple)
	
			wavsol_pkl1 = thnames[index1]
			pdict1      = pickle.load(open(wavsol_pkl1,'r'))
			global1     = pdict1['p1']
			or01        = pdict1['or0']
			orwa1       = pdict1['orwa']
			oro0        = pdict1['oro0']
			wavsol_pkl2 = thnames[index2]
			pdict2      = pickle.load(open(wavsol_pkl2,'r'))
			global2     = pdict2['p1']

			final  = np.zeros( [11, nn, np.shape(obj_S)[2]] )
			finalS = np.zeros( [2, nn, np.shape(obj_Ss)[1]] )
			equis  = np.arange( np.shape(obj_S)[2] ) 
			order  = orwa1
			orre   = or01

			while order < n_last:
				
				m = order + oro0
				chebs = GLOBALutils.Calculate_chebs(equis, m, order0=oro0, ntotal=ntotal, npix=obj_S.shape[2], Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
				WavSol  = lbary_ltopo*(1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(global1,chebs,ncoef_x,ncoef_m)
				WavSol2 = lbary_ltopo*(1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(global2,chebs,ncoef_x,ncoef_m)
				
				final[0,orre,:] = GLOBALutils.ToVacuum(0.5*(WavSol[::-1]+WavSol2[::-1]))
				final[1,orre,:] = obj_S[orre,1,:][::-1]
				final[2,orre,:] = obj_S[orre,2,:][::-1]
				tck = scipy.interpolate.splrep(blaze[0,orre,:],blaze[1,orre,:],k=3,s=0)
				blin = scipy.interpolate.splev(final[0,orre,:],tck,der = 0)
				final[3,orre,:] = final[1,orre,:]/blin
				final[4,orre,:] = final[2,orre,:]*(blin**2)
				
				X = final[0,orre,:]
				Z = final[3,orre,:]
				Z[:150]  = 0.0
				Z[-150:] = 0.0

				I    = np.isnan(Z)
				Z[I] = 0.0
				I    = np.where(Z!=0)[0]
				if len(I)>0:			
					cp    = continuum.NORM_single(X,Z,orden=1)
					ratio = np.polyval(cp,X)
				else:
					ratio = np.ones(len(X))

				ratio = np.polyval(cp,X)
				final[5,orre,:] = final[3,orre,:]/ratio
				Inan = np.where( np.isnan(final[1,orre,:]) == True )[0]
				final[5,orre,Inan] = 1.
				final[6,orre,:] = final[4,orre,:]*(ratio**2)
				final[7,orre,:] = ratio
				final[8,orre,:] = ratio*blin / np.sqrt( ratio *blin / GAIN + (RON/GAIN)**2 )
				
				spl        = scipy.interpolate.splrep(np.arange(len(final[0,orre,:])),final[0,orre,:] ,k=3)
				dlambda_dx = scipy.interpolate.splev(np.arange(len(final[0,orre,:])), spl, der=1)
				NN         = np.average(dlambda_dx)
				dlambda_dx /= NN
				
				LL               = np.where( final[5,orre,:] != 0 )[0]
				medflx           = scipy.signal.medfilt( final[5,orre,LL], 3 )
				fres             = final[5,orre,LL] - medflx
				fdev             = np.sqrt(np.var(fres))
				II               = np.where(final[5,orre,:] > 1. + 5*fdev)[0]
				final[5,orre,II] = 1.

				final[9,orre,:]  = final[5,orre,:] * (dlambda_dx ** 1)
				final[10,orre,:] = final[6,orre,:] / (dlambda_dx ** 2)

				finalS[0,orre,:] = WavSol[::-1]
				finalS[1,orre,:] = obj_Ss[orre,:][::-1]


				order += 1
				orre += 1

			hdu = pyfits.PrimaryHDU( final )
			hdu = GLOBALutils.update_header(hdu,'HIERARCH MJD',scmjd)
			hdu = GLOBALutils.update_header(hdu,'HIERARCH MBJD',scmbjd)
			hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START DATE',hd['DATE-OBS'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START UT',hd['UT-TIME'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH TEXP (S)',hd['EXPTIME'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (KM/S)',bcvel_baryc,'[km/s]')
			hdu = GLOBALutils.update_header(hdu,'HIERARCH (LAMBDA_BARY / LAMBDA_TOPO)',lbary_ltopo)
			hdu = GLOBALutils.update_header(hdu,'HIERARCH TARGET NAME',nombre)
			hdu = GLOBALutils.update_header(hdu,'HIERARCH RA',hd['RA'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC',hd['DEC'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH RA-D',hd['RA-D'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC-D',hd['DEC-D'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH RA BARY',RA)
			hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC BARY',DEC)			
			hdu = GLOBALutils.update_header(hdu,'HIERARCH EQUINOX',hd['EQUINOX'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LATITUDE',hd['SITELAT'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LONGITUDE',hd['SITELONG'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS ALTITUDE',hd['SITEALT'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH TARG AIRMASS',hd['AIRMASS'])
			hdu = GLOBALutils.update_header(hdu,'HIERARCH MOON_VEL',refvel,'[km/s]')
			hdu = GLOBALutils.update_header(hdu,'HIERARCH MOONST',moon_state)
			hdu = GLOBALutils.update_header(hdu,'HIERARCH LUNATION',lunation)
			hdu = GLOBALutils.update_header(hdu,'HIERARCH MOONSEP',moonsep2)
			hdu = GLOBALutils.update_header(hdu,'HIERARCH MOONALT',float(mephem.alt))
			hdu = GLOBALutils.update_header(hdu,'HIERARCH SMOONALT',str(mephem.alt))
		
			if (os.access(dirout+'proc/'+nf,os.F_OK)):
					os.remove( dirout+'proc/'+nf )
			hdu.writeto( dirout+'proc/'+nf )
		
#####################################################################################################################
############################################### Spectral and RV analysis ############################################
#####################################################################################################################

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

print 'Starting with the Post-processing of the spectra...'

if not JustExtract:

    for obj in objects:
        hd     = pyfits.getheader(obj)
        nombre = hd['OBJECT']
        nama   = dirout+'proc/' + nombre+'_'+hd['DATE-OBS']+'_'+hd['UT-TIME'][:2]+'-'+hd['UT-TIME'][3:5]+'-'+hd['UT-TIME'][6:]
        fit    = nama+'_final.fits'
        print '\n\t\tWorking on spectrum:', fit
        know_moon = False
        if fsim.split('/')[-1] in spec_moon:
            I = np.where(obj.split('/')[-1] == spec_moon)[0]
            know_moon = True
            here_moon = use_moon[I]

        hd   = pyfits.getheader(fit)
        spec = pyfits.getdata(fit)
        hdu  = pyfits.open(fit,'update')

        obname      = hd['TARGET NAME']
        refvel      = hd['MOON_VEL']
        lunation    = hd['LUNATION']
        moon_state  = hd['MOONST']
        moonsep     = hd['MOONSEP']
        moon_alt    = hd['MOONALT']
        smoon_alt   = hd['SMOONALT']
        lbary_ltopo = hd['(LAMBDA_BARY / LAMBDA_TOPO)']
        mjd         = hd['MJD']
        TEXP        = hd['TEXP (S)']
        mbjd        = hd['MBJD']

        for i in range(spec.shape[1]):
            if spec[0,i,500] < 5200:
                tuc = i
                break
        SNR_5130 = np.median(spec[8,tuc][800:1200] )

        if DoClass:
            print '\t\tSpectral Analysis:'
            # spectral analysis
            # First, query SIMBAD with the object name
            query_success = False
            query_success,sp_type_query = GLOBALutils.simbad_query_obname(obname)
            # Now, query SIMBAD by coordinates if above not successful
            if (not query_success):
                query_success,sp_type_query = GLOBALutils.simbad_query_coords('12:00:00','00:00:00')
            print "\t\t\tSpectral type returned by SIMBAD query:",sp_type_query

            hdu[0] = GLOBALutils.update_header(hdu[0],'HIERARCH SIMBAD SPTYP', sp_type_query)

            pars_file = dirout + nombre+'_'+hd['HIERARCH SHUTTER START DATE']+'_' \
				+ hd['HIERARCH SHUTTER START UT'][:2]+'-'+hd['HIERARCH SHUTTER START UT'][3:5]+ \
				'-' + hd['HIERARCH SHUTTER START UT'][6:] + '_stellar_pars.txt'

            if os.access(pars_file,os.F_OK) == False or force_stellar_pars:
                print "\t\t\tEstimating atmospheric parameters:"
                spec2 = spec.copy()
                if resolution > 45000:
                    Rx = np.around(1./np.sqrt(1./40000.**2 - 1./resolution**2))
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

        # store the parameters measured for this epoch
        T_eff_epoch = T_eff
        logg_epoch  = logg
        Z_epoch     = Z
        vsini_epoch = vsini
        vel0_epoch  = vel0
	
        hdu[0] = GLOBALutils.update_header(hdu[0],'HIERARCH VEL0', vel0)
        hdu[0] = GLOBALutils.update_header(hdu[0],'HIERARCH TEFF', float(T_eff))
        hdu[0] = GLOBALutils.update_header(hdu[0],'HIERARCH LOGG', float(logg))
        hdu[0] = GLOBALutils.update_header(hdu[0],'HIERARCH Z', Z)
        hdu[0] = GLOBALutils.update_header(hdu[0],'HIERARCH VSINI', vsini)

        print "\t\tRadial Velocity analysis:"
        # assign mask
        sp_type, mask = GLOBALutils.get_mask_reffile(obname,reffile=reffile,base='../data/xc_masks/')
        print "\t\t\tWill use",sp_type,"mask for CCF."

        first_o = 0
        for i in range(spec.shape[1]):
            if spec[0,i,0]<6500:
                first_o = i
                break

        last_o = spec.shape[1]
        for i in range(spec.shape[1]):
            if spec[0,i,0]<4300:
                last_o = i
                break

        spec1 = spec[:,first_o:last_o+1,:]
	
        # Read in mask
        ml, mh, weight = np.loadtxt(mask,unpack=True)
        ml_v = GLOBALutils.ToVacuum( ml )
        mh_v = GLOBALutils.ToVacuum( mh )

        # make mask larger accounting for factor ~3 lower res in DUPONT w/r to HARPS
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
            #first rough correlation to find the minimum
            vels, xc_full, sn, nlines_ccf, W_ccf = GLOBALutils.XCor(spec1, ml_v, mh_v,\
		                            weight, 0, lbary_ltopo, vel_width=300, vel_step=3, start_order=0,\
		                            spec_order=9, iv_order=10, sn_order=8,max_vel_rough=300)
            xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3, Simple=True, start_order=0, W=W_ccf)
            #Normalize the continuum of the CCF robustly with R     
            yy = scipy.signal.medfilt(xc_av,11)
            pred = lowess(yy, vels,frac=0.4,it=10,return_sorted=False)
            tck1 = scipy.interpolate.splrep(vels,pred,k=1)
            xc_av_orig = xc_av.copy()
            xc_av /= pred
            pred_rough = pred.copy()

            vel0_xc = vels[ np.argmin( xc_av ) ] 

            rvels, rxc_av, rpred, rxc_av_orig, rvel0_xc = vels.copy(), \
				xc_av.copy(), pred.copy(), xc_av_orig.copy(), vel0_xc
            xc_av_rough = xc_av
            vels_rough  = vels
		        
            vel_width = np.maximum( 20.0, 6*disp )
            vels, xc_full, sn, nlines_ccf, W_ccf =\
					GLOBALutils.XCor(spec1, ml_v, mh_v, weight, vel0_xc, lbary_ltopo,\
		                        start_order=0, vel_width=vel_width, vel_step=0.1, spec_order=9, \
		                        iv_order=10, sn_order=8,max_vel_rough=300)
		        
            xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3,\
		                        Simple=True, W=W_ccf, start_order=0)

            pred = scipy.interpolate.splev(vels,tck1)
            xc_av /= pred

            if sp_type == 'M5':
                moon_sig = 2.5
            elif sp_type == 'K5':
                moon_sig = 3.3
            else:
                moon_sig = 4.5

            p1,XCmodel,p1gau,XCmodelgau,Ls2 = \
						GLOBALutils.XC_Final_Fit( vels, xc_av , sigma_res=4,\
		                                horder=8, moonv = refvel, moons = moon_sig, moon = False)

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
		         'lunation':lunation,'mephem':mephem,'texp':TEXP}

        pkl_xc = fit[:-4] + obname + '_XC_' + sp_type + '.pkl'
        pickle.dump( xc_dict, open( pkl_xc, 'w' ) )
        ccf_pdf = dirout + 'proc/' + fit.split('/')[-1][:-4] + obname + '_XCs_' + sp_type + '.pdf'
        if not avoid_plot:
            GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)

        airmass  = hdu[0].header['TARG AIRMASS']
        seeing   = -999

        if sp_type == 'G2':
            if T_eff < 6000:
                D = 0.24416
                C = 0.00181
            else:
                D = 0.33491
                C = 0.00113
        elif  sp_type == 'K5':	
            D = 0.20695
            C = 0.00321
        else:	
            D = 0.20695
            C = 0.00321

        RV     = np.around(p1gau_m[1],4)  
        BS     = np.around(SP,4)  
        RVerr2 = 0.400
        BSerr  = np.around(D / float(np.round(SNR_5130)) + C,4)

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
        hdu[0] = GLOBALutils.update_header(hdu[0],'RV', RV)
        hdu[0] = GLOBALutils.update_header(hdu[0],'RV_E', RVerr2)
        hdu[0] = GLOBALutils.update_header(hdu[0],'BS', BS)
        hdu[0] = GLOBALutils.update_header(hdu[0],'BS_E', BSerr)
        hdu[0] = GLOBALutils.update_header(hdu[0],'DISP', disp_epoch)
        hdu[0] = GLOBALutils.update_header(hdu[0],'SNR', SNR_5130)
        hdu[0] = GLOBALutils.update_header(hdu[0],'SNR_R', SNR_5130_R)
        hdu[0] = GLOBALutils.update_header(hdu[0],'INST', 'DUPONT')
        hdu[0] = GLOBALutils.update_header(hdu[0],'RESOL', resolution)
        hdu[0] = GLOBALutils.update_header(hdu[0],'PIPELINE', 'CERES')
        hdu[0] = GLOBALutils.update_header(hdu[0],'XC_MIN', XC_min)
        hdu[0] = GLOBALutils.update_header(hdu[0],'BJD_OUT', bjd_out)

        line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f   dupont   ceres   %8d %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
                      (obname, bjd_out, RV, RVerr2, BS, BSerr, resolution, T_eff_epoch, logg_epoch,\
		       Z_epoch, vsini_epoch, XC_min, disp_epoch, TEXP, SNR_5130_R, ccf_pdf)
        f_res.write(line_out)
        hdu.close()
f_res.close()