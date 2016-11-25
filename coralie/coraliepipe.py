import sys
from pylab import *

base = '../'

sys.path.append(base+"utils/Continuum")
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
import coralieutils
import continuum
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
force_pre_process  = False
force_flat_extract = False
force_thar_extract = False
force_thfp_extract = False
force_tharxc       = False
force_thar_wavcal  = False
force_thfp_wavcal  = False
force_sci_extract  = False
force_spectral_file_build = True
force_stellar_pars = False
dumpargon          = False
minlines_glob_ob   = 700
minlines_glob_co   = 500

Inverse_m          = True
use_cheby          = True
MRMS               = 100   # max rms in m/s, global wav solution

trace_degree       = 4
Marsh_alg          = 0
ext_aperture       = 3
NSigma_Marsh       = 5
NCosmic_Marsh      = 10
S_Marsh            = 0.4
N_Marsh            = 3      # grado polinomio 
min_extract_col    = 50
max_extract_col    = 2000
n_useful           = 70    # up to which order do we care?

# Number of coefficients for the global wavelength solution
ncoef_x            = 4
ncoef_m            = 6
npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2

models_path = base+"data/COELHO_MODELS/R_40000b/"    # path to the synthetic models 
order_dir   = base+"coralie/wavcals/"  # path to reference files for the wavelength solution

#############################

# file containing the log
log = dirout+'night.log'

print "\n\n\tCoralie Euler1.2m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'

# classification of input images according to header info
biases, ob_flats, co_flats, ob_loc, co_loc, ThAr_ref, ThFP_ref,\
        simThAr_sci,sim_FP_sci,ThAr_ref_dates,ThFP_ref_dates,obnames,\
	obnames_FP,exptimes, exptimes_FP = coralieutils.FileClassify(dirin,log)

# Pre-process
if ( (os.access(dirout+'FlatOb.fits',os.F_OK) == False) or \
     (os.access(dirout+'FlatCo.fits',os.F_OK) == False) or \
     (os.access(dirout+'trace.pkl',os.F_OK) == False)  or \
     (os.access(dirout+'MasterBias.fits',os.F_OK) == False)  or \
     (force_pre_process) ):
    print "\tNo previous pre-processing files or found"
    pre_process = 1
else:
    print "\tPre-processing files found, going straight to extraction"
    pre_process = 0

if (pre_process == 1):
    print "\tGenerating Master calibration frames..."
    # median combine Biases
    MasterBias, RO_bias, GA_bias = coralieutils.MedianCombine(biases,ZF=0.)
    hdu = pyfits.PrimaryHDU( MasterBias )
    if (os.access(dirout+'MasterBias.fits',os.F_OK)):
        os.remove(dirout+'MasterBias.fits')
    hdu.writeto(dirout+'MasterBias.fits')
    print "\t\t-> Masterbias: done!"
    # median combine list of ob flats
    Flat_ob, RO_ob, GA_ob = coralieutils.MedianCombine(ob_flats,ZF=MasterBias)
    # save this file for later reference
    hdu = pyfits.PrimaryHDU( Flat_ob )
    if (os.access(dirout+'FlatOb.fits',os.F_OK)):
        os.remove(dirout+'FlatOb.fits')
    hdu.writeto(dirout+'FlatOb.fits')

    # median combine list of co flats
    Flat_co,RO_co,GA_co = coralieutils.MedianCombine(co_flats,ZF=MasterBias)
    hdu = pyfits.PrimaryHDU(Flat_co)
    if (os.access(dirout+'FlatCo.fits',os.F_OK)):
        os.remove(dirout+'FlatCo.fits')
    hdu.writeto(dirout+'FlatCo.fits')
    print "\t\t-> Masterflats: done!"
    
    # Find orders & traces
    print "\tTracing echelle orders..."
    c_ob, nord_ob = GLOBALutils.get_them(Flat_ob, 8, trace_degree,maxords=-1,mode=1)
    c_co, nord_co = GLOBALutils.get_them(Flat_co, 8, trace_degree,maxords=-1,startfrom=300,mode=1)

    trace_dict = {'c_ob':c_ob,
                  'c_co':c_co,
                  'nord_ob':nord_ob, 'nord_co':nord_co,
                  'GA_ob': GA_ob, 'RO_ob': RO_ob,
                  'GA_co': GA_co, 'RO_co': RO_co}
    pickle.dump( trace_dict, open( dirout+"trace.pkl", 'w' ) )

else:
    trace_dict = pickle.load( open( dirout+"trace.pkl", 'r' ) )
    c_co = trace_dict['c_co']
    c_ob = trace_dict['c_ob']
    nord_ob = trace_dict['nord_ob']
    nord_co = trace_dict['nord_co']
    # recover GA*, RO*
    GA_ob = trace_dict['GA_ob']
    RO_ob = trace_dict['RO_ob']
    GA_co = trace_dict['GA_co']
    RO_co = trace_dict['RO_co']
    # recover flats & master bias
    h = pyfits.open(dirout+'FlatOb.fits')
    Flat_ob = h[0].data
    h = pyfits.open(dirout+'FlatCo.fits')
    Flat_co = h[0].data
    h = pyfits.open(dirout+'MasterBias.fits')
    MasterBias = h[0].data

c_all = GLOBALutils.Mesh(c_ob,c_co)
print '\n\tExtraction of Flat calibration frames:'
# Extract flat spectra, object
P_ob_fits = dirout + 'P_ob.fits'
S_flat_ob_fits = dirout +'S_flat_ob.fits'
P_ob = np.zeros( Flat_ob.shape )
S_flat_ob = np.zeros((nord_ob, 3, Flat_ob.shape[1]) )
if ( os.access(P_ob_fits,os.F_OK) == False ) or ( os.access(S_flat_ob_fits,os.F_OK) == False ) or \
   (force_flat_extract):
    print "\t\tNo extracted flat object spectra found or extraction forced, extracting and saving..."
    
    print "\t\t\tWill extract",nord_ob,"orders for object fibre..."
    P_ob = GLOBALutils.obtain_P(Flat_ob,c_ob,ext_aperture,RO_ob,\
                                    GA_ob,NSigma_Marsh, S_Marsh, \
                    N_Marsh, Marsh_alg, min_extract_col,\
                    max_extract_col, npools)

    S_flat_ob  = GLOBALutils.optimal_extraction(Flat_ob,P_ob,c_ob,ext_aperture,\
                                                RO_ob,GA_ob,S_Marsh,NCosmic_Marsh,\
                                                min_extract_col,max_extract_col,npools)

    # write P_on and S_flat_ob as fits files
    if (os.access(P_ob_fits,os.F_OK)):
        os.remove( P_ob_fits )
    if (os.access(S_flat_ob_fits,os.F_OK)):
        os.remove( S_flat_ob_fits )
    
    hdu = pyfits.PrimaryHDU( P_ob )
    hdu.writeto( P_ob_fits )
    hdu = pyfits.PrimaryHDU( S_flat_ob )
    hdu.writeto( S_flat_ob_fits )

   
else:
    print "\t\tExtracted flat object spectra found, loading..."
    P_ob       = pyfits.getdata( P_ob_fits )
    S_flat_ob  = pyfits.getdata( S_flat_ob_fits )

# Extract flat spectra, comparison
P_co_fits = dirout + 'P_co.fits'
S_flat_co_fits = dirout +'S_flat_co.fits'
P_co = np.zeros( Flat_co.shape )
S_flat_co = np.zeros((nord_co, 3, Flat_co.shape[1]) )
if ( os.access(P_co_fits,os.F_OK) == False ) or ( os.access(S_flat_co_fits,os.F_OK) == False ) or (force_flat_extract):
    print "\t\tNo extracted flat comparison spectra found or extraction forced, extracting and saving..."
    
    print "\t\t\tWill extract",nord_co,"orders for comparison fibre"
    P_co = GLOBALutils.obtain_P(Flat_co,c_co,ext_aperture,RO_co,\
                                    GA_co,NSigma_Marsh, S_Marsh, \
                    N_Marsh, Marsh_alg, min_extract_col,\
                    max_extract_col, npools)
    
    S_flat_co  = GLOBALutils.optimal_extraction(Flat_co,P_co,c_co,ext_aperture,RO_co,GA_co,\
                                                S_Marsh,NCosmic_Marsh,min_extract_col,\
                                                max_extract_col,npools) 

    # write P_on and S_flat_co as fits files
    if (os.access(P_co_fits,os.F_OK)):
        os.remove( P_co_fits )
    if (os.access(S_flat_co_fits,os.F_OK)):
        os.remove( S_flat_co_fits )

    hdu = pyfits.PrimaryHDU( P_co )
    hdu.writeto( P_co_fits )
    hdu = pyfits.PrimaryHDU( S_flat_co )
    hdu.writeto( S_flat_co_fits )
   
else:
    print "\t\tExtracted flat comparison spectra found, loading..."
    P_co       = pyfits.getdata( P_co_fits )
    S_flat_co  = pyfits.getdata( S_flat_co_fits )

# Normalize flat field spectra.
S_flat_ob_n, maxvals_ob = GLOBALutils.FlatNormalize_single( S_flat_ob, mid=int(0.5*S_flat_ob.shape[2]))
S_flat_co_n, maxvals_co = GLOBALutils.FlatNormalize_single( S_flat_co, mid=int(0.5*S_flat_co.shape[2]))

print '\n\tExtraction of ThAr calibration frames:'
# Extract all ThAr files
for fsim in ThAr_ref:
    hthar = pyfits.open( fsim )
    dthar = coralieutils.OverscanTrim( pyfits.getdata( fsim ) )
    ron  = hthar[0].header['HIERARCH ESO CORA CCD RON']
    gain = hthar[0].header['HIERARCH ESO CORA CCD GAIN']
    thar_fits_ob = dirout + fsim.split('/')[-1][:-8]+'spec.ob.fits.S'
    thar_fits_co = dirout + fsim.split('/')[-1][:-8]+'spec.co.fits.S'  
    if ( os.access(thar_fits_ob,os.F_OK) == False ) or \
       ( os.access(thar_fits_co,os.F_OK) == False ) or \
       (force_thar_extract):

        print "\t\tNo previous extraction or extraction forced for ThAr file", fsim, "extracting..."
	thar_S_ob  = GLOBALutils.optimal_extraction(dthar,P_ob,c_ob,ext_aperture,ron,gain,\
                                                    S_Marsh,100.,min_extract_col,max_extract_col,npools)
	thar_S_co  = GLOBALutils.optimal_extraction(dthar,P_co,c_co,ext_aperture,ron,gain,\
                                                    S_Marsh,100.,min_extract_col,max_extract_col,npools)           
       
        # save as fits file
        if (os.access(thar_fits_ob,os.F_OK)):
            os.remove( thar_fits_ob )
        if (os.access(thar_fits_co,os.F_OK)):
            os.remove( thar_fits_co )
            
        hdu = pyfits.PrimaryHDU( thar_S_ob )
        hdu.writeto( thar_fits_ob )
        hdu = pyfits.PrimaryHDU( thar_S_co )
        hdu.writeto( thar_fits_co )
    else:
        print "\t\tThAr file", fsim, "all ready extracted, loading..."

print "\n\tWavelength solution of ThAr calibration spectra:"
# compute wavelength calibration files
sorted_ThAr_dates = np.argsort( ThAr_ref_dates )
p0_array = np.zeros( (len(ThAr_ref_dates), npar_wsol) )

for i in range(len(sorted_ThAr_dates)):
    index = sorted_ThAr_dates[i]
    wavsol_pkl   = dirout + ThAr_ref[index].split('/')[-1][:-8]+'wavsolpars.pkl'
    thar_fits_ob = dirout + ThAr_ref[index].split('/')[-1][:-8]+'spec.ob.fits.S'
    thar_fits_co = dirout + ThAr_ref[index].split('/')[-1][:-8]+'spec.co.fits.S'
    if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
        print "\t\tComputing wavelength solution of ThAr file", ThAr_ref[index] 
        hthar = pyfits.open( ThAr_ref[index] )
        mjd, mjd0 = coralieutils.mjd_fromheader( hthar )
        thar_S_ob = pyfits.getdata( thar_fits_ob )
        thar_S_co = pyfits.getdata( thar_fits_co )

        lines_thar_ob  = thar_S_ob[:,1,:]
        iv_thar_ob     = thar_S_ob[:,2,:]
        lines_thar_co  = thar_S_co[:,1,:]
        iv_thar_co     = thar_S_co[:,2,:]

        All_Pixel_Centers = np.array([])
        All_Wavelengths   = np.array([])
        All_Orders        = np.array([])
        All_Centroids     = np.array([])
        All_Sigmas        = np.array([])
        All_Intensities   = np.array([])

        for order in range(n_useful):
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)
            
            thar_order_orig = lines_thar_ob[order,:]
            IV              = iv_thar_ob[order,:]
            wei             = np.sqrt( IV )
            bkg             = GLOBALutils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
            thar_order      = thar_order_orig - bkg

            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths,\
            rms_ms, residuals, centroids, sigmas, intensities \
                = GLOBALutils.Initial_Wav_Calibration(order_dir+'order_'+order_s+'o.iwdat',\
                                                      thar_order,order,wei,rmsmax=5000000,\
                                                      minlines=10,FixEnds=True,Dump_Argon=dumpargon,\
                                                      Dump_AllLines=True, Cheby=use_cheby)
            if (order == 35): 
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

        p0 = np.zeros( npar_wsol )
        p0[0] =  (35+89) * Global_ZP 
        p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
            GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders,\
                                                np.ones(All_Intensities.shape), p0, Cheby=use_cheby,\
                                                maxrms=MRMS, Inv=Inverse_m,minlines=700,order0=89, \
                                                ntotal=n_useful,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

        # Now calibrate COMPARISON orders. Use p1 above as p0
        All_Pixel_Centers_co = np.array([])
        All_Wavelengths_co   = np.array([])
        All_Orders_co        = np.array([])
        All_Centroids_co     = np.array([])
        All_Sigmas_co        = np.array([])
        All_Intensities_co   = np.array([])

        for order in range(22,n_useful):
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)
            
            thar_order_orig = lines_thar_co[order-22,:]
            IV              = iv_thar_co[order-22,:]
            wei             = np.sqrt( IV )
            bkg             = GLOBALutils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
            thar_order      = thar_order_orig - bkg
            
            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
                = GLOBALutils.Initial_Wav_Calibration( order_dir+'order_'+order_s+'o.iwdat', thar_order, order, wei, \
                                                            rmsmax=5000000, minlines=10,FixEnds=True,Dump_Argon=dumpargon, \
                                                            Dump_AllLines=True, Cheby=use_cheby)
          
            All_Pixel_Centers_co = np.append( All_Pixel_Centers_co, pixel_centers )
            All_Wavelengths_co   = np.append( All_Wavelengths_co, wavelengths )
            All_Orders_co        = np.append( All_Orders_co, np.zeros( len(pixel_centers) ) + order )
            All_Centroids_co     = np.append( All_Centroids_co, centroids)
            All_Sigmas_co        = np.append( All_Sigmas_co, sigmas)
            All_Intensities_co   = np.append( All_Intensities_co, intensities )

        p1_co, G_pix_co, G_ord_co, G_wav_co, II_co, rms_ms_co, G_res_co = \
            GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co, All_Wavelengths_co, All_Orders_co,\
                                                     np.ones(All_Intensities_co.shape), p1, Cheby=use_cheby,\
                                                     maxrms=MRMS, Inv=Inverse_m,minlines=500,order0=89,ntotal=n_useful,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
           
        # end COMPARISON orders. 
        pdict = {'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                     'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Orders':All_Orders, 'All_Sigmas':All_Sigmas,
                 'p1_co':p1_co, 'G_pix_co':G_pix_co, 'G_ord_co':G_ord_co, 'G_wav_co':G_wav_co, 'II_co':II_co, 'rms_ms_co':rms_ms_co,\
                     'G_res_co':G_res_co, 'All_Centroids_co':All_Centroids_co}
        pickle.dump( pdict, open( wavsol_pkl, 'w' ) )
        #print "Median sigma:", np.median( All_Sigmas )    
        p0_array[i,:] = p1

    else:
        print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl
        pdict         = pickle.load(open(wavsol_pkl,'r'))
        p0_array[i,:] = pdict['p1']

p0_G  = np.median(p0_array,axis=0)

if len(ThFP_ref) > 0:
    print '\n\tExtraction of Fabry-Perot calibration frames:'
else:
    print '\n\tNo Fabry-Perot calibration images found, moving on'
# Now extract ThAr-FP images
for fsim in ThFP_ref:
    hthfp = pyfits.open( fsim )
    thfp_fits_ob = dirout + fsim.split('/')[-1][:-8]+'spec.ob.fits.S'
    thfp_fits_co = dirout + fsim.split('/')[-1][:-8]+'spec.co.fits.S'
 
    if ( os.access(thfp_fits_ob,os.F_OK) == False ) or \
       ( os.access(thfp_fits_co,os.F_OK) == False ) or \
       (force_thfp_extract):
        print "\t\tNo previous extraction or extraction forced for ThFP file", fsim, "extracting..."
        dthfp = coralieutils.OverscanTrim( pyfits.getdata( fsim ) )
        Centers = np.zeros((len(c_all),dthfp.shape[1]))
        for i in range(c_all.shape[0]):
            Centers[i,:]=scipy.polyval(c_all[i,:],np.arange(len(Centers[i,:])))
        bac = GLOBALutils.get_scat(dthfp,Centers,span=5)
        dthfp -= bac

	thfp_S_ob  = GLOBALutils.optimal_extraction(dthfp,P_ob,c_ob,ext_aperture,\
                                                    hthfp[0].header['HIERARCH ESO CORA CCD RON'],\
                                                    hthfp[0].header['HIERARCH ESO CORA CCD GAIN'],\
                                                    S_Marsh,100.,min_extract_col,max_extract_col,npools)
	thfp_S_co  = GLOBALutils.optimal_extraction(dthfp,P_co,c_co,ext_aperture,\
                                                    hthfp[0].header['HIERARCH ESO CORA CCD RON'],
                                                    hthfp[0].header['HIERARCH ESO CORA CCD GAIN'],
                                                    S_Marsh,100.,min_extract_col,max_extract_col,npools)
        # save as fits file
        if (os.access(thfp_fits_ob,os.F_OK)):
            os.remove( thfp_fits_ob )
        if (os.access(thfp_fits_co,os.F_OK)):
            os.remove( thfp_fits_co )
            
        hdu = pyfits.PrimaryHDU( thfp_S_ob )
        hdu.writeto( thfp_fits_ob )
        hdu = pyfits.PrimaryHDU( thfp_S_co )
        hdu.writeto( thfp_fits_co )
    else:
        print "\t\tFP file", fsim, "all ready extracted, loading..."

# Now calibrate the ThFP spectra with the closest ThAr spectrum
print '\n\tWavelength solution of Fabry-Perot spectra with closest ThAr spectrum:'
for fsim in ThFP_ref:
    hthfp = pyfits.open( fsim )
    mjd, mjd0 = coralieutils.mjd_fromheader(hthfp)
    im = np.argmin(np.absolute(np.array(ThAr_ref_dates) - mjd))
    wavsol_dict = pickle.load(open(dirout + ThAr_ref[im].split('/')[-1][:-8]+'wavsolpars.pkl','r'))
    thfp_fits_ob = dirout + fsim.split('/')[-1][:-8]+'spec.ob.fits.S'
    thfp_fits_co = dirout + fsim.split('/')[-1][:-8]+'spec.co.fits.S'
    wavsol_pkl_fp = dirout + fsim.split('/')[-1][:-8]+'wavsolpars.pkl'
    if ( os.access(wavsol_pkl_fp,os.F_OK) == False ) or (force_thfp_wavcal):
            print '\t\tCalibrating', fsim,'...'
	    thar_fp = pyfits.getdata(thfp_fits_ob)
	    lines_thar_ob  = thar_fp[:,1,:]
	    iv_thar_ob     = thar_fp[:,2,:]

	    All_Pixel_Centers = np.array([])
	    All_Wavelengths   = np.array([])
	    All_Orders        = np.array([])
	    All_Centroids     = np.array([])
	    All_Sigmas        = np.array([])
	    All_Intensities   = np.array([])

	    for order in range(n_useful):
		order_s = str(order)
		if (order < 10):
		    order_s = '0'+str(order)
		thar_order_orig = lines_thar_ob[order,:]
		IV              = iv_thar_ob[order,:]
		wei             = np.sqrt( IV )
		bkg             = GLOBALutils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
		thar_order      = thar_order_orig - bkg

		coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
		            = GLOBALutils.Initial_Wav_Calibration( order_dir+'order_'+order_s+'o.iwdat', thar_order, order, wei, \
		                                                        rmsmax=5000000, minlines=10,FixEnds=True,Dump_Argon=dumpargon, \
		                                                        Dump_AllLines=True, Cheby=use_cheby)

		All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
		All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
		All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order )
		All_Centroids     = np.append( All_Centroids, centroids)
		All_Sigmas        = np.append( All_Sigmas, sigmas)
		All_Intensities   = np.append( All_Intensities, intensities )
		
	    p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
		GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders,\
		                                             np.ones(All_Intensities.shape), p0_G, Cheby=use_cheby,\
		                                             maxrms=100, Inv=Inverse_m, minlines=minlines_glob_ob,\
                                                             order0=89,ntotal=n_useful,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

	    p_shift, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(G_pix, G_wav, G_ord,\
		                                                   np.ones(G_wav.shape), wavsol_dict['p1'],\
		                                                   Cheby=True,Inv=True,maxrms=100,minlines=minlines_glob_ob,\
								   order0=89,ntotal=n_useful,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

	    pdict = {'p1':p1,'p_shift':p_shift,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
		'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Orders':All_Orders, 'All_Sigmas':All_Sigmas, 'p1_co':wavsol_dict['p1_co']}
	    pickle.dump( pdict, open( wavsol_pkl_fp, 'w' ) )
    else:
	print '\t\tFP spectrum', fsim, 'already calibrated, loading...'

### start of science frame reductions ###

new_list = []
new_list_obnames = []
new_list_texp = []
for i in range(len(simThAr_sci)):
    fsim = simThAr_sci[i]
    obname = obnames[i]
    texp = exptimes[i]
    if (object2do == 'all'):
        new_list.append(fsim)
        new_list_obnames.append( obname )
        new_list_texp.append( texp )

    else:
        if (obname == object2do):
            new_list.append(fsim)
            new_list_obnames.append( obname )
            new_list_texp.append( texp )

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

spec_moon = np.array(spec_moon)
use_moon  = np.array(use_moon)

# now extract the images
for nlisti in range(len(new_list)):	
    fsim   = new_list[ nlisti ]
    obname = new_list_obnames[ nlisti ]
    TEXP   = np.around(new_list_texp[ nlisti ])

    know_moon = False
    if fsim.split('/')[-1] in spec_moon:
        I = np.where(fsim.split('/')[-1] == spec_moon)[0]
        know_moon = True
        here_moon = use_moon[I]

    # get header h of image
    h = pyfits.open(fsim)
    print '\n'
    print "\t--> Working on image: ", fsim
    
    # get mjd and mjd0
    mjd,mjd0 = coralieutils.mjd_fromheader(h)

    #  get gain and readnoise of object 
    ronoise = h[0].header['HIERARCH ESO CORA CCD RON']
    gain    = h[0].header['HIERARCH ESO CORA CCD GAIN']
    
    # Object name

    print "\t\tObject name:",obname

    # Open file, trim, overscan subtract and MasterBias subtract
    data = h[0].data
    data = coralieutils.OverscanTrim(data)
    data -= MasterBias

    bacfile = dirout + 'BAC_' + fsim.split('/')[-1][:-4]+'fits'''
    if (os.access(bacfile,os.F_OK))== False:
        Centers = np.zeros((len(c_all),data.shape[1]))
        for i in range(c_all.shape[0]):
            Centers[i,:]=scipy.polyval(c_all[i,:],np.arange(len(Centers[i,:])))
        bac = GLOBALutils.get_scat(data,Centers,span=5)
        hdbac = pyfits.PrimaryHDU( bac )
        hdbac.writeto(bacfile)
    else:
        bac = pyfits.getdata(bacfile)
    data -= bac

    ra,dec = h[0].header['RA'],h[0].header['DEC']

    ra2,dec2 = GLOBALutils.getcoords(obname,mjd,filen=reffile)
    if ra2 !=0 and dec2 != 0:
        ra = ra2
        dec = dec2
    else:
        print '\t\tUsing the coordinates found in the image header.'

    # Find lambda_bary/lambda_topo using JPLEPHEM
    altitude    = h[0].header['HIERARCH ESO OBS GEO ALTITUDE']
    latitude    = h[0].header['HIERARCH ESO OBS GEO LATITU']
    longitude   = h[0].header['HIERARCH ESO OBS GEO LONGIT']
    epoch       = h[0].header['HIERARCH ESO OBS EQUICAT']

    iers                    = GLOBALutils.JPLiers( baryc_dir, mjd-999.0, mjd+999.0 )
    obsradius, R0           = GLOBALutils.JPLR0( latitude, altitude)
    obpos                   = GLOBALutils.obspos( longitude, obsradius, R0 )

    jplephem.set_ephemeris_dir( baryc_dir , ephemeris )
    jplephem.set_observer_coordinates( obpos[0], obpos[1], obpos[2] )

    res         = jplephem.doppler_fraction(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
    lbary_ltopo = 1.0 + res['frac'][0]
    bcvel_baryc = ( lbary_ltopo - 1.0 ) * 2.99792458E5

    print "\t\tBarycentric velocity:", bcvel_baryc

    res = jplephem.pulse_delay(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)   
    mbjd = mjd + res['delay'][0] / (3600.0 * 24.0)

    # Moon Phase Calculations
    gobs      = ephem.Observer()  
    gobs.name = 'Swiss1.2'  
    gobs.lat  = rad(latitude)  # lat/long in decimal degrees  
    gobs.long = rad(longitude)

    DDATE = h[0].header['HIERARCH ESO CORA SHUTTER START DATE']
    HHOUR = h[0].header['HIERARCH ESO CORA SHUTTER START HOUR']

    Mho = str(int(HHOUR))
    if len(Mho)<2:
        Mho = '0'+Mho
    mins = (HHOUR - int(Mho))*60.
    Mmi = str(int(mins))
    if len(Mmi)<2:
        Mmi = '0'+Mmi
    segs = (mins - int(Mmi))*60.
    if segs<10:
        Mse = '0'+str(segs)[:5]
    else:
        Mse = str(segs)[:6]

    gobs.date = str(DDATE[:4]) + '-' +  str(DDATE[4:6]) + '-' + str(DDATE[6:]) + ' ' +  Mho + ':' + Mmi +':' +Mse

    mephem = ephem.Moon()
    mephem.compute(gobs)
    Mcoo = jplephem.object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Mp   = jplephem.barycentric_object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Sp   = jplephem.barycentric_object_track("Sun", int(mjd), float(mjd%1), 1, 0.0)
    res  = jplephem.object_doppler("Moon", int(mjd), mjd%1, 1, 0.0)
    lunation,moon_state,moonsep,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,ra,dec)
    refvel = bcvel_baryc + moonvel
    print '\t\tRadial Velocity of sacttered moonlight:',refvel

    sorted_indices    = np.argsort( np.abs( np.array(ThAr_ref_dates) - mjd ) )
    sorted_indices_FP = np.argsort( np.abs( np.array(ThFP_ref_dates) - mjd ) )

    print '\t\tExtraction:'
    # optimally and simply extract spectra
    sci_fits_ob = dirout + fsim.split('/')[-1][:-8]+'spec.ob.fits.S'
    sci_fits_co = dirout + fsim.split('/')[-1][:-8]+'spec.co.fits.S'
    sci_fits_ob_simple = dirout + fsim.split('/')[-1][:-8]+'spec.simple.ob.fits.S'
    sci_fits_co_simple = dirout + fsim.split('/')[-1][:-8]+'spec.simple.co.fits.S'
    sci_fits_bac = dirout + fsim.split('/')[-1][:-8]+'spec.simple.bac.fits.S'
    if ( os.access(sci_fits_ob,os.F_OK) == False ) or \
       ( os.access(sci_fits_co,os.F_OK) == False ) or \
       ( os.access(sci_fits_ob_simple,os.F_OK) == False ) or \
       ( os.access(sci_fits_co_simple,os.F_OK) == False ) or \
       ( os.access(sci_fits_bac,os.F_OK) == False ) or \
       (force_sci_extract):

        print "\t\t\tNo previous extraction or extraction forced for science file", fsim, "extracting..."
        sci_Ss_ob = GLOBALutils.simple_extraction(data,c_ob,ext_aperture,\
                                                  min_extract_col,max_extract_col,npools)
        sci_Ss_co = GLOBALutils.simple_extraction(data,c_co,ext_aperture,\
                                                  min_extract_col,max_extract_col,npools)
        sci_S_ob  = GLOBALutils.optimal_extraction(data,P_ob,c_ob,ext_aperture,\
                                                   h[0].header['HIERARCH ESO CORA CCD RON'],\
                                                   h[0].header['HIERARCH ESO CORA CCD GAIN'],\
                                                   S_Marsh,NCosmic_Marsh,\
                                                   min_extract_col,max_extract_col,npools)
        sci_S_co  = GLOBALutils.optimal_extraction(data,P_co,c_co,ext_aperture,\
                                                   h[0].header['HIERARCH ESO CORA CCD RON'],\
                                                   h[0].header['HIERARCH ESO CORA CCD GAIN'],\
                                                   S_Marsh,2.*NCosmic_Marsh,\
                                                   min_extract_col,max_extract_col,npools)
        sci_bac =  GLOBALutils.simple_extraction(bac,c_ob,ext_aperture,\
                                                  min_extract_col,max_extract_col,npools)       

        # save as fits file
        if (os.access(sci_fits_ob,os.F_OK)):
            os.remove( sci_fits_ob )
        if (os.access(sci_fits_co,os.F_OK)):
            os.remove( sci_fits_co )
        if (os.access(sci_fits_ob_simple,os.F_OK)):
            os.remove( sci_fits_ob_simple )
        if (os.access(sci_fits_co_simple,os.F_OK)):
            os.remove( sci_fits_co_simple )
        if (os.access(sci_fits_bac,os.F_OK)):
            os.remove( sci_fits_bac )          
        hdu = pyfits.PrimaryHDU( sci_S_ob )
        hdu.writeto( sci_fits_ob )
        hdu = pyfits.PrimaryHDU( sci_S_co )
        hdu.writeto( sci_fits_co )
        hdu = pyfits.PrimaryHDU( sci_Ss_ob )
        hdu.writeto( sci_fits_ob_simple )
        hdu = pyfits.PrimaryHDU( sci_Ss_co )
        hdu.writeto( sci_fits_co_simple )
        hdu = pyfits.PrimaryHDU( sci_bac )
        hdu.writeto( sci_fits_bac )
    else:
        print '\t\t\t'+fsim, "has already been extracted, reading in product fits files..."
        sci_S_ob = pyfits.getdata( sci_fits_ob )
        sci_S_co = pyfits.getdata( sci_fits_co )
        sci_Ss_ob = pyfits.getdata( sci_fits_ob_simple )
        sci_Ss_co = pyfits.getdata( sci_fits_ob_simple )
        sci_bac = pyfits.getdata( sci_fits_bac )
    
    fout = 'proc/'+ obname + '_' + \
        h[0].header['HIERARCH ESO CORA SHUTTER START DATE'] + '_' +\
        'UT' + fsim[-17:-9] + '_' +\
        'sp.fits'

    #Build spectra
    if ( os.access(dirout+fout ,os.F_OK) == False ) or (force_spectral_file_build):
        # initialize file that will have the spectra
        # n_useful should be nord_ob, but we still have not calibrated that bluest order -- TODO
        spec = np.zeros((11, n_useful, data.shape[1]))
        hdu = pyfits.PrimaryHDU( spec )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MJD', mjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MBJD', mbjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START DATE', h[0].header['HIERARCH ESO CORA SHUTTER START DATE'] )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START UT',  h[0].header['HIERARCH ESO CORA SHUTTER START HOUR'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TEXP (S)',h[0].header['HIERARCH ESO OBS TEXP'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH FLUX WEIGHTED MEAN F ',h[0].header['HIERARCH ESO CORA PM FLUX TMMEAN'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (KM/S)', bcvel_baryc)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH (LAMBDA_BARY / LAMBDA_TOPO)', lbary_ltopo)    
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARGET NAME', obname)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA',h[0].header['HIERARCH ESO TEL TARG ALPHA'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC',h[0].header['HIERARCH ESO TEL TARG DELTA'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA DEG',h[0].header['RA'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC DEG',h[0].header['DEC'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA BARY',ra)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC BARY',dec)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH EQUINOX',h[0].header['HIERARCH ESO OBS EQUICAT'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LATITUDE',h[0].header['HIERARCH ESO OBS GEO LATITU'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LONGITUDE',h[0].header['HIERARCH ESO OBS GEO LONGIT'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS ALTITUDE',h[0].header['HIERARCH ESO OBS GEO ALTITUDE'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARG AIRMASS',h[0].header['HIERARCH ESO OBS TARG AIRMASS'])

    print '\t\tWavelength calibration:'
    print '\t\t\tComparision fibre is '+ h[0].header['HIERARCH ESO TPL TYPE']
    if h[0].header['HIERARCH ESO TPL TYPE'] == 'OBTH':
        # get ThAr closest in time
        indice = sorted_indices[0]
        thar_fits_ob = dirout + ThAr_ref[indice].split('/')[-1][:-8]+'spec.ob.fits.S'
        thar_fits_co = dirout + ThAr_ref[indice].split('/')[-1][:-8]+'spec.co.fits.S'
        pkl_wsol = dirout + ThAr_ref[indice].split('/')[-1][:-8]+'wavsolpars.pkl'
        print "\t\t\tUnpickling reference wavelength solution from", pkl_wsol, " ..."
        wsol_dict = pickle.load(open(pkl_wsol,'r'))
        # Extract thAr lines from comparison orders
        lines_thar_co  = sci_S_co[:,1,:]
        iv_thar_co     = sci_S_co[:,2,:]
        All_Pixel_Centers_co = np.array([])
        All_Wavelengths_co   = np.array([])
        All_Orders_co        = np.array([])
        All_Centroids_co     = np.array([])
        All_Sigmas_co        = np.array([])
        All_Intensities_co   = np.array([])
        for order in range(22,n_useful):
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)
            thar_order_orig = lines_thar_co[order-22,:]
            IV              = iv_thar_co[order-22,:]
            wei             = np.sqrt( IV )
            bkg             = GLOBALutils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
            thar_order      = thar_order_orig - bkg
		    
            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
		        = GLOBALutils.Initial_Wav_Calibration( order_dir+'order_'+order_s+'o.iwdat', thar_order, order, wei, \
		                                                    rmsmax=5000000, minlines=10,FixEnds=True,Dump_Argon=dumpargon, \
		                                                    Dump_AllLines=True, Cheby=use_cheby)
		
            All_Pixel_Centers_co = np.append( All_Pixel_Centers_co, pixel_centers )
            All_Wavelengths_co   = np.append( All_Wavelengths_co, wavelengths )
            All_Orders_co        = np.append( All_Orders_co, np.zeros( len(pixel_centers) ) + order )
            All_Centroids_co     = np.append( All_Centroids_co, centroids)
            All_Sigmas_co        = np.append( All_Sigmas_co, sigmas)
            All_Intensities_co   = np.append( All_Intensities_co, intensities )
	    
        # get a global solution for the lines found
        p1_co, G_pix_co, G_ord_co, G_wav_co, II_co, rms_ms_co, G_res_co = \
		    GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co, All_Wavelengths_co, All_Orders_co,\
		                                             np.ones(All_Intensities_co.shape), wsol_dict['p1_co'], Cheby=use_cheby,\
		                                             maxrms=MRMS, Inv=Inverse_m,minlines=minlines_glob_co,\
                                                             order0=89,ntotal=n_useful,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
	    
        # get shift with respect to reference ThAr
        p_shift, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(G_pix_co, G_wav_co, G_ord_co,\
		                                                   np.ones(G_wav_co.shape), wsol_dict['p1_co'],\
		                                                   Cheby=True,Inv=True,maxrms=100,minlines=minlines_glob_co,\
                                                                   order0=89,ntotal=n_useful,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
	    
        precision = rms_ms/np.sqrt(len(I))
        good_quality = True
        if (precision > 10):
            good_quality = False
		
        # Apply new wavelength solution including barycentric correction

    else:
        indice = sorted_indices_FP[0]
        thfp_fits_co = dirout + ThFP_ref[indice].split('/')[-1][:-8]+'spec.co.fits.S'
        pkl_wsol = dirout + ThFP_ref[indice].split('/')[-1][:-8]+'wavsolpars.pkl'
        print "\t\t\tUnpickling reference wavelength solution from", pkl_wsol, " ..."
        wsol_dict = pickle.load(open(pkl_wsol,'r'))
        lines_thar_co = np.zeros(sci_Ss_co.shape)
        lines_thar_co_ref = np.zeros(sci_Ss_co.shape)
        for si in range(S_flat_co_n.shape[0]):
            JI = np.where(S_flat_co_n[si,1]>0)[0]
            lines_thar_co[si,JI] = sci_S_co[si,1,JI] / S_flat_co_n[si,1,JI]
            lines_thar_co_ref[si,JI]  = pyfits.getdata(thfp_fits_co)[si,1,JI] / S_flat_co_n[si,1,JI]
            JI1 = np.where(lines_thar_co[si]<0)[0]
            JI2 = np.where(lines_thar_co_ref[si]<0)[0]
            lines_thar_co[si,JI1] = 0.
            lines_thar_co_ref[si,JI2] = 0.
        #lines_thar_co      = sci_S_co[:,1,:] / S_flat_co_simple_n
        #lines_thar_co_ref  = pyfits.getdata(thfp_fits_co)[:,1,:] / S_flat_co_simple_n
        rv_fps = []
        for order in range(nord_co):
            I = np.where(np.isnan(lines_thar_co[order]))[0]
            lines_thar_co[order][I]=0.
            I = np.where(np.isnan(lines_thar_co_ref[order]))[0]
            lines_thar_co_ref[order][I]=0.
            try:
                tc = GLOBALutils.fp_base(lines_thar_co[order])
                tcr = GLOBALutils.fp_base(lines_thar_co_ref[order])
                IJ1 = np.where(tc!=0)[0]
                IJ2 = np.where(tcr!=0)[0]
                tc /= np.median(tc[IJ1])
                tcr /= np.median(tcr[IJ2])
                rv_fp = GLOBALutils.ccf_fp(tc,tcr,wsol_dict['p1_co'],order+22,order0=89,nx=ncoef_x,nm=ncoef_m,npix=len(tc))
            except:
                rv_fp = -999
            rv_fps.append(rv_fp)
        #plot(rv_fps,'ro')

        rv_fps = np.array(rv_fps)
        I = np.where(rv_fps!=-999)[0]
        rv_fps = rv_fps[I]
        rv_fps = GLOBALutils.sig_cli2(rv_fps,ns=3.)
        #plot(rv_fps,'ro')
        #show()
        #print np.median(rv_fps),np.sqrt(np.var(rv_fps))/np.sqrt(float(len(rv_fps)))
        fp_shift = np.median(rv_fps)
        p_sh = wsol_dict['p_shift'] * 299792458. * 1e-6
        fp_shift += p_sh
        p_shift = 1e6*fp_shift/299792458.
        print '\t\t\tFP shift = ',fp_shift[0],'+-',np.sqrt(np.var(rv_fps))/np.sqrt(float(len(rv_fps))),'m/s'
        good_quality = True

        equis = np.arange( data.shape[1] )        

        for order in range(n_useful):
            m = order + 89
            chebs = GLOBALutils.Calculate_chebs(equis, m, Inverse=Inverse_m,order0=89,ntotal=n_useful,npix=data.shape[1],nx=ncoef_x,nm=ncoef_m)
            if good_quality:
                WavSol = lbary_ltopo * (1.0 + 1.0e-6*p_shift) * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wsol_dict['p1'],chebs,ncoef_x,ncoef_m)
            else:
                WavSol = lbary_ltopo * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wsol_dict['p1'],chebs,ncoef_x,ncoef_m)
            spec[0,order,:] = WavSol
            spec[1,order,:] = sci_S_ob[order,1, :]
            spec[2,order,:] = sci_S_ob[order,2, :]
            # Flat-fielded spectrum
            fn = S_flat_ob_n[order,1,:]
            L  = np.where( fn > 0 )
            spec[3,order,:][L] = sci_S_ob[order,1,:][L] / S_flat_ob_n[order,1,:][L]
            spec[4,order,:][L] = sci_S_ob[order,2,:][L] * ( S_flat_ob_n[order,1,:][L] ** 2 )

        # Continuum normalized spectrum
        ron  = h[0].header['HIERARCH ESO CORA CCD RON']
        gain = h[0].header['HIERARCH ESO CORA CCD GAIN']

        wav_temp, norm_spec = continuum.NORM2( spec[0,:,:],spec[3,:,:])
        for order in range(n_useful):
            L  = np.where( spec[1,order,:] != 0 )
            spec[5,order,:][L] = norm_spec[order][L]
            nJ = np.where(np.isnan(spec[5,order])==True)[0]
            nJ2 = np.where(np.isinf(spec[5,order])==True)[0]
            spec[5,order,nJ] = 1.0
            spec[5,order,nJ2] = 1.0
            ratio              = spec[3,order,:][L] / norm_spec[order][L]  
            spec[6,order,:][L] = spec[4,order,:][L] * (ratio ** 2 )
            spec[7,order,:][L] = ratio
            #spec[8,order,:][L] = ratio * S_flat_ob_n[order,1,:][L] / np.sqrt( ratio * S_flat_ob_n[order,1,:][L] / gain + ext_aperture*2*(ron/gain)**2 + sci_bac[order,:][L] / gain )
            spec[8,order,:][L] = ratio * S_flat_ob_n[order,1,:][L] / np.sqrt( ratio * S_flat_ob_n[order,1,:][L] / gain + (ron/gain)**2 )

            spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
            dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
            NN            = np.average(dlambda_dx)
            dlambda_dx    /= NN
            # clean-up of CRs in continuum-normalized spectrum. Set troublesome pixels to 1
            median_cn_spectrum = np.zeros( spec[5,order,:].shape )
            median_cn_spectrum[L] = scipy.signal.medfilt( spec[5,order,:][L], 7 )
            LK = np.where(spec[8,order] == 0.)[0]
            spec[8,order,LK] = 0.000001
            LL =  np.where(spec[5,order] > 1. + 5./spec[8,order])
            LL2 = np.where(spec[5,order] < - 5./spec[8,order])
            spec[8,order,LK] = 0.
            spec[5,order,:][LL] = 1
            spec[5,order,:][LL2] = 1
            spec[5,order,:][LK] = 0
            spec[6,order,:][LL] = spec[8,order,:][LL] ** 2
            spec[6,order,:][LL2] = spec[8,order,:][LL] ** 2

            spec[9,order,:][L] = spec[5,order,:][L] * (dlambda_dx[L] ** 1) 
            spec[10,order,:][L] = spec[6,order,:][L] / (dlambda_dx[L] ** 2)

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
            #    query_success,sp_type_query = GLOBALutils.simbad_query_coords('12:00:00','00:00:00')
            #print "\t\t\tSpectral type returned by SIMBAD query:",sp_type_query

            hdu = GLOBALutils.update_header(hdu,'HIERARCH SIMBAD SPTYP', sp_type_query)
            pars_file = dirout + fsim.split('/')[-1][:-8]+'_stellar_pars.txt'

            if os.access(pars_file,os.F_OK) == False or force_stellar_pars:
                print "\t\t\tEstimating atmospheric parameters:"
                Rx = np.around(1./np.sqrt(1./40000.**2 - 1./60000**2))
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

            moonsep_cor = h[0].header['HIERARCH ESO OBS MOON SEP']
		
            if sp_type == 'M5':
                moon_sig = 2.5
            elif sp_type == 'K5':
                moon_sig = 3.3
            else:
                moon_sig = 4.5

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
		         'lunation':lunation,'mephem':mephem,'texp':h[0].header['EXPTIME']}

        pkl_xc = dirout + fsim.split('/')[-1][:-8]+obname+'_XC_'+sp_type+'.pkl'
        pickle.dump( xc_dict, open( pkl_xc, 'w' ) )

        ccf_pdf = dirout + 'proc/' + fsim.split('/')[-1][:-4] + obname + '_XCs_' + sp_type + '.pdf'

        if not avoid_plot:
            GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)

        SNR_5130 = np.median(spec[8,30,1000:1101] )
        airmass  = h[0].header['HIERARCH ESO OBS TARG AIRMASS']
        seeing   = h[0].header['HIERARCH ESO OBS AMBI DIMM SEEING']
        
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

        if not good_quality:
            RVerr2 = np.sqrt(0.03**2 + RVerr2**2)

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
        hdu = GLOBALutils.update_header(hdu,'INST', 'CORALIE')
        hdu = GLOBALutils.update_header(hdu,'RESOL', '60000')
        hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')
        hdu = GLOBALutils.update_header(hdu,'XC_MIN', XC_min)
        hdu = GLOBALutils.update_header(hdu,'BJD_OUT', bjd_out)
	    
        line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f   coralie   ceres   60000 %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
                      (obname, bjd_out, RV, RVerr2, BS, BSerr, T_eff_epoch, logg_epoch, Z_epoch, vsini_epoch, XC_min, disp_epoch,\
        TEXP, SNR_5130_R, ccf_pdf)
        f_res.write(line_out)
        if (os.access( dirout + fout,os.F_OK)):
            os.remove( dirout + fout)
        hdu.writeto( dirout + fout )

    else:
        print "\t\tReading spectral file from", fout
        spec = pyfits.getdata( fout )

f_res.close()
