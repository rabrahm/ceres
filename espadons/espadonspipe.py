import sys
import matplotlib
matplotlib.use("Agg") 

from pylab import *

base = '../'
sys.path.append(base+"utils/Correlation")
sys.path.append(base+"utils/GLOBALutils")

baryc_dir= base+'utils/SSEphem/'
sys.path.append(baryc_dir)
ephemeris='DEc403'

import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ceres modules
import espadonsutils
import correlation
import GLOBALutils

# other useful modules
import argparse
import ephem

import jplephem
#from PyAstronomy import pyasl
from math import radians as rad
from astropy.io import fits as pyfits
import pickle
import os
import scipy
import scipy.interpolate
from scipy import interpolate, ndimage

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
parser.add_argument('-mode',default='ss')
parser.add_argument('-amps',default='a')

args   = parser.parse_args()
dirin            = args.directorio
avoid_plot       = args.avoid_plot
dirout           = args.dirout
DoClass          = args.do_class
JustExtract      = args.just_extract
npools           = int(args.npools)
object2do        = args.o2do
reffile          = args.reffile
mode             = args.mode
amps             = args.amps

if dirin[-1] != '/':
    dirin = dirin + '/'

if dirout == 'default':
    dirout = dirin[:-1]+'_red'+mode+'_'+amps+'/'

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
force_thar_extract = False
force_thar_wavcal  = False
force_tharxc       = False
force_sci_extract  = False
force_stellar_pars = False
force_trace    = False
dumpargon          = False
force_spectral_file_build = True

bad_colummn        = True
Inverse_m          = True
use_cheby          = True

MRMS               = 60   # max rms in m/s, global wav solution

trace_degree       = 4
Marsh_alg          = 0
ext_aperture       = 15
NSigma_Marsh       = 10
NCosmic_Marsh      = 10
S_Marsh            = 0.4
N_Marsh            = 4      # grado polinomio 
min_extract_col    = 50
max_extract_col    = 4600

porder_single      = 5
ncoef_x            = 5
ncoef_m            = 7
npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2
nconts = np.array([1,1,1,1,1,1,2,2,2,1,3,3,3,3,2,2,3,3,3,4,3,3,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
lim_iz = np.array([100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,150,200,250,300,350,400,450,500,500,800])

models_path = base+"data/COELHO_MODELS/R_40000b/"
order_dir   = base+"espadons/wavcals/"
 
OO0 = 22

RES = 80000.
if mode != 'so':
    ext_aperture = 7
    RES = 65000.

# file containing the log
log = dirout+'night.log'

print "\n\n\tESPaDOnS CFHT3.6m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'

biases, flats, ThAr_ref, sim_sci, ThAr_ref_dates = espadonsutils.FileClassify(dirin,log,mode=mode,amps=amps)
ThAr_ref = ThAr_ref[:2]
ThAr_ref_dates = ThAr_ref_dates[:2]
if ( (os.access(dirout+'Flat_'+mode+'.fits',os.F_OK) == False)        or \
    (os.access(dirout+'trace_'+mode+'.pkl',os.F_OK) == False)        or \
    (os.access(dirout+'MasterBias_'+mode+'.fits',os.F_OK) == False)  or \
    (force_pre_process) ):
    print "\tNo previous pre-processing files or found"
    pre_process = 1
else:
    print "\tPre-processing files found, going straight to extraction"
    pre_process = 0

if (pre_process == 1):
    print "\t\tGenerating Master calibration frames..."
    # median combine Biases
    MasterBias, RO_bias, GA_bias = espadonsutils.MedianCombine(biases,zero_bo=False)
    hdu = pyfits.PrimaryHDU( MasterBias )
    if (os.access(dirout+'MasterBias_'+mode+'.fits',os.F_OK)):
        os.remove(dirout+'MasterBias_'+mode+'.fits')
    hdu.writeto(dirout+'MasterBias_'+mode+'.fits')
    print "\t\t-> Masterbias: done!"

    # median combine list of flats  
    Flat, RO_flat, GA_flat = espadonsutils.MedianCombine(flats, zero_bo=True, zero=dirout+'MasterBias_'+mode+'.fits')
    hdu = pyfits.PrimaryHDU( Flat )
    if (os.access(dirout+'Flat_'+mode+'.fits',os.F_OK)):
        os.remove(dirout+'Flat_'+mode+'.fits')
    hdu.writeto(dirout+'Flat_'+mode+'.fits')
    print "\t\t-> Masterflat: done!"
    GFlat = Flat.copy()
    if mode == 'so':
        d1 = np.loadtxt('mid_co_so.txt')
    else:
        d1 = np.loadtxt('mid_co_ss.txt')
    d2 = Flat[np.around(.5*Flat.shape[0]).astype('int')]
    d1 = d1/np.median(d1)
    d2 = d2/np.median(d2)
    ejx = np.arange(Flat.shape[1])-int(np.around(.5*Flat.shape[1]))
    corr = scipy.signal.correlate(d1, d2, mode='same')
    Imax = np.argmax(corr)
    traces_shift = ejx[Imax]

    for i in range(Flat.shape[0]):
        GFlat[i] = scipy.ndimage.filters.gaussian_filter(Flat[i],2)

    GFlat = GFlat.T
    print "\tTracing echelle orders..."
    if mode == 'so':
        c_all = espadonsutils.get_them(GFlat, 14, trace_degree,mode=0, shift=traces_shift)
        nord_all = len(c_all)
        print "\t\t"+str(nord_all)+" orders traced..."
        trace_dict = {'c_all':c_all, \
                  'nord_all':nord_all,\
                  'GA_flat':GA_flat,'RO_flat':RO_flat}

    else:
        c_ob, c_co = espadonsutils.get_them(GFlat, 7, trace_degree, mode=1, shift=traces_shift)
        nord_ob, nord_co = len(c_ob), len(c_co)
        c_all = GLOBALutils.Mesh(c_ob,c_co)
        nord_all = len(c_all)
        print "\t\t"+str(nord_ob)+" object orders traced..."
        print "\t\t"+str(nord_co)+" comparison orders traced..."
        trace_dict = {'c_all':c_all, 'c_ob':c_ob, 'c_co':c_co, \
                  'nord_all':nord_all, 'nord_ob':nord_ob, 'nord_co':nord_co, \
                  'GA_flat':GA_flat,'RO_flat':RO_flat}
    Flat = Flat.T
    pickle.dump( trace_dict, open( dirout+'trace_'+mode+'.pkl', 'w' ) )

else:
    trace_dict = pickle.load( open( dirout+'trace_'+mode+'.pkl', 'r' ) )
    c_all = trace_dict['c_all']
    nord_all = trace_dict['nord_all']
    # recover GA*, RO*
    GA_flat = trace_dict['GA_flat']
    RO_flat = trace_dict['RO_flat']
    # recover flats & master bias
    h = pyfits.open(dirout+'Flat_'+mode+'.fits')
    Flat = h[0].data
    Flat = Flat.T
    h = pyfits.open(dirout+'MasterBias_'+mode+'.fits')
    MasterBias = h[0].data
    if mode != 'so':
        nord_co = trace_dict['nord_co']
        nord_ob = trace_dict['nord_ob']
        c_ob, c_co = trace_dict['c_ob'], trace_dict['c_co']

"""
if mode == 'so':
    imshow(Flat,vmax=3000)
    for i in range(len(c_all)):
        plot(np.polyval(c_all[i],np.arange(4600)),'k')
    show()
else:
    imshow(Flat,vmax=3000)
    for i in range(len(c_ob)):
        plot(np.polyval(c_ob[i],np.arange(4600)),'k')
        plot(np.polyval(c_co[i],np.arange(4600)),'g')

    show()
"""

print '\n\tExtraction of Flat calibration frames:'
Flat = Flat.T
P_fits    = dirout + 'P.fits'
P_ob_fits    = dirout + 'P_ob.fits'
P_co_fits    = dirout + 'P_co.fits'

S_flat_fits        = dirout +'S_flat.fits'
S_flat_ob_fits        = dirout +'S_flat_ob.fits'
S_flat_co_fits        = dirout +'S_flat_co.fits'

do_flat = False
if mode == 'so' and (os.access(P_fits,os.F_OK) == False or os.access(S_flat_fits,os.F_OK) == False ):
    do_flat = True
if mode != 'so' and (os.access(P_ob_fits,os.F_OK) == False or os.access(S_flat_ob_fits,os.F_OK) == False or \
    os.access(P_co_fits,os.F_OK) == False or os.access(S_flat_co_fits,os.F_OK)==False):
    do_flat = True

if do_flat or (force_flat_extract):
    "\t\t\tComputing Background..."
    if mode == 'so':
        Centers = np.zeros((len(c_all),Flat.shape[0]))
        for i in range(c_all.shape[0]):
            Centers[i,:]=scipy.polyval(c_all[i],np.arange(len(Centers[i,:])))
    else:
        Centers = np.zeros((len(c_ob),Flat.shape[0]))
        for i in range(c_ob.shape[0]):
            Centers[i,:]=0.5*( scipy.polyval(c_ob[i],np.arange(len(Centers[i,:]))) + \
                              scipy.polyval(c_co[i],np.arange(len(Centers[i,:]))) )

    bac = GLOBALutils.get_scat(Flat.T,Centers,span=15).T
    #plot(Flat[2000])
    #plot(bac[2000])
    #show()
    #print gfd
    Flat -= bac

    if mode == 'so':
        P = GLOBALutils.obtain_P(Flat.T,c_all,ext_aperture,RO_flat, GA_flat, 100*NSigma_Marsh, S_Marsh, \
                    N_Marsh, Marsh_alg, min_extract_col,max_extract_col, npools)
        if ( os.access(P_fits,os.F_OK) ):
            os.remove(P_fits)
        hdu = pyfits.PrimaryHDU( P )
        hdu.writeto( P_fits )

    else:
        P_ob = GLOBALutils.obtain_P(Flat.T,c_ob,ext_aperture,RO_flat, GA_flat, 100*NSigma_Marsh, S_Marsh, \
                    N_Marsh, Marsh_alg, min_extract_col,max_extract_col, npools)
        P_co = GLOBALutils.obtain_P(Flat.T,c_co,ext_aperture,RO_flat, GA_flat, 100*NSigma_Marsh, S_Marsh, \
                    N_Marsh, Marsh_alg, min_extract_col,max_extract_col, npools)
        if ( os.access(P_ob_fits,os.F_OK) ):
            os.remove(P_ob_fits)
        hdu = pyfits.PrimaryHDU( P_ob )
        hdu.writeto( P_ob_fits )
        if ( os.access(P_co_fits,os.F_OK) ):
            os.remove(P_co_fits)
        hdu = pyfits.PrimaryHDU( P_co )
        hdu.writeto( P_co_fits )
        
    print "\t\t\tNo extracted flat spectra found or extraction forced, extracting and saving..."
    if mode == 'so':
        S_flat  = GLOBALutils.optimal_extraction(Flat.T,P,c_all,ext_aperture,\
                                                RO_flat,GA_flat,S_Marsh,NCosmic_Marsh,\
                                                min_extract_col,max_extract_col,npools)
        if (os.access(S_flat_fits,os.F_OK)):
            os.remove( S_flat_fits )              
        hdu = pyfits.PrimaryHDU( S_flat )
        hdu.writeto( S_flat_fits )
    else:
        S_flat_ob  = GLOBALutils.optimal_extraction(Flat.T,P_ob,c_ob,ext_aperture,\
                                                RO_flat,GA_flat,S_Marsh,NCosmic_Marsh,\
                                                min_extract_col,max_extract_col,npools)
        S_flat_co  = GLOBALutils.optimal_extraction(Flat.T,P_co,c_co,ext_aperture,\
                                                RO_flat,GA_flat,S_Marsh,NCosmic_Marsh,\
                                                min_extract_col,max_extract_col,npools)
        # write P_on and S_flat_ob as fits files
        if (os.access(S_flat_ob_fits,os.F_OK)):
            os.remove( S_flat_ob_fits )              
        hdu = pyfits.PrimaryHDU( S_flat_ob )
        hdu.writeto( S_flat_ob_fits )
        if (os.access(S_flat_co_fits,os.F_OK)):
            os.remove( S_flat_co_fits )              
        hdu = pyfits.PrimaryHDU( S_flat_co )
        hdu.writeto( S_flat_co_fits )
           
else:
    print "\t\tExtracted flat spectra found, loading..."
    if mode == 'so':
        P            = pyfits.getdata( P_fits )
        S_flat        = pyfits.getdata( S_flat_fits )
    else:
        P_ob             = pyfits.getdata( P_ob_fits )
        S_flat_ob        = pyfits.getdata( S_flat_ob_fits )        
        P_co             = pyfits.getdata( P_co_fits )
        S_flat_co        = pyfits.getdata( S_flat_co_fits ) 

if mode == 'so':
    S_flat_n, Snorms = GLOBALutils.FlatNormalize_single( S_flat, mid=int(.5*S_flat.shape[2]))
else:
    S_flat_ob_n, Snorms_ob = GLOBALutils.FlatNormalize_single( S_flat_ob, mid=int(.5*S_flat_ob.shape[2]))
    S_flat_co_n, Snorms_co = GLOBALutils.FlatNormalize_single( S_flat_co, mid=int(.5*S_flat_co.shape[2]))

print '\n\tExtraction of ThAr calibration frames:'
# Extract all ThAr+Ne files
for fsim in ThAr_ref:
    print "\t\tWorking on ThAr+Ne file ", fsim, "..."
    hthar = pyfits.open( fsim )

    dthar = espadonsutils.OverscanTrim(pyfits.getdata(fsim))
    dthar = dthar - MasterBias

    thar_fits = dirout + fsim.split('/')[-1][:-4]+'spec.fits.S'
    thar_fits_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.fits.S'
    thar_ob_fits = dirout + fsim.split('/')[-1][:-4]+'spec.ob.fits.S'
    thar_ob_fits_simple = dirout + fsim.split('/')[-1][:-4]+'spec.ob.simple.fits.S'
    thar_co_fits = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
    thar_co_fits_simple = dirout + fsim.split('/')[-1][:-4]+'spec.co.simple.fits.S'
    do_thar = False
    if mode == 'so' and (os.access(thar_fits,os.F_OK) == False or os.access(thar_fits_simple,os.F_OK) == False): 
        do_thar = True
    if mode != 'so' and (os.access(thar_ob_fits,os.F_OK) == False or os.access(thar_ob_fits_simple,os.F_OK) == False or\
        os.access(thar_co_fits,os.F_OK) == False or os.access(thar_co_fits_simple,os.F_OK) == False ):
        do_thar = True
    if do_thar or force_thar_extract:
        if mode == 'so':
            Centers = np.zeros((len(c_all),dthar.shape[0]))
            for i in range(c_all.shape[0]):
                Centers[i,:]=scipy.polyval(c_all[i],np.arange(len(Centers[i,:])))
        else:
            Centers = np.zeros((len(c_ob),dthar.shape[0]))
            for i in range(c_ob.shape[0]):
                Centers[i,:]=0.5*( scipy.polyval(c_ob[i],np.arange(len(Centers[i,:]))) + \
                              scipy.polyval(c_co[i],np.arange(len(Centers[i,:]))) )
        bac = GLOBALutils.get_scat(dthar.T,Centers,span=15).T

        dthar -= bac
        dthar = dthar.T

        RO_thar, GA_thar = hthar[0].header['RDNOISEA'],hthar[0].header['GAINA']
  
        print "\t\tNo previous extraction or extraction forced for ThAr file", fsim, "extracting..."

        if mode == 'so':
            thar_Ss = GLOBALutils.simple_extraction(dthar,c_all,ext_aperture,min_extract_col,max_extract_col,npools)
            thar_S  = GLOBALutils.optimal_extraction(dthar,P,c_all,ext_aperture,RO_thar, GA_thar,S_Marsh,100.,min_extract_col,max_extract_col,npools)
            if (os.access(thar_fits,os.F_OK)):
                os.remove( thar_fits )
            if (os.access(thar_fits_simple,os.F_OK)):
                os.remove( thar_fits_simple )
            
            hdu = pyfits.PrimaryHDU( thar_S )
            hdu.writeto( thar_fits )
            hdu = pyfits.PrimaryHDU( thar_Ss )
            hdu.writeto( thar_fits_simple )
        else:
            thar_Ss_ob = GLOBALutils.simple_extraction(dthar,c_ob,ext_aperture,min_extract_col,max_extract_col,npools)
            thar_S_ob  = GLOBALutils.optimal_extraction(dthar,P_ob,c_ob,ext_aperture,RO_thar, GA_thar,S_Marsh,100.,min_extract_col,max_extract_col,npools)
            if (os.access(thar_ob_fits,os.F_OK)):
                os.remove( thar_ob_fits )
            if (os.access(thar_ob_fits_simple,os.F_OK)):
                os.remove( thar_ob_fits_simple )
            hdu = pyfits.PrimaryHDU( thar_S_ob )
            hdu.writeto( thar_ob_fits )
            hdu = pyfits.PrimaryHDU( thar_Ss_ob )
            hdu.writeto( thar_ob_fits_simple )

            thar_Ss_co = GLOBALutils.simple_extraction(dthar,c_co,ext_aperture,min_extract_col,max_extract_col,npools)
            thar_S_co  = GLOBALutils.optimal_extraction(dthar,P_co,c_co,ext_aperture,RO_thar, GA_thar,S_Marsh,100.,min_extract_col,max_extract_col,npools)
            if (os.access(thar_co_fits,os.F_OK)):
                os.remove( thar_co_fits )
            if (os.access(thar_co_fits_simple,os.F_OK)):
                os.remove( thar_co_fits_simple )
            hdu = pyfits.PrimaryHDU( thar_S_co )
            hdu.writeto( thar_co_fits )
            hdu = pyfits.PrimaryHDU( thar_Ss_co )
            hdu.writeto( thar_co_fits_simple )
    else:
        print "\t\tThAr file", fsim, "all ready extracted, loading..."

sorted_ThAr_dates = np.argsort( ThAr_ref_dates )
print "\n\tWavelength solution of ThAr calibration spectra:"

if mode == 'so':
    for i in range(len(sorted_ThAr_dates)):
        index      = sorted_ThAr_dates[i]  
        hd         = pyfits.getheader(ThAr_ref[index])
        wavsol_pkl = dirout + ThAr_ref[index].split('/')[-1][:-4]+'wavsolpars.pkl'
        
        if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
            print "\t\tComputing wavelength solution of ThAr file", ThAr_ref[index] 

            hthar        = pyfits.open( ThAr_ref[index] )
            mjd, mjd0    = espadonsutils.mjd_fromheader( hthar )

            thar_fits = dirout + ThAr_ref[index].split('/')[-1][:-4]+'spec.fits.S'
            thar_S    = pyfits.getdata( thar_fits )

            lines_thar  = thar_S[:,1,:]
            iv_thar     = thar_S[:,2,:]
        
            All_Pixel_Centers = np.array([])
            All_Wavelengths   = np.array([])
            All_Orders        = np.array([])
            All_Centroids     = np.array([])
            All_Sigmas        = np.array([])
            All_Intensities   = np.array([])
            All_residuals     = np.array([])

            order = 0
            while order < nord_all:
                order_s = str(order)
                if (order < 10):
                    order_s = '0'+str(order)
                
                thar_order_orig = lines_thar[order,:]
                IV              = iv_thar[order,:]
                wei             = np.sqrt( IV )
                #bkg             = ferosutils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)
                thar_order      = thar_order_orig #- bkg

                coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths,\
                    rms_ms, residuals, centroids, sigmas, intensities =\
                        GLOBALutils.Initial_Wav_Calibration(order_dir+'order_'+\
                        order_s+'.iwdat', thar_order, order, wei, rmsmax=100, \
                        minlines=10,FixEnds=False,Dump_Argon=dumpargon,\
                        Dump_AllLines=True, Cheby=use_cheby,porder=5,del_width=4.0)
                #print order
                #plot(pixel_centers, wavelengths,'ro')
                #show()
                if (order == 20): 
                    if (use_cheby):
                        Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, int(.5*len(thar_order)), len(thar_order) )
                    else:
                        Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )
                if order != 0:
                    All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
                    All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
                    All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order )
                    All_Centroids     = np.append( All_Centroids, centroids)
                    All_Sigmas        = np.append( All_Sigmas, sigmas)
                    All_Intensities   = np.append( All_Intensities, intensities )
                    All_residuals     = np.append( All_residuals, residuals )
                order+=1

            p0    = np.zeros( npar_wsol )
            p0[0] =  (20+OO0) * Global_ZP
            p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
                GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders,\
                                                    np.ones(All_Intensities.shape), p0, Cheby=use_cheby,\
                                                    maxrms=MRMS, Inv=Inverse_m,minlines=1200,order0=OO0, \
                                                    ntotal=nord_all,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
            #for i in np.unique(G_ord):
            #    I = np.where(G_ord == i)[0]
            #    plot(G_wav[I],G_res[I],'.')
            #show()

            pdict = {'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                         'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Wavelengths':All_Wavelengths, 'All_Orders':All_Orders, 'All_Pixel_Centers':All_Pixel_Centers, 'All_Sigmas':All_Sigmas}
            pickle.dump( pdict, open( wavsol_pkl, 'w' ) )

        else:
            print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl

else:
    for i in range(len(sorted_ThAr_dates)):
        index      = sorted_ThAr_dates[i]  
        hd         = pyfits.getheader(ThAr_ref[index])
        wavsol_pkl = dirout + ThAr_ref[index].split('/')[-1][:-4]+'wavsolpars.pkl'
        
        if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
            print "\t\tComputing wavelength solution of ThAr file", ThAr_ref[index] 

            hthar        = pyfits.open( ThAr_ref[index] )
            mjd, mjd0    = espadonsutils.mjd_fromheader( hthar )

            thar_fits_ob = dirout + ThAr_ref[index].split('/')[-1][:-4]+'spec.ob.fits.S'
            thar_fits_co = dirout + ThAr_ref[index].split('/')[-1][:-4]+'spec.co.fits.S'

            thar_S_ob    = pyfits.getdata( thar_fits_ob )
            thar_S_co    = pyfits.getdata( thar_fits_co )

            lines_thar_ob  = thar_S_ob[:,1,:]
            iv_thar_ob     = thar_S_ob[:,2,:]
            lines_thar_co  = thar_S_ob[:,1,:]
            iv_thar_co     = thar_S_ob[:,2,:]
        
            All_Pixel_Centers = np.array([])
            All_Wavelengths   = np.array([])
            All_Orders        = np.array([])
            All_Centroids     = np.array([])
            All_Sigmas        = np.array([])
            All_Intensities   = np.array([])
            All_residuals     = np.array([])

            All_Pixel_Centers_co = np.array([])
            All_Wavelengths_co   = np.array([])
            All_Orders_co        = np.array([])
            All_Centroids_co     = np.array([])
            All_Sigmas_co        = np.array([])
            All_Intensities_co   = np.array([])
            All_residuals_co     = np.array([])

            order = 0
            while order < nord_ob:
                order_s = str(order)
                if (order < 10):
                    order_s = '0'+str(order)
                
                thar_order_orig = lines_thar_ob[order,:]
                IV              = iv_thar_ob[order,:]
                wei             = np.sqrt( IV )
                thar_order      = thar_order_orig

                thar_order_orig_co = lines_thar_co[order,:]
                IV_co              = iv_thar_co[order,:]
                wei_co             = np.sqrt( IV_co )
                thar_order_co      = thar_order_orig_co

                coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths,\
                    rms_ms, residuals, centroids, sigmas, intensities =\
                        GLOBALutils.Initial_Wav_Calibration(order_dir+'order_'+\
                        order_s+'.iwdat', thar_order, order, wei, rmsmax=100, \
                        minlines=10,FixEnds=False,Dump_Argon=dumpargon,\
                        Dump_AllLines=True, Cheby=use_cheby,porder=5,del_width=4.0)
                coeffs_pix2wav_co, coeffs_pix2sigma_co, pixel_centers_co, wavelengths_co,\
                    rms_ms_co, residuals_co, centroids_co, sigmas_co, intensities_co =\
                        GLOBALutils.Initial_Wav_Calibration(order_dir+'order_'+\
                        order_s+'.iwdat', thar_order_co, order, wei_co, rmsmax=100, \
                        minlines=10,FixEnds=False,Dump_Argon=dumpargon,\
                        Dump_AllLines=True, Cheby=use_cheby,porder=5,del_width=4.0)
                #plot(GLOBALutils.Cheby_eval( coeffs_pix2wav, np.arange(len(thar_order)), len(thar_order) ), thar_order)
                #print order
                #plot(pixel_centers, wavelengths,'ro')
                #show()
                if (order == 20): 
                    if (use_cheby):
                        Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, int(.5*len(thar_order)), len(thar_order) )
                    else:
                        Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )
                if order != 0:
                    All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
                    All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
                    All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order )
                    All_Centroids     = np.append( All_Centroids, centroids)
                    All_Sigmas        = np.append( All_Sigmas, sigmas)
                    All_Intensities   = np.append( All_Intensities, intensities )
                    All_residuals     = np.append( All_residuals, residuals )

                    All_Pixel_Centers_co = np.append( All_Pixel_Centers_co, pixel_centers_co )
                    All_Wavelengths_co   = np.append( All_Wavelengths_co, wavelengths_co )
                    All_Orders_co        = np.append( All_Orders_co, np.zeros( len(pixel_centers_co) ) + order )
                    All_Centroids_co     = np.append( All_Centroids_co, centroids_co)
                    All_Sigmas_co        = np.append( All_Sigmas_co, sigmas_co)
                    All_Intensities_co   = np.append( All_Intensities_co, intensities_co )
                    All_residuals_co     = np.append( All_residuals_co, residuals_co )
                order+=1
            #show()
            p0    = np.zeros( npar_wsol )
            p0[0] =  (20+OO0) * Global_ZP
            p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
                GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders,\
                                                    np.ones(All_Intensities.shape), p0, Cheby=use_cheby,\
                                                    maxrms=MRMS, Inv=Inverse_m,minlines=1200,order0=OO0, \
                                                    ntotal=nord_ob,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
            p1_co, G_pix_co, G_ord_co, G_wav_co, II_co, rms_ms_co, G_res_co = \
                GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co, All_Wavelengths_co, All_Orders_co,\
                                                    np.ones(All_Intensities_co.shape), p1, Cheby=use_cheby,\
                                                    maxrms=MRMS, Inv=Inverse_m,minlines=1200,order0=OO0, \
                                                    ntotal=nord_co,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
            #for i in np.unique(G_ord):
            #    I = np.where(G_ord == i)[0]
            #    plot(G_wav[I],G_res[I],'.')
            #show()

            pdict = {'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                     'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Wavelengths':All_Wavelengths, 'All_Orders':All_Orders,\
                     'All_Pixel_Centers':All_Pixel_Centers, 'All_Sigmas':All_Sigmas,\
                     'p1_co':p1_co, 'G_pix_co':G_pix_co, 'G_ord_co':G_ord_co, 'G_wav_co':G_wav_co, 'II_co':II_co, 'rms_ms_co':rms_ms_co,\
                     'G_res_co':G_res_co, 'All_Centroids_co':All_Centroids_co, 'All_Wavelengths_co':All_Wavelengths_co, 'All_Orders_co':All_Orders_co,\
                     'All_Pixel_Centers_co':All_Pixel_Centers_co, 'All_Sigmas_co':All_Sigmas_co,\
                     }
            pickle.dump( pdict, open( wavsol_pkl, 'w' ) )

        else:
            print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl

wavsol_pkl = dirout + ThAr_ref[0].split('/')[-1][:-4]+'wavsolpars.pkl'
dct        = pickle.load(open(wavsol_pkl,'r'))
p_ref      = dct['p1']
if mode != 'so':
    p_ref_co      = dct['p1_co']
"""
mjds_thar,shifts = [],[]
print '\n\tDetermination of instrumental drift...'
for i in range(len(ThAr_ref_dates)):
    wavsol_pkl = dirout + ThAr_ref[i].split('/')[-1][:-4]+'wavsolpars.pkl'
    hthar = pyfits.open( ThAr_ref[i] )
    npix = espadonsutils.OverscanTrim(hthar[0].data).shape[0]
    mjd, mjd0 = espadonsutils.mjd_fromheader( hthar )
    mjds_thar.append(mjd)
    dct = pickle.load(open(wavsol_pkl,'r'))
    nord = nord_all
    if mode != 'so':
        nord = nord_ob
    p_shift, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
                GLOBALutils.Global_Wav_Solution_vel_shift(dct['All_Pixel_Centers'],\
                dct['All_Wavelengths'], dct['All_Orders'],np.ones(len(dct['All_Orders'])), p_ref,\
                minlines=1200, maxrms=MRMS,order0=OO0, ntotal=nord,\
                Cheby=use_cheby, Inv=Inverse_m, npix=npix,nx=ncoef_x,nm=ncoef_m)

    shifts.append(p_shift[0])
"""

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
    if object2do in h[0].header['OBJNAME'] or object2do == 'all':
        new_list.append(fsim)
        mjd, mjd0 = espadonsutils.mjd_fromheader( h )
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
    mjd,mjd0      = espadonsutils.mjd_fromheader(h)
    ronoise, gain = float(h[0].header['RDNOISEA']),float(h[0].header['GAINA'])

    # Object name
    obname    = h[0].header['OBJNAME']
    print "\t\tObject name:",obname

    data = espadonsutils.OverscanTrim( h[0].data ) - MasterBias
    data = data.T

    bacfile = dirout + 'BAC_' + fsim.split('/')[-1][:-4]+'fits'
    if os.access(bacfile,os.F_OK)==False:

        if mode == 'so':
            Centers = np.zeros((len(c_all),data.shape[1]))
            for i in range(c_all.shape[0]):
                Centers[i,:]=scipy.polyval(c_all[i],np.arange(len(Centers[i,:])))
        else:
            Centers = np.zeros((len(c_ob),data.shape[1]))
            for i in range(c_ob.shape[0]):
                Centers[i,:]=0.5*( scipy.polyval(c_ob[i],np.arange(len(Centers[i,:]))) + \
                              scipy.polyval(c_co[i],np.arange(len(Centers[i,:]))) )
        bac = GLOBALutils.get_scat(data,Centers,span=15)
        if (os.access(bacfile,os.F_OK)):
            os.remove( bacfile )
        hdbac = pyfits.PrimaryHDU( bac )
        hdbac.writeto(bacfile)
    else:
        bac = pyfits.getdata(bacfile)

    data = data - bac
    ra          = float(h[0].header['RA_DEG'])
    dec         = float(h[0].header['DEC_DEG'])
    altitude    = 4204.
    latitude    = float(h[0].header['LATITUDE'])
    longitude   = float(h[0].header['LONGITUD'])
    epoch       = float(h[0].header['EQUINOX'])

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
    gobs.name = h[0].header['TELESCOP']
    gobs.lat  = rad(latitude)  # lat/long in decimal degrees  
    gobs.long = rad(longitude)
    timeT = h[0].header['UTC-OBS'].split(':')
    if len(timeT[0]) == 1:
        gobs.date = h[0].header['DATE-OBS'][:10] + ' 0' + h[0].header['UTC-OBS']
    else:
        gobs.date = h[0].header['DATE-OBS'][:10] + ' ' + h[0].header['UTC-OBS']
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

    if mode == 'so':
        sci_fits        = dirout + fsim.split('/')[-1][:-4]+'spec.fits.S'
        sci_fits_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.fits.S'

        if ( os.access(sci_fits,os.F_OK) == False ) or ( os.access(sci_fits_simple,os.F_OK) == False ) or \
           ( force_sci_extract ):

            sci_Ss = GLOBALutils.simple_extraction(data,c_all,ext_aperture,\
                                                      min_extract_col,max_extract_col,npools)
            sci_S  = GLOBALutils.optimal_extraction(data,P,c_all,ext_aperture,\
                                                       ronoise,gain,S_Marsh,NCosmic_Marsh,\
                                                       min_extract_col,max_extract_col,npools)
                
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
    else:
        sci_ob_fits        = dirout + fsim.split('/')[-1][:-4]+'spec.ob.fits.S'
        sci_ob_fits_simple = dirout + fsim.split('/')[-1][:-4]+'spec.ob.simple.fits.S'
        sci_co_fits        = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
        sci_co_fits_simple = dirout + fsim.split('/')[-1][:-4]+'spec.co.simple.fits.S'

        if ( os.access(sci_ob_fits,os.F_OK) == False ) or ( os.access(sci_ob_fits_simple,os.F_OK) == False ) or \
           ( os.access(sci_co_fits,os.F_OK) == False ) or ( os.access(sci_co_fits_simple,os.F_OK) == False ) or \
           ( force_sci_extract ):

            sci_ob_Ss = GLOBALutils.simple_extraction(data,c_ob,ext_aperture,\
                                                      min_extract_col,max_extract_col,npools)
            sci_ob_S  = GLOBALutils.optimal_extraction(data,P_ob,c_ob,ext_aperture,\
                                                       ronoise,gain,S_Marsh,NCosmic_Marsh,\
                                                       min_extract_col,max_extract_col,npools)

            sci_co_Ss = GLOBALutils.simple_extraction(data,c_co,ext_aperture,\
                                                      min_extract_col,max_extract_col,npools)
            sci_co_S  = GLOBALutils.optimal_extraction(data,P_co,c_co,ext_aperture,\
                                                       ronoise,gain,S_Marsh,NCosmic_Marsh,\
                                                       min_extract_col,max_extract_col,npools)
                
            if (os.access(sci_ob_fits,os.F_OK)):
                os.remove( sci_ob_fits )
            if (os.access(sci_ob_fits_simple,os.F_OK)):
                os.remove( sci_ob_fits_simple )
            hdu = pyfits.PrimaryHDU( sci_ob_S )
            hdu.writeto( sci_ob_fits )
            hdu = pyfits.PrimaryHDU( sci_ob_Ss )
            hdu.writeto( sci_ob_fits_simple )

            if (os.access(sci_co_fits,os.F_OK)):
                os.remove( sci_co_fits )
            if (os.access(sci_co_fits_simple,os.F_OK)):
                os.remove( sci_co_fits_simple )
            hdu = pyfits.PrimaryHDU( sci_co_S )
            hdu.writeto( sci_co_fits )
            hdu = pyfits.PrimaryHDU( sci_co_Ss )
            hdu.writeto( sci_co_fits_simple )
        

        else:
            print '\t\t\t '+fsim+" has already been extracted, reading in product fits files..."
            if mode == 'so':
                sci_S  = pyfits.getdata( sci_fits )
                sci_Ss = pyfits.getdata( sci_fits_simple )
            else:
                sci_ob_S  = pyfits.getdata( sci_ob_fits )
                sci_ob_Ss = pyfits.getdata( sci_ob_fits_simple )
                sci_co_S  = pyfits.getdata( sci_co_fits )
                sci_co_Ss = pyfits.getdata( sci_co_fits_simple )

    fout    = 'proc/' + obname + '_' + h[0].header['DATE-OBS'][:4] + h[0].header['DATE-OBS'][5:7] + h[0].header['DATE-OBS'][8:10] + '_' +'UT' + h[0].header['UTC-OBS'] + '_sp.fits'
    fout_co = 'proc/' + obname + '_' + h[0].header['DATE-OBS'][:4] + h[0].header['DATE-OBS'][5:7] + h[0].header['DATE-OBS'][8:10] + '_' +'UT' + h[0].header['UTC-OBS'] + '_sp_co.fits'

    if ( os.access(dirout+fout ,os.F_OK) == False ) or (force_spectral_file_build):
        if mode == 'so':
            nord = nord_all
        else:
            nord = nord_ob
        spec = np.zeros((11, nord, data.shape[1]))
        hdu  = pyfits.PrimaryHDU( spec )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MJD', mjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MBJD', mbjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START DATE', h[0].header['DATE-OBS'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START UT',  h[0].header['UTC-OBS'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TEXP (S)',h[0].header['EXPTIME'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (KM/S)', bcvel_baryc)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH (LAMBDA_BARY / LAMBDA_TOPO)', lbary_ltopo)    
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARGET NAME', obname)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA',h[0].header['RA'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC',h[0].header['DEC'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA BARY',ra)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC BARY',dec)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH EQUINOX',h[0].header['EQUINOX'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LATITUDE',latitude)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LONGITUDE',longitude)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS ALTITUDE',altitude)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARG AIRMASS START',h[0].header['AIRMASS'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MOON VEL',refvel)

        if mode == 'po':
            spec_co = np.zeros((11, nord, data.shape[1]))
            hdu_co  = pyfits.PrimaryHDU( spec_co )
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH MJD', mjd)
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH MBJD', mbjd)
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH SHUTTER START DATE', h[0].header['DATE-OBS'])
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH SHUTTER START UT',  h[0].header['UTC-OBS'])
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH TEXP (S)',h[0].header['EXPTIME'])
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH BARYCENTRIC CORRECTION (KM/S)', bcvel_baryc)
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH (LAMBDA_BARY / LAMBDA_TOPO)', lbary_ltopo)    
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH TARGET NAME', obname)
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH RA',h[0].header['RA'])
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH DEC',h[0].header['DEC'])
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH RA BARY',ra)
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH DEC BARY',dec)
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH EQUINOX',h[0].header['EQUINOX'])
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH OBS LATITUDE',latitude)
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH OBS LONGITUDE',longitude)
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH OBS ALTITUDE',altitude)
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH TARG AIRMASS START',h[0].header['AIRMASS'])
            hdu_co = GLOBALutils.update_header(hdu_co,'HIERARCH MOON VEL',refvel)

        if mode == 'so':
            nord = nord_all
        else:
            nord = nord_ob
            sci_S = sci_ob_S.copy()
            S_flat = S_flat_ob.copy()
            S_flat_n = S_flat_ob_n.copy()
            Snorms = Snorms_ob.copy()
        S_flat = S_flat / Snorms.max()
        Snorms = Snorms / Snorms.max()
        print '\t\tWavelength calibration:'
        #print "\t\t\tInstrumental drift:",(1e-6*p_shift)*299792458.0
        # Apply new wavelength solution including barycentric correction
        equis = np.arange( data.shape[1] )
        #p_shift = scipy.interpolate.splev(mjd,tck_shift)
        ind   = 0
        while ind <  nord:
            order = ind + OO0
            m     = order
            chebs  = GLOBALutils.Calculate_chebs(equis, m, npix=data.shape[1], order0=OO0, ntotal=nord, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
            WavSol = GLOBALutils.ToVacuum( lbary_ltopo  * (1.0/float(m)) * \
                GLOBALutils.Joint_Polynomial_Cheby(p_ref,chebs,ncoef_x,ncoef_m) )

            spec[0,ind,:] = WavSol
            spec[1,ind,:] = sci_S[ind,1]
            spec[2,ind,:] = sci_S[ind,2]
            fn  = S_flat[ind,1,:]
            L  = np.where( fn == 0 )[0]
            spec[3,ind,:] = spec[1,ind,:] / S_flat[ind,1,:]
            spec[4,ind,:] = spec[2,ind] * ( S_flat_n[ind,1,:] ** 2 )
            spec[3,ind,L] = 0.
            spec[4,ind,L] = 0.
            #plot(spec[0,ind],spec[3,ind])
            ccoef = GLOBALutils.get_cont_single(spec[0,ind],spec[3,ind],spec[4,ind],ll=1.5,lu=5,nc=nconts[ind])

            L  = np.where( spec[1,ind] != 0 )
            spec[5,ind,:][L] = spec[3,ind][L] / np.polyval(ccoef,spec[0,ind][L])    
            ratio            = np.polyval(ccoef,spec[0,ind][L]) * Snorms[ind]
            spec[6,ind,:][L] = spec[4,ind][L] * (ratio ** 2 )
            spec[7,ind,:][L] = ratio
            spec[8,ind,:][L] = ratio * S_flat_n[ind,1][L] / np.sqrt( ratio * S_flat_n[ind,1][L] / gain + (ronoise/gain)**2 )
            spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
            dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
            NN            = np.average(dlambda_dx)
            dlambda_dx    /= NN
            LL = np.where(spec[5,ind] > 1 + 10. / scipy.signal.medfilt(spec[8,ind],21))[0]
            spec[5,ind,LL] = 1.

            spec[9,ind][L] = spec[5,ind][L] * (dlambda_dx[L] ** 1) 
            spec[10,ind][L] = spec[6,ind][L] / (dlambda_dx[L] ** 2)
            ind +=1

        if os.access(dirout + fout, os.F_OK):
            os.remove(dirout + fout)
        hdu.writeto(dirout + fout)

        if mode == 'po':
            nord = nord_co
            sci_S = sci_co_S.copy()
            S_flat = S_flat_co.copy()
            S_flat_n = S_flat_co_n.copy()
            Snorms = Snorms_co.copy()
            S_flat = S_flat / Snorms.max()
            Snorms = Snorms / Snorms.max()
            ind   = 0
            while ind <  nord:
                order = ind + OO0
                m     = order
                chebs  = GLOBALutils.Calculate_chebs(equis, m, npix=data.shape[1], order0=OO0, ntotal=nord, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
                WavSol = GLOBALutils.ToVacuum( lbary_ltopo  * (1.0/float(m)) * \
                    GLOBALutils.Joint_Polynomial_Cheby(p_ref_co,chebs,ncoef_x,ncoef_m) )

                spec_co[0,ind,:] = WavSol
                spec_co[1,ind,:] = sci_S[ind,1]
                spec_co[2,ind,:] = sci_S[ind,2]
                fn  = S_flat[ind,1,:]
                L  = np.where( fn == 0 )[0]
                spec_co[3,ind,:] = spec_co[1,ind,:] / S_flat[ind,1,:]
                spec_co[4,ind,:] = spec_co[2,ind] * ( S_flat[ind,1,:] ** 2 )
                spec_co[3,ind,L] = 0.
                spec_co[4,ind,L] = 0.
                #plot(spec[0,ind],spec[3,ind])
                ccoef = GLOBALutils.get_cont_single(spec[0,ind],spec[3,ind],spec[4,ind],ll=1.5,lu=5,nc=nconts[ind])

                L  = np.where( spec_co[1,ind] != 0 )
                spec_co[5,ind,:][L] = spec_co[3,ind][L] / np.polyval(ccoef,spec_co[0,ind][L])    
                ratio            = np.polyval(ccoef,spec_co[0,ind][L]) * Snorms[ind]
                spec_co[6,ind,:][L] = spec_co[4,ind][L] * (ratio ** 2 )
                spec_co[7,ind,:][L] = ratio
                spec_co[8,ind,:][L] = ratio * S_flat_n[ind,1][L] / np.sqrt( ratio * S_flat_n[ind,1][L] / gain + (ronoise/gain)**2 )
                spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
                dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
                NN            = np.average(dlambda_dx)
                dlambda_dx    /= NN
                LL = np.where(spec_co[5,ind] > 1 + 10. / scipy.signal.medfilt(spec_co[8,ind],21))[0]
                spec_co[5,ind,LL] = 1.

                spec_co[9,ind][L] = spec_co[5,ind][L] * (dlambda_dx[L] ** 1) 
                spec_co[10,ind][L] = spec_co[6,ind][L] / (dlambda_dx[L] ** 2)
                ind +=1

            if os.access(dirout + fout_co, os.F_OK):
                os.remove(dirout + fout_co)
            hdu_co.writeto(dirout + fout_co)

    if (not JustExtract):

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

            hdu = GLOBALutils.update_header(hdu,'HIERARCH SIMBAD SPTYP', sp_type_query)

            pars_file = dirout + fsim.split('/')[-1][:-4]+'_stellar_pars.txt'

            if os.access(pars_file,os.F_OK) == False or force_stellar_pars:
                print "\t\t\tEstimating atmospheric parameters:"
                Rx = np.around(1./np.sqrt(1./40000.**2 - 1./float(RES)**2))
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
                p1_m,XCmodel_m,p1gau_m,XCmodelgau_m,Ls2_m = GLOBALutils.XC_Final_Fit( vels, \
                    xc_av , sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = True)
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

        pkl_xc = dirout + fsim.split('/')[-1][:-4]+obname+'_XC_'+sp_type+'.pkl'
        pickle.dump( xc_dict, open( pkl_xc, 'w' ) )

        ccf_pdf = dirout + 'proc/' + fsim.split('/')[-1][:-4] + obname + '_XCs_' + sp_type + '.pdf'

        if not avoid_plot:
            GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)

        SNR_5130 = np.median(spec[8,22,1800:2201] )
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

        if (RVerr2 <= 0.010):
            RVerr2 = 0.010

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
        hdu = GLOBALutils.update_header(hdu,'INST', 'ESPADONS')
        hdu = GLOBALutils.update_header(hdu,'RESOL', '80000')
        hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')
        hdu = GLOBALutils.update_header(hdu,'XC_MIN', XC_min)
        hdu = GLOBALutils.update_header(hdu,'BJD_OUT', bjd_out)

        # write to output
        line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f espadons   ceres %9d %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
            (obname, bjd_out, RV, RVerr2, BS, BSerr, RES, T_eff_epoch, logg_epoch, Z_epoch, vsini_epoch, XC_min, disp_epoch,\
            TEXP, SNR_5130_R, ccf_pdf)
        f_res.write(line_out)
    
        if (os.access( dirout + fout,os.F_OK)):
            os.remove( dirout + fout)
        hdu.writeto( dirout + fout )
    else:
        print "Reading spectral file from", fout
        spec = pyfits.getdata( fout )

f_res.close()