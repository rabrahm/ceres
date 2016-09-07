import sys
from pylab import *
base = '../'

sys.path.append(base+"utils/Correlation")
sys.path.append(base+"utils/OptExtract")
sys.path.append(base+"utils/GLOBALutils")

baryc_dir= base+'utils/SSEphemx/'
sys.path.append(baryc_dir)
ephemeris='DEc403'

# ceres modules
import fideosutils
import correlation
import Marsh
import GLOBALutils

# other useful modules
import argparse
import ephem
import glob
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

args = parser.parse_args()
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
## perhaps put into options ##
force_sci_extract  = False
force_thar_extract = False
force_pre_process  = False
force_flat_extract = False
force_tharxc       = False
force_thar_wavcal  = False
force_spectral_file_build = True
dumpargon          = False
minlines_glob_ob   = 700
minlines_glob_co   = 500

Inverse_m          = True
use_cheby          = True
MRMS               = 200   # max rms in m/s, global wav solution

order0		   = 49
ncoef_x		   = 4
ncoef_m		   = 6
npar_wsol          = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2
trace_degree       = 4
Marsh_alg          = 0
ext_aperture_ob    = 7
ext_aperture_co    = 5
NSigma_Marsh       = 10
NCosmic_Marsh      = 10
S_Marsh            = 0.4
N_Marsh            = 3      # grado polinomio 
min_extract_col    = 0
max_extract_col    = 2048
n_useful           = 70    # up to which order do we care?

order_dir   = "wavcals/"

RO_ob, GA_ob = 6.82,1.5
ronoise, gain = 6.82,1.5

print "\n\n\tFIDEOS ESO1.0m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'

if os.access(dirin + 'log.txt', os.F_OK):
	f = open(dirin + 'log.txt','r')
	lines = f.readlines()
	images, types = [],[]
	for line in lines:
		cos = line.split()
		hd = pyfits.getheader(dirin + cos[0])
		print line,hd['EXPTIME']
		images.append(dirin + cos[0])
		types.append(cos[1])
	images, types = np.array(images), np.array(types)
	I = np.where(types=='bias')
	biases = images[I]
	I = np.where(types=='dark')
	darks = images[I]
	I = np.where(types=='thar_thar')
	thars = images[I]
	I = np.where(types=='flat_flat')
	flats = images[I]
	I = np.where(types=='star_thar')
	sci_sim = images[I]
	I = np.where(types=='dark_flat')
	flats_co = images[I]
	I = np.where(types=='star_dark')
	sci_one = images[I]
else:
	biases  = glob.glob(dirin + 'bias*fit')
	flats   = glob.glob(dirin + 'flat*fit')
	thars   = glob.glob(dirin + 'ThAr*fit')
	sci_sim	= glob.glob(dirin + 'test_interfaz_ThAr*fit')
	#flat_ob = glob.glob(dirin + 'FIDEOS_flat_dark*fits')[0]
	#flat_co = glob.glob(dirin + 'FIDEOS_dark_flat*fits')[0]
	#thars   = glob.glob(dirin + 'FIDEOS_thar_thar*fits')[:2]
	#thars   = glob.glob(dirin + 'ThAr_sc*fit')[:2]


if ( (os.access(dirout+'MasterFlat.fits',os.F_OK) == False) or \
     (os.access(dirout+'trace.pkl',os.F_OK) == False)  or \
     (os.access(dirout+'MasterBias.fits',os.F_OK) == False)  or \
     (force_pre_process) ):
    print "\tNo previous pre-processing files or found"
    pre_process = 1
else:
    print "\tPre-processing files found, going straight to extraction"
    pre_process = 0

trace_file = dirout+'trace.pkl'
if (pre_process == 1):
    #c_ob, nord_ob = GLOBALutils.get_them(Flat_ob, 10.,trace_degree,maxords=-1)
    #c_ob, nord_ob = FIDEOSutils.good_orders(c_ob,nord_ob,Flat_ob.shape[0],Flat_ob.shape[1],ext_aperture)	
    #print nord_ob, 'object orders found'
    #c_co, nord_co = GLOBALutils.get_them(Flat_co, 10.,trace_degree,maxords=-1)
    #c_co, nord_co = FIDEOSutils.good_orders(c_co,nord_co,Flat_co.shape[0],Flat_co.shape[1],ext_aperture)
    #print nord_co, 'comparison orders found'
    MasterBias = GLOBALutils.MedianCombine_simple(biases,ZF=0.)
    hdu = pyfits.PrimaryHDU( MasterBias )
    if (os.access(dirout+'MasterBias.fits',os.F_OK)):
        os.remove(dirout+'MasterBias.fits')
    hdu.writeto(dirout+'MasterBias.fits')

    MasterFlat = GLOBALutils.MedianCombine_simple(flats,ZF=MasterBias)
    MasterFlat = np.fliplr(MasterFlat.T)
    hdu = pyfits.PrimaryHDU( MasterFlat )
    if (os.access(dirout+'MasterFlat.fits',os.F_OK)):
        os.remove(dirout+'MasterFlat.fits')
    hdu.writeto(dirout+'MasterFlat.fits')

    MasterFlat_co = GLOBALutils.MedianCombine_simple(flats_co,ZF=MasterBias)
    MasterFlat_co = np.fliplr(MasterFlat_co.T)
    hdu = pyfits.PrimaryHDU( MasterFlat_co )
    if (os.access(dirout+'MasterFlat_co.fits',os.F_OK)):
        os.remove(dirout+'MasterFlat_co.fits')
    hdu.writeto(dirout+'MasterFlat_co.fits')

    #c_co, nord_co = GLOBALutils.get_them(MasterFlat_co.T, 5.,trace_degree,maxords=-1,nsigmas=5)
    #Centers = np.zeros((nord_co,MasterFlat_co.shape[0]))
    #ejx = np.arange(MasterFlat_co.shape[0])
    #for i in range(nord_co):
    #    Centers[i,:]=scipy.polyval(c_co[i],ejx)
    #bac_co = GLOBALutils.get_scat(MasterFlat_co.T,Centers,span = 6)
    #MasterFlat_co = MasterFlat_co - bac_co.T
    #MasterFlat_ob = MasterFlat - MasterFlat_co

    #plot(MasterFlat_ob[1000])
    #show()
    #plot(np.median(MasterFlat[:,2050:],axis=1))
    #show()
    #plot(MasterFlat[1000])
    #plot(MasterFlat_co[1000])
    #show()

    c_all, nord_all = GLOBALutils.get_them(MasterFlat.T, 4.,trace_degree,maxords=-1,nsigmas=1,mode=1)
    print nord_all
    c_ob, c_co, nord_ob, nord_co = fideosutils.clean_orders(c_all,MasterFlat.T)
    print nord_ob, nord_co
    c_co, nord_co = fideosutils.good_orders(c_co,nord_co,MasterFlat.shape[1],MasterFlat.shape[0],ext_aperture_co)
    c_ob, nord_ob = fideosutils.good_orders(c_ob,nord_ob,MasterFlat.shape[1],MasterFlat.shape[0],ext_aperture_ob)
    print nord_ob, nord_co
    print nord_ob, 'object orders found'
    print nord_co, 'comparison orders found'
    trace_dict = {'c_ob':c_ob,'nord_ob':nord_ob,'c_co':c_co,'nord_co':nord_co}
    pickle.dump( trace_dict, open( trace_file, 'w' ) )

else:
    trace_dict = pickle.load( open( trace_file, 'r' ) )
    c_ob    = trace_dict['c_ob']
    nord_ob = len(c_ob)
    c_co    = trace_dict['c_co']
    nord_co = len(c_co)
    MasterBias = pyfits.getdata(dirout+'MasterBias.fits')
    MasterFlat = pyfits.getdata(dirout+'MasterFlat.fits')

if len(darks)>1:
    dtimes = GLOBALutils.get_dark_times(darks,key='EXPTIME')
    for dtime in dtimes:
        tdarks = GLOBALutils.get_tdarks(darks,dtime,key='EXPTIME')
        dark   = GLOBALutils.MedianCombine_simple(tdarks,ZF=MasterBias)
        hdu = pyfits.PrimaryHDU( dark )
        if (os.access(dirout+'Dark_'+str(dtime)+'.fits',os.F_OK)):
            os.remove(dirout+'Dark_'+str(dtime)+'.fits')
        hdu.writeto(dirout+'Dark_'+str(dtime)+'.fits')

c_all = GLOBALutils.Mesh(c_ob,c_co)
nord_all = nord_ob + nord_co

P_ref_ob = np.zeros( MasterFlat.T.shape )
P_ref_co = np.zeros( MasterFlat.T.shape )
P_ref_fits = dirout + 'P_ref.fits'

if ( os.access(P_ref_fits,os.F_OK) == False ) or (force_flat_extract):
	Flat = MasterFlat.T
	P_ref_ob = GLOBALutils.obtain_P(Flat,c_ob,ext_aperture_ob,RO_ob,\
                                    GA_ob,NSigma_Marsh, S_Marsh, \
				    N_Marsh, Marsh_alg, min_extract_col,\
				    max_extract_col, npools)
	P_ref_co = GLOBALutils.obtain_P(Flat,c_co,ext_aperture_co,RO_ob,\
                                    GA_ob,NSigma_Marsh, S_Marsh, \
				    N_Marsh, Marsh_alg, min_extract_col,\
				    max_extract_col, npools)

	P_ref = P_ref_co + P_ref_ob
	if (os.access(P_ref_fits,os.F_OK)):
	    os.remove( P_ref_fits )
	hdu = pyfits.PrimaryHDU( P_ref )
	hdu.writeto( P_ref_fits )
else:
	P_ref = pyfits.getdata(P_ref_fits)

bac = fideosutils.get_scat(MasterFlat, P_ref,1)
MasterFlat -= bac

P_ob_fits = dirout + 'P_ob.fits'
P_ob = np.zeros( MasterFlat.T.shape )
S_flat_ob_fits = dirout +'S_flat_ob.fits'
S_flat_ob_simple_fits = dirout +'S_flat_ob_simple.fits'

if ( os.access(P_ob_fits,os.F_OK) == False ) or ( os.access(S_flat_ob_fits,os.F_OK) == False ) or \
   ( os.access(S_flat_ob_simple_fits,os.F_OK) == False ) or (force_flat_extract):
    print 'Extracting Flat of object fibre...'
    Flat = MasterFlat.T #- bac
    Centers = np.zeros((nord_ob,MasterFlat.shape[0]))
    ejx = np.arange(MasterFlat.shape[0])
    imshow(Flat,vmin=0,vmax = 10000)
    for i in range(nord_ob):
        Centers[i,:]=scipy.polyval(c_ob[i],ejx)
	plot(ejx,Centers[i,:],'r')
    show()

    #bac = FIDEOSutils.scat_flat(MasterFlat.T,c_ob,c_co)


    for i in range(nord_ob):
	print i
	try:
		P_marsh = GLOBALutils.PCoeff( Flat, c_ob[i,:], ext_aperture_ob, RO_ob, GA_ob,\
						      NSigma_Marsh, S_Marsh, N_Marsh, Marsh_alg ,\
						      min_extract_col,max_extract_col )
		P_ob    += P_marsh
	except:
		print 'problem'
	i+=1

    S_flat_ob_simple = GLOBALutils.simple_extraction(Flat,c_ob,ext_aperture_ob,min_extract_col,\
                                                         max_extract_col,npools)
    S_flat_ob  = GLOBALutils.optimal_extraction(Flat,P_ob,c_ob,ext_aperture_ob,RO_ob,GA_ob,S_Marsh,\
                                                    NCosmic_Marsh,min_extract_col,max_extract_col,npools)

    S_flat_ob_simple = GLOBALutils.invert(S_flat_ob_simple)
    S_flat_ob = GLOBALutils.invert(S_flat_ob)

    if (os.access(P_ob_fits,os.F_OK)):
        os.remove( P_ob_fits )
    if (os.access(S_flat_ob_fits,os.F_OK)):
        os.remove( S_flat_ob_fits )
    if (os.access(S_flat_ob_simple_fits,os.F_OK)):
        os.remove( S_flat_ob_simple_fits )

    hdu = pyfits.PrimaryHDU( P_ob )
    hdu.writeto( P_ob_fits )
    hdu = pyfits.PrimaryHDU( S_flat_ob )
    hdu.writeto( S_flat_ob_fits )
    hdu = pyfits.PrimaryHDU( S_flat_ob_simple )
    hdu.writeto( S_flat_ob_simple_fits )
   
else:
    print "Extracted flat comparison spectra found, loading..."
    P_ob       = pyfits.getdata( P_ob_fits )
    S_flat_ob  = pyfits.getdata( S_flat_ob_fits )
    S_flat_ob_simple = pyfits.getdata( S_flat_ob_simple_fits )

P_co_fits = dirout + 'P_co.fits'
P_co = np.zeros( MasterFlat.T.shape )
S_flat_co_fits = dirout +'S_flat_co.fits'
S_flat_co_simple_fits = dirout +'S_flat_co_simple.fits'


if ( os.access(P_co_fits,os.F_OK) == False ) or ( os.access(S_flat_co_fits,os.F_OK) == False ) or \
   ( os.access(S_flat_co_simple_fits,os.F_OK) == False ) or (force_flat_extract):
    print 'Extracting Flat of comparison fibre...'
    Centers = np.zeros((nord_co,MasterFlat.shape[0]))
    Flat = MasterFlat.T
    #ejx = np.arange(MasterFlat.shape[0])
    #for i in range(nord_co):
    #    Centers[i,:]=scipy.polyval(c_co[i],ejx)
    #bac = GLOBALutils.get_scat(Flat_co,Centers,span=5)
    #Flat_co -= bac

    for i in range(nord_co):
	try:
		P_marsh = GLOBALutils.PCoeff( Flat, c_co[i,:], ext_aperture_co, RO_ob, GA_ob,\
						      NSigma_Marsh, S_Marsh, N_Marsh, Marsh_alg ,\
						      min_extract_col,max_extract_col )
		P_co    += P_marsh
	except:
		print 'problem'
	i+=1
    
    S_flat_co_simple = GLOBALutils.simple_extraction(Flat,c_co,ext_aperture_co,min_extract_col,\
                                                         max_extract_col,npools)
    S_flat_co  = GLOBALutils.optimal_extraction(Flat,P_co,c_co,ext_aperture_co,RO_ob,GA_ob,S_Marsh,\
                                                    NCosmic_Marsh,min_extract_col,max_extract_col,npools) 

    S_flat_co_simple = GLOBALutils.invert(S_flat_co_simple)
    S_flat_co = GLOBALutils.invert(S_flat_co)

    if (os.access(P_co_fits,os.F_OK)):
        os.remove( P_co_fits )
    if (os.access(S_flat_co_fits,os.F_OK)):
        os.remove( S_flat_co_fits )
    if (os.access(S_flat_co_simple_fits,os.F_OK)):
        os.remove( S_flat_co_simple_fits )

    hdu = pyfits.PrimaryHDU( P_co )
    hdu.writeto( P_co_fits )
    hdu = pyfits.PrimaryHDU( S_flat_co )
    hdu.writeto( S_flat_co_fits )
    hdu = pyfits.PrimaryHDU( S_flat_co_simple )
    hdu.writeto( S_flat_co_simple_fits )
   
else:
    print "Extracted flat comparison spectra found, loading..."
    P_co       = pyfits.getdata( P_co_fits )
    S_flat_co  = pyfits.getdata( S_flat_co_fits )
    S_flat_co_simple = pyfits.getdata( S_flat_co_simple_fits )

S_flat_ob_n, S_flat_ob_simple_n \
    = GLOBALutils.FlatNormalize( S_flat_ob, S_flat_ob_simple)
S_flat_co_n, S_flat_co_simple_n \
    = GLOBALutils.FlatNormalize( S_flat_co, S_flat_co_simple)

reffiles = []
for i in range(nord_ob):
	order_s = str(i)
	if i < 10:
	    order_s = '0'+str(i)
	reffiles.append(order_dir+'fideos_'+order_s+'.iwdat')
thars = thars[:1]
for fsim in thars:
    hthar = pyfits.open( fsim )
    dthar = pyfits.getdata( fsim )#- MasterBias
    dthar = np.fliplr(dthar.T)
    thar_fits_ob_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.ob.fits.S'
    thar_fits_co_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.co.fits.S'
    bac_fits_thar = dirout + 'BKG_'+ fsim.split('/')[-1][:-4]+'.fits'
    if ( os.access(thar_fits_ob_simple,os.F_OK) == False ) or ( os.access(thar_fits_co_simple,os.F_OK) == False ) or (force_thar_extract):
        print "No previous extraction or extraction forced for ThAr file", fsim, "extracting..."
	#Centers = np.zeros((nord_ob,dthar.shape[1]))
	#ejx = np.arange(dthar.shape[1])
        #for i in range(nord_ob):
        #    Centers[i,:]=scipy.polyval(c_ob[i],ejx)
        #bac = GLOBALutils.get_scat(dthar,Centers,span=10)
        #if (os.access(bac_fits_thar,os.F_OK)):
        #    os.remove( bac_fits_thar )
        #hdu = pyfits.PrimaryHDU( bac )
        #hdu.writeto( bac_fits_thar )
	bac = FIDEOSutils.get_scat(dthar, P_ob + P_co,1)
        dthar -= bac
	thar_Ss_ob_simple  = GLOBALutils.simple_extraction(dthar.T,c_ob,ext_aperture_ob,min_extract_col,max_extract_col,npools)
	thar_Ss_ob         = GLOBALutils.optimal_extraction(dthar.T,P_ob,c_ob,ext_aperture_ob,RO_ob,GA_ob,S_Marsh,\
                                                    NCosmic_Marsh,min_extract_col,max_extract_col,npools)
	thar_Ss_co_simple  = GLOBALutils.simple_extraction(dthar.T,c_co,ext_aperture_co,min_extract_col,max_extract_col,npools)
	thar_Ss_co         = GLOBALutils.optimal_extraction(dthar.T,P_co,c_co,ext_aperture_co,RO_ob,GA_ob,S_Marsh,\
                                                    NCosmic_Marsh,min_extract_col,max_extract_col,npools)
	
	thar_Ss_ob = GLOBALutils.invert(thar_Ss_ob)
	thar_Ss_co = GLOBALutils.invert(thar_Ss_co)

        if (os.access(thar_fits_ob_simple,os.F_OK)):
            os.remove( thar_fits_ob_simple )
        hdu = pyfits.PrimaryHDU( thar_Ss_ob )
        hdu.writeto( thar_fits_ob_simple )
        if (os.access(thar_fits_co_simple,os.F_OK)):
            os.remove( thar_fits_co_simple )
        hdu = pyfits.PrimaryHDU( thar_Ss_co )
        hdu.writeto( thar_fits_co_simple )


for fsim in thars:
    thar_fits_ob_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.ob.fits.S'
    thar_fits_co_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.co.fits.S'
    thar_spec_ob	= dirout + fsim.split('/')[-1][:-4]+'sp.ob.fits'
    thar_spec_co	= dirout + fsim.split('/')[-1][:-4]+'sp.co.fits'
    wavsol_pkl = dirout + fsim.split('/')[-1][:-4]+'wavsolpars.pkl'
    if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
	print "Working on initial ThAr file", fsim
        hthar = pyfits.open( fsim )
        mjd, mjd0 = FIDEOSutils.mjd_fromheader( hthar )
        thar_fits_ob = dirout + fsim.split('/')[-1][:-4]+'spec.ob.fits.S'
        thar_S_ob = pyfits.getdata( thar_fits_ob_simple )
        thar_fits_co = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
        thar_S_co = pyfits.getdata( thar_fits_co_simple )

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
        All_Pixel_Centers_co = np.array([])
        All_Wavelengths_co   = np.array([])
        All_Orders_co        = np.array([])
        All_Centroids_co     = np.array([])
        All_Sigmas_co        = np.array([])
        All_Intensities_co   = np.array([])
	
	offset_ob = GLOBALutils.get_rough_offset(lines_thar_ob,reffiles)
	offset_co = GLOBALutils.get_rough_offset(lines_thar_co,reffiles)

	t_ords,t_wavs = [],[]
        for order in range(nord_ob):
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)
            
            thar_order_orig = lines_thar_ob[order,:]
            IV              = iv_thar_ob[order,:]
            wei             = np.sqrt( IV )
            bkg             = np.zeros(len(thar_order_orig))# - GLOBALutils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
            thar_order      = thar_order_orig - bkg

            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
                = GLOBALutils.Initial_Wav_Calibration( order_dir+'fideos_'+order_s+'.iwdat', thar_order, order, wei, \
                                                   rmsmax=10000, minlines=10,FixEnds=False,Dump_Argon=False,Cheby=True,porder=ncoef_x,rough_shift = offset_ob)
            if (order == 25): 
                Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 1023, len(thar_order) )

            All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
            All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
            All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order )
            All_Centroids     = np.append( All_Centroids, centroids)
            All_Sigmas        = np.append( All_Sigmas, sigmas)
            All_Intensities   = np.append( All_Intensities, intensities )

	    t_ords.append(order)
	    t_wavs.append(GLOBALutils.Cheby_eval( coeffs_pix2wav, 1023, len(thar_order) ))
	t_ords,t_wavs = np.array(t_ords),np.array(t_wavs)
	#GLOBALutils.get_zero_order_number(t_ords,t_wavs)

        p0 = np.zeros( npar_wsol )
        p0[0] =  (25+49) * Global_ZP 
        p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
            GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders,\
                                                np.ones(All_Intensities.shape), p0, Cheby=True,\
                                                maxrms=MRMS, Inv=Inverse_m,minlines=400,order0=49, \
						ntotal=nord_ob,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

	thar_out = np.zeros((2,nord_ob,lines_thar_ob.shape[1]))
        equis = np.arange( lines_thar_ob.shape[1] )        
	order = 0
        while order < nord_ob:
            m   = order + 49
	    chebs = GLOBALutils.Calculate_chebs(equis, m, Inverse=Inverse_m,order0=49,ntotal=nord_ob,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
	    WavSol = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(p1,chebs,ncoef_x,ncoef_m)
            thar_out[0,order,:] = WavSol
            thar_out[1,order,:] = lines_thar_ob[order]
	    order+=1

	if os.access(thar_spec_ob,os.F_OK):
		os.system('rm '+ thar_spec_ob)
	hdu = pyfits.PrimaryHDU(thar_out)
	hdu.writeto(thar_spec_ob)

        for order in range(nord_co):
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)
            
            thar_order_orig = lines_thar_co[order,:]
            IV              = iv_thar_co[order,:]
            wei             = np.sqrt( IV )
            bkg             = np.zeros(len(thar_order_orig))# - GLOBALutils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
            thar_order      = thar_order_orig - bkg

            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
                = GLOBALutils.Initial_Wav_Calibration( order_dir+'fideos_'+order_s+'.iwdat', thar_order, order, wei, \
                                                   rmsmax=10000, minlines=10,FixEnds=False,Dump_Argon=False,Cheby=True,porder=ncoef_x,rough_shift = offset_co)
            if (order == 25): 
                Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 1023, len(thar_order) )

            All_Pixel_Centers_co = np.append( All_Pixel_Centers, pixel_centers )
            All_Wavelengths_co   = np.append( All_Wavelengths, wavelengths )
            All_Orders_co        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order )
            All_Centroids_co     = np.append( All_Centroids, centroids)
            All_Sigmas_co        = np.append( All_Sigmas, sigmas)
            All_Intensities_co   = np.append( All_Intensities, intensities )

        p0 = np.zeros( npar_wsol )
        p0[0] =  (25+49) * Global_ZP 
        p1_co, G_pix_co, G_ord_co, G_wav_co, II_co, rms_ms_co, G_res_co = \
            GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co, All_Wavelengths_co, All_Orders_co,\
                                                np.ones(All_Intensities_co.shape), p0, Cheby=True,\
                                                maxrms=MRMS, Inv=Inverse_m,minlines=400,order0=49, \
						ntotal=nord_co,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

	thar_out_co = np.zeros((2,nord_co,lines_thar_co.shape[1]))
        equis = np.arange( lines_thar_co.shape[1] )        
	order = 0
        while order < nord_co:
            m   = order + 49
	    chebs = GLOBALutils.Calculate_chebs(equis, m, Inverse=Inverse_m,order0=49,ntotal=nord_co,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
	    WavSol = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(p1_co,chebs,ncoef_x,ncoef_m)
            thar_out_co[0,order,:] = WavSol
            thar_out_co[1,order,:] = lines_thar_co[order]
	    order+=1

	if os.access(thar_spec_co,os.F_OK):
		os.system('rm '+ thar_spec_co)
	hdu = pyfits.PrimaryHDU(thar_out_co)
	hdu.writeto(thar_spec_co)
 
        pdict = {'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                 'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Orders':All_Orders, 'All_Sigmas':All_Sigmas,
                 'p1_co':p1_co, 'G_pix_co':G_pix_co, 'G_ord_co':G_ord_co, 'G_wav_co':G_wav_co, 'II_co':II_co, \
		 'rms_ms_co':rms_ms_co,'G_res_co':G_res_co, 'All_Centroids_co':All_Centroids_co \
		}
        pickle.dump( pdict, open( wavsol_pkl, 'w' ) )

    else:
        print "Using previously computed wavelength solution in file",wavsol_pkl

pkl_wsol = dirout + thars[-1].split('/')[-1][:-4]+'wavsolpars.pkl'
print " Unpickling wavelength solution from", pkl_wsol, " ..."
wsol_dict = pickle.load(open(pkl_wsol,'r'))
sci_all = np.hstack((sci_sim,sci_one))
sci_all = sci_one.copy()
lbary_ltopo = 1.
for fsim in sci_all[-2:]:
    hd = pyfits.getheader(fsim)
    exptime = hd['EXPTIME']
    if os.access(dirout+'Dark_'+str(exptime)+'.fits',os.F_OK):
	DARK = pyfits.getdata(dirout+'Dark_'+str(exptime)+'.fits')
    else:
	DARK = np.zeros(MasterBias.shape)
    hdat  = pyfits.open( fsim )
    data  = pyfits.getdata( fsim )- MasterBias - DARK
    data  = np.fliplr(data.T)
    #imshow(data.T,vmin=np.median(data),vmax=np.median(data)+500)
    #for i in range(nord_ob):
    #	ejx = np.arange(2048)
    #	y = np.polyval(c_ob[i],ejx)
    #	plot(ejx,y)
    #show()
    if fsim in sci_sim:
    	bac = FIDEOSutils.get_scat(data, P_ob + P_co,1)
    else:
	bac = FIDEOSutils.get_scat(data, P_ob,1)
    imshow(bac)	
    show()
    plot(data[1000])
    plot(bac[1000])
    show()
    data -= bac
    plot(data[1000])
    show()
    sci_fits_ob_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.ob.fits.S'
    sci_fits_co_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.co.fits.S'
    sci_fits_ob = dirout + fsim.split('/')[-1][:-4]+'spec.ob.fits.S'
    sci_fits_co = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
    bac_fits    = dirout + 'BKG_'+ fsim.split('/')[-1][:-4]+'.fits'

    if ( os.access(sci_fits_ob,os.F_OK) == False )        or ( os.access(sci_fits_co,os.F_OK) == False )        or \
       ( os.access(sci_fits_ob_simple,os.F_OK) == False ) or ( os.access(sci_fits_co_simple,os.F_OK) == False ) or \
       ( force_sci_extract ):
        print "No previous extraction or extraction forced for science file", fsim, "extracting..."

	sci_Ss_ob = GLOBALutils.simple_extraction(data.T,c_ob,ext_aperture_ob,min_extract_col,max_extract_col,npools)
	sci_S_ob  = GLOBALutils.optimal_extraction(data.T,P_ob,c_ob,ext_aperture_ob,ronoise,gain,S_Marsh,NCosmic_Marsh,min_extract_col,max_extract_col,npools)
	sci_Ss_ob = GLOBALutils.invert(sci_Ss_ob)
	sci_S_ob  = GLOBALutils.invert(sci_S_ob)
        if (os.access(sci_fits_ob_simple,os.F_OK)):
            os.remove( sci_fits_ob_simple )
	if (os.access(sci_fits_ob,os.F_OK)):
            os.remove( sci_fits_ob )
	hdu = pyfits.PrimaryHDU( sci_S_ob )
        hdu.writeto( sci_fits_ob )
        hdu = pyfits.PrimaryHDU( sci_Ss_ob )
        hdu.writeto( sci_fits_ob_simple )

	if fsim in sci_sim:
	    sci_Ss_co = GLOBALutils.simple_extraction(data.T,c_co,ext_aperture_co,min_extract_col,max_extract_col,npools)
	    sci_S_co  = GLOBALutils.optimal_extraction(data.T,P_co,c_co,ext_aperture_co,ronoise,gain,S_Marsh,NCosmic_Marsh,min_extract_col,max_extract_col,npools)
	    sci_Ss_co = GLOBALutils.invert(sci_Ss_co)
	    sci_S_co  = GLOBALutils.invert(sci_S_co)
	    if (os.access(sci_fits_co,os.F_OK)):
		os.remove( sci_fits_co )
	    if (os.access(sci_fits_co_simple,os.F_OK)):
		os.remove( sci_fits_co_simple )
	    hdu = pyfits.PrimaryHDU( sci_S_co )
	    hdu.writeto( sci_fits_co )
	    hdu = pyfits.PrimaryHDU( sci_Ss_co )
	    hdu.writeto( sci_fits_co_simple )

    else:
        print fsim, "has already been extracted, reading in product fits files..."
        sci_S_ob  = pyfits.getdata( sci_fits_ob )
        sci_Ss_ob = pyfits.getdata( sci_fits_ob_simple )
    p_shift = 0.
    if fsim in sci_sim:
            sci_S_co  = pyfits.getdata( sci_fits_co )
            sci_Ss_co = pyfits.getdata( sci_fits_co_simple )

	    lines_thar_co  = sci_S_co[:,1,:]
	    iv_thar_co     = sci_S_co[:,2,:]
	    offset_co = GLOBALutils.get_rough_offset(lines_thar_co,reffiles)
	    All_Pixel_Centers_co = np.array([])
	    All_Wavelengths_co   = np.array([])
	    All_Orders_co        = np.array([])
	    All_Centroids_co     = np.array([])
	    All_Sigmas_co        = np.array([])
	    All_Intensities_co   = np.array([])
	    All_residuals_co     = np.array([])

	    for order in range(nord_co):
		    order_s = str(order)
		    if (order < 10):
		        order_s = '0'+str(order)
		    thar_order_orig = lines_thar_co[order,:]
		    IV              = iv_thar_co[order,:]
		    wei             = np.sqrt( IV )
		    #bkg             = FEROSutils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)
		    thar_order      = thar_order_orig #- bkg

		    coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
		        = GLOBALutils.Initial_Wav_Calibration( order_dir+'fideos_'+order_s+'.iwdat', thar_order, order, wei, \
		                                           rmsmax=10000, minlines=10,FixEnds=False,Dump_Argon=False,Cheby=True,porder=ncoef_x,rough_shift = offset_co)
		    #psh, pix_centers, wavelengthss, rms_mss, residualss  = FEROSutils.Wav_Solution_vel_shift(wsol_dict['c_p2w_c'][order-o0], \
		    #								pixel_centers, wavelengths, maxrms=100, minlines=30, Cheby=use_cheby)
		    #shifts.append((1e-6*psh[0])*299792458.0)
		    if (order == 25): 
		        Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 1023, len(thar_order) )

		    All_Pixel_Centers_co = np.append( All_Pixel_Centers_co, pixel_centers )
		    All_Wavelengths_co   = np.append( All_Wavelengths_co, wavelengths )
		    All_Orders_co        = np.append( All_Orders_co, np.zeros( len(pixel_centers) ) + order )
		    All_Centroids_co     = np.append( All_Centroids_co, centroids)
		    All_Sigmas_co        = np.append( All_Sigmas_co, sigmas)
		    All_Intensities_co   = np.append( All_Intensities_co, intensities )
		    All_residuals_co     = np.append( All_residuals_co, residuals )
		    order+=1

	    p0 = np.zeros( npar_wsol )
	    p0[0] =  (25+49) * Global_ZP 
	    p1_co, G_pix_co, G_ord_co, G_wav_co, II_co, rms_ms_co, G_res_co = \
		    GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co, All_Wavelengths_co, All_Orders_co,\
		                                        np.ones(All_Intensities_co.shape), p0, Cheby=True,\
		                                        maxrms=MRMS, Inv=Inverse_m,minlines=400,order0=49, \
							ntotal=nord_co,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

	    #shifts = np.array(shifts)
	    #shifts = FEROSutils.sigma_clip(shifts)
	    #the_sh   = np.around(shifts.mean(),1)
	    #error_sh = np.around(np.sqrt(np.var(shifts)/float(len(shifts)-1)),1)
	    #print 'Shifts (per order):', the_sh, '+-', error_sh

	    p_shift, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
			    GLOBALutils.Global_Wav_Solution_vel_shift(All_Pixel_Centers_co, All_Wavelengths_co, All_Orders_co,\
		                                        np.ones(All_Intensities_co.shape), wsol_dict['p1_co'],\
				                                           Cheby=True,Inv=True,maxrms=MRMS,minlines=minlines_glob_co,\
		                                                           order0=49,ntotal=nord_co,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

    fout = 'proc/'+ fsim.split('/')[-1][:-4]+'sp.fits'
    spec = np.zeros((11, nord_ob, data.shape[0]))
    equis = np.arange( data.shape[0] )        
    for order in range(nord_ob):
            m = order + 49
            chebs = GLOBALutils.Calculate_chebs(equis, m, Inverse=Inverse_m,order0=49,ntotal=nord_ob,npix=data.shape[1],nx=ncoef_x,nm=ncoef_m)
	    WavSol = lbary_ltopo * (1.0 + 1.0e-6*p_shift) * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wsol_dict['p1'],chebs,ncoef_x,ncoef_m)
            spec[0,order,:] = WavSol
            spec[1,order,:] = sci_S_ob[order,1, :]
	    #spec[1,order,:] = sci_Ss_ob[order, :]
            spec[2,order,:] = sci_S_ob[order,2, :]
            # Flat-fielded spectrum
            fn = S_flat_ob_n[order,1,:]
            L  = np.where( fn > 0 )
            spec[3,order,:][L] = sci_S_ob[order,1,:][L] / S_flat_ob_n[order,1,:][L]
            spec[4,order,:][L] = sci_S_ob[order,2,:][L] * ( S_flat_ob_n[order,1,:][L] ** 2 )
	    plot(spec[0,order,:],spec[1,order,:])
    show()

    ccoefs = GLOBALutils.get_cont(spec[0,:,:],spec[3,:,:])
    for order in range(nord_ob):
	L  = np.where( spec[1,order] != 0 )
        spec[5,order,:][L] = spec[3,order][L] / np.polyval(ccoefs[order],spec[0,order][L])
	ratio = np.polyval(ccoefs[order],spec[0,order][L])
        spec[6,order,:][L] = spec[4,order][L] * (ratio ** 2 )
        spec[7,order,:][L] = ratio
        spec[8,order,:][L] = ratio * S_flat_ob_n[order,1][L] / np.sqrt( ratio * S_flat_ob_n[order,1][L] / gain + (ronoise/gain)**2 )
	plot(spec[0,order,:],spec[8,order,:])
    show()	
    hdu = pyfits.PrimaryHDU( spec )
    if (os.access( dirout + fout,os.F_OK)):
        os.remove( dirout + fout)
    hdu.writeto( dirout + fout )
