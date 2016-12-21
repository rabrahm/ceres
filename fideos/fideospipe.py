import sys
from pylab import *
base = '../'

sys.path.append(base+"utils/Correlation")
sys.path.append(base+"utils/OptExtract")
sys.path.append(base+"utils/GLOBALutils")

baryc_dir= base+'utils/SSEphem/'
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
force_stellar_pars = False

dumpargon          = False
minlines_glob_ob   = 400
minlines_glob_co   = 400

Inverse_m          = True
use_cheby          = True
MRMS               = 100   # max rms in m/s, global wav solution

order0		   = 49
ncoef_x		   = 4
ncoef_m		   = 6
npar_wsol          = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2
trace_degree       = 4
Marsh_alg          = 0
ext_aperture_ob    = 7
ext_aperture_co    = 6
NSigma_Marsh       = 10
NCosmic_Marsh      = 10
S_Marsh            = 0.4
N_Marsh            = 3      # grado polinomio 
min_extract_col    = 0
max_extract_col    = 2048
n_useful           = 70    # up to which order do we care?

order_dir   = "wavcals/"
models_path = base+"data/COELHO_MODELS/R_40000b/"    # path to the synthetic models 

RO_ob, GA_ob = 6.82,1.5
ronoise, gain = 6.82,1.5

log = dirout+'night.log'

print "\n\n\tFIDEOS ESO1.0m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'

"""
if os.access(dirin + 'log.txt', os.F_OK) and 1==2:
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
	biases   = glob.glob(dirin + 'bias*fits')
	flats    = glob.glob(dirin + 'flat*2F*fits')
    flats_co = glob.glob(dirin + 'flat*1F*fits')
	thars    = glob.glob(dirin + 'ThAr*fits')
	sci_sim	 = glob.glob(dirin + 'test_interfaz_ThAr*fits')
	#flat_ob = glob.glob(dirin + 'FIDEOS_flat_dark*fits')[0]
	#flat_co = glob.glob(dirin + 'FIDEOS_dark_flat*fits')[0]
	#thars   = glob.glob(dirin + 'FIDEOS_thar_thar*fits')[:2]
	#thars   = glob.glob(dirin + 'ThAr_sc*fit')[:2]
"""

biases, darks, flats, flats_co, thars, thars_co, sim_sci, dar_sci = fideosutils.FileClassify(dirin,log)

if ( (os.access(dirout+'MasterFlat.fits',os.F_OK) == False) or \
     (os.access(dirout+'MasterFlat_co.fits',os.F_OK) == False) or \
     (os.access(dirout+'MasterFlat_ob.fits',os.F_OK) == False) or \
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
    print "\t\t-> Masterbias: done!"

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

    MasterFlat_ob = fideosutils.get_flatOB(MasterFlat,MasterFlat_co)
    hdu = pyfits.PrimaryHDU( MasterFlat_ob )
    if (os.access(dirout+'MasterFlat_ob.fits',os.F_OK)):
        os.remove(dirout+'MasterFlat_ob.fits')
    hdu.writeto(dirout+'MasterFlat_ob.fits')
    print "\t\t-> Masterflats: done!"

    #print gfd
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
    c_co, nord_co = GLOBALutils.get_them(MasterFlat_co.T, 4.,trace_degree,maxords=-1,nsigmas=5,mode=1)
    c_co, nord_co = fideosutils.good_orders(c_co,nord_co,MasterFlat.shape[1],MasterFlat.shape[0],ext_aperture_co)
    print "\t\t-> Comparison orders found:", nord_co
    temp_flat = fideosutils.make_flatOB(MasterFlat, c_co)
    for i in range(temp_flat.shape[1]):
        temp_flat[:,i] = scipy.ndimage.filters.gaussian_filter(temp_flat[:,i],3.)
    c_ob, nord_ob = GLOBALutils.get_them(temp_flat, 10,trace_degree,maxords=-1,nsigmas=3,mode=1)
    c_ob, nord_ob = fideosutils.good_orders(c_ob,nord_ob,MasterFlat.shape[1],MasterFlat.shape[0],ext_aperture_ob)

    print "\t\t-> Object orders found:", nord_ob
    #c_ob, c_co, nord_ob, nord_co = fideosutils.clean_orders(c_all,MasterFlat.T)
    #print nord_ob, nord_co
    #c_co, nord_co = fideosutils.good_orders(c_co,nord_co,MasterFlat.shape[1],MasterFlat.shape[0],ext_aperture_co)
    #c_ob, nord_ob = fideosutils.good_orders(c_ob,nord_ob,MasterFlat.shape[1],MasterFlat.shape[0],ext_aperture_ob)
    #print nord_ob, nord_co
    #print nord_ob, 'object orders found'
    #print nord_co, 'comparison orders found'
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
    MasterFlat_co = pyfits.getdata(dirout+'MasterFlat_co.fits')
    MasterFlat_ob = pyfits.getdata(dirout+'MasterFlat_ob.fits')

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

P_ob_fits = dirout + 'P_ob.fits'
P_co_fits = dirout + 'P_co.fits'
S_flat_ob_fits = dirout +'S_flat_ob.fits'
S_flat_ob_simple_fits = dirout +'S_flat_ob_simple.fits'
S_flat_co_fits = dirout +'S_flat_co.fits'
S_flat_co_simple_fits = dirout +'S_flat_co_simple.fits'

if ( os.access(P_ob_fits,os.F_OK) == False ) or ( os.access(S_flat_ob_fits,os.F_OK) == False ) or \
    ( os.access(S_flat_ob_simple_fits,os.F_OK) == False ) or ( os.access(P_co_fits,os.F_OK) == False ) or \
    ( os.access(S_flat_co_fits,os.F_OK) == False ) or ( os.access(S_flat_co_simple_fits,os.F_OK) == False ) or \
    (force_flat_extract):

    ejx = np.arange(MasterFlat.shape[0])

    Centers = np.zeros((nord_ob,MasterFlat.shape[0]))
    for i in range(nord_ob):
        Centers[i,:]=scipy.polyval(c_ob[i],ejx)
    bac_ob = GLOBALutils.get_scat(MasterFlat_ob.T,Centers,span=8)

    Centers = np.zeros((nord_co,MasterFlat.shape[0]))
    #imshow(MasterFlat_co.T)
    for i in range(nord_co):
        Centers[i,:]=scipy.polyval(c_co[i],ejx)
    bac_co = GLOBALutils.get_scat(MasterFlat_co.T,Centers,span=8)

    Flat_ob = MasterFlat_ob.T - bac_ob
    Flat_co = MasterFlat_co.T - bac_co
    print 'P_ob...'
    P_ob = GLOBALutils.obtain_P(Flat_ob, c_ob, ext_aperture_ob, RO_ob,\
                                    GA_ob,NSigma_Marsh, S_Marsh, \
				    N_Marsh, Marsh_alg, min_extract_col,\
				    max_extract_col, npools)
    print 'P_co...'

    P_co = GLOBALutils.obtain_P(Flat_co, c_co, ext_aperture_co, RO_ob,\
                                    GA_ob,NSigma_Marsh, S_Marsh, \
				    N_Marsh, Marsh_alg, min_extract_col,\
				    max_extract_col, npools)
    print 'S_ob...'

    S_flat_ob_simple = GLOBALutils.simple_extraction(Flat_ob,c_ob,ext_aperture_ob,min_extract_col,\
                                                         max_extract_col,npools)
    S_flat_ob  = GLOBALutils.optimal_extraction(Flat_ob,P_ob,c_ob,ext_aperture_ob,RO_ob,GA_ob,S_Marsh,\
                                                    NCosmic_Marsh,min_extract_col,max_extract_col,npools)
    print 'S_co...'

    S_flat_co_simple = GLOBALutils.simple_extraction(Flat_co,c_co,ext_aperture_co,min_extract_col,\
                                                         max_extract_col,npools)
    S_flat_co  = GLOBALutils.optimal_extraction(Flat_co,P_co,c_co,ext_aperture_co,RO_ob,GA_ob,S_Marsh,\
                                                    NCosmic_Marsh,min_extract_col,max_extract_col,npools)

    S_flat_ob_simple = GLOBALutils.invert(S_flat_ob_simple)
    S_flat_ob = GLOBALutils.invert(S_flat_ob)

    S_flat_co_simple = GLOBALutils.invert(S_flat_co_simple)
    S_flat_co = GLOBALutils.invert(S_flat_co)

    if (os.access(P_ob_fits,os.F_OK)):
        os.remove( P_ob_fits )
    if (os.access(P_co_fits,os.F_OK)):
        os.remove( P_co_fits )

    hdu = pyfits.PrimaryHDU( P_ob )
    hdu.writeto( P_ob_fits )
    hdu = pyfits.PrimaryHDU( P_co )
    hdu.writeto( P_co_fits )

    if (os.access(S_flat_ob_fits,os.F_OK)):
        os.remove( S_flat_ob_fits )
    if (os.access(S_flat_ob_simple_fits,os.F_OK)):
        os.remove( S_flat_ob_simple_fits )

    hdu = pyfits.PrimaryHDU( S_flat_ob )
    hdu.writeto( S_flat_ob_fits )
    hdu = pyfits.PrimaryHDU( S_flat_ob_simple )
    hdu.writeto( S_flat_ob_simple_fits )

    if (os.access(S_flat_co_fits,os.F_OK)):
        os.remove( S_flat_co_fits )
    if (os.access(S_flat_co_simple_fits,os.F_OK)):
        os.remove( S_flat_co_simple_fits )

    hdu = pyfits.PrimaryHDU( S_flat_co )
    hdu.writeto( S_flat_co_fits )
    hdu = pyfits.PrimaryHDU( S_flat_co_simple )
    hdu.writeto( S_flat_co_simple_fits )

else:
    print "\tExtracted flat comparison spectra found, loading...\n"
    P_ob       = pyfits.getdata( P_ob_fits )
    S_flat_ob  = pyfits.getdata( S_flat_ob_fits )
    S_flat_ob_simple = pyfits.getdata( S_flat_ob_simple_fits )
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
#thars = thars[:3]
thar_times = []
for fsim in thars:
    hthar = pyfits.open( fsim )
    dthar = pyfits.getdata( fsim ) - MasterBias
    dthar = np.fliplr(dthar.T)
    mjd, mjd0 = fideosutils.mjd_fromheader( hthar )
    thar_times.append(mjd)
    thar_fits_ob_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.ob.fits.S'
    thar_fits_co_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.co.fits.S'
    thar_fits_ob = dirout + fsim.split('/')[-1][:-4]+'spec.ob.fits.S'
    thar_fits_co = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
    bac_fits_thar = dirout + 'BKG_'+ fsim.split('/')[-1][:-4]+'.fits'
    if ( os.access(thar_fits_ob_simple,os.F_OK) == False ) or ( os.access(thar_fits_co_simple,os.F_OK) == False ) or (force_thar_extract):
        print "No previous extraction or extraction forced for ThAr file", fsim, "extracting..."
        bac = fideosutils.get_scat(dthar.astype('float'), P_ob.astype('float') + P_co.astype('float'),1)
        dthar -= bac
        thar_Ss_ob_simple  = GLOBALutils.simple_extraction(dthar.T,c_ob,ext_aperture_ob,min_extract_col,\
                                                         max_extract_col,npools)
        thar_Ss_ob  = GLOBALutils.optimal_extraction(dthar.T,P_ob,c_ob,ext_aperture_ob,RO_ob,GA_ob,S_Marsh,\
                                                    NCosmic_Marsh,min_extract_col,max_extract_col,npools)
        thar_Ss_co_simple  = GLOBALutils.simple_extraction(dthar.T,c_co,ext_aperture_co,min_extract_col,\
                                                         max_extract_col,npools)
        thar_Ss_co  = GLOBALutils.optimal_extraction(dthar.T,P_co,c_co,ext_aperture_co,RO_ob,GA_ob,S_Marsh,\
                                                    NCosmic_Marsh,min_extract_col,max_extract_col,npools)
	
        thar_Ss_ob_simple = GLOBALutils.invert(thar_Ss_ob_simple)
        thar_Ss_co_simple = GLOBALutils.invert(thar_Ss_co_simple)
        thar_Ss_ob = GLOBALutils.invert(thar_Ss_ob)
        thar_Ss_co = GLOBALutils.invert(thar_Ss_co)

        if (os.access(thar_fits_ob,os.F_OK)):
            os.remove( thar_fits_ob )
        hdu = pyfits.PrimaryHDU( thar_Ss_ob )
        hdu.writeto( thar_fits_ob )
        if (os.access(thar_fits_co,os.F_OK)):
            os.remove( thar_fits_co )
        hdu = pyfits.PrimaryHDU( thar_Ss_co )
        hdu.writeto( thar_fits_co )

        if (os.access(thar_fits_ob_simple,os.F_OK)):
            os.remove( thar_fits_ob_simple )
        hdu = pyfits.PrimaryHDU( thar_Ss_ob_simple )
        hdu.writeto( thar_fits_ob_simple )
        if (os.access(thar_fits_co_simple,os.F_OK)):
            os.remove( thar_fits_co_simple )
        hdu = pyfits.PrimaryHDU( thar_Ss_co_simple )
        hdu.writeto( thar_fits_co_simple )

thar_times = np.array(thar_times)
thars      = thars[np.argsort(thar_times)]
thar_times = thar_times[np.argsort(thar_times)]

thar_co_times = []
for fsim in thars_co:
    hthar = pyfits.open( fsim )
    dthar = pyfits.getdata( fsim ) - MasterBias
    dthar = np.fliplr(dthar.T)
    mjd, mjd0 = fideosutils.mjd_fromheader( hthar )
    thar_co_times.append(mjd)
    thar_fits_co_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.co.fits.S'
    thar_fits_co = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
    bac_fits_thar = dirout + 'BKG_'+ fsim.split('/')[-1][:-4]+'.fits'
    if ( os.access(thar_fits_co_simple,os.F_OK) == False ) or (force_thar_extract):
        print "No previous extraction or extraction forced for ThAr file", fsim, "extracting..."
        bac = fideosutils.get_scat(dthar.astype('float'), P_co.astype('float'),1)
        dthar -= bac
        thar_Ss_co_simple  = GLOBALutils.simple_extraction(dthar.T,c_co,ext_aperture_co,min_extract_col,\
                                                         max_extract_col,npools)
        thar_Ss_co  = GLOBALutils.optimal_extraction(dthar.T,P_co,c_co,ext_aperture_co,RO_ob,GA_ob,S_Marsh,\
                                                    NCosmic_Marsh,min_extract_col,max_extract_col,npools)
    
        thar_Ss_co_simple = GLOBALutils.invert(thar_Ss_co_simple)
        thar_Ss_co = GLOBALutils.invert(thar_Ss_co)

        if (os.access(thar_fits_co,os.F_OK)):
            os.remove( thar_fits_co )
        hdu = pyfits.PrimaryHDU( thar_Ss_co )
        hdu.writeto( thar_fits_co )

        if (os.access(thar_fits_co_simple,os.F_OK)):
            os.remove( thar_fits_co_simple )
        hdu = pyfits.PrimaryHDU( thar_Ss_co_simple )
        hdu.writeto( thar_fits_co_simple )

thar_co_times = np.array(thar_co_times)
thars_co      = thars_co[np.argsort(thar_co_times)]
thar_co_times = thar_co_times[np.argsort(thar_co_times)]

print '\tSelection of echelle orders...'
force_thar_id = False
if os.access(dirout+'order_id.pkl', os.F_OK) == False or force_thar_id:
    print '\t\tComputing CCF...'
    thar_fits_ob = dirout + thars[0].split('/')[-1][:-4]+'spec.ob.fits.S'
    thar_ob = pyfits.getdata(thar_fits_ob)[:,1,:]
    or20_ob, offset_ob = GLOBALutils.identify_order(thar_ob,order_dir+'fideos_20.iwdat',window=100,di=0.5)
    start_ob = or20_ob - 9
    end_ob   = start_ob + 38

    print '\t\tWill start on order', start_ob, 'for object fibre...'
    print '\t\tOffset for object fibre is ', offset_ob, 'pixels'
    thar_fits_co = dirout + thars[0].split('/')[-1][:-4]+'spec.co.fits.S'
    thar_co = pyfits.getdata(thar_fits_co)[:,1,:]
    or20_co, offset_co = GLOBALutils.identify_order(thar_co,order_dir+'fideos_20.iwdat',window=100,di=0.5)
    start_co = or20_co - 9
    end_co   = start_co + 38

    print '\t\tWill start on order', start_co, 'for comparison fibre...'
    print '\t\tOffset for comparison fibre is ', offset_co, 'pixels'
    thar_dict = {'start_ob':start_ob,'end_ob':end_ob,'start_co':start_co,\
      'end_co':end_co, 'offset_ob':offset_ob, 'offset_co':offset_co}
    pickle.dump( thar_dict, open( dirout+'order_id.pkl', 'w' ) )
else:
    print '\t\tValues loaded from file...'
    thar_dict = pickle.load( open( dirout+'order_id.pkl', 'r' ) )
    start_ob    = thar_dict['start_ob']
    end_ob      = thar_dict['end_ob']
    offset_ob   = thar_dict['offset_ob']
    start_co    = thar_dict['start_co']
    end_co      = thar_dict['end_co']
    offset_co   = thar_dict['offset_co']
    print '\t\tWill start on order', start_ob, 'for object fibre...'
    print '\t\tOffset for object fibre is ', offset_ob, 'pixels'
    print '\t\tWill start on order', start_co, 'for comparison fibre...'
    print '\t\tOffset for comparison fibre is ', offset_co, 'pixels'
ntot = end_ob - start_ob +1
ntot2 = end_co - start_co +1
print 'ntot',ntot,ntot2
for fsim in thars:
    thar_fits_ob_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.ob.fits.S'
    thar_fits_co_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.co.fits.S'
    thar_spec_ob	= dirout + fsim.split('/')[-1][:-4]+'sp.ob.fits'
    thar_spec_co	= dirout + fsim.split('/')[-1][:-4]+'sp.co.fits'
    wavsol_pkl = dirout + fsim.split('/')[-1][:-4]+'wavsolpars.pkl'
    if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
        print "Working on initial ThAr file", fsim
        hthar = pyfits.open( fsim )
        mjd, mjd0 = fideosutils.mjd_fromheader( hthar )
        thar_fits_ob = dirout + fsim.split('/')[-1][:-4]+'spec.ob.fits.S'
        thar_S_ob = pyfits.getdata( thar_fits_ob )
        thar_fits_co = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
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
        All_Pixel_Centers_co = np.array([])
        All_Wavelengths_co   = np.array([])
        All_Orders_co        = np.array([])
        All_Centroids_co     = np.array([])
        All_Sigmas_co        = np.array([])
        All_Intensities_co   = np.array([])
	
        #offset_ob = GLOBALutils.get_rough_offset(lines_thar_ob,reffiles,window=300)
        #offset_co = GLOBALutils.get_rough_offset(lines_thar_co,reffiles,window=300)
        t_ords,t_wavs = [],[]
        order = start_ob
        idorder = 11
        while order <= end_ob:
            order_s = str(idorder)
            if (idorder < 10):
                order_s = '0'+str(idorder)
            thar_order_orig = lines_thar_ob[order,:]
            IV              = iv_thar_ob[order,:]
            wei             = np.ones(len(thar_order_orig))#np.sqrt( IV )
            bkg             = scipy.signal.medfilt(thar_order_orig,101)        
            thar_order      = thar_order_orig - bkg

            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
                = GLOBALutils.Initial_Wav_Calibration( order_dir+'fideos_'+order_s+'.iwdat', thar_order, order, wei, \
                                                   rmsmax=400, sigmai=2.2,minlines=10,FixEnds=False,Dump_Argon=False,\
                                                   Cheby=True,porder=ncoef_x,rough_shift=offset_ob,line_width=6,do_xc=False)
            #plot(pixel_centers,residuals,'.')
            #plot([0,2048],[0,0],'k')
            #show()
            if (idorder == 25): 
                Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 1023, len(thar_order) )

            All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
            All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
            All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + idorder )
            All_Centroids     = np.append( All_Centroids, centroids)
            All_Sigmas        = np.append( All_Sigmas, sigmas)
            All_Intensities   = np.append( All_Intensities, intensities )

            t_ords.append(order)
            t_wavs.append(GLOBALutils.Cheby_eval( coeffs_pix2wav, 1023, len(thar_order) ))
            order += 1
            idorder += 1
        t_ords,t_wavs = np.array(t_ords),np.array(t_wavs)

        #GLOBALutils.get_zero_order_number(t_ords,t_wavs)
        p0 = np.zeros( npar_wsol )
        p0[0] =  (25+49) * Global_ZP 
        p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
            GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders,\
                                                np.ones(All_Intensities.shape), p0, Cheby=True,\
                                                maxrms=MRMS, Inv=Inverse_m,minlines=minlines_glob_ob,order0=49, \
						ntotal=ntot,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
        #for o in np.unique(G_ord):
        #    I = np.where(G_ord == o)[0]
        #    plot(G_wav[I],G_res[I],'.')
        #    plot(np.median(G_wav[I]),np.median(G_res[I]),'ko')
        #show()
        thar_out = np.zeros((2,ntot,lines_thar_ob.shape[1]))
        equis = np.arange( lines_thar_ob.shape[1] )        
        order = start_ob
        idorder = 1
        out_order = 0
        while order <= end_ob:
            m   = idorder + 49
            chebs = GLOBALutils.Calculate_chebs(equis, m, Inverse=Inverse_m,order0=49,ntotal=ntot,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
            WavSol = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(p1,chebs,ncoef_x,ncoef_m)
            thar_out[0,out_order,:] = WavSol
            thar_out[1,out_order,:] = lines_thar_ob[order]
            order+=1
            idorder += 1
            out_order += 1
        
        if os.access(thar_spec_ob,os.F_OK):
            os.system('rm '+ thar_spec_ob)
        hdu = pyfits.PrimaryHDU(thar_out)
        hdu.writeto(thar_spec_ob)

        order = start_co
        idorder = 11
        while order <= end_co:
            order_s = str(idorder)
            if (idorder < 10):
                order_s = '0'+str(idorder)
            
            thar_order_orig = lines_thar_co[order,:]
            IV              = iv_thar_co[order,:]
            wei             = np.ones(len(thar_order_orig))#np.sqrt( IV )
            bkg             = scipy.signal.medfilt(thar_order_orig,101)        
            thar_order      = thar_order_orig - bkg
            IO = np.where(G_ord == idorder)[0]
            #print idorder
            #print pixel_centers
            #print wavelengths
            #print order_s
            #wavelengths, pixel_centers, intensities, sigmas, centroids = GLOBALutils.fit_these_lines(G_wav[IO], \
            #    order_dir+'fideos_'+order_s+'.iwdat', thar_order, order,wei, rough_shift=offset_co, \
            #    line_width=6, do_xc=False)

            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
                = GLOBALutils.Initial_Wav_Calibration( order_dir+'fideos_'+order_s+'.iwdat', thar_order, order, wei, \
                                                   rmsmax=400, minlines=10,FixEnds=False,Dump_Argon=False,Cheby=True,porder=ncoef_x,rough_shift = offset_co,do_xc=False,line_width=6)
            #if (idorder == 25): 
            #    Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 1023, len(thar_order) )

            All_Pixel_Centers_co = np.append( All_Pixel_Centers_co, pixel_centers )
            All_Wavelengths_co   = np.append( All_Wavelengths_co, wavelengths )
            All_Orders_co        = np.append( All_Orders_co, np.zeros( len(pixel_centers) ) + idorder )
            All_Centroids_co     = np.append( All_Centroids_co, centroids)
            All_Sigmas_co        = np.append( All_Sigmas_co, sigmas)
            All_Intensities_co   = np.append( All_Intensities_co, intensities )
            order += 1
            idorder += 1

        #p0 = np.zeros( npar_wsol )
        #p0[0] =  (25+49) * Global_ZP 
        p1_co, G_pix_co, G_ord_co, G_wav_co, II_co, rms_ms_co, G_res_co = \
            GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co, All_Wavelengths_co, All_Orders_co,\
                                                np.ones(All_Intensities_co.shape), p1, Cheby=True,\
                                                maxrms=MRMS, Inv=Inverse_m,minlines=minlines_glob_co,order0=49, \
						ntotal=ntot,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
        #for o in np.unique(G_ord_co):
        #    I = np.where(G_ord_co == o)[0]
        #    plot(G_wav_co[I],G_res_co[I],'.')
        #    plot(np.median(G_wav_co[I]),np.median(G_res_co[I]),'ko')
        #show()
        thar_out_co = np.zeros((2,nord_co,lines_thar_co.shape[1]))
        equis = np.arange( lines_thar_co.shape[1] )        
        order = start_co
        idorder = 11
        out_order = 0
        while order <= end_co:
            m   = idorder + 49
            chebs = GLOBALutils.Calculate_chebs(equis, m, Inverse=Inverse_m,order0=49,ntotal=ntot,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
            WavSol = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(p1_co,chebs,ncoef_x,ncoef_m)
            thar_out_co[0,out_order,:] = WavSol
            thar_out_co[1,out_order,:] = lines_thar_co[order]
            order+=1
            idorder += 1
            out_order += 1

        if os.access(thar_spec_co,os.F_OK):
            os.system('rm '+ thar_spec_co)
        hdu = pyfits.PrimaryHDU(thar_out_co)
        hdu.writeto(thar_spec_co)
 
        pdict = {'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                 'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Orders':All_Orders, 'All_Sigmas':All_Sigmas,\
                 'All_Pixel_Centers': All_Pixel_Centers, 'All_Wavelengths':All_Wavelengths, \
                 'p1_co':p1_co, 'G_pix_co':G_pix_co, 'G_ord_co':G_ord_co, 'G_wav_co':G_wav_co, 'II_co':II_co, \
		 'rms_ms_co':rms_ms_co,'G_res_co':G_res_co, 'All_Centroids_co':All_Centroids_co, 'All_Orders_co':All_Orders_co,\
         'All_Sigmas_co':All_Sigmas_co, 'All_Pixel_Centers_co': All_Pixel_Centers_co, 'All_Wavelengths_co':All_Wavelengths_co}
        pickle.dump( pdict, open( wavsol_pkl, 'w' ) )

    else:
        print "Using previously computed wavelength solution in file",wavsol_pkl

print 'eyeeeeeee'
for fsim in thars_co:
    print fsim
    thar_fits_co_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.co.fits.S'
    thar_spec_co    = dirout + fsim.split('/')[-1][:-4]+'sp.co.fits'
    wavsol_pkl = dirout + fsim.split('/')[-1][:-4]+'wavsolpars.pkl'
    if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
        print "Working on initial ThAr file", fsim
        hthar = pyfits.open( fsim )
        mjd, mjd0 = fideosutils.mjd_fromheader( hthar )
        thar_fits_co = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
        thar_S_co = pyfits.getdata( thar_fits_co )

        lines_thar_co  = thar_S_co[:,1,:]
        iv_thar_co     = thar_S_co[:,2,:]

        All_Pixel_Centers_co = np.array([])
        All_Wavelengths_co   = np.array([])
        All_Orders_co        = np.array([])
        All_Centroids_co     = np.array([])
        All_Sigmas_co        = np.array([])
        All_Intensities_co   = np.array([])
    
        #offset_ob = GLOBALutils.get_rough_offset(lines_thar_ob,reffiles,window=300)
        #offset_co = GLOBALutils.get_rough_offset(lines_thar_co,reffiles,window=300)

        order = start_co
        idorder = 11
        while order <= end_co:
            order_s = str(idorder)
            if (idorder < 10):
                order_s = '0'+str(idorder)
            
            thar_order_orig = lines_thar_co[order,:]
            IV              = iv_thar_co[order,:]
            wei             = np.ones(len(thar_order_orig))#np.sqrt( IV )
            bkg             = scipy.signal.medfilt(thar_order_orig,101)        
            thar_order      = thar_order_orig - bkg

            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
                = GLOBALutils.Initial_Wav_Calibration( order_dir+'fideos_'+order_s+'.iwdat', thar_order, order, wei, \
                                                   rmsmax=400, minlines=10,FixEnds=False,Dump_Argon=False,Cheby=True,porder=ncoef_x,rough_shift = offset_co,do_xc=False,line_width=6)
            if (idorder == 25): 
                Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 1023, len(thar_order) )

            All_Pixel_Centers_co = np.append( All_Pixel_Centers_co, pixel_centers )
            All_Wavelengths_co   = np.append( All_Wavelengths_co, wavelengths )
            All_Orders_co        = np.append( All_Orders_co, np.zeros( len(pixel_centers) ) + idorder )
            All_Centroids_co     = np.append( All_Centroids_co, centroids)
            All_Sigmas_co        = np.append( All_Sigmas_co, sigmas)
            All_Intensities_co   = np.append( All_Intensities_co, intensities )
            order += 1
            idorder += 1

        p0 = np.zeros( npar_wsol )
        p0[0] =  (25+49) * Global_ZP 
        p1_co, G_pix_co, G_ord_co, G_wav_co, II_co, rms_ms_co, G_res_co = \
            GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co, All_Wavelengths_co, All_Orders_co,\
                                                np.ones(All_Intensities_co.shape), p0, Cheby=True,\
                                                maxrms=MRMS, Inv=Inverse_m,minlines=minlines_glob_co,order0=49, \
                        ntotal=ntot,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

        thar_out_co = np.zeros((2,nord_co,lines_thar_co.shape[1]))
        equis = np.arange( lines_thar_co.shape[1] )        
        order = start_co
        idorder = 11
        out_order = 0
        while order <= end_co:
            m   = idorder + 49
            chebs = GLOBALutils.Calculate_chebs(equis, m, Inverse=Inverse_m,order0=49,ntotal=ntot,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)
            WavSol = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(p1_co,chebs,ncoef_x,ncoef_m)
            thar_out_co[0,out_order,:] = WavSol
            thar_out_co[1,out_order,:] = lines_thar_co[order]
            order+=1
            idorder += 1
            out_order += 1

        if os.access(thar_spec_co,os.F_OK):
            os.system('rm '+ thar_spec_co)
        hdu = pyfits.PrimaryHDU(thar_out_co)
        hdu.writeto(thar_spec_co)
 
        pdict = {'mjd':mjd, 'p1_co':p1_co, 'G_pix_co':G_pix_co, 'G_ord_co':G_ord_co, \
                 'G_wav_co':G_wav_co, 'II_co':II_co, 'rms_ms_co':rms_ms_co,'G_res_co':G_res_co, \
                 'All_Orders_co':All_Orders_co, 'All_Pixel_Centers_co': All_Pixel_Centers_co, 'All_Wavelengths_co':All_Wavelengths_co }
        pickle.dump( pdict, open( wavsol_pkl, 'w' ) )

    else:
        print "Using previously computed wavelength solution in file",wavsol_pkl

print 'Reference wavelength image will be:', thars[0]
pkl_wsol = dirout + thars[int(.5*len(thars))].split('/')[-1][:-4]+'wavsolpars.pkl'
print " Unpickling wavelength solution from", pkl_wsol, " ..."
wsol_dict = pickle.load(open(pkl_wsol,'r'))
sci_all = np.hstack((sim_sci, dar_sci))
#sci_all = sci_one.copy()
thar_fits_co = dirout + thars[int(.5*len(thars))].split('/')[-1][:-4]+'spec.co.fits.S'
thar_S_co = pyfits.getdata( thar_fits_co )

print '\nComputing instrumental drifts for ThAr + ThAr imeges...'
drift_file = dirout + '/drifts.txt'
force_drifts = True
drift_out = []
for i in range(len(thars)):
        texp = pyfits.getheader(thars[i])['EXPTIME']
        pkl  = dirout + thars[i].split('/')[-1][:-4]+'wavsolpars.pkl'
        wsol = pickle.load(open(pkl,'r'))

        p_shift_co, pix_centers, orders, wavelengths, I, rms_ms_co, residuals_co  = \
                    GLOBALutils.Global_Wav_Solution_vel_shift(wsol['All_Pixel_Centers_co'], wsol['All_Wavelengths_co'],\
                    wsol['All_Orders_co'], np.ones(len(wsol['All_Orders_co'])), wsol_dict['p1_co'], Cheby=True, \
                    Inv=True,maxrms=MRMS,minlines=minlines_glob_co, order0=49,ntotal=ntot,npix=thar_S_co.shape[2],nx=ncoef_x,nm=ncoef_m)
        p_shift_ob, pix_centers, orders, wavelengths, I, rms_ms_ob, residuals_ob  = \
                    GLOBALutils.Global_Wav_Solution_vel_shift(wsol['All_Pixel_Centers'], wsol['All_Wavelengths'],\
                    wsol['All_Orders'], np.ones(len(wsol['All_Orders'])), wsol_dict['p1'], Cheby=True, \
                    Inv=True,maxrms=MRMS,minlines=minlines_glob_ob, order0=49,ntotal=ntot,npix=thar_S_co.shape[2],nx=ncoef_x,nm=ncoef_m)
        print '\n'
        mat = np.array([thar_times[i],p_shift_ob,rms_ms_ob/np.sqrt(float(len(residuals_ob))), p_shift_co, rms_ms_co/np.sqrt(float(len(residuals_co))),texp])
        if len(drift_out)==0:
            drift_out = mat.copy()
        else:
            drift_out = np.vstack((drift_out,mat))
np.savetxt(drift_file,drift_out)

print '\nComputing instrumental drifts for ThAr + Dark imeges...'

for i in range(len(thars_co)):
        texp = pyfits.getheader(thars_co[i])['EXPTIME']
        pkl  = dirout + thars_co[i].split('/')[-1][:-4]+'wavsolpars.pkl'
        wsol = pickle.load(open(pkl,'r'))

        p_shift_co, pix_centers, orders, wavelengths, I, rms_ms_co, residuals_co  = \
                    GLOBALutils.Global_Wav_Solution_vel_shift(wsol['G_pix_co'], wsol['G_wav_co'],\
                    wsol['G_ord_co'], np.ones(len(wsol['G_ord_co'])), wsol_dict['p1_co'], Cheby=True, \
                    Inv=True,maxrms=MRMS,minlines=minlines_glob_co, order0=49,ntotal=ntot,npix=thar_S_co.shape[2],nx=ncoef_x,nm=ncoef_m)
        mat = np.array([thar_times[i],p_shift_co,rms_ms_co/np.sqrt(float(len(residuals_co))), p_shift_co, rms_ms_co/np.sqrt(float(len(residuals_co))),texp])
        if len(drift_out)==0:
            drift_out = mat.copy()
        else:
            drift_out = np.vstack((drift_out,mat))
np.savetxt(drift_file,drift_out)

### start of science frame reductions ###

new_list = []
new_list_obnames = []
new_list_texp = []
for i in range(len(sci_all)):
    fsim = sci_all[i]
    hd = pyfits.getheader(fsim)
    texp = hd['EXPTIME']
    obname = fsim.split('/')[-1]
    obname = obname.split('-')[-1][:-5]

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
for fsim in new_list:
    print '\n'
    print "\t--> Working on image: ", fsim
    hd = pyfits.getheader(fsim)

    hsim = pyfits.open(fsim)
    mjd, mjd0 = fideosutils.mjd_fromheader( hsim )

    exptime = hd['EXPTIME']
    obname = fsim.split('/')[-1]
    obname = obname.split('-')[-1][:-5]
    print "\t\tObject name:",obname

    know_moon = False
    if fsim.split('/')[-1] in spec_moon:
        I = np.where(fsim.split('/')[-1] == spec_moon)[0]
        know_moon = True
        here_moon = use_moon[I]

    altitude    =  2335.
    latitude    = -29.2543
    longitude   = -70.7346
    epoch       =  2000.0

    known_coords = False    
    sp,ra,dec,known_coords = GLOBALutils.simbad_coords(obname,mjd)
    ras,decs = ra,dec
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
    jplephem.set_observer_coordinates( obpos[0], obpos[1], obpos[2] )

    res         = jplephem.doppler_fraction(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
    lbary_ltopo = 1.0 + res['frac'][0]
    bcvel_baryc = ( lbary_ltopo - 1.0 ) * 2.99792458E5

    print "\t\tBarycentric velocity:", bcvel_baryc

    res = jplephem.pulse_delay(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)   
    mbjd = mjd + res['delay'][0] / (3600.0 * 24.0)

    # Moon Phase Calculations
    gobs      = ephem.Observer()  
    gobs.name = 'ESO1.0'  
    gobs.lat  = rad(latitude)  # lat/long in decimal degrees  
    gobs.long = rad(longitude)

    date = hd['DATE-OBS']
    date = datetime.datetime(int(date[:4]),int(date[5:7]),int(date[8:10]),int(date[11:13]),int(date[14:16]),int(date[17:19]))
    #OJO aquiiiiii
    print 'Warning!!! adding 5 hrs to comute MJD due to problem in header! CHECK in future!!'
    new_date = date + datetime.timedelta(hours=5)

    gobs.date = new_date.strftime('%Y-%m-%d %H:%M:%S')

    mephem = ephem.Moon()
    mephem.compute(gobs)
    Mcoo = jplephem.object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Mp   = jplephem.barycentric_object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Sp   = jplephem.barycentric_object_track("Sun", int(mjd), float(mjd%1), 1, 0.0)
    res  = jplephem.object_doppler("Moon", int(mjd), mjd%1, 1, 0.0)
    lunation,moon_state,moonsep,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,ra,dec)
    refvel = bcvel_baryc + moonvel
    print '\t\tRadial Velocity of sacttered moonlight:',refvel

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
    if fsim in sim_sci:
        bac = fideosutils.get_scat(data, P_ob + P_co,1)
    else:
        bac = fideosutils.get_scat(data, P_ob,1)
    #imshow(bac)	
    #show()
    #plot(data[1000])
    #plot(bac[1000])
    #show()
    data -= bac
    #plot(data[1000])
    #show()
    sci_fits_ob_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.ob.fits.S'
    sci_fits_co_simple = dirout + fsim.split('/')[-1][:-4]+'spec.simple.co.fits.S'
    sci_fits_ob = dirout + fsim.split('/')[-1][:-4]+'spec.ob.fits.S'
    sci_fits_co = dirout + fsim.split('/')[-1][:-4]+'spec.co.fits.S'
    bac_fits    = dirout + 'BKG_'+ fsim.split('/')[-1][:-4]+'.fits'

    if ( os.access(sci_fits_ob,os.F_OK) == False ) or ( os.access(sci_fits_co,os.F_OK) == False ) or \
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

        if fsim in sim_sci:
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
    if fsim in sim_sci:
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

        order = start_co
        idorder = 11
        while order <= end_co:
            order_s = str(idorder)
            if (idorder < 10):
                order_s = '0'+str(idorder)
            thar_order_orig = lines_thar_co[order,:]
            IV              = iv_thar_co[order,:]
            wei             = np.sqrt( IV )
            bkg             = scipy.signal.medfilt(thar_order_orig,101) #FEROSutils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)
            thar_order      = thar_order_orig - bkg

            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
		        = GLOBALutils.Initial_Wav_Calibration( order_dir+'fideos_'+order_s+'.iwdat', thar_order, order, wei, \
		                                           rmsmax=400, minlines=10,FixEnds=False,Dump_Argon=False,Cheby=True,porder=ncoef_x,rough_shift = offset_co,do_xc=False,line_width=6)
		    #psh, pix_centers, wavelengthss, rms_mss, residualss  = FEROSutils.Wav_Solution_vel_shift(wsol_dict['c_p2w_c'][order-o0], \
		    #								pixel_centers, wavelengths, maxrms=100, minlines=30, Cheby=use_cheby)
		    #shifts.append((1e-6*psh[0])*299792458.0)
            if (idorder == 25): 
                Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 1023, len(thar_order) )

            All_Pixel_Centers_co = np.append( All_Pixel_Centers_co, pixel_centers )
            All_Wavelengths_co   = np.append( All_Wavelengths_co, wavelengths )
            All_Orders_co        = np.append( All_Orders_co, np.zeros( len(pixel_centers) ) + idorder )
            All_Centroids_co     = np.append( All_Centroids_co, centroids)
            All_Sigmas_co        = np.append( All_Sigmas_co, sigmas)
            All_Intensities_co   = np.append( All_Intensities_co, intensities )
            All_residuals_co     = np.append( All_residuals_co, residuals )
            order+=1
            idorder += 1

        p0 = np.zeros( npar_wsol )
        p0[0] =  (25+49) * Global_ZP 
 
        p1_co, G_pix_co, G_ord_co, G_wav_co, II_co, rms_ms_co, G_res_co = \
		    GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co, All_Wavelengths_co, All_Orders_co,\
		                                        np.ones(All_Intensities_co.shape), p0, Cheby=True,\
		                                        maxrms=MRMS, Inv=Inverse_m,minlines=minlines_glob_co,order0=49, \
							ntotal=ntot,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

        #shifts = np.array(shifts)
        #shifts = FEROSutils.sigma_clip(shifts)
        #the_sh   = np.around(shifts.mean(),1)
        #error_sh = np.around(np.sqrt(np.var(shifts)/float(len(shifts)-1)),1)
        #print 'Shifts (per order):', the_sh, '+-', error_sh

        p_shift, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
			    GLOBALutils.Global_Wav_Solution_vel_shift(All_Pixel_Centers_co, All_Wavelengths_co, All_Orders_co,\
		                                        np.ones(All_Intensities_co.shape), wsol_dict['p1_co'],\
				                                           Cheby=True,Inv=True,maxrms=MRMS,minlines=minlines_glob_co,\
		                                                           order0=49,ntotal=ntot,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

        mat = np.array([mjd,p_shift,rms_ms/np.sqrt(float(len(residuals))), p_shift, rms_ms/np.sqrt(float(len(residuals))),exptime])
        if len(drift_out)==0:
            drift_out = mat.copy()
        else:
            drift_out = np.vstack((drift_out,mat))
    np.savetxt(drift_file,drift_out)

    fout = 'proc/'+ fsim.split('/')[-1][:-4]+'sp.fits'
    spec = np.zeros((11, ntot, data.shape[0]))
    hdu = pyfits.PrimaryHDU( spec )

    hdu = GLOBALutils.update_header(hdu,'HIERARCH MJD', mjd)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH MBJD', mbjd)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START DATE', new_date.strftime('%Y-%m-%d') )
    hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START UT',  new_date.strftime('%H:%M:%S'))
    hdu = GLOBALutils.update_header(hdu,'HIERARCH TEXP (S)',hd['EXPTIME'])
    hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (KM/S)', bcvel_baryc)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH (LAMBDA_BARY / LAMBDA_TOPO)', lbary_ltopo)    
    hdu = GLOBALutils.update_header(hdu,'HIERARCH TARGET NAME', obname)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH RA',ras)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC',decs)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH RA DEG',ras)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC DEG',decs)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH RA BARY',ra)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC BARY',dec)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH EQUINOX',2000)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LATITUDE',latitude)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LONGITUDE',longitude)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS ALTITUDE',altitude)
    hdu = GLOBALutils.update_header(hdu,'HIERARCH TARG AIRMASS',hd['AIRMASSS'])

    equis = np.arange( data.shape[0] )
    out_order = 0
    order = start_ob
    idorder = 11

    while order <= end_ob:
        m = idorder + 49
        chebs = GLOBALutils.Calculate_chebs(equis, m, Inverse=Inverse_m,order0=49,ntotal=ntot,npix=len(equis),nx=ncoef_x,nm=ncoef_m)
        WavSol = lbary_ltopo * (1.0 + 1.0e-6*p_shift) * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wsol_dict['p1'],chebs,ncoef_x,ncoef_m)
        spec[0,out_order,:] = GLOBALutils.ToVacuum(WavSol)
        spec[1,out_order,:] = sci_S_ob[order,1, :]
        #spec[1,order,:] = sci_Ss_ob[order, :]
        spec[2,out_order,:] = sci_S_ob[order,2, :]
        # Flat-fielded spectrum
        fn = S_flat_ob_n[order,1,:]
        L  = np.where( fn > 0 )
        spec[3,out_order,:][L] = sci_S_ob[order,1,:][L] / S_flat_ob_n[order,1,:][L]
        spec[4,out_order,:][L] = sci_S_ob[order,2,:][L] * ( S_flat_ob_n[order,1,:][L] ** 2 )
        ccoef = GLOBALutils.get_cont_single(spec[0,out_order],spec[3,out_order],spec[4,out_order],ll=1.5,lu=5,nc=3)
        L  = np.where( spec[1,out_order] != 0 )
        spec[5,out_order,:][L] = spec[3,out_order][L] / np.polyval(ccoef,spec[0,out_order][L])

        nJ = np.where(np.isnan(spec[5,out_order])==True)[0]
        nJ2 = np.where(np.isinf(spec[5,out_order])==True)[0]
        spec[5,out_order,nJ] = 1.0
        spec[5,out_order,nJ2] = 1.0
        ratio            = spec[3,out_order,:][L] / spec[5,out_order,:][L]
        spec[6,out_order,:][L] = spec[4,out_order,:][L] * (ratio ** 2 )
        spec[7,out_order,:][L] = ratio
        spec[8,out_order,:][L] = ratio * S_flat_ob_n[order,1,:][L] / np.sqrt( ratio * S_flat_ob_n[order,1,:][L] / gain + (ronoise/gain)**2 )

        rI = np.where(spec[5,out_order] > 1. + 8./spec[8,out_order])
        spec[5,out_order,rI] = 1.
        #In = np.where(np.isinf(spec[5,out_order]))
        #print In
        spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
        dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
        NN            = np.average(dlambda_dx)
        dlambda_dx    /= NN

        spec[9,out_order,:][L]  = spec[5,out_order,:][L] * (dlambda_dx[L] ** 1) 
        spec[10,out_order,:][L] = spec[6,out_order,:][L] / (dlambda_dx[L] ** 2)
        #In1 = np.where(np.isinf(spec[8,out_order]))
        #In2 = np.where(np.isinf(spec[9,out_order]))
        #In3 = np.where(np.isinf(spec[10,out_order]))
        #print In1,In2,In3
        #print In

        order += 1
        out_order += 1
        idorder += 1

    if (os.access( dirout + fout,os.F_OK)):
        os.remove( dirout + fout)
    hdu.writeto( dirout + fout )

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
                spec2 = spec.copy()
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
        ml_v -= 1.5*(av_m - ml_v)
        mh_v += 1.5*(mh_v - av_m)
        mask_hw_kms = (GLOBALutils.Constants.c/1e3) * 0.5*(mh_v - ml_v) / av_m

        #sigma_fout = stellar_pars_dir + obname + '_' +'sigma.txt'

        disp = GLOBALutils.get_disp(obname, reffile=reffile)
        if disp == 0:
            known_sigma = False
            if vsini != -999 and vsini > 4.:
                disp = vsini
            else:
                disp = 4.
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
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight, 0, lbary_ltopo, vel_width=200,vel_step=3,\
                                          spec_order=9,iv_order=10,sn_order=8,max_vel_rough=200)
            v1,x1 = vels.copy(), xc_full.copy()
            xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=0.0, Simple=True, W=W_ccf, start_order=2)
            # Normalize the continuum of the CCF robustly with R     
            yy = scipy.signal.medfilt(xc_av,11)
            pred = lowess(yy, vels,frac=0.4,it=10,return_sorted=False)
            lowess_interp = scipy.interpolate.interp1d(vels,yy)
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

            xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=0.0, Simple=True, W=W_ccf, start_order=2)
            pred = lowess_interp(vels)
            xc_av /= pred
            #for i in range(xc_full.shape[1]):
            #    plot(v1,x1[:,i])
            #    plot(vels,xc_full[:,i])
            #    show()

        
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
                if (disp < 4.0): 
                    disp = 4.0
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
                 'lunation':lunation,'mephem':mephem,'texp':hd['EXPTIME']}

        pkl_xc = dirout + fsim.split('/')[-1][:-8]+obname+'_XC_'+sp_type+'.pkl'
        pickle.dump( xc_dict, open( pkl_xc, 'w' ) )

        ccf_pdf = dirout + 'proc/' + fsim.split('/')[-1][:-4] + obname + '_XCs_' + sp_type + '.pdf'

        if not avoid_plot:
            GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)

        SNR_5130 = np.median(spec[8,30,1000:1101] )
        airmass  = hd['AIRMASSS']
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

        #if not good_quality:
        #    RVerr2 = np.sqrt(0.03**2 + RVerr2**2)

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
        hdu = GLOBALutils.update_header(hdu,'INST', 'FIDEOS')
        hdu = GLOBALutils.update_header(hdu,'RESOL', '40000')
        hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')
        hdu = GLOBALutils.update_header(hdu,'XC_MIN', XC_min)
        hdu = GLOBALutils.update_header(hdu,'BJD_OUT', bjd_out)
        
        line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f   fideos   ceres   40000 %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
                      (obname, bjd_out, RV, RVerr2, BS, BSerr, T_eff_epoch, logg_epoch, Z_epoch, vsini_epoch, XC_min, disp_epoch,\
        hd['EXPTIME'], SNR_5130_R, ccf_pdf)
        f_res.write(line_out)
        if (os.access( dirout + fout,os.F_OK)):
            os.remove( dirout + fout)
        hdu.writeto( dirout + fout )

    else:
        print "\t\tReading spectral file from", fout
        spec = pyfits.getdata( fout )

f_res.close()
np.savetxt(drift_file,drift_out)

