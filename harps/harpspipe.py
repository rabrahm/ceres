import sys
from pylab import *

base = '../'

sys.path.append(base+"utils/Correlation")
sys.path.append(base+"utils/OptExtract")
sys.path.append(base+"utils/GLOBALutils")

baryc_dir= base+'utils/SSEphem/'
sys.path.append(baryc_dir)
ephemeris='DEc403'

import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ceres modules
import harpsutils
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
parser.add_argument('-mode', default='HARPS')

args = parser.parse_args()
dirin            = args.directorio
avoid_plot       = args.avoid_plot
dirout           = args.dirout
DoClass          = args.do_class
JustExtract      = args.just_extract
npools           = int(args.npools)
object2do        = args.o2do
reffile          = args.reffile
mode	         = args.mode

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
force_thar_wavcal  = False
force_tharxc       = False
force_sci_extract  = False
force_stellar_pars = False
force_spectral_file_build = True
dumpargon          = False
minlines_glob      = 1000
minlines_initial   = 50

Inverse_m          = True
use_cheby          = True
MRMS_initial       = 50  # max rms in m/s, initial wav solution
MRMS               = 20   # max rms in m/s, global wav solution

trace_degree       = 6
Marsh_alg          = 0
ext_aperture_R     = 5
ext_aperture_B     = 5
NSigma_Marsh       = 5
NCosmic_Marsh      = 5
S_Marsh            = 0.4
N_Marsh            = 3      # grado polinomio 
min_extract_col    = 0
max_extract_col    = 4095

porder      = 4
ncoef_x_B   = 5
ncoef_m_B   = 7
npar_wsol_B = (min(ncoef_x_B,ncoef_m_B) + 1) * (2*max(ncoef_x_B,ncoef_m_B) - min(ncoef_x_B,ncoef_m_B) + 2) / 2
ncoef_x_R   = 5
ncoef_m_R   = 6
npar_wsol_R = (min(ncoef_x_R,ncoef_m_R) + 1) * (2*max(ncoef_x_R,ncoef_m_R) - min(ncoef_x_R,ncoef_m_R) + 2) / 2

or0_R = 89
or0_B = 116

models_path = base+"data/COELHO_MODELS/R_40000b/"
order_dir   = "wavcals/"
final_wav   = '.iwdat'

RESI = 120000.
if mode=='EGGS':
    RESI         = 85000.
    MRMS_initial = 100.
    MRMS         = 50.
    final_wav    = '_eggs.iwdat'

#############################

print "\n\n\tHARPS ESO3.6m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'

# file containing the log
log = dirout+'night.log'

biases, flats, ob_flats, co_flats, ThAr_ref, sim_sci, ThAr_ref_dates, obnames, exptimes, co_types = harpsutils.FileClassify(dirin,log,mode=mode)

print '\tThis in the log of the night:\n'
f = open(log)
flines = f.readlines()
for line in flines:
	print '\t'+line[:-1]
print '\n'

if (     (os.access(dirout+'FlatOb_'+ mode +'.fits',os.F_OK) == False and len(ob_flats)!=0) or \
         (os.access(dirout+'FlatCo_'+ mode +'.fits',os.F_OK) == False and len(co_flats)!=0) or \
	 (os.access(dirout+'Flat_'+ mode +'.fits',  os.F_OK) == False and len(flats)!=0)    or \
         (os.access(dirout+'trace_'+ mode +'.pkl',os.F_OK) == False)  or \
         (os.access(dirout+'MasterBias_'+ mode +'.fits',os.F_OK) == False and len(biases)!=0)  or \
         (force_pre_process) ):
    print "\tNo previous pre-processing files or found"
    pre_process = 1
else:
    print "\tPre-processing files found, going straight to extraction"
    pre_process = 0

if (pre_process == 1):
    # median combine Biases
    print "\t\tGenerating Master calibration frames..."
    MasterBias, RO_bias, GA_bias = harpsutils.MedianCombine(biases)
    hdu = pyfits.PrimaryHDU( MasterBias )
    if (os.access(dirout+'MasterBias_'+mode+'.fits',os.F_OK)):
        os.remove(dirout+'MasterBias_'+mode+'.fits')
    hdu.writeto(dirout+'MasterBias_'+mode+'.fits')
    print "\t\t-> Masterbias: done!"

    # median combine list of ob flats
    Flat_ob, RO_ob, GA_ob = harpsutils.MedianCombine(ob_flats, zero=dirout+'MasterBias_'+mode+'.fits')
    # save this file for later reference
    hdu = pyfits.PrimaryHDU( Flat_ob )
    if (os.access(dirout+'FlatOb_'+ mode +'.fits',os.F_OK)):
        os.remove(dirout+'FlatOb_'+ mode +'.fits')
    hdu.writeto(dirout+'FlatOb_'+ mode +'.fits')

    # median combine list of co flats
    Flat_co, RO_co, GA_co = harpsutils.MedianCombine(co_flats, zero=dirout+'MasterBias_'+mode+'.fits')
    hdu = pyfits.PrimaryHDU(Flat_co)
    if (os.access(dirout+'FlatCo_'+ mode +'.fits',os.F_OK)):
        os.remove(dirout+'FlatCo_'+ mode +'.fits')
    hdu.writeto(dirout+'FlatCo_'+ mode +'.fits')

    # median combine list of flats
    Flat, RO_fl, GA_fl = harpsutils.MedianCombine(flats, zero=dirout+'MasterBias_'+mode+'.fits')
    hdu = pyfits.PrimaryHDU(Flat)
    if (os.access(dirout+'Flat_'+mode+'.fits',os.F_OK)):
        os.remove(dirout+'Flat_'+mode+'.fits')
    hdu.writeto(dirout+'Flat_'+mode+'.fits')
    print "\t\t-> Masterflat: done!"

    print "\tTracing echelle orders..."
    c_all1,nord_all1 = GLOBALutils.get_them(Flat[:,:,0],ext_aperture_B,trace_degree,mode=1)
    c_all1 = c_all1[5:]
    nord_all1 = len(c_all1) 
    I = np.arange(0,nord_all1,2).astype('int')
    c_ob1 = c_all1[I]
    nord_ob1 = len(I)
    I = np.arange(1,nord_all1,2).astype('int')
    c_co1 = c_all1[I]
    nord_co1 = len(I)
    
    c_all2,nord_all2 = GLOBALutils.get_them(Flat[:,:,1],ext_aperture_R,trace_degree,mode=1)
    I = np.arange(0,nord_all2,2).astype('int')
    c_ob2 = c_all2[I]
    nord_ob2 = len(I)
    I = np.arange(1,nord_all2,2).astype('int')
    c_co2 = c_all2[I]
    nord_co2 = len(I)

    print "\t\t"+str(nord_ob1)+" object orders found in blue chip..."
    print "\t\t"+str(nord_co1)+" comparison orders found in blue chip..."
    print "\t\t"+str(nord_ob2)+" object orders found in red chip..."
    print "\t\t"+str(nord_co2)+" comparison orders found in red chip..."
    # pickle traces
    trace_dict = {'c_ob1':c_ob1,
                  'c_co1':c_co1,
		  'c_ob2':c_ob2,
                  'c_co2':c_co2,
		  'c_all1':c_all1,
		  'c_all2':c_all2,
                  'nord_ob1':nord_ob1, 'nord_co1':nord_co1,
		  'nord_ob2':nord_ob2, 'nord_co2':nord_co2,
                  'GA_ob': GA_ob, 'RO_ob': RO_ob,
                  'GA_co': GA_co, 'RO_co': RO_co,
		  'GA_fl': GA_fl, 'RO_fl': RO_fl}
    pickle.dump( trace_dict, open( dirout+"trace_"+mode+".pkl", 'w' ) )

else:
    trace_dict = pickle.load( open( dirout+"trace_"+mode+".pkl", 'r' ) )
    c_co1 = trace_dict['c_co1']
    c_ob1 = trace_dict['c_ob1']
    c_co2 = trace_dict['c_co2']
    c_ob2 = trace_dict['c_ob2']
    nord_ob1 = trace_dict['nord_ob1']
    nord_co1 = trace_dict['nord_co1']
    nord_ob2 = trace_dict['nord_ob2']
    nord_co2 = trace_dict['nord_co2']
    # recover GA*, RO*
    GA_ob = trace_dict['GA_ob']
    RO_ob = trace_dict['RO_ob']
    GA_co = trace_dict['GA_co']
    RO_co = trace_dict['RO_co']
    GA_fl = trace_dict['GA_fl']
    RO_fl = trace_dict['RO_fl']
    # recover flats & master bias
    h = pyfits.open(dirout+'FlatOb_'+ mode +'.fits')
    Flat_ob = h[0].data
    h = pyfits.open(dirout+'FlatCo_'+ mode +'.fits')
    Flat_co = h[0].data
    h = pyfits.open(dirout+'Flat_'+mode+'.fits')
    Flat = h[0].data
    h = pyfits.open(dirout+'MasterBias_'+mode+'.fits')
    MasterBias = h[0].data
    
# mesh all orders
c_all1  = GLOBALutils.Mesh( c_ob1, c_co1)
c_all2  = GLOBALutils.Mesh( c_ob2, c_co2)

# Extract flat spectra, object
print '\n\tExtraction of Flat calibration frames:'
P_ob_B_fits    = dirout + 'P_ob_B_'+mode+'.fits'
P_ob_R_fits    = dirout + 'P_ob_R_'+mode+'.fits'
B_flat_ob_fits = dirout +'B_flat_ob_'+mode+'.fits'
R_flat_ob_fits = dirout +'R_flat_ob_'+mode+'.fits'

P_ob_B    = np.zeros( (Flat_ob.shape[0],Flat_ob.shape[1]) )
P_ob_R    = np.zeros( (Flat_ob.shape[0],Flat_ob.shape[1]) )
B_flat_ob = np.zeros((nord_ob1, 3, Flat_ob.shape[1]) )
R_flat_ob = np.zeros((nord_ob2, 3, Flat_ob.shape[1]) )

if ( os.access(P_ob_B_fits,os.F_OK) == False ) or ( os.access(B_flat_ob_fits,os.F_OK) == False ) or ( os.access(P_ob_R_fits,os.F_OK) == False ) or ( os.access(R_flat_ob_fits,os.F_OK) == False ) or (force_flat_extract):
    print "\t\tNo extracted flat object spectra found or extraction forced, extracting and saving..."

    bacfile = dirout + 'BACR_FLAT_'+mode+'.fits'
    if os.access(bacfile,os.F_OK) == False:
        CentersR = np.zeros((len(c_ob2),Flat_ob[:,:,1].shape[1]))
        for i in range(len(c_ob2)):
            CentersR[i,:]=np.polyval(c_ob2[i],np.arange(Flat_ob[:,:,1].shape[1]))
        bacR = GLOBALutils.get_scat(Flat_ob[:,:,1],CentersR,span=20)
        if (os.access(bacfile,os.F_OK)):
            os.remove( bacfile )
        hdbac = pyfits.PrimaryHDU( bacR )
        hdbac.writeto(bacfile)
    else:
        bacR = pyfits.getdata(bacfile)

    flR = Flat_ob[:,:,1] - bacR
    
    bacfile = dirout + 'BACB_FLAT_'+mode+'.fits'
    if os.access(bacfile,os.F_OK) == False:
        CentersB = np.zeros((len(c_ob1),Flat_ob[:,:,0].shape[1]))
        for i in range(len(c_ob1)):
            CentersB[i,:]=np.polyval(c_ob1[i],np.arange(Flat_ob[:,:,0].shape[1]))
        bacB = GLOBALutils.get_scat(Flat_ob[:,:,0],CentersB,span=10)
        if (os.access(bacfile,os.F_OK)):
            os.remove( bacfile )
        hdbac = pyfits.PrimaryHDU( bacB )
        hdbac.writeto(bacfile)
    else:
        bacB = pyfits.getdata(bacfile)

    flB = Flat_ob[:,:,0] - bacB

    P_ob_B = GLOBALutils.obtain_P(flB,c_ob1,ext_aperture_B,RO_fl[0],\
                                    GA_fl[0],NSigma_Marsh, S_Marsh, \
				    N_Marsh, Marsh_alg, min_extract_col,\
				    max_extract_col, npools)
    P_co_B = GLOBALutils.obtain_P(flB,c_co1,ext_aperture_B,RO_fl[0],\
                                    GA_fl[0],NSigma_Marsh, S_Marsh, \
				    N_Marsh, Marsh_alg, min_extract_col,\
				    max_extract_col, npools)

    P_ob_R = GLOBALutils.obtain_P(flR,c_ob2,ext_aperture_R,RO_fl[1],\
                                    GA_fl[0],NSigma_Marsh, S_Marsh, \
				    N_Marsh, Marsh_alg, min_extract_col,\
				    max_extract_col, npools)
    P_co_R = GLOBALutils.obtain_P(flR,c_co2,ext_aperture_R,RO_fl[1],\
                                    GA_fl[0],NSigma_Marsh, S_Marsh, \
				    N_Marsh, Marsh_alg, min_extract_col,\
				    max_extract_col, npools)


    print "\t\t\tWill extract",nord_ob1,"blue orders"

    B_flat_ob  = GLOBALutils.optimal_extraction(flB,P_ob_B,c_ob1,ext_aperture_B,RO_fl[0],GA_fl[0],\
                                                S_Marsh,NCosmic_Marsh*10,min_extract_col,\
                                                max_extract_col,npools)
    B_flat_ob = GLOBALutils.invert(B_flat_ob)

    print "\t\t\tWill extract",nord_ob2,"red orders"

    R_flat_ob  = GLOBALutils.optimal_extraction(flR,P_ob_R,c_ob2,ext_aperture_R,RO_fl[1],GA_fl[1],\
                                                S_Marsh,NCosmic_Marsh*10,min_extract_col,\
                                                max_extract_col,npools)
    R_flat_ob = GLOBALutils.invert(R_flat_ob)

    B_flat_ob,R_flat_ob = B_flat_ob[::-1],R_flat_ob[::-1]

    # write P_on and S_flat_ob as fits files
    if (os.access(P_ob_B_fits,os.F_OK)):
        os.remove( P_ob_B_fits )
    if (os.access(B_flat_ob_fits,os.F_OK)):
        os.remove( B_flat_ob_fits )
    if (os.access(P_ob_R_fits,os.F_OK)):
        os.remove( P_ob_R_fits )
    if (os.access(R_flat_ob_fits,os.F_OK)):
        os.remove( R_flat_ob_fits )
    
    hdu = pyfits.PrimaryHDU( P_ob_B )
    hdu.writeto( P_ob_B_fits )
    hdu = pyfits.PrimaryHDU( P_ob_R )
    hdu.writeto( P_ob_R_fits )
    hdu = pyfits.PrimaryHDU( B_flat_ob )
    hdu.writeto( B_flat_ob_fits )
    hdu = pyfits.PrimaryHDU( R_flat_ob )
    hdu.writeto( R_flat_ob_fits )

else:
    print "\t\tExtracted flat object spectra found, loading..."
    P_ob_B       = pyfits.getdata( P_ob_B_fits )
    P_ob_R       = pyfits.getdata( P_ob_R_fits )
    B_flat_ob  = pyfits.getdata( B_flat_ob_fits )
    R_flat_ob  = pyfits.getdata( R_flat_ob_fits )

# Extract flat spectra, comparison
P_co_B_fits    = dirout + 'P_co_B_' + mode + '.fits'
P_co_R_fits    = dirout + 'P_co_R_' + mode + '.fits'
B_flat_co_fits = dirout + 'B_flat_co_' + mode + '.fits'
R_flat_co_fits = dirout + 'R_flat_co_' + mode + '.fits'

P_co_B    = np.zeros( (Flat_co.shape[0],Flat_co.shape[1]) )
P_co_R    = np.zeros( (Flat_co.shape[0],Flat_co.shape[1]) )
B_flat_co = np.zeros((nord_co1, 3, Flat_co.shape[1]) )
R_flat_co = np.zeros((nord_co2, 3, Flat_co.shape[1]) )

if ( os.access(P_co_B_fits,os.F_OK) == False ) or ( os.access(B_flat_co_fits,os.F_OK) == False ) or ( os.access(P_co_R_fits,os.F_OK) == False ) or ( os.access(R_flat_co_fits,os.F_OK) == False ) or (force_flat_extract):
    print "\t\tNo extracted flat comparison spectra found or extraction forced, extracting and saving..."
    
    bacfile = dirout + 'BACR_FLAT_CO_'+mode+'.fits'
    if os.access(bacfile,os.F_OK) == False:
        CentersR = np.zeros((len(c_co2),Flat_co[:,:,1].shape[1]))
        for i in range(len(c_co2)):
            CentersR[i,:]=np.polyval(c_co2[i],np.arange(Flat_co[:,:,1].shape[1]))
        bacR = GLOBALutils.get_scat(Flat_co[:,:,1],CentersR,span=20)
        if (os.access(bacfile,os.F_OK)):
            os.remove( bacfile )
        hdbac = pyfits.PrimaryHDU( bacR )
        hdbac.writeto(bacfile)
    else:
        bacR = pyfits.getdata(bacfile)

    flR = Flat_co[:,:,1] - bacR
    
    bacfile = dirout + 'BACB_FLAT_CO_'+mode+'.fits'
    if os.access(bacfile,os.F_OK) == False:
        CentersB = np.zeros((len(c_co1),Flat_co[:,:,0].shape[1]))
        for i in range(len(c_co1)):
            CentersB[i,:]=np.polyval(c_co1[i],np.arange(Flat_co[:,:,0].shape[1]))
        bacB = GLOBALutils.get_scat(Flat_co[:,:,0],CentersB,span=10)
        if (os.access(bacfile,os.F_OK)):
            os.remove( bacfile )
        hdbac = pyfits.PrimaryHDU( bacB )
        hdbac.writeto(bacfile)
    else:
        bacB = pyfits.getdata(bacfile)

    flB = Flat_co[:,:,0] - bacB

    print "\t\t\tWill extract",nord_co1,"blue orders"
    for i in range(nord_co1):
        P_marsh = GLOBALutils.PCoeff( flB, c_co1[i,:], ext_aperture_B, RO_co[0], GA_co[0], NSigma_Marsh,\
				      S_Marsh, N_Marsh, Marsh_alg, min_extract_col,max_extract_col )
        P_co_B    += P_marsh

    B_flat_co  = GLOBALutils.optimal_extraction(Flat[:,:,0],P_co_B,c_co1,ext_aperture_B,RO_fl[0],GA_fl[0],\
                                                S_Marsh,NCosmic_Marsh*10,min_extract_col,\
                                                max_extract_col,npools)
    B_flat_co = GLOBALutils.invert(B_flat_co)

    print "\t\t\tWill extract",nord_co2,"red orders"
    for i in range(nord_co2):
        P_marsh = GLOBALutils.PCoeff( flR, c_co2[i,:], ext_aperture_R, RO_co[1], GA_co[1], NSigma_Marsh,\
				      S_Marsh, N_Marsh, Marsh_alg, min_extract_col,max_extract_col )
        P_co_R    += P_marsh

    R_flat_co  = GLOBALutils.optimal_extraction(Flat[:,:,1],P_co_R,c_co2,ext_aperture_R,RO_fl[1],GA_fl[1],\
                                                S_Marsh,NCosmic_Marsh*10,min_extract_col,\
                                                max_extract_col,npools)
    R_flat_co = GLOBALutils.invert(R_flat_co)

    B_flat_co,R_flat_co = B_flat_co[::-1],R_flat_co[::-1]

    # write P_on and S_flat_co as fits files
    if (os.access(P_co_B_fits,os.F_OK)):
        os.remove( P_co_B_fits )
    if (os.access(B_flat_co_fits,os.F_OK)):
        os.remove( B_flat_co_fits )
    if (os.access(P_co_R_fits,os.F_OK)):
        os.remove( P_co_R_fits )
    if (os.access(R_flat_co_fits,os.F_OK)):
        os.remove( R_flat_co_fits )
    
    hdu = pyfits.PrimaryHDU( P_co_B )
    hdu.writeto( P_co_B_fits )
    hdu = pyfits.PrimaryHDU( P_co_R )
    hdu.writeto( P_co_R_fits )
    hdu = pyfits.PrimaryHDU( B_flat_co )
    hdu.writeto( B_flat_co_fits )
    hdu = pyfits.PrimaryHDU( R_flat_co )
    hdu.writeto( R_flat_co_fits )

else:
    print "\t\tExtracted flat comparison spectra found, loading..."
    P_co_B    = pyfits.getdata( P_co_B_fits )
    P_co_R    = pyfits.getdata( P_co_R_fits )
    B_flat_co = pyfits.getdata( B_flat_co_fits )
    R_flat_co = pyfits.getdata( R_flat_co_fits )


# Normalize flat field spectra.
B_flat_ob_n,Bnorms = GLOBALutils.FlatNormalize_single( B_flat_ob, mid=2048)
R_flat_ob_n,Rnorms = GLOBALutils.FlatNormalize_single( R_flat_ob, mid=2048)

print '\n\tExtraction of ThAr calibration frames:'
# Extract all ThAr files
for fsim in ThAr_ref:
    print "\t\tWorking on ThAr+Ne file ", fsim, "..."
    hthar = pyfits.open( fsim )
    dtharB = harpsutils.OverscanTrim( hthar[1].data ) - MasterBias[:,:,0]
    dtharR = harpsutils.OverscanTrim( hthar[2].data ) - MasterBias[:,:,1]

    bacfile = dirout + 'BACR_' + fsim.split('/')[-1][:-4]+'fits'
    if os.access(bacfile,os.F_OK) == False:
        CentersR = np.zeros((len(c_all2),dtharR.shape[1]))
        for i in range(len(c_all2)):
            CentersR[i,:]=np.polyval(c_all2[i],np.arange(dtharR.shape[1]))
        bacR = GLOBALutils.get_scat(dtharR,CentersR,span=15)

        if (os.access(bacfile,os.F_OK)):
            os.remove( bacfile )
        hdbac = pyfits.PrimaryHDU( bacR )
        hdbac.writeto(bacfile)
    else:
        bacR = pyfits.getdata(bacfile)

    bacfile = dirout + 'BACB_' + fsim.split('/')[-1][:-4]+'fits'
    if os.access(bacfile,os.F_OK) == False:
        CentersB = np.zeros((len(c_all1),dtharB.shape[1]))
        for i in range(len(c_all1)):
            CentersB[i,:]=np.polyval(c_all1[i],np.arange(dtharB.shape[1]))
        bacB = GLOBALutils.get_scat(dtharB,CentersB,span=7)
        if (os.access(bacfile,os.F_OK)):
            os.remove( bacfile )
        hdbac = pyfits.PrimaryHDU( bacB )
        hdbac.writeto(bacfile)
    else:
        bacB = pyfits.getdata(bacfile)

    dtharR = dtharR - bacR
    dtharB = dtharB - bacB

    thar_fits_ob_R = dirout + fsim.split('/')[-1][:-4]+'spec.ob.R.fits.S'
    thar_fits_co_R = dirout + fsim.split('/')[-1][:-4]+'spec.co.R.fits.S'
    thar_fits_ob_B = dirout + fsim.split('/')[-1][:-4]+'spec.ob.B.fits.S'
    thar_fits_co_B = dirout + fsim.split('/')[-1][:-4]+'spec.co.B.fits.S'
    if ( os.access(thar_fits_ob_B,os.F_OK) == False ) or ( os.access(thar_fits_co_B,os.F_OK) == False ) \
    or ( os.access(thar_fits_ob_R,os.F_OK) == False ) or ( os.access(thar_fits_co_R,os.F_OK) == False ) \
    or ( force_thar_extract ):
        print "\t\tNo previous extraction or extraction forced for ThAr file", fsim, "extracting..."

	tR,tG = hthar[1].header['HIERARCH ESO DET OUT1 RON'],hthar[1].header['HIERARCH ESO DET OUT1 GAIN']

	thar_S_ob_B  = GLOBALutils.optimal_extraction(dtharB,P_ob_B,c_ob1,ext_aperture_B,tR,tG,\
                                                    S_Marsh,0*NCosmic_Marsh,min_extract_col,max_extract_col,npools)
        for i in range(nord_ob1):
	    thar_S_ob_B[i,1,:] = thar_S_ob_B[i,1,:][::-1]
	    thar_S_ob_B[i,2,:] = thar_S_ob_B[i,2,:][::-1]

	thar_S_co_B  = GLOBALutils.optimal_extraction(dtharB,P_co_B,c_co1,ext_aperture_B,tR,tG,\
                                                    S_Marsh,0*NCosmic_Marsh,min_extract_col,max_extract_col,npools)
        for i in range(nord_co1):
	    thar_S_co_B[i,1,:] = thar_S_co_B[i,1,:][::-1]
	    thar_S_co_B[i,2,:] = thar_S_co_B[i,2,:][::-1]

	thar_S_ob_B,thar_S_co_B = thar_S_ob_B[::-1],thar_S_co_B[::-1]

	tR,tG = hthar[2].header['HIERARCH ESO DET OUT1 RON'],hthar[2].header['HIERARCH ESO DET OUT1 GAIN']

	thar_S_ob_R  = GLOBALutils.optimal_extraction(dtharR,P_ob_R,c_ob2,ext_aperture_R,tR,tG,\
                                                    S_Marsh,0*NCosmic_Marsh,min_extract_col,max_extract_col,npools)
        for i in range(nord_ob2):
	    thar_S_ob_R[i,1,:] = thar_S_ob_R[i,1,:][::-1]
	    thar_S_ob_R[i,2,:] = thar_S_ob_R[i,2,:][::-1]

	thar_S_co_R  = GLOBALutils.optimal_extraction(dtharR,P_co_R,c_co2,ext_aperture_R,tR,tG,\
                                                    S_Marsh,0*NCosmic_Marsh,min_extract_col,max_extract_col,npools)
        for i in range(nord_co2):
	    thar_S_co_R[i,1,:] = thar_S_co_R[i,1,:][::-1]
	    thar_S_co_R[i,2,:] = thar_S_co_R[i,2,:][::-1]

	thar_S_ob_R,thar_S_co_R = thar_S_ob_R[::-1],thar_S_co_R[::-1]

        # save as fits file
        if (os.access(thar_fits_ob_R,os.F_OK)):
            os.remove( thar_fits_ob_R )
        if (os.access(thar_fits_ob_B,os.F_OK)):
            os.remove( thar_fits_ob_B )
        if (os.access(thar_fits_co_R,os.F_OK)):
            os.remove( thar_fits_co_R )
        if (os.access(thar_fits_co_B,os.F_OK)):
            os.remove( thar_fits_co_B )
            
        hdu = pyfits.PrimaryHDU( thar_S_ob_B )
        hdu.writeto( thar_fits_ob_B )
        hdu = pyfits.PrimaryHDU( thar_S_ob_R )
        hdu.writeto( thar_fits_ob_R )
        hdu = pyfits.PrimaryHDU( thar_S_co_B )
        hdu.writeto( thar_fits_co_B )
        hdu = pyfits.PrimaryHDU( thar_S_co_R )
        hdu.writeto( thar_fits_co_R )
    else:
        print "\t\tThAr file", fsim, "all ready extracted, loading..."

# create wavelength calibration files
print "\n\tWavelength solution of ThAr calibration spectra:"
sorted_ThAr_dates = np.argsort( ThAr_ref_dates )

for i in range(len(ThAr_ref_dates)):
    index = sorted_ThAr_dates[i]
    wavsol_pkl  = dirout + ThAr_ref[index].split('/')[-1][:-4]+'wavsolpars.pkl'
    wavsol_fits = dirout + ThAr_ref[index].split('/')[-1][:-4]+'spec.fits'
    #force_thar_wavcal = True
    if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
        print "\t\tWorking on initial ThAr file", ThAr_ref[index] 
        hthar = pyfits.open( ThAr_ref[index] )
        mjd, mjd0 = harpsutils.mjd_fromheader( hthar )
        thar_fits_ob_B = dirout + ThAr_ref[index].split('/')[-1][:-4]+'spec.ob.B.fits.S'
        thar_fits_co_B = dirout + ThAr_ref[index].split('/')[-1][:-4]+'spec.co.B.fits.S'
        thar_fits_ob_R = dirout + ThAr_ref[index].split('/')[-1][:-4]+'spec.ob.R.fits.S'
        thar_fits_co_R = dirout + ThAr_ref[index].split('/')[-1][:-4]+'spec.co.R.fits.S'

        thar_S_ob_B = pyfits.getdata( thar_fits_ob_B )
        thar_S_co_B = pyfits.getdata( thar_fits_co_B )
        thar_S_ob_R = pyfits.getdata( thar_fits_ob_R )
        thar_S_co_R = pyfits.getdata( thar_fits_co_R )

        lines_thar_ob_B  = thar_S_ob_B[:,1,:]
        iv_thar_ob_B     = thar_S_ob_B[:,2,:]
        lines_thar_co_B  = thar_S_co_B[:,1,:]
        iv_thar_co_B     = thar_S_co_B[:,2,:]
        lines_thar_ob_R  = thar_S_ob_R[:,1,:]
        iv_thar_ob_R     = thar_S_ob_R[:,2,:]
        lines_thar_co_R  = thar_S_co_R[:,1,:]
        iv_thar_co_R     = thar_S_co_R[:,2,:]
	
	c_p2w_ob_B = np.zeros((nord_ob1,porder+1))
	c_p2w_ob_R = np.zeros((nord_ob2,porder+1))
	c_p2w_co_B = np.zeros((nord_co1,porder+1))
	c_p2w_co_R = np.zeros((nord_co2,porder+1))
	spec_thar_ob = np.zeros((2,nord_ob1+nord_ob2,thar_S_ob_B.shape[2]))

	All_Pixel_Centers_R = np.array([])
	All_Wavelengths_R   = np.array([])
	All_Orders_R        = np.array([])
	All_Centroids_R     = np.array([])
	All_Sigmas_R        = np.array([])
	All_Intensities_R   = np.array([])

	counter = 0
	temp_pix = np.array([])
	temp_res = np.array([])
	meds,ords = [],[]
	for order in range(nord_ob2):
            order_s = str(order)
            if (order < 10):
                order_s = '0' + str(order)
            thar_order_orig = lines_thar_ob_R[order,:]
            IV              = iv_thar_ob_R[order,:]
            wei             = np.sqrt( IV )
            #bkg             = CoralieUtils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
            thar_order      = thar_order_orig #- bkg

	    coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms,\
		 residuals, centroids,sigmas, intensities = \
		 GLOBALutils.Initial_Wav_Calibration( order_dir+'R_order_'+order_s+\
		 '.iwdat', thar_order, order, wei, rmsmax=MRMS_initial, minlines=minlines_initial,\
		 FixEnds=False,Dump_Argon=dumpargon, Dump_AllLines=True, Cheby=use_cheby, porder=porder)

            if (order == int(np.around(0.5*nord_ob2))): 
                if (use_cheby):
                    Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, int(np.around(0.5*len(thar_order))), len(thar_order) )
                else:
                    Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

	    All_Pixel_Centers_R = np.append( All_Pixel_Centers_R, pixel_centers )
	    All_Wavelengths_R   = np.append( All_Wavelengths_R, wavelengths )
	    All_Orders_R        = np.append( All_Orders_R, np.zeros( len(pixel_centers) ) + order )
	    All_Centroids_R     = np.append( All_Centroids_R, centroids)
	    All_Sigmas_R        = np.append( All_Sigmas_R, sigmas)
	    All_Intensities_R   = np.append( All_Intensities_R, intensities )

    	    meds.append(GLOBALutils.Cheby_eval(coeffs_pix2wav,0.5*len(thar_order),len(thar_order)))
	    ords.append(order)
	    spec_thar_ob[0,counter] =  GLOBALutils.Cheby_eval(coeffs_pix2wav,np.arange(len(thar_order)),len(thar_order))
	    spec_thar_ob[1,counter] =  thar_order
	    counter += 1
	    pp1 = GLOBALutils.Cheby_eval(coeffs_pix2wav,pixel_centers + sigmas,len(thar_order))
	    pm1 = GLOBALutils.Cheby_eval(coeffs_pix2wav,pixel_centers - sigmas,len(thar_order))
	    wavsigmas = .5*(pp1 - pm1)
	    reses = wavelengths/(wavsigmas*2.355)
	    jji = 0
	    vecp,vecm = [],[]
	    while jji < len(reses):
		if jji + 5 < len(reses):
			argm = np.argmax(reses[jji:jji+5])
			vecp.append(wavelengths[jji+argm])
			vecm.append(reses[jji+argm])
		jji += 5
	    vecp,vecm = np.array(vecp),np.array(vecm)
	    coef_res = np.polyfit(vecp,vecm,2)
	    #plot(wavelengths,np.polyval(coef_res,wavelengths))
	    #plot(vecp,vecm,'o')
	    #print order, rms_ms/np.sqrt(float(len(wavelengths))), rms_ms, len(residuals)
	    
	    c_p2w_ob_R[order] = coeffs_pix2wav
	    isz  = pixel_centers - sigmas
	    der  = pixel_centers + sigmas
	    isz  = GLOBALutils.Cheby_eval( coeffs_pix2wav, isz,len(thar_order))
	    der  = GLOBALutils.Cheby_eval( coeffs_pix2wav, der,len(thar_order))
	    sig  = 0.5*(der-isz)
	    fwhm = 2.35 * sig
	    resol = wavelengths / fwhm
	    temp_pix = np.hstack((temp_pix,pixel_centers))
	    temp_res = np.hstack((temp_res,resol))
	    #plot(pixel_centers,resol,'o')

        p0 = np.zeros( npar_wsol_R )
        p0[0] =  (int(np.around(0.5*nord_ob2))+or0_R) * Global_ZP 
	#GLOBALutils.get_zero_order_number(ords,meds)
	p1_R, G_pix_R, G_ord_R, G_wav_R, II_R, rms_ms_R, G_res_R = \
		GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_R, All_Wavelengths_R, All_Orders_R,\
		                                             np.ones(All_Intensities_R.shape), p0, Cheby=use_cheby,\
		                                             maxrms=MRMS, Inv=Inverse_m, minlines=minlines_glob,\
                                                             order0=or0_R,ntotal=nord_ob2,npix=len(thar_order),nx=ncoef_x_R,nm=ncoef_m_R)


	I = np.argsort(temp_pix)
	temp_pix,temp_res = temp_pix[I],temp_res[I]
	jj = 0
	xx,yy = [],[]
	while jj<len(thar_order):
		I = np.where((temp_pix>=jj) & (temp_pix<jj+50))[0]
		if len(I)>0:
			xt,yt = temp_pix[I],temp_res[I]
			I = np.argmax(yt)
			xx.append(xt[I])
			yy.append(yt[I])
		jj+=50
	xx = np.array(xx)
	yy = np.array(yy)
	coefs_m = np.polyfit(ords,meds,5)
	#plot(ords, meds - np.polyval(coefs_m,ords),'ro')
	#show()
	#plot(xx,yy)
	coef = np.polyfit(xx,yy,3)
	#print 'Resolution coefs',coef
	#plot(np.arange(4096),np.polyval(coef,np.arange(4096)),linewidth=2.0)
	#show()

	All_Pixel_Centers_B = np.array([])
	All_Wavelengths_B   = np.array([])
	All_Orders_B        = np.array([])
	All_Centroids_B     = np.array([])
	All_Sigmas_B        = np.array([])
	All_Intensities_B   = np.array([])

	meds,ords = [],[]
	for order in range(nord_ob1):
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)
            thar_order_orig = lines_thar_ob_B[order,:]
            IV              = iv_thar_ob_B[order,:]
            wei             = np.sqrt( IV )
            #bkg             = CoralieUtils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
            thar_order      = thar_order_orig #- bkg
	    coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, \
			centroids_B, sigmas, intensities = GLOBALutils.Initial_Wav_Calibration( order_dir\
			+'B_order_'+order_s+final_wav, thar_order, order, wei, rmsmax=MRMS_initial,\
			minlines=minlines_initial, FixEnds=False,Dump_Argon=dumpargon, Dump_AllLines=True,\
			Cheby=use_cheby, porder=porder)

            if (order == int(np.around(0.5*nord_ob1))): 
                if (use_cheby):
                    Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, int(np.around(0.5*len(thar_order))), len(thar_order) )
                else:
                    Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

	    All_Pixel_Centers_B = np.append( All_Pixel_Centers_B, pixel_centers )
	    All_Wavelengths_B   = np.append( All_Wavelengths_B, wavelengths )
	    All_Orders_B        = np.append( All_Orders_B, np.zeros( len(pixel_centers) ) + order )
	    All_Centroids_B     = np.append( All_Centroids_B, centroids)
	    All_Sigmas_B        = np.append( All_Sigmas_B, sigmas)
	    All_Intensities_B   = np.append( All_Intensities_B, intensities )

	    meds.append(GLOBALutils.Cheby_eval(coeffs_pix2wav,0.5*len(thar_order),len(thar_order)))
	    ords.append(order)
	    spec_thar_ob[0,counter] =  GLOBALutils.Cheby_eval(coeffs_pix2wav,np.arange(len(thar_order)),len(thar_order))
	    spec_thar_ob[1,counter] =  thar_order
	    counter += 1
	    pp1 = GLOBALutils.Cheby_eval(coeffs_pix2wav,pixel_centers + sigmas,len(thar_order))
	    pm1 = GLOBALutils.Cheby_eval(coeffs_pix2wav,pixel_centers - sigmas,len(thar_order))
	    wavsigmas = .5*(pp1 - pm1)
	    reses = wavelengths/(wavsigmas*2.355)
	    jji = 0
	    vecp,vecm = [],[]
	    while jji < len(reses):
		if jji + 5 < len(reses):
			argm = np.argmax(reses[jji:jji+5])
			vecp.append(wavelengths[jji+argm])
			vecm.append(reses[jji+argm])
		jji += 5
	    vecp,vecm = np.array(vecp),np.array(vecm)
	    coef_res = np.polyfit(vecp,vecm,2)
	    #plot(wavelengths,np.polyval(coef_res,wavelengths))
	    #plot(vecp,vecm,'o')
	    c_p2w_ob_B[order] = coeffs_pix2wav
	    #print order, rms_ms/np.sqrt(float(len(wavelengths))), rms_ms, len(residuals)
	
        p0 = np.zeros( npar_wsol_B )
        p0[0] =  (int(np.around(0.5*nord_ob1))+or0_B) * Global_ZP 
	#GLOBALutils.get_zero_order_number(ords,meds)
	p1_B, G_pix_B, G_ord_B, G_wav_B, II_B, rms_ms_B, G_res_B = \
		GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_B, All_Wavelengths_B, All_Orders_B,\
		                                             np.ones(All_Intensities_B.shape), p0, Cheby=use_cheby,\
		                                             maxrms=MRMS, Inv=Inverse_m, minlines=minlines_glob,\
                                                             order0=or0_B,ntotal=nord_ob1,npix=len(thar_order),nx=ncoef_x_B,nm=ncoef_m_B)


	nhdu = pyfits.PrimaryHDU(spec_thar_ob)
	if os.access(wavsol_fits,os.F_OK):
		os.system('rm '+wavsol_fits)
	nhdu.writeto(wavsol_fits)

	#plot(ords,meds,'ro')
	#coefs_m = np.polyfit(ords,meds,6)
	#plot(ords, meds - np.polyval(coefs_m,ords),'ro')
	#show()
	All_Pixel_Centers_co_R = np.array([])
	All_Wavelengths_co_R   = np.array([])
	All_Orders_co_R        = np.array([])
	All_Centroids_co_R     = np.array([])
	All_Sigmas_co_R        = np.array([])
	All_Intensities_co_R   = np.array([])
	meds,ords = [],[]
        for order in range(nord_co2):
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)
            
            thar_order_orig = lines_thar_co_R[order,:]
            IV              = iv_thar_co_R[order,:]
            wei             = np.sqrt( IV )
            #bkg             = CoralieUtils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
            thar_order      = thar_order_orig #- bkg
            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
		        = GLOBALutils.Initial_Wav_Calibration( order_dir+'R_order_'+order_s+final_wav, thar_order, order, wei, rmsmax=MRMS_initial, minlines=minlines_initial,\
							      FixEnds=True,Dump_Argon=dumpargon, Dump_AllLines=True, Cheby=use_cheby, porder=porder)

	    c_p2w_co_R[order] = coeffs_pix2wav
	    meds.append(GLOBALutils.Cheby_eval(coeffs_pix2wav,0.5*len(thar_order),len(thar_order)))
	    ords.append(order)
            if (order == int(np.around(0.5*nord_co2))): 
                if (use_cheby):
                    Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, int(np.around(0.5*len(thar_order))), len(thar_order) )
                else:
                    Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

	    All_Pixel_Centers_co_R = np.append( All_Pixel_Centers_co_R, pixel_centers )
	    All_Wavelengths_co_R   = np.append( All_Wavelengths_co_R, wavelengths )
	    All_Orders_co_R        = np.append( All_Orders_co_R, np.zeros( len(pixel_centers) ) + order )
	    All_Centroids_co_R     = np.append( All_Centroids_co_R, centroids)
	    All_Sigmas_co_R        = np.append( All_Sigmas_co_R, sigmas)
	    All_Intensities_co_R   = np.append( All_Intensities_co_R, intensities )
	    
	    #print order, rms_ms/np.sqrt(float(len(wavelengths))), rms_ms, len(residuals)
        p0 = np.zeros( npar_wsol_R )
        p0[0] =  (int(np.around(0.5*nord_co2))+or0_R) * Global_ZP 
	#GLOBALutils.get_zero_order_number(ords,meds)
	p1_co_R, G_pix_co_R, G_ord_co_R, G_wav_co_R, II_co_R, rms_ms_co_R, G_res_co_R = \
		GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co_R, All_Wavelengths_co_R, All_Orders_co_R,\
		                                             np.ones(All_Intensities_co_R.shape), p0, Cheby=use_cheby,\
		                                             maxrms=MRMS, Inv=Inverse_m, minlines=minlines_glob,\
                                                             order0=or0_R,ntotal=nord_co2,npix=len(thar_order),nx=ncoef_x_R,nm=ncoef_m_R)
	All_Pixel_Centers_co_B = np.array([])
	All_Wavelengths_co_B   = np.array([])
	All_Orders_co_B        = np.array([])
	All_Centroids_co_B     = np.array([])
	All_Sigmas_co_B        = np.array([])
	All_Intensities_co_B   = np.array([])
	meds,ords = [],[]
        for order in range(nord_co1):
	    order = order + 1
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)
            
            thar_order_orig = lines_thar_co_B[order-1,:]
            IV              = iv_thar_co_B[order-1,:]
            wei             = np.sqrt( IV )
            #bkg             = CoralieUtils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
            thar_order      = thar_order_orig #- bkg
            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals,\
		centroids, sigmas, intensities = GLOBALutils.Initial_Wav_Calibration( order_dir+\
		'B_order_'+order_s+final_wav, thar_order, order, wei, rmsmax=MRMS_initial, \
		minlines=minlines_initial, FixEnds=True,Dump_Argon=dumpargon, Dump_AllLines=True,\
		Cheby=use_cheby, porder=porder)
	    c_p2w_co_B[order-1] = coeffs_pix2wav

	    meds.append(GLOBALutils.Cheby_eval(coeffs_pix2wav,0.5*len(thar_order),len(thar_order)))
	    ords.append(order)
            if (order == int(np.around(0.5*nord_co1))): 
                if (use_cheby):
                    Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, int(np.around(0.5*len(thar_order))), len(thar_order) )
                else:
                    Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

	    All_Pixel_Centers_co_B = np.append( All_Pixel_Centers_co_B, pixel_centers )
	    All_Wavelengths_co_B   = np.append( All_Wavelengths_co_B, wavelengths )
	    All_Orders_co_B        = np.append( All_Orders_co_B, np.zeros( len(pixel_centers) ) + order )
	    All_Centroids_co_B     = np.append( All_Centroids_co_B, centroids)
	    All_Sigmas_co_B        = np.append( All_Sigmas_co_B, sigmas)
	    All_Intensities_co_B   = np.append( All_Intensities_co_B, intensities )

	    #print order, rms_ms/np.sqrt(float(len(wavelengths))), rms_ms, len(residuals)

        p0 = np.zeros( npar_wsol_B )
        p0[0] =  (int(np.around(0.5*nord_co1))+or0_B) * Global_ZP 
	#GLOBALutils.get_zero_order_number(ords,meds)
	p1_co_B, G_pix_co_B, G_ord_co_B, G_wav_co_B, II_co_B, rms_ms_co_B, G_res_co_B = \
		GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers_co_B, All_Wavelengths_co_B, All_Orders_co_B,\
		                                             np.ones(All_Intensities_co_B.shape), p0, Cheby=use_cheby,\
		                                             maxrms=MRMS, Inv=Inverse_m, minlines=minlines_glob,\
                                                             order0=or0_B,ntotal=nord_co1,npix=len(thar_order),nx=ncoef_x_B,nm=ncoef_m_B)

	#pdict = {'p2w_ob_R':c_p2w_ob_R, 'p2w_ob_B':c_p2w_ob_B, 'p2w_co_R':c_p2w_co_R, 'p2w_co_B':c_p2w_co_B}

        pdict = {'p1_R':p1_R,'p1_B':p1_B,'p1_co_R':p1_co_R,'p1_co_B':p1_co_B,\
	         'G_pix_R':G_pix_R, 'G_pix_B':G_pix_B,'G_pix_co_R':G_pix_co_R, 'G_pix_co_B':G_pix_co_B,\
		 'G_ord_R':G_ord_R, 'G_ord_B':G_ord_B,'G_ord_co_R':G_ord_co_R, 'G_ord_co_B':G_ord_co_B,\
		 'G_wav_R':G_wav_R, 'G_wav_B':G_wav_B,'G_wav_co_R':G_wav_co_R, 'G_wav_co_B':G_wav_co_B,\
		 'II_R':II_R,'II_B':II_B,'II_co_R':II_co_R,'II_co_B':II_co_B,\
		 'rms_ms_R':rms_ms_R,'rms_ms_B':rms_ms_B,'rms_ms_co_R':rms_ms_co_R,'rms_ms_co_B':rms_ms_co_B,\
		 'G_res_R':G_res_R, 'G_res_B':G_res_B,'G_res_co_R':G_res_co_R, 'G_res_co_B':G_res_co_B,\
		 'All_Centroids_R':All_Centroids_R,'All_Centroids_B':All_Centroids_B,\
		 'All_Centroids_co_R':All_Centroids_co_R,'All_Centroids_co_B':All_Centroids_co_B,\
		 'All_Orders_R':All_Orders_R,'All_Orders_B':All_Orders_B,\
		 'All_Orders_co_R':All_Orders_co_R,'All_Orders_co_B':All_Orders_co_B,\
		 'All_Sigmas_R':All_Sigmas_R,'All_Sigmas_B':All_Sigmas_B,\
		 'All_Sigmas_co_R':All_Sigmas_co_R,'All_Sigmas_co_B':All_Sigmas_co_B,\
	         'mjd':mjd,'npix':len(thar_order)}

        pickle.dump( pdict, open( wavsol_pkl, 'w' ) )


    else:
        print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl

xc_fout_f = dirout+'ThAr_XCor+DeltaL.dat'
if ( (os.access(xc_fout_f,os.F_OK) == False) or (force_tharxc)):  
    xc_fout = open(xc_fout_f,'w')
    for i in range(len(sorted_ThAr_dates)):
        index = sorted_ThAr_dates[i]

        fsim = ThAr_ref[index]
        hthar = pyfits.open( fsim )
        mjd, mjd0 = harpsutils.mjd_fromheader( hthar )
	pdict1 = pickle.load( open(dirout + fsim.split('/')[-1][:-4]+'wavsolpars.pkl','r' ) )

        for ii in range(len(sorted_ThAr_dates)):
	    index2 = sorted_ThAr_dates[ii]

            fsim2 = ThAr_ref[index2]
            hthar2 = pyfits.open( fsim2 )
            mjd2, mjd02 = harpsutils.mjd_fromheader( hthar2 )
	    pdict2 = pickle.load( open(dirout + fsim2.split('/')[-1][:-4]+'wavsolpars.pkl','r' ) )
                
	    p_shift_R, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(pdict2['G_pix_R'], pdict2['G_wav_R'], pdict2['G_ord_R'],\
		                                                   np.ones(pdict2['G_wav_R'].shape), pdict1['p1_R'],\
		                                                   Cheby=True,Inv=True,maxrms=MRMS,minlines=minlines_glob,\
								   order0=or0_R,ntotal=nord_ob2,npix=pdict2['npix'],nx=ncoef_x_R,nm=ncoef_m_R)

	    p_shift_B, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(pdict2['G_pix_B'], pdict2['G_wav_B'], pdict2['G_ord_B'],\
		                                                   np.ones(pdict2['G_wav_B'].shape), pdict1['p1_B'],\
		                                                   Cheby=True,Inv=True,maxrms=MRMS,minlines=minlines_glob,\
								   order0=or0_B,ntotal=nord_ob1,npix=pdict2['npix'],nx=ncoef_x_B,nm=ncoef_m_B)

	    p_shift_co_R, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(pdict2['G_pix_co_R'], pdict2['G_wav_co_R'], pdict2['G_ord_co_R'],\
		                                                   np.ones(pdict2['G_wav_co_R'].shape), pdict1['p1_co_R'],\
		                                                   Cheby=True,Inv=True,maxrms=MRMS,minlines=minlines_glob,\
								   order0=or0_R,ntotal=nord_co2,npix=pdict2['npix'],nx=ncoef_x_R,nm=ncoef_m_R)

	    p_shift_co_B, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(pdict2['G_pix_co_B'], pdict2['G_wav_co_B'], pdict2['G_ord_co_B'],\
		                                                   np.ones(pdict2['G_wav_co_B'].shape), pdict1['p1_co_B'],\
		                                                   Cheby=True,Inv=True,maxrms=MRMS,minlines=minlines_glob,\
								   order0=or0_B,ntotal=nord_co1,npix=pdict2['npix'],nx=ncoef_x_B,nm=ncoef_m_B)


	    # write out products
	    line_out = "%20.12f %20.12f %12.6f %12.6f %12.6f %12.6f %s %s\n" % (mjd, mjd2,(1e-6*p_shift_R)*299792458.0,(1e-6*p_shift_B)*299792458.0,(1e-6*p_shift_co_R)*299792458.0,(1e-6*p_shift_co_B)*299792458.0,fsim.split('/')[-1][:-4],fsim2.split('/')[-1][:-4])
	    xc_fout.write(line_out)
	    #xc_fout.flush()
    xc_fout.close()


### start of science frame reductions ###
new_list         = []
new_list_obnames = []
new_list_texp    = []
new_list_cotypes = []
for i in range(len(sim_sci)):
    fsim    = sim_sci[i]
    obname  = obnames[i]
    texp    = exptimes[i]
    co_type = co_types[i] 
    if (object2do == 'all'):
        new_list.append(fsim)
        new_list_obnames.append( obname )
        new_list_texp.append( texp )
	new_list_cotypes.append( co_type )
    if object2do == 'new':
	h = pyfits.open(fsim)
	fout = 'proc/'+ obname + '_' + h[0].header['DATE-OBS'] + 'sp.fits'
        if not os.access(dirout + fout,os.F_OK):
            new_list.append(fsim)
            new_list_obnames.append( obname )
            new_list_texp.append( texp )
	    new_list_cotypes.append( co_type )
    else:
        if (obname == object2do):
            new_list.append(fsim)
            new_list_obnames.append( obname )
            new_list_texp.append( texp )
	    new_list_cotypes.append( co_type )

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
use_moon = np.array(use_moon)

print '\n\tThe following targets will be processed:\n'
for nlisti in range(len(new_list)):
	print '\t'+new_list_obnames[nlisti]

for nlisti in range(len(new_list)):
    fsim   = new_list[ nlisti ]
    obname = new_list_obnames[nlisti]
    TEXP   = new_list_texp[nlisti]
    cotype = new_list_cotypes[nlisti]

    know_moon = False
    here_moon = False
    if fsim.split('/')[-1] in spec_moon:
        I = np.where(fsim.split('/')[-1] == spec_moon)[0]
        know_moon = True
        here_moon = use_moon[I]

    h = pyfits.open(fsim)
    print '\n'
    print "\t--> Working on image: ", fsim
    print "\t\tObject name:",obname
    mjd,mjd0 = harpsutils.mjd_fromheader(h)

    # Open file, trim, overscan subtract and MasterBias subtract
    dataB = h[1].data
    dataB = harpsutils.OverscanTrim(dataB)
    dataB -= MasterBias[:,:,0]
    dataR = h[2].data
    dataR = harpsutils.OverscanTrim(dataR)
    dataR -= MasterBias[:,:,1]
    
    if cotype != 'WAVE':
	spanR,spanB = 20,10
    else:
	spanR,spanB = 10,7
	
    bacfile = dirout + 'BACR_' + fsim.split('/')[-1][:-4]+'fits'
    if os.access(bacfile,os.F_OK) == False:
        CentersR = np.zeros((len(c_ob2),dataR.shape[1]))
        for i in range(len(c_ob2)):
            CentersR[i,:]=np.polyval(c_ob2[i],np.arange(dataR.shape[1]))
        bacR = GLOBALutils.get_scat(dataR,CentersR,span=20)
        if (os.access(bacfile,os.F_OK)):
            os.remove( bacfile )
        hdbac = pyfits.PrimaryHDU( bacR )
        hdbac.writeto(bacfile)
    else:
        bacR = pyfits.getdata(bacfile)

    bacfile = dirout + 'BACB_' + fsim.split('/')[-1][:-4]+'fits'
    if os.access(bacfile,os.F_OK) == False:
        CentersB = np.zeros((len(c_ob1),dataB.shape[1]))
        for i in range(len(c_ob1)):
            CentersB[i,:]=np.polyval(c_ob1[i],np.arange(dataB.shape[1]))
        bacB = GLOBALutils.get_scat(dataB,CentersB,span=10)
        if (os.access(bacfile,os.F_OK)):
            os.remove( bacfile )
        hdbac = pyfits.PrimaryHDU( bacB )
        hdbac.writeto(bacfile)
    else:
        bacB = pyfits.getdata(bacfile)

    dataR = dataR - bacR
    dataB = dataB - bacB

    ron1,gain1 = h[1].header['HIERARCH ESO DET OUT1 RON'],h[1].header['HIERARCH ESO DET OUT1 GAIN']
    ron2,gain2 = h[2].header['HIERARCH ESO DET OUT1 RON'],h[2].header['HIERARCH ESO DET OUT1 GAIN']
    halfcounts = h[0].header['HIERARCH ESO INS DET1 TMMEAN']

    # Find lambda_bary/lambda_topo using baryc
    altitude    = h[0].header['HIERARCH ESO TEL GEOELEV']
    latitude    = h[0].header['HIERARCH ESO TEL GEOLAT']
    longitude   = h[0].header['HIERARCH ESO TEL GEOLON']
    ra          = h[0].header['RA']
    dec         = h[0].header['DEC']
    epoch       = h[0].header['HIERARCH ESO TEL TARG EQUINOX']


    ra2,dec2 = GLOBALutils.getcoords(obname,mjd,filen=reffile)
    if ra2 !=0 and dec2 != 0:
	ra = ra2
	dec = dec2
    else:
	print 'Using the coordinates found in the image header.'

    iers          = GLOBALutils.JPLiers( baryc_dir, mjd-999.0, mjd+999.0 )
    obsradius, R0 = GLOBALutils.JPLR0( latitude, altitude)
    obpos         = GLOBALutils.obspos( longitude, obsradius, R0 )
    jplephem.set_ephemeris_dir( baryc_dir , ephemeris )
    jplephem.set_observer_coordinates( obpos[0], obpos[1], obpos[2] )

    res         = jplephem.doppler_fraction(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
    lbary_ltopo = 1.0 + res['frac'][0]
    bcvel_baryc = ( lbary_ltopo - 1.0 ) * 2.99792458E5

    print "Barycentric velocity:", bcvel_baryc, mjd

    res  = jplephem.pulse_delay(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
    mbjd = mjd + res['delay'][0] / (3600.0 * 24.0)

    # Moon Phase Calculations
    gobs = ephem.Observer()  
    gobs.name='Eso3.6'  
    gobs.lat  = rad(latitude)  # lat/long in decimal degrees  
    gobs.long = rad(longitude)
    DDATE     = h[0].header['DATE-OBS'][:10]
    HHOUR     = h[0].header['DATE-OBS'][11:]
    gobs.date = str(DDATE[:4]) + '-' +  str(DDATE[5:7]) + '-' + str(DDATE[8:]) + ' ' +  HHOUR[:2] + ':' + HHOUR[3:5] +':' + str(float(HHOUR[6:]) + halfcounts * TEXP )
    mephem    = ephem.Moon()
    mephem.compute(gobs)

    mephem = ephem.Moon()
    mephem.compute(gobs)
    Mcoo = jplephem.object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Mp = jplephem.barycentric_object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Sp = jplephem.barycentric_object_track("Sun", int(mjd), float(mjd%1), 1, 0.0)
    res  = jplephem.object_doppler("Moon", int(mjd), mjd%1, 1, 0.0)
    lunation,moon_state,moonsep,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,ra,dec)
    refvel = bcvel_baryc + moonvel
    print '\t\tRadial Velocity of sacttered moonlight:',refvel

    sorted_indices = np.argsort( np.abs( np.array(ThAr_ref_dates) - mjd ) )

    # optimally and simply extract spectra
    sci_fits_ob_B = dirout + fsim.split('/')[-1][:-4]+'spec.ob.B.fits.S'
    sci_fits_co_B = dirout + fsim.split('/')[-1][:-4]+'spec.co.B.fits.S'
    sci_fits_ob_R = dirout + fsim.split('/')[-1][:-4]+'spec.ob.R.fits.S'
    sci_fits_co_R = dirout + fsim.split('/')[-1][:-4]+'spec.co.R.fits.S'
    #force_sci_extract = True
    if ( os.access(sci_fits_ob_B,os.F_OK) == False ) or ( os.access(sci_fits_co_B,os.F_OK) == False ) \
    or ( os.access(sci_fits_ob_R,os.F_OK) == False ) or ( os.access(sci_fits_co_R,os.F_OK) == False ) \
    or (force_sci_extract):

        print "No previous extraction or extraction forced for science file", fsim, "extracting..."
        sci_S_ob_B  = np.zeros( (nord_ob1,3,dataB.shape[1]) )
        sci_S_co_B  = np.zeros( (nord_co1,3,dataB.shape[1]) )
        sci_S_ob_R  = np.zeros( (nord_ob2,3,dataR.shape[1]) )
        sci_S_co_R  = np.zeros( (nord_co2,3,dataR.shape[1]) )
        tpars1,tpars2 = [],[]

        sci_S_ob_B  = GLOBALutils.optimal_extraction(dataB,P_ob_B,c_ob1,ext_aperture_B,\
                                                   ron1,gain1,S_Marsh,NCosmic_Marsh,\
                                                   min_extract_col,max_extract_col,npools)

	sci_S_co_B  = GLOBALutils.optimal_extraction(dataB,P_co_B,c_co1,ext_aperture_B,\
                                                   ron1,gain1,S_Marsh,100.*NCosmic_Marsh,\
                                                   min_extract_col,max_extract_col,npools)

        sci_S_ob_R  = GLOBALutils.optimal_extraction(dataR,P_ob_R,c_ob2,ext_aperture_R,\
                                                   ron2,gain2,S_Marsh,NCosmic_Marsh,\
                                                   min_extract_col,max_extract_col,npools)

	sci_S_co_R  = GLOBALutils.optimal_extraction(dataR,P_co_R,c_co2,ext_aperture_R,\
                                                   ron2,gain2,S_Marsh,100.*NCosmic_Marsh,\
                                                   min_extract_col,max_extract_col,npools)

	sci_S_ob_B  = GLOBALutils.invert(sci_S_ob_B)
	sci_S_co_B  = GLOBALutils.invert(sci_S_co_B)
	sci_S_ob_R  = GLOBALutils.invert(sci_S_ob_R)
	sci_S_co_R  = GLOBALutils.invert(sci_S_co_R)

        sci_S_ob_B,sci_S_co_B,sci_S_ob_R,sci_S_co_R = sci_S_ob_B[::-1],sci_S_co_B[::-1],sci_S_ob_R[::-1],sci_S_co_R[::-1]
        # save as fits file
        if (os.access(sci_fits_ob_B,os.F_OK)):
            os.remove( sci_fits_ob_B )
        if (os.access(sci_fits_co_B,os.F_OK)):
            os.remove( sci_fits_co_B )
        if (os.access(sci_fits_ob_R,os.F_OK)):
            os.remove( sci_fits_ob_R )
        if (os.access(sci_fits_co_R,os.F_OK)):
            os.remove( sci_fits_co_R )
            
        hdu = pyfits.PrimaryHDU( sci_S_ob_B )
        hdu.writeto( sci_fits_ob_B )
        hdu = pyfits.PrimaryHDU( sci_S_co_B )
        hdu.writeto( sci_fits_co_B )
        hdu = pyfits.PrimaryHDU( sci_S_ob_R )
        hdu.writeto( sci_fits_ob_R )
        hdu = pyfits.PrimaryHDU( sci_S_co_R )
        hdu.writeto( sci_fits_co_R )

    else:
        print fsim, "has already been extracted, reading in product fits files..."
        sci_S_ob_B = pyfits.getdata( sci_fits_ob_B )
        sci_S_co_B = pyfits.getdata( sci_fits_co_B )
        sci_S_ob_R = pyfits.getdata( sci_fits_ob_R )
        sci_S_co_R = pyfits.getdata( sci_fits_co_R )

    fout = 'proc/'+ obname + '_' + h[0].header['DATE-OBS'] + 'sp.fits'
    dateobs = h[0].header['DATE-OBS'][:4] + h[0].header['DATE-OBS'][5:7] + h[0].header['DATE-OBS'][8:9]

    #Build spectra
    
    if ( os.access(dirout+fout ,os.F_OK) == False ) or (force_spectral_file_build):
        # initialize file that will have the spectra
        spec = np.zeros((11, nord_ob2 + nord_ob1, dataB.shape[1]))
        hdu = pyfits.PrimaryHDU( spec )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MJD', mjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MBJD', mbjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START DATE', dateobs )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START UT',  h[0].header['UTC'] / 3600.)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TEXP (S)',h[0].header['EXPTIME'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH FLUX WEIGHTED MEAN F ',h[0].header['HIERARCH ESO INS DET1 TMMEAN'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (KM/S)', bcvel_baryc)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH (LAMBDA_BARY / LAMBDA_TOPO)', lbary_ltopo)    
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARGET NAME', obname)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA',h[0].header['RA'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC',h[0].header['DEC'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA BARY',ra)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC BARY',dec)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH EQUINOX',h[0].header['HIERARCH ESO TEL TARG EQUINOX'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LATITUDE',h[0].header['HIERARCH ESO TEL GEOLAT'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LONGITUDE',h[0].header['HIERARCH ESO TEL GEOLON'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS ALTITUDE',h[0].header['HIERARCH ESO TEL GEOELEV'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARG AIRMASS',h[0].header['HIERARCH ESO TEL AIRM START'])

	print '\t\tWavelength calibration:'
	print '\t\t\tComparision fibre is '+ cotype

        indice = sorted_indices[0]
        thar_fits_ob_R = dirout + ThAr_ref[indice].split('/')[-1][:-4]+'spec.ob.R.fits.S'
        thar_fits_co_R = dirout + ThAr_ref[indice].split('/')[-1][:-4]+'spec.co.R.fits.S'
        thar_fits_ob_B = dirout + ThAr_ref[indice].split('/')[-1][:-4]+'spec.ob.B.fits.S'
        thar_fits_co_B = dirout + ThAr_ref[indice].split('/')[-1][:-4]+'spec.co.B.fits.S'
        pkl_wsol = dirout + ThAr_ref[indice].split('/')[-1][:-4]+'wavsolpars.pkl'
        print "\t\t\tUnpickling wavelength solution from", pkl_wsol, " ..."
        wsol_dict = pickle.load(open(pkl_wsol,'r'))

	cotype = 'SKY'
	if cotype == 'WAVE':
		# Extract thAr lines from comparison orders
		lines_thar_co_R  = sci_S_co_R[:,1,:]
		iv_thar_co_R     = sci_S_co_R[:,2,:]
		lines_thar_co_B  = sci_S_co_B[:,1,:]
		iv_thar_co_B     = sci_S_co_B[:,2,:]

		Red_Pixel_Centers_co = np.array([])
		Red_Wavelengths_co   = np.array([])
		Red_Orders_co        = np.array([])
		Red_Centroids_co     = np.array([])
		Red_Sigmas_co        = np.array([])
		Red_Intensities_co   = np.array([])

		for order in range(nord_co2):
		    order_s = str(order)
		    if (order < 10):
		        order_s = '0'+str(order)

		    thar_order_orig = lines_thar_co_R[order,:]
		    IV              = iv_thar_co_R[order,:]
		    wei             = np.sqrt( IV )
		    #bkg             = CoralieUtils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
		    thar_order      = thar_order_orig #- bkg
		    
                    coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
		        = GLOBALutils.Initial_Wav_Calibration( order_dir+'R_order_'+order_s+final_wav, thar_order, order, wei, rmsmax=MRMS_initial, minlines=minlines_initial,  \
							      FixEnds=True,Dump_Argon=dumpargon, Dump_AllLines=True, Cheby=use_cheby, porder=porder)
		
		    Red_Pixel_Centers_co = np.append( Red_Pixel_Centers_co, pixel_centers )
		    Red_Wavelengths_co   = np.append( Red_Wavelengths_co, wavelengths )
		    Red_Orders_co        = np.append( Red_Orders_co, np.zeros( len(pixel_centers) ) + order )
		    Red_Centroids_co     = np.append( Red_Centroids_co, centroids)
		    Red_Sigmas_co        = np.append( Red_Sigmas_co, sigmas)
		    Red_Intensities_co   = np.append( Red_Intensities_co, intensities )

                p1_co_R, G_pix_co_R, G_ord_co_R, G_wav_co_R, II_co_R, rms_ms_co_R, G_res_co_R = \
		    GLOBALutils.Fit_Global_Wav_Solution(Red_Pixel_Centers_co, Red_Wavelengths_co, Red_Orders_co,\
		                                             np.ones(Red_Intensities_co.shape), wsol_dict['p1_co_R'], Cheby=use_cheby,\
		                                             maxrms=MRMS, Inv=Inverse_m, minlines=minlines_glob,\
                                                             order0=or0_R,ntotal=nord_co2,npix=len(thar_order),nx=ncoef_x_R,nm=ncoef_m_R)

                p_shift_co_R, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(G_pix_co_R, G_wav_co_R, G_ord_co_R,\
		                                                   np.ones(G_wav_co_R.shape), wsol_dict['p1_co_R'],\
		                                                   Cheby=True,Inv=True,maxrms=MRMS,minlines=minlines_glob,\
								   order0=or0_R,ntotal=nord_co2,npix=len(thar_order),nx=ncoef_x_R,nm=ncoef_m_R)

		weight_R = (np.sqrt(len(orders)) / rms_ms)**2


		Blue_Pixel_Centers_co = np.array([])
		Blue_Wavelengths_co   = np.array([])
		Blue_Orders_co        = np.array([])
		Blue_Centroids_co     = np.array([])
		Blue_Sigmas_co        = np.array([])
		Blue_Intensities_co   = np.array([])

		for order in range(nord_co1):
		    order = order + 1
		    order_s = str(order)
		    if (order < 10):
		        order_s = '0'+str(order)

		    thar_order_orig = lines_thar_co_B[order-1,:]
		    IV              = iv_thar_co_B[order-1,:]
		    wei             = np.sqrt( IV )       
		    thar_order      = thar_order_orig #- bkg
		    
                    coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms, residuals, centroids, sigmas, intensities \
		        = GLOBALutils.Initial_Wav_Calibration( order_dir+'B_order_'+order_s+final_wav, thar_order, order, wei, rmsmax=MRMS_initial, minlines=50,  \
							      FixEnds=True,Dump_Argon=dumpargon, Dump_AllLines=True, Cheby=use_cheby, porder=porder)
	
		    Blue_Pixel_Centers_co = np.append( Blue_Pixel_Centers_co, pixel_centers )
		    Blue_Wavelengths_co   = np.append( Blue_Wavelengths_co, wavelengths )
		    Blue_Orders_co        = np.append( Blue_Orders_co, np.zeros( len(pixel_centers) ) + order )
		    Blue_Centroids_co     = np.append( Blue_Centroids_co, centroids)
		    Blue_Sigmas_co        = np.append( Blue_Sigmas_co, sigmas)
		    Blue_Intensities_co   = np.append( Blue_Intensities_co, intensities )

                p1_co_B, G_pix_co_B, G_ord_co_B, G_wav_co_B, II_co_B, rms_ms_co_B, G_res_co_B = \
		    GLOBALutils.Fit_Global_Wav_Solution(Blue_Pixel_Centers_co, Blue_Wavelengths_co, Blue_Orders_co,\
		                                             np.ones(Blue_Intensities_co.shape), wsol_dict['p1_co_B'], Cheby=use_cheby,\
		                                             maxrms=MRMS, Inv=Inverse_m, minlines=minlines_glob,\
                                                             order0=or0_B,ntotal=nord_co1,npix=len(thar_order),nx=ncoef_x_B,nm=ncoef_m_B)

                p_shift_co_B, pix_centers, orders, wavelengths, I, rms_ms, residuals  = \
		    GLOBALutils.Global_Wav_Solution_vel_shift(G_pix_co_B, G_wav_co_B, G_ord_co_B,\
		                                                   np.ones(G_wav_co_B.shape), wsol_dict['p1_co_B'],\
		                                                   Cheby=True,Inv=True,maxrms=MRMS,minlines=minlines_glob,\
								   order0=or0_B,ntotal=nord_co1,npix=len(thar_order),nx=ncoef_x_B,nm=ncoef_m_B)
		weight_B = (np.sqrt(len(orders)) / rms_ms)**2

		shift = (p_shift_co_R[0]*weight_R + p_shift_co_B[0]*weight_B) / weight_R + weight_B
		print p_shift_co_R,p_shift_co_B,shift
	else:
		p_shift_co_R = [0.]
		p_shift_co_B = [0.]
		p_shift = 0.
		shift = 0.

	good_quality = True
        hdu = GLOBALutils.update_header(hdu,'HIERARCH THAR SHIFT_R',p_shift_co_R[0])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH THAR SHIFT_B',p_shift_co_B[0])
	hdu = GLOBALutils.update_header(hdu,'HIERARCH THAR SHIFT',shift)
        
        # Apply new wavelength solution including barycentric correction
        equis = np.arange( dataB.shape[1] )        

        for order in range(nord_ob2):
            m = order + or0_R
            chebs = GLOBALutils.Calculate_chebs(equis, m, Inverse=Inverse_m,order0=or0_R,ntotal=nord_ob2,npix=len(equis),nx=ncoef_x_R,nm=ncoef_m_R)
	    WavSol = lbary_ltopo * (1.0 + 1.0e-6*shift) * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wsol_dict['p1_R'],chebs,ncoef_x_R,ncoef_m_R)
            spec[0,order,:] = GLOBALutils.ToVacuum(WavSol)
	    spec[1,order,:] = sci_S_ob_R[order,1, :]
            spec[2,order,:] = sci_S_ob_R[order,2, :]
            fn = R_flat_ob[order,1,:]
            L  = np.where( fn > 0 )[0]
            spec[3,order,:][L] = sci_S_ob_R[order,1,:][L] / R_flat_ob[order,1,:][L]
            spec[4,order,:][L] = sci_S_ob_R[order,2,:][L] * ( R_flat_ob[order,1,:][L] ** 2 )

        for order in range(nord_ob1):
            m = order + or0_B
            chebs = GLOBALutils.Calculate_chebs(equis, m, Inverse=Inverse_m,order0=or0_B,ntotal=nord_ob1,npix=len(equis),nx=ncoef_x_B,nm=ncoef_m_B)
	    WavSol = lbary_ltopo * (1.0 + 1.0e-6*shift) * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wsol_dict['p1_B'],chebs,ncoef_x_B,ncoef_m_B)
            spec[0,order + nord_ob2,:] = GLOBALutils.ToVacuum(WavSol)
	    spec[1,order + nord_ob2,:] = sci_S_ob_B[order,1, :]
            spec[2,order + nord_ob2,:] = sci_S_ob_B[order,2, :]
            fn = B_flat_ob[order,1,:]
            L  = np.where( fn > 0 )[0]
            spec[3,order + nord_ob2,:][L] = sci_S_ob_B[order,1,:][L] / B_flat_ob[order,1,:][L]
            spec[4,order + nord_ob2,:][L] = sci_S_ob_B[order,2,:][L] * ( B_flat_ob[order,1,:][L] ** 2 )

	ccoefs = GLOBALutils.get_cont(spec[0],spec[3])

	for order in range(nord_ob2):
            L  = np.where( spec[1,order,:] != 0 )
            spec[5,order,:][L] = spec[3,order,:][L] / np.polyval(ccoefs[order],spec[0,order,:][L])    
            ratio              = np.polyval(ccoefs[order],spec[0,order,:][L])*Rnorms[order]
	    fn = R_flat_ob_n[order,1,:]
            L  = np.where( fn > 0 )[0]
            spec[3,order,:][L] = sci_S_ob_R[order,1,:][L] / R_flat_ob_n[order,1,:][L]
            spec[4,order,:][L] = sci_S_ob_R[order,2,:][L] * ( R_flat_ob_n[order,1,:][L] ** 2 )
            spec[6,order,:][L] = spec[4,order,:][L] * (ratio ** 2 )
            spec[7,order,:][L] = ratio
            spec[8,order,:][L] = ratio * R_flat_ob_n[order,1,:][L] / np.sqrt( ratio * R_flat_ob_n[order,1,:][L] / gain2 + (ron2/gain2)**2 )
            spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
            dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
            NN            = np.average(dlambda_dx)
            dlambda_dx    /= NN
	    LL = np.where(spec[5,order,:] > 1 + 10. / scipy.signal.medfilt(spec[8,order,:],21))[0]
	    spec[5,order,LL] = 1.
            spec[9,order,:][L] = spec[5,order,:][L] * (dlambda_dx[L] ** 1) 
            spec[10,order,:][L] = spec[6,order,:][L] / (dlambda_dx[L] ** 2)

	for order in range(nord_ob1):
            L  = np.where( spec[1,order + nord_ob2,:] != 0 )
            spec[5,order + nord_ob2,:][L] = spec[3,order + nord_ob2,:][L] / np.polyval(ccoefs[order + nord_ob2],spec[0,order + nord_ob2,:][L])	   
            ratio              = np.polyval(ccoefs[order + nord_ob2],spec[0,order + nord_ob2,:][L])*Bnorms[order]
	    fn = B_flat_ob_n[order,1,:]
            L  = np.where( fn > 0 )
            spec[3,order + nord_ob2,:][L] = sci_S_ob_B[order,1,:][L] / B_flat_ob_n[order,1,:][L]
            spec[4,order + nord_ob2,:][L] = sci_S_ob_B[order,2,:][L] * ( B_flat_ob_n[order,1,:][L] ** 2 )
            spec[6,order + nord_ob2,:][L] = spec[4,order + nord_ob2,:][L] * (ratio ** 2 )
            spec[7,order + nord_ob2,:][L] = ratio
            spec[8,order + nord_ob2,:][L] = ratio * B_flat_ob_n[order,1,:][L] / np.sqrt( ratio * B_flat_ob_n[order,1,:][L] / gain1 + (ron1/gain1)**2 )
            spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
            dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
            NN            = np.average(dlambda_dx)
            dlambda_dx    /= NN	   
	    LL = np.where(spec[5,order + nord_ob2,:] > 1 + 20./scipy.signal.medfilt(spec[8,order + nord_ob2,:],21))[0]
	    spec[5,order + nord_ob2,LL] = 1.
            spec[9,order + nord_ob2,:][L] = spec[5,order + nord_ob2,:][L] * (dlambda_dx[L] ** 1) 
            spec[10,order + nord_ob2,:][L] = spec[6,order + nord_ob2,:][L] / (dlambda_dx[L] ** 2)

	JustExtract = False
        if (not JustExtract):
	    if DoClass:
		print '\t\tSpectral Analysis:'
                # spectral analysis
                query_success,sp_type_query = GLOBALutils.simbad_query_obname(obname)
		# Now, query SIMBAD by coordinates if above not successful
		if (not query_success):
		    query_success,sp_type_query = GLOBALutils.simbad_query_coords('12:00:00','00:00:00')
		print "\t\t\tSpectral type returned by SIMBAD query:",sp_type_query

                hdu = GLOBALutils.update_header(hdu,'HIERARCH SIMBAD SPTYP', sp_type_query)

                pars_file = dirout + fsim.split('/')[-1][:-4]+'_stellar_pars.txt'

		if os.access(pars_file,os.F_OK) == False or force_stellar_pars:
			print "\t\t\tEstimating atmospheric parameters:"
			Rx = np.around(1./np.sqrt(1./40000.**2 - 1./RESI**2))
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
            ml, mh, weight = np.loadtxt(mask,unpack=True)
            ml_v = GLOBALutils.ToVacuum( ml )
            mh_v = GLOBALutils.ToVacuum( mh )
	    av_m = 0.5*( ml_v + mh_v )
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

	    if sp_type == 'M5':
		moon_sig = 4.5
	    elif sp_type == 'K5':
		moon_sig = 4.2
	    else:
		moon_sig = 4.0

            while (cond):
                # first rough correlation to find the minimum
                vels, xc_full, sn, nlines_ccf, W_ccf = \
					GLOBALutils.XCor(spec, ml_v, mh_v, weight,\
					0, lbary_ltopo, vel_width=300, vel_step=3,\
					spec_order=9, iv_order=10, sn_order=8,max_vel_rough=300)

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
                
                rvels, rxc_av, rpred, rxc_av_orig, rvel0_xc = \
					vels.copy(), xc_av.copy(), pred.copy(),\
					xc_av_orig.copy(), vel0_xc

		xc_av_rough = xc_av
                vels_rough  = vels
                
                vel_width = np.maximum( 20.0, 6*disp )
                vels, xc_full, sn, nlines_ccf, W_ccf =\
					GLOBALutils.XCor(spec, ml_v, mh_v, weight,\
					vel0_xc, lbary_ltopo, vel_width=vel_width,\
					vel_step=0.1, spec_order=9, iv_order=10, sn_order=8,max_vel_rough=300)

		xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3.0, Simple=True, W=W_ccf)

                pred = np.array( approx(Temp[0],Temp[1],xout=vels, method="linear", rule=2) )[1]
                xc_av /= pred

		p1,XCmodel,p1gau,XCmodelgau,Ls2 = \
					GLOBALutils.XC_Final_Fit( vels, xc_av, sigma_res = 4,\
					 horder=8, moonv=refvel, moons=moon_sig, moon=False)

		moonmatters = False
		if (know_moon and here_moon):
			moonmatters = True
			ismoon = True
			confused = False
			p1_m,XCmodel_m,p1gau_m,XCmodelgau_m,Ls2_m = GLOBALutils.XC_Final_Fit( vels, xc_av, \
				sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = True)
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
			 'lunation':lunation,'mephem':mephem,'texp':h[0].header['EXPTIME']}

	    pkl_xc = dirout + fsim.split('/')[-1][:-4]+obname+'_XC_'+sp_type+'.pkl'
            pickle.dump( xc_dict, open( pkl_xc, 'w' ) )

	    ccf_pdf = dirout + 'proc/' + fsim.split('/')[-1][:-4] + obname + '_XCs_' + sp_type + '.pdf'

	    if not avoid_plot:
	        GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)

            SNR_5130 = np.median(spec[8,28,1900:2101] ) 
            airmass  = h[0].header['HIERARCH ESO TEL AIRM START']
            seeing   = h[0].header['HIERARCH ESO TEL AMBI FWHM START']

	    B,A = -0.00257864,0.07765779
	    RVerr  =  B + ( 1.6 + 0.2 * p1gau[2] ) * A / np.round(SNR_5130)
	    depth_fact = 1. + p1gau[0]/(p1gau[2]*np.sqrt(2*np.pi))
	    if depth_fact < 0.6:
		depth_fact = 0.6
	    depth_fact = (1 - 0.6) / (1 - depth_fact)
	    RVerr *= depth_fact
	    if RVerr < 0.002:
	    	RVerr = .002

	    B,A = -0.00348879, 0.10220848
	    BSerr = B + ( 1.6 + 0.2 * p1gau[2] ) * A / np.round(SNR_5130)
	    if BSerr<0.002:
		BSerr = .002

	    RV     = np.around(p1gau_m[1],4)  
	    BS     = np.around(SP,4)   
            RVerr2 = np.around(RVerr,4)
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
	    SNR_5130_R = np.around(SNR_5130*np.sqrt(3.0))

	    disp_epoch = np.around(p1gau_m[2],1)
            hdu = GLOBALutils.update_header(hdu,'RV', RV)
            hdu = GLOBALutils.update_header(hdu,'RV_E', RVerr2)
            hdu = GLOBALutils.update_header(hdu,'BS', BS)
            hdu = GLOBALutils.update_header(hdu,'BS_E', BSerr)
            hdu = GLOBALutils.update_header(hdu,'DISP', disp_epoch)
            hdu = GLOBALutils.update_header(hdu,'SNR', SNR_5130)
            hdu = GLOBALutils.update_header(hdu,'SNR_R', SNR_5130_R)
	    hdu = GLOBALutils.update_header(hdu,'INST', 'HARPS')
	    hdu = GLOBALutils.update_header(hdu,'RESOL', RESI)
	    hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')
	    hdu = GLOBALutils.update_header(hdu,'XC_MIN', XC_min)
	    hdu = GLOBALutils.update_header(hdu,'BJD_OUT', bjd_out)
    
            line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f     harps   ceres   %8d %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
                      (obname, bjd_out, RV, RVerr2, BS, BSerr, RESI, T_eff_epoch, logg_epoch, Z_epoch, vsini_epoch, XC_min, disp_epoch,\
		       TEXP, SNR_5130_R, ccf_pdf)
	    f_res.write(line_out)
        if (os.access( dirout + fout,os.F_OK)):
            os.remove( dirout + fout)
        hdu.writeto( dirout + fout )

    else:
        print "\t\tReading spectral file from", fout
        spec = pyfits.getdata( fout )

f_res.close()
