import sys
import matplotlib
matplotlib.use("Agg")
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
 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ceres modules
import vbtutils
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
from astropy.io import fits as pyfits
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
parser.add_argument('-resolution',default='60000')
parser.add_argument('-binning', default='1x1')
parser.add_argument('-gangle', default=71.726)


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
binning              = args.binning
binx                 = int(binning.split('x')[0])
biny                 = int(binning.split('x')[1])
gangle               = float(args.gangle)
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

####### GLOBAL VARIABLES #####
force_pre_process  = False
force_bl           = False  #
force_bkg          = False
force_P            = False
force_thar_extract = False
force_thar_wavcal  = True 
force_sci_extract  = False
force_sci_proc     = False  #
force_RV           = False  #
force_stellar_pars = False
force_corr         = False
force_flat_extract = False

bad_colummn        = True

Inverse_m          = True
use_cheby          = True

MRMS              = 80
trace_degree      = 5
Marsh_alg         = 0
ext_aperture      = int(np.around(10. / float(biny)))
NSigma_Marsh      = 5
NCosmic_Marsh     = 5
S_Marsh           = 0.4
N_Marsh           = 3
min_extract_col   = int(np.around(0. / float(binx)))
max_extract_col   = int(np.around(4096. / float(binx)))

#npar_wsol = 27 #number of parameters of wavelength solution
ncoef_x   = 4
ncoef_m   = 8
npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2
n_useful = 50
oro0      = 39
bsec      = [0,50,4146,4196]

order_dir   = base+'vbt/wavcals/'
sufix       = '.iwdat'

models_path = base+'data/COELHO_MODELS/R_40000b/'

print "\n\n\tEchelle Vainu Bappu 2.34m Telescope PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print '\n'

# file containing the log
log = dirout+'night.log'
biases, flats, objects, ThAr_ref, darks = vbtutils.FileClassify(dirin,log)

hd = pyfits.getheader(ThAr_ref[0])
if hd['GRATANGL'] == 71.726:
    sufix = '.2.iwdat'
elif hd['GRATANGL'] == 72.429:
    sufix = '.iwdat'

#ThAr_ref = ThAr_ref[:2]
if dark_substraction == True and len(darks)<3:
    dark_substraction = False
f = open(log,'r')
lines = f.readlines()

print '\tThese are all the images to proccess:'
for bias in biases:
    hd = pyfits.getheader(bias)
    print '\tbias', hd['IMAGETYP'], hd['IMAGETYP'], hd['DATE-OBS'],bias
print '\n'
for dark in darks:
    hd = pyfits.getheader(dark)
    print '\tdark', hd['IMAGETYP'], hd['IMAGETYP'], hd['DATE-OBS'],dark
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

if ( (os.access(dirout+'MasterFlat.fits',os.F_OK) == False) or\
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
    MasterBias, RO_bias, GA_bias = vbtutils.MedianCombine(biases, zero_bo=False, dark_bo=False, flat_bo=False)
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
            DARK, RON, GAIN = vbtutils.MedianCombine(dark_groups[i], zero_bo=True, zero=dirout+'MasterBias.fits',dark_bo=False)
            hdu = pyfits.PrimaryHDU( DARK )
            hdu = GLOBALutils.update_header(hdu,'EXPTIME',dark_times[i])
            if os.access(dirout+'DARK_'+str(int(dark_times[i]))+'s.fits',os.F_OK):
                os.remove(dirout+'DARK_'+str(int(dark_times[i]))+'s.fits')
            hdu.writeto( dirout+'DARK_'+str(int(dark_times[i]))+'s.fits' )
            MDARKS.append(dirout+'DARK_'+str(int(dark_times[i]))+'s.fits')
            i+=1
    print "\t\t-> Masterdarks: done!"

    Flat, RO_flat, GA_flat = vbtutils.MedianCombine(flats, zero_bo=True, dark_bo=False, flat_bo=False,zero=dirout+'MasterBias.fits')
    hdu = pyfits.PrimaryHDU( Flat )

    if (os.access(dirout+'MasterFlat.fits',os.F_OK)):
        os.remove(dirout+'MasterFlat.fits')
    hdu.writeto(dirout+'MasterFlat.fits')
    print "\t\t-> Masterflat: done!"

    # Find orders & traces
    print "\tTracing echelle orders..."
    h = pyfits.open(dirin+stst)[0]
    hth = pyfits.getheader(dirin+stst)
    d = h.data[0]
    d = vbtutils.OverscanTrim(d,bsec)
    d -= MasterBias
    c_all, nord = GLOBALutils.get_them(d,20,trace_degree,mode=1,nsigmas=3,endat=4100)
    print '\t\t'+str(nord)+' orders found...'
    
    trace_dict = {'c_all':c_all, 'nord':nord, 'DARKS':MDARKS, 'dtimes':dark_times, 'RO_flat':RO_flat, 'GA_flat':GA_flat}
    pickle.dump( trace_dict, open( dirout+"trace.pkl", 'w' ) )


else:
    trace_dict = pickle.load( open( dirout+"trace.pkl", 'r' ) )

    c_all = trace_dict['c_all']
    nord = trace_dict['nord']
    RO_flat = trace_dict['RO_flat']
    GA_flat = trace_dict['GA_flat']

    h = pyfits.open(dirout+'MasterBias.fits')
    MasterBias = h[0].data

    MDARKS = trace_dict['DARKS']
    dark_times = trace_dict['dtimes']

    h = pyfits.open(dirout+'MasterFlat.fits')
    Flat = h[0].data

print '\n\tExtraction of Master Flat:'
flat_simple_fits = dirout + 'Flat.spec.simple.fits'
flat_fits        = dirout + 'Flat.spec.fits'
P_fits           = dirout + 'P.fits'
if ( os.access(flat_simple_fits,os.F_OK) == False ) or ( os.access(flat_fits,os.F_OK) == False ) or force_flat_extract: 
    Centers = np.zeros((len(c_all),Flat.shape[1]))
    for i in range(nord):
        Centers[i,:]=scipy.polyval(c_all[i,:],np.arange(len(Centers[i,:])))
    bkg = GLOBALutils.get_scat(Flat,Centers,span=15)
    Flat -= bkg
    flat_simple = GLOBALutils.simple_extraction(Flat,c_all,ext_aperture,min_extract_col,max_extract_col,npools)
    P = GLOBALutils.obtain_P(Flat,c_all,ext_aperture,RO_flat,GA_flat,NSigma_Marsh,S_Marsh,N_Marsh,Marsh_alg,min_extract_col,max_extract_col,npools)
    flat        = GLOBALutils.optimal_extraction(Flat,P,c_all,ext_aperture,RO_flat,GA_flat,S_Marsh,NCosmic_Marsh,min_extract_col,max_extract_col,npools)
    flat_simple = flat_simple[::-1]
    flat        = flat[::-1]
    if (os.access(flat_simple_fits,os.F_OK)):
        os.remove( flat_simple_fits )
    hdu = pyfits.PrimaryHDU( flat_simple )
    hdu.writeto( flat_simple_fits )
    if (os.access(flat_fits,os.F_OK)):
        os.remove( flat_fits )
    hdu = pyfits.PrimaryHDU( flat )
    hdu.writeto( flat_fits )

    if (os.access(P_fits,os.F_OK)):
        os.remove( P_fits )
    hdu = pyfits.PrimaryHDU( P )
    hdu.writeto( P_fits )

else:
    flat_simple = pyfits.getdata(flat_simple_fits)
    flat        = pyfits.getdata(flat_fits)

print '\n\tExtraction of ThAr calibration frames:'

for fsim in ThAr_ref:
    hth = pyfits.getheader(fsim)
    thmjd,mjd0 = vbtutils.mjd_fromheader(hth)
    dth = pyfits.getdata(fsim)[0]
    dth = vbtutils.OverscanTrim(dth,bsec)
    dth -= MasterBias

    thar_fits_simple = dirout+'ThAr_'+hth['DATE-OBS']+'.spec.simple.fits'

    if ( os.access(thar_fits_simple,os.F_OK) == False ) or (force_thar_extract): 
        print "\t\tNo previous extraction or extraction forced for ThAr file", fsim, "extracting..."
        Centers = np.zeros((len(c_all),dth.shape[1]))
        for i in range(nord):
            Centers[i,:]=scipy.polyval(c_all[i,:],np.arange(len(Centers[i,:])))
        bkg = GLOBALutils.get_scat(dth,Centers,span=15)
        dth -= bkg

        thar_Ss = GLOBALutils.simple_extraction(dth,c_all,ext_aperture,min_extract_col,max_extract_col,npools)
        thar_Ss = thar_Ss[::-1]
        if (os.access(thar_fits_simple,os.F_OK)):
            os.remove( thar_fits_simple )
        hdu = pyfits.PrimaryHDU( thar_Ss )
        hdu.writeto( thar_fits_simple )
    else:
        print "\t\tThAr file", fsim, "all ready extracted, loading..."


print "\n\tWavelength solution of ThAr calibration spectra:"
# compute wavelength calibration files

#force_thar_wavcal = True
for thar in ThAr_ref:

    hth = pyfits.getheader(thar)
    thmjd,mjd0 = vbtutils.mjd_fromheader(hth)
    thar_fits_simple = dirout+'ThAr_'+hth['DATE-OBS']+'.spec.simple.fits'
    wavsol_pkl       = dirout+'ThAr_'+hth['DATE-OBS']+'.wavsolpars.pkl'
    thar_Ss = pyfits.getdata(thar_fits_simple)
    RON = hth['RDNOISE']
    GAIN = hth['GAIN']

    if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
        print " \t\tWorking on ThAr file", thar 
        lines_thar = thar_Ss[:,:]

        orders_offset, rough_shift = vbtutils.get_thar_offsets(lines_thar,binning=1,pref='order_',suf=sufix)
        print 'orders_ofset:',orders_offset
        print 'rough_shift:',rough_shift
        
        orderi = 0
        if orders_offset < 0:
            orderi = - orders_offset
        orderf = nord - 1
        if orderf + orders_offset >= n_useful:
            orderf = n_useful - orders_offset - 1

        iv_thar = 1/((lines_thar/GAIN) + (RON**2/GAIN**2))

        All_Pixel_Centers = np.array([])
        All_Wavelengths   = np.array([])
        All_Orders        = np.array([])
        All_Centroids     = np.array([])
        All_Sigmas        = np.array([])
        All_Intensities   = np.array([])
        All_Residuals     = np.array([])
        All_Sigmas        = np.array([])
    
        for order in range(orderi,orderf+1):
            order_s = str(order)
            if (order < 10):
                order_s = '0'+str(order)
            thar_order_orig = lines_thar[order,:]
            IV              = iv_thar[order,:]
            wei             = np.sqrt( IV )
            bkg             = scipy.signal.medfilt(thar_order_orig,201)
            thar_order      = thar_order_orig - thar_order_orig.min()
            f = open(order_dir+'order_'+order_s+sufix)
            lines = f.readlines()
            if len(lines)>5:
                coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths, rms_ms,\
                residuals, centroids, sigmas, intensities = GLOBALutils.Initial_Wav_Calibration( \
                                  order_dir+'order_'+order_s+sufix, thar_order, order, wei,\
                                  rmsmax=200, minlines=10, FixEnds=False, Dump_Argon=False, \
                                  Dump_AllLines=True, Cheby=use_cheby, rough_shift=rough_shift,do_xc=False)
                """
                fwhms_lns = sigmas*2.355
                inis_lns  = pixel_centers - fwhms_lns*0.5
                fins_lns  = pixel_centers + fwhms_lns*0.5           
                inis_wvs  = GLOBALutils.Cheby_eval(coeffs_pix2wav,inis_lns,float(len(thar_order)))
                fins_wvs  = GLOBALutils.Cheby_eval(coeffs_pix2wav,fins_lns,float(len(thar_order)))
                fwhms_wvs = inis_wvs - fins_wvs
                resolution2 = wavelengths / fwhms_wvs

                print "\t\t\tmedian Resolution of order", order, '=', np.around(np.median(resolution2))
                #"""
                if (order == int(0.5*n_useful)): 
                    if (use_cheby):
                        Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 0.5*len(thar_order), len(thar_order) )
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
            
        p0    = np.zeros( npar_wsol )
        p0[0] = int(.5*n_useful) * Global_ZP

        p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
            GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders, \
            np.ones(All_Intensities.shape), p0, npix=lines_thar.shape[1], Cheby=use_cheby, \
            maxrms=MRMS, Inv=Inverse_m, minlines=1000,order0=oro0,nx=ncoef_x,nm=ncoef_m,ntotal=n_useful)
        """
        fl = open('../fies/lovis.txt')
        lls, lts,sts = [],[],[]
        lines = fl.readlines()
        for line in lines:
            cos = line.split()
            lls.append(float(cos[0]))
            lts.append(cos[3])
            sts.append(float(cos[2]))
        lls,lts,sts = np.array(lls),np.array(lts),np.array(sts)
        equis = np.arange(4096)
        for i in np.arange(orderi,orderf+1):
            if i > 14:
                m = i + oro0
                chebs   = GLOBALutils.Calculate_chebs(equis, m, order0=oro0, ntotal=n_useful, npix=len(equis), Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
                WavSol  = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(p1,chebs,ncoef_x,ncoef_m)
                tck = scipy.interpolate.splrep(WavSol,equis,k=3)
                px  = interpolate.splev(lls,tck)
                J = np.where((px>0)&(px<1000))[0]
                print J
                px = px[J]
                lt = lts[J]
                st = sts[J]
                ls = lls[J]
                J = np.where(st>5000)[0]
                px,lt,st, ls = px[J],lt[J],st[J],ls[J]

                if i<10:
                    si = '0'+str(int(i))
                else:
                    si = str(int(i)) 
                fr = open('wavcals/order_'+si+sufix,'r')
                tmplines = fr.readlines()
                f = open('wavcals/order_'+si+sufix,'w')


                imin = False
                vec = []
                for j in np.arange(len(px)-1):
                    if px[j+1]-px[j] > 15:
                        vec.append(j)
                        imin = False
                    else:
                        vec.append(j)
                        imin = True

                    if not imin:
                        line = str(len(vec))+'\t'
                        for k in vec:
                            line+= str(px[k]) + '\t' + str(ls[k]) + '\t'
                        for k in vec:
                            line+= str(lt[k]) + '\t'
                        line+= '\n'
                        f.write(line)
                        imin = False
                        vec = []
                for line in tmplines:
                    f.write(line)
                f.close()
        """


        #    I = np.where(G_ord==i)[0]
        #    plot(G_wav[I],G_res[I],'.')
        #    plot(np.median(G_wav[I]),np.median(G_res[I]),'ko')
        #axhline(0)
        #show()
        pdict = {'p1':p1, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II,\
                     'rms_ms':rms_ms, 'G_res':G_res, 'All_Centroids':All_Centroids,\
                     'All_Sigmas':All_Sigmas, 'orders_offset':orders_offset,\
                     'rough_shift':rough_shift}
        pickle.dump( pdict, open( wavsol_pkl, 'w' ) )
    else:
        print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl
        pdict = pickle.load(open(wavsol_pkl, 'r'))

print '\n\tProcessing of science images:'

new_list         = []
new_list_obnames = []
new_list_texp    = []
for i in range(len(objects)):
    fsim = objects[i]
    hd = pyfits.getheader(objects[i])
    obname = hd['OBJECT']
    texp   = hd['EXPTIME']
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

for obj in objects:

    print '\n'
    print "\t--> Working on image: ", obj

    hd = pyfits.getheader(obj)
    nombre    = hd['OBJECT'].replace(' ','')
    RON, GAIN = hd['RDNOISE'], hd['GAIN']

    print "\t\tObject name:",nombre
    
    nama = nombre + '_' + hd['DATE-OBS'] 

    obj_fits        = dirout + nama + '.spec.fits.S'
    obj_fits_simple = dirout + nama + '.spec.simple.fits.S'
    bkg_obj_fits    = dirout + 'Bkg_' + nama + '.fits'  
    P_fits          = dirout + 'P_' + nama + '.fits'

    if ( os.access(obj_fits,os.F_OK) == False )  or\
       ( os.access(obj_fits_simple,os.F_OK) == False ) or\
       (force_sci_extract) or ( os.access(P_fits,os.F_OK) == False ):

        print "\t\tNo previous extraction or extraction forced for science file", obj, "extracting..."

        dat = pyfits.getdata(obj)[0]
        hdt = pyfits.getheader(obj)
        dat = vbtutils.OverscanTrim(dat,bsec)
        dat -= MasterBias
        c_alls, pshift = GLOBALutils.retrace( dat, c_all,span=30 )
        Centers = np.zeros((len(c_alls),dat.shape[1]))
        for i in range(nord):
            Centers[i,:]=scipy.polyval(c_alls[i,:],np.arange(len(Centers[i,:])))

        #print 'Scatter Light Determination...'

        if ( os.access(bkg_obj_fits,os.F_OK) == False or force_bkg):
            bkg = GLOBALutils.get_scat(dat,Centers,span=15)
            if (os.access(bkg_obj_fits,os.F_OK)):
                os.remove( bkg_obj_fits )
            hdu = pyfits.PrimaryHDU( bkg )
            hdu.writeto( bkg_obj_fits )
        else:
            bkg = pyfits.getdata(bkg_obj_fits)
        dat -= bkg

        if os.access(P_fits,os.F_OK) == False or force_P:
            P = GLOBALutils.obtain_P(dat,c_alls,ext_aperture,RON,\
                                GAIN,NSigma_Marsh,S_Marsh,N_Marsh,Marsh_alg,\
                                min_extract_col,max_extract_col,npools)
            if (os.access(P_fits,os.F_OK)):
                os.remove( P_fits )
            hdu = pyfits.PrimaryHDU( P )
            hdu.writeto( P_fits )
        else:
            P = pyfits.getdata(P_fits)

        obj_S  = GLOBALutils.optimal_extraction(dat,P,c_alls,ext_aperture,RON,\
                                 GAIN,S_Marsh,NCosmic_Marsh,min_extract_col,max_extract_col,npools)[::-1]

        if (os.access(obj_fits,os.F_OK)):
            os.remove( obj_fits )
        hdu = pyfits.PrimaryHDU( obj_S )
        hdu.writeto( obj_fits )

        obj_Ss = GLOBALutils.simple_extraction(dat,c_alls,ext_aperture,min_extract_col,max_extract_col,npools)[::-1]
        
        if (os.access(obj_fits_simple,os.F_OK)):
            os.remove( obj_fits_simple )
        hdu = pyfits.PrimaryHDU( obj_Ss )
        hdu.writeto( obj_fits_simple )
#print gfd

#################################################################################################################
############################################  Final output ######################################################
#################################################################################################################

print "\n\tBuilding the final output spectra..."

for obj in objects:
    hd     = pyfits.getheader(obj)
    nombre = hd['OBJECT'].replace(' ','')
    nama = nombre + '_' + hd['DATE-OBS'] 
    nf     = nama+'_final.fits'
    print "\n\t\t-->Building", nf
    print dirout+'proc/'+nf
    if os.access(dirout+'proc/'+nf,os.F_OK) == False or force_sci_proc:
        # Get observing info from header
        hd = pyfits.getheader(obj)
        nombre    = hd['OBJECT'].replace(' ','')
        exptime   = hd['EXPTIME']
        try:
            RA        = hd['TELRA']
        except:
            RA        = hd['RA']
        try:
            DEC       = hd['TELDEC']
        except:
            DEC       = hd['DEC']
        RON       = hd['RDNOISE']
        GAIN      = hd['GAIN']
        altitude  = 725.
        latitude  = 12.577
        longitude = 78.827
        epoch     = 2000.

        cos = RA.replace(' ','').split(':')
        if float(cos[0])>=0:
            RA = float(cos[0])+float(cos[1])/60.+float(cos[2])/3600.
            RA = RA*360./24.
        else:
            RA = -float(cos[0])+float(cos[1])/60.+float(cos[2])/3600.
            RA = -RA*360./24.

        cos = DEC.replace(' ','').split(':')
        if float(cos[0])>=0:
            DEC = float(cos[0])+float(cos[1])/60.+float(cos[2])/3600.
        else:
            DEC = -float(cos[0])+float(cos[1])/60.+float(cos[2])/3600.

        scmjd,scmjd0 = vbtutils.mjd_fromheader(hd)
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
        bcvel_baryc = ( lbary_ltopo - 1.0 ) * 2.99792458E5  #This in the barycentric velocity
        res         = jplephem.pulse_delay(RA/15.0, DEC, int(scmjd), scmjd%1, 1, 0.0)
        scmbjd      = scmjd + res['delay'][0] / (3600.0 * 24.0) #This is the modified barycentric julian day of the observation

        # set observatory info to retrive info about the moon
        gobs = ephem.Observer()
        gobs.name = 'VBT'
        gobs.lat  = rad(latitude)
        gobs.long = rad(longitude)
        #gobs.date = hd['UT-DATE'] + ' ' + hd['UT-TIME'].replace(':','_')
        gobs.date = hd['DATE-OBS'].replace('T',' ')

        mephem = ephem.Moon()
        mephem.compute(gobs)
        Mcoo = jplephem.object_track("Moon", int(scmjd), float(scmjd%1), 1, 0.0)
        Mp   = jplephem.barycentric_object_track("Moon", int(scmjd), float(scmjd%1), 1, 0.0)
        Sp   = jplephem.barycentric_object_track("Sun", int(scmjd), float(scmjd%1), 1, 0.0)
        res      = jplephem.object_doppler("Moon", int(scmjd), scmjd%1, 1, 0.0)
        lunation,moon_state,moonsep2,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,RA,DEC)
        refvel = bcvel_baryc + moonvel  #This is the velocity of the spectrum of the moon with the applied barycentric correction in the direction of the target. 

        print '\t\t\tBarycentric velocity:',refvel
        
        obj_fits        = dirout+nama+'.spec.fits.S'
        obj_fits_simple = dirout+nama+'.spec.simple.fits.S'
        obj_S           = pyfits.getdata(obj_fits)
        obj_Ss          = pyfits.getdata(obj_fits_simple)

        hth = pyfits.getheader(ThAr_ref[0])
        wavsol_pkl       = dirout+'ThAr_'+hth['DATE-OBS']+'.wavsolpars.pkl'
        pdict      = pickle.load(open(wavsol_pkl,'r'))
        global1     = pdict['p1']

        orders_offset = pdict['orders_offset']
        orderi = 0
        if orders_offset < 0:
            orderi = - orders_offset
        orderf = nord - 1
        if orderf + orders_offset >= n_useful:
            orderf = n_useful - orders_offset - 1


        final  = np.zeros( [11, orderf - orderi + 1, np.shape(obj_S)[2]] )
        equis  = np.arange( obj_S.shape[2] ) 

        for order in range(orderi,orderf+1):
            m = order + oro0
            chebs = GLOBALutils.Calculate_chebs(equis, m, order0=oro0, ntotal=n_useful, npix=obj_S.shape[2], Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
            WavSol  = lbary_ltopo*(1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(global1,chebs,ncoef_x,ncoef_m)

            final[0,order,:] = WavSol
            final[1,order,:] = obj_S[order,1,:]
            final[2,order,:] = obj_S[order,2,:]
            final[3,order,:] = final[1,order,:]/flat[order,1]
            final[4,order,:] = final[2,order,:]*(flat[order,1]**2)

            ccoef = GLOBALutils.get_cont_single(final[0,order],final[3,order],final[4,order],ll=1.5,lu=5,nc=3)
            L  = np.where( final[1,order] != 0 )

            ratio = np.polyval(ccoef,final[0,order])
            final[5,order,:] = final[3,order,:]/ratio
            Inan = np.where( np.isnan(final[1,order,:]) == True )[0]
            final[5,order,Inan] = 1.
            final[6,order,:] = final[4,order,:]*(ratio**2)
            final[7,order,:] = ratio
            final[8,order,:] = ratio*flat[order,1] / np.sqrt( ratio *flat[order,1] / GAIN + (RON/GAIN)**2 )

            spl        = scipy.interpolate.splrep(np.arange(len(final[0,order,:])),final[0,order,:] ,k=3)
            dlambda_dx = scipy.interpolate.splev(np.arange(len(final[0,order,:])), spl, der=1)
            NN         = np.average(dlambda_dx)
            dlambda_dx /= NN

            LL = np.where(final[5,order] > 1 + 10. / scipy.signal.medfilt(final[8,order],21))[0]
            final[5,order,LL] = 1.
            final[9,order][L] = final[5,order][L] * (dlambda_dx[L] ** 1) 
            final[10,order][L] = final[6,order][L] / (dlambda_dx[L] ** 2)

        hdu = pyfits.PrimaryHDU( final )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MJD',scmjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MBJD',scmbjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START DATE',hd['DATE-OBS'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TEXP (S)',hd['EXPTIME'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (KM/S)',bcvel_baryc,'[km/s]')
        hdu = GLOBALutils.update_header(hdu,'HIERARCH (LAMBDA_BARY / LAMBDA_TOPO)',lbary_ltopo)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARGET NAME',nombre)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA',RA)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC',DEC)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA BARY',RA)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC BARY',DEC)            
        hdu = GLOBALutils.update_header(hdu,'HIERARCH EQUINOX',2000)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LATITUDE',latitude)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LONGITUDE',longitude)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS ALTITUDE',altitude)
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
        hdo     = pyfits.getheader(obj)
        nombre = hdo['OBJECT'].replace(' ','')
        nama = nombre + '_' + hdo['DATE-OBS'] 
        fit    = dirout + 'proc/' + nama + '_final.fits'
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
            if spec[0,i,1000] < 5200:
                tuc = i
                break
        SNR_5130 = np.median(spec[8,tuc][1800:2200] )

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

            pars_file = dirout + nombre+'_'+hdo['DATE-OBS']+'_stellar_pars.txt'

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
        hdu[0] = GLOBALutils.update_header(hdu[0],'INST', 'VBT')
        hdu[0] = GLOBALutils.update_header(hdu[0],'RESOL', resolution)
        hdu[0] = GLOBALutils.update_header(hdu[0],'PIPELINE', 'CERES')
        hdu[0] = GLOBALutils.update_header(hdu[0],'XC_MIN', XC_min)
        hdu[0] = GLOBALutils.update_header(hdu[0],'BJD_OUT', bjd_out)

        line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f   vbt   ceres   %8d %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
                      (obname, bjd_out, RV, RVerr2, BS, BSerr, resolution, T_eff_epoch, logg_epoch,\
               Z_epoch, vsini_epoch, XC_min, disp_epoch, TEXP, SNR_5130_R, ccf_pdf)
        f_res.write(line_out)
        hdu.close()
f_res.close()
