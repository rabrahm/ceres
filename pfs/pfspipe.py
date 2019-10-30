import sys
import matplotlib
matplotlib.use("Agg")
from pylab import *
base = '../'

sys.path.append(base+"utils/Correlation")
sys.path.append(base+"utils/GLOBALutils/")

baryc_dir= base+'utils/SSEphem/'
sys.path.append(baryc_dir)
ephemeris='DEc403'
 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ceres modules
import pfsutils
import correlation
import GLOBALutils

# other useful modules
import argparse
import ephem
import jplephem
from math import radians as rad
from astropy.io import fits as pyfits
import pickle
import os
import scipy
import scipy.interpolate
from scipy import interpolate

import statsmodels.api as sm
lowess = sm.nonparametric.lowess

parser = argparse.ArgumentParser()
parser.add_argument('directorio')
parser.add_argument('-avoid_plot', action="store_true", default=False)
parser.add_argument('-dirout',default='default')
parser.add_argument('-do_class', action="store_true", default=False)
parser.add_argument('-just_extract', action="store_true", default=False)
parser.add_argument('-npools', default=1)
parser.add_argument('-o2do',default='all')
parser.add_argument('-reffile',default='default')
parser.add_argument('-ofind', default='last')
parser.add_argument('-binning', default='2x2')

args = parser.parse_args()
dirin            = args.directorio
avoid_plot       = args.avoid_plot
dirout           = args.dirout
DoClass          = args.do_class
JustExtract      = args.just_extract
npools           = int(args.npools)
object2do        = args.o2do
reffile          = args.reffile
stst             = args.ofind
binning          = args.binning 

if dirin[-1] != '/':
    dirin = dirin + '/'

if dirout == 'default':
    dirout = dirin[:-1]+'_red'+binning+'/'

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
use_median_shift = 1
force_pre_process  = False
force_flat     = False
force_flat_extract = False
force_thar_extract = False
force_thar_wavcal  = False
force_tharxc       = False
force_P_det        = False
force_P            = False
force_sci_extract  = False
force_spectral_file_build = True
force_stellar_pars = False
dumpargon          = False

minlines_glob      = 500
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
min_extract_col    = 0
max_extract_col    = 2032
bias_section       = [1,2032,2160]
over_section       = [0,2032,2160]
lsnr1              = 1000
lsnr2              = 1101

fact = 1.
if binning == '1x1':
    ext_aperture *= 2
    min_extract_col    *= 2
    max_extract_col    *= 2
    bias_section[1]    = bias_section[1] * 2
    bias_section[2]    = bias_section[2] * 2
    over_section[1]    = over_section[1] * 2
    over_section[2]    = over_section[2] * 2
    lsnr1              = 1900
    lsnr2              = 2101
    fact = 2.
    

ncoef_x   = 6
ncoef_m   = 7
npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2

models_path = base+"data/COELHO_MODELS/R_40000b/"
order_dir   = "wavcals/"

n_useful = 66    # up to which order do we care?
#############################

log = dirout+'night.log'

print "\n\n\tPFS MAGELLAN6.5m  PIPELINE\n"
print "\tRAW data is in ",dirin
print "\tProducts of reduction will be in",dirout
print "\tWill use only images with binning of", binning
print '\n'

biases, quartz, iodines, science, thars, thar_dates, obnames, exptimes = pfsutils.FileClassify(dirin,log,binning)
nightlog = open(log,'r')
loglines = nightlog.readlines()
print "\tThese are all the images to proccess:"
for line in loglines:
    print '\t'+line[:-1]
if stst == 'last':
    if os.access(dirout+'findstar.txt',os.F_OK):
        fst = open(dirout+'findstar.txt','r')
        stst = fst.readline()
        fst.close()
    else:
        raise ValueError("There is not a previously defined standard star file!!! You have to enter one (i.e. -ofind pfs0001.fits).\n")
else:
    fst = open(dirout+'findstar.txt','w')
    fst.write(stst+'\n')
    fst.close()

if (     (os.access(dirout+'Flat.fits',       os.F_OK) == False)    or \
         (os.access(dirout+'MasterBias.fits', os.F_OK) == False)    or \
         (os.access(dirout+'trace.pkl',       os.F_OK) == False)    or \
         (force_pre_process) ):
    print "\tNo previous pre-processing files or found"
    pre_process = 1
else:
    print "\tPre-processing files found, going straight to extraction"
    pre_process = 0

if (pre_process == 1):
    # median combine Biases
    print "\t\tGenerating Master calibration frames..."

    if len(biases)!=0:
        MasterBias, RO_bias, GA_bias = pfsutils.MedianCombine(biases, bias_section, over_section)
    else:
        MasterBias, RO_bias, GA_bias = np.zeros((max_extract_col,max_extract_col)),0,1
        print "\t\tWarning: 0 biases"
    hdu = pyfits.PrimaryHDU( MasterBias )
    if (os.access(dirout+'MasterBias.fits',os.F_OK)):
        os.remove(dirout+'MasterBias.fits')
    hdu.writeto(dirout+'MasterBias.fits')
    print "\t\t-> Masterbias: done!"
    #MasterBias = pyfits.getdata(dirout+'MasterBias.fits')
    #RO_bias, GA_bias = 1.,1.

    # median combine list of ob flats
    Flat, RO_flat, GA_flat = pfsutils.MedianCombine(quartz, bias_section, over_section, bias = MasterBias)
    # save this file for later reference
    hdu = pyfits.PrimaryHDU( Flat )
    if (os.access(dirout+'Flat.fits',os.F_OK)):
        os.remove(dirout+'Flat.fits')
    hdu.writeto(dirout+'Flat.fits')
    print "\t\t-> Masterflat: done!"
    #Flat = pyfits.getdata(dirout+'Flat.fits')
    #RO_flat, GA_flat = 1.,1.
    
    # Find orders & traces
    print "\tTracing echelle orders..."
    h = pyfits.open(dirin+stst)[0]
    hth = pyfits.getheader(dirin+stst)
    d = h.data
    d = pfsutils.OverscanTrim(d,bias_section, over_section)
    d -= MasterBias
    c_all, nord = GLOBALutils.get_them(d,ext_aperture-1,trace_degree, mode=1, nsigmas=5)

    trace_dict = {'c_all':c_all, 'nord':nord, 'GA_bias': GA_bias, 'RO_bias': RO_bias, 'GA_flat': GA_flat, 'RO_flat': RO_flat}
    pickle.dump( trace_dict, open( dirout+'trace.pkl', 'w' ) )


else:
    print '\t\tLoading Masterbias, Masterflat and traces'
    trace_dict = pickle.load( open( dirout+'trace.pkl', 'r' ) )
    c_all = trace_dict['c_all']
    nord = trace_dict['nord']
    GA_bias = trace_dict['GA_bias']
    RO_bias = trace_dict['RO_bias']
    GA_flat = trace_dict['GA_flat']
    RO_flat = trace_dict['RO_flat']
    # recover flats & master bias
    h = pyfits.open(dirout+'Flat.fits')
    Flat = h[0].data
    h = pyfits.open(dirout+'MasterBias.fits')
    MasterBias = h[0].data
print '\n'

# Make flat
NorFlat_fits = dirout + 'NorFlat.fits'
force_flat = False
if ( os.access(NorFlat_fits,os.F_OK) == False ) or (force_flat):
    print "\tNo previous Normalized Flat computed or determination forced, computing and saving..."
    NorFlat = np.ones(Flat.shape)
    eje = np.arange(Flat.shape[1])
    medianFlat = scipy.signal.medfilt(Flat,[5,9])
    #NorFlat = Flat/medianFlat
    
    vals = []
    ffact = 1.
    if binning == '1x1':
        ffact = 2.
    for order in range(nord):
        y = np.polyval(c_all[order],eje)
        maxval = 0.
        pointsy,pointsx = [],[]
        for xi in range(Flat.shape[1]):
            yi = np.around(y[xi] - 8*ffact).astype('int')
            #print xi,yi
            while yi < Flat.shape[0] and yi < np.around(y[xi] + 8*ffact +1).astype('int') and medianFlat[yi,xi]!=0.:
                pointsy.append(yi)
                pointsx.append(xi)
                if medianFlat[yi,xi] >= maxval:
                    maxval = medianFlat[yi,xi]
                yi+=1
        mat = (np.array(pointsy),np.array(pointsx))
        vals.append(maxval)
        NorFlat[mat] = Flat[mat]/medianFlat[mat]
    I = np.where(np.isnan(NorFlat)==True)
    NorFlat[I] = 1.
    I = np.where(NorFlat==0)
    NorFlat[I] = 1.
    if (os.access(NorFlat_fits,os.F_OK)):
        os.remove( NorFlat_fits )    
    hdu = pyfits.PrimaryHDU( NorFlat )
    hdu.writeto( NorFlat_fits )
else:
    print "\tNormalized Flat found, loading..."
    NorFlat = pyfits.getdata( NorFlat_fits )    

# Compute profile for optimal extraction
P_fits = dirout + 'P.fits'

if ( os.access(P_fits,os.F_OK) == False ) or (force_P_det):
    print "\n\tNo previous P matrix determination or determination forced, computing and saving..."
    
    h = pyfits.open(dirin+stst)[0]
    hth = pyfits.getheader(dirin+stst)
    d = h.data
    d = pfsutils.OverscanTrim(d,bias_section, over_section)
    d -= MasterBias
    d /= NorFlat
    Centers = np.zeros((len(c_all),d.shape[1]))
    for i in range(nord):
        Centers[i,:]=scipy.polyval(c_all[i,:],np.arange(len(Centers[i,:])))

    force_bkg = True
    bkg_obj_fits = dirout + 'BKG_' + h.header['UT-DATE'] + '_' + h.header['UT-TIME'] +'.'+ h.header['OBJECT'] +'.fits'
    if ( os.access(bkg_obj_fits,os.F_OK) == False or force_bkg):
        bkg = GLOBALutils.get_scat(d,Centers,span=ext_aperture)
        if (os.access(bkg_obj_fits,os.F_OK)):
            os.remove( bkg_obj_fits )
        hdu = pyfits.PrimaryHDU( bkg )
        hdu.writeto( bkg_obj_fits )
    else:
        bkg = pyfits.getdata(bkg_obj_fits)

    d -= bkg

    RO,GA = hth['ENOISE'],hth['EGAIN']

    P = GLOBALutils.obtain_P(d,c_all,ext_aperture,RO,\
                                    GA,NSigma_Marsh, S_Marsh, \
                    N_Marsh, Marsh_alg, min_extract_col,\
                    max_extract_col, npools)

    # write P as fits file
    if (os.access(P_fits,os.F_OK)):
        os.remove( P_fits )    
    hdu = pyfits.PrimaryHDU( P )
    hdu.writeto( P_fits )
   
else:
    print "\n\tP matrix found, loading..."
    P      = pyfits.getdata( P_fits )


# Extract Flat
print '\n\tExtraction of Flat calibration frames:'
Flat_spec_fits = dirout + 'Flat_spec.fits'
if ( os.access(Flat_spec_fits,os.F_OK) == False ) or (force_flat_extract):
    print "\t\tNo previous Flat extracted or extraction forced, extracting and saving..."
    flat_S  = np.zeros( (nord,3,Flat.shape[1]) )
    NFlat = Flat/NorFlat
    flat_S = GLOBALutils.simple_extraction(NFlat,c_all,ext_aperture,min_extract_col,max_extract_col,npools)
    #for i in range(nord):
        #    flat_S[i,:,:] = getSpectrum( P, NFlat, c_all[i,:], ext_aperture, RO_flat, GA_flat, S_Marsh, 0*NCosmic_Marsh, min_extract_col,max_extract_col )
    
    #flat_S = PFSutils.FlatNormalize( flat_S, mid = int(np.round(.5*flat_S.shape[2]))-1 )
    if (os.access(Flat_spec_fits,os.F_OK)):
        os.remove( Flat_spec_fits )
    hdu = pyfits.PrimaryHDU( flat_S )
    hdu.writeto( Flat_spec_fits )
else: 
    print "\t\tExtracted flat found, loading..."
    flat_S = pyfits.getdata( Flat_spec_fits )

print '\n\tExtraction of ThAr calibration frames:'
# Extract all ThAr files
for fsim in thars:
    hthar = pyfits.open( fsim )[0]
    dthar = pfsutils.OverscanTrim( hthar.data, bias_section, over_section ) - MasterBias
    dthar /= NorFlat
    hd = pyfits.getheader(fsim)
    thar_fits        = dirout + 'PFS_' + hthar.header['UT-DATE'] + '_' + hthar.header['UT-TIME'] +'.ThAr.spec.fits'
    thar_fits_simple = dirout + 'PFS_' + hthar.header['UT-DATE'] + '_' + hthar.header['UT-TIME'] +'.ThAr.spec.simple.fits'

    if ( os.access(thar_fits,os.F_OK) == False ) or ( os.access(thar_fits_simple,os.F_OK) == False ) or (force_thar_extract):
        print "\t\tNo previous extraction or extraction forced for ThAr file", fsim, "extracting..."

        thar_Ss = GLOBALutils.simple_extraction(dthar,c_all,ext_aperture,min_extract_col,max_extract_col,npools)
        thar_S  = GLOBALutils.optimal_extraction(dthar,P,c_all,ext_aperture,hthar.header['ENOISE'],\
          hthar.header['EGAIN'],S_Marsh,NCosmic_Marsh,min_extract_col,max_extract_col,npools)
            
        # save as fits file
        if (os.access(thar_fits,os.F_OK)):
            os.remove( thar_fits )
        if (os.access(thar_fits_simple,os.F_OK)):
            os.remove( thar_fits_simple )
            
        hdu = pyfits.PrimaryHDU( thar_S )
        #hdu.header = hd
        hdu = GLOBALutils.update_header(hdu,'SLITNAME', hd['SLITNAME'])
        hdu.writeto( thar_fits )

        hdu = pyfits.PrimaryHDU( thar_Ss )
        #hdu.header = hd
        hdu = GLOBALutils.update_header(hdu,'SLITNAME', hd['SLITNAME'])
        hdu.writeto( thar_fits_simple )


print "\n\tWavelength solution of ThAr calibration spectra:"
# Compute wavelength calibration of ThAr
sorted_thar_dates = np.argsort( thar_dates )
p0_array = np.zeros( (5, npar_wsol) )

for i in range(len(thars)):
    index = sorted_thar_dates[i]
    hthar = pyfits.open( thars[index] )
    wavsol_pkl = dirout + 'PFS_' + hthar[0].header['UT-DATE'] + '_' + hthar[0].header['UT-TIME']+'.wavsolpars.pkl'
    if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
        print "\t\tWorking on initial ThAr file", thars[index] 
        
        mjd, mjd0 = pfsutils.mjd_fromheader( hthar )
        thar_fits = dirout + 'PFS_' + hthar[0].header['UT-DATE'] + '_' + hthar[0].header['UT-TIME'] +'.ThAr.spec.fits'
        thar_S = pyfits.getdata( thar_fits )
        hd = pyfits.getheader(thar_fits)
        thar_out = dirout + 'PFS_' + hthar[0].header['UT-DATE'] + '_' + hthar[0].header['UT-TIME'] +'.ThAr.wav.fits'

        thar_data = np.zeros((2,thar_S.shape[0],thar_S.shape[2]))
        thar_data[1]= thar_S[:,1,:]

        lines_thar  = thar_S[:,1,:]
        iv_thar     = thar_S[:,2,:]

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
            
            thar_order_orig = lines_thar[order,:]
            IV              = iv_thar[order,:]
            wei             = np.sqrt( IV )
            #bkg             = PFSutils.Lines_mBack(thar_order_orig, IV,  thres_rel=3)        
            #thar_order      = thar_order_orig - bkg
            thar_order      = thar_order_orig
        
            coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths,\
            rms_ms, residuals, centroids, sigmas, intensities =\
                GLOBALutils.Initial_Wav_Calibration(order_dir+'norder_'+\
                order_s+'o.iwdat', thar_order, order, wei, rmsmax=100, \
                minlines=10,FixEnds=False,Dump_Argon=dumpargon,\
                Dump_AllLines=True, Cheby=use_cheby,porder=ncoef_x,fact=fact,del_width=30.)

            if (order == 32): 
                if (use_cheby):
                    Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, 0.5*len(thar_order), len(thar_order))
                else:
                    Global_ZP = scipy.polyval( coeffs_pix2wav, 0.0 )

            All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
            All_Wavelengths   = np.append( All_Wavelengths, wavelengths )
            All_Orders        = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order )
            All_Centroids     = np.append( All_Centroids, centroids)
            All_Sigmas        = np.append( All_Sigmas, sigmas)
            All_Intensities   = np.append( All_Intensities, intensities )

    #show()
        p0 = np.zeros( npar_wsol )
        p0[0] =  (32+93) * Global_ZP 
        p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = \
            GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths, All_Orders,\
                                                np.ones(All_Intensities.shape), p0, Cheby=use_cheby,\
                                                maxrms=50, Inv=Inverse_m,minlines=minlines_glob,order0=93, \
                                                ntotal=n_useful,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m)

        equis = np.arange( thar_S.shape[2] )    
        for order in range(n_useful):
            m = order + 93
            chebs = GLOBALutils.Calculate_chebs(equis, m, npix=thar_S.shape[2], order0=93, ntotal=n_useful, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
            WavSol = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(p1,chebs,ncoef_x,ncoef_m)   
            thar_data[0,order,:] = WavSol

        if (os.access(thar_out,os.F_OK)):
            os.remove( thar_out )
            
        hdu = pyfits.PrimaryHDU( thar_data )
        hdu = GLOBALutils.update_header(hdu,'SLITNAME', hd['SLITNAME'])
        hdu.writeto( thar_out )

        #for orde in np.unique(G_ord):
        #   I = np.where(G_ord == orde)[0]
        #   plot(G_wav[I],G_res[I],'o')
        #show()
    
        pdict = {'p1':p1,'mjd':mjd, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms,\
                     'G_res':G_res, 'All_Centroids':All_Centroids, 'All_Orders':All_Orders, 'All_Sigmas':All_Sigmas}
        pickle.dump( pdict, open( wavsol_pkl, 'w' ) )
        #print "\tMedian sigma:", np.median( All_Sigmas )

    else:
        print "\t\tUsing previously computed wavelength solution in file",wavsol_pkl
        pdict           = pickle.load(open(wavsol_pkl,'r'))


xc_fout_f = dirout+'ThAr_XCor+DeltaL.dat'
force_tharxc = True
if ( (os.access(xc_fout_f,os.F_OK) == False) or (force_tharxc)):  
    xc_fout = open(xc_fout_f,'w')
    for i in range(len(sorted_thar_dates)):
        fsim = thars[sorted_thar_dates[i]]
        hthar = pyfits.open( fsim )

        thar_fits = dirout + 'PFS_' + hthar[0].header['UT-DATE'] + '_' + hthar[0].header['UT-TIME'] +'.ThAr.spec.fits'
        thar_S = pyfits.getdata( thar_fits )

        mjd, mjd0 = pfsutils.mjd_fromheader( hthar )
        for ii in range(len(sorted_thar_dates)):
            fsim2 = thars[sorted_thar_dates[ii]]
            hthar2 = pyfits.open( fsim2 )

            thar_fits2 = dirout + 'PFS_' + hthar2[0].header['UT-DATE'] + '_' + hthar2[0].header['UT-TIME'] +'.ThAr.spec.fits'
            thar_S2 = pyfits.getdata( thar_fits2 )

            mjd2, mjd02 = pfsutils.mjd_fromheader( hthar2 )
            if (thar_fits != thar_fits2) and (mjd2 > mjd):
                # extract products

                lines_thar  = thar_S[:,1,:]
                lines_thar2 = thar_S2[:,1,:]
                
                pdict1 = pickle.load( open(dirout + 'PFS_' + hthar[0].header['UT-DATE'] + '_' + hthar[0].header['UT-TIME']+'.wavsolpars.pkl','r' ) )
                pdict2 = pickle.load( open(dirout + 'PFS_' + hthar2[0].header['UT-DATE'] + '_' + hthar2[0].header['UT-TIME']+'.wavsolpars.pkl','r' ) )
                
                # cross-correlate
                nob=n_useful
                out_shifts = np.zeros( nob )
                out_dl     = np.zeros( nob )
                for j in range(nob):
                    XCor, shifts = GLOBALutils.XC_ThAr(lines_thar[j,:], lines_thar2[j,:], 20)
                    p1    = GLOBALutils.XC_Gau_Fit(shifts, XCor, back_lag=2, usemin=0)
                    out_shifts[j] = p1[1] 
                    # Now measure Delta \lambda
                    m = j + 93
                    equis = np.arange(50,1950)
                    chebs = GLOBALutils.Calculate_chebs(equis, m, npix=lines_thar.shape[1], order0=93, ntotal=n_useful, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
                    WavSol1 = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(pdict1['p1'],chebs,ncoef_x,ncoef_m)    
                    WavSol2 = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(pdict2['p1'],chebs,ncoef_x,ncoef_m) 
                    out_dl[j]  = np.mean( WavSol1 - WavSol2 ) * 299792458.0 / WavSol1[1000]
                    rep = scipy.interpolate.splrep(equis,WavSol1)
                    dldx  = scipy.interpolate.splev(1000.0,rep,der=1)
                    out_shifts[j] = p1[1] * dldx * 299792458.0 / WavSol1[1000]

                uno = np.median( out_shifts )
                dos = np.mean( out_dl ) 

                # write out products
                line_out = "%20.12f %20.12f %12.6f %14.8f %s %s\n" % (mjd, mjd2, uno, dos,thar_fits.split('/')[-1],thar_fits2.split('/')[-1])
                xc_fout.write(line_out)
                #xc_fout.flush()
    xc_fout.close()

print '\n\tThe following targets will be processed:'

### start of science frame reductions ###
new_list = []
new_list_obnames = []
new_list_texp = []
for i in range(len(science)):
    fsim   = science[i]
    obname = obnames[i]
    texp   = exptimes[i]
    if (object2do == 'all'):
        new_list.append(fsim)
        new_list_obnames.append( obname )
        new_list_texp.append( texp )
        print '\t'+obname
    else:
        if (obname == object2do):
            new_list.append(fsim)
            new_list_obnames.append( obname )
            new_list_texp.append( texp )
        print '\t'+obname

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

#new_list = new_list[1:]           
for nlisti in range(len(new_list)):
    fsim = new_list[ nlisti ]
    obname = new_list_obnames[ nlisti ]
    TEXP =  new_list_texp[ nlisti ]

    print '\n'
    print "\t--> Working on image: ", fsim
    print "\t\tObject name:",obname
    know_moon = False
    if fsim.split('/')[-1] in spec_moon:
        I = np.where(fsim.split('/')[-1] == spec_moon)[0]
        know_moon = True
        here_moon = use_moon[I]

    h = pyfits.open(fsim)
    # get mjd and mjd0
    mjd,mjd0 = pfsutils.mjd_fromheader(h)

    #  get gain and readnoise of object 
    ronoise = h[0].header['ENOISE']
    gain    = h[0].header['EGAIN']
    iodine  = h[0].header['IODINE']
    
    print "\t\tWorking on an observation of object",obname

    # Find lambda_bary/lambda_topo using baryc
    altitude    = h[0].header['SITEALT']
    latitude    = h[0].header['SITELAT']
    longitude   = h[0].header['SITELONG']
    ra          = h[0].header['RA-D']
    dec         = h[0].header['DEC-D']
    epoch       = h[0].header['EQUINOX']

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

    res = jplephem.doppler_fraction(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)
    lbary_ltopo = 1.0 + res['frac'][0]
    bcvel_baryc = ( lbary_ltopo - 1.0 ) * 2.99792458E5

    print "\t\tBarycentric velocity:", bcvel_baryc

    res = jplephem.pulse_delay(ra/15.0, dec, int(mjd), mjd%1, 1, 0.0)  
    mbjd = mjd + res['delay'][0] / (3600.0 * 24.0)

    # Moon Phase Calculations
    gobs = ephem.Observer()  
    gobs.name='Clay_Mag_2'  
    gobs.lat=rad(latitude)  # lat/long in decimal degrees  
    gobs.long=rad(longitude)
    DDATE = h[0].header['UT-DATE']
    HHOUR = h[0].header['UT-TIME']
    Mho = HHOUR[:2]
    Mmi = HHOUR[3:5]
    Mse = HHOUR[6:]
    gobs.date = str(DDATE[:4]) + '-' +  str(DDATE[5:6]) + '-' + str(DDATE[7:]) + ' ' +  Mho + ':' + Mmi +':' +Mse
    mephem = ephem.Moon()
    mephem.compute(gobs)
    Mcoo = jplephem.object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Mp = jplephem.barycentric_object_track("Moon", int(mjd), float(mjd%1), 1, 0.0)
    Sp = jplephem.barycentric_object_track("Sun", int(mjd), float(mjd%1), 1, 0.0)
    res = jplephem.object_doppler("Moon", int(mjd), mjd%1, 1, 0.0)
    lunation,moon_state,moonsep,moonvel = GLOBALutils.get_lunar_props(ephem,gobs,Mcoo,Mp,Sp,res,ra,dec)
    refvel = bcvel_baryc + moonvel
    print '\t\tRadial Velocity of sacttered moonlight:',refvel

    sorted_indices = np.argsort( np.abs( np.array(thar_dates) - mjd ) )

    # optimally and simply extract spectra
    sci_fits = dirout + 'PFS_' + h[0].header['UT-DATE'] + '_' + h[0].header['UT-TIME'] +'.'+ obname +'.spec.fits'
    sci_fits_simple = dirout + 'PFS_' + h[0].header['UT-DATE'] + '_' + h[0].header['UT-TIME'] +'.'+ obname +'.spec.simple.fits'
    P_fits = dirout + 'P_' + h[0].header['UT-DATE'] + '_' + h[0].header['UT-TIME'] +'.'+ obname +'.fits'


    # Open file, trim, overscan subtract and MasterBias subtract
    data = h[0].data
    data = pfsutils.OverscanTrim(data,bias_section, over_section)
    force_sci_extract = True
    if ( os.access(sci_fits,os.F_OK) == False ) or ( os.access(sci_fits_simple,os.F_OK) == False ) or (force_sci_extract):

        data -= MasterBias
        data /= NorFlat
        
        print '\t\t\tRecentering traces...'
        c_alls, pshift = GLOBALutils.retrace( data, c_all )

        Centers = np.zeros((len(c_alls),data.shape[1]))
        for i in range(nord):
            Centers[i,:]=scipy.polyval(c_alls[i,:],np.arange(len(Centers[i,:])))

        force_bkg = False
        bkg_obj_fits = dirout + 'BKG_' + h[0].header['UT-DATE'] + '_' + h[0].header['UT-TIME'] +'.'+ obname +'.fits'
        if ( os.access(bkg_obj_fits,os.F_OK) == False or force_bkg):
            bkg = GLOBALutils.get_scat(data,Centers,span=ext_aperture)
            if (os.access(bkg_obj_fits,os.F_OK)):
                os.remove( bkg_obj_fits )
            hdu = pyfits.PrimaryHDU( bkg )
            hdu.writeto( bkg_obj_fits )
        else:
            bkg = pyfits.getdata(bkg_obj_fits)

        data -= bkg

        if os.access(P_fits,os.F_OK) == False or force_P:
            P = GLOBALutils.obtain_P(data,c_alls,ext_aperture,ronoise,\
                                    gain,NSigma_Marsh, S_Marsh, \
                    N_Marsh, Marsh_alg, min_extract_col,\
                    max_extract_col, npools)

            if (os.access(P_fits,os.F_OK)):
                os.remove( P_fits )
            hdu = pyfits.PrimaryHDU( P )
            hdu.writeto( P_fits )
        else:
            P = pyfits.getdata(P_fits)

        if ( os.access(sci_fits,os.F_OK) == False ) or ( os.access(sci_fits_simple,os.F_OK) == False ) or (force_sci_extract):

            print "\t\tNo previous extraction or extraction forced for science file", fsim, "extracting..."
            sci_Ss = GLOBALutils.simple_extraction(data,c_alls,ext_aperture,min_extract_col,max_extract_col,npools)
            sci_S  = GLOBALutils.optimal_extraction(data,P,c_alls,ext_aperture,ronoise,gain,S_Marsh,NCosmic_Marsh,\
         min_extract_col,max_extract_col,npools)
            
        # save as fits file
            if (os.access(sci_fits,os.F_OK)):
                os.remove( sci_fits )
            if (os.access(sci_fits_simple,os.F_OK)):
                os.remove( sci_fits_simple )
            
            hdu = pyfits.PrimaryHDU( sci_S )
            hdu.writeto( sci_fits )
            hdu = pyfits.PrimaryHDU( sci_Ss )
            hdu.writeto( sci_fits_simple )

        else:
            print '\t\t'+fsim, "has already been extracted, reading in product fits files..."
            sci_S = pyfits.getdata( sci_fits )
            sci_Ss = pyfits.getdata( sci_fits_simple )
    
    
    fout = 'proc/'+ obname + '_' + \
         h[0].header['UT-DATE'] + '_' + \
         h[0].header['UT-TIME'] + '_' + \
         iodine + '_' + 'sp.fits'

    #Build spectra
    
    if ( os.access(dirout+fout ,os.F_OK) == False ) or (force_spectral_file_build):
        # initialize file that will have the spectra
        # n_useful should be nord_ob, but we still have not calibrated that bluest order -- TODO
        spec = np.zeros((11, n_useful, data.shape[1]))
        hdu = pyfits.PrimaryHDU( spec )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MJD', mjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MBJD', mbjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START DATE', h[0].header['UT-DATE'] )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START UT',  h[0].header['UT-TIME'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TEXP (S)',h[0].header['EXPTIME'])
        #hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (km/s) [OBSOLETE]', bcvel)
        #hdu = GLOBALutils.update_header(hdu,'HIERARCH EARTH ROTATION CORRECTION (km/s) [OBSOLETE]', gcvel)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (KM/S)', bcvel_baryc)
        #hdu = GLOBALutils.update_header(hdu,'HIERARCH EARTH ROTATION CORRECTION (km/s)', vel_rot)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH (LAMBDA_BARY / LAMBDA_BARY)', lbary_ltopo)    
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARGET NAME', obname)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA',h[0].header['RA-D'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC',h[0].header['DEC-D'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH RA BARY',ra)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC BARY',dec)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH EQUINOX',h[0].header['EQUINOX'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LATITUDE',h[0].header['SITELAT'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LONGITUDE',h[0].header['SITELONG'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS ALTITUDE',h[0].header['SITEALT'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARG AIRMASS',h[0].header['AIRMASS'])

        # get ThAr closest in time
        indice = sorted_indices[0]
        hthar = pyfits.open(thars[indice])
        thar_fits_ob = dirout + 'PFS_' + hthar[0].header['UT-DATE'] + '_' + hthar[0].header['UT-TIME'] +'.ThAr.spec.fits'
        pkl_wsol = dirout + 'PFS_' + hthar[0].header['UT-DATE'] + '_' + hthar[0].header['UT-TIME'] +'.wavsolpars.pkl'
        print "\t\tUnpickling wavelength solution from", pkl_wsol, " ..."
        wsol_dict = pickle.load(open(pkl_wsol,'r'))

        # Apply new wavelength solution including barycentric correction
        equis = np.arange( data.shape[1] )        

        for order in range(n_useful):
            m = order + 93
            chebs = GLOBALutils.Calculate_chebs(equis, m, npix=data.shape[1], order0=93, ntotal=n_useful, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
            WavSol = lbary_ltopo * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wsol_dict['p1'],chebs,ncoef_x,ncoef_m)   
            spec[0,order,:] = GLOBALutils.ToVacuum(WavSol)
            spec[1,order,:] = sci_S[order,1, :]
            spec[2,order,:] = sci_S[order,2, :]
            spec[3,order,:] = spec[1,order,:] / flat_S[order]
            spec[4,order,:] = spec[2,order,:] * flat_S[order] ** 2
            #cont_coef = continuum.NORM_single(spec[0,order,:],spec[3,order,:],orden=3)
            cont_coef = GLOBALutils.get_cont_single(spec[0,order],spec[3,order],spec[4,order],ll=1.5,lu=5,nc=3)
            ratio = np.polyval(cont_coef, spec[0,order,:])
            L  = np.where( spec[1,order,:] != 0 )
            spec[5,order,:][L] = spec[3,order,:][L] / ratio[L]
            spec[6,order,:][L] = spec[4,order,:][L] * (ratio[L] ** 2 )
            spec[7,order,:][L] = ratio[L]
            spec[8,order,:][L] = ratio[L] * flat_S[order][L] / np.sqrt( ratio[L] * flat_S[order][L] / gain + (ronoise/gain)**2 )
            rI = np.where(spec[5,order] > 1. + 8./spec[8,order])
            spec[5,order,rI] = 1.
            #plot(spec[0,order],spec[5,order])

            spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), WavSol,k=3)
            dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
            NN            = np.average(dlambda_dx)
            dlambda_dx    /= NN

            spec[9,order,:][L] = spec[5,order,:][L] * (dlambda_dx[L] ** 1) 
            spec[10,order,:][L] = spec[6,order,:][L] / (dlambda_dx[L] ** 2)

    #show()
    #JustExtract = False
    if (not JustExtract):
        if h[0].header['SLITNAME'] == '0.5x2.5':
            ref_RES = 76000
        elif h[0].header['SLITNAME'] == '0.7x2.5':
            ref_RES = 54000
        elif h[0].header['SLITNAME'] == '0.3x2.5':
            ref_RES = 127000
	print ref_RES
        if DoClass:
            # spectral analysis
            print "\t\tSpectral Analysis..."
            query_success = False
            # First, query SIMBAD with the object name
            query_success,sp_type_query = GLOBALutils.simbad_query_obname(obname)
            # Now, query SIMBAD by coordinates if above not successful
            if (not query_success):
                query_success,sp_type_query = GLOBALutils.simbad_query_coords('12:00:00','00:00:00')
            print "\t\t\tSpectral type returned by SIMBAD query:",sp_type_query

            hdu = GLOBALutils.update_header(hdu,'HIERARCH SIMBAD SPTYP', sp_type_query)

            pars_file = dirout + 'PFS_' + hthar[0].header['UT-DATE'] + '_' + hthar[0].header['UT-TIME'] +'_stellar_pars.txt'

            if os.access(pars_file,os.F_OK) == False or force_stellar_pars:
                print "\t\t\tEstimating atmospheric parameters:"
                spec2 = spec.copy()
                Rx = np.around(1./np.sqrt(1./40000.**2 - 1./ref_RES**2))
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
        
        velw  = 300
        velsh = 3.
            
        # Read in mask
        ml, mh, weight = np.loadtxt(mask,unpack=True)
        ml_v = GLOBALutils.ToVacuum( ml )
        mh_v = GLOBALutils.ToVacuum( mh )
       
        # make mask larger accounting for factor ~2.5 lower res in FEROS w/r to HARPS
        av_m = 0.5*( ml_v + mh_v )
        ml_v -= 1.*(av_m - ml_v)
        mh_v += 1.*(mh_v - av_m)
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
            # first rough correlation to find the minimum
            vels, xc_full, sn, nlines_ccf, W_ccf = \
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight, 0, lbary_ltopo, vel_width=velw,vel_step=velsh,\
                                          spec_order=9,iv_order=10,sn_order=8,max_vel_rough=velw)

            xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3., Simple=True, W=W_ccf)
            yy     = scipy.signal.medfilt(xc_av,11)
            pred = lowess(yy, vels,frac=0.4,it=10,return_sorted=False)
            tck1 = scipy.interpolate.splrep(vels,pred,k=1)
            xc_av_orig = xc_av.copy()
            xc_av /= pred
            vel0_xc = vels[ np.argmin( xc_av ) ]
            rvels, rxc_av, rpred, rxc_av_orig, rvel0_xc = vels.copy(), xc_av.copy(), pred.copy(), xc_av_orig.copy(), vel0_xc
            xc_av_rough = xc_av
            vels_rough  = vels
                
            vel_width = np.maximum( 20.0, 6*disp )

            vels, xc_full, sn, nlines_ccf, W_ccf =\
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight, vel0_xc, lbary_ltopo, vel_width=vel_width,vel_step=0.2,\
                                          spec_order=9,iv_order=10,sn_order=8,max_vel_rough=velw)
        
            xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3., Simple=True, W=W_ccf)
            pred = scipy.interpolate.splev(vels,tck1)
            xc_av /= pred
        
            if sp_type == 'M5':
                moon_sig = 2.5
            elif sp_type == 'K5':
                moon_sig = 3.3
            else:
                moon_sig = 4.5

            p1,XCmodel,p1gau,XCmodelgau,Ls2 = GLOBALutils.XC_Final_Fit( vels, xc_av , sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = False)

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

            #print p1gau[1]
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

        BSerr = -999.00

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


        SNR_5130 = np.median(spec[8,26,lsnr1:lsnr2] )

        airmass  = h[0].header['AIRMASS']
        seeing   = -999
                    
        if sp_type=='G2':
            if T_eff < 6000:
                D = 0.32815
                C = 0.00453
            else:
                D = 0.32815
                C = 0.00453
        elif  sp_type == 'K5':
            D = 0.27404
            C = 0.00433
        else:
            D = 0.27404
            C = 0.00433

        BSerr = D / float(np.round(SNR_5130)) + C
        RVerr = 0.5

        RV     = np.around(p1gau_m[1],4)  
        BS     = np.around(SP,4)   
        RVerr2 = np.around(RVerr,4)
        BSerr  = np.around(BSerr,4)
        print '\t\t\tRV = '+str(RV)+' +- '+str(RVerr2)
        print '\t\t\tBS = '+str(BS)+' +- '+str(BSerr)

        bjd_out = 2400000.5 + mbjd
        T_eff_err = 100
        logg_err = 0.5
        Z_err = 0.5
        vsini_err = 2
        XC_min = np.abs(np.around(np.min(XCmodel),2))

        SNR_5130 = np.around(SNR_5130)
        SNR_5130_R = np.around(SNR_5130*np.sqrt(2.3))
        # write to output
        disp_epoch = np.around(p1gau_m[2],1)
        hdu = GLOBALutils.update_header(hdu,'RV', RV)
        hdu = GLOBALutils.update_header(hdu,'RV_E', RVerr2)
        hdu = GLOBALutils.update_header(hdu,'BS', BS)
        hdu = GLOBALutils.update_header(hdu,'BS_E', BSerr)
        hdu = GLOBALutils.update_header(hdu,'DISP', disp_epoch)
        hdu = GLOBALutils.update_header(hdu,'SNR', SNR_5130)
        hdu = GLOBALutils.update_header(hdu,'SNR_R', SNR_5130_R)
        hdu = GLOBALutils.update_header(hdu,'INST', 'PFS')
        hdu = GLOBALutils.update_header(hdu,'RESOL', ref_RES)
        hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')
        hdu = GLOBALutils.update_header(hdu,'XC_MIN', XC_min)
        hdu = GLOBALutils.update_header(hdu,'BJD_OUT', bjd_out)
        line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f     pfs   ceres   %8d %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
                      (obname, bjd_out, RV, RVerr2, BS, BSerr, ref_RES, T_eff_epoch, logg_epoch, Z_epoch, vsini_epoch, XC_min, disp_epoch,\
               TEXP, SNR_5130_R, ccf_pdf)
        f_res.write(line_out)
        if (os.access( dirout + fout,os.F_OK)):
            os.remove( dirout + fout)
        hdu.writeto( dirout + fout )

f_res.close()
