import matplotlib
matplotlib.use('Agg')
from pylab import *
import scipy
from scipy import interpolate
from scipy import signal
sys.path.append("../utils/GLOBALutils")
import GLOBALutils
import os 
import glob
from astropy.io import fits as pyfits
from scipy import ndimage
import emcee
import corner

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

def FileClassify(diri, log):

	# define output lists
	biases         = []
	darks		   = []
	flats          = []
	flats_co	   = []
	thar           = []
	thar_co		   = []
	sim_sci        = []
	dar_sci		   = []
	thart = []
	simt = []

	f = open(log,'w')

	#Do not consider the images specified in dir+badfiles.txt
	bad_files = []
	if os.access(diri+'bad_files.txt',os.F_OK):
		bf = open(diri+'bad_files.txt')
		linesbf = bf.readlines()
		for line in linesbf:
			bad_files.append(diri+line[:-1])
		bf.close()
    
	all_files = glob.glob(diri+"/*fits")
	for archivo in all_files:
		#print archivo
		dump = False
		for bf in bad_files:
			if archivo == bf:
				dump = True
				break
		if dump == False:
			h = pyfits.open(archivo)
			hd = pyfits.getheader(archivo)

			if 'bias' in archivo:
				biases.append(archivo)
			elif 'dark' in archivo:
				darks.append(archivo)
			elif 'flat' in archivo:
				if '1F' in archivo:
					flats_co.append(archivo)
				elif '2F' in archivo:
					flats.append(archivo)
			elif 'ThAr' in archivo:
				if '1F' in archivo:
					thar_co.append(archivo)
				elif '2F' in archivo:
					thar.append(archivo)
					thart.append(mjd_fromheader(h)[0])
			else:
				if '1F' in archivo:
					dar_sci.append(archivo)
				else:
					sim_sci.append(archivo)
					simt.append(mjd_fromheader(h)[0])

				#line = "%-15s %10s %10s %8.2f %4.2f %8s %11s %s\n" % (obname, ra, delta, texp, airmass, date, hour, archivo)
				#f.write(line)

	thart = np.array(thart)
	simt = np.array(simt)
	f.close()
	biases, darks, flats, flats_co, thar, thar_co, sim_sci, dar_sci = \
		np.array(biases), np.array(darks), np.array(flats), np.array(flats_co), np.array(thar), \
		np.array(thar_co), np.array(sim_sci), np.array(dar_sci)
	It = np.argsort(thart)
	thar = thar[It]
	It = np.argsort(simt)
	sim_sci = sim_sci[It]
	return biases, darks, flats, flats_co, thar, thar_co, sim_sci, dar_sci

def make_flatOB(MasterFlat, c_co,exap=5):
	nord_co = len(c_co)
	flat = MasterFlat.T
	img_out = flat.copy()
	Centers = np.zeros((len(c_co),flat.shape[1]))
	ejx = np.arange(flat.shape[1])
	for o in range(nord_co):
		Centers[o,:]=scipy.polyval(c_co[o],ejx)
	for x in range(flat.shape[1]):
		baso = np.min(flat[:,x])
		for o in range(nord_co):
			cen = np.around(Centers[o,x])
			try:
				bas = np.min(flat[int(np.around(cen-3*exap)):int(np.around(cen+3*exap)),x])
			except:
				bas = baso
			img_out[int(np.around(cen-exap)):int(np.around(cen+exap+1)),x] = bas
	return img_out

def get_flatOB(MasterFlat, MasterFlat_co):
	d1 = MasterFlat[1000]
	d2 = MasterFlat_co[1000]
	ref2 = np.arange(len(d2))
	tck = scipy.interpolate.splrep(ref2,d2,k=3)
	shfts = np.arange(-6.,6.,0.001)
	ccfs = []
	for shft in shfts:
		ref = ref2 + shft
		t2 = scipy.interpolate.splev(ref,tck)
		ccfs.append(np.sum(t2*d1))
	shft = shfts[np.argmax(ccfs)]

	new_img = MasterFlat_co.copy()
	for x in range(MasterFlat_co.shape[0]):
		d = MasterFlat_co[x]
		ref = np.arange(len(d))
		tck = scipy.interpolate.splrep(ref,d,k=1)
		new_img[x] = MasterFlat[x] - scipy.interpolate.splev(ref+shft,tck)
	return new_img

def get_data(path):
    d = pyfits.getdata(path)
    d = np.fliplr(d.T)
    return d

def good_orders(coef,nord,ny,nx,ext_aperture):
    Centers = np.zeros((nord,nx))
    ejx = np.arange(nx)
    bad_inx = []
    for i in range(nord):
        Centers[i,:]=scipy.polyval(coef[i],ejx)
	I = np.where(Centers[i,:]+ext_aperture>ny)[0]
	if len(I)>0:
	    bad_inx.append(i)
    bad_inx = np.array(bad_inx)
    im = np.min(bad_inx)
    return coef[:im], nord-len(bad_inx)

def clean_orders(c_all,data,exap=5):
    nords = len(c_all)
    medc = int(.5*data.shape[1])
    d = np.median(data[:,medc-exap:medc+exap+1],axis=1)
    #plot(d)
    #show()
    Centers = np.zeros((len(c_all),data.shape[1]))
    ejx = np.arange(data.shape[1])
    for i in range(len(c_all)):
        Centers[i,:]=scipy.polyval(c_all[i],ejx)
    ccens = Centers[:,medc]
    dist1 = ccens[-1]-ccens[-2]
    dist2 = ccens[-2]-ccens[-3]
  
    if dist1 < dist2:
	print 'uno'
	i = nords - 1
    else:
	print 'dos'
	i = nords - 2
    c_co, c_ob = [],[]
    while i > 0:
	print i
        if len(c_co) == 0:
            c_co = c_all[i]
            c_ob = c_all[i-1]
        else:
            c_co = np.vstack((c_all[i],c_co))
            c_ob = np.vstack((c_all[i-1],c_ob))
        i -= 2

	
    #for i in range(len(c_ob)):
    #    Centers[i,:]=scipy.polyval(c_ob[i],ejx)
    #	plot(ejx,Centers[i,:],'r')
    #for i in range(len(c_co)):
    #    Centers[i,:]=scipy.polyval(c_co[i],ejx)
    #	plot(ejx,Centers[i,:],'b')
    #show()
    return c_ob, c_co, len(c_ob), len(c_co)

def mjd_fromheader(h):
    """
    return modified Julian date from header
    """
    datetu = h[0].header['DATE-OBS'][:10]
    ut     = h[0].header['DATE-OBS'][11:]
    mjd0,mjd,i = GLOBALutils.iau_cal2jd(int(datetu[0:4]),int(datetu[5:7]),int(datetu[8:10]))
    ut = float(ut[:2])+ float(ut[3:5])/60. + float(ut[6:])/3600.
    mjd_start = mjd + ut/24.0

    #print 'Warning!!! adding 5 hrs to comute MJD due to problem in header! CHECK in future!!'
    #mjd_start += 5./24.


    secinday = 24*3600.0
    fraction = 0.5
    texp     = h[0].header['EXPTIME'] #sec

    #print 'Warning!!! assuming that the provided date is at the end of the exposure!!! Check in the future!!'
    #mjd = mjd_start - (fraction * texp) / secinday

    mjd = mjd_start + (fraction * texp) / secinday

    return mjd, mjd0

def scat_flat(d,c_ob,c_co):
    Centers_co = np.zeros((len(c_co),d.shape[1]))
    Centers_ob = np.zeros((len(c_ob),d.shape[1]))
    ejx = np.arange(d.shape[1])
    ejy = np.arange(d.shape[0])
    for i in range(len(c_co)):
        Centers_co[i,:]=np.around(scipy.polyval(c_co[i],ejx)).astype('int')
        Centers_ob[i,:]=np.around(scipy.polyval(c_ob[i],ejx)).astype('int')
    i = 0
    while i < d.shape[1]:
        line = d[:,i]
	refmins = []
	valmins = []
	for j in range(len(Centers_co)):
	    p1 = Centers_ob[j,i]
	    p2 = Centers_co[j,i]
	    minv = line[p1-10:p2+10]
	    refv = ejy[p1-10:p2+10]
	    im = np.argmin(minv)
	    refmins.append(refv[im])
	    valmins.append(minv[im])
	refmins,valmis = np.array(refmins),np.array(valmins)
	#plot(line)
	#plot(refmins,valmins,'ro')
	#show()
	#print gfd
	    
def get_scat(sc,ps,binning):
	binning = float(binning)
	#print sc.shape,ps.shape
	ejex = np.arange(sc.shape[1])
	bac = sc.copy()
	for col in range(sc.shape[0]):
		#col = 1000
		#plot(sc[col])
		#show()
		#print gfds
		J = np.where(ps[:,col]==0)[0]
		Z = np.zeros((sc.shape[1]))
		Z[J] = 1.
		#plot(sc[col])
		#plot(sc[col]*Z,linewidth=2.0)
		#ncol = scipy.signal.medfilt(sc[col],11)
		Z2 = np.append(Z[-1],Z[:-1])
		I = np.where((Z!=0) & (Z2==0))[0]
		J = np.where((Z2!=0) & (Z==0))[0]
		J = J[1:]
		I = I[:-1]
		points = []
		vals = []
		for i in range(len(I)):
			#plot(numpy.mean(ejex[I[i]:J[i]]),numpy.median(sc[2000,I[i]:J[i]]),'ro')
			points.append(np.mean(ejex[I[i]:J[i]]))
			vals.append(np.median(sc[col,I[i]:J[i]]))
		points,vals = np.array(points),np.array(vals)
		#pi = 0
		#npoints,nvals = [],[]
		#while pi < len(points)-1:
		#	if vals[pi]>vals[pi+1]:
		#		npoints.append(points[pi+1])
		#		nvals.append(vals[pi+1])
		#	else:
		#		npoints.append(points[pi])
		#		nvals.append(vals[pi])
		#	pi+=2
		#points,vals = np.array(npoints),np.array(nvals)
		
		#vals = scipy.signal.medfilt(vals,3)
		#plot(points,vals,'ro')
		#plot(points,vals,'r')
		#plot(points,scipy.signal.medfilt(vals,3),'k')
		#show()
		F = np.where(vals < 10000)[0]
		vals = vals[F]
		points = points[F]
		tck = interpolate.splrep(points,vals,k=1)
		scat = interpolate.splev(ejex,tck)
		
		scat[:I[0]] = scat[I[0]]
		#scat[J[-1]:] = 0.
		#plot(scat)
		bac[col] = scat
		#plot(ejex,scat,'r')
		#show()
	bacm = signal.medfilt2d(bac,[51,1])
	#plot(bac[:,1000])
	#plot(sc[1000])
	#plot(bacm[1000])
	#show()
	return bacm

def get_name(arch):
	obname = pyfits.getheader(arch)['OBJECT']
	if 'HD10700' in arch:
		obname = 'HD10700'
	return obname

	    
def Initial_Wav_Calibration(filename,spec,order,wei, porder=3, rmsmax=75, minlines=10, FixEnds=True, \
                            Dump_Argon=False, Dump_AllLines=False, Cheby=False, rough_shift = 0.0,del_width=5.0, \
			    binning=1,line_width=4, fact=1,do_xc=True,sigmai=0.7,pixelization=False,fibre_ap=2.3):

	f = open(filename).readlines()
        bad_indices = np.array([])
	pixel_centers = array([])
	wavelengths   = array([])
	sigmas        = array([])
	centroids     = array([])
	intensities   = array([])

	delta=0.
	N_l = 0
	out = []
	for line in f:
	    if line[0]!='#':
		w = line.split()
		nlines = int(w[0])
		pix = float(w[1])*fact/float(binning) + delta + rough_shift
		wav = float(w[2])
		if pix > 20 and pix < len(spec)-20:
			sigma = sigmai * fact / float(binning)
			X = np.arange(int(np.around(pix))-5,int(np.around(pix))+6)
			Y = spec[X]
			guess = np.array([pix-0.5*fibre_ap,np.max(Y),sigmai])
			p1 = FitFideosCompProf(X, Y, guess, fibre_ap)

			pixel_centers = np.append(pixel_centers,p1[0])
			sigmas        = np.append(sigmas,p1[2])
			wavelengths   = np.append(wavelengths,wav)
			intensities   = np.append(intensities,p1[1])
			centroids     = np.append(centroids, p1[0])
			N_l += 1

	#show()
	pixel_centers2 = np.around(pixel_centers).astype('int')
	I = np.where((pixel_centers2>=0) & (pixel_centers2<len(spec)))
	pixel_centers2 = pixel_centers2[I]
	#plot(pixel_centers2,spec[pixel_centers2],'go')
	#show()
	#print gfd

	#I = np.where((pixel_centers>0) & (pixel_centers<2048))[0]
	#plot(np.around(pixel_centers[I]).astype('int'),spec[np.around(pixel_centers[I]).astype('int')],'ro')
	#show()
	I1 = np.where(pixel_centers<50)[0]
	I2 = np.where(pixel_centers>len(spec)-50)[0]
	II = np.hstack((I1,I2))
	bad_indices = np.hstack((np.array(bad_indices),II))
	bad_indices = list(np.unique(bad_indices))
	# now, do the polynomial fit, rejecting some lines until RMS is below rmsmax
	I = range( N_l )

	for bi in bad_indices:
		I.remove( bi )
		N_l -= 1

	if (Cheby):
		coeffs_pix2wav   = Cheby_Fit(pixel_centers[I], wavelengths[I], porder,len(spec))   
		coeffs_pix2sigma = Cheby_Fit(pixel_centers[I], sigmas[I], porder,len(spec))
	else:
		coeffs_pix2wav   = scipy.polyfit(pixel_centers[I], wavelengths[I], porder)    
		coeffs_pix2sigma = scipy.polyfit(pixel_centers[I], sigmas[I], porder)

	rmsms, residuals = rms_ms(coeffs_pix2wav, pixel_centers[I], wavelengths[I], len(spec), Cheby=Cheby)

	if (FixEnds):
		minI = np.min( I ) + 1
		maxI = np.max( I ) - 1
	else:
		minI = np.min( I )
		maxI = np.max( I )
	#if order==26:
	#	    plot(pixel_centers[I],residuals,'ro')
	#	    plot([0,4096],[0,0])
	#plot(np.arange(4096),Cheby_eval(coeffs_pix2wav,np.arange(4096),len(spec)))
	#show()
	#print dfgh
	count = 0
	while ((N_l > minlines) and (rmsms > rmsmax)):
		rmsms, residuals = rms_ms(coeffs_pix2wav, pixel_centers[I], wavelengths[I], len(spec), Cheby=Cheby)
		index_worst = np.argmax( np.absolute(residuals) )
		I.pop( index_worst)
		N_l -= 1
		if (Cheby):
			coeffs_pix2wav   = Cheby_Fit(pixel_centers[I], wavelengths[I], porder,len(spec))   
			coeffs_pix2sigma = Cheby_Fit(pixel_centers[I], sigmas[I], porder,len(spec))
		else:
			coeffs_pix2wav   = scipy.polyfit(pixel_centers[I], wavelengths[I], porder)    
			coeffs_pix2sigma = scipy.polyfit(pixel_centers[I], sigmas[I], porder)
		count +=1

	rmsms, residuals = rms_ms(coeffs_pix2wav, pixel_centers[I], wavelengths[I], len(spec), Cheby=Cheby)    

	pci = np.around(pixel_centers).astype('int')
	#plot(spec)
	#plot(pci,spec[pci],'ro')
	#plot(pci[I],spec[pci[I]],'bo')
	#show()
	#plot(wavelengths[I],residuals,'ro')
	#show()
	#print "RMS is ", rmsms, "using ", N_l, " lines at indices ", I
	#plot(pixel_centers[I],wavelengths[I],'ro')
	#if order == 26:
	#	    plot(pixel_centers[I],residuals-0.1,'bo')
	#	    plot([0,4096],[-0.1,-0.1])
	#	    #plot(np.arange(4096),Cheby_eval(coeffs_pix2wav,np.arange(4096),len(spec)))
	#	    show()
	#print order, len(pixel_centers), len(I)
	return coeffs_pix2wav, coeffs_pix2sigma, pixel_centers[I], wavelengths[I], \
        rmsms, residuals, centroids[I], sigmas[I], intensities[I]

def rms_ms(coeffs_pix2wav, pixel_centers, wavelengths, npix, Cheby=False):
    " Returns rms deviation of best fit in m/s"

    if (Cheby):
        residuals = Cheby_eval(coeffs_pix2wav,pixel_centers,npix) - wavelengths
        central_wav = 0.5 * (Cheby_eval(coeffs_pix2wav,50.,npix) + Cheby_eval(coeffs_pix2wav,npix-50,npix))
    else:
        residuals = scipy.polyval(coeffs_pix2wav,pixel_centers) - wavelengths
        central_wav = 0.5 * (scipy.polyval(coeffs_pix2wav,50.) + scipy.polyval(coeffs_pix2wav,npix-50))

    rms_ms = np.sqrt( np.var( residuals ) ) * 299792458.0 / central_wav
    
    return rms_ms, residuals

def fit2d(sc,trace,filename,spec,order,wei, porder=3, rmsmax=75, minlines=10, FixEnds=True, \
                            Dump_Argon=False, Dump_AllLines=False, Cheby=False, rough_shift = 0.0,del_width=5.0, \
			    binning=1,line_width=4, fact=1,do_xc=True,sigmai=0.8,pixelization=False):

	f = open(filename).readlines()

	pixel_centers = array([])
	wavelengths   = array([])
	sigmas        = array([])
	centroids     = array([])
	intensities   = array([])

	delta=0.
	N_l = 0
	out = []
	for line in f:
	    if line[0]!='#':
		w = line.split()
		nlines = int(w[0])
		pix = float(w[1])*fact/float(binning) + delta + rough_shift
		wav = float(w[2])
		if pix > 20 and pix < len(spec)-20:
			#X = array(range(pix-line_width,pix+line_width+1))
			#Y = spec[pix-line_width:pix+line_width+1]
			sigma = sigmai * fact / float(binning)
			#p1 = FitFideosCompProf(X, Y, B, mu[0], sigma[0])
			#print trace
			y1 = 2064 - np.polyval(trace,2048.-pix)
			#print '#',2048-pix,y1
			
			X = np.arange(int(np.around(2048.-pix))-5,int(np.around(2048-pix))+5)
			Y = np.arange(int(np.around(y1))-5,int(np.around(y1))+6)

			F = sc.T
			F = F[X,:]
			F = F[:,Y]
			B1 = F[0]
			B2 = F[-1]
			B3 = F[:,0]
			B4 = F[:,-1]
			BAC = np.median(np.hstack((B1,B2,B3,B4)))
			F = F - BAC
			Ferr = np.sqrt(np.sqrt(np.absolute(F))**2 + 9.)
			#pars = np.array([1.8128e+03, 1.7167e+03, 7.500e+03, 2.3, 6.15])
			#model = FideosCompProf(pars,X,Y)
			#model2 = FideosCompProf(np.array([2048-pix, y1,F.max()]),X,Y)
			#for i in range(F.shape[1]):
			#	errorbar(X,F[:,i],yerr=Ferr[:,i],fmt='b')
			#	plot(X,model[:,i],'r')
			#	plot(Y,model2[i],'k')
			#show()

			#print X,Y,pix,y1
			#p1 = FitFideos2DProf(X, Y, F, 2048-pix, y1)
			guess = np.array([2048-pix, y1,F.max(),2.,10.])
			ndim = len(guess)
			nwalkers=100
			pos = [guess + np.random.randn(ndim)*np.array([2,2,500.,2.,4.]) for i in range(nwalkers)]
			sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(X,Y,F,Ferr))
			sampler.run_mcmc(pos, 500)
			samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
			#print samples.shape
			print np.median(samples,axis=0)
			if len(out) == 0:
				out = np.median(samples,axis=0)
			else:
				out = np.vstack((out,np.median(samples,axis=0)))
			#fig = corner.corner(samples, labels=["$Xpos$", "$Ypos$", "Amp", 'Fibwidth', 'Inst'])
			#show()
	return out

def rms_ms(coeffs_pix2wav, pixel_centers, wavelengths, npix, Cheby=False):
    " Returns rms deviation of best fit in m/s"

    if (Cheby):
        residuals = Cheby_eval(coeffs_pix2wav,pixel_centers,npix) - wavelengths
        central_wav = 0.5 * (Cheby_eval(coeffs_pix2wav,50.,npix) + Cheby_eval(coeffs_pix2wav,npix-50,npix))
    else:
        residuals = scipy.polyval(coeffs_pix2wav,pixel_centers) - wavelengths
        central_wav = 0.5 * (scipy.polyval(coeffs_pix2wav,50.) + scipy.polyval(coeffs_pix2wav,npix-50))

    rms_ms = np.sqrt( np.var( residuals ) ) * 299792458.0 / central_wav
    
    return rms_ms, residuals

def Cheby_Fit(x,y,order,npix):
    """
    Fits Chebyshev polynomials to y as a function of x

    """
    med = .5*npix
    #normalize x
    x_norm = (x-med) / med

    def fitfunc(p,chebs,order):
        ret_val = 0.0
        for i in range(0,order+1):
            ret_val += p[order-i]*chebs[i]
        return ret_val
    errfunc = lambda p,chebs,y,order: np.ravel( (fitfunc(p,chebs,order)-y) )
    
    def get_chebs(x,order):
        chebs = []
        for i in range(0,order+1):
            chebs.append( scipy.special.chebyt(i)(x) )
        return chebs

    p0 = np.zeros( order + 1 )
    p0[order] = np.mean( y )
    chebs = get_chebs( x_norm, order)
    p1, success = scipy.optimize.leastsq(errfunc, p0, args=(chebs, y, order))
    return p1

def Cheby_eval(p,x,npix):
    """
    evaluates Chebyshev polynomial fit at x given best-fit parameters p
    """
    med = .5*npix
    x_norm = (x-med) / med
    order = len(p) - 1
    ret_val = 0.0
    for i in range(order + 1):
        ret_val += p[order - i]*scipy.special.chebyt(i)(x_norm)
    
    return ret_val

def LineFit_SingleSigma(X, Y, B, mu, sigma, weight,pixelization=False):
    """
    This function fits a series of Gaussians simultaneously, given a set
    of input pixels, sigmas, and intensities
    
    Sigma is the same for all lines
    
    Returns (mu_i, sigm_i, int_i), i=1,.,n where n is the number of components
    """

    # get number of components to be fit
    n = len(mu)

    def fitfunc(p, x, n):

	if pixelization:
	    lxo = len(x)
	    xo = x.copy()
	    x = np.arange(x[0]-0.5,x[-1]+0.5,0.01)

        ret = np.zeros(len(x))
        for i in range(n):
            ret += ( p[i*2+1] * IntGaussian(x,p[i*2+2],p[0]) )

	if pixelization:
	    ret = ret.reshape((lxo,100))
	    ret = np.mean(ret,axis=1)


        return ret
    errfunc = lambda p,x,n,y,weight: np.ravel( (fitfunc(p,x,n)-y)*weight )

    # Build error array
    #err = sqrt((B+Y)/gain + readnoise**2)
    
    # Build input parameters list
    # length of inout parameter list: 3 * n
    p0 = np.zeros(2*n+1)
    p0[0] = sigma[0]
    for i in range(n):
        p0[i*2+1] = Y[int(round(mu[i])-X[0])]
        p0[i*2+2] = mu[i]

    # perform fit
    #plot(X,Y,'b')
    
    p1, success = scipy.optimize.leastsq(errfunc,p0, args=(X,n,Y-B, weight))
    #plot(X,fitfunc(p1,X,n),'r')
    # build output consistent with LineFit
    p_output = np.zeros(3*n)
    for i in range(n):
        p_output[i*3]   = p1[i*2 + 1]
        p_output[i*3+1] = p1[i*2 + 2]
        p_output[i*3+2] = p1[0]

    return p_output, success  

def TruncGaussian(x,y,x0,y0,dev):
	if x > x0:
		return 0
	elif x < x0 - 3*dev:
		return 0
	c1 = (x - x0)**2/(2.*dev*dev)
	c2 = (y - y0)**2/(2.*dev*dev)
	z = np.exp(-(c1+c2))
	return z

def HalfFiber(x,y,x0,y0,lim):
	if x > x0:
		return 0.
	if np.sqrt((x-x0)**2 + (y-y0)**2)<lim:
		return 1.
	else:
		return 0.

def HalfFiberPos(x,y,x0,y0,lim):
	if x < x0:
		return 0.
	if np.sqrt((x-x0)**2 + (y-y0)**2)<lim:
		return 1.
	else:
		return 0.

def FitFideosCompProf(X, Y, p0, fibre_ap):

    def AnaliticProfile1D(p,x,R):
        space = 100.
        xt    = np.arange(x[0],x[-1]+1,1./space)
        x0 = xt - p[0]
        J = np.where((R**2-x0**2>=0) & (x0>=0))[0]
        out = np.zeros(len(x0))
        out[J] = 2*np.sqrt(R**2-x0[J]**2)
        out = scipy.ndimage.filters.gaussian_filter(out, p[2]*space)
        out /= np.max(out)
        out = p[1]*out
        out = np.mean(out.reshape((int(len(out)/space),int(space))),axis=1)
        return out
    errfunc = lambda p,X,Y,fibre_ap: np.ravel( (AnaliticProfile1D(p,X,fibre_ap)-Y) )

    p1, success = scipy.optimize.leastsq(errfunc,p0, args=(X,Y,fibre_ap))
    #plot(X,Y)
    #plot(X,AnaliticProfile1D(p1,X,fibre_ap))
    #plot(X,AnaliticProfile1D(p0,X,fibre_ap))
    return p1

def FideosCompProf(p,X,Y):
	    Xt = np.arange(X[0]-0.5,X[-1]+0.6,0.1)
	    Yt = np.arange(Y[0]-0.5,Y[-1]+0.6,0.1)
	    mat = np.zeros((len(Xt),len(Yt)))
	    for x in np.arange(len(Xt)):
	        for y in  np.arange(len(Yt)):
		    mat[x,y] = HalfFiber(Xt[x],Yt[y],p[0],p[1],p[3])
	    fin2 = scipy.ndimage.filters.gaussian_filter(mat, p[4])
	    fin2 /= np.max(fin2)
	    fin2 = p[2]*fin2
	    
	    mat2 = np.zeros((len(X),len(Y)))
	    for x in np.arange(len(X)):
	        I = np.where((Xt>=X[x]-0.5)&(Xt<X[x]+0.4999))[0]
		temp = fin2[I]
		for y in  np.arange(len(Y)):
			J = np.where((Yt>=Y[y]-0.5)&(Yt<Y[y]+0.4999))[0]
			mat2[x,y] = np.mean(temp[:,J])
	    return mat2


def FitFideos2DProf(X, Y, F, x1, y1):
	def FideosCompProf(p,X,Y):
	    Xt = np.arange(X[0]-0.5,X[-1]+0.6,0.1)
	    Yt = np.arange(Y[0]-0.5,Y[-1]+0.6,0.1)
	    mat = np.zeros((len(Xt),len(Yt)))
	    for x in np.arange(len(Xt)):
	        for y in  np.arange(len(Yt)):
		    mat[x,y] = HalfFiber(Xt[x],Yt[y],p[0],p[1],2.)
	    fin2 = scipy.ndimage.filters.gaussian_filter(mat, 5.)
	    fin2 /= np.max(fin2)
	    fin2 = p[2]*fin2
	    
	    mat2 = np.zeros((len(X),len(Y)))
	    for x in np.arange(len(X)):
	        I = np.where((Xt>=X[x]-0.5)&(Xt<X[x]+0.4999))[0]
		temp = fin2[I]
		for y in  np.arange(len(Y)):
			J = np.where((Yt>=Y[y]-0.5)&(Yt<Y[y]+0.4999))[0]
			mat2[x,y] = np.sum(temp[:,J])
	    return mat2
	errfunc = lambda p,X,Y,F: np.ravel( (FideosCompProf(p,X,Y)-F) )
	
	p0 = np.array([x1,y1,np.max(F)])
        p1, success = scipy.optimize.leastsq(errfunc,p0, args=(X,Y,F))
	print p0
	print p1
	imshow(F)
	show()
	imshow(FideosCompProf(p1,X,Y))
	show()
	return p1


def lnlike(theta, X,Y,F, Ferr):
	model = FideosCompProf(theta,X,Y)
	inv_sigma2 = 1.0/(Ferr**2)
	#ret = -np.log(2*np.pi) + np.log(np.sum(np.exp(-0.5*((y-model)/yerr)**2)/yerr))
	ret = -0.5*(np.sum(inv_sigma2*(F-model)**2 - np.log(inv_sigma2)))
	#for i in range(len(F)):
	#	errorbar(Y,F[i],yerr=Ferr[i],fmt='b')
	#for j in model:
	#	plot(Y,j,'r')
	#show()
	#print theta, ret
	if np.isnan(ret):
		return -np.inf
	else:
		return ret


def lnprior(theta,X,Y):
	if X[0] < theta[0] < X[-1] and Y[0] < theta[1] < Y[-1] and 0 < theta[2] < 70000 and 0.1 < theta[3] < 5 and  0.1 < theta[4] < 20:
	    return 0.0
	return -np.inf

def lnprob(theta, X,Y,F,Ferr):
    lp = lnprior(theta,X,Y)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta,X,Y,F,Ferr)

def AnaliticProfile1D(p,x,R=2.3):
    space = 100.
    xt    = np.arange(x[0],x[-1]+1,1./space)
    x0 = xt - p[0]
    J = np.where((R**2-x0**2>=0) & (x0>=0))[0]
    out = np.zeros(len(x0))
    out[J] = 2*np.sqrt(R**2-x0[J]**2)
    out = scipy.ndimage.filters.gaussian_filter(out, p[2]*space)
    out /= np.max(out)
    out = p[1]*out
    out = np.mean(out.reshape((len(out)/space,space)),axis=1)
    return out


def get_airmass(ra,dec,lat,lon,alt,dateobs):

	target = SkyCoord(ra,dec, unit="deg")  
	obs    = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=alt*u.m)
	#utcoffset = -4*u.hour  # Eastern Daylight Time
	time = Time(dateobs)
	altaz = target.transform_to(AltAz(obstime=time,location=obs))  

	return altaz.secz

def compute_RON(mbias,biases,gain=1.):
	stdevs = []
	for fbias in biases:
		bias = pyfits.getdata(fbias)
		tbias = bias-mbias
		stdevs.append(np.sqrt(np.var(tbias)))
	stdevs = np.array(stdevs)
	return np.median(stdevs)*gain
		
def fill_list_header(hdu,ofiles):
    for i in range(len(ofiles)):
        hdu = GLOBALutils.update_header(hdu,'RAW'+str(i),ofiles[i].split('/')[-1])
    return hdu