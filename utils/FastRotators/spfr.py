from pylab import *
import pyfits
from PyAstronomy import pyasl
import scipy
from scipy import interpolate
from scipy import ndimage
from scipy import signal
import pickle
from matplotlib.backends.backend_pdf import PdfPages

#from pyevolve import G1DList
#from pyevolve import GSimpleGA

from multiprocessing import Pool
import time

def download_models(webpage='http://svo2.cab.inta-csic.es/theory/models/coelho/high/data/',dest='../../data/'):
	
	os.system('mkdir '+dest+'/COELHO2014')
	cwd = os.getcwd()
	os.chdir(dest+'/COELHO2014')

	tf = np.arange(6000,1001,250)
	gf = np.arange(3.0,4.6,0.5)
	zf = np.array([-1.,-0.5,0.0,0.2])

	for t in tf:
		for g in gf:
			for z in zf:
				modname = get_modname(t,g,z)
				os.system('wget ' + webpage+modname+'.fits')
				os.system('wget ' + webpage+modname+'plc.fits')
	os.chdir(cwd)

	return True

def n_Edlen(l):
    sigma = 1e4 / l
    sigma2 = sigma*sigma
    n = 1 + 1e-8 * (8342.13 + 2406030 / (130-sigma2) + 15997/(38.9-sigma2))
    return n

def n_Morton(l):
    sigma = 1e4 / l
    sigma2 = sigma*sigma
    n = 1 + 6.4328e-5 + 2.94981e-2 / (146.-sigma2) + 2.5540e-4/(41.-sigma2)
    return n

def ToAir(l):
    return (l / n_Edlen(l))

def ToVacuum(l):
    cond = 1
    l_prev = l.copy()
    while(cond):
        l_new = n_Edlen(l_prev) * l
        if (max(np.absolute(l_new - l_prev)) < 1e-10): cond = 0
        l_prev = l_new
    return l_prev

def get_modname(t,g,z):
	st = str(int(t))
	if t<10000:
		st = '0'+st
	sg = '+'+str(np.around(g,1))
	if z < 0:
		sz = 'm'
	else:
		sz = 'p'
	z=float(z)
	sz = sz + str(np.around(np.absolute(z),1))
	sz = sz.replace('.','')
	return 't'+st+'_g'+sg+'_'+sz+'p00_hr'

def get_model(t,g,z,model_path='/Users/rabrahm/data/COELHO2014/'):
	modname = model_path + get_modname(t,g,z)
	try:
		out = pyfits.getdata(modname+'.fits')
	except:
		out = pyfits.getdata(modname+'plc.fits')
	return out

def get_near(x,vec):

	if x == vec[0]:
		mmin = vec[0]
		mmax = vec[1]
	elif x == vec[-1]:
		mmin = vec[-2]
		mmax = vec[-1]
	else:
		tvec = vec - x
		In  = np.where(tvec < 0)[0]
		mmin = tvec[In].max() + x
		Ix = np.where(tvec >= 0)[0]
		mmax = tvec[Ix].min() + x
	return mmin,mmax

def trilinear_interpolation(t,g,z,model_path='/Users/rabrahm/data/COELHO2014/'):
	teffs = np.arange(6000,10001,250)
	loggs = np.arange(3.0,4.6,0.5)
	fehs  = np.array([-1.,-0.5,0.0,0.2])
	x0,x1 = get_near(t,teffs)
	y0,y1 = get_near(g,loggs)
	z0,z1 = get_near(z,fehs)
	xd = (t-x0)/(x1-x0)
	yd = (g-y0)/(y1-y0)
	zd = (z-z0)/(z1-z0)

	try:
		hd = pyfits.getheader(model_path+get_modname(x0,y0,z0)+'.fits')
	except:
		hd = pyfits.getheader(model_path+get_modname(x0,y0,z0)+'plc.fits')

	c000 = get_model(x0,y0,z0,model_path)
	c001 = get_model(x0,y0,z1,model_path)
	c010 = get_model(x0,y1,z0,model_path)
	c100 = get_model(x1,y0,z0,model_path)
	c110 = get_model(x1,y1,z0,model_path)
	c101 = get_model(x1,y0,z1,model_path)
	c011 = get_model(x0,y1,z1,model_path)
	c111 = get_model(x1,y1,z1,model_path)

	wav = np.arange(len(c111))*hd['CDELT1'] + hd['CRVAL1']

	c00 = c000*(1-xd) + c100*xd
	c01 = c001*(1-xd) + c101*xd
	c10 = c010*(1-xd) + c110*xd
	c11 = c011*(1-xd) + c111*xd

	c0 = c00*(1-yd) + c10*yd
	c1 = c01*(1-yd) + c11*yd

	c = c0*(1-zd) + c1*zd

	return wav,c 

def normalize_model(w,f):
	ow = w.copy()
	of = f.copy()
	#plot(w,f)
	while True:
		#medflts = scipy.signal.medfilt(f,1001)
		coef = np.polyfit(w,f,6)
		fited = np.polyval(coef,w)
		res = f - fited
		I = np.where(res > -np.sqrt(np.var(res)))[0]
		w,f = w[I],f[I]
		if len(w) < 0.3* len(ow):
			break
	#plot(ow,np.polyval(coef,ow))
	#show()
	return coef

def spec_ccf(sw,sf,mw,mf,vi,vf,dv):
	mf = mf -1
	mf = -mf
	#plot(mw,mf)
	tck = interpolate.splrep(mw,mf,k=1)
	v = vi
	retccf = []
	vels = []
	while v<=vf:
		swt = sw * (1 + v/299792.458)
		mft = interpolate.splev(swt,tck)
		mft -= np.mean(mft)
		sft = sf - np.mean(sf)
		#sft = sf.copy()
		retccf.append(np.sum(mft*sft)/np.sqrt(np.sum(mft**2)*np.sum(sft**2)))
		vels.append(v)
		v+=dv

	return np.array(vels),np.array(retccf)

def ccf_simple(sw,sf,mw,mf,rv):
	mf = mf -1
	mf = -mf
	#plot(mw,mf)
	tck = interpolate.splrep(mw,mf,k=1)
	swt = sw * (1 + rv/299792.458)
	mft = interpolate.splev(swt,tck)
	mft -= np.mean(mft)
	sft = sf - np.mean(sf)
	return np.sum(mft*sft)/np.sqrt(np.sum(mft**2)*np.sum(sft**2))

def clean_strong_lines(mw,sc):
	#""""
	I = np.where((mw>6520)&(mw<6600))[0]
	sc[I] = 1.
	I = np.where((mw>5888)&(mw<5897))[0]
	sc[I] = 1.
	I = np.where((mw>4310)&(mw<4360))[0]
	sc[I] = 1.
	I = np.where((mw>4840)&(mw<4880))[0]
	sc[I] = 1.
	I = np.where((mw>4070)&(mw<4130))[0]
	sc[I] = 1.
	I = np.where((mw>3875)&(mw<3900))[0]
	sc[I] = 1.
	I = np.where((mw>3920)&(mw<3945))[0]
	sc[I] = 1.
	I = np.where((mw>3955)&(mw<3980))[0]
	sc[I] = 1.
	I = np.where(mw<3850)[0]
	sc[I] = 1.
	#"""
	return sc

def RVforFR(wavs,flxs,teff=6700,logg=4.0,feh=-1.0,vsini=100.,model_path='/Users/rabrahm/data/COELHO2014/',vmin=-1000.,vmax=1000.,vstep=10.):
	def fitfunc(p,x):
		ret = p[3] + p[0] * np.exp(-.5*((x-p[1])/p[2])**2)
		return ret
	errfunc = lambda p,x,y: np.ravel( (fitfunc(p,x)-y) )

	#sc = get_model(teff,logg,feh)
	#hd = pyfits.getheader(model_path+get_modname(7000,4.5,0.0)+'.fits')
	#wav = np.arange(len(sc))*hd['CDELT1'] + hd['CRVAL1']
	teff = float(teff)

	try:
		sc = get_model(teff,logg,feh)
		hd = pyfits.getheader(model_path+get_modname(7000,4.5,0.0)+'.fits')
		mw = np.arange(len(sc))*hd['CDELT1'] + hd['CRVAL1']
	except:
		mw,sc = trilinear_interpolation(teff,logg,feh,model_path)

	sc = clean_strong_lines(mw,sc)

	mw = ToVacuum(mw)
	ccftot = []

	II = np.where(sc != 1)[0]
	coef = normalize_model(mw[II],sc[II])
	sc  /= np.polyval(coef,mw)

	for i in range(wavs.shape[0]):
		scf = flxs[i]
		scw = wavs[i]

		J = np.where(scf!=0)[0]
		scw,scf = scw[J],scf[J]
		I = np.where((mw>scw[0]-100) & (mw<scw[-1]+100))
		tmf = pyasl.fastRotBroad(mw[I], sc[I], 0.5, vsini)
		#plot(mw[I],tmf)
		ccv,ccf = spec_ccf(scw,scf,mw[I],tmf,vmin,vmax,vstep)
		ccf = np.array(ccf)
		if len(ccftot)==0:
			ccftot = ccf.copy()
		else:
			ccftot = np.vstack((ccftot,ccf))

	ccftot = np.mean(ccftot,axis=0)

	p0 = [ccftot.min(),ccv[np.argmin(ccftot)],vsini,ccftot[0]]
	p1, success = scipy.optimize.leastsq(errfunc,p0, args=(ccv,ccftot))

	return p1,ccv,ccftot,fitfunc(p1,ccv)
	
def calc_bss2(vels,xc,coef, bot_i=0.15, bot_f=0.4, top_i=0.6, top_f=0.9, dt=0.01):
	try:
		
		I1 = np.where((vels>coef[1]-3*coef[2]) & (vels<coef[1]) )[0]
		I2 = np.where((vels<coef[1]+3*coef[2]) & (vels>coef[1]) )[0]
		I3 = np.where(vels<coef[1]-4*coef[2])[0]
		I4 = np.where(vels>coef[1]+4*coef[2])[0]
		I = np.hstack((I3,I4))
		base = np.median(xc[I])

		xc = base - xc
		xc /= xc.max()


		v1,x1 = vels[I1],xc[I1]
		v2,x2 = vels[I2],xc[I2]
		#plot(v1,x1)
		#plot(v2,x2)
		#show()
		dp = top_f
		vect = []
		while dp >= top_i:
			lb = np.where(x1>dp)[0][0]
			m = (v1[lb] - v1[lb-1])/(x1[lb]-x1[lb-1])
			n = v1[lb] - m*x1[lb]
			bs1 = m*dp+n

			lb = np.where(x2>dp)[0][-1]
			m = (v2[lb] - v2[lb+1])/(x2[lb]-x2[lb+1])
			n = v2[lb] - m*x2[lb]
			bs2 = m*dp+n
			vect.append(0.5*(bs2+bs1))
			dp-=dt
		vect = np.array(vect)

		dp = bot_f
		vecb = []
		while dp >= bot_i:

			lb = np.where(x1>dp)[0][0]
			m = (v1[lb] - v1[lb-1])/(x1[lb]-x1[lb-1])
			n = v1[lb] - m*x1[lb]
			bs1 = m*dp+n

			lb = np.where(x2>dp)[0][-1]
			m = (v2[lb] - v2[lb+1])/(x2[lb]-x2[lb+1])
			n = v2[lb] - m*x2[lb]
			bs2 = m*dp+n
			vecb.append(0.5*(bs2+bs1))
			dp-=dt
		vecb = np.array(vecb)

		return np.median(vecb) - np.median(vect) 
	except:
		return -999.0
"""
def lnlike(theta, W, F, Ferr):
	mw,sc = trilinear_interpolation(int(theta[0]),theta[1],theta[2])
	sct  = clean_strong_lines(mw,sc.copy())
	#plot(mw,sc)
	#show()
	coef = normalize_model(mw,sct)
	sc /= np.polyval(coef,mw)

	#print gfd
	mw = ToVacuum(mw)
	mw *= 1 + theta[3]/299792.458

	totD,totM,totE = np.array([]),np.array([]),np.array([])
	for i in range(W.shape[0]):
		scf = F[i]
		scw = W[i]
		scfe = Ferr[i]

		J = np.where(scf!=0)[0]
		scw,scf,scfe = scw[J],scf[J],scfe[J]

		I = np.where((mw>scw[0]-10) & (mw<scw[-1]+10))
		tmf = pyasl.fastRotBroad(mw[I], sc[I], 0.5, theta[4])
		tck = interpolate.splrep(mw[I],tmf,k=1)
		tmf = interpolate.splev(scw,tck)

		tmf = clean_strong_lines(scw,tmf.copy())
		I = np.where(tmf!=1)[0]
		#plot(scw,tmf)
		#plot(scw[I],tmf[I])
		#plot(scw[I],scf[I])

		#show()
		#print gfd

		tmf = tmf[I]
		scf = scf[I]
		scfe = scfe[I]

		tmf /= np.sum(tmf)
		tsf = scf/np.sum(scf)
		tse = scfe*(np.sum(scf)**2)

		totD = np.hstack((totD,tsf))
		totM = np.hstack((totM,tmf))
		totE = np.hstack((totE,tse))
		#plot(scw[I],tsf)
		#plot(scw[I],tmf)
		#plot(scw[I],tsf + 1./np.sqrt(tse))

	#show()
	#print fds
	#print theta
	#show()
	#print gvfd
	#ret = -np.log(2*np.pi) + np.log(np.sum(np.exp(-0.5*((y-model)/yerr)**2)/yerr))
	#ret = -0.5*(np.sum(inv_sigma2*(F-model)**2 - np.log(inv_sigma2)))
	ret = -0.5*(np.sum(totE*(totD-totM)**2 - np.log(totE)))

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


def lnprior(theta):
	if 6000 < theta[0] < 9000 and 3.0 < theta[1] < 4.5 and -1 < theta[2] < 0.2 and -500 < theta[3] < 500 and  1. < theta[4] < 500.:
	    return 0.0
	return -np.inf

def lnprob(theta, W,F,Ferr):
	lp = lnprior(theta)
	if not np.isfinite(lp):
		return -np.inf
	return lp + lnlike(theta,W,F,Ferr)
"""
def multiccf(pars):
	teff,logg,feh,vsini=pars[0],pars[1],pars[2],pars[3]
	vmin=-500
	vmax=500.
	vstep=20.
	sc = get_model(teff,logg,feh)
	hd = pyfits.getheader(model_path+get_modname(7000,4.5,0.0)+'.fits')
	wav = np.arange(len(sc))*hd['CDELT1'] + hd['CRVAL1']
	try:
		sc = get_model(teff,logg,feh)
		hd = pyfits.getheader(model_path+get_modname(7000,4.5,0.0)+'.fits')
		mw = np.arange(len(sc))*hd['CDELT1'] + hd['CRVAL1']
	except:
		mw,sc = trilinear_interpolation(teff,logg,feh,model_path)

	sc = clean_strong_lines(mw,sc)

	mw = ToVacuum(mw)
	ccftot = []

	II = np.where(sc != 1)[0]
	coef = normalize_model(mw[II],sc[II])
	sc  /= np.polyval(coef,mw)

	for i in range(wavs.shape[0]):
		scf = flxs[i]
		scw = wavs[i]

		J = np.where(scf!=0)[0]
		scw,scf = scw[J],scf[J]
		I = np.where((mw>scw[0]-100) & (mw<scw[-1]+100))
		tmf = pyasl.fastRotBroad(mw[I], sc[I], 0.5, vsini)
		#plot(mw[I],tmf)
		ccv,ccf = spec_ccf(scw,scf,mw[I],tmf,vmin,vmax,vstep)
		ccf = np.array(ccf)
		if len(ccftot)==0:
			ccftot = ccf.copy()
		else:
			ccftot = np.vstack((ccftot,ccf))

	ccftot = np.mean(ccftot,axis=0)
	return ccftot.min()


def get_pars_fr(wavst,flxst,model_patht='/Users/rabrahm/data/COELHO2014/',npools=4):
	t0 = time.time()
	global wavs,flxs
	global model_path
	wavs,flxs=wavst.copy(),flxst.copy()
	model_path=model_patht


	gt = np.array([6000,7000,8000,9000,10000])
	gg = np.array([3.0,3.5,4.0,4.5])
	gz = np.array([-1,-0.5,0.0,0.2])
	gr = np.array([10.,50.,100.,150.,200.,250.,300.])

	#"""
	tr = np.tile(gr,len(gt)*len(gg)*len(gz))
	tg = np.repeat(np.tile(gg,len(gt)),len(gr)*len(gz))
	tz = np.repeat(np.tile(gz,len(gt)*len(gg)),len(gr))
	tt = np.repeat(gt,len(gg)*len(gr)*len(gz))
	tot = np.vstack((tt,tg,tz,tr)).T

	p = Pool(npools)
	vals = np.array((p.map(multiccf, list(tot))))
	p.terminate()
	I = np.argmin(vals)
	best_vals = tot[I]
	bt,bg,bz,br = best_vals[0],best_vals[1],best_vals[2],best_vals[3]
	#"""
	t1 = time.time()
	print bt,bg,bz,br, (t1-t0)/60.,'mins'


	#bt,bg,bz,br = 7000.,4.5, 0.2, 100.0
	gt = np.arange(bt-1000,bt+1001,250)
	I = np.where((gt>=6000) & (gt<=10000))[0]
	gt = gt[I]
	gr = np.arange(br-60.,br+61.,20.)
	I = np.where(gr>=10)[0]
	gr = gr[I]

	tr = np.tile(gr,len(gt)*len(gg)*len(gz))
	tg = np.repeat(np.tile(gg,len(gt)),len(gr)*len(gz))
	tz = np.repeat(np.tile(gz,len(gt)*len(gg)),len(gr))
	tt = np.repeat(gt,len(gg)*len(gr)*len(gz))
	tot = np.vstack((tt,tg,tz,tr)).T

	p = Pool(npools)
	vals = np.array((p.map(multiccf, list(tot))))
	p.terminate()
	I = np.argmin(vals)
	best_vals = tot[I]
	bt,bg,bz,br = best_vals[0],best_vals[1],best_vals[2],best_vals[3]
	t2 = time.time()
	print bt,bg,bz,br, (t2-t1)/60.,'mins'
	#np.savetxt('temp_grid.txt',vals)

	grid = np.reshape(vals,(len(gt),len(gg),len(gz),len(gr)))
	tckt = interpolate.splrep(gt,np.arange(len(gt)),k=1)
	tckg = interpolate.splrep(gg,np.arange(len(gg)),k=1)
	tckz = interpolate.splrep(gz,np.arange(len(gz)),k=1)
	tckr = interpolate.splrep(gr,np.arange(len(gr)),k=1)

	itckt = interpolate.splrep(np.arange(len(gt)),gt,k=1)
	itckg = interpolate.splrep(np.arange(len(gg)),gg,k=1)
	itckz = interpolate.splrep(np.arange(len(gz)),gz,k=1)
	itckr = interpolate.splrep(np.arange(len(gr)),gr,k=1)

	st = np.arange(gt[0],gt[-1]+1,10.)
	sg = np.arange(gg[0],gg[-1]+0.01,0.1)
	sz = np.arange(gz[0],gz[-1]+0.01,0.1)
	sr = np.arange(gr[0],gr[-1]+1.,5.)

	st = interpolate.splev(st,tckt)
	sg = interpolate.splev(sg,tckg)
	sz = interpolate.splev(sz,tckz)
	sr = interpolate.splev(sr,tckr)

	tr2 = np.tile(sr,len(st)*len(sg)*len(sz))
	tg2 = np.repeat(np.tile(sg,len(st)),len(sr)*len(sz))
	tz2 = np.repeat(np.tile(sz,len(st)*len(sg)),len(sr))
	tt2 = np.repeat(st,len(sg)*len(sr)*len(sz))
	tot2 = np.vstack((tt2,tg2,tz2,tr2))

	zi = ndimage.map_coordinates(grid, tot2, order=3, mode='nearest')
	I  = np.argmin(zi)
	minval = tot2[:,I]

	mint = interpolate.splev(minval[0],itckt)
	ming = interpolate.splev(minval[1],itckg)
	minz = interpolate.splev(minval[2],itckz)
	minr = interpolate.splev(minval[3],itckr)

	#d = {'grid':grid, 'zi':zi, 'tot2':tot2, 'gt':gt, 'gg':gg, 'gz':gz, 'gr':gr}
	#pickle.dump(d,open('temp_dict.pkl'))

	return float(mint),float(ming),float(minz),float(minr)

def plot_CCF_FR(xc_dict,path='XC.pdf'):

	vels       = xc_dict['vels']
	xc_av      = xc_dict['xc_av']
	XCmodelgau = xc_dict['XCmodelgau']
	#refvel     = xc_dict['refvel']
	p1gau      = xc_dict['p1gau']

	f1 = figure()

	pp = PdfPages(path)
	ax1 = f1.add_subplot(111)
	ax1.plot(vels, xc_av,'b.', label='CCF')
	ax1.plot(vels, XCmodelgau,'r-',label='Gaussian fit')

	xlabel('Velocity (km/s)')
	ylabel('XC')
	ax1.axvline(p1gau[1],linestyle=':',color='r')
	ax1.axhline(0.0,linestyle='-')
	title('Average Cross-Correlation Function + Fit')
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles[::-1], labels[::-1],prop={'size':6})
	pp.savefig()
	pp.close()
	clf()
	pass

"""
def trans_chromosome(chromosome):
	teff = chromosome[0]*100.+chromosome[1]*10.+chromosome[2]
	m = (10000.- 6000.)/999.
	n = 6000. 
	teff = teff*m + n
	logg = chromosome[3] + chromosome[4]*0.1
	m = (4.5 - 3.0)/9.9
	n = 3. 
	logg = logg*m + n

	feh  = chromosome[5] + chromosome[6]*0.1
	m = (0.2 - -1.)/9.9
	n = -1.
	feh = feh*m + n

	vsini  = chromosome[7]*10. + chromosome[8]
	m = (300. - 10.)/99.
	n = 10.
	vsini = vsini*m + n

	return teff, logg, feh, vsini

global wavs, flxs

def find_pars_GA(wavs,flxs,model_path='/Users/rabrahm/data/COELHO2014/'):

	def eval_func(chromosome):
		print list(chromosome)
		teff, logg, feh, vsini = trans_chromosome(chromosome)
		print teff, logg, feh, vsini
		pt,vels,ccf,mod = RVforFR(wavs,flxs,teff=teff,logg=logg,feh=feh,vsini=vsini,model_path=model_path)
		score = -ccf.min()
		return score

	genome = G1DList.G1DList(9)
	genome.evaluator.set(eval_func)

	ga = GSimpleGA.GSimpleGA(genome, interactiveMode=True)
	ga.setGenerations(40)
	ga.setMutationRate(0.2)
	ga.setPopulationSize(20)
  	#ga.setCrossoverRate(1.0)
	genome.setParams(rangemin=0, rangemax=9)
	#ga.setMultiProcessing(True)
	ga.evolve(freq_stats=10)
	print ga.bestIndividual()
	print trans_chromosome(ga.bestIndividual())

"""




