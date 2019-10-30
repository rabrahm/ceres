try:
    import pyfits
except:
    from astropy.io import fits as pyfits

from pylab import *
import scipy
from scipy import optimize, interpolate
import pickle
import sys
sys.path.append("../utils/GLOBALutils")
import GLOBALutils
from multiprocessing import Pool


def clipp(vec,n=5.):
	I = np.arange(len(vec))
	res = np.zeros(len(vec))
	while True:
		mean = np.mean(vec[I])
		dev = np.sqrt(np.var(vec[I]))
		res[I] = np.absolute(vec[I] - mean)
		J = np.where(res>n*dev)[0]
		if len(J)==0:
			break
		imax = np.argmax(res)
		res[imax] = -999
		Z = np.where(I==imax)[0][0]
		I = np.delete(I,Z)
	return I

def xcor(w,f1,f2,vi=-1000,vf=1000,dv=10):
    vels = np.arange(vi,vf,dv)
    xcors = []
    for v in vels:
	tw = w * ( 1 + v / 299792458. )
        tck = interpolate.splrep(tw,f2,k=3)
        tf2 = interpolate.splev(w,tck)

	N1 = np.sum(f1**2)
	N2 = np.sum(tf2**2)
	CCF = np.sum(f1*tf2)/(np.sqrt(N1*N2))
	xcors.append(CCF)
	
    #plot(vels,xcors)
    #show()
    return vels, np.array(xcors)
    

def ccf_chunks(wave,fp1,fp2,lim1=50,lim2=-50,chsize=50):

    if lim2 < 0:
	lim2 = len(fp1) + lim2
    rvs = []
    wave,fp1,fp2 =  wave[lim1:lim2],fp1[lim1:lim2],fp2[lim1:lim2]
    i=0
    while i < len(fp1):	
	w,f1,f2 = wave[i:i+chsize], fp1[i:i+chsize], fp2[i:i+chsize]

	coef = np.polyfit(w,f1,1)
	f1 = f1 / np.polyval(coef,w)
	f2 = f2 / np.polyval(coef,w)
	
	vels,ccf = xcor(w,f1,f2)
	tvels = np.arange(vels[0],vels[-1],0.01)
	tck = interpolate.splrep(vels,ccf,k=3)
	tccf = interpolate.splev(tvels,tck)
	rv = tvels[np.argmax(tccf)]
	rvs.append(rv)
	if i == len(fp1):
	    break
	elif i + 2*chsize > len(fp1):
	    i = len(fp1)
	else:
	    i += chsize
    rvs = np.array(rvs)
    I = clipp(rvs,n=3)
    rvs = rvs[I]
    return rvs
     

def FitFP(X, Y, p0):
    def AnaliticProfile1D(p,x):
	m = p[0]
	n = p[1]
	pt = p[2:]
        space = 100.
        xt    = np.arange(x[0],x[-1]+1,1./space)
	nlines = len(pt)/3
	base = xt*m + n
	i=0
	while i < len(pt):
		base += pt[i]*np.exp(-0.5*((xt-pt[2+i])/pt[1+i])**2)
		i+=3
        out = np.mean(base.reshape((int(len(base)/space),int(space))),axis=1)
        return out
    errfunc = lambda p,X,Y: np.ravel( (AnaliticProfile1D(p,X)-Y) )

    p1, success = scipy.optimize.leastsq(errfunc,p0, args=(X,Y))

    return p1

def clean(x,vec):
    res  = np.zeros(len(vec))
    while True:
	I = np.where(res!=-999)[0]
	coef = np.polyfit(x[I],vec[I],1)
	y = np.polyval(coef,x[I])
	res[I] = vec[I] - y
	dev = np.sqrt(np.var(res[I]))
	J = np.where((np.absolute(res)>3*dev) & (res != -999))[0]
	if len(J)==0:
	    break
	else:
	    argmax = np.argmax(np.absolute(res[J]))
	    argmax = J[argmax]
	    res[argmax] = -999
    J = np.where(res!=-999)[0]
    return coef

def InitialGuess(path, lim1=50,lim2=-50, oi=0,of=-1):
    sc = pyfits.getdata(path)
    if lim2 < 0:
	lim2 = sc.shape[2] + lim2
    if of < 0:
	of = sc.shape[0] + of

    pxs = {}
    for order in np.arange(oi,of,1):
	fp = sc[order,1]
	ejex = np.arange(len(fp))
	fp1 = np.hstack((fp[-1],fp[0:-1]))
	fp2 = np.hstack((fp[1:],fp[0]))
	fp3 = np.hstack((fp[-2:],fp[0:-2]))
	fp4 = np.hstack((fp[2:],fp[:2]))
	fp5 = np.hstack((fp[-3:],fp[0:-3]))
	fp6 = np.hstack((fp[3:],fp[:3]))
	I1 = np.where((fp1<fp) & (fp2<=fp) & (fp3<fp) & (fp4<fp) & (fp5<fp) & (fp6<fp))[0]
	I2 = np.hstack((I1[-1],I1[0:-1]))
	I = I1-I2
	coef = clean(I1,I)
	I3 = np.where((fp1<fp) & (fp2<=fp) & (fp3<fp) & (fp4<fp) & (fp5<fp) & (fp6<fp) & (ejex>0.5*len(fp)-50) & (ejex<0.5*len(fp)+50))[0][0]
	#print I3
	i= I3
	vec1 = [I3]
	#plot(sc[0,order],sc[1,order])
	while i < lim2:
	     vec1.append(i)
	     i += int(np.around(np.polyval(coef,i)))
	     try:
	        i = i - 4 + np.argmax(fp[i-4:i+4])
	     except:
		break
	i = I3
	while i > lim1:
	     vec1.append(i)
	     i -= int(np.around(np.polyval(coef,i)))
	     try:
	        i = i - 4 + np.argmax(fp[i-4:i+4])
	     except:
		break

	vec1 = np.sort(np.unique(np.around(np.array(vec1)).astype('int')))
	dist = 0.5*(vec1[-2] - vec1[-3])
	if vec1[-1] + dist > lim2:
		vec1 = np.delete(vec1,-1)
	dist = 0.5*(vec1[3] - vec1[2])
	if vec1[0] - dist < lim1:
		vec1 = np.delete(vec1,0)
	#plot(fp)
	#plot(vec1.astype('int'),fp[vec1.astype('int')],'ro')
	#show()
	pxs['order_'+str(int(order))] = vec1.astype('float')

    return pxs

def ParallelFit(i):
    try:
	    mid  = int(np.around(0.5*len(vec1)))
	    dist = int(np.around(0.5*(vec1[mid]-vec1[mid-1])))
	    if i == 0:
		#dist = int(np.around(0.5*(vec1[i+1]-vec1[i])))
		ejx = np.arange( int(np.around(vec1[i]))-dist,int(np.around(vec1[i+1]))+dist+1)
		ejy = sc[order,1,ejx]	
		temp = ejy - np.min(ejy)
		p0 = [0,np.min(ejy),temp[int(vec1[i]-ejx[0])],1.,vec1[i],temp[int(vec1[i+1]-ejx[0])],1.,vec1[i+1]]
	    elif i == len(vec1)-1:
		#dist = int(np.around(0.5*(vec1[i]-vec1[i-1])))
		ejx = np.arange(int(np.around(vec1[i-1]))-dist,int(np.around(vec1[i]))+dist+1)
		ejy = sc[order,1,ejx]
		temp = ejy - np.min(ejy)
		p0 = [0,np.min(ejy),temp[int(vec1[i-1]-ejx[0])],1.,vec1[i-1],temp[int(vec1[i]-ejx[0])],1.,vec1[i]]
	    else:
		#dist = int(np.around(0.5*(vec1[i]-vec1[i-1])))
		ejx = np.arange(int(np.around(vec1[i-1]))-dist,int(np.around(vec1[i+1]))+dist+1)
		ejy = sc[order,1,ejx]
		#print vec1[i-1],vec1[i],vec1[i+1],dist, ejx, ejy
		temp = ejy - np.min(ejy)
		p0 = [0,np.min(ejy),temp[int(vec1[i-1]-ejx[0])],1.,vec1[i-1],temp[int(vec1[i]-ejx[0])],1.,vec1[i],temp[int(vec1[i+1]-ejx[0])],1.,vec1[i+1]]
	    coefs = FitFP(ejx,ejy,p0)
	    pt = coefs[2:]
	    if i == 0:
		#amps.append(pt[0])
		#sigs.append(pt[1])
		#meds.append(pt[2])
		return np.array([pt[2],pt[0]])
	    else:
		#amps.append(pt[3])
		#sigs.append(pt[4])
		#meds.append(pt[5])
		return np.array([pt[5],pt[3]])
    except:
	return [-999,-999]


def GetFPLines(path, vec, lim1=50,lim2=-50,npools=1,oi=0,of=-1):
    global sc,vec1,order
    sc = pyfits.getdata(path)
    if lim2 < 0:
	lim2 = sc.shape[2] + lim2
    if of < 0:
	of = sc.shape[0] + of
    pxs = {}
    for order in np.arange(oi,of,1):
	#print order
	#order=48
	meds,sigs = [],[]
	vec1 = vec['order_'+str(int(order))]
	ies = np.arange(len(vec1))
	#for i in ies:
	#   print i,ParallelFit(i)
	#print gfd
	p = Pool(npools)
	mat = np.array((p.map(ParallelFit, ies)))

	p.terminate()
	meds = mat[:,0]
        amps = mat[:,1]
	"""
	for i in np.arange(len(vec1)):
		if i == 0:
			dist = int(np.around(0.5*(vec1[i+1]-vec1[i])))
			ejx = np.arange( int(np.around(vec1[i]))-dist,int(np.around(vec1[i+1]))+dist+1)
			ejy = sc[order,1,ejx]	
			temp = ejy - np.min(ejy)
			p0 = [0,np.min(ejy),temp[vec1[i]-ejx[0]],1.,vec1[i],temp[vec1[i+1]-ejx[0]],1.,vec1[i+1]]
		elif i == len(vec1)-1:
			dist = int(np.around(0.5*(vec1[i]-vec1[i-1])))
			ejx = np.arange(int(np.around(vec1[i-1]))-dist,int(np.around(vec1[i]))+dist+1)
			ejy = sc[order,1,ejx]
			temp = ejy - np.min(ejy)
			p0 = [0,np.min(ejy),temp[vec1[i-1]-ejx[0]],1.,vec1[i-1],temp[vec1[i]-ejx[0]],1.,vec1[i]]
		else:
			dist = int(np.around(0.5*(vec1[i]-vec1[i-1])))
			ejx = np.arange(int(np.around(vec1[i-1]))-dist,int(np.around(vec1[i+1]))+dist+1)
			ejy = sc[order,1,ejx]
			temp = ejy - np.min(ejy)
			p0 = [0,np.min(ejy),temp[vec1[i-1]-ejx[0]],1.,vec1[i-1],temp[vec1[i]-ejx[0]],1.,vec1[i],temp[vec1[i+1]-ejx[0]],1.,vec1[i+1]]

		coefs = FitFP(ejx,ejy,p0)
		pt = coefs[2:]
		if i == 0:
			sigs.append(pt[1])
			meds.append(pt[2])
		else:
			sigs.append(pt[4])
			meds.append(pt[5])
	    """
	#sigs,meds = np.array(sigs),np.array(meds)

	#m = order + or0 + 22
	#chebs = GLOBALutils.Calculate_chebs(meds, m, Inverse=True,order0=or0,ntotal=n_useful,npix=sc.shape[2],nx=ncoef_x,nm=ncoef_m)
	#WavSol =  (1.0 + 1.0e-6*wavsol_dict['p_shift']) * (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(wavsol_dict['p1_co'],chebs,ncoef_x,ncoef_m)
	I = np.where(meds>sc.shape[2]-1)[0]
	meds[I] = float(sc.shape[2]-2)
	pxs['order_'+str(int(order))] = meds
        pxs['order_amps_'+str(int(order))] = amps
        #pxs['order_'+str(int(order))+'_sigs'] = sigs
	#print meds
	#plot(sc[order,1])
	#plot(meds,sc[order,1,np.around(meds).astype('int')],'ro')
	#show()

    return pxs

