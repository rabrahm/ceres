from pylab import *
import scipy
from scipy import interpolate
from scipy import signal
sys.path.append("../utils/GLOBALutils")
import GLOBALutils
import os 
import glob
import pyfits

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
			elif 'ThAr' in archivo and '300' in archivo:
				if '1F' in archivo:
					thar_co.append(archivo)
				elif '2F' in archivo:
					thar.append(archivo)
			else:
				if '1F' in archivo:
					dar_sci.append(archivo)
				else:
					sim_sci.append(archivo)

				#line = "%-15s %10s %10s %8.2f %4.2f %8s %11s %s\n" % (obname, ra, delta, texp, airmass, date, hour, archivo)
				#f.write(line)

	f.close()
	biases, darks, flats, flats_co, thar, thar_co, sim_sci, dar_sci = \
		np.array(biases), np.array(darks), np.array(flats), np.array(flats_co), np.array(thar), \
		np.array(thar_co), np.array(sim_sci), np.array(dar_sci)
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
				bas = np.min(flat[cen-3*exap:cen+3*exap,x])
			except:
				bas = baso
			img_out[cen-exap:cen+exap+1,x] = bas
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
    plot(d)
    show()
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

	
    for i in range(len(c_ob)):
        Centers[i,:]=scipy.polyval(c_ob[i],ejx)
	plot(ejx,Centers[i,:],'r')
    for i in range(len(c_co)):
        Centers[i,:]=scipy.polyval(c_co[i],ejx)
	plot(ejx,Centers[i,:],'b')
    show()
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

    print 'Warning!!! adding 5 hrs to comute MJD due to problem in header! CHECK in future!!'
    mjd_start += 5./24.


    secinday = 24*3600.0
    fraction = 0.5
    texp     = h[0].header['EXPTIME'] #sec

    print 'Warning!!! assuming that the provided date is at the end of the exposure!!! Check in the future!!'
    mjd = mjd_start - (fraction * texp) / secinday

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
	plot(line)
	plot(refmins,valmins,'ro')
	show()
	print gfd
	    
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

	    
	    
