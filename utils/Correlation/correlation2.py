#!/usr/bin/python
from pylab import *
base = '/data/echelle/'
import numpy
import scipy
import time
import os
import math
import pyfits
import vels
from scipy import optimize
from scipy import interpolate
from scipy import integrate
import copy
from pylab import *

def n_Edlen(l):
    """ 
    Refractive index according to Edlen 1966 

    """

    sigma = 1e4 / l
    sigma2 = sigma*sigma
    n = 1 + 1e-8 * (8342.13 + 2406030 / (130-sigma2) + 15997/(38.9-sigma2))

    return n

def n_Morton(l):
    """
    Refractive index used by JJ: Morton 1991, ApJS, 77, 119
    """
    sigma = 1e4 / l
    sigma2 = sigma*sigma
    n = 1 + 6.4328e-5 + 2.94981e-2 / (146.-sigma2) + 2.5540e-4/(41.-sigma2)

    return n

def ToAir(l):
    """ 

    Transforms input wavelength (A) from Vacuum to Air
    Uses refractive index as per Edlen 1966

    """

    return (l / n_Edlen(l))

def ToVacuum(l):
    """ 

    Transforms input wavelength (A) from Air to Vacuum
    Uses refractive index as per Edlen 1966


    """

    cond = 1
    l_prev = l.copy()
    while(cond):
        l_new = n_Edlen(l_prev) * l
        if (max(abs(l_new - l_prev)) < 1e-10): cond = 0
        l_prev = l_new
    
    return l_prev

def el_stl(wam,fl,SLi,SLf):
	for i in range(len(SLi)):
		if SLi[i]>wam[-1]:
			break
		I = numpy.where((wam >= SLi[i]) & (wam<=SLf[i]))[0]
		fl[I]=1.0
	return fl

def corr(L,F,mwa,FF,ons):
	CCF,NOR = 0,0
	
	o = 0
	while o < L.shape[0]:
		k = 0
		while k<len(ons):
			if ons[k]==o:
				o=o+1
			k=k+1
		
		I  = numpy.where(F[o]!=0.0)[0]
		if I.shape[0] > 0:
			Fl = F[o][I]
			Ll = L[o][I]
			
			I   = numpy.where( (mwa > Ll[0]) & (mwa< Ll[-1]) )[0]
			MML = mwa[I]
			MF  = FF[I]
			#ajustec = continuum.NORM_single(MML, MF, orden=2)
			#plot( MML, MF )
			#plot( MML, numpy.polyval(ajustec, MML))
			
			#NMF = MF / numpy.polyval(ajustec, MML)
			tckM = interpolate.splrep(MML,MF,k=3,s=0)
			NFM  = interpolate.splev(Ll,tckM,der=0)
			
			a    = integrate.simps(Fl[1:-1]*NFM[1:-1],Ll[1:-1])
			b    = integrate.simps(Fl[1:-1]*Fl[1:-1],Ll[1:-1])
			c    = integrate.simps(NFM[1:-1]*NFM[1:-1],Ll[1:-1])
				
			CCF  = CCF+(a/math.sqrt(b*c))*(Ll[-2]-Ll[1])
			NOR = NOR+Ll[-2]-Ll[1]
					
		o += 1
	#show()
	return CCF/NOR

def get_con_coef(MML,MDR,n):
	Io = len(MML)
	coef = numpy.polyfit(MML,MDR,n)
	res = MDR - numpy.polyval(coef,MML)
	dev = numpy.sqrt(numpy.var( res ))
	J = numpy.where( numpy.absolute(res) > 3*dev)[0]
	I = numpy.where( numpy.absolute(res) <= 3*dev)[0]
	cond = True
	if len(J) == 0 or len(I)< 0.3*Io:
		cond = False
	while cond:
		MML,MDR = MML[I],MDR[I]
		coef = numpy.polyfit(MML,MDR,n)
		res = MDR - numpy.polyval(coef,MML)
		dev = numpy.sqrt(numpy.var( res ))
		J = numpy.where( numpy.absolute(res) > 3*dev)[0]
		I = numpy.where( numpy.absolute(res) <= 3*dev)[0]
		if len(J) == 0 or len(I)< 0.3*Io:
			cond = False
	return coef


def get_chis(L,F,snr,mwa,FF,ons):
	res = numpy.array([])
	o = 0
	while o < L.shape[0]:
		k = 0
		while k<len(ons):
			if ons[k]==o:
				o=o+1
			k=k+1
		I  = numpy.where(F[o]!=0.0)[0]
		if I.shape[0] > 0:
			Fl = F[o][I]
			Ll = L[o][I]
			SN = snr[o][I]
			I   = numpy.where( (mwa > L[o][0]) & (mwa< L[o][-1]) )[0]
			MML = mwa[I]
			MF  = FF[I]
			#ajustec = continuum.NORM_single(MML, MF, orden=2)
			#plot( MML, MF )
			#plot( MML, numpy.polyval(ajustec, MML))
			
			#NMF = MF / numpy.polyval(ajustec, MML)
			tckM = interpolate.splrep(MML,MF,k=3,s=0)
			NFM  = interpolate.splev(Ll,tckM,der=0)

			MDR = NFM/Fl
			coef = get_con_coef(Ll,MDR,3)
			Fl = Fl*numpy.polyval(coef,Ll)
			#plot(Ll,NFM,linewidth=2.0)
			#plot(Ll,Fl)
			err = numpy.polyval(coef,Ll)/SN
			res = numpy.hstack((res,(Fl - NFM)/err))
		o+=1
	chis = numpy.add.reduce(res**2)/float(len(res)-5)
	
	#show()
		
	return chis
def bad_orders(L,bl):
	bo = []

	o = 0
	while o < L.shape[0]:
		i=0
		while i<len(bl):
			I = numpy.where((L[o]>bl[i][0]) & (L[o]<bl[i][1]))[0]
			if I.shape[0]!=0:
				bo.append(o)
			i = i+1
		o=o+1
	return bo

def get_name(t,g,z):
	T = str(int(t))
	
	if g == 0.0:
		G = '00'
	elif g == 0.5:
		G = '05'
	else:
		G = str(int(g*10))

	if z < 0:
		siZ = 'm'
		modz = -z
	else:
		siZ = 'p'
		modz = z

	if modz == 0.0:
		Z = '00'
	elif modz == 0.2:
		Z = '02'
	elif modz == 0.5:
		Z = '05'
	else:
		Z = str(int(modz*10))

	return T+'_'+G+'_'+siZ+Z+'p00.ms.fits'

def orden(V,W):
	j=0
	lar = len(V)
	
	while j<lar-1:
		p = V[j]
		pp = W[j]
		i = j+1
		lu = j
		while i < lar:
			if V[i]<p:
				p = V[i]
				pp = W[i]
				lu = i
			i=i+1
		V[lu]=V[j]
		W[lu]=W[j]
		V[j]=p
		W[j]=pp
		j=j+1
	return V,W

def gauss1(params,x):
	C = params[0]
	A = params[1]
	med = params[2]
	sig = params[3]
	g = C+A*numpy.exp(-0.5*(x-med)*(x-med)/(sig*sig))
	return g

def res_gauss1(params,g,x):
	return g-gauss1(params,x)

def CCF(spec, model_path='/dummy/path/',doplot = False, plot_dir = '/home/rabrahm/',plot_name = 'MY_LUP'):
	"""
	This function finds an aproximation to the stellar parameters (Teff, log(g), [Fe/H]) 
	of the input echelle spectrum using a CCF with model spectra. This code also is 
	constructed to find the radial velocity of the star and v*sin(i).
	"""
	
	lux = 299792.458

	L1,F1,SNR1 = spec[0,:,:], spec[5,:,:], spec[8,:,:]

	#width_path = '/home/rabrahm/Desktop/corr2/'
	#slines_path = '/home/rabrahm/Desktop/corr2/'
	width_path = '/data/echelle/ecpipe/Correlation/'
	slines_path = '/data/echelle/ecpipe/Correlation/'
	
	SLi,SLf = numpy.loadtxt(slines_path+'lines2.dat',dtype=None,unpack=True)
	
	SLi = ToVacuum(SLi)
	SLf = ToVacuum(SLf)
	
	AT,AG,AZ,A000,A025,A050,A075,A100,A150,A200,A250,A300,A350,A400,A450,A500 =\
            numpy.loadtxt(width_path+'anchos50000.dat',dtype=None,unpack=True)
	
	vsini = [0.0,2.5,5.0,7.5,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0]
	
	or01 = 0
	for i in range(L1.shape[0]):
		I = numpy.where(L1[i] < 6300.0)[0]
		if len(I) > 0:
			or01 = i
			break

	or02 = L1.shape[0]-1
	for i in range(L1.shape[0]):
		I = numpy.where(L1[i] < 4500.0)[0]
		if len(I) > 0:
			or02 = i
			break

	or03 = 0
	for i in range(L1.shape[0]):
		I = numpy.where(L1[i] < 6250.0)[0]
		if len(I) > 0:
			or03 = i
			break

	or04 = L1.shape[0]-1
	for i in range(L1.shape[0]):
		I = numpy.where(L1[i] < 5500.0)[0]
		if len(I) > 0:
			or04 = i
			break

	or05 = 0
	for i in range(L1.shape[0]):
		I = numpy.where(L1[i] < 5190.0)[0]
		W = numpy.where(L1[i] < 5170.0)[0]
		if len(I) > 0 and len(W) > 0:
			or05 = i
			break
	#print or01,or02,or03,or04,or05
       	guess = [1.0,1.0,1.0]
	
	L = L1[or01:or02]
	F = F1[or01:or02]
	SNR = SNR1[or01:or02]
	Lm = L1[or03:or04]
	Fm = F1[or03:or04]
	SNRm = SNR1[or03:or04]
	Lg = L1[or05]
	Fg = F1[or05]
	SNRg = SNR1[or05]
	
	bad_lines = [[6860,6900],[6550,6580],[6270,6320],[4850,4880]]
	ons = bad_orders(L,bad_lines)
	
	modi = 'R_0.0_5000_30_p00p00.ms.fits'
	#print 'Radial velocity calculation via CCF with: '+modi 
	
	sci = pyfits.getdata(model_path+'vsini_0.0/'+modi)
	hdi = pyfits.getheader(model_path+'vsini_0.0/'+modi)
	
	wam1 = ToVacuum(numpy.arange(len(sci))*hdi['CD1_1']+hdi['CRVAL1'])
	
	Im = numpy.where((wam1 > 5400.0) & (wam1 < 6350.0))[0]
	wam = wam1[Im]
	flm = sci[Im]

	Ig = numpy.where((wam1 > 5100.0) & (wam1 < 5250.0))[0]
	wag = wam1[Ig]

	It = numpy.where((wam1 > 4000.0) & (wam1 < 7500.0))[0]
	wat = wam1[It]

	for i in range(Lm.shape[0]):
		I = numpy.where((Fm[i] != 0.0) & (Fm[i] < 2.0))[0]
		Ls = Lm[i][I]
		Fs = Fm[i][I]
		
		I = numpy.where((wam > Ls[0]-5.0) & (wam < Ls[-1]+5.0))[0]
		MLs = wam[I]
		MFs = flm[I]
		
		vv,cc = vels.CCF(MLs,MFs,Ls,Fs,-200.0,200.0)

		if i == 0:
			ccf = numpy.array(cc)*(Ls[-1]-Ls[0])
			nor = Ls[-1]-Ls[0]
		else:
			ccf = ccf + numpy.array(cc)*(Ls[-1]-Ls[0])
			nor = nor + (Ls[-1]-Ls[0])
		

	ccf = ccf/nor
	vv = numpy.array(vv)

	B = 0.5*(ccf[0]+ccf[-1])
	A = numpy.max(ccf)-B
	med = vv[numpy.where(ccf == numpy.max(ccf))[0]]
	sig = 20.0
	guess1 = [B,A,med,sig]
	ajustep=optimize.leastsq(res_gauss1,guess1,args=(ccf,vv))
	velo = ajustep[0][2]

	#print 'The initial radial velocity is: '+str(velo)+' km/s'
	
	#print 'Determining parameters of the initial model'

	vecti = [3500,4000,4500,5000,5500,6500]
	vecgi = [0.0,1.5, 3.0,4.5]
	#veczi = [-2.0,-1.0,0.0]
	veczi = [-1.0,0.0]
	ccmax = 0
	names = []
	cc = 0
	nor = 0
	for g in vecgi:
		for z in veczi:
			for t in vecti:
				nam = get_name(t,g,z)
				names.append(nam)
	namemax = get_name(5500,4.5,0.0)
	TEI = 5500
	MEI = 0.0
	LGI = 4.5
	for nam in names:
		if ( (os.access(model_path+'vsini_0.0/R_0.0_'+nam,os.F_OK) == True)):
			mod = pyfits.getdata(model_path+'vsini_0.0/R_0.0_'+nam)
			flm = mod[Im]
			cc = 0.0
			nor = 0.0
			for o in range(Lm.shape[0]):
				I = numpy.where(Fm[i] != 0.0)[0]
				Ls = Lm[i][I]
				Fs = Fm[i][I]
				I = numpy.where((wam > Ls[0]-5.0) & (wam < Ls[-1]+5.0))[0]
				MLs = wam[I]*(1+velo/lux)
				MFs = flm[I]
					
				tck = interpolate.splrep(MLs,MFs,k=3,s=0)
				NMFs = interpolate.splev(Ls,tck,der=0)
					
				cc = cc + integrate.simps(NMFs[1:-1]*Fs[1:-1],Ls[1:-1])/math.sqrt(integrate.simps(Fs[1:-1]*Fs[1:-1],Ls[1:-1])*integrate.simps(NMFs[1:-1]*NMFs[1:-1],Ls[1:-1]))*(Ls[-1]-Ls[0])
				nor = nor + (Ls[-1]-Ls[0])
			cc = cc/nor		
			
			if cc >=ccmax:
				ccmax = cc
				namemax = nam

	mod = pyfits.getheader(model_path+'vsini_0.0/R_0.0_'+namemax)
	TEI = mod['TEFF']

	#print 'Teff (initial) = '+str(TEI)+' K'
	
	
	
	if TEI <= 4000:
		rot = 5.0
		LGI = 3.0
		MTI = 0.0
		late = True
                velo2 = velo

	else:
		late = False
		t = TEI
		vecgi  = [1.0,2.0,3.0,4.0]
		#veczi  = [-2.0,-1.0,0.0]
		veczi = [-1.0,0.0]
		dif = 1000.0

		for z in veczi:
			vals = []
			di = 1000
			meg = 3.0
			for g in vecgi:
				nam = get_name(t,g,z)
				if ( (os.access(model_path+'vsini_0.0/R_0.0_'+nam,os.F_OK) == True)):
					mod = pyfits.getdata(model_path+'vsini_0.0/R_0.0_'+nam)
					flm = mod[Im]
					flm = el_stl(wam,flm,SLi,SLf)
					intfl = 0.0
					intflm = 0.0
					di = 1000
					for o in range(Lm.shape[0]):
						I = numpy.where(Fm[i] != 0.0)[0]
						Ls = Lm[i][I]
						Fs = Fm[i][I]
						Fs = el_stl(Ls,Fs,SLi*(1+velo/lux),SLf*(1+velo/lux))
						I = numpy.where((wam > Ls[0]-5.0) & (wam < Ls[-1]+5.0))[0]
						MLs = wam[I]*(1+velo/lux)
						MFs = flm[I]
						tck = interpolate.splrep(MLs,MFs,k=3,s=0)
						NMFs = interpolate.splev(Ls,tck,der=0)

						intfl = intfl + integrate.simps(Fs[1:-1]*Fs[1:-1],Ls[1:-1])
						intflm = intflm + integrate.simps(NMFs[1:-1]*NMFs[1:-1],Ls[1:-1])
					intfl = math.sqrt(intfl)
					intflm = math.sqrt(intflm)	
					dif1 = math.sqrt((intflm-intfl)**2)
					vals.append(dif1)
					if dif1 < di:
						di = dif1
						meg = g
			dif2 = numpy.mean(vals)
			#print z,dif2
			if dif2<dif:
				dif =dif2
				LGI = meg
				MEI = z
				TEI = t
				
	
		#print '[Fe/H] (initial) = '+str(MEI)
	
		vecgi  = [0.5,1.5,2.5,3.5,4.5]
		cfi2 = 0.0
		for g in vecgi:
			nam = get_name(TEI,g,MEI)
			if ( (os.access(model_path+'vsini_0.0/R_0.0_'+nam,os.F_OK) == True)):
				mod = pyfits.getdata(model_path+'vsini_0.0/R_0.0_'+nam)
				
				mflg = mod[Ig]
				
				I = numpy.where(Fg != 0.0)[0]
				Ls = Lg[I]
				Fs = Fg[I]

				I = numpy.where((wag > Ls[0]-5.0) & (wag < Ls[-1]+5.0))[0]
				MLs = wag[I]*(1+velo/lux)
				MFs = mflg[I]
				tck = interpolate.splrep(MLs,MFs,k=3,s=0)
				NMFs = interpolate.splev(Ls,tck,der=0)
				
				cc2 = integrate.simps(Fs[1:-1]*NMFs[1:-1])/math.sqrt(integrate.simps(Fs[1:-1]*Fs[1:-1])*integrate.simps(NMFs[1:-1]*NMFs[1:-1]))
				#print g,cc2
				if cc2 > cfi2:
					cfi2 = cc2
					LGI  = g

		#print 'Log(g) (initial) = '+str(LGI)	
	
	itera = 0
	maximo = 0
	calculated = []
	rotss = []
	vecR = []
	vecT  = []
	vecG  = []
	
	vecZ  = []
	vecCF = []
	while itera < 4:
		
		if late==False:

			MOG = get_name(TEI,LGI,MEI)	
			sc = pyfits.getdata(model_path+'vsini_0.0/R_0.0_'+MOG)

			#print '-Calculating radial shift and v*sin(i) with model: '+MOG
	
			flm = sc[Im]
			ies = []
			flm = el_stl(wam,flm,SLi,SLf)
			for i in range(Lm.shape[0]):
				I = numpy.where(Fm[i] != 0.0)[0]
				Ls = Lm[i][I]
				Fs = Fm[i][I]
				Fs = el_stl(Ls,Fs,SLi*(1+velo/lux),SLf*(1+velo/lux))
				I = numpy.where((wam > Ls[0]-5.0) & (wam < Ls[-1]+5.0))[0]
				ies.append(I)
				MLs = wam[I]
				MFs = flm[I]
				vv,cc = vels.CCF(MLs,MFs,Ls,Fs,-200,200)
				if i == 0:
					ccf = numpy.array(cc)*(Ls[-1]-Ls[0])
					nor = Ls[-1]-Ls[0]
				else:
					ccf = ccf + numpy.array(cc)*(Ls[-1]-Ls[0])
					nor = nor + (Ls[-1]-Ls[0])
			nor2=1.
			ccf = ccf/nor
			vv = numpy.array(vv)
			B = 0.5*(ccf[0]+ccf[-1])
			A = numpy.max(ccf)-B
			med = vv[numpy.where(ccf == numpy.max(ccf))[0]]
			sig = 20.0
			guess1 = [B,A,med,sig]
			ajustep=optimize.leastsq(res_gauss1,guess1,args=(ccf,vv))
			velo2 = ajustep[0][2]	
			sig2 = ajustep[0][3]
			sig2 = math.sqrt(sig2*sig2)
			#plt.plot(vv,ccf,vv,gauss1(ajustep[0],vv))
			#plt.show()
			#print 'radial velocity = '+str(velo2)+' km/s'
			#print 'Sigma = '+str(sig2)+' km/s'
			"""			
				
			vi = 0.0
			difsigmin = 1000
			while vi <= 20.0:
				ai = 0
				if ( (os.access(model_path+'vsini_'+str(vi)+'/R_'+str(vi)+'_'+MOG,os.F_OK) == True)):
					modt = pyfits.getdata(model_path+'vsini_'+str(vi)+'/R_'+str(vi)+'_'+MOG)
					flm2 = modt[Im]
					flm2 = el_stl(wam,flm2,SLi,SLf)
					for I in ies:
						Fs = flm[I]
						Ls = wam[I]
						Fs = el_stl(Ls,Fs,SLi,SLf)
						MFs = flm2[I]
						vv2,cc2 = vels.CCF(Ls,MFs,Ls,Fs,-200,200)
						if ai == 0:
							ccf2 = numpy.array(cc2)*(Ls[-1]-Ls[0])
							nor2 = Ls[-1]-Ls[0]
						else:
							ccf2 = ccf2 + numpy.array(cc2)*(Ls[-1]-Ls[0])
							nor2 = nor2 + (Ls[-1]-Ls[0])
						ai += 1
					cc2 = ccf2/nor2
					vv2 = numpy.array(vv2)
					B3 = 0.5*(cc2[0]+cc2[-1])
					A3 = numpy.max(cc2)-B3
					med3 = 0.0
					sig3 = 20.0
					guess1 = [B3,A3,med3,sig3]
					ajustep=optimize.leastsq(res_gauss1,guess1,args=(cc2,numpy.array(vv2)))
					cte3 = ajustep[0][0]
					no3 = ajustep[0][1]
					med3 = ajustep[0][2]	
					sig3 = ajustep[0][3]
					
					#plt.plot(vv2,cc2)
					
					print vi,sig3
					difsig = math.sqrt((sig2-sig3)**2)
					if difsig < difsigmin:
						difsigmin = difsig
						rot = vi
				ai +=1
				
				if vi <= 7.5:
					vi = vi+2.5
				elif vi < 50.0:
					vi = vi+5.0
				else:
					break
					
			#plt.show()
			"""
			if sig2 >= vsini[0]:
				I = numpy.where((AT == TEI) & (AG == LGI) & (AZ == MEI))[0][0]
				anchos = numpy.array([A000[I],A025[I],A050[I],A075[I],A100[I],A150[I],A200[I],A250[I],A300[I],A350[I],A400[I],A450[I],A500[I]])
				kis = 0
				while kis < len(anchos)-1:
					if anchos[kis]>anchos[kis+1]:
					        portemp = anchos[kis]
						anchos[kis] = anchos[kis+1]
						anchos[kis+1] = portemp
					kis+=1
				tck = interpolate.splrep(anchos,numpy.array(vsini),k=3,s=0)
				calrot = interpolate.splev(sig2,tck,der=0)
				difs = (numpy.array(vsini) - calrot)**2
				AI = numpy.where(difs == numpy.min(difs))[0]
				rot = vsini[AI]
			else:
				rot = vsini[0]
				calrot = 0.0
			
			#rot = 2.5
			#print 'v*sin(i) = '+str(calrot)+' km/s'
	
		RI = numpy.where(numpy.array(rotss) == rot)[0]
		
		if len(RI) > 0:
			break

		else:
			
			model_path1 = model_path+"vsini_"+str(rot)+"/"
	
			""" conociendo la velocidad radial realizo una busqueda gruesa de parametros estelares calculando el maximo de la CCF"""
			#print "-Searching the optimal stellar model"

			#maxT,maxG,maxZ,maxCor,maxvs = 5500,4.0,0.0,0.0,0.0
			maxCor = 0.0
			minCHIS = 100000000000000000000.
			if TEI >= 3500 and TEI <=4250:
				modT = [3500,4000,4500,5000]
			else:
				modT = [4000,4500,5000,5500,6000,6500,7000]
				
			modG = [1.0,2.5,4.0]
			#modZ = [-2.0,-1.0,0.0]
			modZ = [-1.0,0.0]
			MOD = []

			i=0
			while i < len(modZ):
				ig = 0
				while ig < len(modG):
					it = 0			
					while it < len(modT):
						MK = 'R_'+str(rot)+'_'+get_name(modT[it],modG[ig],modZ[i])
						if os.access(model_path1+MK,os.F_OK):
							MOD.append(MK)
						it = it+1
					ig =ig+1
				i = i+1
			
			for m in MOD:
				
				hd = pyfits.getheader(model_path1+m)
				T  = hd['TEFF']
				G  = hd['LOG_G']
				Z  = hd['FEH']
				vs = hd['VSINI']

				sc = pyfits.getdata(model_path1+m)
				FF = sc[It]
				mwa = wat*(1+velo2/lux)	

				#NCCF = corr(L,F,mwa,FF,ons)
				CHIS = get_chis(L,F,SNR,mwa,FF,ons)

				NCCF = 1./CHIS
				if NCCF > maxCor:
					maxCor = NCCF
					maxG   = G
					maxT   = T
					maxZ   = Z
					maxvs  = vs
				elif NCCF<0:
					#print "Problem with spectrum!!! -> Negative value of CCF."
					maxG   = 4.5
					maxT   = 5500
					maxZ   = 0
					maxvs  = 0
				
				

			#print 'maxgrueso',maxT,maxG,maxZ,maxvs

			""" A partir de los parametros encontrados en la busqueda gruesa procedo a constrenir los limites de la busqueda fina"""


			#Tfin,Gfin,Zfin,rotfin,velfin = 5500,4.0,0.0,0.0,0.0
			""" Ahora se buscan los par'ametros estelares optimos mediante una exploracion fina"""

			if maxT == 3500 or maxT == 4000:
				modT = [3500,3750,4000,4250,4500]

			elif maxT == 7000 or maxT == 6500:
				modT = [6000,6250,6500,6750,7000]
			else:
				modT = [maxT-750,maxT-500,maxT-250,maxT,maxT+250,maxT+500,maxT+750]
	
			if maxG == 1.0:
				modG = [0.0,0.5,1.0,1.5,2.0,2.5]
			if maxG == 2.5:
				modG = [1.0,1.5,2.0,2.5,3.0,3.5,4.0]
			if maxG == 4.0:
				modG = [2.5,3.0,3.5,4.0,4.5,5.0]

			if maxZ == -2.0:
				modZ = [-2.5,-2.0,-1.5,-1.0]
			if maxZ == -1.0:
				modZ = [-1.5,-1.0,-0.5,0.0]
			if maxZ == 0.0:
				modZ = [-0.5,0.0,0.2,0.5]

			MOD = []

			for i in modZ:
				for ig in modG:	
					for it in modT:
						MK = 'R_'+str(rot)+'_'+get_name(it,ig,i)
						if os.access(model_path1+MK,os.F_OK):
							MOD.append(MK)
	
			#maxT,maxG,maxZ,maxCor = 5500,4.0,0,0
			maxCor = 0.0
			minCHIS = 10000000000000000.
			for m in MOD:
				calculated.append(m)
				hd = pyfits.getheader(model_path1+m)
				T  = hd['TEFF']
				G  = hd['LOG_G']
				Z  = hd['FEH']

				sc = pyfits.getdata(model_path1+m)
				FF = sc[It]
				mwa = wat*(1+velo2/lux)
		
				#NCCF = corr(L,F,mwa,FF,ons)
				CHIS = get_chis(L,F,SNR,mwa,FF,ons)
				NCCF = 1./CHIS
				
				vecT.append(T)
				vecG.append(G)
				vecZ.append(Z)
				vecR.append(rot)
				vecCF.append(NCCF)
				
				if NCCF > maxCor:
					maxCor = NCCF
					maxT   = T
					maxG   = G
					maxZ   = Z

				elif NCCF < 0:
					print "Problem with spectrum!!! -> Negative value of CCF."
					maxG   = 4.5
					maxT   = 5500
					maxZ   = 0
					maxvs  = 0
			
			#print "maxfino",maxT,maxG,maxZ,maxCor
			#print "maxfino",minT,minG,minZ,maxCor
			
			TEI = maxT
			LGI = maxG
			MEI = maxZ

			ultrot = rot

			if maxCor > maximo:
				maximo = maxCor
				Tfin,Gfin,Zfin,rotfin,velfin = maxT,maxG,maxZ,rot,velo2
				#Tf,Gf,Zf = intert,interg,maZ

			late = False
			rotss.append(rot)
			
			
		itera = itera+1
	if maximo == 0:
		Tfin,Gfin,Zfin,rotfin,velfin = 5500, 4.5,0,0,velo2
	#print 'Pars fase 0:', Tfin, Gfin, Zfin,rotfin,maximo
	mejor = False
	if rotfin == 0.0:
		nrot = [0.0,2.5]
	elif rotfin == 2.5:
		nrot = [0.0,2.5,5.0]
	elif rotfin == 5.0:
		nrot = [2.5,5.0,7.5]
	elif rotfin == 7.5:
		nrot = [5.0,7.5,10.0]
	else:
		nrot = [rotfin-5.0,rotfin,rotfin+5.0]

	if Tfin == 3500:
		nT = [3500,3750,4000]
	elif Tfin == 7000:
		nT = [6500,6750,7000]
	else:
		nT = [Tfin-250,Tfin,Tfin+250]

	if Gfin == 0.0:
		nG = [0.0,0.5,1.0]
	elif Gfin == 0.5:
		nG = [0.0,0.5,1.0]
	elif Gfin == 4.5:
		nG = [4.0,4.5,5.0]
	elif Gfin == 5.0:
		nG = [4.0,4.5,5.0]
	else:
		nG = [Gfin-0.5,Gfin,Gfin+0.5]

	if Zfin == -2.5:
		nZ = [-2.5,-2.0,-1.5]
	elif Zfin == 0.5:
		nZ = [0.0,0.2,0.5]
	elif Zfin == 0.2:
		nZ = [0.0,0.2,0.5]
	elif Zfin == 0.0:
		nZ = [-0.5,0.0,0.2,0.5]
	else:
		nZ = [Zfin-0.5,Zfin,Zfin+0.5]

	for v in nrot:
		model_path2 = model_path+'vsini_'+str(v)+'/'
		names = []
		calc = numpy.array(calculated)

		for t in nT:
			for g in nG:
				for z in nZ:
					nam = 'R_'+str(v)+'_'+get_name(t,g,z)
					I = numpy.where(calc == nam)[0]
					if len(I)==0 and os.access(model_path2+nam,os.F_OK):
						names.append(nam)

		for fits in names:
			calculated.append(fits)
			hd = pyfits.getheader(model_path2+fits)
			T  = hd['TEFF']
			G  = hd['LOG_G']
			Z  = hd['FEH']

			sc = pyfits.getdata(model_path2+fits)
			FF = sc[It]
			mwa = wat*(1+velfin/lux)
		
			NCCF = corr(L,F,mwa,FF,ons)
			CHIS = get_chis(L,F,SNR,mwa,FF,ons)
			#print T,G,Z,vs,CHIS	
			vecT.append(T)
			vecG.append(G)
			vecZ.append(Z)
			vecCF.append(NCCF)
			vecR.append(v)
			#print T,G,Z,NCCF

			if NCCF > maximo:
				mejor = True
				maximo = NCCF
				Tfin   = T
				Gfin   = G
				Zfin   = Z
				rotfin = v

		#print Tfin,Gfin,Zfin,rotfin,maximo
		maximCCF = maximo

	total_inter = 1
	if total_inter == 0:		
		deltaV = 0.5
		deltaG = 0.05
		deltaT = 50.0
		deltaZ = 0.1

		vCF = numpy.array(vecCF)
		vvT = numpy.array(vecT)
		vvG = numpy.array(vecG)
		vvZ = numpy.array(vecZ)
		vvV = numpy.array(vecR)
			
		ejeT = numpy.arange(numpy.min(nT),numpy.max(nT)+deltaT,deltaT)
		ejeV = numpy.arange(numpy.min(nrot),numpy.max(nrot)+deltaV,deltaV)
		ejeG = numpy.arange(numpy.min(nG),numpy.max(nG)+deltaG,deltaG)
		ejeZ = numpy.arange(numpy.min(nZ),numpy.max(nZ)+deltaZ,deltaZ)

		lejeV = len(ejeV)
		lejeG = len(ejeG)
		lejeT = len(ejeT)
		lejeZ = len(ejeZ)
			
		matCCF = numpy.zeros([lejeT,lejeG,lejeZ,lejeV],float)

		for v in nrot:
			posV = int(round((v-ejeV[0])/deltaV))
			for z in nZ:
				posZ = int(round((z-ejeZ[0])/deltaZ))
				for g in nG:
					posG = int(round((g-ejeG[0])/deltaG))
					I = numpy.where((vvV == v) & (vvG ==g) & (vvZ == z))[0]
					if len(I) > 0:
						vTt,vCFt = orden(vvT[I],vCF[I])
				
						if len(I)>3:
							tck  = interpolate.splrep(vTt,vCFt,k=3,s=0)
							ynew = interpolate.splev(ejeT,tck,der=0)
						elif len(I) > 1:
							tck  = interpolate.splrep(vTt,vCFt,k=len(I)-1,s=0)
							ynew = interpolate.splev(ejeT,tck,der=0)
						else:
							ynew = numpy.zeros(len(ejeT),float)+vCFt[0]
		
						matCCF[:,posG,posZ,posV]=ynew
	
		for v in nrot:
			posV = int(round((v-ejeV[0])/deltaV))
			for z in nZ:
				posZ = int(round((z-ejeZ[0])/deltaZ))
				for t in ejeT:
					posT = int(round((t-ejeT[0])/deltaT))
					y1 = matCCF[posT,:,posZ,posV]
					I = numpy.where(y1 != 0.0)[0]
					y1b = y1[I]
					x1b = ejeG[I]
					if len(I) > 0:
						if len(I)>3:
							tck  = interpolate.splrep(x1b,y1b,k=3,s=0)
							ynew = interpolate.splev(ejeG,tck,der=0)
						elif len(I) > 1:
							tck  = interpolate.splrep(x1b,y1b,k=len(I)-1,s=0)
							ynew = interpolate.splev(ejeG,tck,der=0)
						else:
							ynew = numpy.zeros(len(ejeG),float)+y1b[0]
						matCCF[posT,:,posZ,posV]=ynew
		for v in nrot:
			posV = int(round((v-ejeV[0])/deltaV))
			for t in ejeT:
				posT = int(round((t-ejeT[0])/deltaT))
				for g in ejeG:
					posG = int(round((g-ejeG[0])/deltaG))
					y1 = matCCF[posT,posG,:,posV]
					I = numpy.where(y1 != 0.0)[0]
					y1b = y1[I]
					x1b = ejeZ[I]
					if len(I) > 0:
						if len(I)>3:
							tck  = interpolate.splrep(x1b,y1b,k=3,s=0)
							ynew = interpolate.splev(ejeZ,tck,der=0)
						elif len(I) > 1:
							tck  = interpolate.splrep(x1b,y1b,k=len(I)-1,s=0)
							ynew = interpolate.splev(ejeZ,tck,der=0)
						else:
							ynew = numpy.zeros(len(ejeZ),float)+y1b[0]
						matCCF[posT,posG,:,posV]=ynew

		for t in ejeT:
			posT = int(round((t-ejeT[0])/deltaT))
			for g in ejeG:
				posG = int(round((g-ejeG[0])/deltaG))
				for z in ejeZ:
					posZ = int(round((z-ejeZ[0])/deltaZ))
					y1 = matCCF[posT,posG,posZ,:]
					I = numpy.where(y1 != 0.0)[0]
					y1b = y1[I]
					x1b = ejeV[I]
					if len(I) > 0:
						if len(I)>3:
							tck  = interpolate.splrep(x1b,y1b,k=3,s=0)
							ynew = interpolate.splev(ejeV,tck,der=0)
						elif len(I) > 1:
							tck  = interpolate.splrep(x1b,y1b,k=len(I)-1,s=0)
							ynew = interpolate.splev(ejeV,tck,der=0)
						else:
							ynew = numpy.zeros(len(ejeV),float)+y1b[0]
						matCCF[posT,posG,posZ,:]=ynew
	
		I = numpy.where(matCCF == numpy.max(matCCF))

		intert = ejeT[I[0]]		
		interg = ejeG[I[1]]
		interz = ejeZ[I[2]]
		interv = ejeV[I[3]]
		maximCCF = numpy.max(CCF)
		#print 'interp',intert,interg,interz,interv

	elif total_inter == 1:
		if Tfin == 3500:
			nT = [3500,3750,4000,4250]
		elif Tfin == 3750:
			nT = [3500,3750,4000,4250]
		elif Tfin == 7000:
			nT = [6250,6500,6750,7000]
		elif Tfin == 6750:
			nT = [6250,6500,6750,7000]
		else:
			nT = [Tfin-500,Tfin-250,Tfin,Tfin+250,Tfin+500]

		if Gfin == 0.0:
			nG = [0.0,0.5,1.0,1.5]
		elif Gfin == 0.5:
			nG = [0.0,0.5,1.0,1.5]
		elif Gfin == 4.5:
			nG = [3.5,4.0,4.5,5.0]
		elif Gfin == 5.0:
			nG = [3.5,4.0,4.5,5.0]
		else:
			nG = [Gfin-1.0,Gfin-0.5,Gfin,Gfin+0.5,Gfin+1.0]
		calc = numpy.array(calculated)
		model_path2 = model_path+'vsini_'+str(rotfin)+'/'
		names = []

		for t in nT:
			for g in nG:
				nam = 'R_'+str(rotfin)+'_'+get_name(t,g,Zfin)
				I = numpy.where(calc == nam)[0]

				if len(I) == 0 and os.access(model_path2+nam,os.F_OK):
					names.append(nam)

		for fits in names:
			calculated.append(fits)
			hd = pyfits.getheader(model_path2+fits)
			T  = hd['TEFF']
			G  = hd['LOG_G']
			Z  = hd['FEH']

			sc = pyfits.getdata(model_path2+fits)
			FF = sc[It]
			mwa = wat*(1+velfin/lux)
		
			NCCF = corr(L,F,mwa,FF,ons)
			CHIS = get_chis(L,F,SNR,mwa,FF,ons)
			#print T,G,Z,vs,CHIS
			vecZ.append(Z)	
			vecT.append(T)
			vecG.append(G)
			vecR.append(rotfin)
			vecCF.append(NCCF)
			

		VZ = numpy.array(vecZ)
		VR = numpy.array(vecR)
		VT = numpy.array(vecT)
		VG = numpy.array(vecG)
		VF = numpy.array(vecCF)

		I = numpy.where((VZ == Zfin) & (VR == rotfin))[0]
		VT2 = VT[I]
		VG2 = VG[I]
		VF2 = VF[I]
		deltaT = 50.0
		deltaG = 0.05
		
		ejeT = numpy.arange(numpy.min(numpy.array(nT)),numpy.max(numpy.array(nT))+deltaT,deltaT)
		ejeG = numpy.arange(numpy.min(numpy.array(nG)),numpy.max(numpy.array(nG))+deltaG,deltaG)

		lejeT = len(ejeT)
		lejeG = len(ejeG)
			
		matCCF = numpy.zeros([lejeT,lejeG],float)

		for g in nG:
			pos = int(round((g-ejeG[0])/deltaG))
			I = numpy.where(VG2 == g)[0]
			if len(I) > 0:
				
				vTt,vCFt = orden(VT2[I],VF2[I])
				#print vTt,vCFt
				if len(I)>3:
					tck  = interpolate.splrep(vTt,vCFt,k=3,s=0)
					ynew = interpolate.splev(ejeT,tck,der=0)
				elif len(I) > 1:
					tck  = interpolate.splrep(vTt,vCFt,k=len(I)-1,s=0)
					ynew = interpolate.splev(ejeT,tck,der=0)
				else:
					ynew = numpy.zeros(len(ejeT),float)+vCFt[0]
				
				matCCF[:,pos]=ynew
			
		for t in ejeT:
			pos1 = int(round((t-ejeT[0])/deltaT))
			y1 = matCCF[pos1,:]
			I = numpy.where(y1 != 0.0)[0]
			y1b = y1[I]
			x1b = ejeG[I]
			if len(I) > 0:
				if len(I)>3:
					tck  = interpolate.splrep(x1b,y1b,k=3,s=0)
					ynew = interpolate.splev(ejeG,tck,der=0)
				elif len(I) > 1:
					tck  = interpolate.splrep(x1b,y1b,k=len(I)-1,s=0)
					ynew = interpolate.splev(ejeG,tck,der=0)
				else:
					ynew = numpy.zeros(len(ejeG),float)+y1b[0]
			matCCF[pos1,:]=ynew
			
		I = numpy.where(matCCF == numpy.max(matCCF))
			
		intert = ejeT[I[0]][0]
		interg = ejeG[I[1]][0]
		interz = Zfin
		interv = rotfin
		#print intert,interg
		maximCCF = numpy.max(matCCF)

	else:
		intert,interg,interz,interv = Tfin,Gfin,Zfin,rotfin
	#print 'HI'
	if doplot:
		nam = 'R_'+str(rotfin)+'_'+get_name(Tfin,Gfin,Zfin)
		hd = pyfits.getheader(model_path2+nam)
		sc = pyfits.getdata(model_path2+nam)
		PWAV = ToVacuum(numpy.arange(len(sc))*hd['CD1_1']+hd['CRVAL1'])
		PWAV = PWAV*(1+velfin/lux)
		for i in range(L1.shape[0]):
			I = numpy.where((PWAV > L1[i,0]) & (PWAV < L1[i,-1]))[0]
			#print L1[i]
			#print F1[i]
			#print i
			f = figure()
			ylim(0.,1.01)
			plot(L1[i],F1[i])
			plot(PWAV[I],sc[I])
			xlabel('wavelenth [A]')
			ylabel('Continuum Normalized Flux')
			savefig(plot_dir + plot_name + '_' + str(int(L1[i,0])) + '_' + str(int(L1[i,-1])) + '.pdf', format='pdf')
			
			
	return [intert, interg, interz, interv, velfin,maximCCF]
