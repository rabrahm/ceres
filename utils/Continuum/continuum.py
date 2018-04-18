import sys
import math
import pyfits
import numpy
import scipy
from scipy import optimize
import FunNorm
import time
import os
import matplotlib.pyplot as plt
import gc
import warnings

def division(L,F,n):

	x = int(len(L)/n)
	vec = numpy.zeros((2,x),float)

	i = 0
	while i < x:

		if i < x-1:
			vec[0][i] = numpy.mean(L[i*n:(i+1)*n])
			vec[1][i] = numpy.median(F[i*n:(i+1)*n])
		elif i == x-1:
			vec[0][i] = numpy.mean(L[i*n:len(L)])
			vec[1][i] = numpy.median(F[i*n:len(F)])
		
		i = i+1
	
	g = [-0.1,1.0]
	if vec.shape[1] > 2:
		ajuste = optimize.leastsq(res_lin,g,args=(vec[1][:],vec[0][:]))

	else:
		ajuste = optimize.leastsq(res_lin,g,args=(F,L))

	return ajuste

def lin(params,x):
	m = params[0]
	n = params[1]
	y = m*x+n
	return y

def res_lin(params,y,x):
	return y-lin(params,x)

def curva(params,x):
	a_0=params[0]
	a_1=params[1]
	a_2=params[2]
	a_3=params[3]
	a_4=params[4]
	
	y=a_0+a_1*x+a_2*x*x+a_3*(x**3)+a_4*(x**4)
	return y

def res_cur(params,y,x):
	return y-curva(params,x)

def NORM(L,F,orden=2):

	ordenes = L.shape[0]
	pixeles = L.shape[1]

	guess = [-0.1,1.0]
	guess2 = [1,0.1,0.1,0.1,0.1]
	
	div = 4
	rec = 20

	para = numpy.zeros( (2,div),float )
	continuo = numpy.zeros( (2,ordenes),list )
	strl = numpy.zeros( (6,2),float )
	laro = numpy.zeros( ordenes,int )
	
	#strli,strlf = numpy.loadtxt('/home/rabrahm/Desktop/corr2/strong_lines.dat',dtype=None,usecols=(0,1),unpack=True)
	strli,strlf = numpy.loadtxt('../utils/Continuum/strong_lines.dat',dtype=None,usecols=(0,1),unpack=True)

	esto = numpy.zeros( (2,ordenes),list )
	ajent = numpy.zeros( (ordenes,pixeles),float )
	unos = numpy.zeros(pixeles,float)
	NORMALIZADO = numpy.zeros(ordenes,list)
	unos = unos+1.0
	
	Lfilter = numpy.zeros(ordenes,list)
	Ffilter = numpy.zeros(ordenes,list)
	coefs = numpy.zeros([ordenes,orden+1])
	
	o=0
	while o<ordenes:

		Lmien = numpy.zeros(pixeles,float)
		Fmien = numpy.zeros(pixeles,float)

		"""
		elemino flujos nulos
		"""

		I = numpy.where( F[o] != 0.0 )[0]
		Llimpio = L[o][I]
		Flimpio = F[o][I]
		
		"""
		hago filtro de 5 pixeles
		"""
	
		i = 0
		k = 0

		largo_nul = Llimpio.shape[0]
	
		while i < largo_nul:

			if i+9 < largo_nul:
				Lmien[k] = numpy.mean(Llimpio[i:i+5])
				Fmien[k] = numpy.median(Flimpio[i:i+5])
				i = i+5
				k = k+1
			else:
				Lmien[k] = numpy.mean(Llimpio[i:largo_nul])
				Fmien[k] = numpy.median(Flimpio[i:largo_nul])
				i = i+9
				k = k+1

		largo_fil = k	
		Llimpio = Lmien[0:largo_fil]
		Flimpio = Fmien[0:largo_fil]
		
		Lfilter[o] = list(Llimpio)
		Ffilter[o] = list(Flimpio)
	
		esto[0,o] = list(Llimpio)
		esto[1,o] = list(Flimpio)

		"""
		no considero strong lines:
		"""
		rellenado = False
		
		j = 0
		while j<strli.shape[0]:
			
			I = numpy.where( (Llimpio > strli[j]) & (Llimpio < strlf[j]) )[0]
			Flimpio[I]=-1

			if I.shape[0] > 0:
				rellenado == True

			j=j+1
	
		"""
		completo flujo de strong lines con ajuste lineal 
		"""
		if rellenado:
			Flimpio = FunNorm.Rell(Llimpio.astype("double"),Flimpio.astype("double"),Flimpio.shape[0])

		"""
		divido cada orden en 3 y ajusto una recta
		"""
	
		can = int(Llimpio.shape[0]/div)
		mm = numpy.zeros(div,float)
		nn = numpy.zeros(div,float)

		c=0
		while c < div:

			vectorL = Llimpio[c*can:(c+1)*can]
			vectorF = Flimpio[c*can:(c+1)*can]
			para = division(vectorL,vectorF,rec)
			
			mm[c] = para[0][0]
			nn[c] = para[0][1]
			c = c+1

		c = 0
		FF = []
		LL = []

		while c < div:

			m = mm[c]
			n = nn[c]

			if c != div-1:
				segl = Llimpio[c*can:(c+1)*can]
				segf = segl*m+n
			elif c == div-1:
				segl = Llimpio[c*can:Llimpio.shape[0]]
				segf = segl*m+n

			LL = LL+list(segl)
			FF = FF+list(segf)
			c = c+1

		LL = numpy.array(LL)
		FF = numpy.array(FF)

		ajustep = numpy.polyfit(LL,FF,orden)
		FF = numpy.polyfit(ajustep, LL)
		"""
		ajustep=optimize.leastsq(res_cur,guess2,args=(FF,LL))
		a0=ajustep[0][0]	
		a1=ajustep[0][1]
		a2=ajustep[0][2]
		a3=ajustep[0][3]
		a4=ajustep[0][4]
	
		FF=a0+a1*LL+a2*(LL**2)+a3*(LL**3)+a4*(LL**4)
		"""
		"""
		considero los peaks sobre 1 sigma y ajusto polinomio de orden 4
		"""

		JJ = FF.copy()
		RES = Flimpio-FF
		veva = numpy.zeros(pixeles,float)
		AJ = numpy.zeros((Flimpio.shape[0]),float)
		
		I = numpy.where(RES > 0.0)[0]
		veva = RES[I]

		dev = math.sqrt(numpy.mean(veva*veva))
		
		i=0
		while i<FF.shape[0]:

			if RES[i] > dev and RES[i] <= 4*dev:
				AJ[i] = Flimpio[i]
			elif RES[i] > 4*dev:
				AJ[i] = -1
				Flimpio[i] = FF[i]
			else:
				AJ[i] = -1

			i=i+1

		JJ = FF+dev
		
		if AJ[0]== -1:
			AJ[0] = numpy.median(Flimpio[0:20])
		if AJ[-1] == -1:
			AJ[-1] = numpy.median(Flimpio[-21:-1])
	
		FFF = FunNorm.Rell(Llimpio,AJ,Flimpio.shape[0])
		LLL = Llimpio.copy()

		ajustep = numpy.polyfit(LLL,FFF,orden)
		CON = numpy.polyfit(ajustep, LLL)
		"""
		ajustep = optimize.leastsq(res_cur,guess2,args=(FFF,LLL))
		a0 = ajustep[0][0]	
		a1 = ajustep[0][1]
		a2 = ajustep[0][2]
		a3 = ajustep[0][3]
		a4 = ajustep[0][4]
	
		CON = a0+a1*LLL+a2*(LLL**2)+a3*(LLL**3)+a4*(LLL**4)
		"""
		continuo[0,o] = list(LLL)
		continuo[1,o] = list(CON)

		laro[o] = LLL.shape[0]

		esto[0,o] = list(Llimpio)
		esto[1,o] = list(Flimpio)

		o=o+1

	"""
	itero para subir el continuo
	"""

	#print "iterando para subir el continuo"

	j = 0
	while j < 20:

		o = 0
		while o < ordenes:
			
			ya = numpy.array(continuo[1,o])
			RESI = numpy.array(esto[1,o])-ya
			larg = laro[o]
			
			loco=numpy.array(esto[0,o])
			foco=numpy.array(esto[1,o])

			I = numpy.where(RESI > 0)[0]
			
			ya[I] = foco[I]
			

			ajustec	 = numpy.polyfit(loco,ya)
			aji = numpy.polyval(ajustec,loco[0:larg])
			invor = numpy.arange(orden+1)[::-1]
			for orde in invor:
				coefs[o,orde] = invor[orde]
 			"""
			ajustec = optimize.leastsq(res_cur,guess2,args=(ya,loco))
			A0 = ajustec[0][0]
			A1 = ajustec[0][1]	
			A2 = ajustec[0][2]	
			A3 = ajustec[0][3]	
			A4 = ajustec[0][4]	
			coefs[o,0]=A0
			coefs[o,1]=A1
			coefs[o,2]=A2
			coefs[o,3]=A3
			coefs[o,4]=A4
			aji = A0+A1*(loco[0:larg]**1)+A2*(loco[0:larg]**2)+A3*(loco[0:larg]**3)+A4*(loco[0:larg]**4)
			"""
			continuo[1,o] = list(aji)
			
			o=o+1
		j=j+1	
	


	"""
	perfecciono bordes del ajuste con pendiente de orden anterior y posterior
	"""

	#print "arreglando bordes de cada orden"

	o = 0
	while o < ordenes-1:
		
		larg = laro[o]
		longitud = numpy.array(continuo[0,o])
		upL = numpy.array(continuo[0,o+1])
		upF = numpy.array(continuo[1,o+1])
		
		lami = longitud[0]
		lamf = lami+5
	
		I = numpy.where( (upL >= lami) & (upL <= lamf))[0]
		
		if len(I)>4:
			
			veL = upL[I]
			veF = upF[I]
			
			ajuster = optimize.leastsq(res_lin,guess,args=(veF,veL))
			pen = ajuster[0][0]

			j = 0
			
			while j < 5:

				flujo = numpy.array(continuo[1,o])
				I = numpy.where( longitud >= lamf)[0][0]
				x = longitud[I]
				y = flujo[I]
			
				coe = y-pen*x

				flujo[0:I] = longitud[0:I]*pen+coe
				
				ajustec = numpy.polyfit(longitud,flujo,orden)
				tri = numpy.polyval(ajustec, longitud)
		
				"""
				ajustec = optimize.leastsq(res_cur,guess2,args=(flujo,longitud))
				A0 = ajustec[0][0]
				A1 = ajustec[0][1]	
				A2 = ajustec[0][2]	
				A3 = ajustec[0][3]	
				A4 = ajustec[0][4]	
	
				tri=A0+A1*longitud+A2*(longitud**2)+A3*(longitud**3)+A4*(longitud**4)
				"""
				continuo[1,o] = list(tri)

				j = j+1
			
			
			if o == 0:
			
				AJUSTE_FINAL = A0+A1*(L[o]**1)+A2*(L[o]**2)+A3*(L[o]**3)+A4*(L[o]**4)
				NORMALIZADO[o] = list(F[o]/AJUSTE_FINAL)

		else:
			AJUSTE_FINAL = numpy.polyval(coefs[o],L[o])
			#AJUSTE_FINAL = coefs[o,0]+coefs[o,1]*(L[o]**1)+coefs[o,2]*(L[o]**2)+coefs[o,3]*(L[o]**3)+coefs[o,4]*(L[o]**4)
			NORMALIZADO[o] = list(F[o]/AJUSTE_FINAL)
			
		o=o+1	

	o = ordenes-1
	while o > 0:
		
		longitud=numpy.array(continuo[0,o])
		downL=numpy.array(continuo[0,o-1])
		downF=numpy.array(continuo[1,o-1])
		larg=laro[o]
		lamf=longitud[larg-1]
		lami=lamf-5
		
		I = numpy.where((downL>=lami) & (downL<=lamf))[0]
		if len(I)>4:
			veL = downL[I]
			veF = downF[I]
		
			ajuster = optimize.leastsq(res_lin,guess,args=(veF,veL))
			pen = ajuster[0][0]

			j = 0
			while j < 3:

				flujo = numpy.array(continuo[1,o])
				I = numpy.where(longitud<=lami)[0][-1]
				x = longitud[I]
				y = flujo[I]
			
				coe = y-pen*x
			
				flujo[I:] = longitud[I:]*pen+coe
				ajustec = numpy.polyfit(longitud,flujo,orden)
				tri = numpy.polyval(ajustec, longitud)
				"""
				ajustec = optimize.leastsq(res_cur,guess2,args=(flujo,longitud))
				A0 = ajustec[0][0]
				A1 = ajustec[0][1]	
				A2 = ajustec[0][2]	
				A3 = ajustec[0][3]	
				A4 = ajustec[0][4]	
	
				tri=A0+A1*(longitud[:]**1)+A2*(longitud[:]**2)+A3*(longitud[:]**3)+A4*(longitud[:]**4)
				"""
				continuo[1,o]=list(tri)

				j=j+1
			AJUSTE_FINAL = numpy.polyval( ajustec,L[o] )
			#AJUSTE_FINAL = A0+A1*(L[o]**1)+A2*(L[o]**2)+A3*(L[o]**3)+A4*(L[o]**4)
		
			NORMALIZADO[o] = list(F[o]/AJUSTE_FINAL)
		
		else:
			AJUSTE_FINAL = numpy.polyval(coefs[o],L[o])
			#AJUSTE_FINAL = coefs[o,0]+coefs[o,1]*(L[o]**1)+coefs[o,2]*(L[o]**2)+coefs[o,3]*(L[o]**3)+coefs[o,4]*(L[o]**4)
			NORMALIZADO[o] = list(F[o]/AJUSTE_FINAL)
		o=o-1

	AAAA = numpy.zeros(ordenes,list)
	BBBB = numpy.zeros(ordenes,list)

	o = 0
	while o < ordenes:
		
		AAAA[o] = numpy.array(NORMALIZADO[o])	
		BBBB[o]=L[o]
		
		o = o+1;
	
	return BBBB,AAAA;

def NORM2(L,F):

	ordenes = L.shape[0]
	pixeles = L.shape[1]

	guess = [-0.1,1.0]
	guess2 = [1,0.1,0.1,0.1,0.1]
	
	div = 4
	rec = 20

	para = numpy.zeros( (2,div),float )
	continuo = numpy.zeros( (2,ordenes),list )
	strl = numpy.zeros( (6,2),float )
	laro = numpy.zeros( ordenes,int )
	
	#strli,strlf = numpy.loadtxt('/home/rabrahm/Desktop/corr2/strong_lines.dat',dtype=None,usecols=(0,1),unpack=True)
	strli,strlf = numpy.loadtxt('../utils/Continuum/strong_lines.dat',dtype=None,usecols=(0,1),unpack=True)

	esto = numpy.zeros( (2,ordenes),list )
	ajent = numpy.zeros( (ordenes,pixeles),float )
	unos = numpy.zeros(pixeles,float)
	NORMALIZADO = numpy.zeros(ordenes,list)
	unos = unos+1.0
	
	Lfilter = numpy.zeros(ordenes,list)
	Ffilter = numpy.zeros(ordenes,list)
	coefs = numpy.zeros([ordenes,5])
	
	o=0
	while o<ordenes:

		Lmien = numpy.zeros(pixeles,float)
		Fmien = numpy.zeros(pixeles,float)

		"""
		elemino flujos nulos
		"""

		I = numpy.where( F[o] != 0.0 )[0]
		Llimpio = L[o][I]
		Flimpio = F[o][I]
		
		"""
		hago filtro de 5 pixeles
		"""
	
		i = 0
		k = 0

		largo_nul = Llimpio.shape[0]
	
		while i < largo_nul:

			if i+9 < largo_nul:
				Lmien[k] = numpy.mean(Llimpio[i:i+5])
				Fmien[k] = numpy.median(Flimpio[i:i+5])
				i = i+5
				k = k+1
			else:
				Lmien[k] = numpy.mean(Llimpio[i:largo_nul])
				Fmien[k] = numpy.median(Flimpio[i:largo_nul])
				i = i+9
				k = k+1

		largo_fil = k	
		Llimpio = Lmien[0:largo_fil]
		Flimpio = Fmien[0:largo_fil]
		
		Lfilter[o] = list(Llimpio)
		Ffilter[o] = list(Flimpio)
	
		esto[0,o] = list(Llimpio)
		esto[1,o] = list(Flimpio)

		"""
		no considero strong lines:
		"""
		rellenado = False
		
		j = 0
		while j<strli.shape[0]:
			
			I = numpy.where( (Llimpio > strli[j]) & (Llimpio < strlf[j]) )[0]
			Flimpio[I]=-1

			if I.shape[0] > 0:
				rellenado == True

			j=j+1
	
		"""
		completo flujo de strong lines con ajuste lineal 
		"""
		if rellenado:
			Flimpio = FunNorm.Rell(Llimpio.astype("double"),Flimpio.astype("double"),Flimpio.shape[0])

		"""
		divido cada orden en 3 y ajusto una recta
		"""
	
		can = int(Llimpio.shape[0]/div)
		mm = numpy.zeros(div,float)
		nn = numpy.zeros(div,float)

		c=0
		while c < div:

			vectorL = Llimpio[c*can:(c+1)*can]
			vectorF = Flimpio[c*can:(c+1)*can]
			para = division(vectorL,vectorF,rec)
			
			mm[c] = para[0][0]
			nn[c] = para[0][1]
			c = c+1

		c = 0
		FF = []
		LL = []
		
		while c < div:

			m = mm[c]
			n = nn[c]

			if c != div-1:
				segl = Llimpio[c*can:(c+1)*can]
				segf = segl*m+n
			elif c == div-1:
				segl = Llimpio[c*can:Llimpio.shape[0]]
				segf = segl*m+n

			LL = LL+list(segl)
			FF = FF+list(segf)
			c = c+1

		LL = numpy.array(LL)
		FF = numpy.array(FF)

		ajustep=optimize.leastsq(res_cur,guess2,args=(FF,LL))
		a0=ajustep[0][0]	
		a1=ajustep[0][1]
		a2=ajustep[0][2]
		a3=ajustep[0][3]
		a4=ajustep[0][4]
	
		FF=a0+a1*LL+a2*(LL**2)+a3*(LL**3)+a4*(LL**4)

		"""
		considero los peaks sobre 1 sigma y ajusto polinomio de orden 4
		"""

		JJ = FF.copy()
		RES = Flimpio-FF
		veva = numpy.zeros(pixeles,float)
		AJ = numpy.zeros((Flimpio.shape[0]),float)
		
		I = numpy.where(RES > 0.0)[0]
		veva = RES[I]

		dev = math.sqrt(numpy.mean(veva*veva))
		
		i=0
		while i<FF.shape[0]:

			if RES[i] > dev and RES[i] <= 4*dev:
				AJ[i] = Flimpio[i]
			elif RES[i] > 4*dev:
				AJ[i] = -1
				Flimpio[i] = FF[i]
			else:
				AJ[i] = -1

			i=i+1

		JJ = FF+dev
		
		if AJ[0]== -1:
			AJ[0] = numpy.median(Flimpio[0:20])
		if AJ[-1] == -1:
			AJ[-1] = numpy.median(Flimpio[-21:-1])
	
		FFF = FunNorm.Rell(Llimpio,AJ,Flimpio.shape[0])
		LLL = Llimpio.copy()

		ajustep = optimize.leastsq(res_cur,guess2,args=(FFF,LLL))
		a0 = ajustep[0][0]	
		a1 = ajustep[0][1]
		a2 = ajustep[0][2]
		a3 = ajustep[0][3]
		a4 = ajustep[0][4]
	
		CON = a0+a1*LLL+a2*(LLL**2)+a3*(LLL**3)+a4*(LLL**4)

		continuo[0,o] = list(LLL)
		continuo[1,o] = list(CON)

		laro[o] = LLL.shape[0]

		esto[0,o] = list(Llimpio)
		esto[1,o] = list(Flimpio)

		o=o+1

	"""
	itero para subir el continuo
	"""

	#print "iterando para subir el continuo"

	j = 0
	while j < 20:

		o = 0
		while o < ordenes:
			
			ya = numpy.array(continuo[1,o])
			RESI = numpy.array(esto[1,o])-ya
			larg = laro[o]
			
			loco=numpy.array(esto[0,o])
			foco=numpy.array(esto[1,o])

			I = numpy.where(RESI > 0)[0]
			
			ya[I] = foco[I]
									
			ajustec = optimize.leastsq(res_cur,guess2,args=(ya,loco))
			A0 = ajustec[0][0]
			A1 = ajustec[0][1]	
			A2 = ajustec[0][2]	
			A3 = ajustec[0][3]	
			A4 = ajustec[0][4]	
			coefs[o,0]=A0
			coefs[o,1]=A1
			coefs[o,2]=A2
			coefs[o,3]=A3
			coefs[o,4]=A4
			aji = A0+A1*(loco[0:larg]**1)+A2*(loco[0:larg]**2)+A3*(loco[0:larg]**3)+A4*(loco[0:larg]**4)
		
			continuo[1,o] = list(aji)
			
			o=o+1
		j=j+1	
	


	"""
	perfecciono bordes del ajuste con pendiente de orden anterior y posterior
	"""

	#print "arreglando bordes de cada orden"

	o = 0
	while o < ordenes-1:
		
		larg = laro[o]
		longitud = numpy.array(continuo[0,o])
		upL = numpy.array(continuo[0,o+1])
		upF = numpy.array(continuo[1,o+1])
		
		lami = longitud[0]
		lamf = lami+5
	
		I = numpy.where( (upL >= lami) & (upL <= lamf))[0]
		
		if len(I)>4:
			
			veL = upL[I]
			veF = upF[I]
			
			ajuster = optimize.leastsq(res_lin,guess,args=(veF,veL))
			pen = ajuster[0][0]

			j = 0
			
			while j < 5:

				flujo = numpy.array(continuo[1,o])
				I = numpy.where( longitud >= lamf)[0][0]
				x = longitud[I]
				y = flujo[I]
			
				coe = y-pen*x

				flujo[0:I] = longitud[0:I]*pen+coe
					
				ajustec = optimize.leastsq(res_cur,guess2,args=(flujo,longitud))
				A0 = ajustec[0][0]
				A1 = ajustec[0][1]	
				A2 = ajustec[0][2]	
				A3 = ajustec[0][3]	
				A4 = ajustec[0][4]	
	
				tri=A0+A1*longitud+A2*(longitud**2)+A3*(longitud**3)+A4*(longitud**4)
	
				continuo[1,o] = list(tri)

				j = j+1
			
			
			if o == 0:
			
				AJUSTE_FINAL = A0+A1*(L[o]**1)+A2*(L[o]**2)+A3*(L[o]**3)+A4*(L[o]**4)
				NORMALIZADO[o] = list(F[o]/AJUSTE_FINAL)

		else:
			AJUSTE_FINAL = coefs[o,0]+coefs[o,1]*(L[o]**1)+coefs[o,2]*(L[o]**2)+coefs[o,3]*(L[o]**3)+coefs[o,4]*(L[o]**4)
			NORMALIZADO[o] = list(F[o]/AJUSTE_FINAL)
			
		o=o+1	

	o = ordenes-1
	while o > 0:
		
		longitud=numpy.array(continuo[0,o])
		downL=numpy.array(continuo[0,o-1])
		downF=numpy.array(continuo[1,o-1])
		larg=laro[o]
		lamf=longitud[larg-1]
		lami=lamf-5
		
		I = numpy.where((downL>=lami) & (downL<=lamf))[0]
		if len(I)>4:
			veL = downL[I]
			veF = downF[I]
		
			ajuster = optimize.leastsq(res_lin,guess,args=(veF,veL))
			pen = ajuster[0][0]

			j = 0
			while j < 3:

				flujo = numpy.array(continuo[1,o])
				I = numpy.where(longitud<=lami)[0][-1]
				x = longitud[I]
				y = flujo[I]
			
				coe = y-pen*x
			
				flujo[I:] = longitud[I:]*pen+coe

				ajustec = optimize.leastsq(res_cur,guess2,args=(flujo,longitud))
				A0 = ajustec[0][0]
				A1 = ajustec[0][1]	
				A2 = ajustec[0][2]	
				A3 = ajustec[0][3]	
				A4 = ajustec[0][4]	
	
				tri=A0+A1*(longitud[:]**1)+A2*(longitud[:]**2)+A3*(longitud[:]**3)+A4*(longitud[:]**4)
				continuo[1,o]=list(tri)

				j=j+1
		
			AJUSTE_FINAL = A0+A1*(L[o]**1)+A2*(L[o]**2)+A3*(L[o]**3)+A4*(L[o]**4)
		
			NORMALIZADO[o] = list(F[o]/AJUSTE_FINAL)
		
		else:
			AJUSTE_FINAL = coefs[o,0]+coefs[o,1]*(L[o]**1)+coefs[o,2]*(L[o]**2)+coefs[o,3]*(L[o]**3)+coefs[o,4]*(L[o]**4)
			NORMALIZADO[o] = list(F[o]/AJUSTE_FINAL)
		o=o-1

	AAAA = numpy.zeros(ordenes,list)
	BBBB = numpy.zeros(ordenes,list)

	o = 0
	while o < ordenes:
		
		AAAA[o] = numpy.array(NORMALIZADO[o])	
		BBBB[o]=L[o]
		
		o = o+1;
	
	return BBBB,AAAA;

def NORM_single(L, F, orden = 5):

	warnings.simplefilter('ignore', numpy.RankWarning)

	pixeles = L.shape[0]
	guess = [-0.1,1.0]
	guess2 = [1,0.1,0.1,0.1,0.1]
	
	div = 4
	rec = 20
	nnn = orden
	para = numpy.zeros( (2,div),float )
	#strli,strlf = numpy.loadtxt('/data/echelle/ecpipe/Continuum/strong_lines.dat',dtype=None,usecols=(0,1),unpack=True)
	strli = numpy.array([3965.0,3930.0])
	strlf = numpy.array([3975.0,3940.0])
	unos = numpy.zeros(pixeles,float) + 1.0
	
	Lmien = numpy.zeros(pixeles,float)
	Fmien = numpy.zeros(pixeles,float)

	"""
	elemino flujos nulos
	"""

	I = numpy.where( F != 0.0)[0]
	Llimpio = L[I]
	Flimpio = F[I]

	"""
	calculo mediana con 5 pixeles
	"""
	
	i = 0
	k = 0

	largo_nul = Llimpio.shape[0]
	
	while i < largo_nul:
		if i+9 < largo_nul:
			Lmien[k] = numpy.mean(Llimpio[i:i+5])
			Fmien[k] = numpy.median(Flimpio[i:i+5])
			i = i+5
			k = k+1
		else:
			Lmien[k] = numpy.mean(Llimpio[i:largo_nul])
			Fmien[k] = numpy.median(Flimpio[i:largo_nul])
			i = i+9
			k = k+1

	largo_fil = k	
	Llimpio = Lmien[0:largo_fil]
	Flimpio = Fmien[0:largo_fil]
	
	"""
	no considero strong lines:
	"""

	rellenado = False

	j = 0
	while j<strli.shape[0]:
			
		I = numpy.where( (Llimpio > strli[j]) & (Llimpio < strlf[j]) )[0]
		Flimpio[I]=-1

		if I.shape[0]>0:
			rellenado == True

		j=j+1
	
	"""
	completo flujo de strong lines con ajuste lineal 
	"""
	if rellenado:
		#print '!!!!!!!!!!!!!!!!!!!'
		Flimpio = FunNorm.Rell(Llimpio.astype("double"),Flimpio.astype("double"),Flimpio.shape[0])
		#print sys.getrefcount(Flimpio)

	"""
	divido cada orden en 3 y ajusto una recta
	"""
	
	can = int(Llimpio.shape[0]/div)
	mm = numpy.zeros(div,float)
	nn = numpy.zeros(div,float)

	c=0
	while c < div:

		vectorL = Llimpio[c*can:(c+1)*can+1]
		vectorF = Flimpio[c*can:(c+1)*can+1]
		para = division(vectorL,vectorF,rec)
			
		mm[c] = para[0][0]
		nn[c] = para[0][1]
		c = c+1

	c = 0
	FF = []
	LL = []
	
	while c < div:

		m = mm[c]
		n = nn[c]

		if c != div-1:
			segl = Llimpio[c*can:(c+1)*can]
			segf = segl*m+n
		elif c == div-1:
			segl = Llimpio[c*can:Llimpio.shape[0]]
			segf = segl*m+n

		LL = LL+list(segl)
		FF = FF+list(segf)
		c = c+1

	LL = numpy.array(LL)
	FF = numpy.array(FF)
	ajustep = numpy.polyfit(LL,FF,nnn)
	FF = numpy.polyval(ajustep,LL)
	"""
	ajustep=optimize.leastsq(res_cur,guess2,args=(FF,LL))
	a0=ajustep[0][0]	
	a1=ajustep[0][1]
	a2=ajustep[0][2]
	a3=ajustep[0][3]
	a4=ajustep[0][4]
	
	FF=a0+a1*LL+a2*(LL**2)+a3*(LL**3)+a4*(LL**4)
	"""
	"""
	considero los peaks sobre 1 sigma y ajusto polinomio de orden 4
	"""

	JJ = FF.copy()
	RES = Flimpio-FF
	veva = numpy.zeros(pixeles,float)
	AJ = numpy.zeros((Flimpio.shape[0]),float)
	
	WI = numpy.where(RES < 0.0)[0]
	devm = 	math.sqrt(numpy.mean(RES[WI]*RES[WI]))
	
	I = numpy.where(RES > 0.0)[0]
	veva = RES[I]

	dev = math.sqrt(numpy.mean(veva*veva))
	
	i=0
	while i<FF.shape[0]:

		if RES[i] > -devm and RES[i] <= 4*dev:
			AJ[i] = Flimpio[i]
		elif RES[i] > 4*dev:
			AJ[i] = -1
			Flimpio[i] = FF[i]
		else:
			AJ[i] = -1

		i=i+1

	JJ = FF+dev
		
	if AJ[0 ]== -1:
		AJ[0] = numpy.median(Flimpio[0:20])
	if AJ[-1] == -1:
		AJ[-1] = numpy.median(Flimpio[-21:-1])
	
	FFF = FunNorm.Rell(Llimpio,AJ,Flimpio.shape[0])

	LLL = Llimpio.copy()
	ajustep = numpy.polyfit(LLL,FFF,nnn)
	CON = numpy.polyval(ajustep,LLL)
	"""
	ajustep = optimize.leastsq(res_cur,guess2,args=(FFF,LLL))
	a0 = ajustep[0][0]	
	a1 = ajustep[0][1]
	a2 = ajustep[0][2]
	a3 = ajustep[0][3]
	a4 = ajustep[0][4]
	
	CON = a0+a1*LLL+a2*(LLL**2)+a3*(LLL**3)+a4*(LLL**4)
	"""
	"""
	itero para subir el continuo
	"""
	
	#print "iterando para subir el continuo"

	j = 0
	while j < 10:

		
		RESI = Flimpio-CON
		larg = Flimpio.shape[0]
			
		loco = Llimpio.copy()
		foco = Flimpio.copy()
		
		I = numpy.where(RESI > 0.0)
		CON[I] = foco[I]
		
		ajustec = numpy.polyfit(loco,CON,nnn)	
		aji = numpy.polyval(ajustec,loco)
		
		"""					
		ajustec = optimize.leastsq(res_cur,guess2,args=(CON,loco))
		A0 = ajustec[0][0]
		A1 = ajustec[0][1]	
		A2 = ajustec[0][2]	
		A3 = ajustec[0][3]	
		A4 = ajustec[0][4]	
		aji = A0+A1*(loco[0:larg]**1)+A2*(loco[0:larg]**2)+A3*(loco[0:larg]**3)+A4*(loco[0:larg]**4)
		"""

		CON = aji.copy()

		
		j=j+1
	
	I = numpy.where( (Llimpio > Llimpio[0]+5.0) & (Llimpio < Llimpio[0]+10.0) )[0]
	if len(I) > 2:
		ajl = numpy.polyfit(Llimpio[I],aji[I],1)
		pendi = ajl[0]
		"""
		ajl = optimize.leastsq(res_lin,guess,args=(aji[I],Llimpio[I]))
		pendi = ajl[0][0]
		"""
		coefi = aji[I[0]]-pendi*Llimpio[I[0]]	
		aj2 = aji.copy()
		aj2[0:I[0]] =  Llimpio[0:I[0]]*pendi + coefi
	else:
		aj2 = aji.copy()
	
	I = numpy.where( (Llimpio > Llimpio[-1]-10.0) & (Llimpio < Llimpio[-1]-5.0) )[0]
	if len(I) > 2:
		ajl = numpy.polyfit(Llimpio[I],aji[I],1)
		pendi = ajl[0]
		"""
		ajl = optimize.leastsq(res_lin,guess,args=(aji[I],Llimpio[I]))
		pendi = ajl[0][0]
		"""
		coefi = aji[I[-1]]-pendi*Llimpio[I[-1]]
		aj2[I[-1]:] =  Llimpio[I[-1]:]*pendi + coefi
	else:
		aj2 = aji.copy()
	ajustec = numpy.polyfit(Llimpio,aj2,nnn)
	ajf = numpy.polyval(ajustec,L)
	"""
	ajustec = optimize.leastsq(res_cur,guess2,args=(aj2,Llimpio))
	A0 = ajustec[0][0]
	A1 = ajustec[0][1]	
	A2 = ajustec[0][2]	
	A3 = ajustec[0][3]	
	A4 = ajustec[0][4]

	ajf = A0+A1*(L**1)+A2*(L**2)+A3*(L**3)+A4*(L**4)
	"""
	NF = F/ajf
	
	gc.collect()
	#print sys.getrefcount(para)
	#print sys.getrefcount(unos)
	#print sys.getrefcount(Lmien)
	return ajustec