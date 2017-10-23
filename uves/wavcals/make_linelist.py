from pylab import *

refo = 3
while refo < 23:
	srefo = str(refo)
	if refo <10:
		srefo = '0'+srefo

	d = np.loadtxt('all_lines.txt')
	allw = d[:,1]

	d2 = np.loadtxt('rorder_'+srefo+'.iwdat')
	p,w = d2[:,0],d2[:,1]
	c = np.polyfit(w,p,4)
	npx = np.polyval(c,allw)

	I = np.where((npx>10) & (npx<4038))[0]
	npx,nwv = npx[I],allw[I]

	f = open('nrorder_'+srefo+'.iwdat','w')
	for i in range(len(npx)):
		if nwv[i] > w[0] - 100 and nwv[i] < w[-1] + 100:
			#print npx[i],nwv[i]
			f.write('1\t'+str(npx[i])+'\t'+str(nwv[i])+'\tNN\n')
	f.close()
	print p - np.polyval(c,w)
	refo+=1
	#plot(w,p,'ro')
	#plot(nwv,npx,'b')
	#show()