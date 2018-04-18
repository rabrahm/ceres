import pyfits
import argparse
import os
import spfr
import glob
from pylab import *
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('directorio')
parser.add_argument('-pars',default='6000,4.5,0.0,100')
parser.add_argument('-model_path',default='../../data/COELHO2014/')

args   = parser.parse_args()
dirin  = args.directorio
pars   = args.pars
model_path  = args.model_path

teff  = float(pars.split(',')[0])
logg  = float(pars.split(',')[1])
feh   = float(pars.split(',')[2])
vsini = float(pars.split(',')[3])

try:
	sc = pyfits.getdata(dirin)
	fits = [dirin]
	dirin = dirin.split('/')
	dirin = dirin[:-1]
	dt = ''
	for fn in dirin:
		dt += fn+'/'
	dirin = dt
except:
	fits = glob.glob(dirin + '/*fits')

os.system('mkdir '+dirin+'/FR')


for fit in fits:
	print 'RV computation of file ' + fit + ' whith model:' + pars

	sc = pyfits.getdata(fit)
	hd = pyfits.getheader(fit)

	I = np.where(sc[0,:,0]<5200)[0]
	npix = sc.shape[2]
	SNR_5130 = np.median(sc[8,I[0],int(0.5*npix)-100:int(0.5*npix)+100] )
	if SNR_5130 < 1.:
		SNR_5130 = 1.
	pix = (sc[0,0,0]/(sc[0,0,1]-sc[0,0,0]))/float(hd['RESOL'])
	SNR_5130_R = np.around(SNR_5130*np.sqrt(pix))

	p1gau,vels,xc_av,XCmodelgau = spfr.RVforFR(sc[0],sc[5], teff=teff, logg=logg, feh=feh, vsini=vsini, model_path=model_path,vstep=5.)

	A = 0.11081
	B = 0.0016
	D = 0.32815
	C = 0.00453

	RVerr =  B + ( 1.6 + 0.2 * p1gau[2] ) * A / np.round(SNR_5130)
	BSerr =  D / float(np.round(SNR_5130)) + C

	RVerr =  B + (1.6+0.2*p1gau[2]*0.5)*A/np.round(SNR_5130)
	depth_fact = 1. + p1gau[0]/(p1gau[2]*np.sqrt(2*np.pi))
	if depth_fact < 0.6:
		depth_fact = 0.6

	if depth_fact >= 1.:
		RVerr2 = -999.000
	else:
		depth_fact = (1 - 0.6) / (1 - depth_fact)
		RVerr2 = RVerr * depth_fact
 
 	RVerr2 = np.around(RVerr2,4)
	BSerr  = np.around(BSerr,4)

	print '\tRV = ',p1gau[1], '+/-', RVerr2, 'km/s'

	p1gau_m = p1gau
	XCmodel = XCmodelgau
	xc_dict = {'p1gau':p1gau, 'vels':vels, 'xc_av':xc_av,'XCmodelgau':XCmodelgau}
	ccf_pdf = dirin + '/FR/' + fit.split('/')[-1][:-4] + '_XC_FR.pdf'
	pkl_xc  = dirin + '/FR/' + fit.split('/')[-1][:-4] + '_XC_FR.pkl'
	pickle.dump( xc_dict, open( pkl_xc, 'w' ) )
	spfr.plot_CCF_FR(xc_dict,path=ccf_pdf)
	SP2 = spfr.calc_bss2(vels,xc_av,p1gau)

	print '\tBS = ',SP2, '+/-', BSerr,'km/s'

	line_out = "%-15s %18.8f %9.4f %7.4f %9.3f %5.3f   %s   ceres_FR   %8d %6d %5.2f %5.2f %5.1f %4.2f %5.2f %6.1f %4d %s\n"%\
               (hd['HIERARCH TARGET NAME'], 2400000.5 + float(hd['HIERARCH MBJD']), p1gau[1], RVerr2, SP2, BSerr, hd['INST'], int(hd['RESOL']), teff, logg, feh, vsini, xc_av.min(), p1gau[2],\
		       hd['TEXP (s)'], SNR_5130, ccf_pdf)

	i = 0
	isin = False
	if os.access(dirin+'/FR/results_FR.txt',os.F_OK):
		f = open(dirin+'/FR/results_FR.txt','r')
		lines = f.readlines()
		for line in lines:
			cos = line.split()
			if cos[1] == line_out.split()[1]:
				lines[i] = line_out
				isin = True
				break
			i+=1
		if not isin:
			lines.append(line_out)
	else:
		lines = [line_out]

	f = open(dirin+'/FR/results_FR.txt','w')
	for line in lines:
		f.write(line)
	f.close()





