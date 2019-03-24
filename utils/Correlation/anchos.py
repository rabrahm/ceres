from __future__ import print_function
import pyfits
import numpy
import scipy
import math
import os
import toVac
import prueba
from scipy import optimize
import matplotlib.pylab as plt
import readcol
def gauss1(params,x):
    C = params[0]
    A = params[1]
    med = params[2]
    sig = params[3]
    g = C+A*numpy.exp(-0.5*(x-med)*(x-med)/(sig*sig))
    return g

def res_gauss1(params,g,x):
    return g-gauss1(params,x)

def el_stl(wam,fl,SLi,SLf):
    for i in range(len(SLi)):
        if SLi[i]>wam[-1]:
            break
        I = numpy.where((wam >= SLi[i]) & (wam<=SLf[i]))[0]
        fl[I]=1.0
    return fl

SLi,SLf = readcol.readcol('lines2.dat',twod=False)
model_path = '/media/VERBATIM/COHELO_MODELS/RES_MOD/R_60000b/'
model_path = '/data/ajordan/COHELO_MODELS/R_60000b/'

non_rot = os.listdir(model_path+'vsini_0.0/')
rot = [0.0,2.5,5.0,7.5,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0]
i=0
f = open('anchos.dat','w')
verdad = False
non_rot = non_rot[200:]
for fits in non_rot:

    sc = pyfits.getdata(model_path+'vsini_0.0/'+fits)
    hd = pyfits.getheader(model_path+'vsini_0.0/'+fits)
    wa1 = numpy.arange(len(sc))*hd['CD1_1']+hd['CRVAL1']
    sc = el_stl(wa1,sc,SLi,SLf)
    wa = toVac.ToVacuum(wa1)
    I = numpy.where((wa > 5500) & (wa < 6250))[0]
    wa = wa[I]
    fl = sc[I]
    T = hd['TEFF']
    G = hd['LOG_G']
    Z = hd['FEH']
    devs = ''
    if T == 4500 and G == 2.5 and Z == -1.5:
        verdad=True
    #elif T == 4000 and G == 3.0:
    #       break
    if verdad:
        for v in rot:
            if os.access(model_path+'vsini_'+str(v)+'/'+'R_'+str(v)+'_'+fits[-22:],os.F_OK):
                sc2 = pyfits.getdata(model_path+'vsini_'+str(v)+'/'+'R_'+str(v)+'_'+fits[-22:])
                sc2 = el_stl(wa1,sc2,SLi,SLf)
                fl2 = sc2[I]
                vv,cc = prueba.CCF(wa,fl,wa,fl2,-200.0,200.0)
                cc2 = numpy.array(cc)
                vv2 = numpy.array(vv)
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

                devs = devs+str(sig3)+'\t'

            else:
                devs = devs+str(0.0)+'\t'
        f.write(str(T)+'\t'+str(G)+'\t'+str(Z)+'\t'+devs+'\n')

    print(T,G,Z,devs)
