import pyfits
from pylab import *
import scipy
from scipy import interpolate
import sys
base = '../'
sys.path.append(base+"utils/GLOBALutils")
import GLOBALutils

lux = 2.99792458E5
scnd = np.array([  2.14285714e-07,   7.00000000e-04,  -3.04357143e+00])

sc = pyfits.getdata('/Users/rabrahm/data/espadons/20150107_red/1772293c.spec.ob.fits.S')
f = open('/Users/rabrahm/codes/ceres/espadons/thar_list.txt','r')
lines = f.readlines()
tlines = []
names  = []
for line in lines:
    cos = line.split()
    tlines.append(float(cos[0]))
    try:
        names.append(cos[1])
    except:
        names.append('?')
tlines = np.array(tlines)
names = np.array(names)
coefs = np.loadtxt('/Users/rabrahm/codes/ceres/espadons/wcal_ref.dat')
vels = np.arange(-500,500,10.)

nspec = np.ones((11,sc.shape[0],sc.shape[2]))
nspec[5] = sc[:,1,:]

for i in range(nspec.shape[1]):
    npix = np.arange(nspec.shape[2])
    wav  = np.polyval(coefs[i][1:][::-1],npix)*10
    nspec[0,i,:]  = wav
    spl           = interpolate.splrep(np.arange(wav.shape[0]), wav,k=3)
    dlambda_dx    = interpolate.splev(np.arange(wav.shape[0]), spl, der=1)
    NN            = np.average(dlambda_dx)
    dlambda_dx    /= NN
    nspec[9,i,:]  = nspec[5,i,:] * (dlambda_dx ** 1)
    nspec[8,i,:]  = 100.

ml_v = tlines - 0.01
mh_v = tlines + 0.01
weight = np.ones(len(tlines))
vels, xc_full, sn, nlines_ccf, W_ccf = \
                    GLOBALutils.XCor(nspec, ml_v, mh_v, weight, 0, 1., vel_width=300,vel_step=3.,\
                                          spec_order=9,iv_order=10,sn_order=8,max_vel_rough=300)
xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3., Simple=True, W=W_ccf)
#plot(vels,xc_av)
#show()
j = np.median(xc_av)
im = np.argmax(xc_av)
rv = vels[im]

vels, xc_full, sn, nlines_ccf, W_ccf = \
                    GLOBALutils.XCor(nspec, ml_v, mh_v, weight, rv, 1., vel_width=20,vel_step=0.3,\
                                          spec_order=9,iv_order=10,sn_order=8,max_vel_rough=300)
xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=3., Simple=True, W=W_ccf)
cent = np.sum(xc_av*vels) / np.sum(xc_av)
print cent
#axvline(cent)
#plot(vels,xc_av-j)
#show()
#print gfvcdxs
for i in range(nspec.shape[1]):
    si = str(i)
    if i<10:
        si = '0'+si
    f = open('order_'+si+'.iwdat','w')
    npix = np.arange(nspec.shape[2])
    wav  = np.polyval(coefs[i][1:][::-1],npix)*10
    wav  = wav*(1-cent/lux)
    tck = interpolate.splrep(npix,wav,k=3)
    wav = interpolate.splev(npix+np.polyval(scnd,npix),tck)

    I = np.where((tlines>wav[0]) & (tlines<wav[-1]))[0]

    tck = interpolate.splrep(wav,npix,k=3)
    lns = tlines[I]
    nms = names[I]
    plns = interpolate.splev(lns,tck)
    I = np.where((plns>50)&(plns<len(wav)-50))[0]
    lns = lns[I]
    plns = plns[I]
    nms = nms[I]
    #print plns
    #print lns
    if i ==30:
        plot(npix,sc[i,1])
        plot(plns,sc[i,1,plns.astype('int')],'ro')
        show()
    toy = True
    grp1,grp2,grp3 = [lns[0]],[plns[0]],[nms[0]]
    o=1
    while o < len(lns):
        if not toy:
            grp1,grp2,grp3 = [lns[o]],[plns[o]],[nms[o]]
            toy=True
        else:
            if plns[o]-plns[o-1]<3:
                grp1.append(lns[o])
                grp2.append(plns[o])
                grp3.append(nms[o])
            else:
                toy=False
                line = str(len(grp1))+'\t'
                for u in range(len(grp1)):
                    line += str(np.around(grp2[u]))+'\t'+str(grp1[u])+'\t'
                for u in range(len(grp3)):
                    line += str(grp3[u])+'\t'
                line += '\n'
                f.write(line)
                o-=1
        o+=1
    f.close()
