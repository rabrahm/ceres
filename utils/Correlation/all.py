from __future__ import print_function
import pyfits
import numpy
import os
import continuum
import correlation

dire = '/home/rabrahm/Desktop/spec2/'
mods = os.listdir(dire)
t=[]
g = []
z=[]
r=[]
v=[]
nam = []
mods = mods[:-5]
for fits in mods:
    print(fits)
    sc = pyfits.getdata(dire+fits)
    L,F = sc[0],sc[1]
    Ln,Fn = continuum.NORM2(L,F)
    Tfin, Gfin, Zfin, rot, velo2 = correlation.CCF(Ln,Fn)
    t.append(Tfin)
    g.append(Gfin)
    z.append(Zfin)
    r.append(rot)
    v.append(velo2)
    nam.append(fits)

cant = len(t)

o=0
while o< cant:
    print(nam[o],t[o],g[o],z[o],r[o],v[o])
    o+=1
