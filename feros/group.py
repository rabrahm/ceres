from __future__ import print_function
import os
from pylab import *
f = open('/data/echelle/feros/k2_list.txt','r')
lines = f.readlines()
rv,rve,jd,refs = [],[],[],[]
for line in lines:
    path = line[:-1]
    #line = 'python ferospipe.py '+path+' -npools 16 -do_class -avoid_plot -o2do HD157347 -reffile /data/echelle/feros/reffile.txt'
    line = 'python ferospipe.py '+path+' -npools 16 -do_class'
    os.system(line)
    """
    pr = path[:-1]+'_red/proc/results.txt'

    try:
            fo = open(pr,'r')
            lines2 = fo.readlines()
            for line2 in lines2:
                    cos = line2.split()
                    if float(cos[16])>50 and cos[0]=='HD157347':
                            print line2[:-1]
                            rv.append(float(cos[2]))
                            rve.append(float(cos[3]))
                            jd.append(float(cos[1]))
                            st = path[:-1]
                            st = st.split('/')[-1]
                            mo = float(st[4:6])
                            da = float(st[6:])
                            ref = mo + da/30.5
                            refs.append(ref)
    except:
            print 'bad'



jd,rv,refs = np.array(jd),np.array(rv),np.array(refs)
I = np.where(refs<20)[0]
print np.sqrt(np.var(rv[I]))
errorbar(refs,rv,yerr=rve,fmt='ro')
show()
"""
