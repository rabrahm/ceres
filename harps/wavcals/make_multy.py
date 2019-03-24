from __future__ import print_function
import numpy

#inp = 'B_order_00.dat'
ii = 10
while ii < 44:
    inp = 'B_order_'+str(ii)+'.dat'
    print(inp)
    f = open(inp,'r')
    lines = f.readlines()
    f.close()
    pix,wav, elm = [],[],[]
    for line in lines:
        cos = line.split()
        pix.append(int(cos[1]))
        wav.append(float(cos[2]))
        elm.append(cos[3])
    vecs = []
    vec = [0]
    i = 1
    while i < len(pix):
        if pix[i] - pix[i-1] <= 8:
            vec.append(i)
        else:
            vecs.append(vec)
            vec=[i]
        i+=1
    vecs.append(vec)
    f = open(inp[:11]+'iwdat','w')
    for vec in vecs:
        lout = str(len(vec))
        for i in vec:
            lout += '\t' + str(pix[i]) + '\t' + str(wav[i])
        for i in vec:
            lout += '\t' + elm[i]
        lout += '\n'
        f.write(lout)
    f.close()
    ii+=1
