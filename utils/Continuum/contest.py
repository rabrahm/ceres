import sys
import numpy
#import FunNorm
import gc
from scipy import signal
import scipy

def divis(Lx,Fx,nx):


    lt = int(len(Lx)/float(nx))
    vecA,vecB = [],[]

    for i in range(lt):

        if i < lt-1:
            vecA.append(numpy.mean(Lx[int(i*nx):int((i+1)*nx)]))
            vecB.append(numpy.median(Fx[int(i*nx):int((i+1)*nx)]))
        else:
            vecA.append(numpy.mean(Lx[int(i*nx):]))
            vecB.append(numpy.median(Fx[int(i*nx):]))
    vecA = numpy.array(vecA)
    vecB = numpy.array(vecB)

    if len(vecA) > 2:
        ajuste = numpy.polyfit(vecA,vecB,1)
    else:
        ajuste = numpy.polyval(Lx,Fx,1)

    vecA = None
    vecB = None
    Lx = None
    Fx = None

    return ajuste


def NORM_single(L, F, orden = 1):
    ajustep = numpy.polyfit(L,F,orden)


    pixeles = L.shape[0]
    guess = [-0.1,1.0]
    guess2 = [1,0.1,0.1,0.1,0.1]

    div = 4
    rec = 20
    nnn = orden
    para = numpy.zeros( (2,div),float )
    #strli,strlf = numpy.loadtxt('/data/echelle/ecpipe/Continuum/strong_lines.dat',dtype=None,usecols=(0,1),unpack=True)
    strli = numpy.array([3965.0,3930.0,5164.0])
    strlf = numpy.array([3975.0,3940.0,5191.0])

    #Lmien = numpy.zeros(int(pixeles/5.),float)
    #Fmien = numpy.zeros(int(pixeles/5.),float)

    """
    elemino flujos nulos
    """

    I = numpy.where( F != 0.0)[0]
    Llimpio = L[I]
    Flimpio1 = F[I]

    """
    calculo mediana con 5 pixeles
    """
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

    print k, Lmien.shape
    Llimpio = Lmien
    Flimpio = Fmien
    """
    Flimpio = scipy.signal.medfilt(Flimpio1,5)



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
        print '!!!!!!!!!!!!!!!!!!!'
        Flimpio = FunNorm.Rell(Llimpio.astype("double"),Flimpio.astype("double"),Flimpio.shape[0])
        print sys.getrefcount(Flimpio)

    """
    divido cada orden en 3 y ajusto una recta
    """

    can = int(Llimpio.shape[0]/float(div))
    mm,nn=[],[]
    for c in range(div):
        para = divis( Llimpio[int(c*can):int((c+1)*can)], Flimpio[int(c*can):int((c+1)*can)], rec )
        #para = [1,1]
        mm.append(para[0])
        nn.append(para[1])
        para = None
    mm = numpy.array(mm)
    nn = numpy.array(nn)

    """
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

    JJ = FF.copy()
    RES = Flimpio-FF
    veva = numpy.zeros(pixeles,float)
    AJ = numpy.zeros((Flimpio.shape[0]),float)

    WI = numpy.where(RES < 0.0)[0]
    devm =  numpy.sqrt(numpy.mean(RES[WI]*RES[WI]))

    I = numpy.where(RES > 0.0)[0]
    veva = RES[I]

    dev = numpy.sqrt(numpy.mean(veva*veva))

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
    """
    itero para subir el continuo
    """
    """
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

            CON = aji.copy()


            j=j+1
    I = numpy.where( (Llimpio > Llimpio[0]+5.0) & (Llimpio < Llimpio[0]+10.0) )[0]
    if len(I) > 2:
            ajl = numpy.polyfit(Llimpio[I],aji[I],1)
            pendi = ajl[0]
            coefi = aji[I[0]]-pendi*Llimpio[I[0]]
            aj2 = aji.copy()
            aj2[0:I[0]] =  Llimpio[0:I[0]]*pendi + coefi
    else:
            aj2 = aji.copy()

    I = numpy.where( (Llimpio > Llimpio[-1]-10.0) & (Llimpio < Llimpio[-1]-5.0) )[0]
    if len(I) > 2:
            ajl = numpy.polyfit(Llimpio[I],aji[I],1)
            pendi = ajl[0]
            coefi = aji[I[-1]]-pendi*Llimpio[I[-1]]
            aj2[I[-1]:] =  Llimpio[I[-1]:]*pendi + coefi
    else:
            aj2 = aji.copy()
    ajustec = numpy.polyfit(Llimpio,aj2,nnn)
    ajf = numpy.polyval(ajustec,L)

    NF = F/ajf

    para = numpy.empty(0)
    strli = numpy.empty(0)
    strlf = numpy.empty(0)
    unos = numpy.empty(0)
    Lmien = numpy.empty(0)
    Fmien = numpy.empty(0)
    I = numpy.empty(0)
    Llimpio = numpy.empty(0)
    Flimpio = numpy.empty(0)
    mm = numpy.empty(0)
    nn = numpy.empty(0)
    vectorL = numpy.empty(0)
    vectorF = numpy.empty(0)
    segl = numpy.empty(0)
    segf = numpy.empty(0)
    LL = numpy.empty(0)
    FF = numpy.empty(0)
    JJ = numpy.empty(0)
    RES = numpy.empty(0)
    veva = numpy.empty(0)
    AJ = numpy.empty(0)
    WI = numpy.empty(0)
    FFF = numpy.empty(0)
    LLL = numpy.empty(0)
    CON = numpy.empty(0)
    RESI = numpy.empty(0)
    loco = numpy.empty(0)
    foco = numpy.empty(0)
    aji = numpy.empty(0)
    coefi = numpy.empty(0)
    aj2 = numpy.empty(0)
    ajustec = numpy.empty(0)
    ajf = numpy.empty(0)
    L = numpy.empty(0)
    F = numpy.empty(0)
    gc.collect()
    """
    gc.collect()
    return ajustep
