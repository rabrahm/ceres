
# Makefile for Python wrapper and test programs associated with the
# JPL DE403 ephemerides.

PY  = /usr/include/python2.7/
SWIG = /usr/local/Cellar/swig/2.0.12/bin/swig
INUM = /System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/core/include/numpy/
ISOFA = ./SOFA
LSOFA = ./SOFA/lib
CPP = g++

#INCLUDE = -I. -I$(ISLALIB)
INCLUDE = -I. -I$(ISOFA)
#LIBRARIES = -L$(LSLALIB)
LIBRARIES = -L$(LSOFA)

OBJ_FILES = jplephem.o\
	    python_client_fn.o \
	    astrtime.o \
	    delay.o \
	    jpl_eph.o \
	    ptw_aux.o

EXECUTABLES = _jplephem.so \
	      alltimes \
	      asc2bin \
	      closureT \
	      delayT \
	      earthT \
	      ephemT \
	      gbxyz \
	      gc_fn_test \
	      ierspast \
	      jupiterT \
	      lst \
	      makeiers \
	      marsT \
	      mdy2mjd \
	      mercuryT \
	      mjd2mdy \
	      moonT \
	      moonvelT \
	      neptuneT \
	      nutT \
	      plutoT \
	      saturnT \
	      ssobject \
	      sunT \
	      topochk \
	      uranusT \
	      venusT

all: $(EXECUTABLES)

_jplephem.so : jplephem_wrap.o $(OBJ_FILES)
	$(CPP) -lpython2.7 -shared jplephem_wrap.o $(OBJ_FILES) $(LIBRARIES) \
					-lsofa_c -o _jplephem.so

jplephem_wrap.o : jplephem_wrap.c jplephem.h
	$(CPP) jplephem_wrap.c -c -g -I$(PY)/Include -I$(PY) $(INCLUDE) \
	-I$(INUM)

jplephem_wrap.c : jplephem.i
	$(SWIG) -python -I$(LIB) $(INCLUDE) jplephem.i

jplephem.o : jplephem.cc jplephem.h
	$(CPP) jplephem.cc -c -g $(INCLUDE) -I$(PY) -I$(INUM) \
	-I$(PY)/Include -o jplephem.o

alltimes : alltimes.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o alltimes alltimes.o astrtime.o ptw_aux.o $(LIBRARIES) \
	-lsofa_c -lm

alltimes.o : alltimes.cc astrtime.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) alltimes.cc

asc2bin : asc2bin.o
	$(CPP) -o asc2bin asc2bin.o

asc2bin.o : asc2bin.cc jpl_eph.h
	$(CPP) -c asc2bin.cc

astrtime.o : astrtime.cc astrtime.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) astrtime.cc

closureT : closureT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o closureT closureT.o jpl_eph.o astrtime.o ptw_aux.o \
	$(LIBRARIES) -lsofa_c -lm

closureT.o : closureT.cc jpl_eph.h
	$(CPP) -c closureT.cc

delay.o : delay.cc jpl_eph.h astrtime.h delay.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) delay.cc

delayT : delayT.o jpl_eph.o astrtime.o delay.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o delayT delayT.o jpl_eph.o astrtime.o ptw_aux.o delay.o \
	$(LIBRARIES) -lsofa_c -lm

delayT.o : delayT.cc jpl_eph.h astrtime.h delay.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) delayT.cc

earthT : earthT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o earthT earthT.o jpl_eph.o astrtime.o ptw_aux.o $(LIBRARIES) \
	-lsofa_c -lm

earthT.o : earthT.cc jpl_eph.h
	$(CPP) -c earthT.cc

ephemT : ephemT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o ephemT ephemT.o jpl_eph.o astrtime.o ptw_aux.o $(LIBRARIES) \
	-lsofa_c -lm

ephemT.o : ephemT.cc jpl_eph.h
	$(CPP) -c ephemT.cc

gbxyz : gbxyz.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o gbxyz gbxyz.o jpl_eph.o astrtime.o ptw_aux.o $(LIBRARIES) \
	-lsofa_c -lm

gbxyz.o : gbxyz.cc jpl_eph.h astrtime.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) gbxyz.cc

gc_fn_test : gc_fn_test.o jpl_eph.o astrtime.o delay.o python_client_fn.o \
	ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o gc_fn_test gc_fn_test.o jpl_eph.o astrtime.o delay.o \
	ptw_aux.o python_client_fn.o $(LIBRARIES) -lsofa_c -lm

gc_fn_test.o : gc_fn_test.cc python_client.h
	$(CPP) -c $(INCLUDE) gc_fn_test.cc

ierspast : ierspast.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o ierspast ierspast.o astrtime.o ptw_aux.o $(LIBRARIES) \
	-lsofa_c -lm

ierspast.o : ierspast.cc astrtime.h
	$(CPP) -c $(INCLUDE) ierspast.cc

jpl_eph.o : jpl_eph.cc jpl_eph.h astrtime.h
	$(CPP) -c jpl_eph.cc

jupiterT : jupiterT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o jupiterT jupiterT.o jpl_eph.o astrtime.o ptw_aux.o \
	$(LIBRARIES) -lsofa_c -lm

jupiterT.o : jupiterT.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) jupiterT.cc

lst : lst.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o lst lst.o astrtime.o ptw_aux.o $(LIBRARIES) -lsofa_c -lm

lst.o : lst.cc astrtime.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) lst.cc

makeiers : makeiers.o
	$(CPP) -o makeiers makeiers.o -lm

makeiers.o : makeiers.cc
	$(CPP) -c makeiers.cc

marsT : marsT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o marsT marsT.o jpl_eph.o astrtime.o ptw_aux.o $(LIBRARIES) \
	-lsofa_c -lm

marsT.o : marsT.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) marsT.cc

mdy2mjd : mdy2mjd.cc
	$(CPP) -o mdy2mjd mdy2mjd.cc

mercuryT : mercuryT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o mercuryT mercuryT.o jpl_eph.o astrtime.o ptw_aux.o \
	$(LIBRARIES) -lsofa_c -lm

mercuryT.o : mercuryT.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) mercuryT.cc

mjd2mdy : mjd2mdy.cc
	$(CPP) -o mjd2mdy mjd2mdy.cc

moonT : moonT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o moonT moonT.o jpl_eph.o astrtime.o ptw_aux.o $(LIBRARIES) \
	-lsofa_c -lm

moonT.o : moonT.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) moonT.cc

moonvelT : moonvelT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o moonvelT moonvelT.o jpl_eph.o astrtime.o ptw_aux.o \
	$(LIBRARIES) -lsofa_c -lm

moonvelT.o : moonvelT.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) moonvelT.cc

neptuneT : neptuneT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o neptuneT neptuneT.o jpl_eph.o astrtime.o ptw_aux.o \
	$(LIBRARIES) -lsofa_c -lm

neptuneT.o : neptuneT.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) neptuneT.cc

nutT : nutT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o nutT nutT.o jpl_eph.o astrtime.o ptw_aux.o $(LIBRARIES) \
	-lsofa_c -lm

nutT.o : nutT.cc jpl_eph.h
	$(CPP) -c nutT.cc

plutoT : plutoT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o plutoT plutoT.o jpl_eph.o astrtime.o ptw_aux.o $(LIBRARIES) \
	-lsofa_c -lm

plutoT.o : plutoT.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) plutoT.cc

ptw_aux.o : ptw_aux.cc $(ISOFA)/sofa.h $(ISOFA)/sofam.h
	$(CPP) -c $(INCLUDE) ptw_aux.cc

python_client_fn.o : python_client_fn.cc python_client.h \
		jpl_eph.h astrtime.h delay.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) python_client_fn.cc

saturnT : saturnT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o saturnT saturnT.o jpl_eph.o astrtime.o ptw_aux.o \
	$(LIBRARIES) -lsofa_c -lm

saturnT.o : saturnT.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) saturnT.cc

ssobject : ssobject.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o ssobject -Bstatic ssobject.o jpl_eph.o astrtime.o ptw_aux.o -Bdynamic \
	$(LIBRARIES) -lsofa_c -lm

ssobject.o : ssobject.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) ssobject.cc

sunT : sunT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o sunT sunT.o jpl_eph.o astrtime.o ptw_aux.o $(LIBRARIES) \
	-lsofa_c -lm

sunT.o : sunT.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) sunT.cc

tempoT : tempoT.o jpl_eph.o astrtime.o ptw_aux.o delay.o $(LSOFA)/libsofa_c.a
	$(CPP) -o tempoT tempoT.o jpl_eph.o astrtime.o ptw_aux.o delay.o \
	$(LIBRARIES) -lsofa_c -lm

tempoT.o : tempoT.cc jpl_eph.h astrtime.h delay.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) tempoT.cc

topochk : topochk.o astrtime.o ptw_aux.o jpl_eph.o $(LSOFA)/libsofa_c.a
	$(CPP) -o topochk topochk.o astrtime.o ptw_aux.o jpl_eph.o \
	$(LIBRARIES) -lsofa_c -lm

topochk.o : topochk.cc astrtime.h jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) topochk.cc

uranusT : uranusT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o uranusT uranusT.o jpl_eph.o astrtime.o ptw_aux.o \
	$(LIBRARIES) -lsofa_c -lm

uranusT.o : uranusT.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) uranusT.cc

venusT : venusT.o jpl_eph.o astrtime.o ptw_aux.o $(LSOFA)/libsofa_c.a
	$(CPP) -o venusT venusT.o jpl_eph.o astrtime.o ptw_aux.o $(LIBRARIES) \
	-lsofa_c -lm

venusT.o : venusT.cc jpl_eph.h $(ISOFA)/sofa.h
	$(CPP) -c $(INCLUDE) venusT.cc

clean:
	rm -f *.o *.so core*

force:
	rm -f *.o *.so core* jplephem_wrap.c jplephem_wrap.doc jplephem.py

# To test
# >>> from jplephem import *
# >>> from numarray import *

# >>> print pyslaAltaz(0.1, 1.0, 3.14159 * 38.0 / 180.0)
# >>> gb_long = -pyslaDaf2r(79, 50, 23.42)['radians']
# >>> gb_lat = pyslaDaf2r(38, 25, 59.26)['radians']
# >>> print gb_long, gb_lat

