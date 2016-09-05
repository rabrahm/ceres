C        PROGRAM TEST
C
C        REAL*8 BCVEL,HCVEL,GCVEL
C        REAL*8 DJD,DLONG,DLAT,DALT,DRA,DEC,DEQ
C	REAL*8 U,V,TCORB
C
C        DJD=2455085.807413D0
C	DJD=2455208.7033919417D0
C        DLONG=70.733D0
C        DLAT=-29.259D0
C        DALT=2378D0
C        DRA=343.97500610d0/15.0d0
C	DRA=100.0071d0/15.0d0
C        DEC=-26.65888977
C	DEC=-48.54195527d0
C        DEQ=2000.0d0
C
C        CALL BCV(DJD,DLONG,DLAT,DALT,DRA,DEC,DEQ,BCVEL,HCVEL,GCVEL,U,V,TCORB)
C
C        print *, BCVEL, HCVEL, GCVEL, TCORB
C
C        END

c*** File rvsao/Util/bcv.f
c*** January 30, 1997
c*** By G. Torres (1989)
c*** Modified by D. Mink

c  BCV calculates the correction required to reduce observed (topocentric)
c  radial velocities of a given star to the barycenter of the solar system.
c  It includes correction for the effect of the earth's rotation.
c  The maximum error of this routine is not expected to be larger than 0.6 m/sec

	Subroutine BCV (DJD,DLONG,DLAT,DALT,DRA,DEC,
     $     DEQ,BCVEL,HCVEL,GCVEL, U, V, TCORB)

	Real*8 DJD
c				Heliocentric Julian date (days)
	Real*8 DLONG
c				Geodetic longitude (degrees, west positive)
	Real*8 DLAT
c				Geodetic latitude (degrees)
	Real*8 DALT
c				Altitude above sea level (meters)
	Real*8 DRA
c				Right ascension of star (hours)
	Real*8 DEC
c				Declination of star (degrees)
	Real*8 DEQ
c				Mean equator and equinox for coordinates
c				e.g., 1950.0
	Real*8 BCVEL
c				Barycentric correction (km/s) (returned)
	Real*8 HCVEL
c				Heliocentric correction (km/s) (returned)
	Real*8 GCVEL
c				Geocentric correction (km/s) (returned)

	Real*8 DLONGS
c				Geodetic longitude (radians, west+)
	Real*8 DLATS
c				Geodetic latitude (radians)
	Real*8 DRAS
c				Right ascension of star (radians)
	Real*8 DECS
c				Declination of star (radians)

	Real*8 DC(3),DCC(3),DPREMA(3,3),DVELH(3),DVELB(3)
        Real*8 PH(3),PB(3),SLIGHT
	Real*8 DPI,DAUKM,DCTROP,DCBES,DC1900,DCT0,DTR,DST
	Real*8 DEQT,DRA2,DEC2,DHA,DARG
	Integer*4 K

        Real*8 U,V,TCORB

	Data DPI/3.1415926535897932d0/, DAUKM/1.4959787d08/,
     1	     DCTROP/365.24219572d0/, DCBES/0.313d0/, DC1900/1900.0d0/
     2	     DCT0/2415020.0d0/,SLIGHT/299792.458d0/

	DTR = DPI / 180.0d0

	DLONGS = DLONG * DTR
	DLATS = DLAT * DTR
	DRAS = DRA*15.0d0*DTR
	DECS = DEC*DTR
C	Open (15,FILE='testbcv')
C	Write(15,*) 'Julian Date:',DJD
C	Write(15,*) 'Long:',DLONGS,' Lat:',DLATS,DALT
C	Write(15,*) 'RA:',DRAS,' Dec:',DECS,DEQ
c  Calculate local sidereal time

	Call SIDTIM (DJD,DLONGS,DST)

c  Precess R.A. and Dec. to mean equator and equinox of date (DEQT)

	DEQT = (DJD - DCT0 - DCBES)/DCTROP + DC1900
	DC(1) = Dcos(DRAS) * Dcos(DECS)
	DC(2) = Dsin(DRAS) * Dcos(DECS)
	DC(3) =	      Dsin(DECS)

	Call PRE (DEQ,DEQT,DPREMA)
	Do 100 K=1,3
	    DCC(K)=DC(1)*DPREMA(K,1)+DC(2)*DPREMA(K,2)+DC(3)*DPREMA(K,3)
100	    Continue

	If (DCC(1) .ne. 0.0d0) Then
	    DARG = DCC(2) / DCC(1)
	    DRA2 = Datan (DARG)
	    If (DCC(1) .lt. 0.0d0) Then
		DRA2 = DRA2 + DPI
	    Elseif (DCC(2) .lt. 0.0d0) Then
		DRA2 = DRA2 + 2.0d0*DPI
	    Endif
	Else
	    If (DCC(2) .gt. 0.0d0) Then
		DRA2 = DPI/2.0d0
	    Else
		DRA2 = 1.5d0*DPI
	    Endif
	Endif

	DEC2 = DASIN (DCC(3))

c  Calculate hour angle = local sidereal time - R.A.

	DHA = DST - DRA2
c	DHA = Dmod (DHA + 2.0d0*DPI , 2.0d0*DPI)

c	Write(15,*) 'RA2=',DRA2,'  DEC2=',DEC2,'  ALT=',DALT
c	Write(15,*) 'ST=',DST,'  HA=',DHA,'  LAT:',DLATS

c  Calculate observer's geocentric velocity
c  (altitude assumed to be zero)

	Call GEOVEL(DLATS,DALT,DEC2,-DHA,GCVEL,U,V)

c  Calculate components of earth's barycentric velocity,
c  DVELB(I), I=1,2,3  in units of A.U./S

C
C       ORIGINAL CALL!!! -- OJO -- AJ
C
C	Call BARVEL (DJD,DEQT,DVELH,DVELB)

C	Write(15,*) 'BCV: DJD=',DJD,'  DEQT=',DEQT

        CALL BARVEL2(DJD,DEQT,DVELH,DVELB,PH,PB)
c	Write(15,*) 'BCV: DVELH=',DVELH
c	Write(15,*) 'BCV: DVELB=',DVELB
c	Write(15,*) 'BCV: DCC=',DCC
c	Write(15,*) 'BCV: DAUKM=',DAUKM



c  Project barycentric velocity to the direction of the star, and
c  convert to km/s

	BCVEL = 0.0d0
	HCVEL = 0.0d0
        TCORB  = 0.0d0
	Do 200 K=1,3
	    BCVEL = BCVEL + DVELB(K)*DCC(K)*DAUKM
	    HCVEL = HCVEL + DVELH(K)*DCC(K)*DAUKM
            TCORB = TCORB + PB(K)*DCC(K)*DAUKM
200	    Continue

        TCORB  = TCORB / SLIGHT

c	Write(15,*) BCVEL
c	Close (15)
	Return

	End
c------------------------------------------------------------------------------
c
c       Subroutine:     sidereal (DJD,DLONG,DST)
c
c       PURPOSE:	COMPUTES THE MEAN LOCAL sidereal TIME.
c
c       INPUT:	  DJD   = JULIAN DATE
c		       DLONG = OBSERVER'S LONGITUDE (RADIANS)
c
c       OUTPUT:	 DST = MEAN LOCAL sidereal TIME (RADIANS)
c
c       NOTE:	   CONSTANTS TAKEN FROM THE AMERICAN EPHEMERIS
c		       AND NAUTICAL ALMANAC, 1980)
c
c       AUTHOR:	 G. TORRES (1989)
c
c------------------------------------------------------------------------------

	Subroutine SIDTIM (DJD,DLONG,DST)

	Real*8 DJD
c			Julian Date
	Real*8 DLONG
c			Longitude
	Real*8 DST
c			Sidereal TIme (returned)

	Real*8 DTPI, DJD0, DST0, DUT, DT

	Real*8 D1,D2,D3

	Real*8 DPI, DF, DCT0, DCJUL

	DPI = 3.141592653589793d0
	DTPI = 2.0d0*DPI
	DF = 1.00273790934d0
	DCT0 = 2415020.0d0
	DCJUL = 36525.0d0

c   Constants D1,D2,D3 for calculating Greenwich Mean Sidereal Time at 0 UT
	D1 = 1.739935934667999d0
	D2 = 6.283319509909095d02
	D3 = 6.755878646261384d-06

	DJD0 = IDINT(DJD) + 0.5d0
	If (DJD0.GT.DJD)  DJD0 = DJD0 - 1.0d0
	DUT = (DJD - DJD0)*DTPI

	DT = (DJD0 - DCT0)/DCJUL
	DST0 = D1 + D2*DT + D3*DT*DT
	DST0 = Dmod (DST0,DTPI)
	DST = DF*DUT + DST0 - DLONG
	DST = Dmod (DST + 2.0d0*DTPI , DTPI)

	Return

	End
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c
c       Subroutine:     PRE (DEQ1,DEQ2,DPREMA)
c
c       PURPOSE:	CALCULATES THE MATRIX OF GENERAL PRECESSION FROM
c		       DEQ1 TO DEQ2.
c
c       INPUT:	  DEQ1 = INITIAL EPOCH OF MEAN EQUATOR AND EQUINOX
c		       DEQ2 = FINAL EPOCH OF MEAN EQUATOR AND EQUINOX
c
c       OUTPUT:	 DPREMA = 3 X 3 MATRIX OF GENERAL PRECESSION
c
c       NOTE:	   THE PRECESSION ANGLES (DZETA,DZETT,AND DTHET) ARE
c		       COMPUTED FROM THE CONSTANTS (DC1-DC9) CORRESPONDING
c		       TO THE DEFINITIONS IN THE EXPLANATORY SUPPLEMENT
c		       TO THE AMERICAN EPHEMERIS (1961, P.30F).
c
c       AUTHOR:	 P. STUMPFF (IBM-VERSION 1979): ASTRON. ASTROPHYS.
c			SUPPL. SER. 41, 1 (1980)
c		       M. H. SLOVAK (VAX 11/780 IMPLEMENTATION 1986)
c		       G. TORRES (1989)
c
c------------------------------------------------------------------------------

	Subroutine PRE (DEQ1,DEQ2,DPREMA)

	Real*8 DEQ1, DEQ2
	Real*8 DPREMA(3,3)

	Real*8 DT0, DT, DTS, DTC, DZETA, DZETT, DTHET, DSZETA, DCZETA
	Real*8 DSZETT, DCZETT, DSTHET, DCTHET, DA, DB, DC, DD

	Real*8 DCSAR, DC1900, DC1M2, DC1, DC2, DC3, DC4, DC5, DC6
	Real*8 DC7, DC8, DC9

	DCSAR = 4.848136812d-6
	DC1900 = 1900.0d0
	DC1M2 = 0.01d0
	DC1 = 2304.25d0
	DC2 = 1.396d0
	DC3 = 0.302d0
	DC4 = 0.018d0
	DC5 = 0.791d0
	DC6 = 2004.683d0
	DC7 = -0.853d0
	DC8 = -0.426d0
	DC9 = -0.042d0

	DT0 = (DEQ1 - DC1900)*DC1M2
	DT = (DEQ2 - DEQ1)*DC1M2
	DTS = DT * DT
	DTC = DTS * DT
	DZETA = ((DC1+DC2*DT0)*DT+DC3*DTS+DC4*DTC)*DCSAR
	DZETT = DZETA + DC5*DTS*DCSAR
	DTHET = ((DC6+DC7*DT0)*DT+DC8*DTS+DC9*DTC)*DCSAR
	DSZETA = Dsin(DZETA)
	DCZETA = Dcos(DZETA)
	DSZETT = Dsin(DZETT)
	DCZETT = Dcos(DZETT)
	DSTHET = Dsin(DTHET)
	DCTHET = Dcos(DTHET)
	DA = DSZETA * DSZETT
	DB = DCZETA * DSZETT
	DC = DSZETA * DCZETT
	DD = DCZETA * DCZETT

c	Write(15,*) 'DZETA=',DZETA,'  DZETT=',DZETT,'  DTHET=',DTHET
c	Write(15,*) 'CZETA=',DCZETA,'  SZETA=',DSZETA
c	Write(15,*) 'CZETT=',DCZETT,'  SZETT=',DSZETT
c	Write(15,*) 'CTHET=',DCTHET,'  STHET=',DSTHET

	DPREMA(1,1) = DD * DCTHET - DA
	DPREMA(1,2) = -1.d0 * DC * DCTHET - DB
	DPREMA(1,3) = -1.d0 * DSTHET * DCZETT
	DPREMA(2,1) = DB * DCTHET + DC
	DPREMA(2,2) = -1.d0 * DA * DCTHET + DD
	DPREMA(2,3) = -1.d0 * DSTHET * DSZETT
	DPREMA(3,1) = DCZETA * DSTHET
	DPREMA(3,2) = -1.d0 * DSZETA * DSTHET
	DPREMA(3,3) = DCTHET

	Return

	End

c------------------------------------------------------------------------------ 
c
c
c       Subroutine:     GEOVEL (DPHI,DH,DEC,DHA,DVELG,U,V)
c
c       PURPOSE:	CALCULATES THE CORRECTION REQUIRED TO TRANSFORM
c		       THE TOPOCENTRIC RADIAL VELOCITY OF A GIVEN STAR
c		       TO GEOCENTRIC.
c		       - THE MAXIMUM ERROR OF THIS ROUTINE IS NOT EXPECTED
c		       TO BE LARGER THAN 0.1 M/S.
c
c       INPUT:	  DPHI = OBSERVER'S GEODETIC LATITUDE (RADIANS)
c		       DH = OBSERVER'S ALTITUDE ABOVE SEA LEVEL (METERS)
c
c		       DEC = STAR'S DECLINATION (RADIANS) FOR MEAN
c			      EQUATOR AND EQUINOX OF DATE
c		       DHA  = HOUR ANGLE (RADIANS)
c
c       OUTPUT:	 DVELG = GEOCENTRIC CORRECTION (KM/S)
c
c       NOTES:	  VR = R.W.COS(DEC).SIN(HOUR ANGLE), WHERE R =
c		       GEOCENTRIC RADIUS AT OBSERVER'S LATITUDE AND
c		       ALTITUDE, AND W = 2.PI/T, T = LENGTH OF sidereal
c		       DAY (SEC).  THE HOUR ANGLE IS POSITIVE EAST OF
c			THE MERIDIAN.
c			OTHER RELEVANT EQUATIONS FROM E. W. WOOLARD
c		       & G. M. CLEMENCE (1966), SPHERICAL ASTRONOMY,
c			P.45 AND P.48
c
c       AUTHOR:	 G. TORRES (1989)
c
C       MODIFIED: 07.Dec.2009 -- Andres Jordan
C                 Outputs also U,V now, where
C
C                 U: distance from Earth spin axis (km)
C                 V: distance north of equatorial plane (km)
C
c------------------------------------------------------------------------------

	Subroutine GEOVEL(DPHI,DH,DEC,DHA,DVELG,U,V)

	Real*8 DPHI,DH,DEC,DHA,DVELG
        Real*8 U,V

	Real*8 DE2,D1,D2,DR0,DA,DF,DW,DPHIG,DRH

c  Earth's equatorial radius (KM)
	DA = 6378.140d0

c  Polar flattening
	DF = 0.00335281d0

c  Angular rotation rate (2.PI/T)
	DW = 7.2921158554d-05

	DE2 = DF*(2.0d0 - DF)

c  Calculate geocentric radius DR0 at sea level (KM)
	D1 = 1.0d0 - DE2*(2.0d0 - DE2)*Dsin(DPHI)**2
	D2 = 1.0d0 - DE2*Dsin(DPHI)**2
	DR0 = DA * DSQRT(D1/D2)

c  Calculate geocentric latitude DPHIG
	D1 = DE2*Dsin(2.0d0*DPHI)
	D2 = 2.0d0*D2
	DPHIG = DPHI - DataN(D1/D2)

c  Calculate geocentric radius DRH at altitude DH (KM)
	DRH = DR0*Dcos(DPHIG) + DH/1.0D3*Dcos(DPHI)

c  Projected component to star at declination = DEC and
c  at hour angle = DHA, in units of km/s
	DVELG = DW * DRH * Dcos(DEC) * Dsin(DHA)

c  U,V
        U = DRH * Dcos(DPHIG)
        V = DRH * Dsin(DPHIG)
	Return

	End

	Subroutine BARVEL2 (DJE,DEQ,DVELH,DVELB, PH,PB)
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c  
c     This subroutine is identical to BARVEL, but uses IAU routines from SOFA
c     to compute the heliocentric and barycentric velocities. These have maximum
c     deviations from JPL ephemeris of 46 mm/sc, much better than the 42 cm/s of the
c     typical Stumpff subroutine
c
c     Written Dec.07.2009 -- Andres Jordan
c
	Real*8 DJE,DEQ,DEQDAT
	Real*8 DVELH(3)
c				Heliocentric velocity correction (returned)
	Real*8 DVELB(3)
c				Barycentric velocity correction (returned)
     
        Real*8 PH(3), PB(3)
 	Real*8 DCT0,DCBES,DCTROP,DC1900

c       Velocities we will obtain from IAU routiner
        Double Precision PVH(3,2), PVB(3,2)
        INTEGER JSTAT

	Real*8 DPREMA(3,3)

        REAL*8 SECDAY

        SECDAY = 24D0 * 3600D0

	Data  DCT0/2415020.0d0/,DCBES/0.313d0/,DCTROP/365.24219572d0/,
     1       DC1900/1900.0d0/

        Call iau_EPV00 ( DJE, 0.0d0, PVH , PVB, JSTAT )


	DEQDAT = (DJE - DCT0 - DCBES)/DCTROP + DC1900

	Call PRE(DEQDAT,DEQ,DPREMA)

	Do 801 N = 1 , 3
	    DVELH(N) = PVH(1,2)*DPREMA(N,1) + PVH(2,2)*DPREMA(N,2) + 
     1		       PVH(3,2)*DPREMA(N,3)
	    DVELB(N) = PVB(1,2)*DPREMA(N,1) + PVB(2,2)*DPREMA(N,2) + 
     1		       PVB(3,2)*DPREMA(N,3)
            PH(N)    = PVH(1,1)*DPREMA(N,1) + PVH(2,1)*DPREMA(N,2) + 
     1		       PVH(3,1)*DPREMA(N,3)
            PB(N)    = PVB(1,1)*DPREMA(N,1) + PVB(2,1)*DPREMA(N,2) + 
     1		       PVB(3,1)*DPREMA(N,3)

 801     Continue

	Do 901 N = 1 , 3
	    DVELH(N) = DVELH(N) / SECDAY
	    DVELB(N) = DVELB(N) / SECDAY

 901     Continue


	Return

        End

c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c
c       Subroutine:     BARVEL (DJE,DEQ,DVELH,DVELB)
c
c       PURPOSE:	CALCULATES THE HELIOCENTRIC AND BARYCENTRIC
c		       VELOCITY COMPONENTS OF THE EARTH.  THE LARGEST
c		       DEVIATIONS FROM THE JPL-DE96 EPHEMERIS ARE 42 CM/S
c		       FOR BOTH HELIOCENTRIC AND BARYCENTRIC VELOCITY
c		       COMPONENTS.
c
c       INPUT:	  DJE = JULIAN EPHERMERIS DATE
c		       DEQ = EPOCH OF MEAN EQUATOR AND MEAN EQUINOX OF DVELH
c			     AND DVELB.  If DEQ = 0, BOTH VECTORS ARE REFER-
c			     RED TO THE MEAN EQUATOR AND EQUINOX OF DJE.
c
c       OUTPUT:	 DVELH(K) = HELIOCENTRIC VELOCITY COMPONENTS
c		       DVELB(K) = BARYCENTRIC VELOCITY COMPONENTS
c			(DX/DT, DY/DT, DZ/DT, K=1,2,3  A.U./S)
c
c       AUTHOR:	 P. STUMPFF  (IBM-VERSION 1979): ASTRON. ASTROPHYS.
c			SUPPL. SER. 41, 1 (1980)
c		       M. H. SLOVAK (VAX 11/780 IMPLEMENTATION 1986)
c		       G. TORRES (1989)
c
c------------------------------------------------------------------------------

	Subroutine BARVEL (DJE,DEQ,DVELH,DVELB)

	Real*8 DJE,DEQ
	Real*8 DVELH(3)
c				Heliocentric velocity correction (returned)
	Real*8 DVELB(3)
c				Barycentric velocity correction (returned)

	Integer*4 K, N
	Real*8 DCFEL(3,8),DCEPS(3),DCARGS(2,15)
	Real*8 DCARGM(2,3), E, G
	Real*8 CC2PI,CCSEC3,CCKM,CCMLD,CCFDI,T,TSQ,A
	Real*8 DC1MME,CCSGD,DT,DTSQ,DML,DLOCAL,PERTLD,PERTR,PERTRD
	Real*8 COSA,SINA,ESQ,PARAM,DPARAM,TWOE,TWOG,F,SINF,COSF
	Real*8 DYHD,DZHD,B,TL,PERTPD,PERTP,DRLD,DRD,DXHD,PHID,PSID
	Real*8 PLON,POMG,PECC,DYAHD,DZAHD,DYABD,DZABD,DEQDAT
	Real*8 DCSLD,DXBD,DYBD,DZBD
	Real*8 CCAMPS(5,15),CCSEC(3,4)
	Real*8 CCPAMV(4),SN(4),CCSEL(3,17),CCAMPM(4,3)

	Common /BARXYZ/ DPREMA(3,3),DPSI,D1PDRO,DSINLS,DCOSLS,DSINEP,
     1			DCOSEP,FORBEL(7),SORBEL(17),SINLP(4),COSLP(4),
     2			SINLM,COSLM,SIGMA,IDEQ
	Real*8 DPREMA,DPSI,D1PDRO,DSINLS,DCOSLS,DSINEP,DCOSEP
	Real*8 FORBEL,SORBEL,SINLP,COSLP,SINLM,COSLM,SIGMA
	Integer*4 IDEQ

	Equivalence (SORBEL(1),E),(FORBEL(1),G)

	Real*8 DC2PI,DC1,DCT0,DCJUL,DCBES,DCTROP,DC1900,DTL,DEPS
	Real*8 PHI,PERTL

	Data DC2PI/6.2831853071796d0/, CC2PI/6.283185/,
     1	     DC1  /1.0d0/ , DCT0/2415020.0d0/, DCJUL/36525.0d0/,
     2	     DCBES/0.313d0/, DCTROP/365.24219572d0/, DC1900/1900.0d0/

c  Constants DCFEL(I,K) of fast-changing elements

c		       I = 1	     I = 2	   I = 3

	Data DCFEL/ 1.7400353d+00, 6.2833195099091d+02, 5.2796d-06,
     1	      6.2565836d+00, 6.2830194572674d+02,-2.6180d-06,
     1	      4.7199666d+00, 8.3997091449254d+03,-1.9780d-05,
     1	      1.9636505d-01, 8.4334662911720d+03,-5.6044d-05,
     1	      4.1547339d+00, 5.2993466764997d+01, 5.8845d-06,
     1	      4.6524223d+00, 2.1354275911213d+01, 5.6797d-06,
     1	      4.2620486d+00, 7.5025342197656d+00, 5.5317d-06,
     1	      1.4740694d+00, 3.8377331909193d+00, 5.6093d-06/

c   CONSTANTS DCEPS AND CCSEL(I,K) OF SLOWLY CHANGING ELEMENTS

c		       I = 1	I = 2	 I = 3

	Data DCEPS/ 4.093198d-01,-2.271110d-04,-2.860401d-08/

	Data CCSEL/ 1.675104d-02,-4.179579d-05,-1.260516d-07,
     1	      2.220221d-01, 2.809917d-02, 1.852532d-05,
     1	      1.589963d+00, 3.418075d-02, 1.430200d-05,
     1	      2.994089d+00, 2.590824d-02, 4.155840d-06,
     1	      8.155457d-01, 2.486352d-02, 6.836840d-06,
     1	      1.735614d+00, 1.763719d-02, 6.370440d-06,
     1	      1.968564d+00, 1.524020d-02,-2.517152d-06,
     1	      1.282417d+00, 8.703393d-03, 2.289292d-05,
     1	      2.280820d+00, 1.918010d-02, 4.484520d-06,
     1	      4.833473d-02, 1.641773d-04,-4.654200d-07,
     1	      5.589232d-02,-3.455092d-04,-7.388560d-07,
     1	      4.634443d-02,-2.658234d-05, 7.757000d-08,
     1	      8.997041d-03, 6.329728d-06,-1.939256d-09,
     1	      2.284178d-02,-9.941590d-05, 6.787400d-08,
     1	      4.350267d-02,-6.839749d-05,-2.714956d-07,
     1	      1.348204d-02, 1.091504d-05, 6.903760d-07,
     1	      3.106570d-02,-1.665665d-04,-1.590188d-07/

c   CONSTANTS OF THE ARGUMENTS OF THE SHORT-PERIOD PERTURBATIONS BY
c   THE PLANETS:  DCARGS(I,K)

c			I = 1	     I = 2

	Data DCARGS/ 5.0974222d+00,-7.8604195454652d+02,
     1	       3.9584962d+00,-5.7533848094674d+02,
     1	       1.6338070d+00,-1.1506769618935d+03,
     1	       2.5487111d+00,-3.9302097727326d+02,
     1	       4.9255514d+00,-5.8849265665348d+02,
     1	       1.3363463d+00,-5.5076098609303d+02,
     1	       1.6072053d+00,-5.2237501616674d+02,
     1	       1.3629480d+00,-1.1790629318198d+03,
     1	       5.5657014d+00,-1.0977134971135d+03,
     1	       5.0708205d+00,-1.5774000881978d+02,
     1	       3.9318944d+00, 5.2963464780000d+01,
     1	       4.8989497d+00, 3.9809289073258d+01,
     1	       1.3097446d+00, 7.7540959633708d+01,
     1	       3.5147141d+00, 7.9618578146517d+01,
     1	       3.5413158d+00,-5.4868336758022d+02/

c   AMPLITUDES CCAMPS(N,K) OF THE SHORT-PERIOD PERTURBATIONS

c	  N = 1	N = 2	N = 3	N = 4	N = 5

	Data CCAMPS/
     1 -2.279594d-5, 1.407414d-5, 8.273188d-6, 1.340565d-5,-2.490817d-7,
     1 -3.494537d-5, 2.860401d-7, 1.289448d-7, 1.627237d-5,-1.823138d-7,
     1  6.593466d-7, 1.322572d-5, 9.258695d-6,-4.674248d-7,-3.646275d-7,
     1  1.140767d-5,-2.049792d-5,-4.747930d-6,-2.638763d-6,-1.245408d-7,
     1  9.516893d-6,-2.748894d-6,-1.319381d-6,-4.549908d-6,-1.864821d-7,
     1  7.310990d-6,-1.924710d-6,-8.772849d-7,-3.334143d-6,-1.745256d-7,
     1 -2.603449d-6, 7.359472d-6, 3.168357d-6, 1.119056d-6,-1.655307d-7,
     1 -3.228859d-6, 1.308997d-7, 1.013137d-7, 2.403899d-6,-3.736225d-7,
     1  3.442177d-7, 2.671323d-6, 1.832858d-6,-2.394688d-7,-3.478444d-7,
     1  8.702406d-6,-8.421214d-6,-1.372341d-6,-1.455234d-6,-4.998479d-8,
     1 -1.488378d-6,-1.251789d-5, 5.226868d-7,-2.049301d-7, 0.0d0,
     1 -8.043059d-6,-2.991300d-6, 1.473654d-7,-3.154542d-7, 0.0d0,
     1  3.699128d-6,-3.316126d-6, 2.901257d-7, 3.407826d-7, 0.0d0,
     1  2.550120d-6,-1.241123d-6, 9.901116d-8, 2.210482d-7, 0.0d0,
     1 -6.351059d-7, 2.341650d-6, 1.061492d-6, 2.878231d-7, 0.0d0/

c   CONSTANTS OF THE SECULAR PERTURBATIONS IN LONGITUDE CCSEC3 AND
c   CCSEC(N,K)

c			N = 1	N = 2	  N = 3

	Data CCSEC3/-7.757020d-08/

	Data CCSEC/  1.289600d-06, 5.550147d-01, 2.076942d+00,
     1	       3.102810d-05, 4.035027d+00, 3.525565d-01,
     1	       9.124190d-06, 9.990265d-01, 2.622706d+00,
     1	       9.793240d-07, 5.508259d+00, 1.559103d+01/

c   Sidereal RATE DCSLD IN LONGITUDE, RATE CCSGD IN MEAN ANOMALY

	Data DCSLD/ 1.990987d-07/, CCSGD/ 1.990969d-07/

c   SOME CONSTANTS USED IN THE CALCULATION OF THE LUNAR CONTRIBUTION

	Data CCKM/3.122140d-05/, CCMLD/2.661699d-06/, CCFDI/2.399485d-07/

c   CONSTANTS DCARGM(I,K) OF THE ARGUMENTS OF THE PERTURBATIONS OF THE
c   MOTION OF THE MOON

c			 I = 1	     I = 2

	Data DCARGM/  5.1679830d+00, 8.3286911095275d+03,
     1		5.4913150d+00,-7.2140632838100d+03,
     1		5.9598530d+00, 1.5542754389685d+04/

c   AMPLITUDES CCAMPM(N,K) OF THE PERTURBATIONS OF THE MOON

c	   N = 1	 N = 2	 N = 3	 N = 4

	Data CCAMPM/
     1	 1.097594d-01, 2.896773d-07, 5.450474d-02, 1.438491d-07,
     1	-2.223581d-02, 5.083103d-08, 1.002548d-02,-2.291823d-08,
     1	 1.148966d-02, 5.658888d-08, 8.249439d-03, 4.063015d-08/

c  CCPAMV = A*M*DL/DT (PLANETS); DC1MME = 1 - MASS(EARTH+MOON)

	Data CCPAMV/8.326827d-11,1.843484d-11,1.988712d-12,1.881276d-12/,
     1	     DC1MME/0.99999696d0/

c  Program execution begins

c  Control-parameter IDEQ, and time-arguments

	IDEQ = DEQ
	DT = (DJE - DCT0)/DCJUL
	T = DT
	DTSQ = DT * DT
	TSQ = DTSQ

c  Values of all elements for the instant DJE

	Do 100 K = 1 , 8
	    DLOCAL = Dmod (DCFEL(1,K)+DT*DCFEL(2,K)+DTSQ*DCFEL(3,K), DC2PI)
	    If(K.EQ.1) DML = DLOCAL
	    If(K .ne. 1) FORBEL(K-1) = DLOCAL
100	    Continue

	DEPS = Dmod (DCEPS(1)+DT*DCEPS(2)+DTSQ*DCEPS(3), DC2PI)

	Do 200 K = 1 , 17
	    SORBEL(K) = Dmod (CCSEL(1,K)+T*CCSEL(2,K)+TSQ*CCSEL(3,K), DC2PI)
200	    Continue

c  Secular perturbations in longitude

	Do 300 K = 1 , 4
	    A = Dmod (CCSEC(2,K)+T*CCSEC(3,K), DC2PI)
	    SN(K) = SIN(A)
300	    Continue

c  PERIODIC PERTURBATIONS OF THE EMB (EARTH-MOON BARYCENTER)

	PERTL = CCSEC(1,1)	    *SN(1) +CCSEC(1,2)*SN(2)
     1		+(CCSEC(1,3)+T*CCSEC3)*SN(3) +CCSEC(1,4)*SN(4)

	PERTLD = 0.0
	PERTR =  0.0
	PERTRD = 0.0

	DO 400 K = 1 , 15
	    A = Dmod (DCARGS(1,K)+DT*DCARGS(2,K), DC2PI)
	    COSA = COS(A)
	    SINA = SIN(A)
	    PERTL = PERTL + CCAMPS(1,K)*COSA + CCAMPS(2,K)*SINA
	    PERTR = PERTR + CCAMPS(3,K)*COSA + CCAMPS(4,K)*SINA
	    If(K.GE.11) Go to 400
	    PERTLD = PERTLD + (CCAMPS(2,K)*COSA-CCAMPS(1,K)*SINA)*CCAMPS(5,K)
	    PERTRD = PERTRD + (CCAMPS(4,K)*COSA-CCAMPS(3,K)*SINA)*CCAMPS(5,K)
400	    Continue

c  ELLIPTIC PART OF THE MOTION OF THE EMB

	ESQ = E * E
	DPARAM = DC1 - ESQ
	PARAM = DPARAM
	TWOE = E + E
	TWOG = G + G
	PHI = TWOE*((1.d0 - ESQ*(1.d0/8.d0))*SIN(G) + E*(5.d0/8.d0)*SIN(TWOG)
     1	      + ESQ*0.5416667d0*SIN(G+TWOG) )
	F = G + PHI
	SINF = SIN(F)
	COSF = COS(F)
	DPSI = DPARAM/(DC1 + E*COSF)
	PHID = TWOE*CCSGD*((1.d0+ESQ*1.5d0)*COSF+E*(1.25d0-SINF*SINF*0.5d0))
	PSID = CCSGD*E*SINF/SQRT(PARAM)

c  PERTURBED HELIOCENTRIC MOTION OF THE EMB

	D1PDRO = (DC1 + PERTR)
	DRD = D1PDRO*(PSID+DPSI*PERTRD)
	DRLD = D1PDRO*DPSI*(DCSLD+PHID+PERTLD)
	DTL = Dmod (DML+PHI+PERTL, DC2PI)
	DSINLS = Dsin(DTL)
	DCOSLS = Dcos(DTL)
	DXHD = DRD*DCOSLS - DRLD*DSINLS
	DYHD = DRD*DSINLS + DRLD*DCOSLS

c  INFLUENCE OF ECCENTRICITY, EVECTION AND VARIATION OF THE GEOCENTRIC
c  MOTION OF THE MOON

	PERTL =  0.d0
	PERTLD = 0.d0
	PERTP =  0.d0
	PERTPD = 0.d0

	Do 500 K = 1, 3
	    A = Dmod (DCARGM(1,K) + DT*DCARGM(2,K), DC2PI)
	    SINA = SIN(A)
	    COSA = COS(A)
	    PERTL   = PERTL  + CCAMPM(1,K)*SINA
	    PERTLD  = PERTLD + CCAMPM(2,K)*COSA
	    PERTP   = PERTP  + CCAMPM(3,K)*COSA
	    PERTPD  = PERTPD - CCAMPM(4,K)*SINA
500	    Continue

c   HELIOCENTRIC MOTION OF THE EARTH

	TL =  FORBEL(2) + PERTL
	SINLM = SIN(TL)
	COSLM = COS(TL)
	SIGMA = CCKM/(1.0 + PERTP)
	A = SIGMA*(CCMLD+PERTLD)
	B = SIGMA*PERTPD
	DXHD = DXHD + A*SINLM + B*COSLM
	DYHD = DYHD - A*COSLM + B*SINLM
	DZHD =      - SIGMA*CCFDI*COS(FORBEL(3))

c   BARYCENTRIC MOTION OF THE EARTH

	DXBD = DXHD*DC1MME
	DYBD = DYHD*DC1MME
	DZBD = DZHD*DC1MME

	Do 600 K = 1 , 4
	    PLON = FORBEL(K+3)
	    POMG = SORBEL(K+1)
	    PECC = SORBEL(K+9)
	    TL = Dmod (PLON+2.0*PECC*SIN(PLON-POMG), CC2PI)
	    SINLP(K) = SIN(TL)
	    COSLP(K) = COS(TL)
	    DXBD = DXBD + CCPAMV(K)*(SINLP(K) + PECC*SIN(POMG))
	    DYBD = DYBD - CCPAMV(K)*(COSLP(K) + PECC*COS(POMG))
	    DZBD = DZBD - CCPAMV(K)*SORBEL(K+13)*COS(PLON-SORBEL(K+5))
600	    Continue

c   TRANSITION TO MEAN EQUATOR OF DATE

	  DCOSEP = Dcos(DEPS)
	  DSINEP = Dsin(DEPS)
	  DYAHD = DCOSEP*DYHD - DSINEP*DZHD
	  DZAHD = DSINEP*DYHD + DCOSEP*DZHD
	  DYABD = DCOSEP*DYBD - DSINEP*DZBD
	  DZABD = DSINEP*DYBD + DCOSEP*DZBD
c
	  If (IDEQ .ne. 0)  Go to 700
	  DVELH(1) = DXHD
	  DVELH(2) = DYAHD
	  DVELH(3) = DZAHD
	  DVELB(1) = DXBD
	  DVELB(2) = DYABD
	  DVELB(3) = DZABD

	  Return

c   GENERAL PRECESSION FROM EPOCH DJE TO DEQ

700	DEQDAT = (DJE - DCT0 - DCBES)/DCTROP + DC1900

	Call PRE(DEQDAT,DEQ,DPREMA)

	Do 800 N = 1 , 3
	    DVELH(N) = DXHD*DPREMA(N,1) + DYAHD*DPREMA(N,2) + 
     1		       DZAHD*DPREMA(N,3)
	    DVELB(N) = DXBD*DPREMA(N,1) + DYABD*DPREMA(N,2) + 
     1		       DZABD*DPREMA(N,3)
800	    Continue

	Return

	End

c	1990	Return HCV AND BCV

c Apr  1 1992	Use altitude when computing velocity correction
c Aug 14 1992	Return 0 from JULDAY if year is <= 0

c Apr 13 1994	Declare more variables
c Jun 15 1994	Fix parameter descriptions
c Nov 16 1994	Remove ^L's between subroutines
c Nov 16 1994	Shorten three 80-character lines

c Jul 13 1995	Move Julian Date to SPP subroutine

c Jan  3 1997	Note that Julian Date is heliocentric
c Jan 24 1997	Make all constants explicitly double precision
c Jan 30 1997	Fix to match RVCORRECT

      SUBROUTINE iau_CAL2JD ( IY, IM, ID, DJM0, DJM, J )
*+
*  - - - - - - - - - - -
*   i a u _ C A L 2 J D
*  - - - - - - - - - - -
*
*  Gregorian Calendar to Julian Date.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     IY,IM,ID    i     year, month, day in Gregorian calendar (Note 1)
*
*  Returned:
*     DJM0        d     MJD zero-point: always 2400000.5
*     DJM         d     Modified Julian Date for 0 hrs
*     J           i     status:
*                           0 = OK
*                          -1 = bad year   (Note 3: JD not computed)
*                          -2 = bad month  (JD not computed)
*                          -3 = bad day    (JD computed)
*
*  Notes:
*
*  1) The algorithm used is valid from -4800 March 1, but this
*     implementation rejects dates before -4799 January 1.
*
*  2) The Julian Date is returned in two pieces, in the usual SOFA
*     manner, which is designed to preserve time resolution.  The
*     Julian Date is available as a single number by adding DJM0 and
*     DJM.
*
*  3) In early eras the conversion is from the "Proleptic Gregorian
*     Calendar";  no account is taken of the date(s) of adoption of
*     the Gregorian Calendar, nor is the AD/BC numbering convention
*     observed.
*
*  Reference:
*
*     Explanatory Supplement to the Astronomical Almanac,
*     P. Kenneth Seidelmann (ed), University Science Books (1992),
*     Section 12.92 (p604).
*
*  This revision:  2001 September 16
*
*  Copyright (C) 2008 IAU SOFA Review Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER IY, IM, ID
      DOUBLE PRECISION DJM0, DJM
      INTEGER J, MY, IYPMY

*  Earliest year allowed (4800BC)
      INTEGER IYMIN
      PARAMETER ( IYMIN = -4799 )

*  Month lengths in days
      INTEGER MTAB(12)
      DATA MTAB / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Preset status.
      J = 0

*  Validate year.
      IF ( IY.LT.IYMIN ) THEN
         J = -1
      ELSE

*     Validate month.
         IF ( IM.GE.1 .AND. IM.LE.12 ) THEN

*        Allow for leap year.
            IF ( MOD(IY,4) .EQ. 0 ) THEN
               MTAB(2) = 29
            ELSE
               MTAB(2) = 28
            END IF
            IF ( MOD(IY,100).EQ.0 .AND. MOD(IY,400).NE.0 ) MTAB(2) = 28

*        Validate day.
            IF ( ID.LT.1 .OR. ID.GT.MTAB(IM) ) J = -3

*        Result.
            MY = ( IM - 14 ) / 12
            IYPMY = IY + MY
            DJM0 = 2400000.5D0
            DJM = DBLE( ( 1461 * ( IYPMY + 4800 ) ) / 4
     :                + (  367 * ( IM-2 - 12*MY ) ) / 12
     :                - (    3 * ( ( IYPMY + 4900 ) / 100 ) ) / 4
     :                + ID - 2432076)

*        Bad month
         ELSE
            J = -2
         END IF
      END IF

*  Finished.

*+-----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  Standards Of Fundamental Astronomy Review Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Review Board ("the Board").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and 
*     restrictions listed below.
*
*  3. You (the user) may copy and adapt the SOFA software and its 
*     algorithms for your own purposes and you may copy and distribute
*     a resulting "derived work" to others on a world-wide, royalty-free 
*     basis, provided that the derived work complies with the following
*     requirements: 
*
*     a) Your work shall be marked or carry a statement that it (i) uses
*        routines and computations derived by you from software provided 
*        by SOFA under license to you; and (ii) does not contain
*        software provided by SOFA or software that has been distributed
*        by or endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon and/or differs from the
*        original SOFA software.
*
*     c) The name(s) of all routine(s) that you distribute shall differ
*        from the SOFA names, even when the SOFA content has not been
*        otherwise changed.
*
*     d) The routine-naming prefix "iau" shall not be used.
*
*     e) The origin of the SOFA components of your derived work must not
*        be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     f) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have granted 
*        a further right to modify the source code of your derived work.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall acknowledge
*     that the SOFA software was used in obtaining those results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or by
*     inappropriate modification.
*
*  6. The SOFA software is provided "as is" and the Board makes no 
*     warranty as to its use or performance.   The Board does not and 
*     cannot warrant the performance or results which the user may obtain 
*     by using the SOFA software.  The Board makes no warranties, express 
*     or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Board be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Board representative has been advised of such
*     damages, or for any claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms 
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*     Internet email: sofa@rl.ac.uk
*     Postal address: IAU SOFA Center
*                     Rutherford Appleton Laboratory
*                     Chilton, Didcot, Oxon OX11 0QX
*                     United Kingdom
*
*-----------------------------------------------------------------------

      END
      SUBROUTINE iau_DAT ( IY, IM, ID, FD, DELTAT, J )
*+
*  - - - - - - - -
*   i a u _ D A T
*  - - - - - - - -
*
*  For a given UTC date, calculate delta(AT) = TAI-UTC.
*
*     :------------------------------------------:
*     :                                          :
*     :                 IMPORTANT                :
*     :                                          :
*     :  A new version of this routine must be   :
*     :  produced whenever a new leap second is  :
*     :  announced.  There are five items to     :
*     :  change on each such occasion:           :
*     :                                          :
*     :  1) The parameter NDAT must be           :
*     :     increased by 1.                      :
*     :                                          :
*     :  2) A new line must be added to the set  :
*     :     of DATA statements that initialize   :
*     :     the arrays IDATE and DATS.           :
*     :                                          :
*     :  3) The parameter IYV must be set to     :
*     :     the current year.                    :
*     :                                          :
*     :  4) The "Latest leap second" comment     :
*     :     below must be set to the new leap    :
*     :     second date.                         :
*     :                                          :
*     :  5) The "This revision" comment, later,  :
*     :     must be set to the current date.     :
*     :                                          :
*     :  Change (3) must also be carried out     :
*     :  whenever the routine is re-issued,      :
*     :  even if no leap seconds have been       :
*     :  added.                                  :
*     :                                          :
*     :  Latest leap second:  2008 December 31   :
*     :                                          :
*     :__________________________________________:
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     IY       i     UTC:  year (Notes 1 and 2)
*     IM       i           month (Note 2)
*     ID       i           day (Notes 2 and 3)
*     FD       d           fraction of day (Note 4)
*
*  Returned:
*     DELTAT   d     TAI minus UTC, seconds
*     J        i     status (Note 5):
*                       1 = dubious year (Note 1)
*                       0 = OK
*                      -1 = bad year
*                      -2 = bad month
*                      -3 = bad day (Note 3)
*                      -4 = bad fraction (Note 4)
*
*  Notes:
*
*  1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
*     to call the routine with an earlier date.  If this is attempted,
*     zero is returned together with a warning status.
*
*     Because leap seconds cannot, in principle, be predicted in
*     advance, a reliable check for dates beyond the valid range is
*     impossible.  To guard against gross errors, a year five or more
*     after the release year of the present routine (see parameter IYV)
*     is considered dubious.  In this case a warning status is returned
*     but the result is computed in the normal way.
*
*     For both too-early and too-late years, the warning status is J=+1.
*     This is distinct from the error status J=-1, which signifies a
*     year so early that JD could not be computed.
*
*  2) If the specified date is for a day which ends with a leap second,
*     the UTC-TAI value returned is for the period leading up to the
*     leap second.  If the date is for a day which begins as a leap
*     second ends, the UTC-TAI returned is for the period following the
*     leap second.
*
*  3) The day number must be in the normal calendar range, for example
*     1 through 30 for April.  The "almanac" convention of allowing
*     such dates as January 0 and December 32 is not supported in this
*     routine, in order to avoid confusion near leap seconds.
*
*  4) The fraction of day is used only for dates before the introduction
*     of leap seconds, the first of which occurred at the end of 1971.
*     It is tested for validity (zero to less than 1 is the valid range)
*     even if not used;  if invalid, zero is used and status J=-4 is
*     returned.  For many applications, setting FD to zero is
*     acceptable;  the resulting error is always less than 3 ms (and
*     occurs only pre-1972).
*
*  5) The status value returned in the case where there are multiple
*     errors refers to the first error detected.  For example, if the
*     month and day are 13 and 32 respectively, J=-2 (bad month) will be
*     returned.
*
*  6) In cases where a valid result is not available, zero is returned.
*
*  References:
*
*  1) For dates from 1961 January 1 onwards, the expressions from the
*     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
*
*  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
*     the 1992 Explanatory Supplement.
*
*  Called:
*     iau_CAL2JD   Gregorian calendar to Julian Day number
*
*  This revision:  2008 July 5
*
*  Copyright (C) 2008 IAU SOFA Review Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER IY, IM, ID
      DOUBLE PRECISION FD, DELTAT
      INTEGER J

*  Release year for this version of iau_DAT
      INTEGER IYV
      PARAMETER ( IYV = 2009 )

*  Number of Delta(AT) changes (increase by 1 for each new leap second)
      INTEGER NDAT
      PARAMETER ( NDAT = 39 )

*  Number of Delta(AT) expressions before leap seconds were introduced
      INTEGER NERA1
      PARAMETER ( NERA1 = 14 )

*  Dates (year, month) on which new Delta(AT) came into force
      INTEGER IDATE(2,NDAT)

*  New Delta(AT) which came into force on the given dates
      DOUBLE PRECISION DATS(NDAT)

*  Reference dates (MJD) and drift rates (s/day), pre leap seconds
      DOUBLE PRECISION DRIFT(2,NERA1)

*  Miscellaneous local variables
      LOGICAL MORE
      INTEGER JS, I, M, IS
      DOUBLE PRECISION DAT, DJM0, DJM

*  Dates and Delta(AT)s
      DATA (IDATE(I, 1),I=1,2),DATS(1)  / 1960,  1,  1.4178180D0 /
      DATA (IDATE(I, 2),I=1,2),DATS(2)  / 1961,  1,  1.4228180D0 /
      DATA (IDATE(I, 3),I=1,2),DATS(3)  / 1961,  8,  1.3728180D0 /
      DATA (IDATE(I, 4),I=1,2),DATS(4)  / 1962,  1,  1.8458580D0 /
      DATA (IDATE(I, 5),I=1,2),DATS(5)  / 1963, 11,  1.9458580D0 /
      DATA (IDATE(I, 6),I=1,2),DATS(6)  / 1964,  1,  3.2401300D0 /
      DATA (IDATE(I, 7),I=1,2),DATS(7)  / 1964,  4,  3.3401300D0 /
      DATA (IDATE(I, 8),I=1,2),DATS(8)  / 1964,  9,  3.4401300D0 /
      DATA (IDATE(I, 9),I=1,2),DATS(9)  / 1965,  1,  3.5401300D0 /
      DATA (IDATE(I,10),I=1,2),DATS(10) / 1965,  3,  3.6401300D0 /
      DATA (IDATE(I,11),I=1,2),DATS(11) / 1965,  7,  3.7401300D0 /
      DATA (IDATE(I,12),I=1,2),DATS(12) / 1965,  9,  3.8401300D0 /
      DATA (IDATE(I,13),I=1,2),DATS(13) / 1966,  1,  4.3131700D0 /
      DATA (IDATE(I,14),I=1,2),DATS(14) / 1968,  2,  4.2131700D0 /
      DATA (IDATE(I,15),I=1,2),DATS(15) / 1972,  1, 10D0 /
      DATA (IDATE(I,16),I=1,2),DATS(16) / 1972,  7, 11D0 /
      DATA (IDATE(I,17),I=1,2),DATS(17) / 1973,  1, 12D0 /
      DATA (IDATE(I,18),I=1,2),DATS(18) / 1974,  1, 13D0 /
      DATA (IDATE(I,19),I=1,2),DATS(19) / 1975,  1, 14D0 /
      DATA (IDATE(I,20),I=1,2),DATS(20) / 1976,  1, 15D0 /
      DATA (IDATE(I,21),I=1,2),DATS(21) / 1977,  1, 16D0 /
      DATA (IDATE(I,22),I=1,2),DATS(22) / 1978,  1, 17D0 /
      DATA (IDATE(I,23),I=1,2),DATS(23) / 1979,  1, 18D0 /
      DATA (IDATE(I,24),I=1,2),DATS(24) / 1980,  1, 19D0 /
      DATA (IDATE(I,25),I=1,2),DATS(25) / 1981,  7, 20D0 /
      DATA (IDATE(I,26),I=1,2),DATS(26) / 1982,  7, 21D0 /
      DATA (IDATE(I,27),I=1,2),DATS(27) / 1983,  7, 22D0 /
      DATA (IDATE(I,28),I=1,2),DATS(28) / 1985,  7, 23D0 /
      DATA (IDATE(I,29),I=1,2),DATS(29) / 1988,  1, 24D0 /
      DATA (IDATE(I,30),I=1,2),DATS(30) / 1990,  1, 25D0 /
      DATA (IDATE(I,31),I=1,2),DATS(31) / 1991,  1, 26D0 /
      DATA (IDATE(I,32),I=1,2),DATS(32) / 1992,  7, 27D0 /
      DATA (IDATE(I,33),I=1,2),DATS(33) / 1993,  7, 28D0 /
      DATA (IDATE(I,34),I=1,2),DATS(34) / 1994,  7, 29D0 /
      DATA (IDATE(I,35),I=1,2),DATS(35) / 1996,  1, 30D0 /
      DATA (IDATE(I,36),I=1,2),DATS(36) / 1997,  7, 31D0 /
      DATA (IDATE(I,37),I=1,2),DATS(37) / 1999,  1, 32D0 /
      DATA (IDATE(I,38),I=1,2),DATS(38) / 2006,  1, 33D0 /
      DATA (IDATE(I,39),I=1,2),DATS(39) / 2009,  1, 34D0 /

*  Reference dates and drift rates
      DATA (DRIFT(I, 1),I=1,2) / 37300D0, 0.001296D0 /
      DATA (DRIFT(I, 2),I=1,2) / 37300D0, 0.001296D0 /
      DATA (DRIFT(I, 3),I=1,2) / 37300D0, 0.001296D0 /
      DATA (DRIFT(I, 4),I=1,2) / 37665D0, 0.0011232D0 /
      DATA (DRIFT(I, 5),I=1,2) / 37665D0, 0.0011232D0 /
      DATA (DRIFT(I, 6),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I, 7),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I, 8),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I, 9),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I,10),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I,11),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I,12),I=1,2) / 38761D0, 0.001296D0 /
      DATA (DRIFT(I,13),I=1,2) / 39126D0, 0.002592D0 /
      DATA (DRIFT(I,14),I=1,2) / 39126D0, 0.002592D0 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Initialize the result to zero and the status to OK.
      DAT = 0D0
      JS = 0

*  If invalid fraction of a day, set error status and give up.
      IF ( FD.LT.0D0 .OR. FD.GE.1D0 ) THEN
         JS = -4
         GO TO 9000
      END IF

*  Convert the date into an MJD.
      CALL iau_CAL2JD ( IY, IM, ID, DJM0, DJM, JS )

*  If invalid year, month, or day, give up.
      IF ( JS .LT. 0 ) GO TO 9000

*  If pre-UTC year, set warning status and give up.
      IF ( IY .LT. IDATE(1,1) ) THEN
         JS = 1
         GO TO 9000
      END IF

*  If suspiciously late year, set warning status but proceed.
      IF ( IY .GT. IYV+5 ) JS = 1

*  Combine year and month.
      M = 12*IY+IM

*  Prepare to search the tables.
      MORE = .TRUE.

*  Find the most recent table entry.
      DO 1 I=NDAT,1,-1
         IF ( MORE ) THEN
            IS = I
            MORE = M .LT. ( 12*IDATE(1,I) + IDATE(2,I) )
         END IF
 1    CONTINUE

*  Get the Delta(AT).
      DAT = DATS(IS)

*  If pre-1972, adjust for drift.
      IF ( IS .LE. NERA1 ) DAT = DAT +
     :                          ( DJM + FD - DRIFT(1,IS) ) * DRIFT(2,IS)

*  Return the Delta(AT) value and the status.
 9000 CONTINUE
      DELTAT = DAT
      J = JS

*  Finished.

*+-----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  Standards Of Fundamental Astronomy Review Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Review Board ("the Board").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and 
*     restrictions listed below.
*
*  3. You (the user) may copy and adapt the SOFA software and its 
*     algorithms for your own purposes and you may copy and distribute
*     a resulting "derived work" to others on a world-wide, royalty-free 
*     basis, provided that the derived work complies with the following
*     requirements: 
*
*     a) Your work shall be marked or carry a statement that it (i) uses
*        routines and computations derived by you from software provided 
*        by SOFA under license to you; and (ii) does not contain
*        software provided by SOFA or software that has been distributed
*        by or endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon and/or differs from the
*        original SOFA software.
*
*     c) The name(s) of all routine(s) that you distribute shall differ
*        from the SOFA names, even when the SOFA content has not been
*        otherwise changed.
*
*     d) The routine-naming prefix "iau" shall not be used.
*
*     e) The origin of the SOFA components of your derived work must not
*        be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     f) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have granted 
*        a further right to modify the source code of your derived work.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall acknowledge
*     that the SOFA software was used in obtaining those results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or by
*     inappropriate modification.
*
*  6. The SOFA software is provided "as is" and the Board makes no 
*     warranty as to its use or performance.   The Board does not and 
*     cannot warrant the performance or results which the user may obtain 
*     by using the SOFA software.  The Board makes no warranties, express 
*     or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Board be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Board representative has been advised of such
*     damages, or for any claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms 
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*     Internet email: sofa@rl.ac.uk
*     Postal address: IAU SOFA Center
*                     Rutherford Appleton Laboratory
*                     Chilton, Didcot, Oxon OX11 0QX
*                     United Kingdom
*
*-----------------------------------------------------------------------

      END
      DOUBLE PRECISION FUNCTION iau_DTDB ( DATE1, DATE2,
     :                                     UT, ELONG, U, V )
*+
*  - - - - - - - - -
*   i a u _ D T D B
*  - - - - - - - - -
*
*  An approximation to TDB-TT, the difference between barycentric
*  dynamical time and terrestrial time, for an observer on the Earth.
*
*  The different time scales - proper, coordinate and realized - are
*  related to each other:
*
*            TAI             <-  physically realized
*             :
*          offset            <-  observed (nominally +32.184s)
*             :
*            TT              <-  terrestrial time
*             :
*    rate adjustment (L_G)   <-  definition of TT
*             :
*            TCG             <-  time scale for GCRS
*             :
*      "periodic" terms      <-  iau_DTDB is an implementation
*             :
*    rate adjustment (L_C)   <-  function of solar-system ephemeris
*             :
*            TCB             <-  time scale for BCRS
*             :
*    rate adjustment (-L_B)  <-  definition of TDB
*             :
*            TDB             <-  TCB scaled to track TT
*             :
*      "periodic" terms      <-  -iau_DTDB is an approximation
*             :
*            TT              <-  terrestrial time
*
*  Adopted values for the various constants can be found in the IERS
*  Conventions (McCarthy & Petit 2003).
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     DATE1,DATE2     d    date, TDB (Notes 1-3)
*     UT              d    universal time (UT1, fraction of one day)
*     ELONG           d    longitude (east positive, radians)
*     U               d    distance from Earth spin axis (km)
*     V               d    distance north of equatorial plane (km)
*
*  Returned:
*    iau_DTDB         d    TDB-TT (seconds)
*
*  Notes:
*
*  1) The date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the arguments DATE1 and DATE2.  For
*     example, JD(TDB)=2450123.7 could be expressed in any of these
*     ways, among others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in cases
*     where the loss of several decimal digits of resolution is
*     acceptable.  The J2000 method is best matched to the way the
*     argument is handled internally and will deliver the optimum
*     resolution.  The MJD method and the date & time methods are both
*     good compromises between resolution and convenience.
*
*     Although the date is, formally, barycentric dynamical time (TDB),
*     the terrestrial dynamical time (TT) can be used with no practical
*     effect on the accuracy of the prediction.
*
*  2) TT can be regarded as a coordinate time that is realized as an
*     offset of 32.184s from International Atomic Time, TAI.  TT is a
*     specific linear transformation of geocentric coordinate time TCG,
*     which is the time scale for the Geocentric Celestial Reference
*     System, GCRS.
*
*  3) TDB is a coordinate time, and is a specific linear transformation
*     of barycentric coordinate time TCB, which is the time scale for
*     the Barycentric Celestial Reference System, BCRS.
*
*  4) The difference TCG-TCB depends on the masses and positions of the
*     bodies of the solar system and the velocity of the Earth.  It is
*     dominated by a rate difference, the residual being of a periodic
*     character.  The latter, which is modeled by the present routine,
*     comprises a main (annual) sinusoidal term of amplitude
*     approximately 0.00166 seconds, plus planetary terms up to about
*     20 microseconds, and lunar and diurnal terms up to 2 microseconds.
*     These effects come from the changing transverse Doppler effect
*     and gravitational red-shift as the observer (on the Earth's
*     surface) experiences variations in speed (with respect to the
*     BCRS) and gravitational potential.
*
*  5) TDB can be regarded as the same as TCB but with a rate adjustment
*     to keep it close to TT, which is convenient for many applications.
*     The history of successive attempts to define TDB is set out in
*     Resolution 3 adopted by the IAU General Assembly in 2006, which
*     defines a fixed TDB(TCB) transformation that is consistent with
*     contemporary solar-system ephemerides.  Future ephemerides will
*     imply slightly changed transformations between TCG and TCB, which
*     could introduce a linear drift between TDB and TT;  however, any
*     such drift is unlikely to exceed 1 nanosecond per century.
*
*  6) The geocentric TDB-TT model used in the present routine is that of
*     Fairhead & Bretagnon (1990), in its full form.  It was originally
*     supplied by Fairhead (private communications with P.T.Wallace,
*     1990) as a Fortran subroutine.  The present routine contains an
*     adaptation of the Fairhead code.  The numerical results are
*     essentially unaffected by the changes, the differences with
*     respect to the Fairhead & Bretagnon original being at the 1D-20 s
*     level.
*
*     The topocentric part of the model is from Moyer (1981) and
*     Murray (1983), with fundamental arguments adapted from
*     Simon et al. 1994.  It is an approximation to the expression
*     ( v / c ) . ( r / c ), where v is the barycentric velocity of
*     the Earth, r is the geocentric position of the observer and
*     c is the speed of light.
*
*     By supplying zeroes for U and V, the topocentric part of the
*     model can be nullified, and the routine will return the Fairhead
*     & Bretagnon result alone.
*
*  7) During the interval 1950-2050, the absolute accuracy is better
*     than +/- 3 nanoseconds relative to time ephemerides obtained by
*     direct numerical integrations based on the JPL DE405 solar system
*     ephemeris.
*
*  8) It must be stressed that the present routine is merely a model,
*     and that numerical integration of solar-system ephemerides is the
*     definitive method for predicting the relationship between TCG and
*     TCB and hence between TT and TDB.
*
*  References:
*
*     Fairhead, L., & Bretagnon, P., Astron.Astrophys., 229, 240-247
*     (1990).
*
*     IAU 2006 Resolution 3.
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Moyer, T.D., Cel.Mech., 23, 33 (1981).
*
*     Murray, C.A., Vectorial Astrometry, Adam Hilger (1983).
*
*     Seidelmann, P.K. et al., Explanatory Supplement to the
*     Astronomical Almanac, Chapter 2, University Science Books (1992).
*
*     Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G. & Laskar, J., Astron.Astrophys., 282, 663-683 (1994).
*
*  This revision:  2008 May 24
*
*  Copyright (C) 2008 IAU SOFA Review Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, UT, ELONG, U, V

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Degrees to radians
      DOUBLE PRECISION DD2R
      PARAMETER ( DD2R = 1.745329251994329576923691D-2 )

*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian millennium
      DOUBLE PRECISION DJM
      PARAMETER ( DJM = 365250D0 )

      DOUBLE PRECISION T, TSOL, W, ELSUN, EMSUN, D, ELJ, ELS,
     :                 WT, W0, W1, W2, W3, W4, WF, WJ

*
*  =====================
*  Fairhead et al. model
*  =====================
*
*  787 sets of three coefficients.
*
*  Each set is amplitude (microseconds)
*              frequency (radians per Julian millennium since J2000),
*              phase (radians).
*
*  Sets   1-474 are the T**0 terms,
*   "   475-679  "   "  T**1   "
*   "   680-764  "   "  T**2   "
*   "   765-784  "   "  T**3   "
*   "   785-787  "   "  T**4   "  .
*
      DOUBLE PRECISION FAIRHD(3,787)
      INTEGER I,J
      DATA ((FAIRHD(I,J),I=1,3),J=  1, 10) /
     : 1656.674564D-6,    6283.075849991D0, 6.240054195D0,
     :   22.417471D-6,    5753.384884897D0, 4.296977442D0,
     :   13.839792D-6,   12566.151699983D0, 6.196904410D0,
     :    4.770086D-6,     529.690965095D0, 0.444401603D0,
     :    4.676740D-6,    6069.776754553D0, 4.021195093D0,
     :    2.256707D-6,     213.299095438D0, 5.543113262D0,
     :    1.694205D-6,      -3.523118349D0, 5.025132748D0,
     :    1.554905D-6,   77713.771467920D0, 5.198467090D0,
     :    1.276839D-6,    7860.419392439D0, 5.988822341D0,
     :    1.193379D-6,    5223.693919802D0, 3.649823730D0 /
      DATA ((FAIRHD(I,J),I=1,3),J= 11, 20) /
     :    1.115322D-6,    3930.209696220D0, 1.422745069D0,
     :    0.794185D-6,   11506.769769794D0, 2.322313077D0,
     :    0.447061D-6,      26.298319800D0, 3.615796498D0,
     :    0.435206D-6,    -398.149003408D0, 4.349338347D0,
     :    0.600309D-6,    1577.343542448D0, 2.678271909D0,
     :    0.496817D-6,    6208.294251424D0, 5.696701824D0,
     :    0.486306D-6,    5884.926846583D0, 0.520007179D0,
     :    0.432392D-6,      74.781598567D0, 2.435898309D0,
     :    0.468597D-6,    6244.942814354D0, 5.866398759D0,
     :    0.375510D-6,    5507.553238667D0, 4.103476804D0 /
      DATA ((FAIRHD(I,J),I=1,3),J= 21, 30) /
     :    0.243085D-6,    -775.522611324D0, 3.651837925D0,
     :    0.173435D-6,   18849.227549974D0, 6.153743485D0,
     :    0.230685D-6,    5856.477659115D0, 4.773852582D0,
     :    0.203747D-6,   12036.460734888D0, 4.333987818D0,
     :    0.143935D-6,    -796.298006816D0, 5.957517795D0,
     :    0.159080D-6,   10977.078804699D0, 1.890075226D0,
     :    0.119979D-6,      38.133035638D0, 4.551585768D0,
     :    0.118971D-6,    5486.777843175D0, 1.914547226D0,
     :    0.116120D-6,    1059.381930189D0, 0.873504123D0,
     :    0.137927D-6,   11790.629088659D0, 1.135934669D0 /
      DATA ((FAIRHD(I,J),I=1,3),J= 31, 40) /
     :    0.098358D-6,    2544.314419883D0, 0.092793886D0,
     :    0.101868D-6,   -5573.142801634D0, 5.984503847D0,
     :    0.080164D-6,     206.185548437D0, 2.095377709D0,
     :    0.079645D-6,    4694.002954708D0, 2.949233637D0,
     :    0.062617D-6,      20.775395492D0, 2.654394814D0,
     :    0.075019D-6,    2942.463423292D0, 4.980931759D0,
     :    0.064397D-6,    5746.271337896D0, 1.280308748D0,
     :    0.063814D-6,    5760.498431898D0, 4.167901731D0,
     :    0.048042D-6,    2146.165416475D0, 1.495846011D0,
     :    0.048373D-6,     155.420399434D0, 2.251573730D0 /
      DATA ((FAIRHD(I,J),I=1,3),J= 41, 50) /
     :    0.058844D-6,     426.598190876D0, 4.839650148D0,
     :    0.046551D-6,      -0.980321068D0, 0.921573539D0,
     :    0.054139D-6,   17260.154654690D0, 3.411091093D0,
     :    0.042411D-6,    6275.962302991D0, 2.869567043D0,
     :    0.040184D-6,      -7.113547001D0, 3.565975565D0,
     :    0.036564D-6,    5088.628839767D0, 3.324679049D0,
     :    0.040759D-6,   12352.852604545D0, 3.981496998D0,
     :    0.036507D-6,     801.820931124D0, 6.248866009D0,
     :    0.036955D-6,    3154.687084896D0, 5.071801441D0,
     :    0.042732D-6,     632.783739313D0, 5.720622217D0 /
      DATA ((FAIRHD(I,J),I=1,3),J= 51, 60) /
     :    0.042560D-6,  161000.685737473D0, 1.270837679D0,
     :    0.040480D-6,   15720.838784878D0, 2.546610123D0,
     :    0.028244D-6,   -6286.598968340D0, 5.069663519D0,
     :    0.033477D-6,    6062.663207553D0, 4.144987272D0,
     :    0.034867D-6,     522.577418094D0, 5.210064075D0,
     :    0.032438D-6,    6076.890301554D0, 0.749317412D0,
     :    0.030215D-6,    7084.896781115D0, 3.389610345D0,
     :    0.029247D-6,  -71430.695617928D0, 4.183178762D0,
     :    0.033529D-6,    9437.762934887D0, 2.404714239D0,
     :    0.032423D-6,    8827.390269875D0, 5.541473556D0 /
      DATA ((FAIRHD(I,J),I=1,3),J= 61, 70) /
     :    0.027567D-6,    6279.552731642D0, 5.040846034D0,
     :    0.029862D-6,   12139.553509107D0, 1.770181024D0,
     :    0.022509D-6,   10447.387839604D0, 1.460726241D0,
     :    0.020937D-6,    8429.241266467D0, 0.652303414D0,
     :    0.020322D-6,     419.484643875D0, 3.735430632D0,
     :    0.024816D-6,   -1194.447010225D0, 1.087136918D0,
     :    0.025196D-6,    1748.016413067D0, 2.901883301D0,
     :    0.021691D-6,   14143.495242431D0, 5.952658009D0,
     :    0.017673D-6,    6812.766815086D0, 3.186129845D0,
     :    0.022567D-6,    6133.512652857D0, 3.307984806D0 /
      DATA ((FAIRHD(I,J),I=1,3),J= 71, 80) /
     :    0.016155D-6,   10213.285546211D0, 1.331103168D0,
     :    0.014751D-6,    1349.867409659D0, 4.308933301D0,
     :    0.015949D-6,    -220.412642439D0, 4.005298270D0,
     :    0.015974D-6,   -2352.866153772D0, 6.145309371D0,
     :    0.014223D-6,   17789.845619785D0, 2.104551349D0,
     :    0.017806D-6,      73.297125859D0, 3.475975097D0,
     :    0.013671D-6,    -536.804512095D0, 5.971672571D0,
     :    0.011942D-6,    8031.092263058D0, 2.053414715D0,
     :    0.014318D-6,   16730.463689596D0, 3.016058075D0,
     :    0.012462D-6,     103.092774219D0, 1.737438797D0 /
      DATA ((FAIRHD(I,J),I=1,3),J= 81, 90) /
     :    0.010962D-6,       3.590428652D0, 2.196567739D0,
     :    0.015078D-6,   19651.048481098D0, 3.969480770D0,
     :    0.010396D-6,     951.718406251D0, 5.717799605D0,
     :    0.011707D-6,   -4705.732307544D0, 2.654125618D0,
     :    0.010453D-6,    5863.591206116D0, 1.913704550D0,
     :    0.012420D-6,    4690.479836359D0, 4.734090399D0,
     :    0.011847D-6,    5643.178563677D0, 5.489005403D0,
     :    0.008610D-6,    3340.612426700D0, 3.661698944D0,
     :    0.011622D-6,    5120.601145584D0, 4.863931876D0,
     :    0.010825D-6,     553.569402842D0, 0.842715011D0 /
      DATA ((FAIRHD(I,J),I=1,3),J= 91,100) /
     :    0.008666D-6,    -135.065080035D0, 3.293406547D0,
     :    0.009963D-6,     149.563197135D0, 4.870690598D0,
     :    0.009858D-6,    6309.374169791D0, 1.061816410D0,
     :    0.007959D-6,     316.391869657D0, 2.465042647D0,
     :    0.010099D-6,     283.859318865D0, 1.942176992D0,
     :    0.007147D-6,    -242.728603974D0, 3.661486981D0,
     :    0.007505D-6,    5230.807466803D0, 4.920937029D0,
     :    0.008323D-6,   11769.853693166D0, 1.229392026D0,
     :    0.007490D-6,   -6256.777530192D0, 3.658444681D0,
     :    0.009370D-6,  149854.400134205D0, 0.673880395D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=101,110) /
     :    0.007117D-6,      38.027672636D0, 5.294249518D0,
     :    0.007857D-6,   12168.002696575D0, 0.525733528D0,
     :    0.007019D-6,    6206.809778716D0, 0.837688810D0,
     :    0.006056D-6,     955.599741609D0, 4.194535082D0,
     :    0.008107D-6,   13367.972631107D0, 3.793235253D0,
     :    0.006731D-6,    5650.292110678D0, 5.639906583D0,
     :    0.007332D-6,      36.648562930D0, 0.114858677D0,
     :    0.006366D-6,    4164.311989613D0, 2.262081818D0,
     :    0.006858D-6,    5216.580372801D0, 0.642063318D0,
     :    0.006919D-6,    6681.224853400D0, 6.018501522D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=111,120) /
     :    0.006826D-6,    7632.943259650D0, 3.458654112D0,
     :    0.005308D-6,   -1592.596013633D0, 2.500382359D0,
     :    0.005096D-6,   11371.704689758D0, 2.547107806D0,
     :    0.004841D-6,    5333.900241022D0, 0.437078094D0,
     :    0.005582D-6,    5966.683980335D0, 2.246174308D0,
     :    0.006304D-6,   11926.254413669D0, 2.512929171D0,
     :    0.006603D-6,   23581.258177318D0, 5.393136889D0,
     :    0.005123D-6,      -1.484472708D0, 2.999641028D0,
     :    0.004648D-6,    1589.072895284D0, 1.275847090D0,
     :    0.005119D-6,    6438.496249426D0, 1.486539246D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=121,130) /
     :    0.004521D-6,    4292.330832950D0, 6.140635794D0,
     :    0.005680D-6,   23013.539539587D0, 4.557814849D0,
     :    0.005488D-6,      -3.455808046D0, 0.090675389D0,
     :    0.004193D-6,    7234.794256242D0, 4.869091389D0,
     :    0.003742D-6,    7238.675591600D0, 4.691976180D0,
     :    0.004148D-6,    -110.206321219D0, 3.016173439D0,
     :    0.004553D-6,   11499.656222793D0, 5.554998314D0,
     :    0.004892D-6,    5436.993015240D0, 1.475415597D0,
     :    0.004044D-6,    4732.030627343D0, 1.398784824D0,
     :    0.004164D-6,   12491.370101415D0, 5.650931916D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=131,140) /
     :    0.004349D-6,   11513.883316794D0, 2.181745369D0,
     :    0.003919D-6,   12528.018664345D0, 5.823319737D0,
     :    0.003129D-6,    6836.645252834D0, 0.003844094D0,
     :    0.004080D-6,   -7058.598461315D0, 3.690360123D0,
     :    0.003270D-6,      76.266071276D0, 1.517189902D0,
     :    0.002954D-6,    6283.143160294D0, 4.447203799D0,
     :    0.002872D-6,      28.449187468D0, 1.158692983D0,
     :    0.002881D-6,     735.876513532D0, 0.349250250D0,
     :    0.003279D-6,    5849.364112115D0, 4.893384368D0,
     :    0.003625D-6,    6209.778724132D0, 1.473760578D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=141,150) /
     :    0.003074D-6,     949.175608970D0, 5.185878737D0,
     :    0.002775D-6,    9917.696874510D0, 1.030026325D0,
     :    0.002646D-6,   10973.555686350D0, 3.918259169D0,
     :    0.002575D-6,   25132.303399966D0, 6.109659023D0,
     :    0.003500D-6,     263.083923373D0, 1.892100742D0,
     :    0.002740D-6,   18319.536584880D0, 4.320519510D0,
     :    0.002464D-6,     202.253395174D0, 4.698203059D0,
     :    0.002409D-6,       2.542797281D0, 5.325009315D0,
     :    0.003354D-6,  -90955.551694697D0, 1.942656623D0,
     :    0.002296D-6,    6496.374945429D0, 5.061810696D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=151,160) /
     :    0.003002D-6,    6172.869528772D0, 2.797822767D0,
     :    0.003202D-6,   27511.467873537D0, 0.531673101D0,
     :    0.002954D-6,   -6283.008539689D0, 4.533471191D0,
     :    0.002353D-6,     639.897286314D0, 3.734548088D0,
     :    0.002401D-6,   16200.772724501D0, 2.605547070D0,
     :    0.003053D-6,  233141.314403759D0, 3.029030662D0,
     :    0.003024D-6,   83286.914269554D0, 2.355556099D0,
     :    0.002863D-6,   17298.182327326D0, 5.240963796D0,
     :    0.002103D-6,   -7079.373856808D0, 5.756641637D0,
     :    0.002303D-6,   83996.847317911D0, 2.013686814D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=161,170) /
     :    0.002303D-6,   18073.704938650D0, 1.089100410D0,
     :    0.002381D-6,      63.735898303D0, 0.759188178D0,
     :    0.002493D-6,    6386.168624210D0, 0.645026535D0,
     :    0.002366D-6,       3.932153263D0, 6.215885448D0,
     :    0.002169D-6,   11015.106477335D0, 4.845297676D0,
     :    0.002397D-6,    6243.458341645D0, 3.809290043D0,
     :    0.002183D-6,    1162.474704408D0, 6.179611691D0,
     :    0.002353D-6,    6246.427287062D0, 4.781719760D0,
     :    0.002199D-6,    -245.831646229D0, 5.956152284D0,
     :    0.001729D-6,    3894.181829542D0, 1.264976635D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=171,180) /
     :    0.001896D-6,   -3128.388765096D0, 4.914231596D0,
     :    0.002085D-6,      35.164090221D0, 1.405158503D0,
     :    0.002024D-6,   14712.317116458D0, 2.752035928D0,
     :    0.001737D-6,    6290.189396992D0, 5.280820144D0,
     :    0.002229D-6,     491.557929457D0, 1.571007057D0,
     :    0.001602D-6,   14314.168113050D0, 4.203664806D0,
     :    0.002186D-6,     454.909366527D0, 1.402101526D0,
     :    0.001897D-6,   22483.848574493D0, 4.167932508D0,
     :    0.001825D-6,   -3738.761430108D0, 0.545828785D0,
     :    0.001894D-6,    1052.268383188D0, 5.817167450D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=181,190) /
     :    0.001421D-6,      20.355319399D0, 2.419886601D0,
     :    0.001408D-6,   10984.192351700D0, 2.732084787D0,
     :    0.001847D-6,   10873.986030480D0, 2.903477885D0,
     :    0.001391D-6,   -8635.942003763D0, 0.593891500D0,
     :    0.001388D-6,      -7.046236698D0, 1.166145902D0,
     :    0.001810D-6,  -88860.057071188D0, 0.487355242D0,
     :    0.001288D-6,   -1990.745017041D0, 3.913022880D0,
     :    0.001297D-6,   23543.230504682D0, 3.063805171D0,
     :    0.001335D-6,    -266.607041722D0, 3.995764039D0,
     :    0.001376D-6,   10969.965257698D0, 5.152914309D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=191,200) /
     :    0.001745D-6,  244287.600007027D0, 3.626395673D0,
     :    0.001649D-6,   31441.677569757D0, 1.952049260D0,
     :    0.001416D-6,    9225.539273283D0, 4.996408389D0,
     :    0.001238D-6,    4804.209275927D0, 5.503379738D0,
     :    0.001472D-6,    4590.910180489D0, 4.164913291D0,
     :    0.001169D-6,    6040.347246017D0, 5.841719038D0,
     :    0.001039D-6,    5540.085789459D0, 2.769753519D0,
     :    0.001004D-6,    -170.672870619D0, 0.755008103D0,
     :    0.001284D-6,   10575.406682942D0, 5.306538209D0,
     :    0.001278D-6,      71.812653151D0, 4.713486491D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=201,210) /
     :    0.001321D-6,   18209.330263660D0, 2.624866359D0,
     :    0.001297D-6,   21228.392023546D0, 0.382603541D0,
     :    0.000954D-6,    6282.095528923D0, 0.882213514D0,
     :    0.001145D-6,    6058.731054289D0, 1.169483931D0,
     :    0.000979D-6,    5547.199336460D0, 5.448375984D0,
     :    0.000987D-6,   -6262.300454499D0, 2.656486959D0,
     :    0.001070D-6, -154717.609887482D0, 1.827624012D0,
     :    0.000991D-6,    4701.116501708D0, 4.387001801D0,
     :    0.001155D-6,     -14.227094002D0, 3.042700750D0,
     :    0.001176D-6,     277.034993741D0, 3.335519004D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=211,220) /
     :    0.000890D-6,   13916.019109642D0, 5.601498297D0,
     :    0.000884D-6,   -1551.045222648D0, 1.088831705D0,
     :    0.000876D-6,    5017.508371365D0, 3.969902609D0,
     :    0.000806D-6,   15110.466119866D0, 5.142876744D0,
     :    0.000773D-6,   -4136.910433516D0, 0.022067765D0,
     :    0.001077D-6,     175.166059800D0, 1.844913056D0,
     :    0.000954D-6,   -6284.056171060D0, 0.968480906D0,
     :    0.000737D-6,    5326.786694021D0, 4.923831588D0,
     :    0.000845D-6,    -433.711737877D0, 4.749245231D0,
     :    0.000819D-6,    8662.240323563D0, 5.991247817D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=221,230) /
     :    0.000852D-6,     199.072001436D0, 2.189604979D0,
     :    0.000723D-6,   17256.631536341D0, 6.068719637D0,
     :    0.000940D-6,    6037.244203762D0, 6.197428148D0,
     :    0.000885D-6,   11712.955318231D0, 3.280414875D0,
     :    0.000706D-6,   12559.038152982D0, 2.824848947D0,
     :    0.000732D-6,    2379.164473572D0, 2.501813417D0,
     :    0.000764D-6,   -6127.655450557D0, 2.236346329D0,
     :    0.000908D-6,     131.541961686D0, 2.521257490D0,
     :    0.000907D-6,   35371.887265976D0, 3.370195967D0,
     :    0.000673D-6,    1066.495477190D0, 3.876512374D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=231,240) /
     :    0.000814D-6,   17654.780539750D0, 4.627122566D0,
     :    0.000630D-6,      36.027866677D0, 0.156368499D0,
     :    0.000798D-6,     515.463871093D0, 5.151962502D0,
     :    0.000798D-6,     148.078724426D0, 5.909225055D0,
     :    0.000806D-6,     309.278322656D0, 6.054064447D0,
     :    0.000607D-6,     -39.617508346D0, 2.839021623D0,
     :    0.000601D-6,     412.371096874D0, 3.984225404D0,
     :    0.000646D-6,   11403.676995575D0, 3.852959484D0,
     :    0.000704D-6,   13521.751441591D0, 2.300991267D0,
     :    0.000603D-6,  -65147.619767937D0, 4.140083146D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=241,250) /
     :    0.000609D-6,   10177.257679534D0, 0.437122327D0,
     :    0.000631D-6,    5767.611978898D0, 4.026532329D0,
     :    0.000576D-6,   11087.285125918D0, 4.760293101D0,
     :    0.000674D-6,   14945.316173554D0, 6.270510511D0,
     :    0.000726D-6,    5429.879468239D0, 6.039606892D0,
     :    0.000710D-6,   28766.924424484D0, 5.672617711D0,
     :    0.000647D-6,   11856.218651625D0, 3.397132627D0,
     :    0.000678D-6,   -5481.254918868D0, 6.249666675D0,
     :    0.000618D-6,   22003.914634870D0, 2.466427018D0,
     :    0.000738D-6,    6134.997125565D0, 2.242668890D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=251,260) /
     :    0.000660D-6,     625.670192312D0, 5.864091907D0,
     :    0.000694D-6,    3496.032826134D0, 2.668309141D0,
     :    0.000531D-6,    6489.261398429D0, 1.681888780D0,
     :    0.000611D-6, -143571.324284214D0, 2.424978312D0,
     :    0.000575D-6,   12043.574281889D0, 4.216492400D0,
     :    0.000553D-6,   12416.588502848D0, 4.772158039D0,
     :    0.000689D-6,    4686.889407707D0, 6.224271088D0,
     :    0.000495D-6,    7342.457780181D0, 3.817285811D0,
     :    0.000567D-6,    3634.621024518D0, 1.649264690D0,
     :    0.000515D-6,   18635.928454536D0, 3.945345892D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=261,270) /
     :    0.000486D-6,    -323.505416657D0, 4.061673868D0,
     :    0.000662D-6,   25158.601719765D0, 1.794058369D0,
     :    0.000509D-6,     846.082834751D0, 3.053874588D0,
     :    0.000472D-6,  -12569.674818332D0, 5.112133338D0,
     :    0.000461D-6,    6179.983075773D0, 0.513669325D0,
     :    0.000641D-6,   83467.156352816D0, 3.210727723D0,
     :    0.000520D-6,   10344.295065386D0, 2.445597761D0,
     :    0.000493D-6,   18422.629359098D0, 1.676939306D0,
     :    0.000478D-6,    1265.567478626D0, 5.487314569D0,
     :    0.000472D-6,     -18.159247265D0, 1.999707589D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=271,280) /
     :    0.000559D-6,   11190.377900137D0, 5.783236356D0,
     :    0.000494D-6,    9623.688276691D0, 3.022645053D0,
     :    0.000463D-6,    5739.157790895D0, 1.411223013D0,
     :    0.000432D-6,   16858.482532933D0, 1.179256434D0,
     :    0.000574D-6,   72140.628666286D0, 1.758191830D0,
     :    0.000484D-6,   17267.268201691D0, 3.290589143D0,
     :    0.000550D-6,    4907.302050146D0, 0.864024298D0,
     :    0.000399D-6,      14.977853527D0, 2.094441910D0,
     :    0.000491D-6,     224.344795702D0, 0.878372791D0,
     :    0.000432D-6,   20426.571092422D0, 6.003829241D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=281,290) /
     :    0.000481D-6,    5749.452731634D0, 4.309591964D0,
     :    0.000480D-6,    5757.317038160D0, 1.142348571D0,
     :    0.000485D-6,    6702.560493867D0, 0.210580917D0,
     :    0.000426D-6,    6055.549660552D0, 4.274476529D0,
     :    0.000480D-6,    5959.570433334D0, 5.031351030D0,
     :    0.000466D-6,   12562.628581634D0, 4.959581597D0,
     :    0.000520D-6,   39302.096962196D0, 4.788002889D0,
     :    0.000458D-6,   12132.439962106D0, 1.880103788D0,
     :    0.000470D-6,   12029.347187887D0, 1.405611197D0,
     :    0.000416D-6,   -7477.522860216D0, 1.082356330D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=291,300) /
     :    0.000449D-6,   11609.862544012D0, 4.179989585D0,
     :    0.000465D-6,   17253.041107690D0, 0.353496295D0,
     :    0.000362D-6,   -4535.059436924D0, 1.583849576D0,
     :    0.000383D-6,   21954.157609398D0, 3.747376371D0,
     :    0.000389D-6,      17.252277143D0, 1.395753179D0,
     :    0.000331D-6,   18052.929543158D0, 0.566790582D0,
     :    0.000430D-6,   13517.870106233D0, 0.685827538D0,
     :    0.000368D-6,   -5756.908003246D0, 0.731374317D0,
     :    0.000330D-6,   10557.594160824D0, 3.710043680D0,
     :    0.000332D-6,   20199.094959633D0, 1.652901407D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=301,310) /
     :    0.000384D-6,   11933.367960670D0, 5.827781531D0,
     :    0.000387D-6,   10454.501386605D0, 2.541182564D0,
     :    0.000325D-6,   15671.081759407D0, 2.178850542D0,
     :    0.000318D-6,     138.517496871D0, 2.253253037D0,
     :    0.000305D-6,    9388.005909415D0, 0.578340206D0,
     :    0.000352D-6,    5749.861766548D0, 3.000297967D0,
     :    0.000311D-6,    6915.859589305D0, 1.693574249D0,
     :    0.000297D-6,   24072.921469776D0, 1.997249392D0,
     :    0.000363D-6,    -640.877607382D0, 5.071820966D0,
     :    0.000323D-6,   12592.450019783D0, 1.072262823D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=311,320) /
     :    0.000341D-6,   12146.667056108D0, 4.700657997D0,
     :    0.000290D-6,    9779.108676125D0, 1.812320441D0,
     :    0.000342D-6,    6132.028180148D0, 4.322238614D0,
     :    0.000329D-6,    6268.848755990D0, 3.033827743D0,
     :    0.000374D-6,   17996.031168222D0, 3.388716544D0,
     :    0.000285D-6,    -533.214083444D0, 4.687313233D0,
     :    0.000338D-6,    6065.844601290D0, 0.877776108D0,
     :    0.000276D-6,      24.298513841D0, 0.770299429D0,
     :    0.000336D-6,   -2388.894020449D0, 5.353796034D0,
     :    0.000290D-6,    3097.883822726D0, 4.075291557D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=321,330) /
     :    0.000318D-6,     709.933048357D0, 5.941207518D0,
     :    0.000271D-6,   13095.842665077D0, 3.208912203D0,
     :    0.000331D-6,    6073.708907816D0, 4.007881169D0,
     :    0.000292D-6,     742.990060533D0, 2.714333592D0,
     :    0.000362D-6,   29088.811415985D0, 3.215977013D0,
     :    0.000280D-6,   12359.966151546D0, 0.710872502D0,
     :    0.000267D-6,   10440.274292604D0, 4.730108488D0,
     :    0.000262D-6,     838.969287750D0, 1.327720272D0,
     :    0.000250D-6,   16496.361396202D0, 0.898769761D0,
     :    0.000325D-6,   20597.243963041D0, 0.180044365D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=331,340) /
     :    0.000268D-6,    6148.010769956D0, 5.152666276D0,
     :    0.000284D-6,    5636.065016677D0, 5.655385808D0,
     :    0.000301D-6,    6080.822454817D0, 2.135396205D0,
     :    0.000294D-6,    -377.373607916D0, 3.708784168D0,
     :    0.000236D-6,    2118.763860378D0, 1.733578756D0,
     :    0.000234D-6,    5867.523359379D0, 5.575209112D0,
     :    0.000268D-6, -226858.238553767D0, 0.069432392D0,
     :    0.000265D-6,  167283.761587465D0, 4.369302826D0,
     :    0.000280D-6,   28237.233459389D0, 5.304829118D0,
     :    0.000292D-6,   12345.739057544D0, 4.096094132D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=341,350) /
     :    0.000223D-6,   19800.945956225D0, 3.069327406D0,
     :    0.000301D-6,   43232.306658416D0, 6.205311188D0,
     :    0.000264D-6,   18875.525869774D0, 1.417263408D0,
     :    0.000304D-6,   -1823.175188677D0, 3.409035232D0,
     :    0.000301D-6,     109.945688789D0, 0.510922054D0,
     :    0.000260D-6,     813.550283960D0, 2.389438934D0,
     :    0.000299D-6,  316428.228673312D0, 5.384595078D0,
     :    0.000211D-6,    5756.566278634D0, 3.789392838D0,
     :    0.000209D-6,    5750.203491159D0, 1.661943545D0,
     :    0.000240D-6,   12489.885628707D0, 5.684549045D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=351,360) /
     :    0.000216D-6,    6303.851245484D0, 3.862942261D0,
     :    0.000203D-6,    1581.959348283D0, 5.549853589D0,
     :    0.000200D-6,    5642.198242609D0, 1.016115785D0,
     :    0.000197D-6,     -70.849445304D0, 4.690702525D0,
     :    0.000227D-6,    6287.008003254D0, 2.911891613D0,
     :    0.000197D-6,     533.623118358D0, 1.048982898D0,
     :    0.000205D-6,   -6279.485421340D0, 1.829362730D0,
     :    0.000209D-6,  -10988.808157535D0, 2.636140084D0,
     :    0.000208D-6,    -227.526189440D0, 4.127883842D0,
     :    0.000191D-6,     415.552490612D0, 4.401165650D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=361,370) /
     :    0.000190D-6,   29296.615389579D0, 4.175658539D0,
     :    0.000264D-6,   66567.485864652D0, 4.601102551D0,
     :    0.000256D-6,   -3646.350377354D0, 0.506364778D0,
     :    0.000188D-6,   13119.721102825D0, 2.032195842D0,
     :    0.000185D-6,    -209.366942175D0, 4.694756586D0,
     :    0.000198D-6,   25934.124331089D0, 3.832703118D0,
     :    0.000195D-6,    4061.219215394D0, 3.308463427D0,
     :    0.000234D-6,    5113.487598583D0, 1.716090661D0,
     :    0.000188D-6,    1478.866574064D0, 5.686865780D0,
     :    0.000222D-6,   11823.161639450D0, 1.942386641D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=371,380) /
     :    0.000181D-6,   10770.893256262D0, 1.999482059D0,
     :    0.000171D-6,    6546.159773364D0, 1.182807992D0,
     :    0.000206D-6,      70.328180442D0, 5.934076062D0,
     :    0.000169D-6,   20995.392966449D0, 2.169080622D0,
     :    0.000191D-6,   10660.686935042D0, 5.405515999D0,
     :    0.000228D-6,   33019.021112205D0, 4.656985514D0,
     :    0.000184D-6,   -4933.208440333D0, 3.327476868D0,
     :    0.000220D-6,    -135.625325010D0, 1.765430262D0,
     :    0.000166D-6,   23141.558382925D0, 3.454132746D0,
     :    0.000191D-6,    6144.558353121D0, 5.020393445D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=381,390) /
     :    0.000180D-6,    6084.003848555D0, 0.602182191D0,
     :    0.000163D-6,   17782.732072784D0, 4.960593133D0,
     :    0.000225D-6,   16460.333529525D0, 2.596451817D0,
     :    0.000222D-6,    5905.702242076D0, 3.731990323D0,
     :    0.000204D-6,     227.476132789D0, 5.636192701D0,
     :    0.000159D-6,   16737.577236597D0, 3.600691544D0,
     :    0.000200D-6,    6805.653268085D0, 0.868220961D0,
     :    0.000187D-6,   11919.140866668D0, 2.629456641D0,
     :    0.000161D-6,     127.471796607D0, 2.862574720D0,
     :    0.000205D-6,    6286.666278643D0, 1.742882331D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=391,400) /
     :    0.000189D-6,     153.778810485D0, 4.812372643D0,
     :    0.000168D-6,   16723.350142595D0, 0.027860588D0,
     :    0.000149D-6,   11720.068865232D0, 0.659721876D0,
     :    0.000189D-6,    5237.921013804D0, 5.245313000D0,
     :    0.000143D-6,    6709.674040867D0, 4.317625647D0,
     :    0.000146D-6,    4487.817406270D0, 4.815297007D0,
     :    0.000144D-6,    -664.756045130D0, 5.381366880D0,
     :    0.000175D-6,    5127.714692584D0, 4.728443327D0,
     :    0.000162D-6,    6254.626662524D0, 1.435132069D0,
     :    0.000187D-6,   47162.516354635D0, 1.354371923D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=401,410) /
     :    0.000146D-6,   11080.171578918D0, 3.369695406D0,
     :    0.000180D-6,    -348.924420448D0, 2.490902145D0,
     :    0.000148D-6,     151.047669843D0, 3.799109588D0,
     :    0.000157D-6,    6197.248551160D0, 1.284375887D0,
     :    0.000167D-6,     146.594251718D0, 0.759969109D0,
     :    0.000133D-6,   -5331.357443741D0, 5.409701889D0,
     :    0.000154D-6,      95.979227218D0, 3.366890614D0,
     :    0.000148D-6,   -6418.140930027D0, 3.384104996D0,
     :    0.000128D-6,   -6525.804453965D0, 3.803419985D0,
     :    0.000130D-6,   11293.470674356D0, 0.939039445D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=411,420) /
     :    0.000152D-6,   -5729.506447149D0, 0.734117523D0,
     :    0.000138D-6,     210.117701700D0, 2.564216078D0,
     :    0.000123D-6,    6066.595360816D0, 4.517099537D0,
     :    0.000140D-6,   18451.078546566D0, 0.642049130D0,
     :    0.000126D-6,   11300.584221356D0, 3.485280663D0,
     :    0.000119D-6,   10027.903195729D0, 3.217431161D0,
     :    0.000151D-6,    4274.518310832D0, 4.404359108D0,
     :    0.000117D-6,    6072.958148291D0, 0.366324650D0,
     :    0.000165D-6,   -7668.637425143D0, 4.298212528D0,
     :    0.000117D-6,   -6245.048177356D0, 5.379518958D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=421,430) /
     :    0.000130D-6,   -5888.449964932D0, 4.527681115D0,
     :    0.000121D-6,    -543.918059096D0, 6.109429504D0,
     :    0.000162D-6,    9683.594581116D0, 5.720092446D0,
     :    0.000141D-6,    6219.339951688D0, 0.679068671D0,
     :    0.000118D-6,   22743.409379516D0, 4.881123092D0,
     :    0.000129D-6,    1692.165669502D0, 0.351407289D0,
     :    0.000126D-6,    5657.405657679D0, 5.146592349D0,
     :    0.000114D-6,     728.762966531D0, 0.520791814D0,
     :    0.000120D-6,      52.596639600D0, 0.948516300D0,
     :    0.000115D-6,      65.220371012D0, 3.504914846D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=431,440) /
     :    0.000126D-6,    5881.403728234D0, 5.577502482D0,
     :    0.000158D-6,  163096.180360983D0, 2.957128968D0,
     :    0.000134D-6,   12341.806904281D0, 2.598576764D0,
     :    0.000151D-6,   16627.370915377D0, 3.985702050D0,
     :    0.000109D-6,    1368.660252845D0, 0.014730471D0,
     :    0.000131D-6,    6211.263196841D0, 0.085077024D0,
     :    0.000146D-6,    5792.741760812D0, 0.708426604D0,
     :    0.000146D-6,     -77.750543984D0, 3.121576600D0,
     :    0.000107D-6,    5341.013788022D0, 0.288231904D0,
     :    0.000138D-6,    6281.591377283D0, 2.797450317D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=441,450) /
     :    0.000113D-6,   -6277.552925684D0, 2.788904128D0,
     :    0.000115D-6,    -525.758811831D0, 5.895222200D0,
     :    0.000138D-6,    6016.468808270D0, 6.096188999D0,
     :    0.000139D-6,   23539.707386333D0, 2.028195445D0,
     :    0.000146D-6,   -4176.041342449D0, 4.660008502D0,
     :    0.000107D-6,   16062.184526117D0, 4.066520001D0,
     :    0.000142D-6,   83783.548222473D0, 2.936315115D0,
     :    0.000128D-6,    9380.959672717D0, 3.223844306D0,
     :    0.000135D-6,    6205.325306007D0, 1.638054048D0,
     :    0.000101D-6,    2699.734819318D0, 5.481603249D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=451,460) /
     :    0.000104D-6,    -568.821874027D0, 2.205734493D0,
     :    0.000103D-6,    6321.103522627D0, 2.440421099D0,
     :    0.000119D-6,    6321.208885629D0, 2.547496264D0,
     :    0.000138D-6,    1975.492545856D0, 2.314608466D0,
     :    0.000121D-6,     137.033024162D0, 4.539108237D0,
     :    0.000123D-6,   19402.796952817D0, 4.538074405D0,
     :    0.000119D-6,   22805.735565994D0, 2.869040566D0,
     :    0.000133D-6,   64471.991241142D0, 6.056405489D0,
     :    0.000129D-6,     -85.827298831D0, 2.540635083D0,
     :    0.000131D-6,   13613.804277336D0, 4.005732868D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=461,470) /
     :    0.000104D-6,    9814.604100291D0, 1.959967212D0,
     :    0.000112D-6,   16097.679950283D0, 3.589026260D0,
     :    0.000123D-6,    2107.034507542D0, 1.728627253D0,
     :    0.000121D-6,   36949.230808424D0, 6.072332087D0,
     :    0.000108D-6,  -12539.853380183D0, 3.716133846D0,
     :    0.000113D-6,   -7875.671863624D0, 2.725771122D0,
     :    0.000109D-6,    4171.425536614D0, 4.033338079D0,
     :    0.000101D-6,    6247.911759770D0, 3.441347021D0,
     :    0.000113D-6,    7330.728427345D0, 0.656372122D0,
     :    0.000113D-6,   51092.726050855D0, 2.791483066D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=471,480) /
     :    0.000106D-6,    5621.842923210D0, 1.815323326D0,
     :    0.000101D-6,     111.430161497D0, 5.711033677D0,
     :    0.000103D-6,     909.818733055D0, 2.812745443D0,
     :    0.000101D-6,    1790.642637886D0, 1.965746028D0,

*  T
     :  102.156724D-6,    6283.075849991D0, 4.249032005D0,
     :    1.706807D-6,   12566.151699983D0, 4.205904248D0,
     :    0.269668D-6,     213.299095438D0, 3.400290479D0,
     :    0.265919D-6,     529.690965095D0, 5.836047367D0,
     :    0.210568D-6,      -3.523118349D0, 6.262738348D0,
     :    0.077996D-6,    5223.693919802D0, 4.670344204D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=481,490) /
     :    0.054764D-6,    1577.343542448D0, 4.534800170D0,
     :    0.059146D-6,      26.298319800D0, 1.083044735D0,
     :    0.034420D-6,    -398.149003408D0, 5.980077351D0,
     :    0.032088D-6,   18849.227549974D0, 4.162913471D0,
     :    0.033595D-6,    5507.553238667D0, 5.980162321D0,
     :    0.029198D-6,    5856.477659115D0, 0.623811863D0,
     :    0.027764D-6,     155.420399434D0, 3.745318113D0,
     :    0.025190D-6,    5746.271337896D0, 2.980330535D0,
     :    0.022997D-6,    -796.298006816D0, 1.174411803D0,
     :    0.024976D-6,    5760.498431898D0, 2.467913690D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=491,500) /
     :    0.021774D-6,     206.185548437D0, 3.854787540D0,
     :    0.017925D-6,    -775.522611324D0, 1.092065955D0,
     :    0.013794D-6,     426.598190876D0, 2.699831988D0,
     :    0.013276D-6,    6062.663207553D0, 5.845801920D0,
     :    0.011774D-6,   12036.460734888D0, 2.292832062D0,
     :    0.012869D-6,    6076.890301554D0, 5.333425680D0,
     :    0.012152D-6,    1059.381930189D0, 6.222874454D0,
     :    0.011081D-6,      -7.113547001D0, 5.154724984D0,
     :    0.010143D-6,    4694.002954708D0, 4.044013795D0,
     :    0.009357D-6,    5486.777843175D0, 3.416081409D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=501,510) /
     :    0.010084D-6,     522.577418094D0, 0.749320262D0,
     :    0.008587D-6,   10977.078804699D0, 2.777152598D0,
     :    0.008628D-6,    6275.962302991D0, 4.562060226D0,
     :    0.008158D-6,    -220.412642439D0, 5.806891533D0,
     :    0.007746D-6,    2544.314419883D0, 1.603197066D0,
     :    0.007670D-6,    2146.165416475D0, 3.000200440D0,
     :    0.007098D-6,      74.781598567D0, 0.443725817D0,
     :    0.006180D-6,    -536.804512095D0, 1.302642751D0,
     :    0.005818D-6,    5088.628839767D0, 4.827723531D0,
     :    0.004945D-6,   -6286.598968340D0, 0.268305170D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=511,520) /
     :    0.004774D-6,    1349.867409659D0, 5.808636673D0,
     :    0.004687D-6,    -242.728603974D0, 5.154890570D0,
     :    0.006089D-6,    1748.016413067D0, 4.403765209D0,
     :    0.005975D-6,   -1194.447010225D0, 2.583472591D0,
     :    0.004229D-6,     951.718406251D0, 0.931172179D0,
     :    0.005264D-6,     553.569402842D0, 2.336107252D0,
     :    0.003049D-6,    5643.178563677D0, 1.362634430D0,
     :    0.002974D-6,    6812.766815086D0, 1.583012668D0,
     :    0.003403D-6,   -2352.866153772D0, 2.552189886D0,
     :    0.003030D-6,     419.484643875D0, 5.286473844D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=521,530) /
     :    0.003210D-6,      -7.046236698D0, 1.863796539D0,
     :    0.003058D-6,    9437.762934887D0, 4.226420633D0,
     :    0.002589D-6,   12352.852604545D0, 1.991935820D0,
     :    0.002927D-6,    5216.580372801D0, 2.319951253D0,
     :    0.002425D-6,    5230.807466803D0, 3.084752833D0,
     :    0.002656D-6,    3154.687084896D0, 2.487447866D0,
     :    0.002445D-6,   10447.387839604D0, 2.347139160D0,
     :    0.002990D-6,    4690.479836359D0, 6.235872050D0,
     :    0.002890D-6,    5863.591206116D0, 0.095197563D0,
     :    0.002498D-6,    6438.496249426D0, 2.994779800D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=531,540) /
     :    0.001889D-6,    8031.092263058D0, 3.569003717D0,
     :    0.002567D-6,     801.820931124D0, 3.425611498D0,
     :    0.001803D-6,  -71430.695617928D0, 2.192295512D0,
     :    0.001782D-6,       3.932153263D0, 5.180433689D0,
     :    0.001694D-6,   -4705.732307544D0, 4.641779174D0,
     :    0.001704D-6,   -1592.596013633D0, 3.997097652D0,
     :    0.001735D-6,    5849.364112115D0, 0.417558428D0,
     :    0.001643D-6,    8429.241266467D0, 2.180619584D0,
     :    0.001680D-6,      38.133035638D0, 4.164529426D0,
     :    0.002045D-6,    7084.896781115D0, 0.526323854D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=541,550) /
     :    0.001458D-6,    4292.330832950D0, 1.356098141D0,
     :    0.001437D-6,      20.355319399D0, 3.895439360D0,
     :    0.001738D-6,    6279.552731642D0, 0.087484036D0,
     :    0.001367D-6,   14143.495242431D0, 3.987576591D0,
     :    0.001344D-6,    7234.794256242D0, 0.090454338D0,
     :    0.001438D-6,   11499.656222793D0, 0.974387904D0,
     :    0.001257D-6,    6836.645252834D0, 1.509069366D0,
     :    0.001358D-6,   11513.883316794D0, 0.495572260D0,
     :    0.001628D-6,    7632.943259650D0, 4.968445721D0,
     :    0.001169D-6,     103.092774219D0, 2.838496795D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=551,560) /
     :    0.001162D-6,    4164.311989613D0, 3.408387778D0,
     :    0.001092D-6,    6069.776754553D0, 3.617942651D0,
     :    0.001008D-6,   17789.845619785D0, 0.286350174D0,
     :    0.001008D-6,     639.897286314D0, 1.610762073D0,
     :    0.000918D-6,   10213.285546211D0, 5.532798067D0,
     :    0.001011D-6,   -6256.777530192D0, 0.661826484D0,
     :    0.000753D-6,   16730.463689596D0, 3.905030235D0,
     :    0.000737D-6,   11926.254413669D0, 4.641956361D0,
     :    0.000694D-6,    3340.612426700D0, 2.111120332D0,
     :    0.000701D-6,    3894.181829542D0, 2.760823491D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=561,570) /
     :    0.000689D-6,    -135.065080035D0, 4.768800780D0,
     :    0.000700D-6,   13367.972631107D0, 5.760439898D0,
     :    0.000664D-6,    6040.347246017D0, 1.051215840D0,
     :    0.000654D-6,    5650.292110678D0, 4.911332503D0,
     :    0.000788D-6,    6681.224853400D0, 4.699648011D0,
     :    0.000628D-6,    5333.900241022D0, 5.024608847D0,
     :    0.000755D-6,    -110.206321219D0, 4.370971253D0,
     :    0.000628D-6,    6290.189396992D0, 3.660478857D0,
     :    0.000635D-6,   25132.303399966D0, 4.121051532D0,
     :    0.000534D-6,    5966.683980335D0, 1.173284524D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=571,580) /
     :    0.000543D-6,    -433.711737877D0, 0.345585464D0,
     :    0.000517D-6,   -1990.745017041D0, 5.414571768D0,
     :    0.000504D-6,    5767.611978898D0, 2.328281115D0,
     :    0.000485D-6,    5753.384884897D0, 1.685874771D0,
     :    0.000463D-6,    7860.419392439D0, 5.297703006D0,
     :    0.000604D-6,     515.463871093D0, 0.591998446D0,
     :    0.000443D-6,   12168.002696575D0, 4.830881244D0,
     :    0.000570D-6,     199.072001436D0, 3.899190272D0,
     :    0.000465D-6,   10969.965257698D0, 0.476681802D0,
     :    0.000424D-6,   -7079.373856808D0, 1.112242763D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=581,590) /
     :    0.000427D-6,     735.876513532D0, 1.994214480D0,
     :    0.000478D-6,   -6127.655450557D0, 3.778025483D0,
     :    0.000414D-6,   10973.555686350D0, 5.441088327D0,
     :    0.000512D-6,    1589.072895284D0, 0.107123853D0,
     :    0.000378D-6,   10984.192351700D0, 0.915087231D0,
     :    0.000402D-6,   11371.704689758D0, 4.107281715D0,
     :    0.000453D-6,    9917.696874510D0, 1.917490952D0,
     :    0.000395D-6,     149.563197135D0, 2.763124165D0,
     :    0.000371D-6,    5739.157790895D0, 3.112111866D0,
     :    0.000350D-6,   11790.629088659D0, 0.440639857D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=591,600) /
     :    0.000356D-6,    6133.512652857D0, 5.444568842D0,
     :    0.000344D-6,     412.371096874D0, 5.676832684D0,
     :    0.000383D-6,     955.599741609D0, 5.559734846D0,
     :    0.000333D-6,    6496.374945429D0, 0.261537984D0,
     :    0.000340D-6,    6055.549660552D0, 5.975534987D0,
     :    0.000334D-6,    1066.495477190D0, 2.335063907D0,
     :    0.000399D-6,   11506.769769794D0, 5.321230910D0,
     :    0.000314D-6,   18319.536584880D0, 2.313312404D0,
     :    0.000424D-6,    1052.268383188D0, 1.211961766D0,
     :    0.000307D-6,      63.735898303D0, 3.169551388D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=601,610) /
     :    0.000329D-6,      29.821438149D0, 6.106912080D0,
     :    0.000357D-6,    6309.374169791D0, 4.223760346D0,
     :    0.000312D-6,   -3738.761430108D0, 2.180556645D0,
     :    0.000301D-6,     309.278322656D0, 1.499984572D0,
     :    0.000268D-6,   12043.574281889D0, 2.447520648D0,
     :    0.000257D-6,   12491.370101415D0, 3.662331761D0,
     :    0.000290D-6,     625.670192312D0, 1.272834584D0,
     :    0.000256D-6,    5429.879468239D0, 1.913426912D0,
     :    0.000339D-6,    3496.032826134D0, 4.165930011D0,
     :    0.000283D-6,    3930.209696220D0, 4.325565754D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=611,620) /
     :    0.000241D-6,   12528.018664345D0, 3.832324536D0,
     :    0.000304D-6,    4686.889407707D0, 1.612348468D0,
     :    0.000259D-6,   16200.772724501D0, 3.470173146D0,
     :    0.000238D-6,   12139.553509107D0, 1.147977842D0,
     :    0.000236D-6,    6172.869528772D0, 3.776271728D0,
     :    0.000296D-6,   -7058.598461315D0, 0.460368852D0,
     :    0.000306D-6,   10575.406682942D0, 0.554749016D0,
     :    0.000251D-6,   17298.182327326D0, 0.834332510D0,
     :    0.000290D-6,    4732.030627343D0, 4.759564091D0,
     :    0.000261D-6,    5884.926846583D0, 0.298259862D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=621,630) /
     :    0.000249D-6,    5547.199336460D0, 3.749366406D0,
     :    0.000213D-6,   11712.955318231D0, 5.415666119D0,
     :    0.000223D-6,    4701.116501708D0, 2.703203558D0,
     :    0.000268D-6,    -640.877607382D0, 0.283670793D0,
     :    0.000209D-6,    5636.065016677D0, 1.238477199D0,
     :    0.000193D-6,   10177.257679534D0, 1.943251340D0,
     :    0.000182D-6,    6283.143160294D0, 2.456157599D0,
     :    0.000184D-6,    -227.526189440D0, 5.888038582D0,
     :    0.000182D-6,   -6283.008539689D0, 0.241332086D0,
     :    0.000228D-6,   -6284.056171060D0, 2.657323816D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=631,640) /
     :    0.000166D-6,    7238.675591600D0, 5.930629110D0,
     :    0.000167D-6,    3097.883822726D0, 5.570955333D0,
     :    0.000159D-6,    -323.505416657D0, 5.786670700D0,
     :    0.000154D-6,   -4136.910433516D0, 1.517805532D0,
     :    0.000176D-6,   12029.347187887D0, 3.139266834D0,
     :    0.000167D-6,   12132.439962106D0, 3.556352289D0,
     :    0.000153D-6,     202.253395174D0, 1.463313961D0,
     :    0.000157D-6,   17267.268201691D0, 1.586837396D0,
     :    0.000142D-6,   83996.847317911D0, 0.022670115D0,
     :    0.000152D-6,   17260.154654690D0, 0.708528947D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=641,650) /
     :    0.000144D-6,    6084.003848555D0, 5.187075177D0,
     :    0.000135D-6,    5756.566278634D0, 1.993229262D0,
     :    0.000134D-6,    5750.203491159D0, 3.457197134D0,
     :    0.000144D-6,    5326.786694021D0, 6.066193291D0,
     :    0.000160D-6,   11015.106477335D0, 1.710431974D0,
     :    0.000133D-6,    3634.621024518D0, 2.836451652D0,
     :    0.000134D-6,   18073.704938650D0, 5.453106665D0,
     :    0.000134D-6,    1162.474704408D0, 5.326898811D0,
     :    0.000128D-6,    5642.198242609D0, 2.511652591D0,
     :    0.000160D-6,     632.783739313D0, 5.628785365D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=651,660) /
     :    0.000132D-6,   13916.019109642D0, 0.819294053D0,
     :    0.000122D-6,   14314.168113050D0, 5.677408071D0,
     :    0.000125D-6,   12359.966151546D0, 5.251984735D0,
     :    0.000121D-6,    5749.452731634D0, 2.210924603D0,
     :    0.000136D-6,    -245.831646229D0, 1.646502367D0,
     :    0.000120D-6,    5757.317038160D0, 3.240883049D0,
     :    0.000134D-6,   12146.667056108D0, 3.059480037D0,
     :    0.000137D-6,    6206.809778716D0, 1.867105418D0,
     :    0.000141D-6,   17253.041107690D0, 2.069217456D0,
     :    0.000129D-6,   -7477.522860216D0, 2.781469314D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=661,670) /
     :    0.000116D-6,    5540.085789459D0, 4.281176991D0,
     :    0.000116D-6,    9779.108676125D0, 3.320925381D0,
     :    0.000129D-6,    5237.921013804D0, 3.497704076D0,
     :    0.000113D-6,    5959.570433334D0, 0.983210840D0,
     :    0.000122D-6,    6282.095528923D0, 2.674938860D0,
     :    0.000140D-6,     -11.045700264D0, 4.957936982D0,
     :    0.000108D-6,   23543.230504682D0, 1.390113589D0,
     :    0.000106D-6,  -12569.674818332D0, 0.429631317D0,
     :    0.000110D-6,    -266.607041722D0, 5.501340197D0,
     :    0.000115D-6,   12559.038152982D0, 4.691456618D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=671,680) /
     :    0.000134D-6,   -2388.894020449D0, 0.577313584D0,
     :    0.000109D-6,   10440.274292604D0, 6.218148717D0,
     :    0.000102D-6,    -543.918059096D0, 1.477842615D0,
     :    0.000108D-6,   21228.392023546D0, 2.237753948D0,
     :    0.000101D-6,   -4535.059436924D0, 3.100492232D0,
     :    0.000103D-6,      76.266071276D0, 5.594294322D0,
     :    0.000104D-6,     949.175608970D0, 5.674287810D0,
     :    0.000101D-6,   13517.870106233D0, 2.196632348D0,
     :    0.000100D-6,   11933.367960670D0, 4.056084160D0,

*  T^2
     :    4.322990D-6,    6283.075849991D0, 2.642893748D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=681,690) /
     :    0.406495D-6,       0.000000000D0, 4.712388980D0,
     :    0.122605D-6,   12566.151699983D0, 2.438140634D0,
     :    0.019476D-6,     213.299095438D0, 1.642186981D0,
     :    0.016916D-6,     529.690965095D0, 4.510959344D0,
     :    0.013374D-6,      -3.523118349D0, 1.502210314D0,
     :    0.008042D-6,      26.298319800D0, 0.478549024D0,
     :    0.007824D-6,     155.420399434D0, 5.254710405D0,
     :    0.004894D-6,    5746.271337896D0, 4.683210850D0,
     :    0.004875D-6,    5760.498431898D0, 0.759507698D0,
     :    0.004416D-6,    5223.693919802D0, 6.028853166D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=691,700) /
     :    0.004088D-6,      -7.113547001D0, 0.060926389D0,
     :    0.004433D-6,   77713.771467920D0, 3.627734103D0,
     :    0.003277D-6,   18849.227549974D0, 2.327912542D0,
     :    0.002703D-6,    6062.663207553D0, 1.271941729D0,
     :    0.003435D-6,    -775.522611324D0, 0.747446224D0,
     :    0.002618D-6,    6076.890301554D0, 3.633715689D0,
     :    0.003146D-6,     206.185548437D0, 5.647874613D0,
     :    0.002544D-6,    1577.343542448D0, 6.232904270D0,
     :    0.002218D-6,    -220.412642439D0, 1.309509946D0,
     :    0.002197D-6,    5856.477659115D0, 2.407212349D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=701,710) /
     :    0.002897D-6,    5753.384884897D0, 5.863842246D0,
     :    0.001766D-6,     426.598190876D0, 0.754113147D0,
     :    0.001738D-6,    -796.298006816D0, 2.714942671D0,
     :    0.001695D-6,     522.577418094D0, 2.629369842D0,
     :    0.001584D-6,    5507.553238667D0, 1.341138229D0,
     :    0.001503D-6,    -242.728603974D0, 0.377699736D0,
     :    0.001552D-6,    -536.804512095D0, 2.904684667D0,
     :    0.001370D-6,    -398.149003408D0, 1.265599125D0,
     :    0.001889D-6,   -5573.142801634D0, 4.413514859D0,
     :    0.001722D-6,    6069.776754553D0, 2.445966339D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=711,720) /
     :    0.001124D-6,    1059.381930189D0, 5.041799657D0,
     :    0.001258D-6,     553.569402842D0, 3.849557278D0,
     :    0.000831D-6,     951.718406251D0, 2.471094709D0,
     :    0.000767D-6,    4694.002954708D0, 5.363125422D0,
     :    0.000756D-6,    1349.867409659D0, 1.046195744D0,
     :    0.000775D-6,     -11.045700264D0, 0.245548001D0,
     :    0.000597D-6,    2146.165416475D0, 4.543268798D0,
     :    0.000568D-6,    5216.580372801D0, 4.178853144D0,
     :    0.000711D-6,    1748.016413067D0, 5.934271972D0,
     :    0.000499D-6,   12036.460734888D0, 0.624434410D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=721,730) /
     :    0.000671D-6,   -1194.447010225D0, 4.136047594D0,
     :    0.000488D-6,    5849.364112115D0, 2.209679987D0,
     :    0.000621D-6,    6438.496249426D0, 4.518860804D0,
     :    0.000495D-6,   -6286.598968340D0, 1.868201275D0,
     :    0.000456D-6,    5230.807466803D0, 1.271231591D0,
     :    0.000451D-6,    5088.628839767D0, 0.084060889D0,
     :    0.000435D-6,    5643.178563677D0, 3.324456609D0,
     :    0.000387D-6,   10977.078804699D0, 4.052488477D0,
     :    0.000547D-6,  161000.685737473D0, 2.841633844D0,
     :    0.000522D-6,    3154.687084896D0, 2.171979966D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=731,740) /
     :    0.000375D-6,    5486.777843175D0, 4.983027306D0,
     :    0.000421D-6,    5863.591206116D0, 4.546432249D0,
     :    0.000439D-6,    7084.896781115D0, 0.522967921D0,
     :    0.000309D-6,    2544.314419883D0, 3.172606705D0,
     :    0.000347D-6,    4690.479836359D0, 1.479586566D0,
     :    0.000317D-6,     801.820931124D0, 3.553088096D0,
     :    0.000262D-6,     419.484643875D0, 0.606635550D0,
     :    0.000248D-6,    6836.645252834D0, 3.014082064D0,
     :    0.000245D-6,   -1592.596013633D0, 5.519526220D0,
     :    0.000225D-6,    4292.330832950D0, 2.877956536D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=741,750) /
     :    0.000214D-6,    7234.794256242D0, 1.605227587D0,
     :    0.000205D-6,    5767.611978898D0, 0.625804796D0,
     :    0.000180D-6,   10447.387839604D0, 3.499954526D0,
     :    0.000229D-6,     199.072001436D0, 5.632304604D0,
     :    0.000214D-6,     639.897286314D0, 5.960227667D0,
     :    0.000175D-6,    -433.711737877D0, 2.162417992D0,
     :    0.000209D-6,     515.463871093D0, 2.322150893D0,
     :    0.000173D-6,    6040.347246017D0, 2.556183691D0,
     :    0.000184D-6,    6309.374169791D0, 4.732296790D0,
     :    0.000227D-6,  149854.400134205D0, 5.385812217D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=751,760) /
     :    0.000154D-6,    8031.092263058D0, 5.120720920D0,
     :    0.000151D-6,    5739.157790895D0, 4.815000443D0,
     :    0.000197D-6,    7632.943259650D0, 0.222827271D0,
     :    0.000197D-6,      74.781598567D0, 3.910456770D0,
     :    0.000138D-6,    6055.549660552D0, 1.397484253D0,
     :    0.000149D-6,   -6127.655450557D0, 5.333727496D0,
     :    0.000137D-6,    3894.181829542D0, 4.281749907D0,
     :    0.000135D-6,    9437.762934887D0, 5.979971885D0,
     :    0.000139D-6,   -2352.866153772D0, 4.715630782D0,
     :    0.000142D-6,    6812.766815086D0, 0.513330157D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=761,770) /
     :    0.000120D-6,   -4705.732307544D0, 0.194160689D0,
     :    0.000131D-6,  -71430.695617928D0, 0.000379226D0,
     :    0.000124D-6,    6279.552731642D0, 2.122264908D0,
     :    0.000108D-6,   -6256.777530192D0, 0.883445696D0,

*  T^3
     :    0.143388D-6,    6283.075849991D0, 1.131453581D0,
     :    0.006671D-6,   12566.151699983D0, 0.775148887D0,
     :    0.001480D-6,     155.420399434D0, 0.480016880D0,
     :    0.000934D-6,     213.299095438D0, 6.144453084D0,
     :    0.000795D-6,     529.690965095D0, 2.941595619D0,
     :    0.000673D-6,    5746.271337896D0, 0.120415406D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=771,780) /
     :    0.000672D-6,    5760.498431898D0, 5.317009738D0,
     :    0.000389D-6,    -220.412642439D0, 3.090323467D0,
     :    0.000373D-6,    6062.663207553D0, 3.003551964D0,
     :    0.000360D-6,    6076.890301554D0, 1.918913041D0,
     :    0.000316D-6,     -21.340641002D0, 5.545798121D0,
     :    0.000315D-6,    -242.728603974D0, 1.884932563D0,
     :    0.000278D-6,     206.185548437D0, 1.266254859D0,
     :    0.000238D-6,    -536.804512095D0, 4.532664830D0,
     :    0.000185D-6,     522.577418094D0, 4.578313856D0,
     :    0.000245D-6,   18849.227549974D0, 0.587467082D0 /
      DATA ((FAIRHD(I,J),I=1,3),J=781,787) /
     :    0.000180D-6,     426.598190876D0, 5.151178553D0,
     :    0.000200D-6,     553.569402842D0, 5.355983739D0,
     :    0.000141D-6,    5223.693919802D0, 1.336556009D0,
     :    0.000104D-6,    5856.477659115D0, 4.239842759D0,

*  T^4
     :    0.003826D-6,    6283.075849991D0, 5.705257275D0,
     :    0.000303D-6,   12566.151699983D0, 5.407132842D0,
     :    0.000209D-6,     155.420399434D0, 1.989815753D0 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Time since J2000.0 in Julian millennia.
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJM

*  =================
*  Topocentric terms
*  =================

*  Convert UT to local solar time in radians.
      TSOL = MOD(UT,1D0) * D2PI + ELONG

*  FUNDAMENTAL ARGUMENTS:  Simon et al. 1994.

*  Combine time argument (millennia) with deg/arcsec factor.
      W = T / 3600D0

*  Sun Mean Longitude.
      ELSUN = MOD(280.46645683D0 + 1296027711.03429D0 * W, 360D0) * DD2R

*  Sun Mean Anomaly.
      EMSUN = MOD(357.52910918D0 + 1295965810.481D0 * W, 360D0) * DD2R

*  Mean Elongation of Moon from Sun.
      D = MOD(297.85019547D0 + 16029616012.090D0 * W, 360D0) * DD2R

*  Mean Longitude of Jupiter.
      ELJ = MOD(34.35151874D0 + 109306899.89453D0 * W, 360D0) * DD2R

*  Mean Longitude of Saturn.
      ELS = MOD(50.07744430D0 + 44046398.47038D0 * W, 360D0) * DD2R

*  TOPOCENTRIC TERMS:  Moyer 1981 and Murray 1983.
      WT =  + 0.00029D-10 * U * SIN(TSOL + ELSUN - ELS)
     :      + 0.00100D-10 * U * SIN(TSOL - 2D0*EMSUN)
     :      + 0.00133D-10 * U * SIN(TSOL - D)
     :      + 0.00133D-10 * U * SIN(TSOL + ELSUN - ELJ)
     :      - 0.00229D-10 * U * SIN(TSOL + 2D0*ELSUN + EMSUN)
     :      - 0.0220 D-10 * V * COS(ELSUN + EMSUN)
     :      + 0.05312D-10 * U * SIN(TSOL - EMSUN)
     :      - 0.13677D-10 * U * SIN(TSOL + 2D0*ELSUN)
     :      - 1.3184 D-10 * V * COS(ELSUN)
     :      + 3.17679D-10 * U * SIN(TSOL)

*  =====================
*  Fairhead et al. model
*  =====================

*  T**0
      W0 = 0D0
      DO 10 J=474,1,-1
         W0 = W0 + FAIRHD(1,J) * SIN(FAIRHD(2,J)*T + FAIRHD(3,J))
 10   CONTINUE

*  T**1
      W1 = 0D0
      DO 11 J=679,475,-1
         W1 = W1 + FAIRHD(1,J) * SIN(FAIRHD(2,J)*T + FAIRHD(3,J))
 11   CONTINUE

*  T**2
      W2 = 0D0
      DO 12 J=764,680,-1
         W2 = W2 + FAIRHD(1,J) * SIN(FAIRHD(2,J)*T + FAIRHD(3,J))
 12   CONTINUE

*  T**3
      W3 = 0D0
      DO 13 J=784,765,-1
         W3 = W3 + FAIRHD(1,J) * SIN(FAIRHD(2,J)*T + FAIRHD(3,J))
 13   CONTINUE

*  T**4
      W4 = 0D0
      DO 14 J=787,785,-1
         W4 = W4 + FAIRHD(1,J) * SIN(FAIRHD(2,J)*T + FAIRHD(3,J))
 14   CONTINUE

*  Multiply by powers of T and combine.
      WF = T * ( T * ( T * ( T * W4 + W3 ) + W2 ) + W1 ) + W0

*  Adjustments to use JPL planetary masses instead of IAU.
      WJ =   0.00065D-6 * SIN(6069.776754D0*T + 4.021194D0) +
     :       0.00033D-6 * SIN( 213.299095D0*T + 5.543132D0) +
     :    ( -0.00196D-6 * SIN(6208.294251D0*T + 5.696701D0) ) +
     :    ( -0.00173D-6 * SIN(  74.781599D0*T + 2.435900D0) ) +
     :       0.03638D-6 * T * T

*  ============
*  Final result
*  ============

*  TDB-TT in seconds.
      iau_DTDB = WT + WF + WJ

*  Finished.

*+-----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  Standards Of Fundamental Astronomy Review Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Review Board ("the Board").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and 
*     restrictions listed below.
*
*  3. You (the user) may copy and adapt the SOFA software and its 
*     algorithms for your own purposes and you may copy and distribute
*     a resulting "derived work" to others on a world-wide, royalty-free 
*     basis, provided that the derived work complies with the following
*     requirements: 
*
*     a) Your work shall be marked or carry a statement that it (i) uses
*        routines and computations derived by you from software provided 
*        by SOFA under license to you; and (ii) does not contain
*        software provided by SOFA or software that has been distributed
*        by or endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon and/or differs from the
*        original SOFA software.
*
*     c) The name(s) of all routine(s) that you distribute shall differ
*        from the SOFA names, even when the SOFA content has not been
*        otherwise changed.
*
*     d) The routine-naming prefix "iau" shall not be used.
*
*     e) The origin of the SOFA components of your derived work must not
*        be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     f) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have granted 
*        a further right to modify the source code of your derived work.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall acknowledge
*     that the SOFA software was used in obtaining those results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or by
*     inappropriate modification.
*
*  6. The SOFA software is provided "as is" and the Board makes no 
*     warranty as to its use or performance.   The Board does not and 
*     cannot warrant the performance or results which the user may obtain 
*     by using the SOFA software.  The Board makes no warranties, express 
*     or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Board be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Board representative has been advised of such
*     damages, or for any claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms 
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*     Internet email: sofa@rl.ac.uk
*     Postal address: IAU SOFA Center
*                     Rutherford Appleton Laboratory
*                     Chilton, Didcot, Oxon OX11 0QX
*                     United Kingdom
*
*-----------------------------------------------------------------------

      END
      SUBROUTINE iau_EPV00 ( DATE1, DATE2, PVH, PVB, JSTAT )
*+
*  - - - - - - - - - -
*   i a u _ E P V 0 0
*  - - - - - - - - - -
*
*  Earth position and velocity, heliocentric and barycentric, with
*  respect to the Barycentric Celestial Reference System.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     DATE1    d         TDB date part A (Note 1)
*     DATE2    d         TDB date part B (Note 1)
*
*  Returned:
*     PVH      d(3,2)    heliocentric Earth position/velocity (AU,AU/day)
*     PVB      d(3,2)    barycentric Earth position/velocity (AU,AU/day)
*     JSTAT    i         status: 0 = OK
*                               +1 = warning: date outside 1900-2100 AD
*
*  Notes:
*
*  1) The epoch EPOCH1+EPOCH2 is a Julian Date, apportioned in
*     any convenient way between the two arguments.  For example,
*     JD(TDB)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*              EPOCH1        EPOCH2
*
*           2450123.7D0        0D0        (JD method)
*            2451545D0      -1421.3D0     (J2000 method)
*           2400000.5D0     50123.2D0     (MJD method)
*           2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*     However, the accuracy of the result is more likely to be
*     limited by the algorithm itself than the way the epoch has been
*     expressed.
*
*  2) On return, the arrays PVH and PVB contain the following:
*
*        PVH(1,1)  x       }
*        PVH(2,1)  y       } heliocentric position, AU
*        PVH(3,1)  z       }
*
*        PVH(1,2)  xdot    }
*        PVH(2,2)  ydot    } heliocentric velocity, AU/d
*        PVH(3,2)  zdot    }
*
*        PVB(1,1)  x       }
*        PVB(2,1)  y       } barycentric position, AU
*        PVB(3,1)  z       }
*
*        PVB(1,2)  xdot    }
*        PVB(2,2)  ydot    } barycentric velocity, AU/d
*        PVB(3,2)  zdot    }
*
*     The vectors are with respect to the Barycentric Celestial
*     Reference System.  The time unit is one day in TDB.
*
*  3) The routine is a SIMPLIFIED SOLUTION from the planetary theory
*     VSOP2000 (X. Moisson, P. Bretagnon, 2001, Celes. Mechanics &
*     Dyn. Astron., 80, 3/4, 205-213) and is an adaptation of original
*     Fortran code supplied by P. Bretagnon (private comm., 2000).
*
*  4) Comparisons over the time span 1900-2100 with this simplified
*     solution and the JPL DE405 ephemeris give the following results:
*
*                                RMS    max
*           Heliocentric:
*              position error    3.7   11.2   km
*              velocity error    1.4    5.0   mm/s
*
*           Barycentric:
*              position error    4.6   13.4   km
*              velocity error    1.4    4.9   mm/s
*
*     Comparisons with the JPL DE406 ephemeris show that by 1800 and
*     2200 the position errors are approximately double their 1900-2100
*     size.  By 1500 and 2500 the deterioration is a factor of 10 and by
*     1000 and 3000 a factor of 60.  The velocity accuracy falls off at
*     about half that rate.
*
*  This revision:  2008 May 24
*
*  Copyright (C) 2008 IAU SOFA Review Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, PVH(3,2), PVB(3,2)
      INTEGER JSTAT

      DOUBLE PRECISION T, T2, XYZ, XYZD, A, B, C, CT, P, CP,
     :                 PH(3), VH(3), PB(3), VB(3), X, Y, Z

      INTEGER I, J, K

*  Days per Julian year
      DOUBLE PRECISION DJY
      PARAMETER ( DJY = 365.25D0 )

*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )
*
*  Matrix elements for orienting the analytical model to DE405/ICRF.
*
*  The corresponding Euler angles are:
*
*                        d  '  "
*    1st rotation    -  23 26 21.4091 about the x-axis  (obliquity)
*    2nd rotation    +         0.0475 about the z-axis  (RA offset)
*
*  These were obtained empirically, by comparisons with DE405 over
*  1900-2100.
*
      DOUBLE PRECISION AM12, AM13, AM21, AM22, AM23, AM32, AM33
      PARAMETER ( AM12 = +0.000000211284D0,
     :            AM13 = -0.000000091603D0,
     :            AM21 = -0.000000230286D0,
     :            AM22 = +0.917482137087D0,
     :            AM23 = -0.397776982902D0,
     :            AM32 = +0.397776982902D0,
     :            AM33 = +0.917482137087D0 )

*  ----------------------
*  Ephemeris Coefficients
*  ----------------------
*
*  The coefficients are stored in arrays of dimension (3,n,3).  There
*  are separate sets of arrays for (i) the Sun to Earth vector and
*  (ii) the Solar-System barycenter to Sun vector.  Each of these two
*  sets contains separate arrays for the terms (n in number) in each
*  power of time (in Julian years since J2000): T^0, T^1 and T^2.
*  Within each array, all the Cartesian x-components, elements (i,j,1),
*  appear first, followed by all the y-components, elements (i,j,2) and
*  finally all the z-components, elements (i,j,3).  At the lowest level
*  are groups of three coefficients.  The first coefficient in each
*  group, element (1,j,k), is the amplitude of the term, the second,
*  element (2,j,k), is the phase and the third, element (3,j,k), is the
*  frequency.
*
*  The naming scheme is such that a block
*
*     DOUBLE PRECISION bn(3,Mbn,3)
*
*  applies to body b and time exponent n:
*
*    . b can be either E (Earth with respect to Sun) or S (Sun with
*      respect to Solar-System Barycenter)
*
*    . n can be 0, 1 or 2, for T^0, T^1 or T^2
*
*  For example, array E2(3,ME2,3) contains the coefficients for
*  the T^2 terms for the Sun-to-Earth vector.
*
*  There is no requirement for the X, Y and Z models for a particular
*  block to use the same number of coefficients.  The number actually
*  used is parameterized, the number of terms being used called NbnX,
*  NbnY, and NbnZ respectively.  The parameter Mbn is the biggest of
*  the three, and defines the array size.  Unused elements are not
*  initialized and are never accessed.
*

      INTEGER NE0(3), NE0X, NE0Y, NE0Z, ME0,
     :        NE1(3), NE1X, NE1Y, NE1Z, ME1,
     :        NE2(3), NE2X, NE2Y, NE2Z, ME2,
     :        NS0(3), NS0X, NS0Y, NS0Z, MS0,
     :        NS1(3), NS1X, NS1Y, NS1Z, MS1,
     :        NS2(3), NS2X, NS2Y, NS2Z, MS2

      PARAMETER ( NE0X = 501, NE0Y = 501, NE0Z = 137, ME0 = NE0X,
     :            NE1X =  79, NE1Y =  80, NE1Z =  12, ME1 = NE1Y,
     :            NE2X =   5, NE2Y =   5, NE2Z =   3, ME2 = NE2X,
     :            NS0X = 212, NS0Y = 213, NS0Z =  69, MS0 = NS0Y,
     :            NS1X =  50, NS1Y =  50, NS1Z =  14, MS1 = NS1X,
     :            NS2X =   9, NS2Y =   9, NS2Z =   2, MS2 = NS2X )

      DOUBLE PRECISION E0(3,ME0,3), E1(3,ME1,3), E2(3,ME2,3),
     :                 S0(3,MS0,3), S1(3,MS1,3), S2(3,MS2,3)

      DATA NE0 / NE0X, NE0Y, NE0Z /
      DATA NE1 / NE1X, NE1Y, NE1Z /
      DATA NE2 / NE2X, NE2Y, NE2Z /
      DATA NS0 / NS0X, NS0Y, NS0Z /
      DATA NS1 / NS1X, NS1Y, NS1Z /
      DATA NS2 / NS2X, NS2Y, NS2Z /

*  Sun-to-Earth, T^0, X
      DATA ((E0(I,J,1),I=1,3),J=  1, 10) /
     :  0.9998292878132D+00, 0.1753485171504D+01, 0.6283075850446D+01,
     :  0.8352579567414D-02, 0.1710344404582D+01, 0.1256615170089D+02,
     :  0.5611445335148D-02, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.1046664295572D-03, 0.1667225416770D+01, 0.1884922755134D+02,
     :  0.3110842534677D-04, 0.6687513390251D+00, 0.8399684731857D+02,
     :  0.2552413503550D-04, 0.5830637358413D+00, 0.5296909721118D+00,
     :  0.2137207845781D-04, 0.1092330954011D+01, 0.1577343543434D+01,
     :  0.1680240182951D-04, 0.4955366134987D+00, 0.6279552690824D+01,
     :  0.1679012370795D-04, 0.6153014091901D+01, 0.6286599010068D+01,
     :  0.1445526946777D-04, 0.3472744100492D+01, 0.2352866153506D+01 /
      DATA ((E0(I,J,1),I=1,3),J= 11, 20) /
     :  0.1091038246184D-04, 0.3689845786119D+01, 0.5223693906222D+01,
     :  0.9344399733932D-05, 0.6073934645672D+01, 0.1203646072878D+02,
     :  0.8993182910652D-05, 0.3175705249069D+01, 0.1021328554739D+02,
     :  0.5665546034116D-05, 0.2152484672246D+01, 0.1059381944224D+01,
     :  0.6844146703035D-05, 0.1306964099750D+01, 0.5753384878334D+01,
     :  0.7346610905565D-05, 0.4354980070466D+01, 0.3981490189893D+00,
     :  0.6815396474414D-05, 0.2218229211267D+01, 0.4705732307012D+01,
     :  0.6112787253053D-05, 0.5384788425458D+01, 0.6812766822558D+01,
     :  0.4518120711239D-05, 0.6087604012291D+01, 0.5884926831456D+01,
     :  0.4521963430706D-05, 0.1279424524906D+01, 0.6256777527156D+01 /
      DATA ((E0(I,J,1),I=1,3),J= 21, 30) /
     :  0.4497426764085D-05, 0.5369129144266D+01, 0.6309374173736D+01,
     :  0.4062190566959D-05, 0.5436473303367D+00, 0.6681224869435D+01,
     :  0.5412193480192D-05, 0.7867838528395D+00, 0.7755226100720D+00,
     :  0.5469839049386D-05, 0.1461440311134D+01, 0.1414349524433D+02,
     :  0.5205264083477D-05, 0.4432944696116D+01, 0.7860419393880D+01,
     :  0.2149759935455D-05, 0.4502237496846D+01, 0.1150676975667D+02,
     :  0.2279109618501D-05, 0.1239441308815D+01, 0.7058598460518D+01,
     :  0.2259282939683D-05, 0.3272430985331D+01, 0.4694002934110D+01,
     :  0.2558950271319D-05, 0.2265471086404D+01, 0.1216800268190D+02,
     :  0.2561581447555D-05, 0.1454740653245D+01, 0.7099330490126D+00 /
      DATA ((E0(I,J,1),I=1,3),J= 31, 40) /
     :  0.1781441115440D-05, 0.2962068630206D+01, 0.7962980379786D+00,
     :  0.1612005874644D-05, 0.1473255041006D+01, 0.5486777812467D+01,
     :  0.1818630667105D-05, 0.3743903293447D+00, 0.6283008715021D+01,
     :  0.1818601377529D-05, 0.6274174354554D+01, 0.6283142985870D+01,
     :  0.1554475925257D-05, 0.1624110906816D+01, 0.2513230340178D+02,
     :  0.2090948029241D-05, 0.5852052276256D+01, 0.1179062909082D+02,
     :  0.2000176345460D-05, 0.4072093298513D+01, 0.1778984560711D+02,
     :  0.1289535917759D-05, 0.5217019331069D+01, 0.7079373888424D+01,
     :  0.1281135307881D-05, 0.4802054538934D+01, 0.3738761453707D+01,
     :  0.1518229005692D-05, 0.8691914742502D+00, 0.2132990797783D+00 /
      DATA ((E0(I,J,1),I=1,3),J= 41, 50) /
     :  0.9450128579027D-06, 0.4601859529950D+01, 0.1097707878456D+02,
     :  0.7781119494996D-06, 0.1844352816694D+01, 0.8827390247185D+01,
     :  0.7733407759912D-06, 0.3582790154750D+01, 0.5507553240374D+01,
     :  0.7350644318120D-06, 0.2695277788230D+01, 0.1589072916335D+01,
     :  0.6535928827023D-06, 0.3651327986142D+01, 0.1176985366291D+02,
     :  0.6324624183656D-06, 0.2241302375862D+01, 0.6262300422539D+01,
     :  0.6298565300557D-06, 0.4407122406081D+01, 0.6303851278352D+01,
     :  0.8587037089179D-06, 0.3024307223119D+01, 0.1672837615881D+03,
     :  0.8299954491035D-06, 0.6192539428237D+01, 0.3340612434717D+01,
     :  0.6311263503401D-06, 0.2014758795416D+01, 0.7113454667900D-02 /
      DATA ((E0(I,J,1),I=1,3),J= 51, 60) /
     :  0.6005646745452D-06, 0.3399500503397D+01, 0.4136910472696D+01,
     :  0.7917715109929D-06, 0.2493386877837D+01, 0.6069776770667D+01,
     :  0.7556958099685D-06, 0.4159491740143D+01, 0.6496374930224D+01,
     :  0.6773228244949D-06, 0.4034162934230D+01, 0.9437762937313D+01,
     :  0.5370708577847D-06, 0.1562219163734D+01, 0.1194447056968D+01,
     :  0.5710804266203D-06, 0.2662730803386D+01, 0.6282095334605D+01,
     :  0.5709824583726D-06, 0.3985828430833D+01, 0.6284056366286D+01,
     :  0.5143950896447D-06, 0.1308144688689D+01, 0.6290189305114D+01,
     :  0.5088010604546D-06, 0.5352817214804D+01, 0.6275962395778D+01,
     :  0.4960369085172D-06, 0.2644267922349D+01, 0.6127655567643D+01 /
      DATA ((E0(I,J,1),I=1,3),J= 61, 70) /
     :  0.4803137891183D-06, 0.4008844192080D+01, 0.6438496133249D+01,
     :  0.5731747768225D-06, 0.3794550174597D+01, 0.3154687086868D+01,
     :  0.4735947960579D-06, 0.6107118308982D+01, 0.3128388763578D+01,
     :  0.4808348796625D-06, 0.4771458618163D+01, 0.8018209333619D+00,
     :  0.4115073743137D-06, 0.3327111335159D+01, 0.8429241228195D+01,
     :  0.5230575889287D-06, 0.5305708551694D+01, 0.1336797263425D+02,
     :  0.5133977889215D-06, 0.5784230738814D+01, 0.1235285262111D+02,
     :  0.5065815825327D-06, 0.2052064793679D+01, 0.1185621865188D+02,
     :  0.4339831593868D-06, 0.3644994195830D+01, 0.1726015463500D+02,
     :  0.3952928638953D-06, 0.4930376436758D+01, 0.5481254917084D+01 /
      DATA ((E0(I,J,1),I=1,3),J= 71, 80) /
     :  0.4898498111942D-06, 0.4542084219731D+00, 0.9225539266174D+01,
     :  0.4757490209328D-06, 0.3161126388878D+01, 0.5856477690889D+01,
     :  0.4727701669749D-06, 0.6214993845446D+00, 0.2544314396739D+01,
     :  0.3800966681863D-06, 0.3040132339297D+01, 0.4265981595566D+00,
     :  0.3257301077939D-06, 0.8064977360087D+00, 0.3930209696940D+01,
     :  0.3255810528674D-06, 0.1974147981034D+01, 0.2146165377750D+01,
     :  0.3252029748187D-06, 0.2845924913135D+01, 0.4164311961999D+01,
     :  0.3255505635308D-06, 0.3017900824120D+01, 0.5088628793478D+01,
     :  0.2801345211990D-06, 0.6109717793179D+01, 0.1256967486051D+02,
     :  0.3688987740970D-06, 0.2911550235289D+01, 0.1807370494127D+02 /
      DATA ((E0(I,J,1),I=1,3),J= 81, 90) /
     :  0.2475153429458D-06, 0.2179146025856D+01, 0.2629832328990D-01,
     :  0.3033457749150D-06, 0.1994161050744D+01, 0.4535059491685D+01,
     :  0.2186743763110D-06, 0.5125687237936D+01, 0.1137170464392D+02,
     :  0.2764777032774D-06, 0.4822646860252D+00, 0.1256262854127D+02,
     :  0.2199028768592D-06, 0.4637633293831D+01, 0.1255903824622D+02,
     :  0.2046482824760D-06, 0.1467038733093D+01, 0.7084896783808D+01,
     :  0.2611209147507D-06, 0.3044718783485D+00, 0.7143069561767D+02,
     :  0.2286079656818D-06, 0.4764220356805D+01, 0.8031092209206D+01,
     :  0.1855071202587D-06, 0.3383637774428D+01, 0.1748016358760D+01,
     :  0.2324669506784D-06, 0.6189088449251D+01, 0.1831953657923D+02 /
      DATA ((E0(I,J,1),I=1,3),J= 91,100) /
     :  0.1709528015688D-06, 0.5874966729774D+00, 0.4933208510675D+01,
     :  0.2168156875828D-06, 0.4302994009132D+01, 0.1044738781244D+02,
     :  0.2106675556535D-06, 0.3800475419891D+01, 0.7477522907414D+01,
     :  0.1430213830465D-06, 0.1294660846502D+01, 0.2942463415728D+01,
     :  0.1388396901944D-06, 0.4594797202114D+01, 0.8635942003952D+01,
     :  0.1922258844190D-06, 0.4943044543591D+00, 0.1729818233119D+02,
     :  0.1888460058292D-06, 0.2426943912028D+01, 0.1561374759853D+03,
     :  0.1789449386107D-06, 0.1582973303499D+00, 0.1592596075957D+01,
     :  0.1360803685374D-06, 0.5197240440504D+01, 0.1309584267300D+02,
     :  0.1504038014709D-06, 0.3120360916217D+01, 0.1649636139783D+02 /
      DATA ((E0(I,J,1),I=1,3),J=101,110) /
     :  0.1382769533389D-06, 0.6164702888205D+01, 0.7632943190217D+01,
     :  0.1438059769079D-06, 0.1437423770979D+01, 0.2042657109477D+02,
     :  0.1326303260037D-06, 0.3609688799679D+01, 0.1213955354133D+02,
     :  0.1159244950540D-06, 0.5463018167225D+01, 0.5331357529664D+01,
     :  0.1433118149136D-06, 0.6028909912097D+01, 0.7342457794669D+01,
     :  0.1234623148594D-06, 0.3109645574997D+01, 0.6279485555400D+01,
     :  0.1233949875344D-06, 0.3539359332866D+01, 0.6286666145492D+01,
     :  0.9927196061299D-07, 0.1259321569772D+01, 0.7234794171227D+01,
     :  0.1242302191316D-06, 0.1065949392609D+01, 0.1511046609763D+02,
     :  0.1098402195201D-06, 0.2192508743837D+01, 0.1098880815746D+02 /
      DATA ((E0(I,J,1),I=1,3),J=111,120) /
     :  0.1158191395315D-06, 0.4054411278650D+01, 0.5729506548653D+01,
     :  0.9048475596241D-07, 0.5429764748518D+01, 0.9623688285163D+01,
     :  0.8889853269023D-07, 0.5046586206575D+01, 0.6148010737701D+01,
     :  0.1048694242164D-06, 0.2628858030806D+01, 0.6836645152238D+01,
     :  0.1112308378646D-06, 0.4177292719907D+01, 0.1572083878776D+02,
     :  0.8631729709901D-07, 0.1601345232557D+01, 0.6418140963190D+01,
     :  0.8527816951664D-07, 0.2463888997513D+01, 0.1471231707864D+02,
     :  0.7892139456991D-07, 0.3154022088718D+01, 0.2118763888447D+01,
     :  0.1051782905236D-06, 0.4795035816088D+01, 0.1349867339771D+01,
     :  0.1048219943164D-06, 0.2952983395230D+01, 0.5999216516294D+01 /
      DATA ((E0(I,J,1),I=1,3),J=121,130) /
     :  0.7435760775143D-07, 0.5420547991464D+01, 0.6040347114260D+01,
     :  0.9869574106949D-07, 0.3695646753667D+01, 0.6566935184597D+01,
     :  0.9156886364226D-07, 0.3922675306609D+01, 0.5643178611111D+01,
     :  0.7006834356188D-07, 0.1233968624861D+01, 0.6525804586632D+01,
     :  0.9806170182601D-07, 0.1919542280684D+01, 0.2122839202813D+02,
     :  0.9052289673607D-07, 0.4615902724369D+01, 0.4690479774488D+01,
     :  0.7554200867893D-07, 0.1236863719072D+01, 0.1253985337760D+02,
     :  0.8215741286498D-07, 0.3286800101559D+00, 0.1097355562493D+02,
     :  0.7185178575397D-07, 0.5880942158367D+01, 0.6245048154254D+01,
     :  0.7130726476180D-07, 0.7674871987661D+00, 0.6321103546637D+01 /
      DATA ((E0(I,J,1),I=1,3),J=131,140) /
     :  0.6650894461162D-07, 0.6987129150116D+00, 0.5327476111629D+01,
     :  0.7396888823688D-07, 0.3576824794443D+01, 0.5368044267797D+00,
     :  0.7420588884775D-07, 0.5033615245369D+01, 0.2354323048545D+02,
     :  0.6141181642908D-07, 0.9449927045673D+00, 0.1296430071988D+02,
     :  0.6373557924058D-07, 0.6206342280341D+01, 0.9517183207817D+00,
     :  0.6359474329261D-07, 0.5036079095757D+01, 0.1990745094947D+01,
     :  0.5740173582646D-07, 0.6105106371350D+01, 0.9555997388169D+00,
     :  0.7019864084602D-07, 0.7237747359018D+00, 0.5225775174439D+00,
     :  0.6398054487042D-07, 0.3976367969666D+01, 0.2407292145756D+02,
     :  0.7797092650498D-07, 0.4305423910623D+01, 0.2200391463820D+02 /
      DATA ((E0(I,J,1),I=1,3),J=141,150) /
     :  0.6466760000900D-07, 0.3500136825200D+01, 0.5230807360890D+01,
     :  0.7529417043890D-07, 0.3514779246100D+01, 0.1842262939178D+02,
     :  0.6924571140892D-07, 0.2743457928679D+01, 0.1554202828031D+00,
     :  0.6220798650222D-07, 0.2242598118209D+01, 0.1845107853235D+02,
     :  0.5870209391853D-07, 0.2332832707527D+01, 0.6398972393349D+00,
     :  0.6263953473888D-07, 0.2191105358956D+01, 0.6277552955062D+01,
     :  0.6257781390012D-07, 0.4457559396698D+01, 0.6288598745829D+01,
     :  0.5697304945123D-07, 0.3499234761404D+01, 0.1551045220144D+01,
     :  0.6335438746791D-07, 0.6441691079251D+00, 0.5216580451554D+01,
     :  0.6377258441152D-07, 0.2252599151092D+01, 0.5650292065779D+01 /
      DATA ((E0(I,J,1),I=1,3),J=151,160) /
     :  0.6484841818165D-07, 0.1992812417646D+01, 0.1030928125552D+00,
     :  0.4735551485250D-07, 0.3744672082942D+01, 0.1431416805965D+02,
     :  0.4628595996170D-07, 0.1334226211745D+01, 0.5535693017924D+00,
     :  0.6258152336933D-07, 0.4395836159154D+01, 0.2608790314060D+02,
     :  0.6196171366594D-07, 0.2587043007997D+01, 0.8467247584405D+02,
     :  0.6159556952126D-07, 0.4782499769128D+01, 0.2394243902548D+03,
     :  0.4987741172394D-07, 0.7312257619924D+00, 0.7771377146812D+02,
     :  0.5459280703142D-07, 0.3001376372532D+01, 0.6179983037890D+01,
     :  0.4863461189999D-07, 0.3767222128541D+01, 0.9027992316901D+02,
     :  0.5349912093158D-07, 0.3663594450273D+01, 0.6386168663001D+01 /
      DATA ((E0(I,J,1),I=1,3),J=161,170) /
     :  0.5673725607806D-07, 0.4331187919049D+01, 0.6915859635113D+01,
     :  0.4745485060512D-07, 0.5816195745518D+01, 0.6282970628506D+01,
     :  0.4745379005326D-07, 0.8323672435672D+00, 0.6283181072386D+01,
     :  0.4049002796321D-07, 0.3785023976293D+01, 0.6254626709878D+01,
     :  0.4247084014515D-07, 0.2378220728783D+01, 0.7875671926403D+01,
     :  0.4026912363055D-07, 0.2864103423269D+01, 0.6311524991013D+01,
     :  0.4062935011774D-07, 0.2415408595975D+01, 0.3634620989887D+01,
     :  0.5347771048509D-07, 0.3343479309801D+01, 0.2515860172507D+02,
     :  0.4829494136505D-07, 0.2821742398262D+01, 0.5760498333002D+01,
     :  0.4342554404599D-07, 0.5624662458712D+01, 0.7238675589263D+01 /
      DATA ((E0(I,J,1),I=1,3),J=171,180) /
     :  0.4021599184361D-07, 0.5557250275009D+00, 0.1101510648075D+02,
     :  0.4104900474558D-07, 0.3296691780005D+01, 0.6709674010002D+01,
     :  0.4376532905131D-07, 0.3814443999443D+01, 0.6805653367890D+01,
     :  0.3314590480650D-07, 0.3560229189250D+01, 0.1259245002418D+02,
     :  0.3232421839643D-07, 0.5185389180568D+01, 0.1066495398892D+01,
     :  0.3541176318876D-07, 0.3921381909679D+01, 0.9917696840332D+01,
     :  0.3689831242681D-07, 0.4190658955386D+01, 0.1192625446156D+02,
     :  0.3890605376774D-07, 0.5546023371097D+01, 0.7478166569050D-01,
     :  0.3038559339780D-07, 0.6231032794494D+01, 0.1256621883632D+02,
     :  0.3137083969782D-07, 0.6207063419190D+01, 0.4292330755499D+01 /
      DATA ((E0(I,J,1),I=1,3),J=181,190) /
     :  0.4024004081854D-07, 0.1195257375713D+01, 0.1334167431096D+02,
     :  0.3300234879283D-07, 0.1804694240998D+01, 0.1057540660594D+02,
     :  0.3635399155575D-07, 0.5597811343500D+01, 0.6208294184755D+01,
     :  0.3032668691356D-07, 0.3191059366530D+01, 0.1805292951336D+02,
     :  0.2809652069058D-07, 0.4094348032570D+01, 0.3523159621801D-02,
     :  0.3696955383823D-07, 0.5219282738794D+01, 0.5966683958112D+01,
     :  0.3562894142503D-07, 0.1037247544554D+01, 0.6357857516136D+01,
     :  0.3510598524148D-07, 0.1430020816116D+01, 0.6599467742779D+01,
     :  0.3617736142953D-07, 0.3002911403677D+01, 0.6019991944201D+01,
     :  0.2624524910730D-07, 0.2437046757292D+01, 0.6702560555334D+01 /
      DATA ((E0(I,J,1),I=1,3),J=191,200) /
     :  0.2535824204490D-07, 0.1581594689647D+01, 0.3141537925223D+02,
     :  0.3519787226257D-07, 0.5379863121521D+01, 0.2505706758577D+03,
     :  0.2578406709982D-07, 0.4904222639329D+01, 0.1673046366289D+02,
     :  0.3423887981473D-07, 0.3646448997315D+01, 0.6546159756691D+01,
     :  0.2776083886467D-07, 0.3307829300144D+01, 0.1272157198369D+02,
     :  0.3379592818379D-07, 0.1747541251125D+01, 0.1494531617769D+02,
     :  0.3050255426284D-07, 0.1784689432607D-01, 0.4732030630302D+01,
     :  0.2652378350236D-07, 0.4420055276260D+01, 0.5863591145557D+01,
     :  0.2374498173768D-07, 0.3629773929208D+01, 0.2388894113936D+01,
     :  0.2716451255140D-07, 0.3079623706780D+01, 0.1202934727411D+02 /
      DATA ((E0(I,J,1),I=1,3),J=201,210) /
     :  0.3038583699229D-07, 0.3312487903507D+00, 0.1256608456547D+02,
     :  0.2220681228760D-07, 0.5265520401774D+01, 0.1336244973887D+02,
     :  0.3044156540912D-07, 0.4766664081250D+01, 0.2908881142201D+02,
     :  0.2731859923561D-07, 0.5069146530691D+01, 0.1391601904066D+02,
     :  0.2285603018171D-07, 0.5954935112271D+01, 0.6076890225335D+01,
     :  0.2025006454555D-07, 0.4061789589267D+01, 0.4701116388778D+01,
     :  0.2012597519804D-07, 0.2485047705241D+01, 0.6262720680387D+01,
     :  0.2003406962258D-07, 0.4163779209320D+01, 0.6303431020504D+01,
     :  0.2207863441371D-07, 0.6923839133828D+00, 0.6489261475556D+01,
     :  0.2481374305624D-07, 0.5944173595676D+01, 0.1204357418345D+02 /
      DATA ((E0(I,J,1),I=1,3),J=211,220) /
     :  0.2130923288870D-07, 0.4641013671967D+01, 0.5746271423666D+01,
     :  0.2446370543391D-07, 0.6125796518757D+01, 0.1495633313810D+00,
     :  0.1932492759052D-07, 0.2234572324504D+00, 0.1352175143971D+02,
     :  0.2600122568049D-07, 0.4281012405440D+01, 0.4590910121555D+01,
     :  0.2431754047488D-07, 0.1429943874870D+00, 0.1162474756779D+01,
     :  0.1875902869209D-07, 0.9781803816948D+00, 0.6279194432410D+01,
     :  0.1874381139426D-07, 0.5670368130173D+01, 0.6286957268481D+01,
     :  0.2156696047173D-07, 0.2008985006833D+01, 0.1813929450232D+02,
     :  0.1965076182484D-07, 0.2566186202453D+00, 0.4686889479442D+01,
     :  0.2334816372359D-07, 0.4408121891493D+01, 0.1002183730415D+02 /
      DATA ((E0(I,J,1),I=1,3),J=221,230) /
     :  0.1869937408802D-07, 0.5272745038656D+01, 0.2427287361862D+00,
     :  0.2436236460883D-07, 0.4407720479029D+01, 0.9514313292143D+02,
     :  0.1761365216611D-07, 0.1943892315074D+00, 0.1351787002167D+02,
     :  0.2156289480503D-07, 0.1418570924545D+01, 0.6037244212485D+01,
     :  0.2164748979255D-07, 0.4724603439430D+01, 0.2301353951334D+02,
     :  0.2222286670853D-07, 0.2400266874598D+01, 0.1266924451345D+02,
     :  0.2070901414929D-07, 0.5230348028732D+01, 0.6528907488406D+01,
     :  0.1792745177020D-07, 0.2099190328945D+01, 0.6819880277225D+01,
     :  0.1841802068445D-07, 0.3467527844848D+00, 0.6514761976723D+02,
     :  0.1578401631718D-07, 0.7098642356340D+00, 0.2077542790660D-01 /
      DATA ((E0(I,J,1),I=1,3),J=231,240) /
     :  0.1561690152531D-07, 0.5943349620372D+01, 0.6272439236156D+01,
     :  0.1558591045463D-07, 0.7040653478980D+00, 0.6293712464735D+01,
     :  0.1737356469576D-07, 0.4487064760345D+01, 0.1765478049437D+02,
     :  0.1434755619991D-07, 0.2993391570995D+01, 0.1102062672231D+00,
     :  0.1482187806654D-07, 0.2278049198251D+01, 0.1052268489556D+01,
     :  0.1424812827089D-07, 0.1682114725827D+01, 0.1311972100268D+02,
     :  0.1380282448623D-07, 0.3262668602579D+01, 0.1017725758696D+02,
     :  0.1811481244566D-07, 0.3187771221777D+01, 0.1887552587463D+02,
     :  0.1504446185696D-07, 0.5650162308647D+01, 0.7626583626240D-01,
     :  0.1740776154137D-07, 0.5487068607507D+01, 0.1965104848470D+02 /
      DATA ((E0(I,J,1),I=1,3),J=241,250) /
     :  0.1374339536251D-07, 0.5745688172201D+01, 0.6016468784579D+01,
     :  0.1761377477704D-07, 0.5748060203659D+01, 0.2593412433514D+02,
     :  0.1535138225795D-07, 0.6226848505790D+01, 0.9411464614024D+01,
     :  0.1788140543676D-07, 0.6189318878563D+01, 0.3301902111895D+02,
     :  0.1375002807996D-07, 0.5371812884394D+01, 0.6327837846670D+00,
     :  0.1242115758632D-07, 0.1471687569712D+01, 0.3894181736510D+01,
     :  0.1450977333938D-07, 0.4143836662127D+01, 0.1277945078067D+02,
     :  0.1297579575023D-07, 0.9003477661957D+00, 0.6549682916313D+01,
     :  0.1462667934821D-07, 0.5760505536428D+01, 0.1863592847156D+02,
     :  0.1381774374799D-07, 0.1085471729463D+01, 0.2379164476796D+01 /
      DATA ((E0(I,J,1),I=1,3),J=251,260) /
     :  0.1682333169307D-07, 0.5409870870133D+01, 0.1620077269078D+02,
     :  0.1190812918837D-07, 0.1397205174601D+01, 0.1149965630200D+02,
     :  0.1221434762106D-07, 0.9001804809095D+00, 0.1257326515556D+02,
     :  0.1549934644860D-07, 0.4262528275544D+01, 0.1820933031200D+02,
     :  0.1252138953050D-07, 0.1411642012027D+01, 0.6993008899458D+01,
     :  0.1237078905387D-07, 0.2844472403615D+01, 0.2435678079171D+02,
     :  0.1446953389615D-07, 0.5295835522223D+01, 0.3813291813120D-01,
     :  0.1388446457170D-07, 0.4969428135497D+01, 0.2458316379602D+00,
     :  0.1019339179228D-07, 0.2491369561806D+01, 0.6112403035119D+01,
     :  0.1258880815343D-07, 0.4679426248976D+01, 0.5429879531333D+01 /
      DATA ((E0(I,J,1),I=1,3),J=261,270) /
     :  0.1297768238261D-07, 0.1074509953328D+01, 0.1249137003520D+02,
     :  0.9913505718094D-08, 0.4735097918224D+01, 0.6247047890016D+01,
     :  0.9830453155969D-08, 0.4158649187338D+01, 0.6453748665772D+01,
     :  0.1192615865309D-07, 0.3438208613699D+01, 0.6290122169689D+01,
     :  0.9835874798277D-08, 0.1913300781229D+01, 0.6319103810876D+01,
     :  0.9639087569277D-08, 0.9487683644125D+00, 0.8273820945392D+01,
     :  0.1175716107001D-07, 0.3228141664287D+01, 0.6276029531202D+01,
     :  0.1018926508678D-07, 0.2216607854300D+01, 0.1254537627298D+02,
     :  0.9500087869225D-08, 0.2625116459733D+01, 0.1256517118505D+02,
     :  0.9664192916575D-08, 0.5860562449214D+01, 0.6259197520765D+01 /
      DATA ((E0(I,J,1),I=1,3),J=271,280) /
     :  0.9612858712203D-08, 0.7885682917381D+00, 0.6306954180126D+01,
     :  0.1117645675413D-07, 0.3932148831189D+01, 0.1779695906178D+02,
     :  0.1158864052160D-07, 0.9995605521691D+00, 0.1778273215245D+02,
     :  0.9021043467028D-08, 0.5263769742673D+01, 0.6172869583223D+01,
     :  0.8836134773563D-08, 0.1496843220365D+01, 0.1692165728891D+01,
     :  0.1045872200691D-07, 0.7009039517214D+00, 0.2204125344462D+00,
     :  0.1211463487798D-07, 0.4041544938511D+01, 0.8257698122054D+02,
     :  0.8541990804094D-08, 0.1447586692316D+01, 0.6393282117669D+01,
     :  0.1038720703636D-07, 0.4594249718112D+00, 0.1550861511662D+02,
     :  0.1126722351445D-07, 0.3925550579036D+01, 0.2061856251104D+00 /
      DATA ((E0(I,J,1),I=1,3),J=281,290) /
     :  0.8697373859631D-08, 0.4411341856037D+01, 0.9491756770005D+00,
     :  0.8869380028441D-08, 0.2402659724813D+01, 0.3903911373650D+01,
     :  0.9247014693258D-08, 0.1401579743423D+01, 0.6267823317922D+01,
     :  0.9205062930950D-08, 0.5245978000814D+01, 0.6298328382969D+01,
     :  0.8000745038049D-08, 0.3590803356945D+01, 0.2648454860559D+01,
     :  0.9168973650819D-08, 0.2470150501679D+01, 0.1498544001348D+03,
     :  0.1075444949238D-07, 0.1328606161230D+01, 0.3694923081589D+02,
     :  0.7817298525817D-08, 0.6162256225998D+01, 0.4804209201333D+01,
     :  0.9541469226356D-08, 0.3942568967039D+01, 0.1256713221673D+02,
     :  0.9821910122027D-08, 0.2360246287233D+00, 0.1140367694411D+02 /
      DATA ((E0(I,J,1),I=1,3),J=291,300) /
     :  0.9897822023777D-08, 0.4619805634280D+01, 0.2280573557157D+02,
     :  0.7737289283765D-08, 0.3784727847451D+01, 0.7834121070590D+01,
     :  0.9260204034710D-08, 0.2223352487601D+01, 0.2787043132925D+01,
     :  0.7320252888486D-08, 0.1288694636874D+01, 0.6282655592598D+01,
     :  0.7319785780946D-08, 0.5359869567774D+01, 0.6283496108294D+01,
     :  0.7147219933778D-08, 0.5516616675856D+01, 0.1725663147538D+02,
     :  0.7946502829878D-08, 0.2630459984567D+01, 0.1241073141809D+02,
     :  0.9001711808932D-08, 0.2849815827227D+01, 0.6281591679874D+01,
     :  0.8994041507257D-08, 0.3795244450750D+01, 0.6284560021018D+01,
     :  0.8298582787358D-08, 0.5236413127363D+00, 0.1241658836951D+02 /
      DATA ((E0(I,J,1),I=1,3),J=301,310) /
     :  0.8526596520710D-08, 0.4794605424426D+01, 0.1098419223922D+02,
     :  0.8209822103197D-08, 0.1578752370328D+01, 0.1096996532989D+02,
     :  0.6357049861094D-08, 0.5708926113761D+01, 0.1596186371003D+01,
     :  0.7370473179049D-08, 0.3842402530241D+01, 0.4061219149443D+01,
     :  0.7232154664726D-08, 0.3067548981535D+01, 0.1610006857377D+03,
     :  0.6328765494903D-08, 0.1313930030069D+01, 0.1193336791622D+02,
     :  0.8030064908595D-08, 0.3488500408886D+01, 0.8460828644453D+00,
     :  0.6275464259232D-08, 0.1532061626198D+01, 0.8531963191132D+00,
     :  0.7051897446325D-08, 0.3285859929993D+01, 0.5849364236221D+01,
     :  0.6161593705428D-08, 0.1477341999464D+01, 0.5573142801433D+01 /
      DATA ((E0(I,J,1),I=1,3),J=311,320) /
     :  0.7754683957278D-08, 0.1586118663096D+01, 0.8662240327241D+01,
     :  0.5889928990701D-08, 0.1304887868803D+01, 0.1232342296471D+02,
     :  0.5705756047075D-08, 0.4555333589350D+01, 0.1258692712880D+02,
     :  0.5964178808332D-08, 0.3001762842062D+01, 0.5333900173445D+01,
     :  0.6712446027467D-08, 0.4886780007595D+01, 0.1171295538178D+02,
     :  0.5941809275464D-08, 0.4701509603824D+01, 0.9779108567966D+01,
     :  0.5466993627395D-08, 0.4588357817278D+01, 0.1884211409667D+02,
     :  0.6340512090980D-08, 0.1164543038893D+01, 0.5217580628120D+02,
     :  0.6325505710045D-08, 0.3919171259645D+01, 0.1041998632314D+02,
     :  0.6164789509685D-08, 0.2143828253542D+01, 0.6151533897323D+01 /
      DATA ((E0(I,J,1),I=1,3),J=321,330) /
     :  0.5263330812430D-08, 0.6066564434241D+01, 0.1885275071096D+02,
     :  0.5597087780221D-08, 0.2926316429472D+01, 0.4337116142245D+00,
     :  0.5396556236817D-08, 0.3244303591505D+01, 0.6286362197481D+01,
     :  0.5396615148223D-08, 0.3404304703662D+01, 0.6279789503410D+01,
     :  0.7091832443341D-08, 0.8532377803192D+00, 0.4907302013889D+01,
     :  0.6572352589782D-08, 0.4901966774419D+01, 0.1176433076753D+02,
     :  0.5960236060795D-08, 0.1874672315797D+01, 0.1422690933580D-01,
     :  0.5125480043511D-08, 0.3735726064334D+01, 0.1245594543367D+02,
     :  0.5928241866410D-08, 0.4502033899935D+01, 0.6414617803568D+01,
     :  0.5249600357424D-08, 0.4372334799878D+01, 0.1151388321134D+02 /
      DATA ((E0(I,J,1),I=1,3),J=331,340) /
     :  0.6059171276087D-08, 0.2581617302908D+01, 0.6062663316000D+01,
     :  0.5295235081662D-08, 0.2974811513158D+01, 0.3496032717521D+01,
     :  0.5820561875933D-08, 0.1796073748244D+00, 0.2838593341516D+00,
     :  0.4754696606440D-08, 0.1981998136973D+01, 0.3104930017775D+01,
     :  0.6385053548955D-08, 0.2559174171605D+00, 0.6133512519065D+01,
     :  0.6589828273941D-08, 0.2750967106776D+01, 0.4087944051283D+02,
     :  0.5383376567189D-08, 0.6325947523578D+00, 0.2248384854122D+02,
     :  0.5928941683538D-08, 0.1672304519067D+01, 0.1581959461667D+01,
     :  0.4816060709794D-08, 0.3512566172575D+01, 0.9388005868221D+01,
     :  0.6003381586512D-08, 0.5610932219189D+01, 0.5326786718777D+01 /
      DATA ((E0(I,J,1),I=1,3),J=341,350) /
     :  0.5504225393105D-08, 0.4037501131256D+01, 0.6503488384892D+01,
     :  0.5353772620129D-08, 0.6122774968240D+01, 0.1735668374386D+03,
     :  0.5786253768544D-08, 0.5527984999515D+01, 0.1350651127443D+00,
     :  0.5065706702002D-08, 0.9980765573624D+00, 0.1248988586463D+02,
     :  0.5972838885276D-08, 0.6044489493203D+01, 0.2673594526851D+02,
     :  0.5323585877961D-08, 0.3924265998147D+01, 0.4171425416666D+01,
     :  0.5210772682858D-08, 0.6220111376901D+01, 0.2460261242967D+02,
     :  0.4726549040535D-08, 0.3716043206862D+01, 0.7232251527446D+01,
     :  0.6029425105059D-08, 0.8548704071116D+00, 0.3227113045244D+03,
     :  0.4481542826513D-08, 0.1426925072829D+01, 0.5547199253223D+01 /
      DATA ((E0(I,J,1),I=1,3),J=351,360) /
     :  0.5836024505068D-08, 0.7135651752625D-01, 0.7285056171570D+02,
     :  0.4137046613272D-08, 0.5330767643283D+01, 0.1087398597200D+02,
     :  0.5171977473924D-08, 0.4494262335353D+00, 0.1884570439172D+02,
     :  0.5694429833732D-08, 0.2952369582215D+01, 0.9723862754494D+02,
     :  0.4009158925298D-08, 0.3500003416535D+01, 0.6244942932314D+01,
     :  0.4784939596873D-08, 0.6196709413181D+01, 0.2929661536378D+02,
     :  0.3983725022610D-08, 0.5103690031897D+01, 0.4274518229222D+01,
     :  0.3870535232462D-08, 0.3187569587401D+01, 0.6321208768577D+01,
     :  0.5140501213951D-08, 0.1668924357457D+01, 0.1232032006293D+02,
     :  0.3849034819355D-08, 0.4445722510309D+01, 0.1726726808967D+02 /
      DATA ((E0(I,J,1),I=1,3),J=361,370) /
     :  0.4002383075060D-08, 0.5226224152423D+01, 0.7018952447668D+01,
     :  0.3890719543549D-08, 0.4371166550274D+01, 0.1491901785440D+02,
     :  0.4887084607881D-08, 0.5973556689693D+01, 0.1478866649112D+01,
     :  0.3739939287592D-08, 0.2089084714600D+01, 0.6922973089781D+01,
     :  0.5031925918209D-08, 0.4658371936827D+01, 0.1715706182245D+02,
     :  0.4387748764954D-08, 0.4825580552819D+01, 0.2331413144044D+03,
     :  0.4147398098865D-08, 0.3739003524998D+01, 0.1376059875786D+02,
     :  0.3719089993586D-08, 0.1148941386536D+01, 0.6297302759782D+01,
     :  0.3934238461056D-08, 0.1559893008343D+01, 0.7872148766781D+01,
     :  0.3672471375622D-08, 0.5516145383612D+01, 0.6268848941110D+01 /
      DATA ((E0(I,J,1),I=1,3),J=371,380) /
     :  0.3768911277583D-08, 0.6116053700563D+01, 0.4157198507331D+01,
     :  0.4033388417295D-08, 0.5076821746017D+01, 0.1567108171867D+02,
     :  0.3764194617832D-08, 0.8164676232075D+00, 0.3185192151914D+01,
     :  0.4840628226284D-08, 0.1360479453671D+01, 0.1252801878276D+02,
     :  0.4949443923785D-08, 0.2725622229926D+01, 0.1617106187867D+03,
     :  0.4117393089971D-08, 0.6054459628492D+00, 0.5642198095270D+01,
     :  0.3925754020428D-08, 0.8570462135210D+00, 0.2139354194808D+02,
     :  0.3630551757923D-08, 0.3552067338279D+01, 0.6294805223347D+01,
     :  0.3627274802357D-08, 0.3096565085313D+01, 0.6271346477544D+01,
     :  0.3806143885093D-08, 0.6367751709777D+00, 0.1725304118033D+02 /
      DATA ((E0(I,J,1),I=1,3),J=381,390) /
     :  0.4433254641565D-08, 0.4848461503937D+01, 0.7445550607224D+01,
     :  0.3712319846576D-08, 0.1331950643655D+01, 0.4194847048887D+00,
     :  0.3849847534783D-08, 0.4958368297746D+00, 0.9562891316684D+00,
     :  0.3483955430165D-08, 0.2237215515707D+01, 0.1161697602389D+02,
     :  0.3961912730982D-08, 0.3332402188575D+01, 0.2277943724828D+02,
     :  0.3419978244481D-08, 0.5785600576016D+01, 0.1362553364512D+02,
     :  0.3329417758177D-08, 0.9812676559709D-01, 0.1685848245639D+02,
     :  0.4207206893193D-08, 0.9494780468236D+00, 0.2986433403208D+02,
     :  0.3268548976410D-08, 0.1739332095686D+00, 0.5749861718712D+01,
     :  0.3321880082685D-08, 0.1423354800666D+01, 0.6279143387820D+01 /
      DATA ((E0(I,J,1),I=1,3),J=391,400) /
     :  0.4503173010852D-08, 0.2314972675293D+00, 0.1385561574497D+01,
     :  0.4316599090954D-08, 0.1012646782616D+00, 0.4176041334900D+01,
     :  0.3283493323850D-08, 0.5233306881265D+01, 0.6287008313071D+01,
     :  0.3164033542343D-08, 0.4005597257511D+01, 0.2099539292909D+02,
     :  0.4159720956725D-08, 0.5365676242020D+01, 0.5905702259363D+01,
     :  0.3565176892217D-08, 0.4284440620612D+01, 0.3932462625300D-02,
     :  0.3514440950221D-08, 0.4270562636575D+01, 0.7335344340001D+01,
     :  0.3540596871909D-08, 0.5953553201060D+01, 0.1234573916645D+02,
     :  0.2960769905118D-08, 0.1115180417718D+01, 0.2670964694522D+02,
     :  0.2962213739684D-08, 0.3863811918186D+01, 0.6408777551755D+00 /
      DATA ((E0(I,J,1),I=1,3),J=401,410) /
     :  0.3883556700251D-08, 0.1268617928302D+01, 0.6660449441528D+01,
     :  0.2919225516346D-08, 0.4908605223265D+01, 0.1375773836557D+01,
     :  0.3115158863370D-08, 0.3744519976885D+01, 0.3802769619140D-01,
     :  0.4099438144212D-08, 0.4173244670532D+01, 0.4480965020977D+02,
     :  0.2899531858964D-08, 0.5910601428850D+01, 0.2059724391010D+02,
     :  0.3289733429855D-08, 0.2488050078239D+01, 0.1081813534213D+02,
     :  0.3933075612875D-08, 0.1122363652883D+01, 0.3773735910827D+00,
     :  0.3021403764467D-08, 0.4951973724904D+01, 0.2982630633589D+02,
     :  0.2798598949757D-08, 0.5117057845513D+01, 0.1937891852345D+02,
     :  0.3397421302707D-08, 0.6104159180476D+01, 0.6923953605621D+01 /
      DATA ((E0(I,J,1),I=1,3),J=411,420) /
     :  0.3720398002179D-08, 0.1184933429829D+01, 0.3066615496545D+02,
     :  0.3598484186267D-08, 0.3505282086105D+01, 0.6147450479709D+01,
     :  0.3694594027310D-08, 0.2286651088141D+01, 0.2636725487657D+01,
     :  0.2680444152969D-08, 0.1871816775482D+00, 0.6816289982179D+01,
     :  0.3497574865641D-08, 0.3143251755431D+01, 0.6418701221183D+01,
     :  0.3130274129494D-08, 0.2462167316018D+01, 0.1235996607578D+02,
     :  0.3241119069551D-08, 0.4256374004686D+01, 0.1652265972112D+02,
     :  0.2601960842061D-08, 0.4970362941425D+01, 0.1045450126711D+02,
     :  0.2690601527504D-08, 0.2372657824898D+01, 0.3163918923335D+00,
     :  0.2908688152664D-08, 0.4232652627721D+01, 0.2828699048865D+02 /
      DATA ((E0(I,J,1),I=1,3),J=421,430) /
     :  0.3120456131875D-08, 0.3925747001137D+00, 0.2195415756911D+02,
     :  0.3148855423384D-08, 0.3093478330445D+01, 0.1172006883645D+02,
     :  0.3051044261017D-08, 0.5560948248212D+01, 0.6055599646783D+01,
     :  0.2826006876660D-08, 0.5072790310072D+01, 0.5120601093667D+01,
     :  0.3100034191711D-08, 0.4998530231096D+01, 0.1799603123222D+02,
     :  0.2398771640101D-08, 0.2561739802176D+01, 0.6255674361143D+01,
     :  0.2384002842728D-08, 0.4087420284111D+01, 0.6310477339748D+01,
     :  0.2842146517568D-08, 0.2515048217955D+01, 0.5469525544182D+01,
     :  0.2847674371340D-08, 0.5235326497443D+01, 0.1034429499989D+02,
     :  0.2903722140764D-08, 0.1088200795797D+01, 0.6510552054109D+01 /
      DATA ((E0(I,J,1),I=1,3),J=431,440) /
     :  0.3187610710605D-08, 0.4710624424816D+01, 0.1693792562116D+03,
     :  0.3048869992813D-08, 0.2857975896445D+00, 0.8390110365991D+01,
     :  0.2860216950984D-08, 0.2241619020815D+01, 0.2243449970715D+00,
     :  0.2701117683113D-08, 0.6651573305272D-01, 0.6129297044991D+01,
     :  0.2509891590152D-08, 0.1285135324585D+01, 0.1044027435778D+02,
     :  0.2623200252223D-08, 0.2981229834530D+00, 0.6436854655901D+01,
     :  0.2622541669202D-08, 0.6122470726189D+01, 0.9380959548977D+01,
     :  0.2818435667099D-08, 0.4251087148947D+01, 0.5934151399930D+01,
     :  0.2365196797465D-08, 0.3465070460790D+01, 0.2470570524223D+02,
     :  0.2358704646143D-08, 0.5791603815350D+01, 0.8671969964381D+01 /
      DATA ((E0(I,J,1),I=1,3),J=441,450) /
     :  0.2388299481390D-08, 0.4142483772941D+01, 0.7096626156709D+01,
     :  0.1996041217224D-08, 0.2101901889496D+01, 0.1727188400790D+02,
     :  0.2687593060336D-08, 0.1526689456959D+01, 0.7075506709219D+02,
     :  0.2618913670810D-08, 0.2397684236095D+01, 0.6632000300961D+01,
     :  0.2571523050364D-08, 0.5751929456787D+00, 0.6206810014183D+01,
     :  0.2582135006946D-08, 0.5595464352926D+01, 0.4873985990671D+02,
     :  0.2372530190361D-08, 0.5092689490655D+01, 0.1590676413561D+02,
     :  0.2357178484712D-08, 0.4444363527851D+01, 0.3097883698531D+01,
     :  0.2451590394723D-08, 0.3108251687661D+01, 0.6612329252343D+00,
     :  0.2370045949608D-08, 0.2608133861079D+01, 0.3459636466239D+02 /
      DATA ((E0(I,J,1),I=1,3),J=451,460) /
     :  0.2268997267358D-08, 0.3639717753384D+01, 0.2844914056730D-01,
     :  0.1731432137906D-08, 0.1741898445707D+00, 0.2019909489111D+02,
     :  0.1629869741622D-08, 0.3902225646724D+01, 0.3035599730800D+02,
     :  0.2206215801974D-08, 0.4971131250731D+01, 0.6281667977667D+01,
     :  0.2205469554680D-08, 0.1677462357110D+01, 0.6284483723224D+01,
     :  0.2148792362509D-08, 0.4236259604006D+01, 0.1980482729015D+02,
     :  0.1873733657847D-08, 0.5926814998687D+01, 0.2876692439167D+02,
     :  0.2026573758959D-08, 0.4349643351962D+01, 0.2449240616245D+02,
     :  0.1807770325110D-08, 0.5700940482701D+01, 0.2045286941806D+02,
     :  0.1881174408581D-08, 0.6601286363430D+00, 0.2358125818164D+02 /
      DATA ((E0(I,J,1),I=1,3),J=461,470) /
     :  0.1368023671690D-08, 0.2211098592752D+01, 0.2473415438279D+02,
     :  0.1720017916280D-08, 0.4942488551129D+01, 0.1679593901136D+03,
     :  0.1702427665131D-08, 0.1452233856386D+01, 0.3338575901272D+03,
     :  0.1414032510054D-08, 0.5525357721439D+01, 0.1624205518357D+03,
     :  0.1652626045364D-08, 0.4108794283624D+01, 0.8956999012000D+02,
     :  0.1642957769686D-08, 0.7344335209984D+00, 0.5267006960365D+02,
     :  0.1614952403624D-08, 0.3541213951363D+01, 0.3332657872986D+02,
     :  0.1535988291188D-08, 0.4031094072151D+01, 0.3852657435933D+02,
     :  0.1593193738177D-08, 0.4185136203609D+01, 0.2282781046519D+03,
     :  0.1074569126382D-08, 0.1720485636868D+01, 0.8397383534231D+02 /
      DATA ((E0(I,J,1),I=1,3),J=471,480) /
     :  0.1074408214509D-08, 0.2758613420318D+01, 0.8401985929482D+02,
     :  0.9700199670465D-09, 0.4216686842097D+01, 0.7826370942180D+02,
     :  0.1258433517061D-08, 0.2575068876639D+00, 0.3115650189215D+03,
     :  0.1240303229539D-08, 0.4800844956756D+00, 0.1784300471910D+03,
     :  0.9018345948127D-09, 0.3896756361552D+00, 0.5886454391678D+02,
     :  0.1135301432805D-08, 0.3700805023550D+00, 0.7842370451713D+02,
     :  0.9215887951370D-09, 0.4364579276638D+01, 0.1014262087719D+03,
     :  0.1055401054147D-08, 0.2156564222111D+01, 0.5660027930059D+02,
     :  0.1008725979831D-08, 0.5454015785234D+01, 0.4245678405627D+02,
     :  0.7217398104321D-09, 0.1597772562175D+01, 0.2457074661053D+03 /
      DATA ((E0(I,J,1),I=1,3),J=481,490) /
     :  0.6912033134447D-09, 0.5824090621461D+01, 0.1679936946371D+03,
     :  0.6833881523549D-09, 0.3578778482835D+01, 0.6053048899753D+02,
     :  0.4887304205142D-09, 0.3724362812423D+01, 0.9656299901946D+02,
     :  0.5173709754788D-09, 0.5422427507933D+01, 0.2442876000072D+03,
     :  0.4671353097145D-09, 0.2396106924439D+01, 0.1435713242844D+03,
     :  0.5652608439480D-09, 0.2804028838685D+01, 0.8365903305582D+02,
     :  0.5604061331253D-09, 0.1638816006247D+01, 0.8433466158131D+02,
     :  0.4712723365400D-09, 0.8979003224474D+00, 0.3164282286739D+03,
     :  0.4909967465112D-09, 0.3210426725516D+01, 0.4059982187939D+03,
     :  0.4771358267658D-09, 0.5308027211629D+01, 0.1805255418145D+03 /
      DATA ((E0(I,J,1),I=1,3),J=491,500) /
     :  0.3943451445989D-09, 0.2195145341074D+01, 0.2568537517081D+03,
     :  0.3952109120244D-09, 0.5081189491586D+01, 0.2449975330562D+03,
     :  0.3788134594789D-09, 0.4345171264441D+01, 0.1568131045107D+03,
     :  0.3738330190479D-09, 0.2613062847997D+01, 0.3948519331910D+03,
     :  0.3099866678136D-09, 0.2846760817689D+01, 0.1547176098872D+03,
     :  0.2002962716768D-09, 0.4921360989412D+01, 0.2268582385539D+03,
     :  0.2198291338754D-09, 0.1130360117454D+00, 0.1658638954901D+03,
     :  0.1491958330784D-09, 0.4228195232278D+01, 0.2219950288015D+03,
     :  0.1475384076173D-09, 0.3005721811604D+00, 0.3052819430710D+03,
     :  0.1661626624624D-09, 0.7830125621203D+00, 0.2526661704812D+03 /
      DATA ((E0(I,J,1),I=1,3),J=501,NE0X) /
     :  0.9015823460025D-10, 0.3807792942715D+01, 0.4171445043968D+03 /

*  Sun-to-Earth, T^1, X
      DATA ((E1(I,J,1),I=1,3),J=  1, 10) /
     :  0.1234046326004D-05, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.5150068824701D-06, 0.6002664557501D+01, 0.1256615170089D+02,
     :  0.1290743923245D-07, 0.5959437664199D+01, 0.1884922755134D+02,
     :  0.1068615564952D-07, 0.2015529654209D+01, 0.6283075850446D+01,
     :  0.2079619142538D-08, 0.1732960531432D+01, 0.6279552690824D+01,
     :  0.2078009243969D-08, 0.4915604476996D+01, 0.6286599010068D+01,
     :  0.6206330058856D-09, 0.3616457953824D+00, 0.4705732307012D+01,
     :  0.5989335313746D-09, 0.3802607304474D+01, 0.6256777527156D+01,
     :  0.5958495663840D-09, 0.2845866560031D+01, 0.6309374173736D+01,
     :  0.4866923261539D-09, 0.5213203771824D+01, 0.7755226100720D+00 /
      DATA ((E1(I,J,1),I=1,3),J= 11, 20) /
     :  0.4267785823142D-09, 0.4368189727818D+00, 0.1059381944224D+01,
     :  0.4610675141648D-09, 0.1837249181372D-01, 0.7860419393880D+01,
     :  0.3626989993973D-09, 0.2161590545326D+01, 0.5753384878334D+01,
     :  0.3563071194389D-09, 0.1452631954746D+01, 0.5884926831456D+01,
     :  0.3557015642807D-09, 0.4470593393054D+01, 0.6812766822558D+01,
     :  0.3210412089122D-09, 0.5195926078314D+01, 0.6681224869435D+01,
     :  0.2875473577986D-09, 0.5916256610193D+01, 0.2513230340178D+02,
     :  0.2842913681629D-09, 0.1149902426047D+01, 0.6127655567643D+01,
     :  0.2751248215916D-09, 0.5502088574662D+01, 0.6438496133249D+01,
     :  0.2481432881127D-09, 0.2921989846637D+01, 0.5486777812467D+01 /
      DATA ((E1(I,J,1),I=1,3),J= 21, 30) /
     :  0.2059885976560D-09, 0.3718070376585D+01, 0.7079373888424D+01,
     :  0.2015522342591D-09, 0.5979395259740D+01, 0.6290189305114D+01,
     :  0.1995364084253D-09, 0.6772087985494D+00, 0.6275962395778D+01,
     :  0.1957436436943D-09, 0.2899210654665D+01, 0.5507553240374D+01,
     :  0.1651609818948D-09, 0.6228206482192D+01, 0.1150676975667D+02,
     :  0.1822980550699D-09, 0.1469348746179D+01, 0.1179062909082D+02,
     :  0.1675223159760D-09, 0.3813910555688D+01, 0.7058598460518D+01,
     :  0.1706491764745D-09, 0.3004380506684D+00, 0.7113454667900D-02,
     :  0.1392952362615D-09, 0.1440393973406D+01, 0.7962980379786D+00,
     :  0.1209868266342D-09, 0.4150425791727D+01, 0.4694002934110D+01 /
      DATA ((E1(I,J,1),I=1,3),J= 31, 40) /
     :  0.1009827202611D-09, 0.3290040429843D+01, 0.3738761453707D+01,
     :  0.1047261388602D-09, 0.4229590090227D+01, 0.6282095334605D+01,
     :  0.1047006652004D-09, 0.2418967680575D+01, 0.6284056366286D+01,
     :  0.9609993143095D-10, 0.4627943659201D+01, 0.6069776770667D+01,
     :  0.9590900593873D-10, 0.1894393939924D+01, 0.4136910472696D+01,
     :  0.9146249188071D-10, 0.2010647519562D+01, 0.6496374930224D+01,
     :  0.8545274480290D-10, 0.5529846956226D-01, 0.1194447056968D+01,
     :  0.8224377881194D-10, 0.1254304102174D+01, 0.1589072916335D+01,
     :  0.6183529510410D-10, 0.3360862168815D+01, 0.8827390247185D+01,
     :  0.6259255147141D-10, 0.4755628243179D+01, 0.8429241228195D+01 /
      DATA ((E1(I,J,1),I=1,3),J= 41, 50) /
     :  0.5539291694151D-10, 0.5371746955142D+01, 0.4933208510675D+01,
     :  0.7328259466314D-10, 0.4927699613906D+00, 0.4535059491685D+01,
     :  0.6017835843560D-10, 0.5776682001734D-01, 0.1255903824622D+02,
     :  0.7079827775243D-10, 0.4395059432251D+01, 0.5088628793478D+01,
     :  0.5170358878213D-10, 0.5154062619954D+01, 0.1176985366291D+02,
     :  0.4872301838682D-10, 0.6289611648973D+00, 0.6040347114260D+01,
     :  0.5249869411058D-10, 0.5617272046949D+01, 0.3154687086868D+01,
     :  0.4716172354411D-10, 0.3965901800877D+01, 0.5331357529664D+01,
     :  0.4871214940964D-10, 0.4627507050093D+01, 0.1256967486051D+02,
     :  0.4598076850751D-10, 0.6023631226459D+01, 0.6525804586632D+01 /
      DATA ((E1(I,J,1),I=1,3),J= 51, 60) /
     :  0.4562196089485D-10, 0.4138562084068D+01, 0.3930209696940D+01,
     :  0.4325493872224D-10, 0.1330845906564D+01, 0.7632943190217D+01,
     :  0.5673781176748D-10, 0.2558752615657D+01, 0.5729506548653D+01,
     :  0.3961436642503D-10, 0.2728071734630D+01, 0.7234794171227D+01,
     :  0.5101868209058D-10, 0.4113444965144D+01, 0.6836645152238D+01,
     :  0.5257043167676D-10, 0.6195089830590D+01, 0.8031092209206D+01,
     :  0.5076613989393D-10, 0.2305124132918D+01, 0.7477522907414D+01,
     :  0.3342169352778D-10, 0.5415998155071D+01, 0.1097707878456D+02,
     :  0.3545881983591D-10, 0.3727160564574D+01, 0.4164311961999D+01,
     :  0.3364063738599D-10, 0.2901121049204D+00, 0.1137170464392D+02 /
      DATA ((E1(I,J,1),I=1,3),J= 61, 70) /
     :  0.3357039670776D-10, 0.1652229354331D+01, 0.5223693906222D+01,
     :  0.4307412268687D-10, 0.4938909587445D+01, 0.1592596075957D+01,
     :  0.3405769115435D-10, 0.2408890766511D+01, 0.3128388763578D+01,
     :  0.3001926198480D-10, 0.4862239006386D+01, 0.1748016358760D+01,
     :  0.2778264787325D-10, 0.5241168661353D+01, 0.7342457794669D+01,
     :  0.2676159480666D-10, 0.3423593942199D+01, 0.2146165377750D+01,
     :  0.2954273399939D-10, 0.1881721265406D+01, 0.5368044267797D+00,
     :  0.3309362888795D-10, 0.1931525677349D+01, 0.8018209333619D+00,
     :  0.2810283608438D-10, 0.2414659495050D+01, 0.5225775174439D+00,
     :  0.3378045637764D-10, 0.4238019163430D+01, 0.1554202828031D+00 /
      DATA ((E1(I,J,1),I=1,3),J= 71,NE1X) /
     :  0.2558134979840D-10, 0.1828225235805D+01, 0.5230807360890D+01,
     :  0.2273755578447D-10, 0.5858184283998D+01, 0.7084896783808D+01,
     :  0.2294176037690D-10, 0.4514589779057D+01, 0.1726015463500D+02,
     :  0.2533506099435D-10, 0.2355717851551D+01, 0.5216580451554D+01,
     :  0.2716685375812D-10, 0.2221003625100D+01, 0.8635942003952D+01,
     :  0.2419043435198D-10, 0.5955704951635D+01, 0.4690479774488D+01,
     :  0.2521232544812D-10, 0.1395676848521D+01, 0.5481254917084D+01,
     :  0.2630195021491D-10, 0.5727468918743D+01, 0.2629832328990D-01,
     :  0.2548395840944D-10, 0.2628351859400D-03, 0.1349867339771D+01 /

*  Sun-to-Earth, T^2, X
      DATA ((E2(I,J,1),I=1,3),J=  1,NE2X) /
     : -0.4143818297913D-10, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.2171497694435D-10, 0.4398225628264D+01, 0.1256615170089D+02,
     :  0.9845398442516D-11, 0.2079720838384D+00, 0.6283075850446D+01,
     :  0.9256833552682D-12, 0.4191264694361D+01, 0.1884922755134D+02,
     :  0.1022049384115D-12, 0.5381133195658D+01, 0.8399684731857D+02 /

*  Sun-to-Earth, T^0, Y
      DATA ((E0(I,J,2),I=1,3),J=  1, 10) /
     :  0.9998921098898D+00, 0.1826583913846D+00, 0.6283075850446D+01,
     : -0.2442700893735D-01, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.8352929742915D-02, 0.1395277998680D+00, 0.1256615170089D+02,
     :  0.1046697300177D-03, 0.9641423109763D-01, 0.1884922755134D+02,
     :  0.3110841876663D-04, 0.5381140401712D+01, 0.8399684731857D+02,
     :  0.2570269094593D-04, 0.5301016407128D+01, 0.5296909721118D+00,
     :  0.2147389623610D-04, 0.2662510869850D+01, 0.1577343543434D+01,
     :  0.1680344384050D-04, 0.5207904119704D+01, 0.6279552690824D+01,
     :  0.1679117312193D-04, 0.4582187486968D+01, 0.6286599010068D+01,
     :  0.1440512068440D-04, 0.1900688517726D+01, 0.2352866153506D+01 /
      DATA ((E0(I,J,2),I=1,3),J= 11, 20) /
     :  0.1135139664999D-04, 0.5273108538556D+01, 0.5223693906222D+01,
     :  0.9345482571018D-05, 0.4503047687738D+01, 0.1203646072878D+02,
     :  0.9007418719568D-05, 0.1605621059637D+01, 0.1021328554739D+02,
     :  0.5671536712314D-05, 0.5812849070861D+00, 0.1059381944224D+01,
     :  0.7451401861666D-05, 0.2807346794836D+01, 0.3981490189893D+00,
     :  0.6393470057114D-05, 0.6029224133855D+01, 0.5753384878334D+01,
     :  0.6814275881697D-05, 0.6472990145974D+00, 0.4705732307012D+01,
     :  0.6113705628887D-05, 0.3813843419700D+01, 0.6812766822558D+01,
     :  0.4503851367273D-05, 0.4527804370996D+01, 0.5884926831456D+01,
     :  0.4522249141926D-05, 0.5991783029224D+01, 0.6256777527156D+01 /
      DATA ((E0(I,J,2),I=1,3),J= 21, 30) /
     :  0.4501794307018D-05, 0.3798703844397D+01, 0.6309374173736D+01,
     :  0.5514927480180D-05, 0.3961257833388D+01, 0.5507553240374D+01,
     :  0.4062862799995D-05, 0.5256247296369D+01, 0.6681224869435D+01,
     :  0.5414900429712D-05, 0.5499032014097D+01, 0.7755226100720D+00,
     :  0.5463153987424D-05, 0.6173092454097D+01, 0.1414349524433D+02,
     :  0.5071611859329D-05, 0.2870244247651D+01, 0.7860419393880D+01,
     :  0.2195112094455D-05, 0.2952338617201D+01, 0.1150676975667D+02,
     :  0.2279139233919D-05, 0.5951775132933D+01, 0.7058598460518D+01,
     :  0.2278386100876D-05, 0.4845456398785D+01, 0.4694002934110D+01,
     :  0.2559088003308D-05, 0.6945321117311D+00, 0.1216800268190D+02 /
      DATA ((E0(I,J,2),I=1,3),J= 31, 40) /
     :  0.2561079286856D-05, 0.6167224608301D+01, 0.7099330490126D+00,
     :  0.1792755796387D-05, 0.1400122509632D+01, 0.7962980379786D+00,
     :  0.1818715656502D-05, 0.4703347611830D+01, 0.6283142985870D+01,
     :  0.1818744924791D-05, 0.5086748900237D+01, 0.6283008715021D+01,
     :  0.1554518791390D-05, 0.5331008042713D-01, 0.2513230340178D+02,
     :  0.2063265737239D-05, 0.4283680484178D+01, 0.1179062909082D+02,
     :  0.1497613520041D-05, 0.6074207826073D+01, 0.5486777812467D+01,
     :  0.2000617940427D-05, 0.2501426281450D+01, 0.1778984560711D+02,
     :  0.1289731195580D-05, 0.3646340599536D+01, 0.7079373888424D+01,
     :  0.1282657998934D-05, 0.3232864804902D+01, 0.3738761453707D+01 /
      DATA ((E0(I,J,2),I=1,3),J= 41, 50) /
     :  0.1528915968658D-05, 0.5581433416669D+01, 0.2132990797783D+00,
     :  0.1187304098432D-05, 0.5453576453694D+01, 0.9437762937313D+01,
     :  0.7842782928118D-06, 0.2823953922273D+00, 0.8827390247185D+01,
     :  0.7352892280868D-06, 0.1124369580175D+01, 0.1589072916335D+01,
     :  0.6570189360797D-06, 0.2089154042840D+01, 0.1176985366291D+02,
     :  0.6324967590410D-06, 0.6704855581230D+00, 0.6262300422539D+01,
     :  0.6298289872283D-06, 0.2836414855840D+01, 0.6303851278352D+01,
     :  0.6476686465855D-06, 0.4852433866467D+00, 0.7113454667900D-02,
     :  0.8587034651234D-06, 0.1453511005668D+01, 0.1672837615881D+03,
     :  0.8068948788113D-06, 0.9224087798609D+00, 0.6069776770667D+01 /
      DATA ((E0(I,J,2),I=1,3),J= 51, 60) /
     :  0.8353786011661D-06, 0.4631707184895D+01, 0.3340612434717D+01,
     :  0.6009324532132D-06, 0.1829498827726D+01, 0.4136910472696D+01,
     :  0.7558158559566D-06, 0.2588596800317D+01, 0.6496374930224D+01,
     :  0.5809279504503D-06, 0.5516818853476D+00, 0.1097707878456D+02,
     :  0.5374131950254D-06, 0.6275674734960D+01, 0.1194447056968D+01,
     :  0.5711160507326D-06, 0.1091905956872D+01, 0.6282095334605D+01,
     :  0.5710183170746D-06, 0.2415001635090D+01, 0.6284056366286D+01,
     :  0.5144373590610D-06, 0.6020336443438D+01, 0.6290189305114D+01,
     :  0.5103108927267D-06, 0.3775634564605D+01, 0.6275962395778D+01,
     :  0.4960654697891D-06, 0.1073450946756D+01, 0.6127655567643D+01 /
      DATA ((E0(I,J,2),I=1,3),J= 61, 70) /
     :  0.4786385689280D-06, 0.2431178012310D+01, 0.6438496133249D+01,
     :  0.6109911263665D-06, 0.5343356157914D+01, 0.3154687086868D+01,
     :  0.4839898944024D-06, 0.5830833594047D-01, 0.8018209333619D+00,
     :  0.4734822623919D-06, 0.4536080134821D+01, 0.3128388763578D+01,
     :  0.4834741473290D-06, 0.2585090489754D+00, 0.7084896783808D+01,
     :  0.5134858581156D-06, 0.4213317172603D+01, 0.1235285262111D+02,
     :  0.5064004264978D-06, 0.4814418806478D+00, 0.1185621865188D+02,
     :  0.3753476772761D-06, 0.1599953399788D+01, 0.8429241228195D+01,
     :  0.4935264014283D-06, 0.2157417556873D+01, 0.2544314396739D+01,
     :  0.3950929600897D-06, 0.3359394184254D+01, 0.5481254917084D+01 /
      DATA ((E0(I,J,2),I=1,3),J= 71, 80) /
     :  0.4895849789777D-06, 0.5165704376558D+01, 0.9225539266174D+01,
     :  0.4215241688886D-06, 0.2065368800993D+01, 0.1726015463500D+02,
     :  0.3796773731132D-06, 0.1468606346612D+01, 0.4265981595566D+00,
     :  0.3114178142515D-06, 0.3615638079474D+01, 0.2146165377750D+01,
     :  0.3260664220838D-06, 0.4417134922435D+01, 0.4164311961999D+01,
     :  0.3976996123008D-06, 0.4700866883004D+01, 0.5856477690889D+01,
     :  0.2801459672924D-06, 0.4538902060922D+01, 0.1256967486051D+02,
     :  0.3638931868861D-06, 0.1334197991475D+01, 0.1807370494127D+02,
     :  0.2487013269476D-06, 0.3749275558275D+01, 0.2629832328990D-01,
     :  0.3034165481994D-06, 0.4236622030873D+00, 0.4535059491685D+01 /
      DATA ((E0(I,J,2),I=1,3),J= 81, 90) /
     :  0.2676278825586D-06, 0.5970848007811D+01, 0.3930209696940D+01,
     :  0.2764903818918D-06, 0.5194636754501D+01, 0.1256262854127D+02,
     :  0.2485149930507D-06, 0.1002434207846D+01, 0.5088628793478D+01,
     :  0.2199305540941D-06, 0.3066773098403D+01, 0.1255903824622D+02,
     :  0.2571106500435D-06, 0.7588312459063D+00, 0.1336797263425D+02,
     :  0.2049751817158D-06, 0.3444977434856D+01, 0.1137170464392D+02,
     :  0.2599707296297D-06, 0.1873128542205D+01, 0.7143069561767D+02,
     :  0.1785018072217D-06, 0.5015891306615D+01, 0.1748016358760D+01,
     :  0.2324833891115D-06, 0.4618271239730D+01, 0.1831953657923D+02,
     :  0.1709711119545D-06, 0.5300003455669D+01, 0.4933208510675D+01 /
      DATA ((E0(I,J,2),I=1,3),J= 91,100) /
     :  0.2107159351716D-06, 0.2229819815115D+01, 0.7477522907414D+01,
     :  0.1750333080295D-06, 0.6161485880008D+01, 0.1044738781244D+02,
     :  0.2000598210339D-06, 0.2967357299999D+01, 0.8031092209206D+01,
     :  0.1380920248681D-06, 0.3027007923917D+01, 0.8635942003952D+01,
     :  0.1412460470299D-06, 0.6037597163798D+01, 0.2942463415728D+01,
     :  0.1888459803001D-06, 0.8561476243374D+00, 0.1561374759853D+03,
     :  0.1788370542585D-06, 0.4869736290209D+01, 0.1592596075957D+01,
     :  0.1360893296167D-06, 0.3626411886436D+01, 0.1309584267300D+02,
     :  0.1506846530160D-06, 0.1550975377427D+01, 0.1649636139783D+02,
     :  0.1800913376176D-06, 0.2075826033190D+01, 0.1729818233119D+02 /
      DATA ((E0(I,J,2),I=1,3),J=101,110) /
     :  0.1436261390649D-06, 0.6148876420255D+01, 0.2042657109477D+02,
     :  0.1220227114151D-06, 0.4382583879906D+01, 0.7632943190217D+01,
     :  0.1337883603592D-06, 0.2036644327361D+01, 0.1213955354133D+02,
     :  0.1159326650738D-06, 0.3892276994687D+01, 0.5331357529664D+01,
     :  0.1352853128569D-06, 0.1447950649744D+01, 0.1673046366289D+02,
     :  0.1433408296083D-06, 0.4457854692961D+01, 0.7342457794669D+01,
     :  0.1234701666518D-06, 0.1538818147151D+01, 0.6279485555400D+01,
     :  0.1234027192007D-06, 0.1968523220760D+01, 0.6286666145492D+01,
     :  0.1244024091797D-06, 0.5779803499985D+01, 0.1511046609763D+02,
     :  0.1097934945516D-06, 0.6210975221388D+00, 0.1098880815746D+02 /
      DATA ((E0(I,J,2),I=1,3),J=111,120) /
     :  0.1254611329856D-06, 0.2591963807998D+01, 0.1572083878776D+02,
     :  0.1158247286784D-06, 0.2483612812670D+01, 0.5729506548653D+01,
     :  0.9039078252960D-07, 0.3857554579796D+01, 0.9623688285163D+01,
     :  0.9108024978836D-07, 0.5826368512984D+01, 0.7234794171227D+01,
     :  0.8887068108436D-07, 0.3475694573987D+01, 0.6148010737701D+01,
     :  0.8632374035438D-07, 0.3059070488983D-01, 0.6418140963190D+01,
     :  0.7893186992967D-07, 0.1583194837728D+01, 0.2118763888447D+01,
     :  0.8297650201172D-07, 0.8519770534637D+00, 0.1471231707864D+02,
     :  0.1019759578988D-06, 0.1319598738732D+00, 0.1349867339771D+01,
     :  0.1010037696236D-06, 0.9937860115618D+00, 0.6836645152238D+01 /
      DATA ((E0(I,J,2),I=1,3),J=121,130) /
     :  0.1047727548266D-06, 0.1382138405399D+01, 0.5999216516294D+01,
     :  0.7351993881086D-07, 0.3833397851735D+01, 0.6040347114260D+01,
     :  0.9868771092341D-07, 0.2124913814390D+01, 0.6566935184597D+01,
     :  0.7007321959390D-07, 0.5946305343763D+01, 0.6525804586632D+01,
     :  0.6861411679709D-07, 0.4574654977089D+01, 0.7238675589263D+01,
     :  0.7554519809614D-07, 0.5949232686844D+01, 0.1253985337760D+02,
     :  0.9541880448335D-07, 0.3495242990564D+01, 0.2122839202813D+02,
     :  0.7185606722155D-07, 0.4310113471661D+01, 0.6245048154254D+01,
     :  0.7131360871710D-07, 0.5480309323650D+01, 0.6321103546637D+01,
     :  0.6651142021039D-07, 0.5411097713654D+01, 0.5327476111629D+01 /
      DATA ((E0(I,J,2),I=1,3),J=131,140) /
     :  0.8538618213667D-07, 0.1827849973951D+01, 0.1101510648075D+02,
     :  0.8634954288044D-07, 0.5443584943349D+01, 0.5643178611111D+01,
     :  0.7449415051484D-07, 0.2011535459060D+01, 0.5368044267797D+00,
     :  0.7421047599169D-07, 0.3464562529249D+01, 0.2354323048545D+02,
     :  0.6140694354424D-07, 0.5657556228815D+01, 0.1296430071988D+02,
     :  0.6353525143033D-07, 0.3463816593821D+01, 0.1990745094947D+01,
     :  0.6221964013447D-07, 0.1532259498697D+01, 0.9517183207817D+00,
     :  0.5852480257244D-07, 0.1375396598875D+01, 0.9555997388169D+00,
     :  0.6398637498911D-07, 0.2405645801972D+01, 0.2407292145756D+02,
     :  0.7039744069878D-07, 0.5397541799027D+01, 0.5225775174439D+00 /
      DATA ((E0(I,J,2),I=1,3),J=141,150) /
     :  0.6977997694382D-07, 0.4762347105419D+01, 0.1097355562493D+02,
     :  0.7460629558396D-07, 0.2711944692164D+01, 0.2200391463820D+02,
     :  0.5376577536101D-07, 0.2352980430239D+01, 0.1431416805965D+02,
     :  0.7530607893556D-07, 0.1943940180699D+01, 0.1842262939178D+02,
     :  0.6822928971605D-07, 0.4337651846959D+01, 0.1554202828031D+00,
     :  0.6220772380094D-07, 0.6716871369278D+00, 0.1845107853235D+02,
     :  0.6586950799043D-07, 0.2229714460505D+01, 0.5216580451554D+01,
     :  0.5873800565771D-07, 0.7627013920580D+00, 0.6398972393349D+00,
     :  0.6264346929745D-07, 0.6202785478961D+00, 0.6277552955062D+01,
     :  0.6257929115669D-07, 0.2886775596668D+01, 0.6288598745829D+01 /
      DATA ((E0(I,J,2),I=1,3),J=151,160) /
     :  0.5343536033409D-07, 0.1977241012051D+01, 0.4690479774488D+01,
     :  0.5587849781714D-07, 0.1922923484825D+01, 0.1551045220144D+01,
     :  0.6905100845603D-07, 0.3570757164631D+01, 0.1030928125552D+00,
     :  0.6178957066649D-07, 0.5197558947765D+01, 0.5230807360890D+01,
     :  0.6187270224331D-07, 0.8193497368922D+00, 0.5650292065779D+01,
     :  0.5385664291426D-07, 0.5406336665586D+01, 0.7771377146812D+02,
     :  0.6329363917926D-07, 0.2837760654536D+01, 0.2608790314060D+02,
     :  0.4546018761604D-07, 0.2933580297050D+01, 0.5535693017924D+00,
     :  0.6196091049375D-07, 0.4157871494377D+01, 0.8467247584405D+02,
     :  0.6159555108218D-07, 0.3211703561703D+01, 0.2394243902548D+03 /
      DATA ((E0(I,J,2),I=1,3),J=161,170) /
     :  0.4995340539317D-07, 0.1459098102922D+01, 0.4732030630302D+01,
     :  0.5457031243572D-07, 0.1430457676136D+01, 0.6179983037890D+01,
     :  0.4863461418397D-07, 0.2196425916730D+01, 0.9027992316901D+02,
     :  0.5342947626870D-07, 0.2086612890268D+01, 0.6386168663001D+01,
     :  0.5674296648439D-07, 0.2760204966535D+01, 0.6915859635113D+01,
     :  0.4745783120161D-07, 0.4245368971862D+01, 0.6282970628506D+01,
     :  0.4745676961198D-07, 0.5544725787016D+01, 0.6283181072386D+01,
     :  0.4049796869973D-07, 0.2213984363586D+01, 0.6254626709878D+01,
     :  0.4248333596940D-07, 0.8075781952896D+00, 0.7875671926403D+01,
     :  0.4027178070205D-07, 0.1293268540378D+01, 0.6311524991013D+01 /
      DATA ((E0(I,J,2),I=1,3),J=171,180) /
     :  0.4066543943476D-07, 0.3986141175804D+01, 0.3634620989887D+01,
     :  0.4858863787880D-07, 0.1276112738231D+01, 0.5760498333002D+01,
     :  0.5277398263530D-07, 0.4916111741527D+01, 0.2515860172507D+02,
     :  0.4105635656559D-07, 0.1725805864426D+01, 0.6709674010002D+01,
     :  0.4376781925772D-07, 0.2243642442106D+01, 0.6805653367890D+01,
     :  0.3235827894693D-07, 0.3614135118271D+01, 0.1066495398892D+01,
     :  0.3073244740308D-07, 0.2460873393460D+01, 0.5863591145557D+01,
     :  0.3088609271373D-07, 0.5678431771790D+01, 0.9917696840332D+01,
     :  0.3393022279836D-07, 0.3814017477291D+01, 0.1391601904066D+02,
     :  0.3038686508802D-07, 0.4660216229171D+01, 0.1256621883632D+02 /
      DATA ((E0(I,J,2),I=1,3),J=181,190) /
     :  0.4019677752497D-07, 0.5906906243735D+01, 0.1334167431096D+02,
     :  0.3288834998232D-07, 0.9536146445882D+00, 0.1620077269078D+02,
     :  0.3889973794631D-07, 0.3942205097644D+01, 0.7478166569050D-01,
     :  0.3050438987141D-07, 0.1624810271286D+01, 0.1805292951336D+02,
     :  0.3601142564638D-07, 0.4030467142575D+01, 0.6208294184755D+01,
     :  0.3689015557141D-07, 0.3648878818694D+01, 0.5966683958112D+01,
     :  0.3563471893565D-07, 0.5749584017096D+01, 0.6357857516136D+01,
     :  0.2776183170667D-07, 0.2630124187070D+01, 0.3523159621801D-02,
     :  0.2922350530341D-07, 0.1790346403629D+01, 0.1272157198369D+02,
     :  0.3511076917302D-07, 0.6142198301611D+01, 0.6599467742779D+01 /
      DATA ((E0(I,J,2),I=1,3),J=191,200) /
     :  0.3619351007632D-07, 0.1432421386492D+01, 0.6019991944201D+01,
     :  0.2561254711098D-07, 0.2302822475792D+01, 0.1259245002418D+02,
     :  0.2626903942920D-07, 0.8660470994571D+00, 0.6702560555334D+01,
     :  0.2550187397083D-07, 0.6069721995383D+01, 0.1057540660594D+02,
     :  0.2535873526138D-07, 0.1079020331795D-01, 0.3141537925223D+02,
     :  0.3519786153847D-07, 0.3809066902283D+01, 0.2505706758577D+03,
     :  0.3424651492873D-07, 0.2075435114417D+01, 0.6546159756691D+01,
     :  0.2372676630861D-07, 0.2057803120154D+01, 0.2388894113936D+01,
     :  0.2710980779541D-07, 0.1510068488010D+01, 0.1202934727411D+02,
     :  0.3038710889704D-07, 0.5043617528901D+01, 0.1256608456547D+02 /
      DATA ((E0(I,J,2),I=1,3),J=201,210) /
     :  0.2220364130585D-07, 0.3694793218205D+01, 0.1336244973887D+02,
     :  0.3025880825460D-07, 0.5450618999049D-01, 0.2908881142201D+02,
     :  0.2784493486864D-07, 0.3381164084502D+01, 0.1494531617769D+02,
     :  0.2294414142438D-07, 0.4382309025210D+01, 0.6076890225335D+01,
     :  0.2012723294724D-07, 0.9142212256518D+00, 0.6262720680387D+01,
     :  0.2036357831958D-07, 0.5676172293154D+01, 0.4701116388778D+01,
     :  0.2003474823288D-07, 0.2592767977625D+01, 0.6303431020504D+01,
     :  0.2207144900109D-07, 0.5404976271180D+01, 0.6489261475556D+01,
     :  0.2481664905135D-07, 0.4373284587027D+01, 0.1204357418345D+02,
     :  0.2674949182295D-07, 0.5859182188482D+01, 0.4590910121555D+01 /
      DATA ((E0(I,J,2),I=1,3),J=211,220) /
     :  0.2450554720322D-07, 0.4555381557451D+01, 0.1495633313810D+00,
     :  0.2601975986457D-07, 0.3933165584959D+01, 0.1965104848470D+02,
     :  0.2199860022848D-07, 0.5227977189087D+01, 0.1351787002167D+02,
     :  0.2448121172316D-07, 0.4858060353949D+01, 0.1162474756779D+01,
     :  0.1876014864049D-07, 0.5690546553605D+01, 0.6279194432410D+01,
     :  0.1874513219396D-07, 0.4099539297446D+01, 0.6286957268481D+01,
     :  0.2156380842559D-07, 0.4382594769913D+00, 0.1813929450232D+02,
     :  0.1981691240061D-07, 0.1829784152444D+01, 0.4686889479442D+01,
     :  0.2329992648539D-07, 0.2836254278973D+01, 0.1002183730415D+02,
     :  0.1765184135302D-07, 0.2803494925833D+01, 0.4292330755499D+01 /
      DATA ((E0(I,J,2),I=1,3),J=221,230) /
     :  0.2436368366085D-07, 0.2836897959677D+01, 0.9514313292143D+02,
     :  0.2164089203889D-07, 0.6127522446024D+01, 0.6037244212485D+01,
     :  0.1847755034221D-07, 0.3683163635008D+01, 0.2427287361862D+00,
     :  0.1674798769966D-07, 0.3316993867246D+00, 0.1311972100268D+02,
     :  0.2222542124356D-07, 0.8294097805480D+00, 0.1266924451345D+02,
     :  0.2071074505925D-07, 0.3659492220261D+01, 0.6528907488406D+01,
     :  0.1608224471835D-07, 0.4774492067182D+01, 0.1352175143971D+02,
     :  0.1857583439071D-07, 0.2873120597682D+01, 0.8662240327241D+01,
     :  0.1793018836159D-07, 0.5282441177929D+00, 0.6819880277225D+01,
     :  0.1575391221692D-07, 0.1320789654258D+01, 0.1102062672231D+00 /
      DATA ((E0(I,J,2),I=1,3),J=231,240) /
     :  0.1840132009557D-07, 0.1917110916256D+01, 0.6514761976723D+02,
     :  0.1760917288281D-07, 0.2972635937132D+01, 0.5746271423666D+01,
     :  0.1561779518516D-07, 0.4372569261981D+01, 0.6272439236156D+01,
     :  0.1558687885205D-07, 0.5416424926425D+01, 0.6293712464735D+01,
     :  0.1951359382579D-07, 0.3094448898752D+01, 0.2301353951334D+02,
     :  0.1569144275614D-07, 0.2802103689808D+01, 0.1765478049437D+02,
     :  0.1479130389462D-07, 0.2136435020467D+01, 0.2077542790660D-01,
     :  0.1467828510764D-07, 0.7072627435674D+00, 0.1052268489556D+01,
     :  0.1627627337440D-07, 0.3947607143237D+01, 0.6327837846670D+00,
     :  0.1503498479758D-07, 0.4079248909190D+01, 0.7626583626240D-01 /
      DATA ((E0(I,J,2),I=1,3),J=241,250) /
     :  0.1297967708237D-07, 0.6269637122840D+01, 0.1149965630200D+02,
     :  0.1374416896634D-07, 0.4175657970702D+01, 0.6016468784579D+01,
     :  0.1783812325219D-07, 0.1476540547560D+01, 0.3301902111895D+02,
     :  0.1525884228756D-07, 0.4653477715241D+01, 0.9411464614024D+01,
     :  0.1451067396763D-07, 0.2573001128225D+01, 0.1277945078067D+02,
     :  0.1297713111950D-07, 0.5612799618771D+01, 0.6549682916313D+01,
     :  0.1462784012820D-07, 0.4189661623870D+01, 0.1863592847156D+02,
     :  0.1384185980007D-07, 0.2656915472196D+01, 0.2379164476796D+01,
     :  0.1221497599801D-07, 0.5612515760138D+01, 0.1257326515556D+02,
     :  0.1560574525896D-07, 0.4783414317919D+01, 0.1887552587463D+02 /
      DATA ((E0(I,J,2),I=1,3),J=251,260) /
     :  0.1544598372036D-07, 0.2694431138063D+01, 0.1820933031200D+02,
     :  0.1531678928696D-07, 0.4105103489666D+01, 0.2593412433514D+02,
     :  0.1349321503795D-07, 0.3082437194015D+00, 0.5120601093667D+01,
     :  0.1252030290917D-07, 0.6124072334087D+01, 0.6993008899458D+01,
     :  0.1459243816687D-07, 0.3733103981697D+01, 0.3813291813120D-01,
     :  0.1226103625262D-07, 0.1267127706817D+01, 0.2435678079171D+02,
     :  0.1019449641504D-07, 0.4367790112269D+01, 0.1725663147538D+02,
     :  0.1380789433607D-07, 0.3387201768700D+01, 0.2458316379602D+00,
     :  0.1019453421658D-07, 0.9204143073737D+00, 0.6112403035119D+01,
     :  0.1297929434405D-07, 0.5786874896426D+01, 0.1249137003520D+02 /
      DATA ((E0(I,J,2),I=1,3),J=261,270) /
     :  0.9912677786097D-08, 0.3164232870746D+01, 0.6247047890016D+01,
     :  0.9829386098599D-08, 0.2586762413351D+01, 0.6453748665772D+01,
     :  0.1226807746104D-07, 0.6239068436607D+01, 0.5429879531333D+01,
     :  0.1192691755997D-07, 0.1867380051424D+01, 0.6290122169689D+01,
     :  0.9836499227081D-08, 0.3424716293727D+00, 0.6319103810876D+01,
     :  0.9642862564285D-08, 0.5661372990657D+01, 0.8273820945392D+01,
     :  0.1165184404862D-07, 0.5768367239093D+01, 0.1778273215245D+02,
     :  0.1175794418818D-07, 0.1657351222943D+01, 0.6276029531202D+01,
     :  0.1018948635601D-07, 0.6458292350865D+00, 0.1254537627298D+02,
     :  0.9500383606676D-08, 0.1054306140741D+01, 0.1256517118505D+02 /
      DATA ((E0(I,J,2),I=1,3),J=271,280) /
     :  0.1227512202906D-07, 0.2505278379114D+01, 0.2248384854122D+02,
     :  0.9664792009993D-08, 0.4289737277000D+01, 0.6259197520765D+01,
     :  0.9613285666331D-08, 0.5500597673141D+01, 0.6306954180126D+01,
     :  0.1117906736211D-07, 0.2361405953468D+01, 0.1779695906178D+02,
     :  0.9611378640782D-08, 0.2851310576269D+01, 0.2061856251104D+00,
     :  0.8845354852370D-08, 0.6208777705343D+01, 0.1692165728891D+01,
     :  0.1054046966600D-07, 0.5413091423934D+01, 0.2204125344462D+00,
     :  0.1215539124483D-07, 0.5613969479755D+01, 0.8257698122054D+02,
     :  0.9932460955209D-08, 0.1106124877015D+01, 0.1017725758696D+02,
     :  0.8785804715043D-08, 0.2869224476477D+01, 0.9491756770005D+00 /
      DATA ((E0(I,J,2),I=1,3),J=281,290) /
     :  0.8538084097562D-08, 0.6159640899344D+01, 0.6393282117669D+01,
     :  0.8648994369529D-08, 0.1374901198784D+01, 0.4804209201333D+01,
     :  0.1039063219067D-07, 0.5171080641327D+01, 0.1550861511662D+02,
     :  0.8867983926439D-08, 0.8317320304902D+00, 0.3903911373650D+01,
     :  0.8327495955244D-08, 0.3605591969180D+01, 0.6172869583223D+01,
     :  0.9243088356133D-08, 0.6114299196843D+01, 0.6267823317922D+01,
     :  0.9205657357835D-08, 0.3675153683737D+01, 0.6298328382969D+01,
     :  0.1033269714606D-07, 0.3313328813024D+01, 0.5573142801433D+01,
     :  0.8001706275552D-08, 0.2019980960053D+01, 0.2648454860559D+01,
     :  0.9171858254191D-08, 0.8992015524177D+00, 0.1498544001348D+03 /
      DATA ((E0(I,J,2),I=1,3),J=291,300) /
     :  0.1075327150242D-07, 0.2898669963648D+01, 0.3694923081589D+02,
     :  0.9884866689828D-08, 0.4946715904478D+01, 0.1140367694411D+02,
     :  0.9541835576677D-08, 0.2371787888469D+01, 0.1256713221673D+02,
     :  0.7739903376237D-08, 0.2213775190612D+01, 0.7834121070590D+01,
     :  0.7311962684106D-08, 0.3429378787739D+01, 0.1192625446156D+02,
     :  0.9724904869624D-08, 0.6195878564404D+01, 0.2280573557157D+02,
     :  0.9251628983612D-08, 0.6511509527390D+00, 0.2787043132925D+01,
     :  0.7320763787842D-08, 0.6001083639421D+01, 0.6282655592598D+01,
     :  0.7320296650962D-08, 0.3789073265087D+01, 0.6283496108294D+01,
     :  0.7947032271039D-08, 0.1059659582204D+01, 0.1241073141809D+02 /
      DATA ((E0(I,J,2),I=1,3),J=301,310) /
     :  0.9005277053115D-08, 0.1280315624361D+01, 0.6281591679874D+01,
     :  0.8995601652048D-08, 0.2224439106766D+01, 0.6284560021018D+01,
     :  0.8288040568796D-08, 0.5234914433867D+01, 0.1241658836951D+02,
     :  0.6359381347255D-08, 0.4137989441490D+01, 0.1596186371003D+01,
     :  0.8699572228626D-08, 0.1758411009497D+01, 0.6133512519065D+01,
     :  0.6456797542736D-08, 0.5919285089994D+01, 0.1685848245639D+02,
     :  0.7424573475452D-08, 0.5414616938827D+01, 0.4061219149443D+01,
     :  0.7235671196168D-08, 0.1496516557134D+01, 0.1610006857377D+03,
     :  0.8104015182733D-08, 0.1919918242764D+01, 0.8460828644453D+00,
     :  0.8098576535937D-08, 0.3819615855458D+01, 0.3894181736510D+01 /
      DATA ((E0(I,J,2),I=1,3),J=311,320) /
     :  0.6275292346625D-08, 0.6244264115141D+01, 0.8531963191132D+00,
     :  0.6052432989112D-08, 0.5037731872610D+00, 0.1567108171867D+02,
     :  0.5705651535817D-08, 0.2984557271995D+01, 0.1258692712880D+02,
     :  0.5789650115138D-08, 0.6087038140697D+01, 0.1193336791622D+02,
     :  0.5512132153377D-08, 0.5855668994076D+01, 0.1232342296471D+02,
     :  0.7388890819102D-08, 0.2443128574740D+01, 0.4907302013889D+01,
     :  0.5467593991798D-08, 0.3017561234194D+01, 0.1884211409667D+02,
     :  0.6388519802999D-08, 0.5887386712935D+01, 0.5217580628120D+02,
     :  0.6106777149944D-08, 0.3483461059895D+00, 0.1422690933580D-01,
     :  0.7383420275489D-08, 0.5417387056707D+01, 0.2358125818164D+02 /
      DATA ((E0(I,J,2),I=1,3),J=321,330) /
     :  0.5505208141738D-08, 0.2848193644783D+01, 0.1151388321134D+02,
     :  0.6310757462877D-08, 0.2349882520828D+01, 0.1041998632314D+02,
     :  0.6166904929691D-08, 0.5728575944077D+00, 0.6151533897323D+01,
     :  0.5263442042754D-08, 0.4495796125937D+01, 0.1885275071096D+02,
     :  0.5591828082629D-08, 0.1355441967677D+01, 0.4337116142245D+00,
     :  0.5397051680497D-08, 0.1673422864307D+01, 0.6286362197481D+01,
     :  0.5396992745159D-08, 0.1833502206373D+01, 0.6279789503410D+01,
     :  0.6572913000726D-08, 0.3331122065824D+01, 0.1176433076753D+02,
     :  0.5123421866413D-08, 0.2165327142679D+01, 0.1245594543367D+02,
     :  0.5930495725999D-08, 0.2931146089284D+01, 0.6414617803568D+01 /
      DATA ((E0(I,J,2),I=1,3),J=331,340) /
     :  0.6431797403933D-08, 0.4134407994088D+01, 0.1350651127443D+00,
     :  0.5003182207604D-08, 0.3805420303749D+01, 0.1096996532989D+02,
     :  0.5587731032504D-08, 0.1082469260599D+01, 0.6062663316000D+01,
     :  0.5935263407816D-08, 0.8384333678401D+00, 0.5326786718777D+01,
     :  0.4756019827760D-08, 0.3552588749309D+01, 0.3104930017775D+01,
     :  0.6599951172637D-08, 0.4320826409528D+01, 0.4087944051283D+02,
     :  0.5902606868464D-08, 0.4811879454445D+01, 0.5849364236221D+01,
     :  0.5921147809031D-08, 0.9942628922396D-01, 0.1581959461667D+01,
     :  0.5505382581266D-08, 0.2466557607764D+01, 0.6503488384892D+01,
     :  0.5353771071862D-08, 0.4551978748683D+01, 0.1735668374386D+03 /
      DATA ((E0(I,J,2),I=1,3),J=341,350) /
     :  0.5063282210946D-08, 0.5710812312425D+01, 0.1248988586463D+02,
     :  0.5926120403383D-08, 0.1333998428358D+01, 0.2673594526851D+02,
     :  0.5211016176149D-08, 0.4649315360760D+01, 0.2460261242967D+02,
     :  0.5347075084894D-08, 0.5512754081205D+01, 0.4171425416666D+01,
     :  0.4872609773574D-08, 0.1308025299938D+01, 0.5333900173445D+01,
     :  0.4727711321420D-08, 0.2144908368062D+01, 0.7232251527446D+01,
     :  0.6029426018652D-08, 0.5567259412084D+01, 0.3227113045244D+03,
     :  0.4321485284369D-08, 0.5230667156451D+01, 0.9388005868221D+01,
     :  0.4476406760553D-08, 0.6134081115303D+01, 0.5547199253223D+01,
     :  0.5835268277420D-08, 0.4783808492071D+01, 0.7285056171570D+02 /
      DATA ((E0(I,J,2),I=1,3),J=351,360) /
     :  0.5172183602748D-08, 0.5161817911099D+01, 0.1884570439172D+02,
     :  0.5693571465184D-08, 0.1381646203111D+01, 0.9723862754494D+02,
     :  0.4060634965349D-08, 0.3876705259495D+00, 0.4274518229222D+01,
     :  0.3967398770473D-08, 0.5029491776223D+01, 0.3496032717521D+01,
     :  0.3943754005255D-08, 0.1923162955490D+01, 0.6244942932314D+01,
     :  0.4781323427824D-08, 0.4633332586423D+01, 0.2929661536378D+02,
     :  0.3871483781204D-08, 0.1616650009743D+01, 0.6321208768577D+01,
     :  0.5141741733997D-08, 0.9817316704659D-01, 0.1232032006293D+02,
     :  0.4002385978497D-08, 0.3656161212139D+01, 0.7018952447668D+01,
     :  0.4901092604097D-08, 0.4404098713092D+01, 0.1478866649112D+01 /
      DATA ((E0(I,J,2),I=1,3),J=361,370) /
     :  0.3740932630345D-08, 0.5181188732639D+00, 0.6922973089781D+01,
     :  0.4387283718538D-08, 0.3254859566869D+01, 0.2331413144044D+03,
     :  0.5019197802033D-08, 0.3086773224677D+01, 0.1715706182245D+02,
     :  0.3834931695175D-08, 0.2797882673542D+01, 0.1491901785440D+02,
     :  0.3760413942497D-08, 0.2892676280217D+01, 0.1726726808967D+02,
     :  0.3719717204628D-08, 0.5861046025739D+01, 0.6297302759782D+01,
     :  0.4145623530149D-08, 0.2168239627033D+01, 0.1376059875786D+02,
     :  0.3932788425380D-08, 0.6271811124181D+01, 0.7872148766781D+01,
     :  0.3686377476857D-08, 0.3936853151404D+01, 0.6268848941110D+01,
     :  0.3779077950339D-08, 0.1404148734043D+01, 0.4157198507331D+01 /
      DATA ((E0(I,J,2),I=1,3),J=371,380) /
     :  0.4091334550598D-08, 0.2452436180854D+01, 0.9779108567966D+01,
     :  0.3926694536146D-08, 0.6102292739040D+01, 0.1098419223922D+02,
     :  0.4841000253289D-08, 0.6072760457276D+01, 0.1252801878276D+02,
     :  0.4949340130240D-08, 0.1154832815171D+01, 0.1617106187867D+03,
     :  0.3761557737360D-08, 0.5527545321897D+01, 0.3185192151914D+01,
     :  0.3647396268188D-08, 0.1525035688629D+01, 0.6271346477544D+01,
     :  0.3932405074189D-08, 0.5570681040569D+01, 0.2139354194808D+02,
     :  0.3631322501141D-08, 0.1981240601160D+01, 0.6294805223347D+01,
     :  0.4130007425139D-08, 0.2050060880201D+01, 0.2195415756911D+02,
     :  0.4433905965176D-08, 0.3277477970321D+01, 0.7445550607224D+01 /
      DATA ((E0(I,J,2),I=1,3),J=381,390) /
     :  0.3851814176947D-08, 0.5210690074886D+01, 0.9562891316684D+00,
     :  0.3485807052785D-08, 0.6653274904611D+00, 0.1161697602389D+02,
     :  0.3979772816991D-08, 0.1767941436148D+01, 0.2277943724828D+02,
     :  0.3402607460500D-08, 0.3421746306465D+01, 0.1087398597200D+02,
     :  0.4049993000926D-08, 0.1127144787547D+01, 0.3163918923335D+00,
     :  0.3420511182382D-08, 0.4214794779161D+01, 0.1362553364512D+02,
     :  0.3640772365012D-08, 0.5324905497687D+01, 0.1725304118033D+02,
     :  0.3323037987501D-08, 0.6135761838271D+01, 0.6279143387820D+01,
     :  0.4503141663637D-08, 0.1802305450666D+01, 0.1385561574497D+01,
     :  0.4314560055588D-08, 0.4812299731574D+01, 0.4176041334900D+01 /
      DATA ((E0(I,J,2),I=1,3),J=391,400) /
     :  0.3294226949110D-08, 0.3657547059723D+01, 0.6287008313071D+01,
     :  0.3215657197281D-08, 0.4866676894425D+01, 0.5749861718712D+01,
     :  0.4129362656266D-08, 0.3809342558906D+01, 0.5905702259363D+01,
     :  0.3137762976388D-08, 0.2494635174443D+01, 0.2099539292909D+02,
     :  0.3514010952384D-08, 0.2699961831678D+01, 0.7335344340001D+01,
     :  0.3327607571530D-08, 0.3318457714816D+01, 0.5436992986000D+01,
     :  0.3541066946675D-08, 0.4382703582466D+01, 0.1234573916645D+02,
     :  0.3216179847052D-08, 0.5271066317054D+01, 0.3802769619140D-01,
     :  0.2959045059570D-08, 0.5819591585302D+01, 0.2670964694522D+02,
     :  0.3884040326665D-08, 0.5980934960428D+01, 0.6660449441528D+01 /
      DATA ((E0(I,J,2),I=1,3),J=401,410) /
     :  0.2922027539886D-08, 0.3337290282483D+01, 0.1375773836557D+01,
     :  0.4110846382042D-08, 0.5742978187327D+01, 0.4480965020977D+02,
     :  0.2934508411032D-08, 0.2278075804200D+01, 0.6408777551755D+00,
     :  0.3966896193000D-08, 0.5835747858477D+01, 0.3773735910827D+00,
     :  0.3286695827610D-08, 0.5838898193902D+01, 0.3932462625300D-02,
     :  0.3720643094196D-08, 0.1122212337858D+01, 0.1646033343740D+02,
     :  0.3285508906174D-08, 0.9182250996416D+00, 0.1081813534213D+02,
     :  0.3753880575973D-08, 0.5174761973266D+01, 0.5642198095270D+01,
     :  0.3022129385587D-08, 0.3381611020639D+01, 0.2982630633589D+02,
     :  0.2798569205621D-08, 0.3546193723922D+01, 0.1937891852345D+02 /
      DATA ((E0(I,J,2),I=1,3),J=411,420) /
     :  0.3397872070505D-08, 0.4533203197934D+01, 0.6923953605621D+01,
     :  0.3708099772977D-08, 0.2756168198616D+01, 0.3066615496545D+02,
     :  0.3599283541510D-08, 0.1934395469918D+01, 0.6147450479709D+01,
     :  0.3688702753059D-08, 0.7149920971109D+00, 0.2636725487657D+01,
     :  0.2681084724003D-08, 0.4899819493154D+01, 0.6816289982179D+01,
     :  0.3495993460759D-08, 0.1572418915115D+01, 0.6418701221183D+01,
     :  0.3130770324995D-08, 0.8912190180489D+00, 0.1235996607578D+02,
     :  0.2744353821941D-08, 0.3800821940055D+01, 0.2059724391010D+02,
     :  0.2842732906341D-08, 0.2644717440029D+01, 0.2828699048865D+02,
     :  0.3046882682154D-08, 0.3987793020179D+01, 0.6055599646783D+01 /
      DATA ((E0(I,J,2),I=1,3),J=421,430) /
     :  0.2399072455143D-08, 0.9908826440764D+00, 0.6255674361143D+01,
     :  0.2384306274204D-08, 0.2516149752220D+01, 0.6310477339748D+01,
     :  0.2977324500559D-08, 0.5849195642118D+01, 0.1652265972112D+02,
     :  0.3062835258972D-08, 0.1681660100162D+01, 0.1172006883645D+02,
     :  0.3109682589231D-08, 0.5804143987737D+00, 0.2751146787858D+02,
     :  0.2903920355299D-08, 0.5800768280123D+01, 0.6510552054109D+01,
     :  0.2823221989212D-08, 0.9241118370216D+00, 0.5469525544182D+01,
     :  0.3187949696649D-08, 0.3139776445735D+01, 0.1693792562116D+03,
     :  0.2922559771655D-08, 0.3549440782984D+01, 0.2630839062450D+00,
     :  0.2436302066603D-08, 0.4735540696319D+01, 0.3946258593675D+00 /
      DATA ((E0(I,J,2),I=1,3),J=431,440) /
     :  0.3049473043606D-08, 0.4998289124561D+01, 0.8390110365991D+01,
     :  0.2863682575784D-08, 0.6709515671102D+00, 0.2243449970715D+00,
     :  0.2641750517966D-08, 0.5410978257284D+01, 0.2986433403208D+02,
     :  0.2704093466243D-08, 0.4778317207821D+01, 0.6129297044991D+01,
     :  0.2445522177011D-08, 0.6009020662222D+01, 0.1171295538178D+02,
     :  0.2623608810230D-08, 0.5010449777147D+01, 0.6436854655901D+01,
     :  0.2079259704053D-08, 0.5980943768809D+01, 0.2019909489111D+02,
     :  0.2820225596771D-08, 0.2679965110468D+01, 0.5934151399930D+01,
     :  0.2365221950927D-08, 0.1894231148810D+01, 0.2470570524223D+02,
     :  0.2359682077149D-08, 0.4220752950780D+01, 0.8671969964381D+01 /
      DATA ((E0(I,J,2),I=1,3),J=441,450) /
     :  0.2387577137206D-08, 0.2571783940617D+01, 0.7096626156709D+01,
     :  0.1982102089816D-08, 0.5169765997119D+00, 0.1727188400790D+02,
     :  0.2687502389925D-08, 0.6239078264579D+01, 0.7075506709219D+02,
     :  0.2207751669135D-08, 0.2031184412677D+01, 0.4377611041777D+01,
     :  0.2618370214274D-08, 0.8266079985979D+00, 0.6632000300961D+01,
     :  0.2591951887361D-08, 0.8819350522008D+00, 0.4873985990671D+02,
     :  0.2375055656248D-08, 0.3520944177789D+01, 0.1590676413561D+02,
     :  0.2472019978911D-08, 0.1551431908671D+01, 0.6612329252343D+00,
     :  0.2368157127199D-08, 0.4178610147412D+01, 0.3459636466239D+02,
     :  0.1764846605693D-08, 0.1506764000157D+01, 0.1980094587212D+02 /
      DATA ((E0(I,J,2),I=1,3),J=451,460) /
     :  0.2291769608798D-08, 0.2118250611782D+01, 0.2844914056730D-01,
     :  0.2209997316943D-08, 0.3363255261678D+01, 0.2666070658668D+00,
     :  0.2292699097923D-08, 0.4200423956460D+00, 0.1484170571900D-02,
     :  0.1629683015329D-08, 0.2331362582487D+01, 0.3035599730800D+02,
     :  0.2206492862426D-08, 0.3400274026992D+01, 0.6281667977667D+01,
     :  0.2205746568257D-08, 0.1066051230724D+00, 0.6284483723224D+01,
     :  0.2026310767991D-08, 0.2779066487979D+01, 0.2449240616245D+02,
     :  0.1762977622163D-08, 0.9951450691840D+00, 0.2045286941806D+02,
     :  0.1368535049606D-08, 0.6402447365817D+00, 0.2473415438279D+02,
     :  0.1720598775450D-08, 0.2303524214705D+00, 0.1679593901136D+03 /
      DATA ((E0(I,J,2),I=1,3),J=461,470) /
     :  0.1702429015449D-08, 0.6164622655048D+01, 0.3338575901272D+03,
     :  0.1414033197685D-08, 0.3954561185580D+01, 0.1624205518357D+03,
     :  0.1573768958043D-08, 0.2028286308984D+01, 0.3144167757552D+02,
     :  0.1650705184447D-08, 0.2304040666128D+01, 0.5267006960365D+02,
     :  0.1651087618855D-08, 0.2538461057280D+01, 0.8956999012000D+02,
     :  0.1616409518983D-08, 0.5111054348152D+01, 0.3332657872986D+02,
     :  0.1537175173581D-08, 0.5601130666603D+01, 0.3852657435933D+02,
     :  0.1593191980553D-08, 0.2614340453411D+01, 0.2282781046519D+03,
     :  0.1499480170643D-08, 0.3624721577264D+01, 0.2823723341956D+02,
     :  0.1493807843235D-08, 0.4214569879008D+01, 0.2876692439167D+02 /
      DATA ((E0(I,J,2),I=1,3),J=471,480) /
     :  0.1074571199328D-08, 0.1496911744704D+00, 0.8397383534231D+02,
     :  0.1074406983417D-08, 0.1187817671922D+01, 0.8401985929482D+02,
     :  0.9757576855851D-09, 0.2655703035858D+01, 0.7826370942180D+02,
     :  0.1258432887565D-08, 0.4969896184844D+01, 0.3115650189215D+03,
     :  0.1240336343282D-08, 0.5192460776926D+01, 0.1784300471910D+03,
     :  0.9016107005164D-09, 0.1960356923057D+01, 0.5886454391678D+02,
     :  0.1135392360918D-08, 0.5082427809068D+01, 0.7842370451713D+02,
     :  0.9216046089565D-09, 0.2793775037273D+01, 0.1014262087719D+03,
     :  0.1061276615030D-08, 0.3726144311409D+01, 0.5660027930059D+02,
     :  0.1010110596263D-08, 0.7404080708937D+00, 0.4245678405627D+02 /
      DATA ((E0(I,J,2),I=1,3),J=481,490) /
     :  0.7217424756199D-09, 0.2697449980577D-01, 0.2457074661053D+03,
     :  0.6912003846756D-09, 0.4253296276335D+01, 0.1679936946371D+03,
     :  0.6871814664847D-09, 0.5148072412354D+01, 0.6053048899753D+02,
     :  0.4887158016343D-09, 0.2153581148294D+01, 0.9656299901946D+02,
     :  0.5161802866314D-09, 0.3852750634351D+01, 0.2442876000072D+03,
     :  0.5652599559057D-09, 0.1233233356270D+01, 0.8365903305582D+02,
     :  0.4710812608586D-09, 0.5610486976767D+01, 0.3164282286739D+03,
     :  0.4909977500324D-09, 0.1639629524123D+01, 0.4059982187939D+03,
     :  0.4772641839378D-09, 0.3737100368583D+01, 0.1805255418145D+03,
     :  0.4487562567153D-09, 0.1158417054478D+00, 0.8433466158131D+02 /
      DATA ((E0(I,J,2),I=1,3),J=491,500) /
     :  0.3943441230497D-09, 0.6243502862796D+00, 0.2568537517081D+03,
     :  0.3952236913598D-09, 0.3510377382385D+01, 0.2449975330562D+03,
     :  0.3788898363417D-09, 0.5916128302299D+01, 0.1568131045107D+03,
     :  0.3738329328831D-09, 0.1042266763456D+01, 0.3948519331910D+03,
     :  0.2451199165151D-09, 0.1166788435700D+01, 0.1435713242844D+03,
     :  0.2436734402904D-09, 0.3254726114901D+01, 0.2268582385539D+03,
     :  0.2213605274325D-09, 0.1687210598530D+01, 0.1658638954901D+03,
     :  0.1491521204829D-09, 0.2657541786794D+01, 0.2219950288015D+03,
     :  0.1474995329744D-09, 0.5013089805819D+01, 0.3052819430710D+03,
     :  0.1661939475656D-09, 0.5495315428418D+01, 0.2526661704812D+03 /
      DATA ((E0(I,J,2),I=1,3),J=501,NE0Y) /
     :  0.9015946748003D-10, 0.2236989966505D+01, 0.4171445043968D+03 /

*  Sun-to-Earth, T^1, Y
      DATA ((E1(I,J,2),I=1,3),J=  1, 10) /
     :  0.9304690546528D-06, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.5150715570663D-06, 0.4431807116294D+01, 0.1256615170089D+02,
     :  0.1290825411056D-07, 0.4388610039678D+01, 0.1884922755134D+02,
     :  0.4645466665386D-08, 0.5827263376034D+01, 0.6283075850446D+01,
     :  0.2079625310718D-08, 0.1621698662282D+00, 0.6279552690824D+01,
     :  0.2078189850907D-08, 0.3344713435140D+01, 0.6286599010068D+01,
     :  0.6207190138027D-09, 0.5074049319576D+01, 0.4705732307012D+01,
     :  0.5989826532569D-09, 0.2231842216620D+01, 0.6256777527156D+01,
     :  0.5961360812618D-09, 0.1274975769045D+01, 0.6309374173736D+01,
     :  0.4874165471016D-09, 0.3642277426779D+01, 0.7755226100720D+00 /
      DATA ((E1(I,J,2),I=1,3),J= 11, 20) /
     :  0.4283834034360D-09, 0.5148765510106D+01, 0.1059381944224D+01,
     :  0.4652389287529D-09, 0.4715794792175D+01, 0.7860419393880D+01,
     :  0.3751707476401D-09, 0.6617207370325D+00, 0.5753384878334D+01,
     :  0.3559998806198D-09, 0.6155548875404D+01, 0.5884926831456D+01,
     :  0.3558447558857D-09, 0.2898827297664D+01, 0.6812766822558D+01,
     :  0.3211116927106D-09, 0.3625813502509D+01, 0.6681224869435D+01,
     :  0.2875609914672D-09, 0.4345435813134D+01, 0.2513230340178D+02,
     :  0.2843109704069D-09, 0.5862263940038D+01, 0.6127655567643D+01,
     :  0.2744676468427D-09, 0.3926419475089D+01, 0.6438496133249D+01,
     :  0.2481285237789D-09, 0.1351976572828D+01, 0.5486777812467D+01 /
      DATA ((E1(I,J,2),I=1,3),J= 21, 30) /
     :  0.2060338481033D-09, 0.2147556998591D+01, 0.7079373888424D+01,
     :  0.2015822358331D-09, 0.4408358972216D+01, 0.6290189305114D+01,
     :  0.2001195944195D-09, 0.5385829822531D+01, 0.6275962395778D+01,
     :  0.1953667642377D-09, 0.1304933746120D+01, 0.5507553240374D+01,
     :  0.1839744078713D-09, 0.6173567228835D+01, 0.1179062909082D+02,
     :  0.1643334294845D-09, 0.4635942997523D+01, 0.1150676975667D+02,
     :  0.1768051018652D-09, 0.5086283558874D+01, 0.7113454667900D-02,
     :  0.1674874205489D-09, 0.2243332137241D+01, 0.7058598460518D+01,
     :  0.1421445397609D-09, 0.6186899771515D+01, 0.7962980379786D+00,
     :  0.1255163958267D-09, 0.5730238465658D+01, 0.4694002934110D+01 /
      DATA ((E1(I,J,2),I=1,3),J= 31, 40) /
     :  0.1013945281961D-09, 0.1726055228402D+01, 0.3738761453707D+01,
     :  0.1047294335852D-09, 0.2658801228129D+01, 0.6282095334605D+01,
     :  0.1047103879392D-09, 0.8481047835035D+00, 0.6284056366286D+01,
     :  0.9530343962826D-10, 0.3079267149859D+01, 0.6069776770667D+01,
     :  0.9604637611690D-10, 0.3258679792918D+00, 0.4136910472696D+01,
     :  0.9153518537177D-10, 0.4398599886584D+00, 0.6496374930224D+01,
     :  0.8562458214922D-10, 0.4772686794145D+01, 0.1194447056968D+01,
     :  0.8232525360654D-10, 0.5966220721679D+01, 0.1589072916335D+01,
     :  0.6150223411438D-10, 0.1780985591923D+01, 0.8827390247185D+01,
     :  0.6272087858000D-10, 0.3184305429012D+01, 0.8429241228195D+01 /
      DATA ((E1(I,J,2),I=1,3),J= 41, 50) /
     :  0.5540476311040D-10, 0.3801260595433D+01, 0.4933208510675D+01,
     :  0.7331901699361D-10, 0.5205948591865D+01, 0.4535059491685D+01,
     :  0.6018528702791D-10, 0.4770139083623D+01, 0.1255903824622D+02,
     :  0.5150530724804D-10, 0.3574796899585D+01, 0.1176985366291D+02,
     :  0.6471933741811D-10, 0.2679787266521D+01, 0.5088628793478D+01,
     :  0.5317460644174D-10, 0.9528763345494D+00, 0.3154687086868D+01,
     :  0.4832187748783D-10, 0.5329322498232D+01, 0.6040347114260D+01,
     :  0.4716763555110D-10, 0.2395235316466D+01, 0.5331357529664D+01,
     :  0.4871509139861D-10, 0.3056663648823D+01, 0.1256967486051D+02,
     :  0.4598417696768D-10, 0.4452762609019D+01, 0.6525804586632D+01 /
      DATA ((E1(I,J,2),I=1,3),J= 51, 60) /
     :  0.5674189533175D-10, 0.9879680872193D+00, 0.5729506548653D+01,
     :  0.4073560328195D-10, 0.5939127696986D+01, 0.7632943190217D+01,
     :  0.5040994945359D-10, 0.4549875824510D+01, 0.8031092209206D+01,
     :  0.5078185134679D-10, 0.7346659893982D+00, 0.7477522907414D+01,
     :  0.3769343537061D-10, 0.1071317188367D+01, 0.7234794171227D+01,
     :  0.4980331365299D-10, 0.2500345341784D+01, 0.6836645152238D+01,
     :  0.3458236594757D-10, 0.3825159450711D+01, 0.1097707878456D+02,
     :  0.3578859493602D-10, 0.5299664791549D+01, 0.4164311961999D+01,
     :  0.3370504646419D-10, 0.5002316301593D+01, 0.1137170464392D+02,
     :  0.3299873338428D-10, 0.2526123275282D+01, 0.3930209696940D+01 /
      DATA ((E1(I,J,2),I=1,3),J= 61, 70) /
     :  0.4304917318409D-10, 0.3368078557132D+01, 0.1592596075957D+01,
     :  0.3402418753455D-10, 0.8385495425800D+00, 0.3128388763578D+01,
     :  0.2778460572146D-10, 0.3669905203240D+01, 0.7342457794669D+01,
     :  0.2782710128902D-10, 0.2691664812170D+00, 0.1748016358760D+01,
     :  0.2711725179646D-10, 0.4707487217718D+01, 0.5296909721118D+00,
     :  0.2981760946340D-10, 0.3190260867816D+00, 0.5368044267797D+00,
     :  0.2811672977772D-10, 0.3196532315372D+01, 0.7084896783808D+01,
     :  0.2863454474467D-10, 0.2263240324780D+00, 0.5223693906222D+01,
     :  0.3333464634051D-10, 0.3498451685065D+01, 0.8018209333619D+00,
     :  0.3312991747609D-10, 0.5839154477412D+01, 0.1554202828031D+00 /
      DATA ((E1(I,J,2),I=1,3),J= 71,NE1Y) /
     :  0.2813255564006D-10, 0.8268044346621D+00, 0.5225775174439D+00,
     :  0.2665098083966D-10, 0.3934021725360D+01, 0.5216580451554D+01,
     :  0.2349795705216D-10, 0.5197620913779D+01, 0.2146165377750D+01,
     :  0.2330352293961D-10, 0.2984999231807D+01, 0.1726015463500D+02,
     :  0.2728001683419D-10, 0.6521679638544D+00, 0.8635942003952D+01,
     :  0.2484061007669D-10, 0.3468955561097D+01, 0.5230807360890D+01,
     :  0.2646328768427D-10, 0.1013724533516D+01, 0.2629832328990D-01,
     :  0.2518630264831D-10, 0.6108081057122D+01, 0.5481254917084D+01,
     :  0.2421901455384D-10, 0.1651097776260D+01, 0.1349867339771D+01,
     :  0.6348533267831D-11, 0.3220226560321D+01, 0.8433466158131D+02 /

*  Sun-to-Earth, T^2, Y
      DATA ((E2(I,J,2),I=1,3),J=  1,NE2Y) /
     :  0.5063375872532D-10, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.2173815785980D-10, 0.2827805833053D+01, 0.1256615170089D+02,
     :  0.1010231999920D-10, 0.4634612377133D+01, 0.6283075850446D+01,
     :  0.9259745317636D-12, 0.2620612076189D+01, 0.1884922755134D+02,
     :  0.1022202095812D-12, 0.3809562326066D+01, 0.8399684731857D+02 /

*  Sun-to-Earth, T^0, Z
      DATA ((E0(I,J,3),I=1,3),J=  1, 10) /
     :  0.2796207639075D-05, 0.3198701560209D+01, 0.8433466158131D+02,
     :  0.1016042198142D-05, 0.5422360395913D+01, 0.5507553240374D+01,
     :  0.8044305033647D-06, 0.3880222866652D+01, 0.5223693906222D+01,
     :  0.4385347909274D-06, 0.3704369937468D+01, 0.2352866153506D+01,
     :  0.3186156414906D-06, 0.3999639363235D+01, 0.1577343543434D+01,
     :  0.2272412285792D-06, 0.3984738315952D+01, 0.1047747311755D+01,
     :  0.1645620103007D-06, 0.3565412516841D+01, 0.5856477690889D+01,
     :  0.1815836921166D-06, 0.4984507059020D+01, 0.6283075850446D+01,
     :  0.1447461676364D-06, 0.3702753570108D+01, 0.9437762937313D+01,
     :  0.1430760876382D-06, 0.3409658712357D+01, 0.1021328554739D+02 /
      DATA ((E0(I,J,3),I=1,3),J= 11, 20) /
     :  0.1120445753226D-06, 0.4829561570246D+01, 0.1414349524433D+02,
     :  0.1090232840797D-06, 0.2080729178066D+01, 0.6812766822558D+01,
     :  0.9715727346551D-07, 0.3476295881948D+01, 0.4694002934110D+01,
     :  0.1036267136217D-06, 0.4056639536648D+01, 0.7109288135493D+02,
     :  0.8752665271340D-07, 0.4448159519911D+01, 0.5753384878334D+01,
     :  0.8331864956004D-07, 0.4991704044208D+01, 0.7084896783808D+01,
     :  0.6901658670245D-07, 0.4325358994219D+01, 0.6275962395778D+01,
     :  0.9144536848998D-07, 0.1141826375363D+01, 0.6620890113188D+01,
     :  0.7205085037435D-07, 0.3624344170143D+01, 0.5296909721118D+00,
     :  0.7697874654176D-07, 0.5554257458998D+01, 0.1676215758509D+03 /
      DATA ((E0(I,J,3),I=1,3),J= 21, 30) /
     :  0.5197545738384D-07, 0.6251760961735D+01, 0.1807370494127D+02,
     :  0.5031345378608D-07, 0.2497341091913D+01, 0.4705732307012D+01,
     :  0.4527110205840D-07, 0.2335079920992D+01, 0.6309374173736D+01,
     :  0.4753355798089D-07, 0.7094148987474D+00, 0.5884926831456D+01,
     :  0.4296951977516D-07, 0.1101916352091D+01, 0.6681224869435D+01,
     :  0.3855341568387D-07, 0.1825495405486D+01, 0.5486777812467D+01,
     :  0.5253930970990D-07, 0.4424740687208D+01, 0.7860419393880D+01,
     :  0.4024630496471D-07, 0.5120498157053D+01, 0.1336797263425D+02,
     :  0.4061069791453D-07, 0.6029771435451D+01, 0.3930209696940D+01,
     :  0.3797883804205D-07, 0.4435193600836D+00, 0.3154687086868D+01 /
      DATA ((E0(I,J,3),I=1,3),J= 31, 40) /
     :  0.2933033225587D-07, 0.5124157356507D+01, 0.1059381944224D+01,
     :  0.3503000930426D-07, 0.5421830162065D+01, 0.6069776770667D+01,
     :  0.3670096214050D-07, 0.4582101667297D+01, 0.1219403291462D+02,
     :  0.2905609437008D-07, 0.1926566420072D+01, 0.1097707878456D+02,
     :  0.2466827821713D-07, 0.6090174539834D+00, 0.6496374930224D+01,
     :  0.2691647295332D-07, 0.1393432595077D+01, 0.2200391463820D+02,
     :  0.2150554667946D-07, 0.4308671715951D+01, 0.5643178611111D+01,
     :  0.2237481922680D-07, 0.8133968269414D+00, 0.8635942003952D+01,
     :  0.1817741038157D-07, 0.3755205127454D+01, 0.3340612434717D+01,
     :  0.2227820762132D-07, 0.2759558596664D+01, 0.1203646072878D+02 /
      DATA ((E0(I,J,3),I=1,3),J= 41, 50) /
     :  0.1944713772307D-07, 0.5699645869121D+01, 0.1179062909082D+02,
     :  0.1527340520662D-07, 0.1986749091746D+01, 0.3981490189893D+00,
     :  0.1577282574914D-07, 0.3205017217983D+01, 0.5088628793478D+01,
     :  0.1424738825424D-07, 0.6256747903666D+01, 0.2544314396739D+01,
     :  0.1616563121701D-07, 0.2601671259394D+00, 0.1729818233119D+02,
     :  0.1401210391692D-07, 0.4686939173506D+01, 0.7058598460518D+01,
     :  0.1488726974214D-07, 0.2815862451372D+01, 0.2593412433514D+02,
     :  0.1692626442388D-07, 0.4956894109797D+01, 0.1564752902480D+03,
     :  0.1123571582910D-07, 0.2381192697696D+01, 0.3738761453707D+01,
     :  0.9903308606317D-08, 0.4294851657684D+01, 0.9225539266174D+01 /
      DATA ((E0(I,J,3),I=1,3),J= 51, 60) /
     :  0.9174533187191D-08, 0.3075171510642D+01, 0.4164311961999D+01,
     :  0.8645985631457D-08, 0.5477534821633D+00, 0.8429241228195D+01,
     : -0.1085876492688D-07, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.9264309077815D-08, 0.5968571670097D+01, 0.7079373888424D+01,
     :  0.8243116984954D-08, 0.1489098777643D+01, 0.1044738781244D+02,
     :  0.8268102113708D-08, 0.3512977691983D+01, 0.1150676975667D+02,
     :  0.9043613988227D-08, 0.1290704408221D+00, 0.1101510648075D+02,
     :  0.7432912038789D-08, 0.1991086893337D+01, 0.2608790314060D+02,
     :  0.8586233727285D-08, 0.4238357924414D+01, 0.2986433403208D+02,
     :  0.7612230060131D-08, 0.2911090150166D+01, 0.4732030630302D+01 /
      DATA ((E0(I,J,3),I=1,3),J= 61, 70) /
     :  0.7097787751408D-08, 0.1908938392390D+01, 0.8031092209206D+01,
     :  0.7640237040175D-08, 0.6129219000168D+00, 0.7962980379786D+00,
     :  0.7070445688081D-08, 0.1380417036651D+01, 0.2146165377750D+01,
     :  0.7690770957702D-08, 0.1680504249084D+01, 0.2122839202813D+02,
     :  0.8051292542594D-08, 0.5127423484511D+01, 0.2942463415728D+01,
     :  0.5902709104515D-08, 0.2020274190917D+01, 0.7755226100720D+00,
     :  0.5134567496462D-08, 0.2606778676418D+01, 0.1256615170089D+02,
     :  0.5525802046102D-08, 0.1613011769663D+01, 0.8018209333619D+00,
     :  0.5880724784221D-08, 0.4604483417236D+01, 0.4690479774488D+01,
     :  0.5211699081370D-08, 0.5718964114193D+01, 0.8827390247185D+01 /
      DATA ((E0(I,J,3),I=1,3),J= 71, 80) /
     :  0.4891849573562D-08, 0.3689658932196D+01, 0.2132990797783D+00,
     :  0.5150246069997D-08, 0.4099769855122D+01, 0.6480980550449D+02,
     :  0.5102434319633D-08, 0.5660834602509D+01, 0.3379454372902D+02,
     :  0.5083405254252D-08, 0.9842221218974D+00, 0.4136910472696D+01,
     :  0.4206562585682D-08, 0.1341363634163D+00, 0.3128388763578D+01,
     :  0.4663249683579D-08, 0.8130132735866D+00, 0.5216580451554D+01,
     :  0.4099474416530D-08, 0.5791497770644D+01, 0.4265981595566D+00,
     :  0.4628251220767D-08, 0.1249802769331D+01, 0.1572083878776D+02,
     :  0.5024068728142D-08, 0.4795684802743D+01, 0.6290189305114D+01,
     :  0.5120234327758D-08, 0.3810420387208D+01, 0.5230807360890D+01 /
      DATA ((E0(I,J,3),I=1,3),J= 81, 90) /
     :  0.5524029815280D-08, 0.1029264714351D+01, 0.2397622045175D+03,
     :  0.4757415718860D-08, 0.3528044781779D+01, 0.1649636139783D+02,
     :  0.3915786131127D-08, 0.5593889282646D+01, 0.1589072916335D+01,
     :  0.4869053149991D-08, 0.3299636454433D+01, 0.7632943190217D+01,
     :  0.3649365703729D-08, 0.1286049002584D+01, 0.6206810014183D+01,
     :  0.3992493949002D-08, 0.3100307589464D+01, 0.2515860172507D+02,
     :  0.3320247477418D-08, 0.6212683940807D+01, 0.1216800268190D+02,
     :  0.3287123739696D-08, 0.4699118445928D+01, 0.7234794171227D+01,
     :  0.3472776811103D-08, 0.2630507142004D+01, 0.7342457794669D+01,
     :  0.3423253294767D-08, 0.2946432844305D+01, 0.9623688285163D+01 /
      DATA ((E0(I,J,3),I=1,3),J= 91,100) /
     :  0.3896173898244D-08, 0.1224834179264D+01, 0.6438496133249D+01,
     :  0.3388455337924D-08, 0.1543807616351D+01, 0.1494531617769D+02,
     :  0.3062704716523D-08, 0.1191777572310D+01, 0.8662240327241D+01,
     :  0.3270075600400D-08, 0.5483498767737D+01, 0.1194447056968D+01,
     :  0.3101209215259D-08, 0.8000833804348D+00, 0.3772475342596D+02,
     :  0.2780883347311D-08, 0.4077980721888D+00, 0.5863591145557D+01,
     :  0.2903605931824D-08, 0.2617490302147D+01, 0.1965104848470D+02,
     :  0.2682014743119D-08, 0.2634703158290D+01, 0.7238675589263D+01,
     :  0.2534360108492D-08, 0.6102446114873D+01, 0.6836645152238D+01,
     :  0.2392564882509D-08, 0.3681820208691D+01, 0.5849364236221D+01 /
      DATA ((E0(I,J,3),I=1,3),J=101,110) /
     :  0.2656667254856D-08, 0.6216045388886D+01, 0.6133512519065D+01,
     :  0.2331242096773D-08, 0.5864949777744D+01, 0.4535059491685D+01,
     :  0.2287898363668D-08, 0.4566628532802D+01, 0.7477522907414D+01,
     :  0.2336944521306D-08, 0.2442722126930D+01, 0.1137170464392D+02,
     :  0.3156632236269D-08, 0.1626628050682D+01, 0.2509084901204D+03,
     :  0.2982612402766D-08, 0.2803604512609D+01, 0.1748016358760D+01,
     :  0.2774031674807D-08, 0.4654002897158D+01, 0.8223916695780D+02,
     :  0.2295236548638D-08, 0.4326518333253D+01, 0.3378142627421D+00,
     :  0.2190714699873D-08, 0.4519614578328D+01, 0.2908881142201D+02,
     :  0.2191495845045D-08, 0.3012626912549D+01, 0.1673046366289D+02 /
      DATA ((E0(I,J,3),I=1,3),J=111,120) /
     :  0.2492901628386D-08, 0.1290101424052D+00, 0.1543797956245D+03,
     :  0.1993778064319D-08, 0.3864046799414D+01, 0.1778984560711D+02,
     :  0.1898146479022D-08, 0.5053777235891D+01, 0.2042657109477D+02,
     :  0.1918280127634D-08, 0.2222470192548D+01, 0.4165496312290D+02,
     :  0.1916351061607D-08, 0.8719067257774D+00, 0.7737595720538D+02,
     :  0.1834720181466D-08, 0.4031491098040D+01, 0.2358125818164D+02,
     :  0.1249201523806D-08, 0.5938379466835D+01, 0.3301902111895D+02,
     :  0.1477304050539D-08, 0.6544722606797D+00, 0.9548094718417D+02,
     :  0.1264316431249D-08, 0.2059072853236D+01, 0.8399684731857D+02,
     :  0.1203526495039D-08, 0.3644813532605D+01, 0.4558517281984D+02 /
      DATA ((E0(I,J,3),I=1,3),J=121,130) /
     :  0.9221681059831D-09, 0.3241815055602D+01, 0.7805158573086D+02,
     :  0.7849278367646D-09, 0.5043812342457D+01, 0.5217580628120D+02,
     :  0.7983392077387D-09, 0.5000024502753D+01, 0.1501922143975D+03,
     :  0.7925395431654D-09, 0.1398734871821D-01, 0.9061773743175D+02,
     :  0.7640473285886D-09, 0.5067111723130D+01, 0.4951538251678D+02,
     :  0.5398937754482D-09, 0.5597382200075D+01, 0.1613385000004D+03,
     :  0.5626247550193D-09, 0.2601338209422D+01, 0.7318837597844D+02,
     :  0.5525197197855D-09, 0.5814832109256D+01, 0.1432335100216D+03,
     :  0.5407629837898D-09, 0.3384820609076D+01, 0.3230491187871D+03,
     :  0.3856739119801D-09, 0.1072391840473D+01, 0.2334791286671D+03 /
      DATA ((E0(I,J,3),I=1,3),J=131,NE0Z) /
     :  0.3856425239987D-09, 0.2369540393327D+01, 0.1739046517013D+03,
     :  0.4350867755983D-09, 0.5255575751082D+01, 0.1620484330494D+03,
     :  0.3844113924996D-09, 0.5482356246182D+01, 0.9757644180768D+02,
     :  0.2854869155431D-09, 0.9573634763143D+00, 0.1697170704744D+03,
     :  0.1719227671416D-09, 0.1887203025202D+01, 0.2265204242912D+03,
     :  0.1527846879755D-09, 0.3982183931157D+01, 0.3341954043900D+03,
     :  0.1128229264847D-09, 0.2787457156298D+01, 0.3119028331842D+03 /

*  Sun-to-Earth, T^1, Z
      DATA ((E1(I,J,3),I=1,3),J=  1, 10) /
     :  0.2278290449966D-05, 0.3413716033863D+01, 0.6283075850446D+01,
     :  0.5429458209830D-07, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.1903240492525D-07, 0.3370592358297D+01, 0.1256615170089D+02,
     :  0.2385409276743D-09, 0.3327914718416D+01, 0.1884922755134D+02,
     :  0.8676928342573D-10, 0.1824006811264D+01, 0.5223693906222D+01,
     :  0.7765442593544D-10, 0.3888564279247D+01, 0.5507553240374D+01,
     :  0.7066158332715D-10, 0.5194267231944D+01, 0.2352866153506D+01,
     :  0.7092175288657D-10, 0.2333246960021D+01, 0.8399684731857D+02,
     :  0.5357582213535D-10, 0.2224031176619D+01, 0.5296909721118D+00,
     :  0.3828035865021D-10, 0.2156710933584D+01, 0.6279552690824D+01 /
      DATA ((E1(I,J,3),I=1,3),J= 11,NE1Z) /
     :  0.3824857220427D-10, 0.1529755219915D+01, 0.6286599010068D+01,
     :  0.3286995181628D-10, 0.4879512900483D+01, 0.1021328554739D+02 /

*  Sun-to-Earth, T^2, Z
      DATA ((E2(I,J,3),I=1,3),J=  1,NE2Z) /
     :  0.9722666114891D-10, 0.5152219582658D+01, 0.6283075850446D+01,
     : -0.3494819171909D-11, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.6713034376076D-12, 0.6440188750495D+00, 0.1256615170089D+02 /

*  SSB-to-Sun, T^0, X
      DATA ((S0(I,J,1),I=1,3),J=  1, 10) /
     :  0.4956757536410D-02, 0.3741073751789D+01, 0.5296909721118D+00,
     :  0.2718490072522D-02, 0.4016011511425D+01, 0.2132990797783D+00,
     :  0.1546493974344D-02, 0.2170528330642D+01, 0.3813291813120D-01,
     :  0.8366855276341D-03, 0.2339614075294D+01, 0.7478166569050D-01,
     :  0.2936777942117D-03, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.1201317439469D-03, 0.4090736353305D+01, 0.1059381944224D+01,
     :  0.7578550887230D-04, 0.3241518088140D+01, 0.4265981595566D+00,
     :  0.1941787367773D-04, 0.1012202064330D+01, 0.2061856251104D+00,
     :  0.1889227765991D-04, 0.3892520416440D+01, 0.2204125344462D+00,
     :  0.1937896968613D-04, 0.4797779441161D+01, 0.1495633313810D+00 /
      DATA ((S0(I,J,1),I=1,3),J= 11, 20) /
     :  0.1434506110873D-04, 0.3868960697933D+01, 0.5225775174439D+00,
     :  0.1406659911580D-04, 0.4759766557397D+00, 0.5368044267797D+00,
     :  0.1179022300202D-04, 0.7774961520598D+00, 0.7626583626240D-01,
     :  0.8085864460959D-05, 0.3254654471465D+01, 0.3664874755930D-01,
     :  0.7622752967615D-05, 0.4227633103489D+01, 0.3961708870310D-01,
     :  0.6209171139066D-05, 0.2791828325711D+00, 0.7329749511860D-01,
     :  0.4366435633970D-05, 0.4440454875925D+01, 0.1589072916335D+01,
     :  0.3792124889348D-05, 0.5156393842356D+01, 0.7113454667900D-02,
     :  0.3154548963402D-05, 0.6157005730093D+01, 0.4194847048887D+00,
     :  0.3088359882942D-05, 0.2494567553163D+01, 0.6398972393349D+00 /
      DATA ((S0(I,J,1),I=1,3),J= 21, 30) /
     :  0.2788440902136D-05, 0.4934318747989D+01, 0.1102062672231D+00,
     :  0.3039928456376D-05, 0.4895077702640D+01, 0.6283075850446D+01,
     :  0.2272258457679D-05, 0.5278394064764D+01, 0.1030928125552D+00,
     :  0.2162007057957D-05, 0.5802978019099D+01, 0.3163918923335D+00,
     :  0.1767632855737D-05, 0.3415346595193D-01, 0.1021328554739D+02,
     :  0.1349413459362D-05, 0.2001643230755D+01, 0.1484170571900D-02,
     :  0.1170141900476D-05, 0.2424750491620D+01, 0.6327837846670D+00,
     :  0.1054355266820D-05, 0.3123311487576D+01, 0.4337116142245D+00,
     :  0.9800822461610D-06, 0.3026258088130D+01, 0.1052268489556D+01,
     :  0.1091203749931D-05, 0.3157811670347D+01, 0.1162474756779D+01 /
      DATA ((S0(I,J,1),I=1,3),J= 31, 40) /
     :  0.6960236715913D-06, 0.8219570542313D+00, 0.1066495398892D+01,
     :  0.5689257296909D-06, 0.1323052375236D+01, 0.9491756770005D+00,
     :  0.6613172135802D-06, 0.2765348881598D+00, 0.8460828644453D+00,
     :  0.6277702517571D-06, 0.5794064466382D+01, 0.1480791608091D+00,
     :  0.6304884066699D-06, 0.7323555380787D+00, 0.2243449970715D+00,
     :  0.4897850467382D-06, 0.3062464235399D+01, 0.3340612434717D+01,
     :  0.3759148598786D-06, 0.4588290469664D+01, 0.3516457698740D-01,
     :  0.3110520548195D-06, 0.1374299536572D+01, 0.6373574839730D-01,
     :  0.3064708359780D-06, 0.4222267485047D+01, 0.1104591729320D-01,
     :  0.2856347168241D-06, 0.3714202944973D+01, 0.1510475019529D+00 /
      DATA ((S0(I,J,1),I=1,3),J= 41, 50) /
     :  0.2840945514288D-06, 0.2847972875882D+01, 0.4110125927500D-01,
     :  0.2378951599405D-06, 0.3762072563388D+01, 0.2275259891141D+00,
     :  0.2714229481417D-06, 0.1036049980031D+01, 0.2535050500000D-01,
     :  0.2323551717307D-06, 0.4682388599076D+00, 0.8582758298370D-01,
     :  0.1881790512219D-06, 0.4790565425418D+01, 0.2118763888447D+01,
     :  0.2261353968371D-06, 0.1669144912212D+01, 0.7181332454670D-01,
     :  0.2214546389848D-06, 0.3937717281614D+01, 0.2968341143800D-02,
     :  0.2184915594933D-06, 0.1129169845099D+00, 0.7775000683430D-01,
     :  0.2000164937936D-06, 0.4030009638488D+01, 0.2093666171530D+00,
     :  0.1966105136719D-06, 0.8745955786834D+00, 0.2172315424036D+00 /
      DATA ((S0(I,J,1),I=1,3),J= 51, 60) /
     :  0.1904742332624D-06, 0.5919743598964D+01, 0.2022531624851D+00,
     :  0.1657399705031D-06, 0.2549141484884D+01, 0.7358765972222D+00,
     :  0.1574070533987D-06, 0.5277533020230D+01, 0.7429900518901D+00,
     :  0.1832261651039D-06, 0.3064688127777D+01, 0.3235053470014D+00,
     :  0.1733615346569D-06, 0.3011432799094D+01, 0.1385174140878D+00,
     :  0.1549124014496D-06, 0.4005569132359D+01, 0.5154640627760D+00,
     :  0.1637044713838D-06, 0.1831375966632D+01, 0.8531963191132D+00,
     :  0.1123420082383D-06, 0.1180270407578D+01, 0.1990721704425D+00,
     :  0.1083754165740D-06, 0.3414101320863D+00, 0.5439178814476D+00,
     :  0.1156638012655D-06, 0.6130479452594D+00, 0.5257585094865D+00 /
      DATA ((S0(I,J,1),I=1,3),J= 61, 70) /
     :  0.1142548785134D-06, 0.3724761948846D+01, 0.5336234347371D+00,
     :  0.7921463895965D-07, 0.2435425589361D+01, 0.1478866649112D+01,
     :  0.7428600285231D-07, 0.3542144398753D+01, 0.2164800718209D+00,
     :  0.8323211246747D-07, 0.3525058072354D+01, 0.1692165728891D+01,
     :  0.7257595116312D-07, 0.1364299431982D+01, 0.2101180877357D+00,
     :  0.7111185833236D-07, 0.2460478875808D+01, 0.4155522422634D+00,
     :  0.6868090383716D-07, 0.4397327670704D+01, 0.1173197218910D+00,
     :  0.7226419974175D-07, 0.4042647308905D+01, 0.1265567569334D+01,
     :  0.6955642383177D-07, 0.2865047906085D+01, 0.9562891316684D+00,
     :  0.7492139296331D-07, 0.5014278994215D+01, 0.1422690933580D-01 /
      DATA ((S0(I,J,1),I=1,3),J= 71, 80) /
     :  0.6598363128857D-07, 0.2376730020492D+01, 0.6470106940028D+00,
     :  0.7381147293385D-07, 0.3272990384244D+01, 0.1581959461667D+01,
     :  0.6402909624032D-07, 0.5302290955138D+01, 0.9597935788730D-01,
     :  0.6237454263857D-07, 0.5444144425332D+01, 0.7084920306520D-01,
     :  0.5241198544016D-07, 0.4215359579205D+01, 0.5265099800692D+00,
     :  0.5144463853918D-07, 0.1218916689916D+00, 0.5328719641544D+00,
     :  0.5868164772299D-07, 0.2369402002213D+01, 0.7871412831580D-01,
     :  0.6233195669151D-07, 0.1254922242403D+01, 0.2608790314060D+02,
     :  0.6068463791422D-07, 0.5679713760431D+01, 0.1114304132498D+00,
     :  0.4359361135065D-07, 0.6097219641646D+00, 0.1375773836557D+01 /
      DATA ((S0(I,J,1),I=1,3),J= 81, 90) /
     :  0.4686510366826D-07, 0.4786231041431D+01, 0.1143987543936D+00,
     :  0.3758977287225D-07, 0.1167368068139D+01, 0.1596186371003D+01,
     :  0.4282051974778D-07, 0.1519471064319D+01, 0.2770348281756D+00,
     :  0.5153765386113D-07, 0.1860532322984D+01, 0.2228608264996D+00,
     :  0.4575129387188D-07, 0.7632857887158D+00, 0.1465949902372D+00,
     :  0.3326844933286D-07, 0.1298219485285D+01, 0.5070101000000D-01,
     :  0.3748617450984D-07, 0.1046510321062D+01, 0.4903339079539D+00,
     :  0.2816756661499D-07, 0.3434522346190D+01, 0.2991266627620D+00,
     :  0.3412750405039D-07, 0.2523766270318D+01, 0.3518164938661D+00,
     :  0.2655796761776D-07, 0.2904422260194D+01, 0.6256703299991D+00 /
      DATA ((S0(I,J,1),I=1,3),J= 91,100) /
     :  0.2963597929458D-07, 0.5923900431149D+00, 0.1099462426779D+00,
     :  0.2539523734781D-07, 0.4851947722567D+01, 0.1256615170089D+02,
     :  0.2283087914139D-07, 0.3400498595496D+01, 0.6681224869435D+01,
     :  0.2321309799331D-07, 0.5789099148673D+01, 0.3368040641550D-01,
     :  0.2549657649750D-07, 0.3991856479792D-01, 0.1169588211447D+01,
     :  0.2290462303977D-07, 0.2788567577052D+01, 0.1045155034888D+01,
     :  0.1945398522914D-07, 0.3290896998176D+01, 0.1155361302111D+01,
     :  0.1849171512638D-07, 0.2698060129367D+01, 0.4452511715700D-02,
     :  0.1647199834254D-07, 0.3016735644085D+01, 0.4408250688924D+00,
     :  0.1529530765273D-07, 0.5573043116178D+01, 0.6521991896920D-01 /
      DATA ((S0(I,J,1),I=1,3),J=101,110) /
     :  0.1433199339978D-07, 0.1481192356147D+01, 0.9420622223326D+00,
     :  0.1729134193602D-07, 0.1422817538933D+01, 0.2108507877249D+00,
     :  0.1716463931346D-07, 0.3469468901855D+01, 0.2157473718317D+00,
     :  0.1391206061378D-07, 0.6122436220547D+01, 0.4123712502208D+00,
     :  0.1404746661924D-07, 0.1647765641936D+01, 0.4258542984690D-01,
     :  0.1410452399455D-07, 0.5989729161964D+01, 0.2258291676434D+00,
     :  0.1089828772168D-07, 0.2833705509371D+01, 0.4226656969313D+00,
     :  0.1047374564948D-07, 0.5090690007331D+00, 0.3092784376656D+00,
     :  0.1358279126532D-07, 0.5128990262836D+01, 0.7923417740620D-01,
     :  0.1020456476148D-07, 0.9632772880808D+00, 0.1456308687557D+00 /
      DATA ((S0(I,J,1),I=1,3),J=111,120) /
     :  0.1033428735328D-07, 0.3223779318418D+01, 0.1795258541446D+01,
     :  0.1412435841540D-07, 0.2410271572721D+01, 0.1525316725248D+00,
     :  0.9722759371574D-08, 0.2333531395690D+01, 0.8434341241180D-01,
     :  0.9657334084704D-08, 0.6199270974168D+01, 0.1272681024002D+01,
     :  0.1083641148690D-07, 0.2864222292929D+01, 0.7032915397480D-01,
     :  0.1067318403838D-07, 0.5833458866568D+00, 0.2123349582968D+00,
     :  0.1062366201976D-07, 0.4307753989494D+01, 0.2142632012598D+00,
     :  0.1236364149266D-07, 0.2873917870593D+01, 0.1847279083684D+00,
     :  0.1092759489593D-07, 0.2959887266733D+01, 0.1370332435159D+00,
     :  0.8912069362899D-08, 0.5141213702562D+01, 0.2648454860559D+01 /
      DATA ((S0(I,J,1),I=1,3),J=121,130) /
     :  0.9656467707970D-08, 0.4532182462323D+01, 0.4376440768498D+00,
     :  0.8098386150135D-08, 0.2268906338379D+01, 0.2880807454688D+00,
     :  0.7857714675000D-08, 0.4055544260745D+01, 0.2037373330570D+00,
     :  0.7288455940646D-08, 0.5357901655142D+01, 0.1129145838217D+00,
     :  0.9450595950552D-08, 0.4264926963939D+01, 0.5272426800584D+00,
     :  0.9381718247537D-08, 0.7489366976576D-01, 0.5321392641652D+00,
     :  0.7079052646038D-08, 0.1923311052874D+01, 0.6288513220417D+00,
     :  0.9259004415344D-08, 0.2970256853438D+01, 0.1606092486742D+00,
     :  0.8259801499742D-08, 0.3327056314697D+01, 0.8389694097774D+00,
     :  0.6476334355779D-08, 0.2954925505727D+01, 0.2008557621224D+01 /
      DATA ((S0(I,J,1),I=1,3),J=131,140) /
     :  0.5984021492007D-08, 0.9138753105829D+00, 0.2042657109477D+02,
     :  0.5989546863181D-08, 0.3244464082031D+01, 0.2111650433779D+01,
     :  0.6233108606023D-08, 0.4995232638403D+00, 0.4305306221819D+00,
     :  0.6877299149965D-08, 0.2834987233449D+01, 0.9561746721300D-02,
     :  0.8311234227190D-08, 0.2202951835758D+01, 0.3801276407308D+00,
     :  0.6599472832414D-08, 0.4478581462618D+01, 0.1063314406849D+01,
     :  0.6160491096549D-08, 0.5145858696411D+01, 0.1368660381889D+01,
     :  0.6164772043891D-08, 0.3762976697911D+00, 0.4234171675140D+00,
     :  0.6363248684450D-08, 0.3162246718685D+01, 0.1253008786510D-01,
     :  0.6448587520999D-08, 0.3442693302119D+01, 0.5287268506303D+00 /
      DATA ((S0(I,J,1),I=1,3),J=141,150) /
     :  0.6431662283977D-08, 0.8977549136606D+00, 0.5306550935933D+00,
     :  0.6351223158474D-08, 0.4306447410369D+01, 0.5217580628120D+02,
     :  0.5476721393451D-08, 0.3888529177855D+01, 0.2221856701002D+01,
     :  0.5341772572619D-08, 0.2655560662512D+01, 0.7466759693650D-01,
     :  0.5337055758302D-08, 0.5164990735946D+01, 0.7489573444450D-01,
     :  0.5373120816787D-08, 0.6041214553456D+01, 0.1274714967946D+00,
     :  0.5392351705426D-08, 0.9177763485932D+00, 0.1055449481598D+01,
     :  0.6688495850205D-08, 0.3089608126937D+01, 0.2213766559277D+00,
     :  0.5072003660362D-08, 0.4311316541553D+01, 0.2132517061319D+00,
     :  0.5070726650455D-08, 0.5790675464444D+00, 0.2133464534247D+00 /
      DATA ((S0(I,J,1),I=1,3),J=151,160) /
     :  0.5658012950032D-08, 0.2703945510675D+01, 0.7287631425543D+00,
     :  0.4835509924854D-08, 0.2975422976065D+01, 0.7160067364790D-01,
     :  0.6479821978012D-08, 0.1324168733114D+01, 0.2209183458640D-01,
     :  0.6230636494980D-08, 0.2860103632836D+01, 0.3306188016693D+00,
     :  0.4649239516213D-08, 0.4832259763403D+01, 0.7796265773310D-01,
     :  0.6487325792700D-08, 0.2726165825042D+01, 0.3884652414254D+00,
     :  0.4682823682770D-08, 0.6966602455408D+00, 0.1073608853559D+01,
     :  0.5704230804976D-08, 0.5669634104606D+01, 0.8731175355560D-01,
     :  0.6125413585489D-08, 0.1513386538915D+01, 0.7605151500000D-01,
     :  0.6035825038187D-08, 0.1983509168227D+01, 0.9846002785331D+00 /
      DATA ((S0(I,J,1),I=1,3),J=161,170) /
     :  0.4331123462303D-08, 0.2782892992807D+01, 0.4297791515992D+00,
     :  0.4681107685143D-08, 0.5337232886836D+01, 0.2127790306879D+00,
     :  0.4669105829655D-08, 0.5837133792160D+01, 0.2138191288687D+00,
     :  0.5138823602365D-08, 0.3080560200507D+01, 0.7233337363710D-01,
     :  0.4615856664534D-08, 0.1661747897471D+01, 0.8603097737811D+00,
     :  0.4496916702197D-08, 0.2112508027068D+01, 0.7381754420900D-01,
     :  0.4278479042945D-08, 0.5716528462627D+01, 0.7574578717200D-01,
     :  0.3840525503932D-08, 0.6424172726492D+00, 0.3407705765729D+00,
     :  0.4866636509685D-08, 0.4919244697715D+01, 0.7722995774390D-01,
     :  0.3526100639296D-08, 0.2550821052734D+01, 0.6225157782540D-01 /
      DATA ((S0(I,J,1),I=1,3),J=171,180) /
     :  0.3939558488075D-08, 0.3939331491710D+01, 0.5268983110410D-01,
     :  0.4041268772576D-08, 0.2275337571218D+01, 0.3503323232942D+00,
     :  0.3948761842853D-08, 0.1999324200790D+01, 0.1451108196653D+00,
     :  0.3258394550029D-08, 0.9121001378200D+00, 0.5296435984654D+00,
     :  0.3257897048761D-08, 0.3428428660869D+01, 0.5297383457582D+00,
     :  0.3842559031298D-08, 0.6132927720035D+01, 0.9098186128426D+00,
     :  0.3109920095448D-08, 0.7693650193003D+00, 0.3932462625300D-02,
     :  0.3132237775119D-08, 0.3621293854908D+01, 0.2346394437820D+00,
     :  0.3942189421510D-08, 0.4841863659733D+01, 0.3180992042600D-02,
     :  0.3796972285340D-08, 0.1814174994268D+01, 0.1862120789403D+00 /
      DATA ((S0(I,J,1),I=1,3),J=181,190) /
     :  0.3995640233688D-08, 0.1386990406091D+01, 0.4549093064213D+00,
     :  0.2875013727414D-08, 0.9178318587177D+00, 0.1905464808669D+01,
     :  0.3073719932844D-08, 0.2688923811835D+01, 0.3628624111593D+00,
     :  0.2731016580075D-08, 0.1188259127584D+01, 0.2131850110243D+00,
     :  0.2729549896546D-08, 0.3702160634273D+01, 0.2134131485323D+00,
     :  0.3339372892449D-08, 0.7199163960331D+00, 0.2007689919132D+00,
     :  0.2898833764204D-08, 0.1916709364999D+01, 0.5291709230214D+00,
     :  0.2894536549362D-08, 0.2424043195547D+01, 0.5302110212022D+00,
     :  0.3096872473843D-08, 0.4445894977497D+01, 0.2976424921901D+00,
     :  0.2635672326810D-08, 0.3814366984117D+01, 0.1485980103780D+01 /
      DATA ((S0(I,J,1),I=1,3),J=191,200) /
     :  0.3649302697001D-08, 0.2924200596084D+01, 0.6044726378023D+00,
     :  0.3127954585895D-08, 0.1842251648327D+01, 0.1084620721060D+00,
     :  0.2616040173947D-08, 0.4155841921984D+01, 0.1258454114666D+01,
     :  0.2597395859860D-08, 0.1158045978874D+00, 0.2103781122809D+00,
     :  0.2593286172210D-08, 0.4771850408691D+01, 0.2162200472757D+00,
     :  0.2481823585747D-08, 0.4608842558889D+00, 0.1062562936266D+01,
     :  0.2742219550725D-08, 0.1538781127028D+01, 0.5651155736444D+00,
     :  0.3199558469610D-08, 0.3226647822878D+00, 0.7036329877322D+00,
     :  0.2666088542957D-08, 0.1967991731219D+00, 0.1400015846597D+00,
     :  0.2397067430580D-08, 0.3707036669873D+01, 0.2125476091956D+00 /
      DATA ((S0(I,J,1),I=1,3),J=201,210) /
     :  0.2376570772738D-08, 0.1182086628042D+01, 0.2140505503610D+00,
     :  0.2547228007887D-08, 0.4906256820629D+01, 0.1534957940063D+00,
     :  0.2265575594114D-08, 0.3414949866857D+01, 0.2235935264888D+00,
     :  0.2464381430585D-08, 0.4599122275378D+01, 0.2091065926078D+00,
     :  0.2433408527044D-08, 0.2830751145445D+00, 0.2174915669488D+00,
     :  0.2443605509076D-08, 0.4212046432538D+01, 0.1739420156204D+00,
     :  0.2319779262465D-08, 0.9881978408630D+00, 0.7530171478090D-01,
     :  0.2284622835465D-08, 0.5565347331588D+00, 0.7426161660010D-01,
     :  0.2467268750783D-08, 0.5655708150766D+00, 0.2526561439362D+00,
     :  0.2808513492782D-08, 0.1418405053408D+01, 0.5636314030725D+00 /
      DATA ((S0(I,J,1),I=1,3),J=211,NS0X) /
     :  0.2329528932532D-08, 0.4069557545675D+01, 0.1056200952181D+01,
     :  0.9698639532817D-09, 0.1074134313634D+01, 0.7826370942180D+02 /

*  SSB-to-Sun, T^1, X
      DATA ((S1(I,J,1),I=1,3),J=  1, 10) /
     : -0.1296310361520D-07, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.8975769009438D-08, 0.1128891609250D+01, 0.4265981595566D+00,
     :  0.7771113441307D-08, 0.2706039877077D+01, 0.2061856251104D+00,
     :  0.7538303866642D-08, 0.2191281289498D+01, 0.2204125344462D+00,
     :  0.6061384579336D-08, 0.3248167319958D+01, 0.1059381944224D+01,
     :  0.5726994235594D-08, 0.5569981398610D+01, 0.5225775174439D+00,
     :  0.5616492836424D-08, 0.5057386614909D+01, 0.5368044267797D+00,
     :  0.1010881584769D-08, 0.3473577116095D+01, 0.7113454667900D-02,
     :  0.7259606157626D-09, 0.3651858593665D+00, 0.6398972393349D+00,
     :  0.8755095026935D-09, 0.1662835408338D+01, 0.4194847048887D+00 /
      DATA ((S1(I,J,1),I=1,3),J= 11, 20) /
     :  0.5370491182812D-09, 0.1327673878077D+01, 0.4337116142245D+00,
     :  0.5743773887665D-09, 0.4250200846687D+01, 0.2132990797783D+00,
     :  0.4408103140300D-09, 0.3598752574277D+01, 0.1589072916335D+01,
     :  0.3101892374445D-09, 0.4887822983319D+01, 0.1052268489556D+01,
     :  0.3209453713578D-09, 0.9702272295114D+00, 0.5296909721118D+00,
     :  0.3017228286064D-09, 0.5484462275949D+01, 0.1066495398892D+01,
     :  0.3200700038601D-09, 0.2846613338643D+01, 0.1495633313810D+00,
     :  0.2137637279911D-09, 0.5692163292729D+00, 0.3163918923335D+00,
     :  0.1899686386727D-09, 0.2061077157189D+01, 0.2275259891141D+00,
     :  0.1401994545308D-09, 0.4177771136967D+01, 0.1102062672231D+00 /
      DATA ((S1(I,J,1),I=1,3),J= 21, 30) /
     :  0.1578057810499D-09, 0.5782460597335D+01, 0.7626583626240D-01,
     :  0.1237713253351D-09, 0.5705900866881D+01, 0.5154640627760D+00,
     :  0.1313076837395D-09, 0.5163438179576D+01, 0.3664874755930D-01,
     :  0.1184963304860D-09, 0.3054804427242D+01, 0.6327837846670D+00,
     :  0.1238130878565D-09, 0.2317292575962D+01, 0.3961708870310D-01,
     :  0.1015959527736D-09, 0.2194643645526D+01, 0.7329749511860D-01,
     :  0.9017954423714D-10, 0.2868603545435D+01, 0.1990721704425D+00,
     :  0.8668024955603D-10, 0.4923849675082D+01, 0.5439178814476D+00,
     :  0.7756083930103D-10, 0.3014334135200D+01, 0.9491756770005D+00,
     :  0.7536503401741D-10, 0.2704886279769D+01, 0.1030928125552D+00 /
      DATA ((S1(I,J,1),I=1,3),J= 31, 40) /
     :  0.5483308679332D-10, 0.6010983673799D+01, 0.8531963191132D+00,
     :  0.5184339620428D-10, 0.1952704573291D+01, 0.2093666171530D+00,
     :  0.5108658712030D-10, 0.2958575786649D+01, 0.2172315424036D+00,
     :  0.5019424524650D-10, 0.1736317621318D+01, 0.2164800718209D+00,
     :  0.4909312625978D-10, 0.3167216416257D+01, 0.2101180877357D+00,
     :  0.4456638901107D-10, 0.7697579923471D+00, 0.3235053470014D+00,
     :  0.4227030350925D-10, 0.3490910137928D+01, 0.6373574839730D-01,
     :  0.4095456040093D-10, 0.5178888984491D+00, 0.6470106940028D+00,
     :  0.4990537041422D-10, 0.3323887668974D+01, 0.1422690933580D-01,
     :  0.4321170010845D-10, 0.4288484987118D+01, 0.7358765972222D+00 /
      DATA ((S1(I,J,1),I=1,3),J= 41,NS1X) /
     :  0.3544072091802D-10, 0.6021051579251D+01, 0.5265099800692D+00,
     :  0.3480198638687D-10, 0.4600027054714D+01, 0.5328719641544D+00,
     :  0.3440287244435D-10, 0.4349525970742D+01, 0.8582758298370D-01,
     :  0.3330628322713D-10, 0.2347391505082D+01, 0.1104591729320D-01,
     :  0.2973060707184D-10, 0.4789409286400D+01, 0.5257585094865D+00,
     :  0.2932606766089D-10, 0.5831693799927D+01, 0.5336234347371D+00,
     :  0.2876972310953D-10, 0.2692638514771D+01, 0.1173197218910D+00,
     :  0.2827488278556D-10, 0.2056052487960D+01, 0.2022531624851D+00,
     :  0.2515028239756D-10, 0.7411863262449D+00, 0.9597935788730D-01,
     :  0.2853033744415D-10, 0.3948481024894D+01, 0.2118763888447D+01 /

*  SSB-to-Sun, T^2, X
      DATA ((S2(I,J,1),I=1,3),J=  1,NS2X) /
     :  0.1603551636587D-11, 0.4404109410481D+01, 0.2061856251104D+00,
     :  0.1556935889384D-11, 0.4818040873603D+00, 0.2204125344462D+00,
     :  0.1182594414915D-11, 0.9935762734472D+00, 0.5225775174439D+00,
     :  0.1158794583180D-11, 0.3353180966450D+01, 0.5368044267797D+00,
     :  0.9597358943932D-12, 0.5567045358298D+01, 0.2132990797783D+00,
     :  0.6511516579605D-12, 0.5630872420788D+01, 0.4265981595566D+00,
     :  0.7419792747688D-12, 0.2156188581957D+01, 0.5296909721118D+00,
     :  0.3951972655848D-12, 0.1981022541805D+01, 0.1059381944224D+01,
     :  0.4478223877045D-12, 0.0000000000000D+00, 0.0000000000000D+00 /

*  SSB-to-Sun, T^0, Y
      DATA ((S0(I,J,2),I=1,3),J=  1, 10) /
     :  0.4955392320126D-02, 0.2170467313679D+01, 0.5296909721118D+00,
     :  0.2722325167392D-02, 0.2444433682196D+01, 0.2132990797783D+00,
     :  0.1546579925346D-02, 0.5992779281546D+00, 0.3813291813120D-01,
     :  0.8363140252966D-03, 0.7687356310801D+00, 0.7478166569050D-01,
     :  0.3385792683603D-03, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.1201192221613D-03, 0.2520035601514D+01, 0.1059381944224D+01,
     :  0.7587125720554D-04, 0.1669954006449D+01, 0.4265981595566D+00,
     :  0.1964155361250D-04, 0.5707743963343D+01, 0.2061856251104D+00,
     :  0.1891900364909D-04, 0.2320960679937D+01, 0.2204125344462D+00,
     :  0.1937373433356D-04, 0.3226940689555D+01, 0.1495633313810D+00 /
      DATA ((S0(I,J,2),I=1,3),J= 11, 20) /
     :  0.1437139941351D-04, 0.2301626908096D+01, 0.5225775174439D+00,
     :  0.1406267683099D-04, 0.5188579265542D+01, 0.5368044267797D+00,
     :  0.1178703080346D-04, 0.5489483248476D+01, 0.7626583626240D-01,
     :  0.8079835186041D-05, 0.1683751835264D+01, 0.3664874755930D-01,
     :  0.7623253594652D-05, 0.2656400462961D+01, 0.3961708870310D-01,
     :  0.6248667483971D-05, 0.4992775362055D+01, 0.7329749511860D-01,
     :  0.4366353695038D-05, 0.2869706279678D+01, 0.1589072916335D+01,
     :  0.3829101568895D-05, 0.3572131359950D+01, 0.7113454667900D-02,
     :  0.3175733773908D-05, 0.4535372530045D+01, 0.4194847048887D+00,
     :  0.3092437902159D-05, 0.9230153317909D+00, 0.6398972393349D+00 /
      DATA ((S0(I,J,2),I=1,3),J= 21, 30) /
     :  0.2874168812154D-05, 0.3363143761101D+01, 0.1102062672231D+00,
     :  0.3040119321826D-05, 0.3324250895675D+01, 0.6283075850446D+01,
     :  0.2699723308006D-05, 0.2917882441928D+00, 0.1030928125552D+00,
     :  0.2134832683534D-05, 0.4220997202487D+01, 0.3163918923335D+00,
     :  0.1770412139433D-05, 0.4747318496462D+01, 0.1021328554739D+02,
     :  0.1377264209373D-05, 0.4305058462401D+00, 0.1484170571900D-02,
     :  0.1127814538960D-05, 0.8538177240740D+00, 0.6327837846670D+00,
     :  0.1055608090130D-05, 0.1551800742580D+01, 0.4337116142245D+00,
     :  0.9802673861420D-06, 0.1459646735377D+01, 0.1052268489556D+01,
     :  0.1090329461951D-05, 0.1587351228711D+01, 0.1162474756779D+01 /
      DATA ((S0(I,J,2),I=1,3),J= 31, 40) /
     :  0.6959590025090D-06, 0.5534442628766D+01, 0.1066495398892D+01,
     :  0.5664914529542D-06, 0.6030673003297D+01, 0.9491756770005D+00,
     :  0.6607787763599D-06, 0.4989507233927D+01, 0.8460828644453D+00,
     :  0.6269725742838D-06, 0.4222951804572D+01, 0.1480791608091D+00,
     :  0.6301889697863D-06, 0.5444316669126D+01, 0.2243449970715D+00,
     :  0.4891042662861D-06, 0.1490552839784D+01, 0.3340612434717D+01,
     :  0.3457083123290D-06, 0.3030475486049D+01, 0.3516457698740D-01,
     :  0.3032559967314D-06, 0.2652038793632D+01, 0.1104591729320D-01,
     :  0.2841133988903D-06, 0.1276744786829D+01, 0.4110125927500D-01,
     :  0.2855564444432D-06, 0.2143368674733D+01, 0.1510475019529D+00 /
      DATA ((S0(I,J,2),I=1,3),J= 41, 50) /
     :  0.2765157135038D-06, 0.5444186109077D+01, 0.6373574839730D-01,
     :  0.2382312465034D-06, 0.2190521137593D+01, 0.2275259891141D+00,
     :  0.2808060365077D-06, 0.5735195064841D+01, 0.2535050500000D-01,
     :  0.2332175234405D-06, 0.9481985524859D-01, 0.7181332454670D-01,
     :  0.2322488199659D-06, 0.5180499361533D+01, 0.8582758298370D-01,
     :  0.1881850258423D-06, 0.3219788273885D+01, 0.2118763888447D+01,
     :  0.2196111392808D-06, 0.2366941159761D+01, 0.2968341143800D-02,
     :  0.2183810335519D-06, 0.4825445110915D+01, 0.7775000683430D-01,
     :  0.2002733093326D-06, 0.2457148995307D+01, 0.2093666171530D+00,
     :  0.1967111767229D-06, 0.5586291545459D+01, 0.2172315424036D+00 /
      DATA ((S0(I,J,2),I=1,3),J= 51, 60) /
     :  0.1568473250543D-06, 0.3708003123320D+01, 0.7429900518901D+00,
     :  0.1852528314300D-06, 0.4310638151560D+01, 0.2022531624851D+00,
     :  0.1832111226447D-06, 0.1494665322656D+01, 0.3235053470014D+00,
     :  0.1746805502310D-06, 0.1451378500784D+01, 0.1385174140878D+00,
     :  0.1555730966650D-06, 0.1068040418198D+01, 0.7358765972222D+00,
     :  0.1554883462559D-06, 0.2442579035461D+01, 0.5154640627760D+00,
     :  0.1638380568746D-06, 0.2597913420625D+00, 0.8531963191132D+00,
     :  0.1159938593640D-06, 0.5834512021280D+01, 0.1990721704425D+00,
     :  0.1083427965695D-06, 0.5054033177950D+01, 0.5439178814476D+00,
     :  0.1156480369431D-06, 0.5325677432457D+01, 0.5257585094865D+00 /
      DATA ((S0(I,J,2),I=1,3),J= 61, 70) /
     :  0.1141308860095D-06, 0.2153403923857D+01, 0.5336234347371D+00,
     :  0.7913146470946D-07, 0.8642846847027D+00, 0.1478866649112D+01,
     :  0.7439752463733D-07, 0.1970628496213D+01, 0.2164800718209D+00,
     :  0.7280277104079D-07, 0.6073307250609D+01, 0.2101180877357D+00,
     :  0.8319567719136D-07, 0.1954371928334D+01, 0.1692165728891D+01,
     :  0.7137705549290D-07, 0.8904989440909D+00, 0.4155522422634D+00,
     :  0.6900825396225D-07, 0.2825717714977D+01, 0.1173197218910D+00,
     :  0.7245757216635D-07, 0.2481677513331D+01, 0.1265567569334D+01,
     :  0.6961165696255D-07, 0.1292955312978D+01, 0.9562891316684D+00,
     :  0.7571804456890D-07, 0.3427517575069D+01, 0.1422690933580D-01 /
      DATA ((S0(I,J,2),I=1,3),J= 71, 80) /
     :  0.6605425721904D-07, 0.8052192701492D+00, 0.6470106940028D+00,
     :  0.7375477357248D-07, 0.1705076390088D+01, 0.1581959461667D+01,
     :  0.7041664951470D-07, 0.4848356967891D+00, 0.9597935788730D-01,
     :  0.6322199535763D-07, 0.3878069473909D+01, 0.7084920306520D-01,
     :  0.5244380279191D-07, 0.2645560544125D+01, 0.5265099800692D+00,
     :  0.5143125704988D-07, 0.4834486101370D+01, 0.5328719641544D+00,
     :  0.5871866319373D-07, 0.7981472548900D+00, 0.7871412831580D-01,
     :  0.6300822573871D-07, 0.5979398788281D+01, 0.2608790314060D+02,
     :  0.6062154271548D-07, 0.4108655402756D+01, 0.1114304132498D+00,
     :  0.4361912339976D-07, 0.5322624319280D+01, 0.1375773836557D+01 /
      DATA ((S0(I,J,2),I=1,3),J= 81, 90) /
     :  0.4417005920067D-07, 0.6240817359284D+01, 0.2770348281756D+00,
     :  0.4686806749936D-07, 0.3214977301156D+01, 0.1143987543936D+00,
     :  0.3758892132305D-07, 0.5879809634765D+01, 0.1596186371003D+01,
     :  0.5151351332319D-07, 0.2893377688007D+00, 0.2228608264996D+00,
     :  0.4554683578572D-07, 0.5475427144122D+01, 0.1465949902372D+00,
     :  0.3442381385338D-07, 0.5992034796640D+01, 0.5070101000000D-01,
     :  0.2831093954933D-07, 0.5367350273914D+01, 0.3092784376656D+00,
     :  0.3756267090084D-07, 0.5758171285420D+01, 0.4903339079539D+00,
     :  0.2816374679892D-07, 0.1863718700923D+01, 0.2991266627620D+00,
     :  0.3419307025569D-07, 0.9524347534130D+00, 0.3518164938661D+00 /
      DATA ((S0(I,J,2),I=1,3),J= 91,100) /
     :  0.2904250494239D-07, 0.5304471615602D+01, 0.1099462426779D+00,
     :  0.2471734511206D-07, 0.1297069793530D+01, 0.6256703299991D+00,
     :  0.2539620831872D-07, 0.3281126083375D+01, 0.1256615170089D+02,
     :  0.2281017868007D-07, 0.1829122133165D+01, 0.6681224869435D+01,
     :  0.2275319473335D-07, 0.5797198160181D+01, 0.3932462625300D-02,
     :  0.2547755368442D-07, 0.4752697708330D+01, 0.1169588211447D+01,
     :  0.2285979669317D-07, 0.1223205292886D+01, 0.1045155034888D+01,
     :  0.1913386560994D-07, 0.1757532993389D+01, 0.1155361302111D+01,
     :  0.1809020525147D-07, 0.4246116108791D+01, 0.3368040641550D-01,
     :  0.1649213300201D-07, 0.1445162890627D+01, 0.4408250688924D+00 /
      DATA ((S0(I,J,2),I=1,3),J=101,110) /
     :  0.1834972793932D-07, 0.1126917567225D+01, 0.4452511715700D-02,
     :  0.1439550648138D-07, 0.6160756834764D+01, 0.9420622223326D+00,
     :  0.1487645457041D-07, 0.4358761931792D+01, 0.4123712502208D+00,
     :  0.1731729516660D-07, 0.6134456753344D+01, 0.2108507877249D+00,
     :  0.1717747163567D-07, 0.1898186084455D+01, 0.2157473718317D+00,
     :  0.1418190430374D-07, 0.4180286741266D+01, 0.6521991896920D-01,
     :  0.1404844134873D-07, 0.7654053565412D-01, 0.4258542984690D-01,
     :  0.1409842846538D-07, 0.4418612420312D+01, 0.2258291676434D+00,
     :  0.1090948346291D-07, 0.1260615686131D+01, 0.4226656969313D+00,
     :  0.1357577323612D-07, 0.3558248818690D+01, 0.7923417740620D-01 /
      DATA ((S0(I,J,2),I=1,3),J=111,120) /
     :  0.1018154061960D-07, 0.5676087241256D+01, 0.1456308687557D+00,
     :  0.1412073972109D-07, 0.8394392632422D+00, 0.1525316725248D+00,
     :  0.1030938326496D-07, 0.1653593274064D+01, 0.1795258541446D+01,
     :  0.1180081567104D-07, 0.1285802592036D+01, 0.7032915397480D-01,
     :  0.9708510575650D-08, 0.7631889488106D+00, 0.8434341241180D-01,
     :  0.9637689663447D-08, 0.4630642649176D+01, 0.1272681024002D+01,
     :  0.1068910429389D-07, 0.5294934032165D+01, 0.2123349582968D+00,
     :  0.1063716179336D-07, 0.2736266800832D+01, 0.2142632012598D+00,
     :  0.1234858713814D-07, 0.1302891146570D+01, 0.1847279083684D+00,
     :  0.8912631189738D-08, 0.3570415993621D+01, 0.2648454860559D+01 /
      DATA ((S0(I,J,2),I=1,3),J=121,130) /
     :  0.1036378285534D-07, 0.4236693440949D+01, 0.1370332435159D+00,
     :  0.9667798501561D-08, 0.2960768892398D+01, 0.4376440768498D+00,
     :  0.8108314201902D-08, 0.6987781646841D+00, 0.2880807454688D+00,
     :  0.7648364324628D-08, 0.2499017863863D+01, 0.2037373330570D+00,
     :  0.7286136828406D-08, 0.3787426951665D+01, 0.1129145838217D+00,
     :  0.9448237743913D-08, 0.2694354332983D+01, 0.5272426800584D+00,
     :  0.9374276106428D-08, 0.4787121277064D+01, 0.5321392641652D+00,
     :  0.7100226287462D-08, 0.3530238792101D+00, 0.6288513220417D+00,
     :  0.9253056659571D-08, 0.1399478925664D+01, 0.1606092486742D+00,
     :  0.6636432145504D-08, 0.3479575438447D+01, 0.1368660381889D+01 /
      DATA ((S0(I,J,2),I=1,3),J=131,140) /
     :  0.6469975312932D-08, 0.1383669964800D+01, 0.2008557621224D+01,
     :  0.7335849729765D-08, 0.1243698166898D+01, 0.9561746721300D-02,
     :  0.8743421205855D-08, 0.3776164289301D+01, 0.3801276407308D+00,
     :  0.5993635744494D-08, 0.5627122113596D+01, 0.2042657109477D+02,
     :  0.5981008479693D-08, 0.1674336636752D+01, 0.2111650433779D+01,
     :  0.6188535145838D-08, 0.5214925208672D+01, 0.4305306221819D+00,
     :  0.6596074017566D-08, 0.2907653268124D+01, 0.1063314406849D+01,
     :  0.6630815126226D-08, 0.2127643669658D+01, 0.8389694097774D+00,
     :  0.6156772830040D-08, 0.5082160803295D+01, 0.4234171675140D+00,
     :  0.6446960563014D-08, 0.1872100916905D+01, 0.5287268506303D+00 /
      DATA ((S0(I,J,2),I=1,3),J=141,150) /
     :  0.6429324424668D-08, 0.5610276103577D+01, 0.5306550935933D+00,
     :  0.6302232396465D-08, 0.1592152049607D+01, 0.1253008786510D-01,
     :  0.6399244436159D-08, 0.2746214421532D+01, 0.5217580628120D+02,
     :  0.5474965172558D-08, 0.2317666374383D+01, 0.2221856701002D+01,
     :  0.5339293190692D-08, 0.1084724961156D+01, 0.7466759693650D-01,
     :  0.5334733683389D-08, 0.3594106067745D+01, 0.7489573444450D-01,
     :  0.5392665782110D-08, 0.5630254365606D+01, 0.1055449481598D+01,
     :  0.6682075673789D-08, 0.1518480041732D+01, 0.2213766559277D+00,
     :  0.5079130495960D-08, 0.2739765115711D+01, 0.2132517061319D+00,
     :  0.5077759793261D-08, 0.5290711290094D+01, 0.2133464534247D+00 /
      DATA ((S0(I,J,2),I=1,3),J=151,160) /
     :  0.4832037368310D-08, 0.1404473217200D+01, 0.7160067364790D-01,
     :  0.6463279674802D-08, 0.6038381695210D+01, 0.2209183458640D-01,
     :  0.6240592771560D-08, 0.1290170653666D+01, 0.3306188016693D+00,
     :  0.4672013521493D-08, 0.3261895939677D+01, 0.7796265773310D-01,
     :  0.6500650750348D-08, 0.1154522312095D+01, 0.3884652414254D+00,
     :  0.6344161389053D-08, 0.6206111545062D+01, 0.7605151500000D-01,
     :  0.4682518370646D-08, 0.5409118796685D+01, 0.1073608853559D+01,
     :  0.5329460015591D-08, 0.1202985784864D+01, 0.7287631425543D+00,
     :  0.5701588675898D-08, 0.4098715257064D+01, 0.8731175355560D-01,
     :  0.6030690867211D-08, 0.4132033218460D+00, 0.9846002785331D+00 /
      DATA ((S0(I,J,2),I=1,3),J=161,170) /
     :  0.4336256312655D-08, 0.1211415991827D+01, 0.4297791515992D+00,
     :  0.4688498808975D-08, 0.3765479072409D+01, 0.2127790306879D+00,
     :  0.4675578609335D-08, 0.4265540037226D+01, 0.2138191288687D+00,
     :  0.4225578112158D-08, 0.5237566010676D+01, 0.3407705765729D+00,
     :  0.5139422230028D-08, 0.1507173079513D+01, 0.7233337363710D-01,
     :  0.4619995093571D-08, 0.9023957449848D-01, 0.8603097737811D+00,
     :  0.4494776255461D-08, 0.5414930552139D+00, 0.7381754420900D-01,
     :  0.4274026276788D-08, 0.4145735303659D+01, 0.7574578717200D-01,
     :  0.5018141789353D-08, 0.3344408829055D+01, 0.3180992042600D-02,
     :  0.4866163952181D-08, 0.3348534657607D+01, 0.7722995774390D-01 /
      DATA ((S0(I,J,2),I=1,3),J=171,180) /
     :  0.4111986020501D-08, 0.4198823597220D+00, 0.1451108196653D+00,
     :  0.3356142784950D-08, 0.5609144747180D+01, 0.1274714967946D+00,
     :  0.4070575554551D-08, 0.7028411059224D+00, 0.3503323232942D+00,
     :  0.3257451857278D-08, 0.5624697983086D+01, 0.5296435984654D+00,
     :  0.3256973703026D-08, 0.1857842076707D+01, 0.5297383457582D+00,
     :  0.3830771508640D-08, 0.4562887279931D+01, 0.9098186128426D+00,
     :  0.3725024005962D-08, 0.2358058692652D+00, 0.1084620721060D+00,
     :  0.3136763921756D-08, 0.2049731526845D+01, 0.2346394437820D+00,
     :  0.3795147256194D-08, 0.2432356296933D+00, 0.1862120789403D+00,
     :  0.2877342229911D-08, 0.5631101279387D+01, 0.1905464808669D+01 /
      DATA ((S0(I,J,2),I=1,3),J=181,190) /
     :  0.3076931798805D-08, 0.1117615737392D+01, 0.3628624111593D+00,
     :  0.2734765945273D-08, 0.5899826516955D+01, 0.2131850110243D+00,
     :  0.2733405296885D-08, 0.2130562964070D+01, 0.2134131485323D+00,
     :  0.2898552353410D-08, 0.3462387048225D+00, 0.5291709230214D+00,
     :  0.2893736103681D-08, 0.8534352781543D+00, 0.5302110212022D+00,
     :  0.3095717734137D-08, 0.2875061429041D+01, 0.2976424921901D+00,
     :  0.2636190425832D-08, 0.2242512846659D+01, 0.1485980103780D+01,
     :  0.3645512095537D-08, 0.1354016903958D+01, 0.6044726378023D+00,
     :  0.2808173547723D-08, 0.6705114365631D-01, 0.6225157782540D-01,
     :  0.2625012866888D-08, 0.4775705748482D+01, 0.5268983110410D-01 /
      DATA ((S0(I,J,2),I=1,3),J=191,200) /
     :  0.2572233995651D-08, 0.2638924216139D+01, 0.1258454114666D+01,
     :  0.2604238824792D-08, 0.4826358927373D+01, 0.2103781122809D+00,
     :  0.2596886385239D-08, 0.3200388483118D+01, 0.2162200472757D+00,
     :  0.3228057304264D-08, 0.5384848409563D+01, 0.2007689919132D+00,
     :  0.2481601798252D-08, 0.5173373487744D+01, 0.1062562936266D+01,
     :  0.2745977498864D-08, 0.6250966149853D+01, 0.5651155736444D+00,
     :  0.2669878833811D-08, 0.4906001352499D+01, 0.1400015846597D+00,
     :  0.3203986611711D-08, 0.5034333010005D+01, 0.7036329877322D+00,
     :  0.3354961227212D-08, 0.6108262423137D+01, 0.4549093064213D+00,
     :  0.2400407324558D-08, 0.2135399294955D+01, 0.2125476091956D+00 /
      DATA ((S0(I,J,2),I=1,3),J=201,210) /
     :  0.2379905859802D-08, 0.5893721933961D+01, 0.2140505503610D+00,
     :  0.2550844302187D-08, 0.3331940762063D+01, 0.1534957940063D+00,
     :  0.2268824211001D-08, 0.1843418461035D+01, 0.2235935264888D+00,
     :  0.2464700891204D-08, 0.3029548547230D+01, 0.2091065926078D+00,
     :  0.2436814726024D-08, 0.4994717970364D+01, 0.2174915669488D+00,
     :  0.2443623894745D-08, 0.2645102591375D+01, 0.1739420156204D+00,
     :  0.2318701783838D-08, 0.5700547397897D+01, 0.7530171478090D-01,
     :  0.2284448700256D-08, 0.5268898905872D+01, 0.7426161660010D-01,
     :  0.2468848123510D-08, 0.5276280575078D+01, 0.2526561439362D+00,
     :  0.2814052350303D-08, 0.6130168623475D+01, 0.5636314030725D+00 /
      DATA ((S0(I,J,2),I=1,3),J=211,NS0Y) /
     :  0.2243662755220D-08, 0.6631692457995D+00, 0.8886590321940D-01,
     :  0.2330795855941D-08, 0.2499435487702D+01, 0.1056200952181D+01,
     :  0.9757679038404D-09, 0.5796846023126D+01, 0.7826370942180D+02 /

*  SSB-to-Sun, T^1, Y
      DATA ((S1(I,J,2),I=1,3),J=  1, 10) /
     :  0.8989047573576D-08, 0.5840593672122D+01, 0.4265981595566D+00,
     :  0.7815938401048D-08, 0.1129664707133D+01, 0.2061856251104D+00,
     :  0.7550926713280D-08, 0.6196589104845D+00, 0.2204125344462D+00,
     :  0.6056556925895D-08, 0.1677494667846D+01, 0.1059381944224D+01,
     :  0.5734142698204D-08, 0.4000920852962D+01, 0.5225775174439D+00,
     :  0.5614341822459D-08, 0.3486722577328D+01, 0.5368044267797D+00,
     :  0.1028678147656D-08, 0.1877141024787D+01, 0.7113454667900D-02,
     :  0.7270792075266D-09, 0.5077167301739D+01, 0.6398972393349D+00,
     :  0.8734141726040D-09, 0.9069550282609D-01, 0.4194847048887D+00,
     :  0.5377371402113D-09, 0.6039381844671D+01, 0.4337116142245D+00 /
      DATA ((S1(I,J,2),I=1,3),J= 11, 20) /
     :  0.4729719431571D-09, 0.2153086311760D+01, 0.2132990797783D+00,
     :  0.4458052820973D-09, 0.5059830025565D+01, 0.5296909721118D+00,
     :  0.4406855467908D-09, 0.2027971692630D+01, 0.1589072916335D+01,
     :  0.3101659310977D-09, 0.3317677981860D+01, 0.1052268489556D+01,
     :  0.3016749232545D-09, 0.3913703482532D+01, 0.1066495398892D+01,
     :  0.3198541352656D-09, 0.1275513098525D+01, 0.1495633313810D+00,
     :  0.2142065389871D-09, 0.5301351614597D+01, 0.3163918923335D+00,
     :  0.1902615247592D-09, 0.4894943352736D+00, 0.2275259891141D+00,
     :  0.1613410990871D-09, 0.2449891130437D+01, 0.1102062672231D+00,
     :  0.1576992165097D-09, 0.4211421447633D+01, 0.7626583626240D-01 /
      DATA ((S1(I,J,2),I=1,3),J= 21, 30) /
     :  0.1241637259894D-09, 0.4140803368133D+01, 0.5154640627760D+00,
     :  0.1313974830355D-09, 0.3591920305503D+01, 0.3664874755930D-01,
     :  0.1181697118258D-09, 0.1506314382788D+01, 0.6327837846670D+00,
     :  0.1238239742779D-09, 0.7461405378404D+00, 0.3961708870310D-01,
     :  0.1010107068241D-09, 0.6271010795475D+00, 0.7329749511860D-01,
     :  0.9226316616509D-10, 0.1259158839583D+01, 0.1990721704425D+00,
     :  0.8664946419555D-10, 0.3353244696934D+01, 0.5439178814476D+00,
     :  0.7757230468978D-10, 0.1447677295196D+01, 0.9491756770005D+00,
     :  0.7693168628139D-10, 0.1120509896721D+01, 0.1030928125552D+00,
     :  0.5487897454612D-10, 0.4439380426795D+01, 0.8531963191132D+00 /
      DATA ((S1(I,J,2),I=1,3),J= 31, 40) /
     :  0.5196118677218D-10, 0.3788856619137D+00, 0.2093666171530D+00,
     :  0.5110853339935D-10, 0.1386879372016D+01, 0.2172315424036D+00,
     :  0.5027804534813D-10, 0.1647881805466D+00, 0.2164800718209D+00,
     :  0.4922485922674D-10, 0.1594315079862D+01, 0.2101180877357D+00,
     :  0.6155599524400D-10, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.4447147832161D-10, 0.5480720918976D+01, 0.3235053470014D+00,
     :  0.4144691276422D-10, 0.1931371033660D+01, 0.6373574839730D-01,
     :  0.4099950625452D-10, 0.5229611294335D+01, 0.6470106940028D+00,
     :  0.5060541682953D-10, 0.1731112486298D+01, 0.1422690933580D-01,
     :  0.4293615946300D-10, 0.2714571038925D+01, 0.7358765972222D+00 /
      DATA ((S1(I,J,2),I=1,3),J= 41,NS1Y) /
     :  0.3545659845763D-10, 0.4451041444634D+01, 0.5265099800692D+00,
     :  0.3479112041196D-10, 0.3029385448081D+01, 0.5328719641544D+00,
     :  0.3438516493570D-10, 0.2778507143731D+01, 0.8582758298370D-01,
     :  0.3297341285033D-10, 0.7898709807584D+00, 0.1104591729320D-01,
     :  0.2972585818015D-10, 0.3218785316973D+01, 0.5257585094865D+00,
     :  0.2931707295017D-10, 0.4260731012098D+01, 0.5336234347371D+00,
     :  0.2897198149403D-10, 0.1120753978101D+01, 0.1173197218910D+00,
     :  0.2832293240878D-10, 0.4597682717827D+00, 0.2022531624851D+00,
     :  0.2864348326612D-10, 0.2169939928448D+01, 0.9597935788730D-01,
     :  0.2852714675471D-10, 0.2377659870578D+01, 0.2118763888447D+01 /

*  SSB-to-Sun, T^2, Y
      DATA ((S2(I,J,2),I=1,3),J=  1,NS2Y) /
     :  0.1609114495091D-11, 0.2831096993481D+01, 0.2061856251104D+00,
     :  0.1560330784946D-11, 0.5193058213906D+01, 0.2204125344462D+00,
     :  0.1183535479202D-11, 0.5707003443890D+01, 0.5225775174439D+00,
     :  0.1158183066182D-11, 0.1782400404928D+01, 0.5368044267797D+00,
     :  0.1032868027407D-11, 0.4036925452011D+01, 0.2132990797783D+00,
     :  0.6540142847741D-12, 0.4058241056717D+01, 0.4265981595566D+00,
     :  0.7305236491596D-12, 0.6175401942957D+00, 0.5296909721118D+00,
     : -0.5580725052968D-12, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.3946122651015D-12, 0.4108265279171D+00, 0.1059381944224D+01 /

*  SSB-to-Sun, T^0, Z
      DATA ((S0(I,J,3),I=1,3),J=  1, 10) /
     :  0.1181255122986D-03, 0.4607918989164D+00, 0.2132990797783D+00,
     :  0.1127777651095D-03, 0.4169146331296D+00, 0.5296909721118D+00,
     :  0.4777754401806D-04, 0.4582657007130D+01, 0.3813291813120D-01,
     :  0.1129354285772D-04, 0.5758735142480D+01, 0.7478166569050D-01,
     : -0.1149543637123D-04, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.3298730512306D-05, 0.5978801994625D+01, 0.4265981595566D+00,
     :  0.2733376706079D-05, 0.7665413691040D+00, 0.1059381944224D+01,
     :  0.9426389657270D-06, 0.3710201265838D+01, 0.2061856251104D+00,
     :  0.8187517749552D-06, 0.3390675605802D+00, 0.2204125344462D+00,
     :  0.4080447871819D-06, 0.4552296640088D+00, 0.5225775174439D+00 /
      DATA ((S0(I,J,3),I=1,3),J= 11, 20) /
     :  0.3169973017028D-06, 0.3445455899321D+01, 0.5368044267797D+00,
     :  0.2438098615549D-06, 0.5664675150648D+01, 0.3664874755930D-01,
     :  0.2601897517235D-06, 0.1931894095697D+01, 0.1495633313810D+00,
     :  0.2314558080079D-06, 0.3666319115574D+00, 0.3961708870310D-01,
     :  0.1962549548002D-06, 0.3167411699020D+01, 0.7626583626240D-01,
     :  0.2180518287925D-06, 0.1544420746580D+01, 0.7113454667900D-02,
     :  0.1451382442868D-06, 0.1583756740070D+01, 0.1102062672231D+00,
     :  0.1358439007389D-06, 0.5239941758280D+01, 0.6398972393349D+00,
     :  0.1050585898028D-06, 0.2266958352859D+01, 0.3163918923335D+00,
     :  0.1050029870186D-06, 0.2711495250354D+01, 0.4194847048887D+00 /
      DATA ((S0(I,J,3),I=1,3),J= 21, 30) /
     :  0.9934920679800D-07, 0.1116208151396D+01, 0.1589072916335D+01,
     :  0.1048395331560D-06, 0.3408619600206D+01, 0.1021328554739D+02,
     :  0.8370147196668D-07, 0.3810459401087D+01, 0.2535050500000D-01,
     :  0.7989856510998D-07, 0.3769910473647D+01, 0.7329749511860D-01,
     :  0.5441221655233D-07, 0.2416994903374D+01, 0.1030928125552D+00,
     :  0.4610812906784D-07, 0.5858503336994D+01, 0.4337116142245D+00,
     :  0.3923022803444D-07, 0.3354170010125D+00, 0.1484170571900D-02,
     :  0.2610725582128D-07, 0.5410600646324D+01, 0.6327837846670D+00,
     :  0.2455279767721D-07, 0.6120216681403D+01, 0.1162474756779D+01,
     :  0.2375530706525D-07, 0.6055443426143D+01, 0.1052268489556D+01 /
      DATA ((S0(I,J,3),I=1,3),J= 31, 40) /
     :  0.1782967577553D-07, 0.3146108708004D+01, 0.8460828644453D+00,
     :  0.1581687095238D-07, 0.6255496089819D+00, 0.3340612434717D+01,
     :  0.1594657672461D-07, 0.3782604300261D+01, 0.1066495398892D+01,
     :  0.1563448615040D-07, 0.1997775733196D+01, 0.2022531624851D+00,
     :  0.1463624258525D-07, 0.1736316792088D+00, 0.3516457698740D-01,
     :  0.1331585056673D-07, 0.4331941830747D+01, 0.9491756770005D+00,
     :  0.1130634557637D-07, 0.6152017751825D+01, 0.2968341143800D-02,
     :  0.1028949607145D-07, 0.2101792614637D+00, 0.2275259891141D+00,
     :  0.1024074971618D-07, 0.4071833211074D+01, 0.5070101000000D-01,
     :  0.8826956060303D-08, 0.4861633688145D+00, 0.2093666171530D+00 /
      DATA ((S0(I,J,3),I=1,3),J= 41, 50) /
     :  0.8572230171541D-08, 0.5268190724302D+01, 0.4110125927500D-01,
     :  0.7649332643544D-08, 0.5134543417106D+01, 0.2608790314060D+02,
     :  0.8581673291033D-08, 0.2920218146681D+01, 0.1480791608091D+00,
     :  0.8430589300938D-08, 0.3604576619108D+01, 0.2172315424036D+00,
     :  0.7776165501012D-08, 0.3772942249792D+01, 0.6373574839730D-01,
     :  0.8311070234408D-08, 0.6200412329888D+01, 0.3235053470014D+00,
     :  0.6927365212582D-08, 0.4543353113437D+01, 0.8531963191132D+00,
     :  0.6791574208598D-08, 0.2882188406238D+01, 0.7181332454670D-01,
     :  0.5593100811839D-08, 0.1776646892780D+01, 0.7429900518901D+00,
     :  0.4553381853021D-08, 0.3949617611240D+01, 0.7775000683430D-01 /
      DATA ((S0(I,J,3),I=1,3),J= 51, 60) /
     :  0.5758000450068D-08, 0.3859251775075D+01, 0.1990721704425D+00,
     :  0.4281283457133D-08, 0.1466294631206D+01, 0.2118763888447D+01,
     :  0.4206935661097D-08, 0.5421776011706D+01, 0.1104591729320D-01,
     :  0.4213751641837D-08, 0.3412048993322D+01, 0.2243449970715D+00,
     :  0.5310506239878D-08, 0.5421641370995D+00, 0.5154640627760D+00,
     :  0.3827450341320D-08, 0.8887314524995D+00, 0.1510475019529D+00,
     :  0.4292435241187D-08, 0.1405043757194D+01, 0.1422690933580D-01,
     :  0.3189780702289D-08, 0.1060049293445D+01, 0.1173197218910D+00,
     :  0.3226611928069D-08, 0.6270858897442D+01, 0.2164800718209D+00,
     :  0.2893897608830D-08, 0.5117563223301D+01, 0.6470106940028D+00 /
      DATA ((S0(I,J,3),I=1,3),J= 61,NS0Z) /
     :  0.3239852024578D-08, 0.4079092237983D+01, 0.2101180877357D+00,
     :  0.2956892222200D-08, 0.1594917021704D+01, 0.3092784376656D+00,
     :  0.2980177912437D-08, 0.5258787667564D+01, 0.4155522422634D+00,
     :  0.3163725690776D-08, 0.3854589225479D+01, 0.8582758298370D-01,
     :  0.2662262399118D-08, 0.3561326430187D+01, 0.5257585094865D+00,
     :  0.2766689135729D-08, 0.3180732086830D+00, 0.1385174140878D+00,
     :  0.2411600278464D-08, 0.3324798335058D+01, 0.5439178814476D+00,
     :  0.2483527695131D-08, 0.4169069291947D+00, 0.5336234347371D+00,
     :  0.7788777276590D-09, 0.1900569908215D+01, 0.5217580628120D+02 /

*  SSB-to-Sun, T^1, Z
      DATA ((S1(I,J,3),I=1,3),J=  1, 10) /
     :  0.5444220475678D-08, 0.1803825509310D+01, 0.2132990797783D+00,
     :  0.3883412695596D-08, 0.4668616389392D+01, 0.5296909721118D+00,
     :  0.1334341434551D-08, 0.0000000000000D+00, 0.0000000000000D+00,
     :  0.3730001266883D-09, 0.5401405918943D+01, 0.2061856251104D+00,
     :  0.2894929197956D-09, 0.4932415609852D+01, 0.2204125344462D+00,
     :  0.2857950357701D-09, 0.3154625362131D+01, 0.7478166569050D-01,
     :  0.2499226432292D-09, 0.3657486128988D+01, 0.4265981595566D+00,
     :  0.1937705443593D-09, 0.5740434679002D+01, 0.1059381944224D+01,
     :  0.1374894396320D-09, 0.1712857366891D+01, 0.5368044267797D+00,
     :  0.1217248678408D-09, 0.2312090870932D+01, 0.5225775174439D+00 /
      DATA ((S1(I,J,3),I=1,3),J= 11,NS1Z) /
     :  0.7961052740870D-10, 0.5283368554163D+01, 0.3813291813120D-01,
     :  0.4979225949689D-10, 0.4298290471860D+01, 0.4194847048887D+00,
     :  0.4388552286597D-10, 0.6145515047406D+01, 0.7113454667900D-02,
     :  0.2586835212560D-10, 0.3019448001809D+01, 0.6398972393349D+00 /

*  SSB-to-Sun, T^2, Z
      DATA ((S2(I,J,3),I=1,3),J=  1,NS2Z) /
     :  0.3749920358054D-12, 0.3230285558668D+01, 0.2132990797783D+00,
     :  0.2735037220939D-12, 0.6154322683046D+01, 0.5296909721118D+00 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Time since reference epoch, years.
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJY
      T2 = T*T

*  Set status.
      IF ( ABS(T) .LE. 100D0 ) THEN
         JSTAT = 0
      ELSE
         JSTAT = 1
      END IF

*  X then Y then Z.
      DO 7 K=1,3

*     Initialize position and velocity component.
         XYZ = 0D0
         XYZD = 0D0

*     ------------------------------------------------
*     Obtain component of Sun to Earth ecliptic vector
*     ------------------------------------------------

*     Sun to Earth, T^0 terms.
         DO 1 J=1,NE0(K)
            A = E0(1,J,K)
            B = E0(2,J,K)
            C = E0(3,J,K)
            P = B + C*T
            XYZ  = XYZ  + A*COS(P)
            XYZD = XYZD - A*C*SIN(P)
 1       CONTINUE

*     Sun to Earth, T^1 terms.
         DO 2 J=1,NE1(K)
            A = E1(1,J,K)
            B = E1(2,J,K)
            C = E1(3,J,K)
            CT = C*T
            P = B + CT
            CP = COS(P)
            XYZ  = XYZ  + A*T*CP
            XYZD = XYZD + A*(CP-CT*SIN(P))
 2       CONTINUE

*     Sun to Earth, T^2 terms.
         DO 3 J=1,NE2(K)
            A = E2(1,J,K)
            B = E2(2,J,K)
            C = E2(3,J,K)
            CT = C*T
            P = B + CT
            CP = COS(P)
            XYZ  = XYZ  + A*T2*CP
            XYZD = XYZD + A*T*(2D0*CP-CT*SIN(P))
 3       CONTINUE

*     Heliocentric Earth position and velocity component.
         PH(K) = XYZ
         VH(K) = XYZD / DJY

*     ------------------------------------------------
*     Obtain component of SSB to Earth ecliptic vector
*     ------------------------------------------------

*     SSB to Sun, T^0 terms.
         DO 4 J=1,NS0(K)
            A = S0(1,J,K)
            B = S0(2,J,K)
            C = S0(3,J,K)
            P = B + C*T
            XYZ  = XYZ  + A*COS(P)
            XYZD = XYZD - A*C*SIN(P)
 4       CONTINUE

*     SSB to Sun, T^1 terms.
         DO 5 J=1,NS1(K)
            A = S1(1,J,K)
            B = S1(2,J,K)
            C = S1(3,J,K)
            CT = C*T
            P = B + CT
            CP = COS(P)
            XYZ  = XYZ  + A*T*CP
            XYZD = XYZD + A*(CP-CT*SIN(P))
 5       CONTINUE

*     SSB to Sun, T^2 terms.
         DO 6 J=1,NS2(K)
            A = S2(1,J,K)
            B = S2(2,J,K)
            C = S2(3,J,K)
            CT = C*T
            P = B + CT
            CP = COS(P)
            XYZ  = XYZ  + A*T2*CP
            XYZD = XYZD + A*T*(2D0*CP-CT*SIN(P))
 6       CONTINUE

*     Barycentric Earth position and velocity component.
         PB(K) = XYZ
         VB(K) = XYZD / DJY

*     Next Cartesian component.
 7    CONTINUE

*  Rotate from ecliptic to ICRF coordinates and return the results.
      X = PH(1)
      Y = PH(2)
      Z = PH(3)
      PVH(1,1) =      X + AM12*Y + AM13*Z
      PVH(2,1) = AM21*X + AM22*Y + AM23*Z
      PVH(3,1) =          AM32*Y + AM33*Z
      X = VH(1)
      Y = VH(2)
      Z = VH(3)
      PVH(1,2) =      X + AM12*Y + AM13*Z
      PVH(2,2) = AM21*X + AM22*Y + AM23*Z
      PVH(3,2) =          AM32*Y + AM33*Z
      X = PB(1)
      Y = PB(2)
      Z = PB(3)
      PVB(1,1) =      X + AM12*Y + AM13*Z
      PVB(2,1) = AM21*X + AM22*Y + AM23*Z
      PVB(3,1) =          AM32*Y + AM33*Z
      X = VB(1)
      Y = VB(2)
      Z = VB(3)
      PVB(1,2) =      X + AM12*Y + AM13*Z
      PVB(2,2) = AM21*X + AM22*Y + AM23*Z
      PVB(3,2) =          AM32*Y + AM33*Z

*  Finished.

*+-----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  Standards Of Fundamental Astronomy Review Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Review Board ("the Board").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and 
*     restrictions listed below.
*
*  3. You (the user) may copy and adapt the SOFA software and its 
*     algorithms for your own purposes and you may copy and distribute
*     a resulting "derived work" to others on a world-wide, royalty-free 
*     basis, provided that the derived work complies with the following
*     requirements: 
*
*     a) Your work shall be marked or carry a statement that it (i) uses
*        routines and computations derived by you from software provided 
*        by SOFA under license to you; and (ii) does not contain
*        software provided by SOFA or software that has been distributed
*        by or endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon and/or differs from the
*        original SOFA software.
*
*     c) The name(s) of all routine(s) that you distribute shall differ
*        from the SOFA names, even when the SOFA content has not been
*        otherwise changed.
*
*     d) The routine-naming prefix "iau" shall not be used.
*
*     e) The origin of the SOFA components of your derived work must not
*        be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     f) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have granted 
*        a further right to modify the source code of your derived work.
*
*  4. In any published work or commercial products which includes
*     results achieved by using the SOFA software, you shall acknowledge
*     that the SOFA software was used in obtaining those results.
*
*  5. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or by
*     inappropriate modification.
*
*  6. The SOFA software is provided "as is" and the Board makes no 
*     warranty as to its use or performance.   The Board does not and 
*     cannot warrant the performance or results which the user may obtain 
*     by using the SOFA software.  The Board makes no warranties, express 
*     or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Board be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Board representative has been advised of such
*     damages, or for any claim by any third party.
*
*  7. The provision of any version of the SOFA software under the terms 
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*     Internet email: sofa@rl.ac.uk
*     Postal address: IAU SOFA Center
*                     Rutherford Appleton Laboratory
*                     Chilton, Didcot, Oxon OX11 0QX
*                     United Kingdom
*
*-----------------------------------------------------------------------

      END
