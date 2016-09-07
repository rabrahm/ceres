#include <sofa.h>
#include <sofam.h>
#include <string.h>
#include "math.h"

void ptwGeoc ( double p, double h, double *r, double *z )
/*
**  - - - - - - - -
**   p t w G e o c
**  - - - - - - - -
**
**  Convert geodetic position to geocentric.
**
**  (double precision)
**
**  Given:
**     p     double     latitude (geodetic, radians)
**     h     double     height above reference spheroid (geodetic, metres)
**
**  Returned:
**     *r    double     distance from Earth axis (AU)
**     *z    double     distance from plane of Earth equator (AU)
**
**  Notes:
**
**  1  Geocentric latitude can be obtained by evaluating atan2(z,r).
**
**  2  IAU 1976 constants are used.
**
**  Reference:
**
**     Green,R.M., Spherical Astronomy, CUP 1985, p98.
**
**  Last revision:   22 July 2004
**
**  Author P.T.Wallace.  Free software, no restrictions.
*/
{
   double sp, cp, c, s;


/* Earth equatorial radius (metres) */
   static const double a0 = 6378140.0;

/* Reference spheroid flattening factor and useful function thereof */
   static const double f = 1.0 / 298.257;
   double b = ( 1.0 - f ) * ( 1.0 - f );

/* Astronomical unit in metres */
   static const double au = 1.49597870e11;

/* Geodetic to geocentric conversion */
   sp = sin ( p );
   cp = cos ( p );
   c = 1.0 / sqrt ( cp * cp + b * sp * sp );
   s = b * c;
   *r = ( a0 * c + h ) * cp / au;
   *z = ( a0 * s + h ) * sp / au;
}

void ptwEqecl ( double dr, double dd, double date,
                double *dl, double *db )
/*
**  - - - - - - - - -
**   p t w E q e c l
**  - - - - - - - - -
**
**  Transformation from J2000.0 equatorial coordinates to ecliptic
**  coordinates (IAU 1976/80).
**
**  (double precision)
**
**  Given:
**     dr,dd       double      J2000.0 mean RA,Dec (radians)
**     date        double      TDB (loosely ET) as Modified Julian Date
**                                              (JD-2400000.5)
**  Returned:
**     *dl,*db     double      ecliptic longitude and latitude
**                             (mean of date, IAU 1980 theory, radians)
**
**  Called:
**     iauS2c, iauPmat76, iauRxp, iauIr, iauObl80, iauRx, slaDcc2s,
**     iauAnp, iauAnpm
**
**
**  Last revision:   2009 December 29
**
**  Author P.T.Wallace.
*/
{
   double rmat[3][3], v1[3], v2[3];


/* Spherical to Cartesian */
   iauS2c ( dr, dd, v1 );

/* Mean J2000.0 to mean of date */
   iauPmat76 ( 2400000.5, date, rmat );
   iauRxp ( rmat, v1, v2 );

/* Equatorial to ecliptic */
   iauIr ( rmat );
   iauRx ( iauObl80 ( 2400000.5, date ), rmat );
   iauRxp ( rmat, v2, v1 );

/* Cartesian to spherical */
   iauC2s ( v1, dl, db );

/* Express in conventional ranges */
   *dl = iauAnp ( *dl );
   *db = iauAnpm ( *db );
}

void ptwFk45z ( double r1950, double d1950, double bepoch,
                double *r2000, double *d2000 )
/*
**  - - - - - - - - -
**   p t w F k 4 5 z
**  - - - - - - - - -
**
**  Convert B1950.0 FK4 star data to J2000.0 FK5 assuming zero
**  proper motion in the FK5 frame (double precision)
**
**  This function converts stars from the old, Bessel-Newcomb, FK4
**  system to the newer, IAU 1976, FK5, Fricke system, in such a
**  way that the FK5 proper motion is zero.  Because such a star
**  has, in general, a non-zero proper motion in the FK4 system,
**  the function requires the epoch at which the position in the
**  FK4 system was determined.
**
**  The method is from Appendix 2 of Ref 1, but using the constants
**  of Ref 4.
**
**  Given:
**     r1950,d1950     double   B1950.0 FK4 RA,Dec at epoch (rad)
**     bepoch          double   Besselian epoch (e.g. 1979.3)
**
**  Returned:
**     *r2000,*d2000   double   J2000.0 FK5 RA,Dec (rad)
**
**  Notes:
**
**  1)  The epoch BEPOCH is strictly speaking Besselian, but
**      if a Julian epoch is supplied the result will be
**      affected only to a negligible extent.
**
**  2)  Conversion from Besselian epoch 1950.0 to Julian epoch
**      2000.0 only is provided for.  Conversions involving other
**      epochs will require use of the appropriate precession,
**      proper motion, and E-terms functions before and/or
**      after FK45Z is called.
**
**  3)  In the FK4 catalogue the proper motions of stars within
**      10 degrees of the poles do not embody the differential
**      E-term effect and should, strictly speaking, be handled
**      in a different manner from stars outside these regions.
**      However, given the general lack of homogeneity of the star
**      data available for routine astrometry, the difficulties of
**      handling positions that may have been determined from
**      astrometric fields spanning the polar and non-polar regions,
**      the likelihood that the differential E-terms effect was not
**      taken into account when allowing for proper motion in past
**      astrometry, and the undesirability of a discontinuity in
**      the algorithm, the decision has been made in this function to
**      include the effect of differential E-terms on the proper
**      motions for all stars, whether polar or not.  At epoch 2000,
**      and measuring on the sky rather than in terms of dRA, the
**      errors resulting from this simplification are less than
**      1 milliarcsecond in position and 1 milliarcsecond per
**      century in proper motion.
**
**  References:
**
**     1  Aoki,S., et al, 1983, Astron. Astrophys., 128, 263
**
**     2  Smith, C.A. et al, 1989, "The transformation of astrometric
**        catalog systems to the equinox J2000.0", Astron.J. 97, 265
**
**     3  Yallop, B.D. et al, 1989, "Transformation of mean star places
**        from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space",
**        Astron.J. 97, 274
**
**     4  Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to
**        the Astronomical Almanac", ISBN 0-935702-68-7
**
**  Called:  iauS2c, iauEpj, iauEpb2jd, iauC2s, iauAnp
**
**  Defined in sofam.h:  D2PI
**
**  Last revision:   30 December 2009
**
**  Author:   P.T.Wallace.
*/
{
   double w, date1, date2;
   int i, j;

/* Position and position+velocity vectors */
   double r0[3], a1[3], v1[3], v2[6];

/* Radians per year to arcsec per century */
   static double pmf = 100.0 * 60.0 * 60.0 * 360.0 / D2PI;

/*
** Canonical constants  (see references)
*/

/* vectors a and adot, and matrix m (only half of which is needed here) */
   static double a[3]  = { -1.62557e-6,  -0.31919e-6, -0.13843e-6 };
   static double ad[3] = {  1.245e-3,    -1.580e-3,   -0.659e-3 };
   static double em[6][3] =
   {
     {  0.9999256782, -0.0111820611, -0.0048579477 },
     {  0.0111820610,  0.9999374784, -0.0000271765 },
     {  0.0048579479, -0.0000271474,  0.9999881997 },
     { -0.000551,     -0.238565,      0.435739     },
     {  0.238514,     -0.002667,     -0.008541     },
     { -0.435623,      0.012254,      0.002117     }
   };

/* Spherical to Cartesian */
   iauS2c ( r1950, d1950, r0 );

/* Adjust vector a to give zero proper motion in FK5 */
   w = ( bepoch - 1950.0 ) / pmf;
   for ( i = 0; i < 3; i++ ) {
      a1[i] = a[i] + w * ad[i];
   }

/* Remove e-terms */
   w = r0[0] * a1[0] + r0[1] * a1[1] + r0[2] * a1[2];
   for ( i = 0; i < 3; i++ ) {
      v1[i] = r0[i] - a1[i] + w * r0[i];
   }

/* Convert position vector to Fricke system */
   for ( i = 0; i < 6; i++ ) {
      w = 0.0;
      for ( j = 0; j < 3; j++ ) {
         w += em[i][j] * v1[j];
      }
      v2[i] = w;
   }

/* Allow for fictitious proper motion in FK4 */
   iauEpb2jd ( bepoch, &date1, &date2 );
   w = ( iauEpj ( date1, date2 ) - 2000.0 ) / pmf;
   for ( i = 0; i < 3; i++ ) {
      v2[i] += w * v2[i+3];
   }

/* Revert to spherical coordinates */
   iauC2s ( v2, &w, d2000 );
   *r2000 = iauAnp ( w );
}

/*--------------------------------------------------------------------*/

void ptwFk425 ( double r1950, double d1950, double dr1950,
                double dd1950, double p1950, double v1950,
                double *r2000, double *d2000, double *dr2000,
                double *dd2000, double *p2000, double *v2000 )
/*
**  - - - - - - - - -
**   p t w F k 4 2 5
**  - - - - - - - - -
**
**  Convert B1950.0 FK4 star data to J2000.0 FK5.
**
**  (double precision)
**
**  This function converts stars from the old, Bessel-Newcomb, FK4
**  system to the newer, IAU 1976, FK5, Fricke system.  The precepts
**  of Smith et al (Ref 1) are followed, using the implementation
**  by Yallop et al (Ref 2) of a matrix method due to Standish.
**  Kinoshita's development of Andoyer's post-Newcomb precession is
**  used.  The numerical constants from Seidelmann et al (Ref 3) are
**  used canonically.
**
**  Given:  (all B1950.0,FK4)
**
**     r1950,d1950     double    B1950.0 RA,dec (rad)
**     dr1950,dd1950   double    B1950.0 proper motions (rad/trop.yr)
**     p1950           double    parallax (arcsec)
**     v1950           double    radial velocity (km/s, +ve = moving away)
**
**  Returned:  (all J2000.0,FK5)
**
**     *r2000,*d2000   double    J2000.0 RA,dec (rad)
**     *dr2000,*dd2000 double    J2000.0 proper motions (rad/jul.yr)
**     *p2000          double    parallax (arcsec)
**     *v2000          double    radial velocity (km/s, +ve = moving away)
**
**  Notes:
**
**  1)  The proper motions in RA are dRA/dt rather than
**      cos(Dec)*dRA/dt, and are per year rather than per century.
**
**  2)  Conversion from Besselian epoch 1950.0 to Julian epoch
**      2000.0 only is provided for.  Conversions involving other
**      epochs will require use of the appropriate precession,
**      proper motion, and E-terms functions before and/or
**      after FK425 is called.
**
**  3)  In the FK4 catalogue the proper motions of stars within
**      10 degrees of the poles do not embody the differential
**      E-term effect and should, strictly speaking, be handled
**      in a different manner from stars outside these regions.
**      However, given the general lack of homogeneity of the star
**      data available for routine astrometry, the difficulties of
**      handling positions that may have been determined from
**      astrometric fields spanning the polar and non-polar regions,
**      the likelihood that the differential E-terms effect was not
**      taken into account when allowing for proper motion in past
**      astrometry, and the undesirability of a discontinuity in
**      the algorithm, the decision has been made in this function to
**      include the effect of differential E-terms on the proper
**      motions for all stars, whether polar or not.  At epoch 2000,
**      and measuring on the sky rather than in terms of dRA, the
**      errors resulting from this simplification are less than
**      1 milliarcsecond in position and 1 milliarcsecond per
**      century in proper motion.
**
**  References:
**
**     1  Smith, C.A. et al, 1989, "The transformation of astrometric
**        catalog systems to the equinox J2000.0", Astron.J. 97, 265
**
**     2  Yallop, B.D. et al, 1989, "Transformation of mean star places
**        from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space",
**        Astron.J. 97, 274
**
**     3  Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to
**        the Astronomical Almanac", ISBN 0-935702-68-7
**
**  Defined in sofam.h:  D2PI
**
**  Last revision:   30 December 2009
**
**  Author P.T.Wallace.
*/
{
   double r, d, ur, ud, px, rv, sr, cr, sd, cd, w, wd,
          x, y, z, xd, yd, zd,
          rxysq, rxyzsq, rxy, rxyz, spxy, spxyz;
   int i, j;

/* Star position and velocity vectors */
   double r0[3], rd0[3];

/* Combined position and velocity vectors */
   double v1[6], v2[6];

/* Radians per year to arcsec per century */
   static double pmf = 100.0 * 60.0 * 60.0 * 360.0 / D2PI;

/* Small number to avoid arithmetic problems */
   double tiny = 1.0e-30;

/*
** Canonical constants  (see references)
*/

/*
** Km per sec to AU per tropical century
** = 86400 * 36524.2198782 / 1.49597870e8
*/
   double vf = 21.095;

/* Constant vector and matrix (by rows) */
   static double a[3]  = { -1.62557e-6,  -0.31919e-6, -0.13843e-6 };
   static double ad[3] = {  1.245e-3,     -1.580e-3,   -0.659e-3  };
   static double em[6][6] =
   {
     {  0.9999256782,              /* em[0][0] */
       -0.0111820611,              /* em[0][1] */
       -0.0048579477,              /* em[0][2] */
        0.00000242395018,          /* em[0][3] */
       -0.00000002710663,          /* em[0][4] */
       -0.00000001177656 },        /* em[0][5] */

     {  0.0111820610,              /* em[1][0] */
        0.9999374784,              /* em[1][1] */
       -0.0000271765,              /* em[1][2] */
        0.00000002710663,          /* em[1][3] */
        0.00000242397878,          /* em[1][4] */
       -0.00000000006587 },        /* em[1][5] */

     {  0.0048579479,              /* em[2][0] */
       -0.0000271474,              /* em[2][1] */
        0.9999881997,              /* em[2][2] */
        0.00000001177656,          /* em[2][3] */
       -0.00000000006582,          /* em[2][4] */
        0.00000242410173 },        /* em[2][5] */

     { -0.000551,                  /* em[3][0] */
       -0.238565,                  /* em[3][1] */
        0.435739,                  /* em[3][2] */
        0.99994704,                /* em[3][3] */
       -0.01118251,                /* em[3][4] */
       -0.00485767 },              /* em[3][5] */

     {  0.238514,                  /* em[4][0] */
       -0.002667,                  /* em[4][1] */
       -0.008541,                  /* em[4][2] */
        0.01118251,                /* em[4][3] */
        0.99995883,                /* em[4][4] */
       -0.00002718 },              /* em[4][5] */

     { -0.435623,                  /* em[5][0] */
        0.012254,                  /* em[5][1] */
        0.002117,                  /* em[5][2] */
        0.00485767,                /* em[5][3] */
       -0.00002714,                /* em[5][4] */
        1.00000956 }               /* em[5][5] */
   };

/* Pick up B1950 data (units radians and arcsec/tc) */
   r = r1950;
   d = d1950;
   ur = dr1950 * pmf;
   ud = dd1950 * pmf;
   px = p1950;
   rv = v1950;

/* Spherical to Cartesian */
   sr = sin ( r );
   cr = cos ( r );
   sd = sin ( d );
   cd = cos ( d );

   r0[0] = cr * cd;
   r0[1] = sr * cd;
   r0[2] = sd;

   w = vf * rv * px;

   rd0[0] = ( -sr * cd * ur ) - ( cr * sd * ud ) + ( w * r0[0] );
   rd0[1] = ( cr * cd * ur ) - ( sr * sd * ud ) + ( w * r0[1] );
   rd0[2] = ( cd * ud ) + ( w * r0[2] );

/* Allow for e-terms and express as position+velocity 6-vector */
   w = ( r0[0] * a[0] ) + ( r0[1] * a[1] ) + ( r0[2] * a[2] );
   wd = ( r0[0] * ad[0] ) + ( r0[1] * ad[1] ) + ( r0[2] * ad[2] );

   for ( i = 0; i < 3; i++ ) {
      v1[i] = r0[i]  - a[i]  + w * r0[i];
      v1[i+3] = rd0[i] - ad[i] + wd * r0[i];
   }

/* Convert position+velocity vector to Fricke system */
   for ( i = 0; i < 6; i++ ) {
      w = 0.0;
      for ( j = 0; j < 6; j++ ) {
         w += em[i][j] * v1[j];
      }
      v2[i] = w;
   }

/* Revert to spherical coordinates */
   x = v2[0];
   y = v2[1];
   z = v2[2];
   xd = v2[3];
   yd = v2[4];
   zd = v2[5];

   rxysq = ( x * x ) + ( y * y );
   rxyzsq = ( rxysq ) + ( z * z );
   rxy = sqrt ( rxysq );
   rxyz = sqrt (  rxyzsq );

   spxy = ( x * xd ) + ( y * yd );
   spxyz = spxy + ( z * zd );

   r = ( x != 0.0 || y != 0.0 ) ? atan2 ( y, x ) : 0.0;
   if ( r < 0.0 ) r += D2PI;
   d = atan2 ( z, rxy );

   if ( rxy > tiny ) {
      ur = ( ( x * yd ) - ( y * xd ) ) / rxysq;
      ud = ( ( zd * rxysq ) - ( z * spxy ) ) / ( rxyzsq * rxy );
   }

   if ( px > tiny ) {
      rv = spxyz / ( px * rxyz * vf );
      px = px / rxyz;
   }

/* Return results */
   *r2000 = r;
   *d2000 = d;
   *dr2000 = ur / pmf;
   *dd2000 = ud / pmf;
   *v2000 = rv;
   *p2000 = px;
}

/*--------------------------------------------------------------------*/

void ptwObs ( int n, char *c, char *name, double *w, double *p, double *h )
/*
**  - - - - - - -
**   p t w O b s
**  - - - - - - -
**
**  Parameters of selected groundbased observing stations.
**
**  Given:
**     n       int     number specifying observing station
**
**  Either given or returned
**    *c       char    identifier specifying observing station
**
**  Returned:
**    *name    char    name of specified observing station
**    *w       double  longitude (radians, West +ve)
**    *p       double  geodetic latitude (radians, North +ve)
**    *h       double  height above sea level (metres)
**
**  Notes:
**
**     Station identifiers may be up to 10 characters long
**     (plus string terminator), and station names may be up to
**     40 characters long (plus string terminator).  Leading or
**     trailing spaces are not supported.
**
**     c and n are alternative ways of specifying the observing
**     station.  The c option, which is the most generally useful,
**     may be selected by specifying an n value of zero or less.
**     If n is 1 or more, the parameters of the nth station
**     in the currently supported list are interrogated (n=1
**     meaning the first station in the list), and the station
**     identifier c is returned as well as name, w, p and h.
**
**     If the station parameters are not available, either because
**     the station identifier c is not recognized, or because an
**     n value greater than the number of stations supported is
**     given, a name of "?" is returned and c, w, p and h are left
**     in their current states.
**
**     Programs can obtain a list of all currently supported
**     stations by calling the function repeatedly, with n=1,2,3...
**     When name="?" is seen, the list of stations has been
**     exhausted.
**
**     Station numbers, identifiers, names and other details are
**     subject to change and should not be hardwired into
**     application programs.
**
**     All station identifiers c are uppercase only;  lowercase
**     characters must be converted to uppercase by the calling
**     program.  The station names returned may contain both upper-
**     and lowercase.  All characters up to the first space are
**     checked;  thus an abbreviated ID will return the parameters
**     for the first station in the list which matches the
**     abbreviation supplied, and no station in the list will ever
**     contain embedded spaces.  c must not have leading spaces.
**
**     IMPORTANT -- BEWARE OF THE LONGITUDE SIGN CONVENTION.  The
**     longitude returned by slaObs is west-positive in accordance
**     with astronomical usage.  However, this sign convention is
**     left-handed and is the opposite of the one used by geographers;
**     elsewhere in slalib the preferable east-positive convention is
**     used.  In particular, note that for use in slaAop, slaAoppa
**     and slaOap the sign of the longitude must be reversed.
**
**     Users are urged to inform the author of any improvements
**     they would like to see made.  For example:
**
**         typographical corrections
**         more accurate parameters
**         better station identifiers or names
**         additional stations
**
**  Defined in sofam.h:  DAS2R
**
**  Last revision:   2009 December 29
**
**  Author P.T.Wallace.
*/

#define WEST(id,iam,as) ( DAS2R * \
            ( (double) (60L * ( (long) (60 * (id)) +(iam) ) ) + (as) ) )
#define NORTH(id,iam,as) ( WEST(id,iam,as) )
#define EAST(id,iam,as) ( -WEST(id,iam,as) )
#define SOUTH(id,iam,as) ( -WEST(id,iam,as) )

{
/* ------------------- Table of station parameters ------------------- */

   static struct station {
      char *id;                 /* identifier */
      char *na;                 /* name */
      double wlong;             /* longitude (west) */
      double phi;               /* latitude */
      double hm;                /* height ASL (metres) */
   } statab[] = {

/* AAT (Observer's Guide) */
      {
         "AAT",
         "Anglo-Australian 3.9m Telescope",
         EAST(149, 3,57.91),
         SOUTH(31,16,37.34),
         1164.0
      },

/* WHT (Gemini, April 1987) */
      {
         "LPO4.2",
         "William Herschel 4.2m Telescope",
         WEST(17,52,53.9),
         NORTH(28,45,38.1),
         2332.0
      },

/* INT (Gemini, April 1987) */
      {
         "LPO2.5",
         "Isaac Newton 2.5m Telescope",
         WEST(17,52,39.5),
         NORTH(28,45,43.2),
         2336.0
      },

/* JKT (Gemini, April 1987) */
      {
         "LPO1",
         "Jacobus Kapteyn 1m Telescope",
         WEST(17,52,41.2),
         NORTH(28,45,39.9),
         2364.0
      },

/* Lick 120" (S.L.Allen, private communication, 2002) */
      {
         "LICK120",
         "Lick 120 inch",
         WEST(121,38,13.689),
         NORTH(37,20,34.931),
         1286.0
      },

/* MMT 6.5m conversion (MMT Observatory website) */
      {
         "MMT",
         "MMT 6.5m, Mt Hopkins",
         WEST(110,53, 4.4),
         NORTH(31,41,19.6),
         2608.0
      },

/* Victoria B.C. 1.85m (1984 Almanac) */
      {
         "VICBC",
         "Victoria B.C. 1.85 metre",
         WEST(123,25, 1.18),
         NORTH(48,31,11.9),
         238.0
      },

/* Las Campanas (1983 Almanac) */
      {
         "DUPONT",
         "Du Pont 2.5m Telescope, Las Campanas",
         WEST(70,42,9.),
         SOUTH(29, 0,11.),
         2280.0
      },

/* Mt Hopkins 1.5m (1983 Almanac) */
      {
         "MTHOP1.5",
         "Mt Hopkins 1.5 metre",
         WEST(110,52,39.00),
         NORTH(31,40,51.4),
         2344.0
      },

/* Mt Stromlo 74" (1983 Almanac) */
      {
         "STROMLO74",
         "Mount Stromlo 74 inch",
         EAST(149, 0,27.59),
         SOUTH(35,19,14.3),
         767.0
      },

/* ANU 2.3m, SSO (Gary Hovey) */
      {
         "ANU2.3",
         "Siding Spring 2.3 metre",
         EAST(149, 3,40.3),
         SOUTH(31,16,24.1),
         1149.0
      },

/* Greenbank 140' (1983 Almanac) */
      {
         "GBVA140",
         "Greenbank 140 foot",
         WEST(79,50, 9.61),
         NORTH(38,26,15.4),
         881.0
      },

/* Cerro Tololo 4m (1982 Almanac) */
      {
         "TOLOLO4M",
         "Cerro Tololo 4 metre",
         WEST(70,48,53.6),
         SOUTH(30, 9,57.8),
         2235.0
      },

/* Cerro Tololo 1.5m (1982 Almanac) */
      {
         "TOLOLO1.5M",
         "Cerro Tololo 1.5 metre",
         WEST(70,48,54.5),
         SOUTH(30, 9,56.3),
         2225.0
      },

/* Tidbinbilla 64m (1982 Almanac) */
      {
         "TIDBINBLA",
         "Tidbinbilla 64 metre",
         EAST(148,58,48.20),
         SOUTH(35,24,14.3),
         670.0
      },

/* Bloemfontein 1.52m (1981 Almanac) */
      {
         "BLOEMF",
         "Bloemfontein 1.52 metre",
         EAST(26,24,18.),
         SOUTH(29, 2,18.),
         1387.0
      },

/* Bosque Alegre 1.54m (1981 Almanac) */
      {
         "BOSQALEGRE",
         "Bosque Alegre 1.54 metre",
         WEST(64,32,48.0),
         SOUTH(31,35,53.),
         1250.0
      },

/* USNO 61" astrographic reflector, Flagstaff (1981 Almanac) */
      {
         "FLAGSTF61",
         "USNO 61 inch astrograph, Flagstaff",
         WEST(111,44,23.6),
         NORTH(35,11, 2.5),
         2316.0
      },

/* Lowell 72" (1981 Almanac) */
      {
         "LOWELL72",
         "Perkins 72 inch, Lowell",
         WEST(111,32, 9.3),
         NORTH(35, 5,48.6),
         2198.0
      },

/* Harvard 1.55m (1981 Almanac) */
      {
         "HARVARD",
         "Harvard College Observatory 1.55m",
         WEST(71,33,29.32),
         NORTH(42,30,19.0),
         185.0
      },

/* Okayama 1.88m (1981 Almanac) */
      {
         "OKAYAMA",
         "Okayama 1.88 metre",
         EAST(133,35,47.29),
         NORTH(34,34,26.1),
         372.0
      },

/* Kitt Peak Mayall 4m (1981 Almanac) */
      {
         "KPNO158",
         "Kitt Peak 158 inch",
         WEST(111,35,57.61),
         NORTH(31,57,50.3),
         2120.0
      },

/* Kitt Peak 90 inch (1981 Almanac) */
      {
         "KPNO90",
         "Kitt Peak 90 inch",
         WEST(111,35,58.24),
         NORTH(31,57,46.9),
         2071.0
      },

/* Kitt Peak 84 inch (1981 Almanac) */
      {
         "KPNO84",
         "Kitt Peak 84 inch",
         WEST(111,35,51.56),
         NORTH(31,57,29.2),
         2096.0
      },

/* Kitt Peak 36 foot (1981 Almanac) */
      {
         "KPNO36FT",
         "Kitt Peak 36 foot",
         WEST(111,36,51.12),
         NORTH(31,57,12.1),
         1939.0
      },

/* Kottamia 74" (1981 Almanac) */
      {
         "KOTTAMIA",
         "Kottamia 74 inch",
         EAST(31,49,30.),
         NORTH(29,55,54.),
         476.0
      },

/* La Silla 3.6m (1981 Almanac) */
      {
         "ESO3.6",
         "ESO 3.6 metre",
         WEST(70,43,36.),
         SOUTH(29,15,36.),
         2428.0
      },

/* Mauna Kea 88 inch
   (IfA website, Richard Wainscoat) */
      {
         "MAUNAK88",
         "Mauna Kea 88 inch",
         WEST(155,28, 9.96),
         NORTH(19,49,22.77),
         4213.6
      },

/* UKIRT
   (Ifa website, Richard Wainscoat) */
      {
         "UKIRT",
         "UK Infra Red Telescope",
         WEST(155,28,13.18),
         NORTH(19,49,20.75),
         4198.5
      },

/* Quebec 1.6m (1981 Almanac) */
      {
         "QUEBEC1.6",
         "Quebec 1.6 metre",
         WEST(71, 9, 9.7),
         NORTH(45,27,20.6),
         1114.0
      },

/* Mt Ekar 1.82m (1981 Almanac) */
      {
         "MTEKAR",
         "Mt Ekar 1.82 metre",
         EAST(11,34,15.),
         NORTH(45,50,48.),
         1365.0
      },

/* Mt Lemmon 60" (1981 Almanac) */
      {
         "MTLEMMON60",
         "Mt Lemmon 60 inch",
         WEST(110,42,16.9),
         NORTH(32,26,33.9),
         2790.0
      },

/* Mt Locke 2.7m (1981 Almanac) */
      {
         "MCDONLD2.7",
         "McDonald 2.7 metre",
         WEST(104, 1,17.60),
         NORTH(30,40,17.7),
         2075.0
      },

/* Mt Locke 2.1m (1981 Almanac) */
      {
         "MCDONLD2.1",
         "McDonald 2.1 metre",
         WEST(104, 1,20.10),
         NORTH(30,40,17.7),
         2075.0
      },

/* Palomar 200" (1981 Almanac) */
      {
         "PALOMAR200",
         "Palomar 200 inch",
         WEST(116,51,50.),
         NORTH(33,21,22.),
         1706.0
      },

/* Palomar 60" (1981 Almanac) */
      {
         "PALOMAR60",
         "Palomar 60 inch",
         WEST(116,51,31.),
         NORTH(33,20,56.),
         1706.0
      },

/* David Dunlap 74" (1981 Almanac) */
      {
         "DUNLAP74",
         "David Dunlap 74 inch",
         WEST(79,25,20.),
         NORTH(43,51,46.),
         244.0
      },

/* Haute Provence 1.93m (1981 Almanac) */
      {
         "HPROV1.93",
         "Haute Provence 1.93 metre",
         EAST(5,42,46.75),
         NORTH(43,55,53.3),
         665.0
      },

/* Haute Provence 1.52m (1981 Almanac) */
      {
         "HPROV1.52",
         "Haute Provence 1.52 metre",
         EAST(5,42,43.82),
         NORTH(43,56, 0.2),
         667.0
      },

/* San Pedro Martir 83" (1981 Almanac) */
      {
         "SANPM83",
         "San Pedro Martir 83 inch",
         WEST(115,27,47.),
         NORTH(31, 2,38.),
         2830.0
      },

/* Sutherland 74" (1981 Almanac) */
      {
         "SAAO74",
         "Sutherland 74 inch",
         EAST(20,48,44.3),
         SOUTH(32,22,43.4),
         1771.0
      },

/* Tautenburg 2m (1981 Almanac) */
      {
         "TAUTNBG",
         "Tautenburg 2 metre",
         EAST(11,42,45.),
         NORTH(50,58,51.),
         331.0
      },

/* Catalina 61" (1981 Almanac) */
      {
         "CATALINA61",
         "Catalina 61 inch",
         WEST(110,43,55.1),
         NORTH(32,25, 0.7),
         2510.0
      },

/* Steward 90" (1981 Almanac) */
      {
         "STEWARD90",
         "Steward 90 inch",
         WEST(111,35,58.24),
         NORTH(31,57,46.9),
         2071.0
      },

/* Russian 6m (1981 Almanac) */
      {
         "USSR6",
         "USSR 6 metre",
         EAST(41,26,30.0),
         NORTH(43,39,12.),
         2100.0
      },

/* Arecibo 1000' (1981 Almanac) */
      {
         "ARECIBO",
         "Arecibo 1000 foot",
         WEST(66,45,11.1),
         NORTH(18,20,36.6),
         496.0
      },

/* Cambridge 5km (1981 Almanac) */
      {
         "CAMB5KM",
         "Cambridge 5km",
         EAST(0, 2,37.23),
         NORTH(52,10,12.2),
         17.0
      },

/* Cambridge 1 mile (1981 Almanac) */
      {
         "CAMB1MILE",
         "Cambridge 1 mile",
         EAST(0, 2,21.64),
         NORTH(52, 9,47.3),
         17.0
      },

/* Bonn 100m (1981 Almanac) */
      {
         "EFFELSBERG",
         "Effelsberg 100 metre",
         EAST(6,53, 1.5),
         NORTH(50,31,28.6),
         366.0
      },

/* Greenbank 300' (1981 Almanac - defunct) */
      {
         "GBVA300",
         "Greenbank 300 foot",
         WEST(79,50,56.36),
         NORTH(38,25,46.3),
         894.0
      },

/* Jodrell Bank Mk 1 (1981 Almanac) */
      {
         "JODRELL1",
         "Jodrell Bank 250 foot",
         WEST(2,18,25.),
         NORTH(53,14,10.5),
         78.0
      },

/* Australia Telescope Parkes Observatory
   (Peter te Lintel Hekkert) */
      {
         "PARKES",
         "Parkes 64 metre",
         EAST(148,15,44.3591),
         SOUTH(32,59,59.8657),
         391.79
      },

/* VLA (1981 Almanac) */
      {
         "VLA",
         "Very Large Array",
         WEST(107,37, 3.82),
         NORTH(34, 4,43.5),
         2124.0
      },

/* Sugar Grove 150' (1981 Almanac) */
      {
         "SUGARGROVE",
         "Sugar Grove 150 foot",
         WEST(79,16,23.),
         NORTH(38,31,14.),
         705.0
      },

/* Russian 600' (1981 Almanac) */
      {
         "USSR600",
         "USSR 600 foot",
         EAST(41,35,25.5),
         NORTH(43,49,32.),
         973.0
      },

/* Nobeyama 45 metre mm dish (based on 1981 Almanac entry) */
      {
         "NOBEYAMA",
         "Nobeyama 45 metre",
         EAST(138,29,12.),
         NORTH(35,56,19.),
         1350.0
      },

/* James Clerk Maxwell 15 metre mm telescope, Mauna Kea
   (IfA website, Richard Wainscoat, height from I.Coulson) */
      {
         "JCMT",
         "JCMT 15 metre",
         WEST(155,28,37.20),
         NORTH(19,49,22.11),
         4111.0
      },

/* ESO 3.5 metre NTT, La Silla (K.Wirenstrand) */
      {
         "ESONTT",
         "ESO 3.5 metre NTT",
         WEST(70,43, 7.),
         SOUTH(29,15,30.),
         2377.0
      },

/* St Andrews University Observatory (1982 Almanac) */
      {
         "ST.ANDREWS",
         "St Andrews",
         WEST(2,48,52.5),
         NORTH(56,20,12.),
         30.0
      },

/* Apache Point 3.5 metre (R.Owen) */
      {
         "APO3.5",
         "Apache Point 3.5m",
         WEST(105,49,11.56),
         NORTH(32,46,48.96),
         2809.0
      },

/* W.M.Keck Observatory, Telescope 1 (site survey)
   (William Lupton) */
      {
         "KECK1",
         "Keck 10m Telescope #1",
         WEST(155,28,28.99),
         NORTH(19,49,33.41),
         4160.0
      },

/* Tautenberg Schmidt (1983 Almanac) */
      {
         "TAUTSCHM",
         "Tautenberg 1.34 metre Schmidt",
         EAST(11,42,45.0),
         NORTH(50,58,51.0),
         331.0
      },

/* Palomar Schmidt (1981 Almanac) */
      {
         "PALOMAR48",
         "Palomar 48-inch Schmidt",
         WEST(116,51,32.0),
         NORTH(33,21,26.0),
         1706.0
      },

/* UK Schmidt, Siding Spring (1983 Almanac) */
      {
         "UKST",
         "UK 1.2 metre Schmidt, Siding Spring",
         EAST(149, 4,12.8),
         SOUTH(31,16,27.8),
         1145.0
      },

/* Kiso Schmidt, Japan (1981 Almanac) */
      {
         "KISO",
         "Kiso 1.05 metre Schmidt, Japan",
         EAST(137,37,42.2),
         NORTH(35,47,38.7),
         1130.0
      },

/* ESO Schmidt, La Silla (1981 Almanac) */
      {
         "ESOSCHM",
         "ESO 1 metre Schmidt, La Silla",
         WEST(70,43,46.5),
         SOUTH(29,15,25.8),
         2347.0
      },

/* Australia Telescope Compact Array (WGS84 coordinates of Station 35,
   Mark Calabretta) */
      {
         "ATCA",
         "Australia Telescope Compact Array",
         EAST(149,33, 0.500),
         SOUTH(30,18,46.385),
         236.9
      },

/* Australia Telescope Mopra Observatory
   (Peter te Lintel Hekkert) */
      {
         "MOPRA",
         "ATNF Mopra Observatory",
         EAST(149, 5,58.732),
         SOUTH(31,16, 4.451),
         850.0
      },

/* Subaru telescope, Mauna Kea
   (IfA website, Richard Wainscoat) */
      {
         "SUBARU",
         "Subaru 8m telescope",
         WEST(155,28,33.67),
         NORTH(19,49,31.81),
         4163.0
      },

/* Canada-France-Hawaii Telescope, Mauna Kea
   (IfA website, Richard Wainscoat) */
      {
         "CFHT",
         "Canada-France-Hawaii 3.6m Telescope",
         WEST(155,28, 7.95),
         NORTH(19,49,30.91),
         4204.1
      },

/* W.M.Keck Observatory, Telescope 2
   (William Lupton) */
      {
         "KECK2",
         "Keck 10m Telescope #2",
         WEST(155,28,27.24),
         NORTH(19,49,35.62),
         4159.6
      },

/* Gemini North, Mauna Kea
   (IfA website, Richard Wainscoat) */
      {
         "GEMININ",
         "Gemini North 8-m telescope",
         WEST(155,28, 8.57),
         NORTH(19,49,25.69),
         4213.4
      },

/* Five College Radio Astronomy Observatory
   (Tim Jenness) */
      {
         "FCRAO",
         "Five College Radio Astronomy Obs",
         WEST(72,20,42.0),
         NORTH(42,23,30.0),
         314.0
      },

/* NASA Infra Red Telescope Facility
   (IfA website, Richard Wainscoat) */
      {
         "IRTF",
         "NASA IR Telescope Facility, Mauna Kea",
         WEST(155,28,19.20),
         NORTH(19,49,34.39),
         4168.1
      },

/* Caltech Submillimeter Observatory
   (IfA website, Richard Wainscoat; height estimated) */
      {
         "CSO",
         "Caltech Sub-mm Observatory, Mauna Kea",
         WEST(155,28,31.79),
         NORTH(19,49,20.78),
         4080.0
      },

/* ESO VLT, UT1
   (ESO website, VLT Whitebook Chapter 2) */
      {
         "VLT1",
         "ESO VLT, Paranal, Chile: UT1",
         WEST(70,24,11.642),
         SOUTH(24,37,33.117),
         2635.43
      },

/* ESO VLT, UT2
   (ESO website, VLT Whitebook Chapter 2) */
      {
         "VLT2",
         "ESO VLT, Paranal, Chile: UT2",
         WEST(70,24,10.855),
         SOUTH(24,37,31.465),
         2635.43
      },

/* ESO VLT, UT3
   (ESO website, VLT Whitebook Chapter 2) */
      {
         "VLT3",
         "ESO VLT, Paranal, Chile: UT3",
         WEST(70,24, 9.896),
         SOUTH(24,37,30.300),
         2635.43
      },

/* ESO VLT, UT4
   (ESO website, VLT Whitebook Chapter 2) */
      {
         "VLT4",
         "ESO VLT, Paranal, Chile: UT4",
         WEST(70,24, 8.000),
         SOUTH(24,37,31.000),
         2635.43
      },

/* Gemini South, Cerro Pachon
   (GPS readings by Patrick Wallace) */
      {
         "GEMINIS",
         "Gemini South 8-m telescope",
         WEST(70,44,11.5),
         SOUTH(30,14,26.7),
         2738.0
      },

/* Cologne Observatory for Submillimeter Astronomy (KOSMA)
   (Holger Jakob) */
      {
         "KOSMA3M",
         "KOSMA 3m telescope, Gornergrat",
         EAST(7,47,3.48),
         NORTH(45,58,59.772),
         3141.0
      },

/* Magellan 1, 6.5m telescope at Las Campanas
   (Skip Schaller) */
      {
         "MAGELLAN1",
         "Magellan 1, 6.5m, Las Campanas",
         WEST(70,41,31.9),
         SOUTH(29, 0,51.7),
         2408.0
      },

/* Magellan 2, 6.5m telescope at Las Campanas
   (Skip Schaller) */
      {
         "MAGELLAN2",
         "Magellan 2, 6.5m, Las Campanas",
         WEST(70,41,33.5),
         SOUTH(29, 0,50.3),
         2408.0
      },

/* Hobby-Eberly Telescope (www.as.utexas.edu/mcdonald/het/het) */
      {
         "HET",
         "Hobby-Eberly Telescope",
         WEST(104, 0,53.0),
         NORTH(30,40,53.2),
         2026.0
      }
   };

   static int NMAX = ( sizeof statab / sizeof ( struct station ) );

/* ------------------------------------------------------------------- */

   int m, i, ic;


/* Station specified by number or identifier? */
   if ( n > 0 ) {

   /* Station specified by number */
      m = n - 1;
      if ( m < NMAX ) {
         strcpy ( c, statab[m].id );
      }
   } else {

   /* Station specified by identifier:  determine corresponding number */
      for ( i = 0; i < 10; i++ ) {
         ic = (int) c[i];
         if ( ic == ' ' || ic == '\0' ) break;
      }
      for ( m = 0; m < NMAX; m++ ) {
         if ( ! strncmp ( c, statab[m].id, i ) ) {
            break;
         }
      }
   }

/* Return parameters of mth station */
   if ( m < NMAX ) {
      strncpy ( name, statab[m].na, 40 );
      *w = statab[m].wlong;
      *p = statab[m].phi;
      *h = statab[m].hm;
   } else {
      strcpy ( name, "?" );
   }
}

/*--------------------------------------------------------------------*/

void ptwCldj ( int iy, int im, int id, double *djm, int *j )
/*
**  - - - - - - - -
**   p t w C l d j
**  - - - - - - - -
**
**  Gregorian calendar to Modified Julian Date.
**
**  Given:
**     iy,im,id     int    year, month, day in Gregorian calendar
**
**  Returned:
**     *djm         double Modified Julian Date (JD-2400000.5) for 0 hrs
**     *j           int    status:
**                           0 = OK
**                           1 = bad year   (MJD not computed)
**                           2 = bad month  (MJD not computed)
**                           3 = bad day    (MJD computed)
**
**  The year must be -4699 (i.e. 4700BC) or later.
**
**  The algorithm is derived from that of Hatcher 1984 (QJRAS 25, 53-55).
**
**  Last revision:   2009 December 31
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   long iyL, imL;

/* Month lengths in days */
   static int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };



/* Validate year */
   if ( iy < -4699 ) { *j = 1; return; }

/* Validate month */
   if ( ( im < 1 ) || ( im > 12 ) ) { *j = 2; return; }

/* Allow for leap year */
   mtab[1] = ( ( ( iy % 4 ) == 0 ) &&
             ( ( ( iy % 100 ) != 0 ) || ( ( iy % 400 ) == 0 ) ) ) ?
             29 : 28;

/* Validate day */
   *j = ( id < 1 || id > mtab[im-1] ) ? 3 : 0;

/* Lengthen year and month numbers to avoid overflow */
   iyL = (long) iy;
   imL = (long) im;

/* Perform the conversion */
   *djm = (double)
        ( ( 1461L * ( iyL - ( 12L - imL ) / 10L + 4712L ) ) / 4L
        + ( 306L * ( ( imL + 9L ) % 12L ) + 5L ) / 10L
        - ( 3L * ( ( iyL - ( 12L - imL ) / 10L + 4900L ) / 100L ) ) / 4L
        + (long) id - 2399904L );
}

/*--------------------------------------------------------------------*/

void ptwCaldj ( int iy, int im, int id, double *djm, int *j )
/*
**  - - - - - - - - -
**   p t w C a l d j
**  - - - - - - - - -
**
**  Gregorian calendar to Modified Julian Date.
**
**  (Includes century default feature:  use slaCldj for years
**   before 100AD.)
**
**  Given:
**     iy,im,id   int      year, month, day in Gregorian calendar
**
**  Returned:
**     *djm       double   Modified Julian Date (JD-2400000.5) for 0 hrs
**     *j         int      status:
**                           0 = ok
**                           1 = bad year   (MJD not computed)
**                           2 = bad month  (MJD not computed)
**                           3 = bad day    (MJD computed)
**
**  Acceptable years are 00-49, interpreted as 2000-2049,
**                       50-99,     "       "  1950-1999,
**                       100 upwards, interpreted literally.
**
**  Called:  ptwCldj
**
**  Last revision:   2009 December 31
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int ny;

/* Default century if appropriate */
   if ( ( iy >= 0 ) && ( iy <= 49 ) )
      ny = iy + 2000;
   else if ( ( iy >= 50 ) && ( iy <= 99 ) )
      ny = iy + 1900;
   else
      ny = iy;

/* Modified Julian Date */
   ptwCldj ( ny, im, id, djm, j );
}
