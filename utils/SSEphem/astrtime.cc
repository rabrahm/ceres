// This file is part of SSEphem.
//
//  SSEphem is a collection of C/C++ source files that provide an interface
//  to the JPL Solar System Ephemeris.
//  Copyright (C) 2009 Associated Universities, Inc. Washington DC, USA
//
//  SSEphem is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  SSEphem is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with SSEphem.  If not, see <http://www.gnu.org/licenses/>.

//               ASTRONOMICAL TIME CONVERSION ROUTINES

// The functions in this file are written to provide conversions to a number
// of astronomical times starting with the generally distributed time, UTC
// or Coordinated Universal Time.  Since the offset and rate of drift of
// between the earth rotation times (UT1, GMST, GAST, LMST, and LAST) and
// the times derived from atomic standardsd (UTC, TAI, TT, and TDB) are
// not entirely predictable for long periods of time, these routines depend
// on two data files that must be kept up to date with information from the
// International Earth Rotation Service.  The file 'leapsec.tab' contains
// the Modified Julian dates of the leap second additions along with the new
// TAI-UTC value that applies after 0 hours of that date.  The file
// 'iers.tab' contains the UT1-UTC values at 5-day intervals along with
// polar motion offsets for those dates.  New IERS data is published every
// Friday so, in principle, this file should be updated weekly.  With each
// report the IERS publishes a prediction algorithm for some months in the
// future, but the error in this algorithm can build up to a significant
// fraction of an arcsecond in UT1-UTC, so frequent maintenance of the
// 'iers.tab' file is important for telescope pointing.

// Data structures and function prototypes are kept in astrtime.h.

#include <stdio.h>
#include <string.h>
#include "sofa.h"
#include "astrtime.h"
#include "ptw_aux.h"

leap_secT ls[MAX_NUM_LEAP_SEC];              // data from leap second file
int num_leap_sec = 0;                        // number of leap seconds found
                                             // in the file
char leap_sec_file[120] = { "leapsec.tab" }; // leap second table file name
FILE *ls_stream = 0;                         // leap second file stream

iersT iers[MAX_NUM_IERS_OFFSETS];            // data from IERS data file
int num_iers_offsets = 0;                    // number of UT1 and pole offsets
                                             // found in the IERS file
char iers_file[120] = { "iers.tab" };        // IERS table file name
FILE *iers_stream = 0;                       // IERS file stream

// VLBI geocentric station coordinates for the Green Bank 140-ft telescope
//observatoryT obs = { 882880.0208, -4924482.4385, 3944130.6438 };

// geocentric station coordinates for the GBT derived from aips++
// measures longitude, latitude, and elevation.
//observatoryT obs = { 880310.5151, -4912152.1552, 3960048.6809 };

// GBT position from 20-meter to GBT Ties VLBI run equivalent to
// geodetic coordinates using the GRS80 ellipsoid and local geoid height
// of -31.10 meters to get: W  79d50'23.419",  N  38d25'59.265", h=854.83 m
observatoryT obs = { 882589.661,  -4924872.388,  3943729.411 };

// UTC to International Atomic Time (TAI).  By definition, UTC and TAI have
// the same rate, but UTC is offset by an integer number of seconds, called
// leap seconds.  The offset is changed as needed to keep UTC within about
// 0.7 seconds of earth rotation time, UT1.  The ascii file 'leapsec.tab'
// contains the dates (MJD) when the offset changes.
//      UTC = TAI - (number of leap seconds)
// The input parameters are mjd and utc (in fraction of a day), and the
// output arguments are tai_mjd and tai (in fraction of a day).  The return
// values is 1, if successful, or 0, if not.

int utc2tai (long mjd, double utc, long &tai_mjd, double &tai)
{
    tai_mjd = mjd;
    tai = utc;
    tai += get_leap_sec(mjd, utc) / 86400.0;
    if (tai >= 1.0) {
        tai -= 1.0;
        tai_mjd++;
    }
    return 1;
}

// ...and this routine goes the other way.  There's a one-second no man's
// land during the leap second that is not accounted for here.

int tai2utc (long tai_mjd, double tai, long &mjd, double &utc)
{
    mjd = tai_mjd;
    double ls = get_leap_sec(tai_mjd, tai);
    utc = tai - ls / 86400.0;
    if (utc < 0.0) {
	utc += 1.0;
	mjd--;
    }
    if (ls > (get_leap_sec(mjd, utc) + 0.00001)) {
	utc += 1.0 / 86400.0;
	if (utc >= 1.0) {
	    utc -= 1.0;
	    mjd++;
	}
    }
    return 1;
}

// Before atomic clocks, Ephemeris Time (ET) was the closest available
// approximation to a uniform time for planetary motion calculations.  This
// has now been replaced by Terrestrial Dynamic Time (TT), which is tied
// directly to atomic time by a constant offset of 32.184 seconds, exactly.
// The purpose of this offset is to maintain continuity between ET and TT.
// TAI, hence TT, extends back only to 1956.
//     TT = TAI + 32.184 = UTC + (number of leap seconds) + 32.184
// The input parameters are mjd and utc (in fraction of a day), and the
// output arguments are tt_mjd and tt (in fraction of a day).  The return
// values is 1, if successful, or 0, if not.

int utc2tt (long mjd, double utc, long &tt_mjd, double &tt)
{
    long tai_mjd;             // International Atomic Time MJD
    double tai;               // International Atomic Time

    if (!utc2tai (mjd, utc, tai_mjd, tai)) {
        return 0;
    }
    tt = tai + 32.184 / 86400.0;
    tt_mjd = tai_mjd;
    if (tt >= 1.0) {
        tt -= 1.0;
        tt_mjd++;
    }
    return 1;
}

// ...and this routine goes the other way, TT -> UTC.

int tt2utc (long tt_mjd, double tt, long &mjd, double &utc)
{
    long tai_mjd;             // International Atomic Time MJD
    double tai;               // International Atomic Time

    tai = tt - 32.184 / 86400.0;
    tai_mjd = tt_mjd;
    if (tai < 0.0) {
        tai += 1.0;
        tai_mjd--;
    }
    if (!tai2utc (tai_mjd, tai, mjd, utc)) {
        return 0;
    }
    return 1;
}

// Barycentric Dynamic Time (TDB) is the same as as Terrestrial Dynamic Time
// (TT) except for relativistic corrections to move the origin to the solar
// system barycenter.  These corrections amount to as much as about 1.6
// millisends and are periodic with an average of zero over many years.  The
// dominant terms in this correction are have annula and semi-annual
// periods, as given by the equation in the Astronomical almanac.  This
// routine uses a more accurate calculation implemented in the Starlink
// library, which claims an accuracy of a few nanoseconds.
// The input parameters are mjd and utc (in fraction of a day), and the
// output arguments are tdb_mjd and tdb (in fraction of a day).  The return
// values is 1, if successful, or 0, if not.

int utc2tdb (long mjd, double utc, long &tdb_mjd, double &tdb)
{
    long tt_mjd;             // Terrestrial coordinate time MJD
    double tt;               // Terrestrial coordinate time
    if (!utc2tt (mjd, utc, tt_mjd, tt)) {
        return 0;
    }
    long ut1_mjd;            // actual earth rotation time MJD
    double ut1;              // actual earth rotation time
    utc2ut1 (mjd, utc, ut1_mjd, ut1);
    
    double obs_long = get_observatory_longitude();
    double obs_r = get_observatory_radius() / 1000.0;
    double obs_z = get_observatory_z() / 1000.0;
    tdb = tt + iauDtdb(2400000.5, (double)tt_mjd + tt, ut1, obs_long,
		       obs_r, obs_z) / 86400.0;
    tdb_mjd = tt_mjd;
    if (tdb > 1.0) {
        tdb -= 1.0;
        tdb_mjd++;
    }
    if (tdb < 0.0) {
        tdb += 1.0;
        tdb_mjd--;
    }
    return 1;
}

// ...and this routine goes the other way, TDB -> UTC.

int tdb2utc (long tdb_mjd, double tdb, long &mjd, double &utc)
{
    // first an approximate conversion to get UT1
    if (!tt2utc (tdb_mjd, tdb, mjd, utc)) {
        return 0;
    }
    long ut1_mjd;            // actual earth rotation time MJD
    double ut1;              // actual earth rotation time
    utc2ut1 (mjd, utc, ut1_mjd, ut1);
    
    double obs_long = get_observatory_longitude();
    double obs_r = get_observatory_radius() / 1000.0;
    double obs_z = get_observatory_z() / 1000.0;
    long tt_mjd;             // Terrestrial coordinate time MJD
    double tt;               // Terrestrial coordinate time
    // tt = tdb - slaRcc((double)tdb_mjd + tdb, ut1, obs_long, obs_r, obs_z) /
    // 86400.0;
    printf("astrtime.cc: obs_long = %f\n", obs_long);
    tt = tdb - iauDtdb(2400000.5, (double)tdb_mjd + tdb, ut1,
		       -obs_long, obs_r, obs_z);
    tt_mjd = tdb_mjd;
    if (tt > 1.0) {
        tt -= 1.0;
        tt_mjd++;
    }
    if (tt < 0.0) {
        tt += 1.0;
        tt_mjd--;
    }
    if (!tt2utc (tt_mjd, tt, mjd, utc)) {
        return 0;
    }
    return 1;
}

// Universal Time (UT1) is directly proportional to the actual rotation of
// the earth, independent of observing location.  It is the observed
// rotation of the earth (UT0) corrected for the observer's shift in
// longitude due to polar motion.  A day in UT1 is very close to a mean
// solar day.  Since the earth's rotation is not uniform, the rate of UT1 is
// not constant, and its offset from atomic time is continually changing in
// an, as yet, unpredictable way.  As of December 1995, UT1 is drifting
// about 0.8 seconds per year with respect to atomic time (TAI or UTC).
// Since UTC is intentionally incremented by integer seconds (leap seconds)
// to stay within 0.7 seconds of UT1, the difference between UT1 and UTC is
// never greater than this.  The difference, DUT1 = UT1 - UTC is monitored
// by the International Earth Rotation Service and published weekly along
// with predictions for a number of months into the future.  For the purpose
// of this routine the pblished and predicted values are kept in the ascii
// file 'iers.tab'.  UT1 is linked directly to the mean siderial time (GMST)
// through and equation as given in the Astronomical Almanac.  UT1 and,
// hence, GMST are obtained with earth rotation measures using VLBI and
// about 23 extragalactic radio sources as primary position standards.  The
// results are published as DUT1 = UT1-UTC and coordinates to the earths
// rotational pole with respect to its nominal position as defined by the
// IERS coordinate frame.
//               UT1 = UTC + DUT1 (from the 'iers.tab' file)
// The input parameters are mjd and utc (in fraction of a day), and the
// output arguments are ut1_mjd and ut1 (in fraction of a day).

void utc2ut1 (long mjd, double utc, long &ut1_mjd, double &ut1)
{
    ut1_mjd = mjd;
    ut1 = utc + get_ut1_offset(mjd, utc) / 86400.0;
    if (ut1 > 1.0) {
        ut1 -= 1.0;
        ut1_mjd++;
    }
    if (ut1 < 0.0) {
        ut1 += 1.0;
        ut1_mjd--;
    }
}

// UT0 is an observatory-specific version of UT1 in the sense that UT0
// contains the effect of polar motion on the observed rotation of the earth
// with respect to a non-rotating reference as defined by distand radio
// sources.  Normally UT1 is derived from UT0 by making the polar motion
// correction, but in the context of this set of routines we can infer our
// local UT0 from UTC and the published UT1-UTC offset by adding the polar
// offset in reverse.  From page 253 of the Explanatory Supplement to the
// Astronomical Almanac (1992)
//     UT0 = UT1 + tan(lat) * (x * sin(long) + y * cos(long))
// where x and y are the published pole offsets, and lat and long are the
// observatory's station coordinates.

void utc2ut0 (long mjd, double utc, long &ut0_mjd, double &ut0)
{
    long ut1_mjd = ut0_mjd = mjd;
    double ut1;
    utc2ut1(mjd, utc, ut1_mjd, ut1);
    
    double x, y;
    if (!get_pole_offsets(mjd, utc, x, y)) {
	ut0 = ut1;
	return;
    }
    x /= 1296000.0;		// convert arcseconds into days
    y /= 1296000.0;
    double lon = get_observatory_longitude();
    double lat = get_observatory_geoc_latitude();
    ut0 = ut1 + tan(lat) * (x * sin(lon) + y * cos(lon));
    ut0_mjd = ut1_mjd;
    if (ut0 < 0.0) {
	ut0 += 1.0;
	ut0_mjd--;
    }
    if (ut0 > 1.0) {
	ut0 -= 1.0;
	ut0_mjd++;
    }
}

// UT2 appears to be of mostly historical interest.  Before 1972 the time
// broadcast services kept their time signals within 0.1 seconds of UT2,
// which is UT1 with annual and semiannual variations in the earth's
// rotation removed.  The formal relation between UT1 and UT2 is
//         UT2 = UT1 + 0.022 * sin(2 * Pi * t) - 0.012 * cos(2 * Pi * t)
//                   - 0.006 * sin(4 * Pi * t) + 0.007 * cos(4 * Pi * t)
// where t = 2000.0 + (MJD - 51544.03) / 365.2422, the Besselian day
// fraction.  See the Explanatory Supplement to the Astronomical Almanac
// (1992) page 85 and the Explanatory Supplement to IERS Bulletins A and B,
// March 1995.

void utc2ut2 (long mjd, double utc, long &ut2_mjd, double &ut2)
{
    long ut1_mjd;
    double ut1;
    utc2ut1(mjd, utc, ut1_mjd, ut1);
    double arg = TWOPI * (2000.0 + ((double)ut1_mjd + ut1 - 51544.03))
								/ 365.2422;
    ut2 = ut1 + (0.022 * sin(arg) - 0.012 * cos(arg)
		 - 0.006 * sin(2.0 * arg) + 0.007 * cos(2.0 * arg)) / 86400.0;
    ut2_mjd = ut1_mjd;
    if (ut2 < 0.0) {
	ut2 += 1.0;
	ut2_mjd--;
    }
    if (ut2 > 1.0) {
	ut2 -= 1.0;
	ut2_mjd++;
    }
}

// Greenwich Mean Sidereal Time is the hour angle of the average position of
// the vernal equinox from the Greenwich meridian.  The "average" position
// of the vernal equixox is defined by the equation on page 50 of the
// Explanatory Supplement to the Astronomical Almanac (1992) linking UT1 to
// GMST.  This equation is coded in the SOFA routine iauGmst82().
//    GMST (in seconds at UT1=0) = 24110.54841 + 8640184.812866 * T
//				   + 0.093104 * T^2 - 0.0000062 * T^3
//    where T is in Julian centuries from 2000 Jan. 1 12h UT1
//	    T = d / 36525
//	    d = JD - 2451545.0
//  The ratio of the rates of UT1 and GMST is the ratio on the earth
// rotation period as measured with respect to the mean solar position and
// with respect to distant stars plus the average drift rate of the vernal
// equinox (precession).

double utc2gmst (long mjd, double utc)
{
    long ut1_mjd;
    double ut1;
    utc2ut1 (mjd, utc, ut1_mjd, ut1);
    return iauGmst82 (2400000.5, (double)ut1_mjd + ut1) / TWOPI;
}

// Greenwich Apparent Sidereal Time is GMST corrected for the shift in the
// position of the vernal equinox due to nutation.  Nutation is the
// mathematically predictable change in the direction of the earth rotation
// axis due to changing external torques from the sun, moon and planets.
// The smoothly varying part of the change in the earth's rotation
// (precession) is already accounted for in GMST.  The currently standard
// nutation theory is composed of 106 non-harmonically-related sine and
// cosine components plus 85 planetary correction terms.  The 4 dominant
// periods are 6798.4, 182.6, 13.7, and 3399.2 days with ampltudes of
// -17.1996, -1.3187, -0.2274, and 0.2062 arcseconds in ecliptic longitude,
// respectively.  The right ascension component of nutation is called the
// "equation of the equinoxes," which is coded in the SOFA routine
// iauEqeq94().
//	GAST = GMST + (equation of the equinoxes)

double utc2gast (long mjd, double utc)
{
    double gmst = utc2gmst(mjd, utc);
// ------- using the equation of the equinoxes -------
    // first get the Ephemeris Time (actually Terrestrial Dynamic Time)
    long tt_mjd;
    double tt;
    utc2tt (mjd, utc, tt_mjd, tt);
    // then compute and add the equation of the equinoxes
    return gmst + iauEqeq94(2400000.5, (double)tt_mjd + tt) / TWOPI;
// ---------------------------------------------------
/* a longer alternative using full nutation and taking the R.A. component
    double vector[3];
    iauS2c (gmst * TWOPI, 0.0, vector);

    long tt_mjd;
    double tt;
    utc2tt(mjd, utc, tt_mjd, tt);

    double nut_mat[3][3];
    iauNum00a(2400000.5, (double)tt_mjd + tt, nut_mat);

    double new_vector[3];
    iauRxp(nut_mat, vector, new_vector);

    double ra, dec;
    iauC2s(new_vector, &ra, &dec);
    ra /= TWOPI;
    if (ra < 0.0) ra += 1.0;
    if (ra >= 1.0) ra -= 1.0;
    return ra;
*/
}

// Local Mean Sidereal time is just GMST plus the observer's longitude
// measured positive to the east of Greenwich.

double utc2lmst (long mjd, double utc)
{
    double lmst = utc2gmst (mjd, utc) + get_observatory_longitude() / TWOPI;
    if (lmst < 0.0) {
        lmst += 1.0;
    }
    if (lmst >= 1.0) {
        lmst -= 1.0;
    }
    return lmst;
}

// Local Apparent Sidereal time is just GAST plus the observer's longitude
// measured positive to the east of Greenwich.

double utc2last (long mjd, double utc)
{
    double last = utc2gast (mjd, utc) + get_observatory_longitude() / TWOPI;
    if (last < 0.0) {
        last += 1.0;
    }
    if (last >= 1.0) {
        last -= 1.0;
    }
    return last;
}

// This routine saves the leap second data file path name and calls routines
// to close the old file, if open, open the named file, load the data
// into the leap second data structure, ls[], and close the file.  If
// successful, it returns the number of leap second entries, otherwise it
// returns zero.

int set_leap_sec_file_name (char *file_name)
{
    strcpy(leap_sec_file, file_name);
    int nn;
    if (open_leap_sec_file()) {
        nn = load_leap_sec_data();
	close_leap_sec_file();
	return nn;
    } else {
        return 0;
    }
}

// This routine opens the leap second ascii data file named in the global
// character string 'leap_sec_file'.  The file stream data structure pointer
// for use in fread() is stored in the global variable 'ls_stream'.  If
// successful, 1 is returned, otherwise 0 is returned.

int open_leap_sec_file()
{
    if (ls_stream != 0) {
	close_leap_sec_file();
    }
    if ((ls_stream = fopen(leap_sec_file, "r")) == NULL) {
        printf ("Connot open leap second tablefile '%s'\n", leap_sec_file);
        ls_stream = 0;
        return 0;
    }
    return 1;
}

// This routine closes the leap second data file, if it is open.

void close_leap_sec_file()
{
    if (ls_stream != 0) {
        fclose(ls_stream);
        ls_stream = 0;
    }
}

// This routine loads the leap second data structure, ls[], from the ascii
// leap second file.  If successful, it returns the number of leap second
// entries, otherwise it returns 0.  The ascii file format is
//    41498	1972 July 1    TAI-UTC = 11.0
//    41683	1973 January 1 TAI-UTC = 12.0
//    42048	1974 January 1 TAI-UTC = 13.0
//      :               :                 :
// The MJD must be the first field on each line, and the TAI-UTC value must
// be preceded by the '=' character.  The three year month day field are
// also expected.  Only the MJD and TAI-UTC values are actually used by the
// time routines, but the  year month day information is kept in the ls[]
// data structure.

int load_leap_sec_data()
{
    char line[100];
    int ln = 0;
    while ((fgets(line, 99, ls_stream) != NULL) && (ln < MAX_NUM_LEAP_SEC)) {
        // 50083.0	1996 January 1 TAI-UTC = 30.0
        if (sscanf(line, "%ld %s %s %s",
                   &ls[ln].mjd, &ls[ln].yr, &ls[ln].month, &ls[ln].day) < 4) {
            printf ("Error reading dates in leap second file '%s', line %d\n",
                    leap_sec_file, ln + 1);
            return 0;
        }
        // skip to the '='
	int i;
        for (i = 0; i < strlen(line); i++) {
            if (line[i] == '=') {
                break;
            }
        }
        if (sscanf(&line[i+1], "%lf", &ls[ln].utc_offset) != 1) {
            printf (
              "Error reading utc_offset in leap second file '%s', line %d\n",
                    leap_sec_file, ln + 1);
            return 0;
        }
        ln++;
    }
    if (ln == MAX_NUM_LEAP_SEC) {
        printf ("May have exceeded the leap second file length limit\n");
        printf (" of %d lines.  Recompile 'astrtime.cc'n", MAX_NUM_LEAP_SEC);
    }
    num_leap_sec = ln;
    return ln;
}

// This routine returns the TAI-UTC offset in effect for the MJD/UTC
// specified.  If the leap second data structure has not already been
// filled, the routines are called to open the file and load the leap second
// data structure, ls[].

double get_leap_sec(long mjd, double utc)
{
    if (num_leap_sec == 0) {
        if (open_leap_sec_file()) {
            if(load_leap_sec_data() == 0) {
		close_leap_sec_file();
                return 0.0;
            }
        } else {
            return 0.0;
        }
	close_leap_sec_file();
    }
    // move to the right integer MJD in case UTC is not between 0.0 and 1.0
    int safety = 1000;
    while ((utc >= 1.0) && (safety-- > 0)) {
        utc -= 1.0;
        mjd++;
    }
    while ((utc < 1.0) && (safety-- > 0)) {
        utc += 1.0;
        mjd--;
    }
    if (safety < 1) {
        printf ("get_leap_sec: utc must be > -1000 and < 1000\n");
        return 0.0;
    }
    if (mjd < ls[0].mjd) {
        printf (
      "get_leap_sec: requested mjd is before beginning of leap second file\n");
        return ls[0].utc_offset - 1.0;
    }
    // search backwards through the data until finding the date preceding
    // the specified MJD/UTC.  Leap seconds are added at 0 hours of the MJD
    // in the table.
    int i;
    for (i = num_leap_sec - 1; i >= 0; i--) {
        if (mjd >= ls[i].mjd) {
            return ls[i].utc_offset;
        }
    }
    // should never reach this point
    printf ("get_leap_sec: Error!\n");
    return 0.0;
}
// This routine saves the IERS data file path name and calls routines
// to close the old file, if open, open the named file, load the data
// into the IERS data structure, ls[], and close the file.  If
// successful, it returns the number of IERS entries, otherwise it
// returns zero.

int set_iers_file_name (char *file_name)
{
    strcpy(iers_file, file_name);
    close_iers_file();
    int nn;
    if (open_iers_file()) {
        nn = load_iers_data();
	close_iers_file();
	return nn;
    } else {
        return 0;
    }
}

// This routine opens the IERS ascii data file named in the global
// character string 'leap_sec_file'.  The file stream data structure pointer
// for use in fread() is stored in the global variable 'iers_stream'.  If
// successful, 1 is returned, otherwise 0 is returned.

int open_iers_file()
{
    if (iers_stream != 0) {
	close_iers_file();
    }
    if ((iers_stream = fopen(iers_file, "r")) == NULL) {
        printf ("Connot open IERS tablefile '%s'\n", iers_file);
        iers_stream = 0;
        return 0;
    }
    return 1;
}

// This routine closes the IERS data file, if it is open.

void close_iers_file()
{
    if (iers_stream != 0) {
        fclose(iers_stream);
        iers_stream = 0;
    }
}

// This routine loads the IERS data structure, iers[], from the ascii
// IERS file.  If successful, it returns the number of entries, otherwise it
// returns 0.  The ascii file format is
// 50070 -0.1566  0.1483 -0.40088
// 50075 -0.1709  0.1621 -0.41367
// 50080 -0.1841  0.1770 -0.42633
// 50083 -0.1882  0.1839  0.56874 <<<<<< leap second
// 50085 -0.1961  0.1929  0.56115
// 50090 -0.2068  0.2096  0.54873
//      :               :                 :
// The MJD must be the first field on each line, and the x and y pole
// offsets are next, and the fourth column is UT1-UTC.  Comments after the
// fourth column are ignored.  Note that there must be an entry for the MJD
// of a leap second to allow a proper transition across this discontinuity.
// All values are for 0 hours of the MJD given.  The entry interval shown is
// 5 days, but this may be shorter or longer as required by the accuracy of
// a linear interpolation between two entries.

int load_iers_data()
{
    char line[100];
    int ln = 0;
    while ((fgets(line, 99, iers_stream) != NULL) &&
                                          (ln < MAX_NUM_IERS_OFFSETS)) {
        //  MJD      x       y     UT1 - UTC
        // 50060  -.12288  .12426  -.389102
        if (sscanf(line,"%ld %lf %lf %lf",
                   &iers[ln].mjd, &iers[ln].x, &iers[ln].y,
                   &iers[ln].ut1_offset) < 4) {
            printf ("Error reading dates in IERS file '%s', line %d\n",
                    iers_file, ln + 1);
            return 0;
        }
        ln++;
    }
    if (ln == MAX_NUM_IERS_OFFSETS) {
        printf ("May have exceeded the IERS file length limit\n");
        printf (" of %d lines.  Recompile 'astrtime.cc'n",
                MAX_NUM_IERS_OFFSETS);
    }
    num_iers_offsets = ln;
    return ln;
}

// This routine returns the UT1-UTC offset in effect for the MJD/UTC
// specified.  If the IERS data structure has not already been
// filled, the routines are called to open the file and load the IERS
// data structure, iers[].

double get_ut1_offset(long mjd, double utc)
{
    if (num_iers_offsets == 0) {
        if (open_iers_file()) {
            if(load_iers_data() == 0) {
		close_iers_file();
                return 0.0;
            }
        } else {
            return 0.0;
        }
	close_iers_file();
    }
    // move to the right integer MJD in case UTC is not between 0.0 and 1.0
    int safety = 1000;
    while ((utc >= 1.0) && (safety-- > 0)) {
        utc -= 1.0;
        mjd++;
    }
    while ((utc < 1.0) && (safety-- > 0)) {
        utc += 1.0;
        mjd--;
    }
    if (safety < 1) {
        printf ("get_ut1_offset: utc must be > -1000 and < 1000\n");
        return 0.0;
    }
    if (mjd < iers[0].mjd) {
        printf (
         "get_ut1_offset: requested mjd is before beginning of IERS file %ld\n",
          iers[0].mjd);
        return iers[0].ut1_offset;
    }
    if (mjd > iers[num_iers_offsets - 1].mjd) {
        printf (
      "get_ut1_offset: requested mjd is beyond end of IERS file\n");
        return iers[num_iers_offsets - 1].ut1_offset;
    }
    // search backwards through the data until finding the date preceding
    // the specified MJD/UTC.
    int i;
    for (i = num_iers_offsets - 2; i >= 0; i--) {
        if (mjd >= iers[i].mjd) {
	    // interpolate to the specified MJD/UTC
            double interval = iers[i+1].mjd - iers[i].mjd;
            double change = iers[i+1].ut1_offset - iers[i].ut1_offset;
            // watch for a leapsecond jump
            if (change > 0.3) {
                change -= 1.0;
            } else if (change < -0.3) {
                change += 1.0;
            }
            double interp = ((double)(mjd - iers[i].mjd) + utc) / interval;
            return iers[i].ut1_offset + interp * change;
        }
    }
    // should never reach this point
    printf ("get_ut1_offset: Error!\n");
    return 0.0;
}

// This routine returns the pole offsets in effect for the MJD/UTC
// specified.  If the IERS data structure has not already been
// filled, the routines are called to open the file and load the IERS
// data structure, iers[].  The pole offsets in x and y are in arcseconds.
// x is in the direction of the Greenwich meridian, and y is toward 90
// degrees west longitude.

int get_pole_offsets(long mjd, double utc, double &x, double &y)
{
    x = y = 0;
    if (num_iers_offsets == 0) {
        if (open_iers_file()) {
            if(load_iers_data() == 0) {
		close_iers_file();
                return 0;
            }
        } else {
            return 0;
        }
	close_iers_file();
    }
    // move to the right integer MJD in case UTC is not between 0.0 and 1.0
    int safety = 1000;
    while ((utc >= 1.0) && (safety-- > 0)) {
        utc -= 1.0;
        mjd++;
    }
    while ((utc < 1.0) && (safety-- > 0)) {
        utc += 1.0;
        mjd--;
    }
    if (safety < 1) {
        printf ("get_ut1_offset: utc must be > -1000 and < 1000\n");
        return 0;
    }
    if (mjd < iers[0].mjd) {
        printf (
         "get_ut1_offset: requested mjd is before beginning of IERS file %ld\n",
          iers[0].mjd);
        return 0;
    }
    if (mjd > iers[num_iers_offsets - 1].mjd) {
        printf (
      "get_ut1_offset: requested mjd is beyond end of IERS file\n");
        return 0;
    }
    // search backwards through the data until finding the date preceding
    // the specified MJD/UTC.
    int i;
    for (i = num_iers_offsets - 2; i >= 0; i--) {
        if (mjd >= iers[i].mjd) {
	    // interpolate to the specified MJD/UTC
            double interval = iers[i+1].mjd - iers[i].mjd;
            double change_x = iers[i+1].x - iers[i].x;
            double change_y = iers[i+1].y - iers[i].y;
            double interp = ((double)(mjd - iers[i].mjd) + utc) / interval;
	    x = iers[i].x + interp * change_x;
	    y = iers[i].y + interp * change_y;
            return 1;
        }
    }
    // should never reach this point
    printf ("get_pole_offsets: Error!\n");
    return 0;
}

// This routine converts the observatory's specified geodetic longitude,
// latitude, and height above the reference geoid into geocentric
// rectangular coordinates and saves them for topocentric calculations.
// The height above the reference geoid is generally withing 100 meters or
// so of the more commonly published height above mean sea level.  Geodetic
// longitude and latitude (measured to the normal to the reference geoid)
// are generally within a minute of arc of the geographic or astronomical
// longitude and latitude (measured to the local gravity vector).  Longitude
// is positive east and it and latitude are in degrees.  The height in
// meters.
// In most cases, specifying the observatory rectangular coordinates
// directly with 'set_observatory_xyz()' will be more accurate.

void set_latlong2xyz (double geodetic_long, double geodetic_lat,
                      double geoid_height)
{
    double lon = geodetic_long * DEG2RAD;
    double lat = geodetic_lat * DEG2RAD;
    double radius, z;
    ptwGeoc(lat, geoid_height, &radius, &z);
    radius *= AU2METERS;
    obs.x = radius * cos(lon);
    obs.y = radius * sin(lon);
    obs.z = z * AU2METERS;
}

// This routine saves the observatory's geocentric rectangular coordinates
// for topocentric calculations.  The x axis is toward the equator/Greenwich
// meridian intersection, y is toward the equator 90 degrees east, and z is
// toward the north pole.  Coordinates must be in meters.

void set_observatory_xyz (double x, double y, double z)
{
    obs.x = x;
    obs.y = y;
    obs.z = z;
}

// This routine converts observatory rectangular coordinates to longitude.
// Geocentric and geodetic longitude are the same as long a an
// earth-centered geoid is assumed.

double get_observatory_longitude()
{
    if ((obs.y == 0.0) && (obs.x == 0.0)) {
        return 0.0;
    }
    return atan2(obs.y, obs.x);
}

// This routine converts the obsevatory's rectangular coordinates into
// geocentric latitude.  The more common latitudes in use or geodetic or
// geographic, but conversion to geodetic is quite complex, and conversion
// to geographic is local gravity dependent.  The difference between
// geocentric and geodetic latitude can be as much as about 10 arcminutes at
// mid latitudes.

double get_observatory_geoc_latitude()
{
    double xy = sqrt(obs.x * obs.x + obs.y * obs.y);
    return atan(obs.z / xy);
}

// This routine returns the observatory's distance to the center of the
// earth, in meters, from the geocentric rectangular coordinates.

double get_observatory_radius()
{
    return sqrt(obs.x * obs.x + obs.y * obs.y + obs.z * obs.z);
}

double get_observatory_x()
{
    return obs.x;
}

double get_observatory_y()
{
    return obs.y;
}

double get_observatory_z()
{
    return obs.z;
}

// This routine returns the observatory's geocentric rectangular coordinate
// position, in km, and velocity, in km/day, with respect to the geocentric
// J2000 equator and equinox in a manner compatible with the JPL planetary
// ephemerides routines.  The observatory's coordinates are uncorrected for
// polar motion (latitude/longitude error) and offset of the Celestial
// Ephemeris Pole (nutation error), which could amount to an error of as
// much as 25 meters or so (80 ns in delay).  Note that the time is
// Terrestrial Dynamic Time.

int get_observatory_posn_vel(double pv[6], double jd, double day, int vel_flag)
{
    double obs_vector[3];
    obs_vector[0] = obs.x;
    obs_vector[1] = obs.y;
    obs_vector[2] = obs.z;
    // ------------------------------------------------
    // This is where to correct for polar motion offset
    // in the observatory's coordinates.
    // ------------------------------------------------
    // with respect to the current rotation axis:
    long tt_mjd = (long)(jd + day - 2400000.5);
    double tt = (jd + day - 2400000.5) - (double)tt_mjd;
    long mjd;
    double utc;
    tt2utc(tt_mjd, tt, mjd, utc);
    double arg = utc2last(mjd, utc) * TWOPI;
    double s = sin(arg);
    double c = cos(arg);
    double xy = sqrt(obs_vector[0] * obs_vector[0] +
		     obs_vector[1] * obs_vector[1]);
    pv[0] = xy * c;		// position
    pv[1] = xy * s;
    pv[2] = obs_vector[2];
    if (vel_flag) {
	pv[3] = -xy * s * 366.25 / 365.25;	// velocity in meters
	pv[4] = xy * c * 366.25 / 365.25;	// per UTC day
	pv[5] = 0.0;
    } else {
	int i;
	for (i = 3; i < 6; i++) pv[i] = 0.0;
    }
    // ----------------------------------------------------------
    // This is where to correct for offset of the Celestial
    // Ephemeris pole.  This offset is the bit of unpredictable
    // nutation not in the theory of nutation used by iauPnm00a()
    // ----------------------------------------------------------
    // then rotate into the J2000 frame
    double pn_mat[3][3];
    iauPnm00a (2400000.5, jd + day - 2400000.5, pn_mat);

    double vector[3], new_vector[3];
    int i;
    for (i = 0; i < 3; i++) vector[i] = pv[i];
    iauTrxp(pn_mat, vector, new_vector);
    for (i = 0; i < 3; i++) pv[i] = new_vector[i] / 1000.0;

    if (vel_flag) {
	for (i = 0; i < 3; i++) vector[i] = pv[i+3];
	iauTrxp(pn_mat, vector, new_vector);
	// convert to km/day
	for (i = 0; i < 3; i++) pv[i+3] = new_vector[i] * TWOPI / 1000.0;
    }
    return 1;
}
