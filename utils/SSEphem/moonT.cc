//This file is part of SSEphem.

//  SSEphem is a collection of C/C++ source files that provide an interface
//  to the JPL Solar System Ephemeris.
//  Copyright (C) 2009 Associated Universities, Inc. Washington DC, USA

//  SSEphem is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  SSEphem is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.

//  You should have received a copy of the GNU Lesser General Public License
//  along with SSEphem.  If not, see <http://www.gnu.org/licenses/>.

// JPL DExxx Moon's position test main program

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "jpl_eph.h"
#include "sofa.h"
#include "astrtime.h"

char *deg2dms_str(double deg)
{
    char *str = new char[50];

    int sign = 1;
    if (deg < 0.0) sign = -1;
    deg = fabs(deg);
    int d = (int)deg;
    deg = 60.0 * (deg - (double)d);
    int min = (int)deg;
    deg = 60.0 * (deg - (double)min);
    if (sign > 0) {
	sprintf(str, "%3d %2d %5.2f", d, min, deg);
    } else {
	sprintf(str, "-%2d %2d %5.2f", d, min, deg);
    }
    return str;
}

char *deg2hms_str(double deg)
{
    char *str = new char[50];

    int sign = 1;
    if (deg < 0.0) sign = -1;
    deg = fabs(deg) / 15.0;
    int d = (int)deg;
    deg = 60.0 * (deg - (double)d);
    int min = (int)deg;
    deg = 60.0 * (deg - (double)min);
    sprintf(str, "%2d %2d %6.3f", sign * d, min, deg);
    return str;
}

int main (int argc, char *argv[])
{
    if (argc < 2) {
        printf ("Usage: moonT ephemeris_file_name\n");
        return 1;
    } else {
	printf ("\nJPL Moon position test.  Compare with\n");
	printf ("1996 Astronomical Almanac even # pages D6-20\n\n");
    }
    if (!open_ephemeris(argv[1])) {
        return 1;
    }
    double jd = 2450082.5;    // Jan 0, 1996
    double tptr[6];

    // print position for 0h TDB on Jan 1, Apr 10, July 19, and Oct. 27, 1996
    static char gdate[4][20] = { "Jan 1, 1996", "Apr 10, 1996",
				 "July 19, 1996", "Oct. 27, 1996" };
    int i = 0;
    double day;
    for (day = 1.0; day < 365; day += 100.0) {
        if (get_moon_posn_vel(tptr, Geocentric, jd, day)) {
	    double d, dec, ra;
	    d = sqrt(tptr[0] * tptr[0] + tptr[1] * tptr[1] + tptr[2] * tptr[2]);
	    // compute and print J2000 coordinates
            double posn[3];
	    int jxx;
            for (jxx = 0; jxx < 3; jxx++) posn[jxx] = tptr[jxx];
	    iauC2s (posn, &ra, &dec);
	    dec *= RAD2DEG;
	    ra *= RAD2DEG;
	    if (ra < 0.0) ra += 360.0;

	    // This position assumes no light travel time so is not
	    // one that could be used to point a telescope
	    printf ("TDB JD %6.3f   %s\n", jd + day, gdate[i]);
	    i++;
	    printf ("  J2000 RA %s  Dec %s, Dist. %10.8f au\n",
		    deg2hms_str(ra), deg2dms_str(dec), d / get_au());

	    // compute and print apparent coordinates
	    // first recompute for light travel advanced time
	    // since the earth's position is computed for the earlier time,
	    // too, this accounts for aberration.
	    double day_light = day - d / (C_LIGHT * 86400.0);
	    if (!get_moon_posn_vel(tptr, Geocentric, jd, day_light))
		return 1;
	    // add precession and nutation
	    double rmatpn[3][3], new_vec[3];
	    double date = jd + day - 2400000.5;
	    iauPnm00a (2400000.5, date, rmatpn);
            for (jxx = 0; jxx < 3; jxx++) posn[jxx] = tptr[jxx];
	    iauRxp (rmatpn, posn, new_vec);
	    double ara, adec;
	    iauC2s (new_vec, &ara, &adec);
	    adec *= RAD2DEG;
	    ara *= RAD2DEG;
	    if (ara < 0.0) ara += 360.0;

	    // compare this with Astron. Almanac even # pages D6-20
	    printf ("  Appar RA %s  Dec %s\n",
		    deg2hms_str(ara), deg2dms_str(adec));

	    // Now do topocentric calculation
	    if (!get_moon_posn_vel(tptr, Topocentric, jd, day))
		return 1;
	    double td = sqrt(tptr[0] * tptr[0] + tptr[1] * tptr[1]
			     + tptr[2] * tptr[2]);
	    // recompute for light travel advanced time
	    // since the earth's position is computed for the earlier time,
	    // too, this accounts for aberration.
//  Lines deleted for comparison with geoentric + topo correction below
//	    day_light = day - d / (C_LIGHT * 86400.0);
//	    if (!get_moon_posn_vel(tptr, Topocentric, jd, day_light))
//		return 1;
	    // add precession and nutation
	    rmatpn[3][3], new_vec[3];
	    date = jd + day - 2400000.5;
	    iauPnm00a (2400000.5, date, rmatpn);
            for (jxx = 0; jxx < 3; jxx++) posn[jxx] = tptr[jxx];
	    iauRxp (rmatpn, posn, new_vec);
	    double tra, tdec;
	    iauC2s (new_vec, &tra, &tdec);
	    tdec *= RAD2DEG;
	    tra *= RAD2DEG;
	    if (tra < 0.0) tra += 360.0;
	    printf ("  Topoc RA %s  Dec %s  Dist. %8.2f km\n",
		    deg2hms_str(tra), deg2dms_str(tdec), td);

	    // compare this with topocentric correction equations given
	    // an page D3 of the Astronomical Almanac.
	    if (!get_moon_posn_vel(tptr, Geocentric, jd, day))
		return 1;
	    d = sqrt(tptr[0] * tptr[0] + tptr[1] * tptr[1]
			     + tptr[2] * tptr[2]);
	    // add precession and nutation
	    date = jd + day - 2400000.5;
	    iauPnm00a (2400000.5, date, rmatpn);
            for (jxx = 0; jxx < 3; jxx++) posn[jxx] = tptr[jxx];
	    iauRxp (rmatpn, posn, new_vec);
	    double gra, gdec;
	    iauC2s (new_vec, &gra, &gdec);
	    if (gra < 0.0) gra += TWOPI;

	    double rho = get_observatory_radius() / 1000.0;
	    double lat = get_observatory_geoc_latitude();
	    long tdb_mjd = (long)(jd + day - 2400000.5);
	    double tdb = (jd + day - 2400000.5) - (double)tdb_mjd;
	    long mjd;
	    double utc;
	    tdb2utc(tdb_mjd, tdb, mjd, utc);
	    double last = utc2last(mjd, utc);
	    double xp = d * cos(gdec) * cos(gra)
		                - rho * cos(lat) * cos(TWOPI * last);
	    double yp = d * cos(gdec) * sin(gra)
				- rho * cos(lat) * sin(TWOPI * last);
	    double zp = d * sin(gdec) - rho * sin(lat);
	    double dp = sqrt(xp * xp + yp * yp + zp * zp);

	    double rap = atan2(yp, xp) * RAD2DEG;
	    if (rap < 0.0) rap += 360.0;
	    double decp = asin(zp / dp) * RAD2DEG;
	    printf ("From Almanac topocentric correction:\n");
	    printf ("  Topoc RA %s  Dec %s  Dist. %8.2f km\n\n",
		    deg2hms_str(rap), deg2dms_str(decp), dp);

//	    printf (
//	      "x,y,z    %10.7f %10.7f %10.7f\nxv,yv,zv %10.7f %10.7f %10.7f\n",
//		   tptr[0] / get_au(), tptr[1] / get_au(), tptr[2] / get_au(),
//		   tptr[3] / get_au(), tptr[4] / get_au(), tptr[5] / get_au());
        } else {
            printf ("Moon position calculation failed\n");
        }
    }
    return 0;
}
