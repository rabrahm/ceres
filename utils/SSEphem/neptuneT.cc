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

// JPL DExxx Neptune's position test main program

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "jpl_eph.h"
#include "sofa.h"
#include "ptw_aux.h"

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
	sprintf(str, "%3d %2d %4.1f", d, min, deg);
    } else {
	sprintf(str, "-%2d %2d %4.1f", d, min, deg);
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
    sprintf(str, "%2d %2d %5.2f", sign * d, min, deg);
    return str;
}

int main (int argc, char *argv[])
{
    if (argc < 2) {
        printf ("Usage: neptuneT ephemeris_file_name\n");
        return 1;
    } else {
        printf ("\nJPL Neptune position test.  Compare with\n");
	printf ("1996 Astronomical Almanac pages E13 and E38-41\n\n");
    }
    if (!open_ephemeris(argv[1])) {
        return 1;
    }
    double jd = 2450082.5;    // Jan 0, 1996
    double tptr[6];

    // print position for 0h TT on Feb 7, Apr 27, Aug 25, and Nov 13, 1995
    static double dx[4] = { 38.0, 118.0, 238.0, 318.0 };
    static char *gdate[4] = { "Feb 7, 1995", "Apr 27, 1995", "Aug 25, 1995",
			      "Nov 13, 1995" };
    int ix;
    for (ix = 0; ix < 4; ix++) {
	double day = dx[ix];
	// start with heliocentric positions to be compared with
	// Astronomical Almanac pages E13
        if (get_planet_posn_vel(tptr, Neptune, Heliocentric, jd, day)) {
	    double d, lon, lat, xra, xdec;
	    d = sqrt(tptr[0] * tptr[0] + tptr[1] * tptr[1] + tptr[2] * tptr[2]);
	    // compute J2000 coordinates
            double posn[3];
	    int jxx;
            for (jxx = 0; jxx < 3; jxx++) posn[jxx] = tptr[jxx];
	    iauC2s (posn, &xra, &xdec);

	    // convert to ecliptic coordinates of date
	    double date = jd + day - 2400000.5;
	    ptwEqecl(xra, xdec, date, &lon, &lat);
	    lat *= RAD2DEG;
	    lon *= RAD2DEG;
	    if (lon < 0.0) lon += 360.0;
	    printf ("TT JD %6.3f  %s\n", jd + day, gdate[ix]);
	    printf ("  Longitude %s  Latitude %s, Dist. %10.8f\n",
		    deg2dms_str(lon), deg2dms_str(lat), d / get_au());

	    // now get geocentric positions to be compared with
	    // Astronomical Almanac pages E38-41
	    //  first get the light travel time
	    if (!get_planet_posn_vel(tptr, Neptune, Geocentric, jd, day))
		return 1;
	    d = sqrt(tptr[0] * tptr[0] + tptr[1] * tptr[1] + tptr[2] * tptr[2]);

	    // then recompute for light travel advanced time
	    // since the earth's position is computed for the earlier time,
	    // too, this accounts for aberration.
	    double day_light = day - d / (C_LIGHT * 86400.0);
	    if (!get_planet_posn_vel(tptr, Neptune, Geocentric, jd, day_light))
		return 1;
	    // add precession and nutation
	    double rmatpn[3][3], new_vec[3], ra, dec;
	    iauPnm00a (2400000.5, date, rmatpn);
            for (jxx = 0; jxx < 3; jxx++) posn[jxx] = tptr[jxx];
	    iauRxp (rmatpn, posn, new_vec);
	    iauC2s (new_vec, &ra, &dec);
	    dec *= RAD2DEG;
	    ra *= RAD2DEG;
	    if (ra < 0.0) ra += 360.0;
	    printf ("  Appar  RA  %s      Dec %s, Dist. %10.8f\n\n",
		    deg2hms_str(ra), deg2dms_str(dec), d / get_au());
        } else {
            printf ("Neptune position calculation failed\n");
        }
    }
    return 0;
}
