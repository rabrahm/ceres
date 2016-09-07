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

// This routine returns the rectangular coordinate components of the
// Green Bank 140-ft telescope for the specified date and time.  The
// reference frame is defined by the solar system ephemeris used.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "jpl_eph.h"
#include "astrtime.h"
#include "sofa.h"

int main (int argc, char *argv[])
{
    if (argc < 4) {
        printf ("Usage: gbxyz ephemeris_file_name MJD UTC(day fraction)\n");
        return 1;
    }
    if (!open_ephemeris(argv[1])) {
        return 1;
    }
    long mjd = atol(argv[2]);	// Specified Modified Julian Date
    double utc = atof(argv[3]);	// Fraction of a day in UTC

    double pv[6];	// Green Bank's barycentric position and velocity
    double pvt[6];	// Green Bank's geocentric position and velocity
    long tt_mjd;	// Date for Terrestrial Dynamic Time
    double tt;		// Fraction of TDT day
    long tdb_mjd;	// Date for Barycentric Dynamic Time
    double tdb;		// Fraction of TDB day

    utc2tt(mjd, utc, tt_mjd, tt);
    double jd = (double)tt_mjd + 2400000.5;	// TDT Julian date
    utc2tdb(mjd, utc, tdb_mjd, tt);
    double jd_tdb = (double)tdb_mjd + 2400000.5;	// TDB Julian date

    int status = get_earth_posn_vel(pv, SS_Barycentric, jd_tdb, tdb, 0);
    get_observatory_posn_vel(pvt, jd, tt, 0);
    int i;
    for (i = 0; i < 6; i++) {
	pvt[i] += pv[i];
    }
    printf ("Green Bank 140-ft Barycentric position:\n");
    if (status) {
	printf ("   X = %15.3f  kilometers\n", pv[0]);
	printf ("   Y = %15.3f\n", pv[1]);
	printf ("   Z = %15.3f\n", pv[2]);
    } else {
	printf ("   Error getting earth's position\n");
    }
    close_ephemeris();
    return 0;
}
