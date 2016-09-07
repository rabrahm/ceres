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

// This routine returns the geometric speed-of-light delay, in seconds,
// between two solar system reference frames for a given J2000 right
// ascension and declination and time.  The reference frame selections are
// given by the RefFrame enum in jpl_eph.h, and ra2000 and dec2000 are in
// radians.  The delay is positive when ref2 is farther from the source at
// (ra, dec) than is ref1.

#include <stdio.h>
#include "jpl_eph.h"
#include "astrtime.h"
#include "delay.h"
#include "sofa.h"

double get_geometric_delay (RefFrame ref1, RefFrame ref2,
			    long mjd, double utc,
			    double ra2000, double dec2000)
{
    int i;
    double pv1[6], pv2[6], pv[6];
    long tdb_mjd, tt_mjd;
    double tdb, tt;

    utc2tdb(mjd, utc, tdb_mjd, tdb);
    utc2tt(mjd, utc, tt_mjd, tt);
    double jd = (double)tdb_mjd + 2400000.5;
    double jd_tt = (double)tt_mjd + 2400000.5;
    int status;
    switch (ref1) {
	case Heliocentric:
	    status = get_posn_vel(pv1, Sun, jd, tdb, 0);
	    break;
	case SS_Barycentric:
	    for (i = 0; i < 6; i++) pv1[i] = 0.0;
	    break;
	case EM_Barycentric:
	    status = get_posn_vel(pv1, EM_Bary, jd, tdb, 0);
	    break;
	case Geocentric:
	    status = get_earth_posn_vel(pv1, SS_Barycentric, jd, tdb, 0);
	    break;
	case Topocentric:
	    status = get_earth_posn_vel(pv, SS_Barycentric, jd, tdb, 0);
	    get_observatory_posn_vel(pv1, jd_tt, tt, 0);
	    for (i = 0; i < 6; i++) pv1[i] += pv[i];
	    break;
	default:
	    printf(
	      "get_geometric_delay: First reference frame not recognized\n");
	    break;
    }
    if (!status) return 0.0;
    switch (ref2) {
	case Heliocentric:
	    status = get_posn_vel(pv2, Sun, jd, tdb, 0);
	    break;
	case SS_Barycentric:
	    for (i = 0; i < 6; i++) pv2[i] = 0.0;
	    break;
	case EM_Barycentric:
	    status = get_posn_vel(pv2, EM_Bary, jd, tdb, 0);
	    break;
	case Geocentric:
	    status = get_earth_posn_vel(pv2, SS_Barycentric, jd, tdb, 0);
	    break;
	case Topocentric:
	    status = get_earth_posn_vel(pv, SS_Barycentric, jd, tdb, 0);
	    get_observatory_posn_vel(pv2, jd_tt, tt, 0);
	    for (i = 0; i < 6; i++) pv2[i] += pv[i];
	    break;
	default:
	    printf(
	      "get_geometric_delay: Second reference frame not recognized\n");
	    break;
    }
    if (!status) return 0.0;
    double p[3];
    for (i = 0; i < 3; i++) {
	p[i] = pv1[i] - pv2[i];
    }
    double vector[3];
    iauS2c(ra2000, dec2000, vector);
    return iauPdp(vector, p) / C_LIGHT;
}
