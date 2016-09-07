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

// JPL DExxx closure tests through various paths around solar system

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "jpl_eph.h"

int main (int argc, char *argv[])
{
    if (argc < 2) {
        printf ("Usage: closureT ephemeris_file_name\n");
        return 1;
    } else {
	printf ("\nJPL Solar System Ephemeris closure tests\n\n");
    }
    if (!open_ephemeris(argv[1])) {
        return 1;
    }
    int i, flag;
    double jd = 2450085.5;    // Jan 3, 1996
    double day = 0.3;
    double pv1[6], pv2[6], pv3[6], d1, d2, d3, d;
    
    if (!get_earth_posn_vel(pv1, Topocentric, jd, day)) {
        return 1;
    }
    d1 = sqrt(pv1[0] * pv1[0] + pv1[1] * pv1[1] + pv1[2] * pv1[2]);
    if (!get_sun_posn_vel(pv2, Geocentric, jd, day)) {
        return 1;
    }
    d2 = sqrt(pv2[0] * pv2[0] + pv2[1] * pv2[1] + pv2[2] * pv2[2]);
    if (!get_sun_posn_vel(pv3, Topocentric, jd, day)) {
        return 1;
    }
    d3 = sqrt(pv3[0] * pv3[0] + pv3[1] * pv3[1] + pv3[2] * pv3[2]);
    flag = 0;
    for (i = 0; i < 6; i++) {
	if (fabs(pv1[i] + pv2[i] - pv3[i]) > 1.0e-6) {
	    flag++;
	}
    }
    if ((d1 < 1.0e-5) || (d2 < 1.0e-10) || (d3 < 1.0e-10)) {
	flag++;
    }
    d = (d1 + d2 + d3) / get_au();
    printf ("Observatory - Geocenter - Sun    (%7.4f au) ", d);
    if (flag) {
	printf ("FAILED\n");
    } else {
	printf ("Passed\n");
    }
    
    if (! get_planet_posn_vel(pv1, Mercury, Topocentric, jd, day)) {
        return 1;
    }
    d1 = sqrt(pv1[0] * pv1[0] + pv1[1] * pv1[1] + pv1[2] * pv1[2]);
    if (!get_planet_posn_vel(pv2, Mercury, Heliocentric, jd, day)) {
        return 1;
    }
    d2 = sqrt(pv2[0] * pv2[0] + pv2[1] * pv2[1] + pv2[2] * pv2[2]);
    if (!get_sun_posn_vel(pv3, Topocentric, jd, day)) {
        return 1;
    }
    d3 = sqrt(pv3[0] * pv3[0] + pv3[1] * pv3[1] + pv3[2] * pv3[2]);
    flag = 0;
    for (i = 0; i < 6; i++) {
	if (fabs(pv1[i] - pv2[i] - pv3[i]) > 1.0e-6) {
	    printf ("diff %d = %e\n", i, fabs(pv1[i] - pv2[i] - pv3[i]));
	    flag++;
	}
    }
    if ((d1 < 1.0e-5) || (d2 < 1.0e-10) || (d3 < 1.0e-10)) {
	flag++;
    }
    d = (d1 + d2 + d3) / get_au();
    printf ("Observatory - Mercury - Sun      (%7.4f au) ", d);
    if (flag) {
	printf ("FAILED\n");
    } else {
	printf ("Passed\n");
    }
    
    if (!get_earth_posn_vel(pv1, Heliocentric, jd, day)) {
        return 1;
    }
    d1 = sqrt(pv1[0] * pv1[0] + pv1[1] * pv1[1] + pv1[2] * pv1[2]);
    if (!get_moon_posn_vel(pv2, Geocentric, jd, day)) {
        return 1;
    }
    d2 = sqrt(pv2[0] * pv2[0] + pv2[1] * pv2[1] + pv2[2] * pv2[2]);
    if (!get_moon_posn_vel(pv3, Heliocentric, jd, day)) {
        return 1;
    }
    d3 = sqrt(pv3[0] * pv3[0] + pv3[1] * pv3[1] + pv3[2] * pv3[2]);
    flag = 0;
    for (i = 0; i < 6; i++) {
	if (fabs(pv1[i] + pv2[i] - pv3[i]) > 1.0e-6) {
	    printf ("diff %d = %e\n", i, fabs(pv1[i] + pv2[i] - pv3[i]));
	    flag++;
	}
    }
    if ((d1 < 1.0e-5) || (d2 < 1.0e-10) || (d3 < 1.0e-10)) {
	flag++;
    }
    d = (d1 + d2 + d3) / get_au();
    printf ("Geocenter - Moon - Sun           (%7.4f au) ", d);
    if (flag) {
	printf ("FAILED\n");
    } else {
	printf ("Passed\n");
    }
    
    if (!get_earth_posn_vel(pv1, EM_Barycentric, jd, day)) {
        return 1;
    }
    d1 = sqrt(pv1[0] * pv1[0] + pv1[1] * pv1[1] + pv1[2] * pv1[2]);
    if (!get_moon_posn_vel(pv2, Geocentric, jd, day)) {
        return 1;
    }
    d2 = sqrt(pv2[0] * pv2[0] + pv2[1] * pv2[1] + pv2[2] * pv2[2]);
    if (!get_moon_posn_vel(pv3, EM_Barycentric, jd, day)) {
        return 1;
    }
    d3 = sqrt(pv3[0] * pv3[0] + pv3[1] * pv3[1] + pv3[2] * pv3[2]);
    flag = 0;
    for (i = 0; i < 6; i++) {
	if (fabs(pv1[i] + pv2[i] - pv3[i]) > 1.0e-6) {
	    printf ("diff %d = %e\n", i, fabs(pv1[i] + pv2[i] - pv3[i]));
	    flag++;
	}
    }
    if ((d1 < 1.0e-5) || (d2 < 1.0e-10) || (d3 < 1.0e-10)) {
	flag++;
    }
    d = (d1 + d2 + d3) / get_au();
    printf ("EM_Barycenter - Geocenter - Moon (%7.4f au) ", d);
    if (flag) {
	printf ("FAILED\n");
    } else {
	printf ("Passed\n");
    }
    
    if (!get_earth_posn_vel(pv1, EM_Barycentric, jd, day)) {
        return 1;
    }
    d1 = sqrt(pv1[0] * pv1[0] + pv1[1] * pv1[1] + pv1[2] * pv1[2]);
    if (!get_sun_posn_vel(pv2, Geocentric, jd, day)) {
        return 1;
    }
    d2 = sqrt(pv2[0] * pv2[0] + pv2[1] * pv2[1] + pv2[2] * pv2[2]);
    if (!get_sun_posn_vel(pv3, EM_Barycentric, jd, day)) {
        return 1;
    }
    d3 = sqrt(pv3[0] * pv3[0] + pv3[1] * pv3[1] + pv3[2] * pv3[2]);
    flag = 0;
    for (i = 0; i < 6; i++) {
	if (fabs(pv1[i] + pv2[i] - pv3[i]) > 1.0e-6) {
	    printf ("diff %d = %e\n", i, fabs(pv1[i] + pv2[i] - pv3[i]));
	    flag++;
	}
    }
    if ((d1 < 1.0e-5) || (d2 < 1.0e-10) || (d3 < 1.0e-10)) {
	flag++;
    }
    d = (d1 + d2 + d3) / get_au();
    printf ("EM_Barycenter - Geocenter - Sun  (%7.4f au) ", d);
    if (flag) {
	printf ("FAILED\n");
    } else {
	printf ("Passed\n");
    }
    
    if (!get_planet_posn_vel(pv1, Pluto, SS_Barycentric, jd, day)) {
        return 1;
    }
    d1 = sqrt(pv1[0] * pv1[0] + pv1[1] * pv1[1] + pv1[2] * pv1[2]);
    if (!get_planet_posn_vel(pv2, Pluto, Heliocentric, jd, day)) {
        return 1;
    }
    d2 = sqrt(pv2[0] * pv2[0] + pv2[1] * pv2[1] + pv2[2] * pv2[2]);
    if (!get_sun_posn_vel(pv3, SS_Barycentric, jd, day)) {
        return 1;
    }
    d3 = sqrt(pv3[0] * pv3[0] + pv3[1] * pv3[1] + pv3[2] * pv3[2]);
    flag = 0;
    for (i = 0; i < 6; i++) {
	if (fabs(pv1[i] - pv2[i] - pv3[i]) > 1.0e-6) {
	    printf ("diff %d = %e\n", i, fabs(pv1[i] - pv2[i] - pv3[i]));
	    flag++;
	}
    }
    if ((d1 < 1.0e-5) || (d2 < 1.0e-10) || (d3 < 1.0e-10)) {
	flag++;
    }
    d = (d1 + d2 + d3) / get_au();
    printf ("Pluto - SS_Barycenter - Sun      (%7.4f au) ", d);
    if (flag) {
	printf ("FAILED\n");
    } else {
	printf ("Passed\n");
    }

    if (!get_planet_posn_vel(pv1, Venus, Heliocentric, jd, day)) {
        return 1;
    }
    d1 = sqrt(pv1[0] * pv1[0] + pv1[1] * pv1[1] + pv1[2] * pv1[2]);
    if (!get_planet_posn_vel(pv2, Venus, EM_Barycentric, jd, day)) {
        return 1;
    }
    d2 = sqrt(pv2[0] * pv2[0] + pv2[1] * pv2[1] + pv2[2] * pv2[2]);
    if (!get_sun_posn_vel(pv3, EM_Barycentric, jd, day)) {
        return 1;
    }
    d3 = sqrt(pv3[0] * pv3[0] + pv3[1] * pv3[1] + pv3[2] * pv3[2]);
    flag = 0;
    for (i = 0; i < 6; i++) {
	if (fabs(pv1[i] - pv2[i] + pv3[i]) > 1.0e-6) {
	    printf ("diff %d = %e\n", i, fabs(pv1[i] - pv2[i] + pv3[i]));
	    flag++;
	}
    }
    if ((d1 < 1.0e-5) || (d2 < 1.0e-10) || (d3 < 1.0e-10)) {
	flag++;
    }
    d = (d1 + d2 + d3) / get_au();
    printf ("Venus - EM_Barycenter - Sun      (%7.4f au) ", d);
    if (flag) {
	printf ("FAILED\n");
    } else {
	printf ("Passed\n");
    }

    if (!get_moon_posn_vel(pv1, Topocentric, jd, day)) {
        return 1;
    }
    d1 = sqrt(pv1[0] * pv1[0] + pv1[1] * pv1[1] + pv1[2] * pv1[2]);
    if (!get_moon_posn_vel(pv2, Heliocentric, jd, day)) {
        return 1;
    }
    d2 = sqrt(pv2[0] * pv2[0] + pv2[1] * pv2[1] + pv2[2] * pv2[2]);
    if (!get_sun_posn_vel(pv3, Topocentric, jd, day)) {
        return 1;
    }
    d3 = sqrt(pv3[0] * pv3[0] + pv3[1] * pv3[1] + pv3[2] * pv3[2]);
    flag = 0;
    for (i = 0; i < 6; i++) {
	if (fabs(pv1[i] - pv2[i] - pv3[i]) > 1.0e-6) {
	    printf ("diff %d = %e\n", i, fabs(pv1[i] - pv2[i] - pv3[i]));
	    flag++;
	}
    }
    if ((d1 < 1.0e-5) || (d2 < 1.0e-10) || (d3 < 1.0e-10)) {
	flag++;
    }
    d = (d1 + d2 + d3) / get_au();
    printf ("Observatory - Moon - Sun         (%7.4f au) ", d);
    if (flag) {
	printf ("FAILED\n");
    } else {
	printf ("Passed\n");
    }
    
    if (!get_planet_posn_vel(pv1, Jupiter, Topocentric, jd, day)) {
        return 1;
    }
    d1 = sqrt(pv1[0] * pv1[0] + pv1[1] * pv1[1] + pv1[2] * pv1[2]);
    if (!get_planet_posn_vel(pv2, Jupiter, Heliocentric, jd, day)) {
        return 1;
    }
    d2 = sqrt(pv2[0] * pv2[0] + pv2[1] * pv2[1] + pv2[2] * pv2[2]);
    if (!get_sun_posn_vel(pv3, Topocentric, jd, day)) {
        return 1;
    }
    d3 = sqrt(pv3[0] * pv3[0] + pv3[1] * pv3[1] + pv3[2] * pv3[2]);
    flag = 0;
    for (i = 0; i < 6; i++) {
	if (fabs(pv1[i] - pv2[i] - pv3[i]) > 1.0e-6) {
	    printf ("diff %d = %e\n", i, fabs(pv1[i] - pv2[i] - pv3[i]));
	    flag++;
	}
    }
    if ((d1 < 1.0e-5) || (d2 < 1.0e-10) || (d3 < 1.0e-10)) {
	flag++;
    }
    d = (d1 + d2 + d3) / get_au();
    printf ("Observatory - Jupiter - Sun      (%7.4f au) ", d);
    if (flag) {
	printf ("FAILED\n");
    } else {
	printf ("Passed\n");
    }
    
    if (!get_earth_posn_vel(pv1, EM_Barycentric, jd, day)) {
        return 1;
    }
    d1 = sqrt(pv1[0] * pv1[0] + pv1[1] * pv1[1] + pv1[2] * pv1[2]);
    if (!get_sun_posn_vel(pv2, Geocentric, jd, day)) {
        return 1;
    }
    d2 = sqrt(pv2[0] * pv2[0] + pv2[1] * pv2[1] + pv2[2] * pv2[2]);
    if (!get_sun_posn_vel(pv3, EM_Barycentric, jd, day)) {
        return 1;
    }
    d3 = sqrt(pv3[0] * pv3[0] + pv3[1] * pv3[1] + pv3[2] * pv3[2]);
    flag = 0;
    for (i = 0; i < 6; i++) {
	if (fabs(pv1[i] + pv2[i] - pv3[i]) > 1.0e-6) {
	    printf ("diff %d = %e\n", i, fabs(pv1[i] + pv2[i] - pv3[i]));
	    flag++;
	}
    }
    if ((d1 < 1.0e-5) || (d2 < 1.0e-10) || (d3 < 1.0e-10)) {
	flag++;
    }
    d = (d1 + d2 + d3) / get_au();
    printf ("EM_Barycenter - Geocenter - Sun  (%7.4f au) ", d);
    if (flag) {
	printf ("FAILED\n");
    } else {
	printf ("Passed\n");
    }
    printf("\n");
    close_ephemeris();
    return 0;
}
