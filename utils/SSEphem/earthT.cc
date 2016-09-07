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

// JPL DExxx earth's position test main program
// Compare with Astronomical Almanac even # pages B44 - B58

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "jpl_eph.h"

int main (int argc, char *argv[])
{
    if (argc < 2) {
        printf ("Usage: earthT ephemeris_file_name\n");
        return 1;
    } else {
	printf ("\nJPL Earth position test.  Compare with\n");
	printf ("1996 Astronomical Almanac even # pages B44 - B58\n\n");
    }
    if (!open_ephemeris(argv[1])) {
        return 1;
    }
    double jd = 2450082.5;    // Jan 0, 1996
    double earth_pv[6];

    // print position for 0h TT on Jan 1, Apr 10, July 19, and Oct. 27, 1996
    static char gdate[4][20] = { "Jan 1, 1996", "Apr 10, 1996",
				 "July 19, 1996", "Oct. 27, 1996" };
    int i = 0;
    double day;
    for (day = 1.0; day < 365; day += 100.0) {
        if (get_earth_posn_vel(earth_pv, SS_Barycentric, jd, day)) {
	    printf ("TT JD %5.1f   %s\n", jd + day, gdate[i]);
	    printf (
	     "x,y,z    %12.9f %12.9f %12.9f\nxv,yv,zv %12.0f %12.0f %12.0f\n\n",
		earth_pv[0] / get_au(), earth_pv[1] / get_au(),
		earth_pv[2] / get_au(), 1.0e9 * earth_pv[3] / get_au(),
		1.0e9 * earth_pv[4] / get_au(), 1.0e9 * earth_pv[5] / get_au());
        } else {
            printf ("Earth position calculation failed\n");
        }
	i++;
    }
    return 0;
}
