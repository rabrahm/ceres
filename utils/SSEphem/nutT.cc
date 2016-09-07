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

// JPL DE403 nutation test main program

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "jpl_eph.h"

int main (int argc, char *argv[])
{
    if (argc < 2) {
        printf ("Usage: nutT ephemeris_file_name\n");
        return 1;
    } else {
        printf ("\nJPL Ephemerides nutation test.  Compare with\n");
	printf ("1996 Astronomical Almanac pages B24-31\n\n");
    }
    if (!open_ephemeris(argv[1])) {
        return 1;
    }
    // try a set of nutation calculations
    double jd = 2450082.5;    // Jan 0, 1996
    double tptr[4];

    // print nutation values for 0h TT on Jan 1, Apr 10, July 19,
    // and Oct. 27, 1996
    static char gdate[4][20] = { "Jan 1, 1996", "Apr 10, 1996",
				 "July 19, 1996", "Oct. 27, 1996" };
    int i = 0;
    double day;
    for (day = 1.0; day < 365; day += 100.0) {
        if (get_nutation(tptr, jd, day)) {
	    int ij;
	    for (ij = 0; ij < 4; ij++) {
		tptr[ij] *= 206264.81;
	    }
	    printf ("TT JD %6.3f   %s\n", jd + day, gdate[i]);
	    printf (
		"Long %6.3f, Obliq %6.3f, Long vel %6.3f, Obliq vel %6.3f\n\n",
                    tptr[0], tptr[1], tptr[2], tptr[3]);
        } else {
            printf ("Nutation calculation failed\n");
	    return 1;
        }
	i++;
    }
    return 0;
}
