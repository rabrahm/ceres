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

// Green Bank station coordinates check

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "jpl_eph.h"
#include "sofa.h"
#include "astrtime.h"
#include "ptw_aux.h"

// X, Y, Z coordinates for the Green Bank telescopes accirate to about 0.003 m
// NRAO 140-ft
//  882880.0208   -4924482.4385    3944130.6438
// NRAO 85-1
//  883555.6841   -4924491.0244    3943961.9293
// NRAO 85-3
//  882325.6940   -4925138.1167    3943397.6444

typedef struct {
    double x;
    double y;
    double z;
} topo;

topo  n140;
topo  n85_1;
topo  n85_3;

int main (int argc, char *argv[])
{
    if (argc < 2) {
	printf ("Usage: topochk ephemeris_file\n");
	return 1;
    }
    if (!open_ephemeris(argv[1])) {
        return 1;
    }
    n140.x =   882880.0208;
    n140.y = -4924482.4385;
    n140.z =  3944130.6438;

    n85_1.x =   883555.6841;
    n85_1.y = -4924491.0244;
    n85_1.z =  3943961.9293;

    n85_3.x =   882325.6940;
    n85_3.y = -4925138.1167;
    n85_3.z =  3943397.6444;

    char name[30];
    double west_long;
    double latitude;
    double height;

    ptwObs (0, "GBVA140", name, &west_long, &latitude, &height);
    printf ("%s, wlong = %f, lat = %f, height = %f\n", name, west_long,
	    latitude, height);
    double radius, z;
    ptwGeoc (latitude, height, &radius, &z);
    printf ("R = %f, Z = %f  (from Starlink coordinates)\n",
	    get_au() * radius, get_au() * z);
    radius = sqrt(n140.x * n140.x + n140.y * n140.y) / 1000.0;
    z = n140.z / 1000.0;
    printf ("R = %f, Z = %f  (from VLBI)\n", radius, z);
    ptwObs (0, "GBVA300", name, &west_long, &latitude, &height);
    printf ("%s, wlong = %f, lat = %f, height = %f\n", name, west_long,
	    latitude, height);
    ptwGeoc (latitude, height, &radius, &z);
    printf ("R = %f, Z = %f  (from Starlink coordinates)\n\n",
	    get_au() * radius, get_au() * z);

    double pv[6];
//    double jd = 2450348.5;	// Sept 22, 1996
    double jd = 2452261.5;	// Dec 18, 2001
    double utc = 5.20 / 24.0;	// roughly 0h local apparent solar time
    get_observatory_posn_vel(pv, jd, utc);
    printf ("Obs:   %13.2f %13.2f %13.2f   %8.4f %8.4f %8.4f\n",
	    pv[0], pv[1], pv[2],
	    pv[3] / 86400.0, pv[4] / 86400.0, pv[5] / 86400.0);

    get_earth_posn_vel(pv, Heliocentric, jd, utc);
    printf ("Earth: %13.2f %13.2f %13.2f   %8.4f %8.4f %8.4f\n",
	    pv[0], pv[1], pv[2],
	    pv[3] / 86400.0, pv[4] / 86400.0, pv[5] / 86400.0);
    return 0;
}
