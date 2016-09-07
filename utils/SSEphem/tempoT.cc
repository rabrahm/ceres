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

// This program generates an ascii file with pulse arrival times in a format
// readable by the 'tempo' program.  The arrival times are predicted times
// using a given P, P_dot, and epoch and the solar system delay calculations
// from the 'jpl_eph.cc' and 'astrtime.cc' routines.  This permits an accurate
// comparison between the 'tempo' and these routines.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "jpl_eph.h"
#include "astrtime.h"
#include "delay.h"
#include "sofa.h"
#include "ptw_aux.h"

//#define RA	192952.042000
//#define DEC	105304.48000

//#define RA	000000.0
//#define DEC	100000.0

#define RA	080000.0
#define DEC	100000.0

//#define RA	180000.0
//#define DEC	664000.0

#define PERIOD	0.2265171544090938
#define P_DOT	1.1552965958e-15
#define EPOCH	41703.500000
#define DM	3.176000
#define P_EPOCH	EPOCH		// position epoch

#define BJ_NUM	0		// 0=B1950, 1=J2000
#define EPHEM	0		// 0=PEP740R, 1=DE200, 2=DE202, 3=DE211

//#define J2000_CHK

void get_2000(double &ra2000, double &dec2000)
{
    double hrs = (int)(RA / 10000.0);
    double ra = RA - hrs * 10000.0;
    double min = (int)(ra / 100.0);
    double sec = ra - min * 100.0;
    ra = TWOPI * (hrs + (min + sec / 60.0) / 60.0) / 24.0;

    double sign = 1.0;
    if (DEC < 0.0) sign = -1.0;
    double dec = fabs(DEC);
    double deg = (int)(dec / 10000.0);
    dec -= deg * 10000.0;
    min = (int)(dec / 100.0);
    sec = dec - min * 100.0;
    dec = sign * TWOPI * (deg + (min + sec / 60.0) / 60.0) / 360.0;

    double b_epoch = iauEpb(2400000.5, P_EPOCH);
    ptwFk45z(ra, dec, b_epoch, &ra2000, &dec2000);
    if (BJ_NUM == 1) {
	ra2000 = ra;
	dec2000 = dec;
    }
#ifdef J2000_CHK	// alternate routine for getting B1950 -> J2000
    double rax, decx, drax, ddecx, px, vx;
    slaFk425(ra,    dec,   0.0,   0.0,   0.0, 0.0,
	     &rax, &decx, &drax, &ddecx, &px, &vx);
    printf ("RA %12.4lf  DEC %12.4lf\n",
	    ra2000 * 206264.8, dec2000 * 206264.8);
    printf ("RA %12.4lf  DEC %12.4lf\n",
	    rax * 206264.8, decx * 206264.8);
    printf ("PMra %8.6lf  PMdec %8.6lf  PAR %8.6lf,  RV %8.6lf\n",
	    drax * 206264.8, ddecx * 206264.8, px, vx);
#endif
}

double create_toa(long mjd, double utc)
{
    double p_freq = 1.0 / PERIOD;
    double p_f_dot = -P_DOT / (PERIOD * PERIOD);
    double ra2000, dec2000;
    get_2000(ra2000, dec2000);

    long tdb_mjd;
    double tdb;
    utc2tdb(mjd, utc, tdb_mjd, tdb);
    double t = (double)tdb_mjd - EPOCH;
    t = (t + tdb) * 86400.0;
    double pulse_phase = fmod(p_freq * t + p_f_dot * t * t / 2.0, 1.0);
    double pulse_freq = p_freq + p_f_dot * t;
    double time = utc - (pulse_phase / pulse_freq) / 86400.0;
    time += fmod (get_geometric_delay (SS_Barycentric, Topocentric,
				       mjd, utc, ra2000, dec2000),
		  1.0 / pulse_freq) / 86400.0;
    long time_mjd = mjd;
    if (time < 0.0) {
	time += 1.0;
	time_mjd--;
    }
    if (time > 1.0) {
	time -= 1.0;
	time_mjd++;
    }
    return (double)time_mjd + time;
}

int main (int argc, char *argv[])
{
    if (argc < 3) {
        printf ("Usage: tempoT ephemeris_file_name TOA_file_name\n");
        return 1;
    } else {
        printf ("\n'tempo' timing comparison\n");
    }
    if (!open_ephemeris(argv[1])) {
        return 1;
    }
    FILE *tempo = fopen(argv[2], "wt");
    if (tempo == NULL) {
	printf ("Cannot open file %s\n", argv[2]);
	return 1;
    }
    fprintf (tempo, "000 1100       0                             %d%d\n",
	     EPHEM, BJ_NUM);
    fprintf (tempo, "TESTpsr%33.6f%20.5f\n", RA, DEC);
    fprintf (tempo, "P%18.16f%19.10f%15.6f%17.1f\n",
	     PERIOD, 1.0e15 * P_DOT, EPOCH, 0.0);
    fprintf (tempo, "%21.6f\n", DM);

//set_latlong2xyz (-79.8360028, 38.437614, 880.9);
    double utc = 0.5;
    for (long mjd = 49465; mjd < 50040; mjd += 20) {
//long mjd = 49555;
//for (utc = 0.0; utc < 1.0; utc += 0.05) {
	double toax = create_toa(mjd, utc);
	fprintf (tempo, "a%22.3f%18.11f%12.1f%18.1f\n",
		 1410.000, toax, 1.0, 0.0);
    }
    close_ephemeris();
    return 0;
}
/*
1929+10                    192952.042000        105304.48000
P0.2265171544090938       1.1552965958   41703.500000              0.0
             3.176000
a              1410.000 49465.38767536036      1769.7               0.0        
a              1410.000 49481.42702670563      1769.7               0.0        
a              1410.000 49526.22971275058      1769.7               0.0        
a              1410.000 49554.18576569262      1769.7               0.0        
a              1410.000 49581.08224824225      1769.7               0.0        
a              1410.000 49581.08834663018      1769.7               0.0        
a              1410.000 49581.09288502779      1769.7               0.0        
a              1410.000 49635.05378722316      1769.7               0.0        
a              1410.000 49827.45728159377      1769.7               0.0        
a              1410.000 49864.25403042885      1769.7               0.0        
a              1410.000 49881.37570751099      1769.7               0.0        
a              1410.000 49885.39707029672      1769.7               0.0        
a              1410.000 49918.26572821527      1769.7               0.0        
a              1410.000 49919.25491029026      1769.7               0.0        
a              1410.000 49919.25787810236      1769.7               0.0        
a              1410.000 49979.11528950181      1769.7               0.0        
a              1410.000 50001.89377003431      1769.7               0.0        
a              1410.000 50038.90196675479      1769.7               0.0        
a              1410.001 50052.88851247383      1769.7               0.0        
*/
