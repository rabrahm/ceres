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

// JPL DExxx Moon's velocity test main program

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
        printf ("Usage: moon_test ephemeris_file_name\n");
        return 1;
    } else {
	printf ("\nJPL Moon velocity test.\n");
    }
    if (!open_ephemeris(argv[1])) {
        return 1;
    }
    double jd = 2450082.5;    // Jan 0, 1996
    double tptr[6];
    double tptrx[6];

    // print position for 0h TT on Jan 1, Apr 10, July 19, and Oct. 27, 1996
    static char gdate[4][20] = { "Jan 1, 1996", "Apr 10, 1996",
				 "July 19, 1996", "Oct. 27, 1996" };
    int i = 0;
    double day;
//    for (day = 1.0; day < 365; day += 100.0) {
    for (day = 1.0; day < 5; day += 100.0) {
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
	    printf ("TT JD %6.3f   %s\n", jd + day, gdate[i]);
	    i++;
	    printf ("  J2000 RA %s  Dec %s, Dist. %10.8f au\n",
		    deg2hms_str(ra), deg2dms_str(dec), d / get_au());

            // compute again for 60 seconds later
            double dayx = day + 60.0 / 86400.0;
            get_moon_posn_vel(tptrx, Geocentric, jd, dayx);
            for (jxx = 0; jxx < 3; jxx++) posn[jxx] = tptrx[jxx];
	    double decx, rax;
	    iauC2s (posn, &rax, &decx);
	    decx *= RAD2DEG;
	    rax *= RAD2DEG;
	    if (rax < 0.0) rax += 360.0;
	    printf ("  J2000 RA %s  Dec %s\n",
		    deg2hms_str(rax), deg2dms_str(decx));

            double kpd2mpm = 60.0 * 360.0 / (TWOPI * d * 1440.0);
            printf ("Vx = %f, Vy = %f, Vz = %f\n",
                    kpd2mpm * tptr[3], kpd2mpm * tptr[4], kpd2mpm * tptr[5]);
            printf ("Delta RA = %10.7f,  Delta Dec = %10.7f  arcmin/min\n",
                    (rax - ra) * 60.0, (decx - dec) * 60.0);

            double dd = kpd2mpm * (tptr[5] * cos(dec / RAD2DEG) -
                                   (tptr[3] * cos(ra / RAD2DEG) +
                         tptr[4] * sin(ra / RAD2DEG)) * sin(dec / RAD2DEG));
//            printf ("Delta Dec = %10.7f  arcmin/min\n", dd);
            double ddx = kpd2mpm * (tptrx[5] * cos(decx / RAD2DEG) -
                                   (tptrx[3] * cos(rax / RAD2DEG) +
                         tptrx[4] * sin(rax / RAD2DEG)) * sin(decx / RAD2DEG));
//            printf ("Delta Dec = %10.7f  arcmin/min\n", ddx);
//            printf ("Delta Dec = %10.7f  arcmin/min\n", (ddx + dd) / 2.0);
            double dr = kpd2mpm * (tptr[4] * cos(ra / RAD2DEG) -
                         tptr[3] * sin(ra / RAD2DEG)) / cos(dec / RAD2DEG);
//            printf ("Delta RA = %10.7f  arcmin/min\n", dr);
            double drx = kpd2mpm * (tptrx[4] * cos(rax / RAD2DEG) -
                         tptrx[3] * sin(rax / RAD2DEG)) / cos(decx / RAD2DEG);
//            printf ("Delta RA = %10.7f  arcmin/min\n", drx);
//            printf ("Delta RA = %10.7f  arcmin/min\n", (drx + dr) / 2.0);
            printf ("Delta RA = %10.7f,  Delta Dec = %10.7f  arcmin/min\n",
                    (drx + dr) / 2.0, (ddx + dd) / 2.0);

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
	    printf ("\n  Appar RA %s  Dec %s\n",
		    deg2hms_str(ara), deg2dms_str(adec));

            // do it again for 60 seconds later
            double day_lightx = day_light + 60.0 / 86400.0;
	    get_moon_posn_vel(tptrx, Geocentric, jd, day_lightx);
	    // add precession and nutation
	    double rmatpnx[3][3], new_vecx[3];
	    date = jd + day + (60.0 / 86400.0) - 2400000.5;
	    iauPnm00a (2400000.5, date, rmatpnx);
            for (jxx = 0; jxx < 3; jxx++) posn[jxx] = tptrx[jxx];
	    iauRxp (rmatpnx, posn, new_vecx);
	    double arax, adecx;
	    iauC2s (new_vecx, &arax, &adecx);
	    adecx *= RAD2DEG;
	    arax *= RAD2DEG;
	    if (arax < 0.0) arax += 360.0;
	    printf ("  Appar RA %s  Dec %s\n",
		    deg2hms_str(arax), deg2dms_str(adecx));
            double vel[3], new_vel[3];
            for (jxx = 0; jxx < 3; jxx++) vel[jxx] = tptr[jxx + 3];
	    iauRxp (rmatpn, vel, new_vel);
            printf ("Vx = %f, Vy = %f, Vz = %f\n",
                    kpd2mpm * new_vel[0], kpd2mpm * new_vel[1],
                    kpd2mpm * new_vel[2]);
            printf ("Delta RA = %10.7f,  Delta Dec = %10.7f  arcmin/min\n",
                    (arax - ara) * 60.0, (adecx - adec) * 60.0);
            dd = kpd2mpm * (new_vel[2] * cos(adec / RAD2DEG) -
                                   (new_vel[0] * cos(ara / RAD2DEG) +
                         new_vel[1] * sin(ara / RAD2DEG)) * sin(adec / RAD2DEG));
//            printf ("Delta Dec = %10.7f  arcmin/min\n", dd);
            double velx[3], new_velx[3];
            for (jxx = 0; jxx < 3; jxx++) velx[jxx] = tptrx[jxx + 3];
	    iauRxp (rmatpnx, velx, new_velx);
            ddx = kpd2mpm * (new_velx[2] * cos(adecx / RAD2DEG) -
                                   (new_velx[0] * cos(arax / RAD2DEG) +
                    new_velx[1] * sin(arax / RAD2DEG)) * sin(adecx / RAD2DEG));
//            printf ("Delta Dec = %10.7f  arcmin/min\n", ddx);
//            printf ("Delta Dec = %10.7f  arcmin/min\n", (ddx + dd) / 2.0);
            dr = kpd2mpm * (new_vel[1] * cos(ara / RAD2DEG) -
                       new_vel[0] * sin(ara / RAD2DEG)) / cos(adec / RAD2DEG);
//            printf ("Delta RA = %10.7f  arcmin/min\n", dr);
            drx = kpd2mpm * (new_velx[1] * cos(arax / RAD2DEG) -
                    new_velx[0] * sin(arax / RAD2DEG)) / cos(adecx / RAD2DEG);
//            printf ("Delta RA = %10.7f  arcmin/min\n", drx);
//            printf ("Delta RA = %10.7f  arcmin/min\n", (drx + dr) / 2.0);
            printf ("Delta RA = %10.7f,  Delta Dec = %10.7f  arcmin/min\n",
                    (drx + dr) / 2.0, (ddx + dd) / 2.0);

	    // Now do topocentric calculation
	    if (!get_moon_posn_vel(tptr, Topocentric, jd, day))
		return 1;
	    double td = sqrt(tptr[0] * tptr[0] + tptr[1] * tptr[1]
			     + tptr[2] * tptr[2]);
            kpd2mpm = 60.0 * 360.0 / (TWOPI * td * 1440.0);
	    // recompute for light travel advanced time
	    // since the earth's position is computed for the earlier time,
	    // too, this accounts for aberration.
	    day_light = day - d / (C_LIGHT * 86400.0);
	    if (!get_moon_posn_vel(tptr, Topocentric, jd, day_light))
		return 1;
	    td = sqrt(tptr[0] * tptr[0] + tptr[1] * tptr[1]
			     + tptr[2] * tptr[2]);
	    // add precession and nutation
	    date = jd + day - 2400000.5;
	    iauPnm00a (2400000.5, date, rmatpn);
            for (jxx = 0; jxx < 3; jxx++) posn[jxx] = tptr[jxx];
	    iauRxp (rmatpn, posn, new_vec);
	    double tra, tdec;
	    iauC2s (new_vec, &tra, &tdec);
	    tdec *= RAD2DEG;
	    tra *= RAD2DEG;
	    if (tra < 0.0) tra += 360.0;
	    printf ("\n  Topoc RA %s  Dec %s  Dist. %8.2f km\n",
		    deg2hms_str(tra), deg2dms_str(tdec), td);

            // do it again for 60 seconds later
            day_lightx = day_light + 60.0 / 86400.0;
            get_moon_posn_vel(tptrx, Topocentric, jd, day_lightx);
	    double tdx = sqrt(tptrx[0] * tptrx[0] + tptrx[1] * tptrx[1]
			     + tptrx[2] * tptrx[2]);
            double kpd2mpmx = 60.0 * 360.0 / (TWOPI * tdx * 1440.0);
	    // add precession and nutation
	    date = jd + day + (60.0 / 86400.0) - 2400000.5;
	    iauPnm00a (2400000.5, date, rmatpnx);
            for (jxx = 0; jxx < 3; jxx++) posn[jxx] = tptrx[jxx];
	    iauRxp (rmatpnx, posn, new_vecx);
	    double trax, tdecx;
	    iauC2s (new_vecx, &trax, &tdecx);
	    tdecx *= RAD2DEG;
	    trax *= RAD2DEG;
	    if (trax < 0.0) trax += 360.0;
	    printf ("  Topoc RA %s  Dec %s\n",
		    deg2hms_str(trax), deg2dms_str(tdecx));
            for (jxx = 0; jxx < 3; jxx++) vel[jxx] = tptr[jxx + 3];
	    iauRxp (rmatpn, vel, new_vel);
            printf ("Vx = %f, Vy = %f, Vz = %f\n",
                    kpd2mpm * new_vel[0], kpd2mpm * new_vel[1],
                    kpd2mpm * new_vel[2]);
            printf ("Delta RA = %10.7f,  Delta Dec = %10.7f  arcmin/min\n",
                    (trax - tra) * 60.0, (tdecx - tdec) * 60.0);
            dd = kpd2mpm * (new_vel[2] * cos(tdec / RAD2DEG) -
                                   (new_vel[0] * cos(tra / RAD2DEG) +
                      new_vel[1] * sin(tra / RAD2DEG)) * sin(tdec / RAD2DEG));
//            printf ("Delta Dec = %10.7f  arcmin/min\n", dd);
            for (jxx = 0; jxx < 3; jxx++) velx[jxx] = tptrx[jxx + 3];
	    iauRxp (rmatpnx, velx, new_velx);
            ddx = kpd2mpmx * (new_velx[2] * cos(tdecx / RAD2DEG) -
                                   (new_velx[0] * cos(trax / RAD2DEG) +
                    new_velx[1] * sin(trax / RAD2DEG)) * sin(tdecx / RAD2DEG));
//            printf ("Delta Dec = %10.7f  arcmin/min\n", ddx);
//            printf ("Delta Dec = %10.7f  arcmin/min\n", (ddx + dd) / 2.0);
            dr = kpd2mpm * (new_vel[1] * cos(tra / RAD2DEG) -
                       new_vel[0] * sin(tra / RAD2DEG)) / cos(tdec / RAD2DEG);
//            printf ("Delta RA = %10.7f  arcmin/min\n", dr);
            drx = kpd2mpmx * (new_velx[1] * cos(trax / RAD2DEG) -
                    new_velx[0] * sin(trax / RAD2DEG)) / cos(tdecx / RAD2DEG);
//            printf ("Delta RA = %10.7f  arcmin/min\n", drx);
//            printf ("Delta RA = %10.7f  arcmin/min\n", (drx + dr) / 2.0);
            printf ("Delta RA = %10.7f,  Delta Dec = %10.7f  arcmin/min\n",
                    (drx + dr) / 2.0, (ddx + dd) / 2.0);
        } else {
            printf ("Moon position calculation failed\n");
        }
    }
    return 0;
}
