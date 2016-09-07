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

// This program is a temporary one to generate an IERS table of pole
// offsets and UT1-UTC values from the extrapolation given in the IERS
// Bulletin A of 27 June, 1996

#include <stdio.h>
#include <math.h>

#define PI                   3.14159265358979323846
#define TWOPI                (2.0 * PI)
#define FOURPI               (4.0 * PI)

int main(void)
{
    FILE *iers = fopen("iers.tmp", "w");
    double last_x, last_y, last_dut1, dx, dy, ddut1;
    last_x = last_y = last_dut1 = dx = dy = ddut1 = 0.0;

/*
    // Oct. 10, 1995 -> Dec 31, 1996
    long mjd;
    for (mjd = 50200; mjd < 50650; mjd += 5) {
        double by = 1995.0 + (double)(mjd - 49718) / 365.25;
        double dut2 =  0.022 * sin(TWOPI * by)
                     - 0.012 * cos(TWOPI * by)
                     - 0.006 * sin(FOURPI * by)
                     + 0.007 * cos(FOURPI * by);
        double a = TWOPI * (double)(mjd - 50260) / 365.25;
        double c = TWOPI * (double)(mjd - 50260) / 435.  ;
        double x = 0.0319 + 0.0019 * cos(a) + 0.0685 * sin(a) +
                            0.1172 * cos(c) + 0.1740 * sin(c);
        double y = 0.3342 + 0.0617 * cos(a) - 0.0109 * sin(a) +
                            0.1740 * cos(c) - 0.1172 * sin(c);
        double dut1 = 0.2274 - 0.00207 * (double)(mjd - 50259) - dut2;

    // Sept. 10, 1998 -> Dec 31, 1998
    long mjd;
    for (mjd = 51067; mjd < 51178; mjd += 5) {
        double by = 1995.0 + (double)(mjd - 49718) / 365.25;
        double dut2 =  0.022 * sin(TWOPI * by)
                     - 0.012 * cos(TWOPI * by)
                     - 0.006 * sin(FOURPI * by)
                     + 0.007 * cos(FOURPI * by);
        double a = TWOPI * (double)(mjd - 50701) / 365.25;
        double c = TWOPI * (double)(mjd - 50701) / 435.  ;
        double x = 0.0432 + 0.0267 * cos(a) - 0.0207 * sin(a) +
                            0.1095 * cos(c) + 0.1589 * sin(c);
        double y = 0.3339 - 0.0164 * cos(a) - 0.0507 * sin(a) +
                            0.1589 * cos(c) - 0.1095 * sin(c);
        double dut1 = 0.3848 - 0.00181 * (double)(mjd - 50701) - dut2;

    // March 4, 1999 -> Dec 31, 1999
    long mjd;
    for (mjd = 51242; mjd < 51543; mjd += 5) {
        double by = 1995.0 + (double)(mjd - 49718) / 365.25;
        double dut2 =  0.022 * sin(TWOPI * by)
                     - 0.012 * cos(TWOPI * by)
                     - 0.006 * sin(FOURPI * by)
                     + 0.007 * cos(FOURPI * by);
        double a = TWOPI * (double)(mjd - 50876) / 365.25;
        double c = TWOPI * (double)(mjd - 50876) / 435.  ;
        double x = 0.0026 - 0.0546 * cos(a) + 0.0285 * sin(a) +
                            0.0131 * cos(c) - 0.1755 * sin(c);
        double y = 0.3272 + 0.0239 * cos(a) + 0.0742 * sin(a) -
                            0.1755 * cos(c) - 0.0131 * sin(c);
        double dut1 = 0.1094 - 0.00184 * (double)(mjd - 50876) - dut2;
*/
    // January 19, 2001 -> Dec 31, 2001
    long mjd;
    for (mjd = 51926; mjd < 52276; mjd += 5) {
        double by = 1995.0 + (double)(mjd - 49718) / 365.25;
        double dut2 =  0.022 * sin(TWOPI * by)
                     - 0.012 * cos(TWOPI * by)
                     - 0.006 * sin(FOURPI * by)
                     + 0.007 * cos(FOURPI * by);
        double a = TWOPI * (double)(mjd - 51561) / 365.25;
        double c = TWOPI * (double)(mjd - 51561) / 435.  ;
        double x = 0.0395 - 0.0564 * cos(a) - 0.0755 * sin(a) +
                            0.0649 * cos(c) + 0.1150 * sin(c);
        double y = 0.3335 - 0.0693 * cos(a) + 0.0461 * sin(a) +
                            0.1150 * cos(c) - 0.0649 * sin(c);
        double dut1 = 0.3365 - 0.00097 * (double)(mjd - 51561) - dut2;
/*
        if (mjd > 50082) dut1 += 1.0;
        // write an entry for the leap second epoch
        if (mjd == 50085) {
            fprintf (iers, "%6ld %7.4lf %7.4lf %8.5lf <<<<<< leap second\n",
                       50083, x - 0.6 * dx, y - 0.6 * dy, dut1 - 0.6 * ddut1);
        }
*/
        fprintf (iers, "%6ld %7.4lf %7.4lf %8.5lf\n", mjd, x, y, dut1);

        dx = x - last_x;
        dy = y - last_y;
        ddut1 = dut1 - last_dut1;
        last_x = x;
        last_y = y;
        last_dut1 = dut1;
    }
    fclose(iers);
    return 0;
}
