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

// This program is a temporary one to generate IERS table entries of
// UT1-UTC values from a table of TAI-UT1 values in the TEMPO directory.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "astrtime.h"

//#define PI                   3.14159265358979323846
#define TWOPI                (2.0 * PI)
#define FOURPI               (4.0 * PI)

int main(void)
{
    FILE *iers = fopen("iers.tmp", "w");
    FILE *tai = fopen("ut1.dat", "rt");

    char line[85];
    long mjdc = 0;
    int ln = 0;
    // skip first 2 lines
    fgets (line, 84, tai);
    fgets (line, 84, tai);
    while (fgets (line, 84, tai) != NULL) {
        line[strlen(line)-1] = '\0';
        ln++;
        long mjd;
        double x[6];
        int num = sscanf(line, "%ld %lf %lf %lf %lf %lf %lf",
                   &mjd, &x[0], &x[1], &x[2], &x[3], &x[4], &x[5]) - 1;
        if ((mjdc != 0) && (mjdc != mjd)) {
            printf ("Following line sync error, line %d\n", ln);
            break;
        } else {
            mjdc = mjd;
        }
        double xp = 0.0;
        double yp = 0.0;
        if (mjd > 49400) {
	    int i;
            for (i = 0; i < num; i++) {
                x[i] /= 10000.0;
                x[i] -= get_leap_sec(mjdc, 0.01);
                fprintf (iers, "%6ld %7.4lf %7.4lf %8.5lf\n",
                                                     mjdc, xp, yp, -x[i]);
                mjdc += 5;
            }
        } else {
            num = 6;
            mjdc = 0;
        }
        if (num != 6) break;
    }
    printf ("  %d lines processed\n", ln);
    fclose(iers);
    fclose(tai);
    return 0;
}
