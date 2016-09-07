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

#include <stdio.h>
#include <stdlib.h>

int main ( int argc, char *argv[] )
{
    if (argc < 4) {
	printf ("Usage: mjd2mdy month day year ( 10 23 1996 )\n");
	return 1;
    }
    int month = atoi(argv[1]);
    if (month < 1 || month > 12) {
	printf ("month out of range\n");
	return 1;
    }
    int day   = atoi(argv[2]);
    if (day < 1 || day > 31) {
	printf ("day out of range\n");
	return 1;
    }
    int year  = atoi(argv[3]);
    if (year <= 50) {
	year += 2000;
    } else if (year < 100) {
	year += 1900;
    }
    if (year < 1994 || year > 2050) {
	printf ("Year must be between 1994 and 2050\n");
	return 1;
    }
    // days in the years up to 2050
    static
    const int ydays[] = {      365, 365, 366, 365, 365, 365, 366,
                          365, 365, 365, 366, 365, 365, 365, 366,
                          365, 365, 365, 366, 365, 365, 365, 366,
                          365, 365, 365, 366, 365, 365, 365, 366,
                          365, 365, 365, 366, 365, 365, 365, 366,
                          365, 365, 365, 366, 365, 365, 365, 366,
                          365, 365, 365, 366, 365, 365, 365, 366,
                          365, 365 };
 
    static int mdays[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
 
    // MJD for 0h Jan 0, 1994 = 49352
    long mjd = 49352;
 
    int yr = year;
    for (int y = 0; yr > 1994; y++, yr--) {
	mjd += ydays[y];
    }
    for (int m = 0; m < month - 1; m++) {
	mjd += mdays[m];
    }
    mjd += day;
    if (ydays[year - 1994] == 366 && month > 2) mjd++;

    printf ("%02d/%02d/%02d == MJD %d\n", month, day, year, mjd);
    return 0;
}
