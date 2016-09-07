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
    int day, month, year;

    if (argc < 2) {
	printf ("Usage: mjd2mdy MJD\n");
	return 1;
    }
    long mjd = atol(argv[1]);

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
    long days_since = mjd - 49352;
 
    year = 0;
    for (int y = 0; y < sizeof(ydays) / sizeof(int); y++) {
        if (days_since <= ydays[y]) {
            year = 1994 + y;
            if (ydays[y] == 366) {
                mdays[1] = 29;
            }
            break;
        }
        days_since -= ydays[y];
    }
    if (year == 0) {
	printf ("MJD year appears to be out of range ( > 2050 )\n");
	return 1;
    }
    for (int m = 0; m < 12; m++) {
        if (days_since <= mdays[m]) {
            month = m + 1;
            day = (int)days_since;
            break;
        }
        days_since -= mdays[m];
    }
    printf ("MJD %d == %02d/%02d/%02d (mm/dd/yy)\n", mjd, month, day, year);
    return 0;
}
