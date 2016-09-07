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
// This program converts the given MJD and UTC and prints the corresponding
// sidereal time.  The hard-coded observatory coordinates are used.

#include <stdio.h>
#include <string.h>
#include "sofa.h"
#include "astrtime.h"
#include "ptw_aux.h"

double str2days (char *str)
{
    char wstr[20];
    strcpy (wstr, str);
    // convert delimiters to spaces
    int i;
    for (i = 0; i < strlen(wstr); i++) {
        if (wstr[i] == ':') {
            wstr[i] = ' ';
        }
    }
    double hrs, min, sec;
    if (sscanf (wstr, "%lf %lf %lf", &hrs, &min, &sec) != 3) {
        printf ("Error converting time string '%s'\n", str);
        return 0.0;
    }
    return hrs / 24.0 + min / 1440.0 + sec /86400.0;
}

void str2ymd(char *str, int &yr, int &mo, int &da)
{
    char wstr[20];
    strcpy (wstr, str);
    // convert delimiters to spaces
    int i;
    for (i = 0; i < strlen(wstr); i++) {
        if (wstr[i] == '/') {
            wstr[i] = ' ';
        }
    }
    if (sscanf (wstr, "%d %d %d", &mo, &da, &yr) != 3) {
        printf ("Error converting time string '%s'\n", str);
    }
}

char *mo2str(int mo)
{
    if (mo <= 0) {
        return "";
    }
    switch (mo) {
        case 1: return "Jan";
        case 2: return "Feb";
        case 3: return "Mar";
        case 4: return "Apr";
        case 5: return "May";
        case 6: return "Jun";
        case 7: return "Jul";
        case 8: return "Aug";
        case 9: return "Sep";
        case 10: return "Oct";
        case 11: return "Nov";
        case 12: return "Dec";
        default: return "";
    }
}

char *lst2str(double lst)
{
    static char str[20];
    lst *= 24.0;
    int hr = (int)lst;
    lst -= (double)hr;
    lst *= 60.0;
    int min = (int)lst;
    lst -= (double)min;
    lst *= 60.0;
    if (lst >= 59.995) {
        lst = 0.0;
        min++;
        if (min >= 60) {
            min--;
            hr++;
        }
    }
    sprintf (str, "%d %d %5.2lf", hr, min, lst);
    return str;
}

int main (int argc, char *argv[])
{
    if (argc != 3) {
        printf ("Usage: lst hh:mm:ss mm/dd/yyyy\n");
        printf ("       time and date are UTC\n");
        return 1;
    }
    double utc = str2days(argv[1]);
    int yr, mo, da, status;
    str2ymd(argv[2], yr, mo, da);

    double dmjd;
    ptwCaldj (yr, mo, da, &dmjd, &status);

    long mjd = (long)(dmjd + 0.1);
    double lst = utc2lmst (mjd, utc);
    printf ("LMST for %s UTC  %s %d, %d is  %s\n",
            argv[1], mo2str(mo), da, yr, lst2str(lst));
    return 0;
}
