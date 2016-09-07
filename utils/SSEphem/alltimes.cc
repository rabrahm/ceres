// This file is part of SSEphem.

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

// This program converts the given MJD and UTC into all of the times computed
// in the astrtime.cc routines.  The hard-coded observatory coordinates for
// the 140-ft are used.

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

char *day2time_str(double time, char *str)
{
    time *= 24.0;
    int hr = (int)time;
    time -= (double)hr;
    time *= 60.0;
    int min = (int)time;
    time -= (double)min;
    time *= 60.0;
    if (time >= 59.9995) {
        time = 0.0;
        min++;
        if (min >= 60) {
            min = 0;
            hr++;
        }
    }
    sprintf (str, "%2d %2d %6.3lf", hr, min, time);
    return str;
}

char *rad2dms_str(double rad, char *str)
{
    int sign = 1;
    if (rad < 0.0) sign = -1;
    rad = fabs(rad * 360.0 / TWOPI);
    int d = (int)rad;
    rad -= (double)d;
    rad *= 60.0;
    int m = (int)rad;
    rad -= (double)m;
    rad *= 60.0;
    if (rad >= 59.995) {
        rad = 0.0;
        m++;
        if (m >= 60) {
            m = 0;
            d++;
        }
    }
    if (sign > 0) {
	sprintf (str, "%2d %2d %5.2lf", d, m, rad);
    } else {
	sprintf (str, "-%2d %2d %5.2lf", d, m, rad);
    }
    return str;
}

int main (int argc, char *argv[])
{
    if (argc != 3) {
        printf ("Usage: alltimes hh:mm:ss mm/dd/yyyy\n");
        printf ("       time and date are UTC\n");
        return 1;
    }
    set_leap_sec_file_name ("leapsec.tab");
    set_iers_file_name ("iers.tab");
    double utc = str2days(argv[1]);
    int yr, mo, da, status;
    str2ymd(argv[2], yr, mo, da);

    double dmjd;
    ptwCaldj (yr, mo, da, &dmjd, &status);

    char str1[20], str2[20];		// places to hold the strings that
					// don't go out of scope
    long mjd = (long)(dmjd + 0.1);
    printf ("\n  At %s UTC, %s %d, %d  MJD %ld\n",
	    argv[1], mo2str(mo), da, yr, mjd);
    printf ("     Geocentric station coordinates:\n");
    printf ("          x,y,z (km) %9.3f %9.3f %9.3f\n",
	    get_observatory_x() / 1000.0, get_observatory_y() / 1000.0,
	    get_observatory_z() / 1000.0);
    printf ("          Long.  %s   Lat.  %s\n",
	    rad2dms_str(get_observatory_longitude(), str1),
	    rad2dms_str(get_observatory_geoc_latitude(), str2));
    double px, py;
    get_pole_offsets(mjd, utc, px, py);
    printf ("     Assumed polar motion offsets: x = %6.4f,  = %6.4f arcsec\n",
	    px, py);
    printf ("     Assumed UT1-UTC: %6.3f seconds\n",
	    get_ut1_offset(mjd, utc));

    printf ("                                                           MJD\n");
    printf ("  Coordinated Universal Time (UTC)......... %s  %ld\n",
	    day2time_str(utc, str1), mjd);

    double est = utc - 5.0 / 24.0;
    if (est < 0.0) est += 1.0;
    printf ("  Eastern Standard Time (EST).............. %s\n",
	    day2time_str(est, str1));

    long xmjd;
    double ut1;
    utc2ut1 (mjd, utc, xmjd, ut1);
    printf ("  Universal Time (UT1)..................... %s  %ld\n",
	    day2time_str(ut1, str1), xmjd);

    double ut0;
    utc2ut0 (mjd, utc, xmjd, ut0);
    printf ("  Universal Time (UT0)..................... %s  %ld\n",
	    day2time_str(ut0, str1), xmjd);

    double ut2;
    utc2ut2 (mjd, utc, xmjd, ut2);
    printf ("  Universal Time (UT2)..................... %s  %ld\n",
	    day2time_str(ut2, str1), xmjd);

    double tai;
    utc2tai (mjd, utc, xmjd, tai);
    printf ("  International Atomic Time (TAI).......... %s  %ld\n",
	    day2time_str(tai, str1), xmjd);

    double tt;
    utc2tt (mjd, utc, xmjd, tt);
    printf ("  Terrestrial Dynamic Time (TT)............ %s  %ld\n",
	    day2time_str(tt, str1), xmjd);

    double tdb;
    utc2tdb (mjd, utc, xmjd, tdb);
    printf ("  Barycentric Dynamic Time (TDB)........... %s  %ld\n",
	    day2time_str(tdb, str1), xmjd);

    double gmst = utc2gmst (mjd, utc);
    printf ("  Greenwich Mean Sidereal Time (GMST)...... %s\n",
	    day2time_str(gmst, str1));

    double gast = utc2gast (mjd, utc);
    printf ("  Greenwich Apparent Sidereal Time (GAST).. %s\n",
	    day2time_str(gast, str1));

    double lmst = utc2lmst (mjd, utc);
    printf ("  Local Mean Sidereal Time (LMST).......... %s\n",
	    day2time_str(lmst, str1));

    double last = utc2last (mjd, utc);
    printf ("  Local Apparent Sidereal Time (LAST)...... %s\n",
	    day2time_str(last, str1));
    printf("\n");

    return 0;
}
