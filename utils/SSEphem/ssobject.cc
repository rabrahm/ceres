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

//   SSOBJECT
// This program returns the topocentric J2000 (Astrometric) position
// (RA, Dec), position derivatives, and radial velocity of the
// specified major solar system object.  This position is corrected
// for light travel time in the J2000 coordinate system, but not for
// precession, nutation, or aberration.  Hence, the position returned
// may be used to locate the object as it would appear on a background
// of stars on a J2000 coordinate grid.

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "sofa.h"
#include "astrtime.h"
#include "jpl_eph.h"
#include "ptw_aux.h"

int set_obs_location();
long get_mjd(char *date);
double get_utc (char *str);
int get_object_posn_vel(double tptr[6], char *object, RefFrame ref,
			double jd, double tt);
char *decimal2str(double d, int num_decimal, char *str);
void make_lower(char *str);

int main (int argc, char *argv[])
{
    if (argc < 4) {
	printf("\n Usage: ssobject object MJD|mm/dd/yyyy UTC\n");
	printf("    object = sun     | moon   | mercury | venus   | mars |\n");
	printf("             jupiter | saturn | uranus  | neptune | pluto\n");
	printf(" Examples:\n");
	printf("            ssobject moon 50086 23:13:26.2\n");
	printf("            ssobject mars 12/10/1995 08:11:53\n");
	printf("\n This program returns the topocentric J2000 (Astrometric)");
	printf(" position\n (RA, Dec), position derivatives, and radial ");
	printf("velocity of the\n specified major solar system object.  ");
	printf("This position is corrected\n for light travel time in the ");
	printf("J2000 coordinate system, but not for\n precession, nutation,");
	printf(" or aberration.  Hence, the position returned\n may be used ");
	printf("to locate the object as it would appear on a background\n ");
	printf("of stars on a J2000 coordinate grid.\n ");
	printf("The location of the 140-ft telescope is assumed, but this ");
	printf("may be\n overridden with the geodetic coordinate environment");
	printf(" variables:\n   export OBS_LONG=-79.123456  (east longitude ");
	printf("in degrees)\n   export OBS_LAT=38.654321\n   export OBS_HGHT");
	printf("=830.0       (height above geoid in meters)\n");
	printf(" The ephemeris, leap second, and UT1 files path must be set,");
	printf(" e.g.,\n");
	printf("   export EPH_PATH=/hyades1/rfisher/SSEphem/DE200\n\n");
	return 1;
    }
    if (set_obs_location()) {
	return 1;
    }
    char object_lc[20];
    strcpy(object_lc, argv[1]);
    make_lower(object_lc);
    long mjd = get_mjd(argv[2]);
    if (mjd < 49000 || mjd > 51000) {
	printf ("Unrecognized date or MJD, %ld\n", mjd);
	printf ("  Must be either mm/dd/yyyy or MJD, e.g., 50086 format\n");
	return 1;
    }
    double utc = get_utc(argv[3]);
    if (utc < 0.0 || utc > 1.0) {
	printf ("Unrecognized UTC, %lf\n", utc);
	printf ("  Must be in hh:mm:ss.s format\n");
	return 1;
    }
    char *eph_path;
    if ((eph_path = getenv("EPH_PATH")) == NULL) {
	printf ("Ephemeris, leap second, and UT1 files path not specified");
	printf (" with 'export EPH_PATH'\n");
	return 1;
    }
    char file_path[120];
    strcpy(file_path, eph_path);
    strcat(file_path, "/DEc200");
    if (!open_ephemeris(file_path)) {
        return 1;
    }
    strcpy(file_path, eph_path);
    strcat(file_path, "/leapsec.tab");
    if (!set_leap_sec_file_name (file_path)) {
	return 1;
    }
    strcpy(file_path, eph_path);
    strcat(file_path, "/iers.tab");
    if (!set_iers_file_name (file_path)) {
        return 1;
    }
    
    // convert UTC to Barycentric Dynamic Time used by the DE200 solar
    // system ephemeris
    long tdb_mjd;
    double tdb;
    utc2tdb(mjd, utc, tdb_mjd, tdb);
    double jd = (double)tdb_mjd + 2400000.5;	// Julian date at 0h TDB
printf ("JD = %lf, tdb = %12.9lf\n", jd, tdb);

    // first get the distance to the object which will give us the
    // light travel time
    double object_posn[6];	// x,y,z,xv,yv,zv
    if (!get_object_posn_vel(object_posn, object_lc, Topocentric, jd, tdb))
	return 1;
    // distance in km
    double d = sqrt(object_posn[0] * object_posn[0] +
		    object_posn[1] * object_posn[1] +
		    object_posn[2] * object_posn[2]);

    // then compute the barycentric position and velocity of the object at
    // the light-travel-advanced time
    RefFrame interm_ref = SS_Barycentric;
    if (!strcmp(object_lc, "moon")) {
	interm_ref = EM_Barycentric;	// more accurate for the moon
    }
    double tdb_light = tdb - d / (C_LIGHT * 86400.0);
    if (!get_object_posn_vel(object_posn, object_lc, interm_ref,
			     jd, tdb_light))
	return 1;

    // but compute the observer's position for the given time and subtract
    // from the object's advanced position
    double earth_posn[6];
    if (!get_earth_posn_vel(earth_posn, interm_ref, jd, tdb))
	return 1;
    double observer_posn[6];
    if(!get_observatory_posn_vel(observer_posn, jd, tdb))
	return 1;
    int i;
    for (i = 0; i < 6; i++) {
	object_posn[i] -= earth_posn[i];// + observer_posn[i];
    }

    // now we can convert the 3-D position and velocity, relative to the
    // observer, into sky coordinates and radial velocity
    double ra, dec, dist, ra_vel, dec_vel, radial_vel;
    //slaDc62s(object_posn, &ra, &dec, &dist, &ra_vel, &dec_vel, &radial_vel);
    //iauPv2s(pv, &theta, &phi, &rx, &td, &pd, &rdx);
    double pv[2][3];
    pv[0][0] = object_posn[0];
    pv[0][1] = object_posn[1];
    pv[0][2] = object_posn[2];
    pv[1][0] = object_posn[3];
    pv[1][1] = object_posn[4];
    pv[1][2] = object_posn[5];
    iauPv2s(pv, &ra, &dec, &dist, &ra_vel, &dec_vel, &radial_vel);

    // convert to hours, degrees, km, arcsec/sec, arcsec/sec, and km/sec
    ra *= RAD2DEG / 15.0;
    if (ra < 0.0) ra += 24.0;
    dec *= RAD2DEG;
    ra_vel *= RAD2DEG / 24.0;
    dec_vel *= RAD2DEG / 24.0;
    radial_vel /= 86400.0;
    char str1[40], str2[40];
    printf ("RA %s  Dec %s  Dist %1.0lf km\n",
	    decimal2str(ra, 4, str1), decimal2str(dec, 3, str2), dist);
    printf ("RAv %9.7lf \"/s  Decv %9.7lf \"/s  Rv %5.3lf km/s\n",
	    ra_vel, dec_vel, radial_vel);
    return 0;
}

int set_obs_location()
{
    char *obs_long = getenv("OBS_LONG");
    char *obs_lat  = getenv("OBS_LAT");
    char *obs_hght = getenv("OBS_HGHT");
    if (obs_long == NULL && obs_lat == NULL && obs_hght == NULL) {
	// use the default observatory location
	return 0;
    } else if ((obs_long != NULL && obs_lat != NULL && obs_hght != NULL) &&
	       (strlen(obs_long) > 1 && strlen(obs_lat) > 1 &&
		strlen(obs_hght) > 0)) {
	double lon = atof(obs_long);
	if (lon < -360.0 || lon > 360.0) {
	    printf (
	      "Observatory longitude, %lf, out of range (-360 to 360)\n", lon);
	    return 1;
	}
	double lat = atof(obs_lat);
	if (lat < -360.0 || lat > 360.0) {
	    printf (
		"Observatory latitude, %lf, out of range (-90 to 90)\n", lat);
	    return 1;
	}
	double hght = atof(obs_hght);
	if (hght > 10000.0 || hght < -1000.0) {
	    printf (
	  "Observatory height, %lf, out of range (-1,000 to 10,000 meters)\n",
		hght);
	    return 1;
	}
	set_latlong2xyz(lon, lat, hght);
    } else {
	printf ("You must specify no or all three environment variables\n");
	printf ("  OBS_LONG  OBS_LAT  OBS_HGHT\n");
	return 1;
    }
    return 0;
}

long get_mjd(char *date)
{
    if (strchr(date, (int)'/') == NULL) {
	return atol(date);
    }
    char wstr[20];
    strcpy (wstr, date);
    // convert delimiters to spaces
    int i;
    for (i = 0; i < strlen(wstr); i++) {
        if (wstr[i] == '/') {
            wstr[i] = ' ';
        }
    }
    int mo, da, yr;
    if (sscanf (wstr, "%d %d %d", &mo, &da, &yr) != 3) {
        printf ("Error converting time string '%s'\n", date);
	return 0;
    }
    double jdm;
    int status;
    ptwCaldj(yr, mo, da, &jdm, &status);
    switch (status) {
	case 1:
	    printf("Bad year in date, %s\n", date);
	    return 0;
	case 2:
	    printf("Bad month in date, %s\n", date);
	    return 0;
	case 3:
	    printf("Bad day in date, %s\n", date);
	    return 0;
	default:
	    return (long)(jdm + 0.01);
    }
}

double get_utc (char *str)
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
        return -1.0;
    }
    return hrs / 24.0 + min / 1440.0 + sec /86400.0;
}

int get_object_posn_vel(double tptr[6], char *object, RefFrame ref,
			double jd, double tdb)
{
    if (!strcmp(object, "sun")) {
	return get_sun_posn_vel(tptr, ref, jd, tdb);
    } else if (!strcmp(object, "moon")) {
	return get_moon_posn_vel(tptr, ref, jd, tdb);
    } else if (!strcmp(object, "mercury")) {
	return get_planet_posn_vel(tptr, Mercury, ref, jd, tdb);
    } else if (!strcmp(object, "venus")) {
	return get_planet_posn_vel(tptr, Venus, ref, jd, tdb);
    } else if (!strcmp(object, "mars")) {
	return get_planet_posn_vel(tptr, Mars, ref, jd, tdb);
    } else if (!strcmp(object, "jupiter")) {
	return get_planet_posn_vel(tptr, Jupiter, ref, jd, tdb);
    } else if (!strcmp(object, "saturn")) {
	return get_planet_posn_vel(tptr, Saturn, ref, jd, tdb);
    } else if (!strcmp(object, "uranus")) {
	return get_planet_posn_vel(tptr, Uranus, ref, jd, tdb);
    } else if (!strcmp(object, "neptune")) {
	return get_planet_posn_vel(tptr, Neptune, ref, jd, tdb);
    } else if (!strcmp(object, "pluto")) {
	return get_planet_posn_vel(tptr, Pluto, ref, jd, tdb);
    } else {
	printf ("Object name not recognized, %s\n", object);
	return 0;
    }
}

char *decimal2str(double d, int num_decimal, char *str)
{
    str[0] = '\0';
    if (d < 0.0) {
	strcpy(str, "-");
    }
    d = fabs(d);
    int hd = (int)d;
    d -= (double)hd;
    d *= 60.0;
    int min = (int)d;
    d -= (double)min;
    double sec = d * 60.0;

    double tol = 0.5;
    int i;
    for (i = 0; i < num_decimal; i++) {
	tol *= 0.1;
    }
    if (sec > (60.0 - tol)) {
	sec = 0.0;
	min++;
	if (min >= 60) {
	    min = 0;
	    hd++;
	}
    }
    char xstr[40];
    sprintf (xstr, "%02d:%02d:%012.9lf", hd, min, sec);
    if (num_decimal > 9) {
	num_decimal = 9;
    } else if (num_decimal < 0) {
	num_decimal = 0;
    }
    char *ptr = strchr(xstr, '.');
    if (num_decimal == 0) {
	*ptr = '\0';
    } else if (strlen(ptr) > num_decimal+1) {
	ptr[num_decimal+1] = '\0';
    }
    strcat(str, xstr);
    return str;
}

void make_lower(char *str)
{
    int i;
    for (i = 0; i < strlen(str); i++) {
	str[i] = tolower(str[i]);
    }
}
