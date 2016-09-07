// This file contains the functions that connect a python client to the JPL
// ephemeris functions and selected StarLink routines.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "jpl_eph.h"
#include "astrtime.h"
#include "delay.h"
#include "sofa.h"
#include "python_client.h"

// SS_OBJECT_DIRECTION
//      Returns topocentric astrometric J2000 coordinates for the specified
// object and time(s).  The astrometric position is the direction of the line
// between the observer's position at the specified time and the position of
// the object at the time the light began its journey, i.e. the object's
// position at the specified time less the light travel time.
// 
// 	Inputs:	Object name (we'll need a list of recognized names)
// 		MJD
// 		UTC in fraction of a day for the first position
// 		Number of positions wanted
// 		Interval between positions in UTC seconds
// 
// 	Output:	RA (fractional hours), Dec (fractional degrees),
// 		RA rate (arcmin/min), Dec rate (arcmin/min)
// 		for each time requested
//	Return: Number of positions found

int ss_object_direction ( char *object, long mjd, double utc,
			  int number_of_positions, double interval,
			  double *ra, double *dec,
			  double *ra_rate, double *dec_rate )
{
    int number_of_positions_done = 0;

    if (number_of_positions <= 0) {
	return 0;
    }
    long tdb_mjd;
    double tdb;
    long mjd_temp = mjd;
    double utc_temp = utc;
    int i;
    for (i = 0; i < number_of_positions; i++) {
	utc2tdb(mjd_temp, utc_temp, tdb_mjd, tdb);
	double jd = (double)tdb_mjd + 2400000.5;
	double pv[6];

	// get the astrometric topocentric position
	//  first get the light travel time
	double tptr[6];
	if (!get_object_posn_vel(tptr, object, Topocentric, jd, tdb))
	    return number_of_positions_done;
	double d = sqrt(tptr[0] * tptr[0] + tptr[1] * tptr[1] +
			tptr[2] * tptr[2]);

	// then recompute object for light travel advanced time
	double day_light = tdb - d / (C_LIGHT * 86400.0);
	double object_posn[6];
	if (!get_object_posn_vel(object_posn, object, Topocentric,
				 jd, day_light)) {
	    return number_of_positions_done;
	}
	// convert from km/day to km/min
	int iy;
	for (iy = 3; iy < 6; iy++) {
	    object_posn[iy] /= 1440.0;
	}
	double ra_rad, dec_rad, radius, ra_rate_rad, dec_rate_rad, rv;
	//slaDc62s (object_posn, &ra_rad, &dec_rad, &radius,
	//	  &ra_rate_rad, &dec_rate_rad, &rv);
	double pvx[2][3];
	pvx[0][0] = object_posn[0];
	pvx[0][1] = object_posn[1];
	pvx[0][2] = object_posn[2];
	pvx[1][0] = object_posn[3];
	pvx[1][1] = object_posn[4];
	pvx[1][2] = object_posn[5];
	iauPv2s(pvx, &ra_rad, &dec_rad, &radius, &ra_rate_rad, &dec_rate_rad,
		&rv);
	double rax = ra_rad * 24.0 / TWOPI;
	double decx = dec_rad * 360.0 / TWOPI;
	PositionT pos = ss_remove_aberration (mjd, utc, rax, decx);
	if (pos.ra < 0.0) pos.ra += 24.0;
	if (pos.ra >= 24.0) pos.ra -= 24.0;
	ra[number_of_positions_done] = pos.ra;
	dec[number_of_positions_done] = pos.dec;
	ra_rate[number_of_positions_done] =
					ra_rate_rad * 60.0 * 360.0 / TWOPI;
	dec_rate[number_of_positions_done] =
					dec_rate_rad * 60.0 * 360.0 / TWOPI;

	utc_temp += interval / 86400.0;
	if (utc_temp >= 1.0) {
	    mjd_temp++;
	    utc_temp -= 1.0;
	}
	number_of_positions_done++;
    }
    return number_of_positions_done;
}

int get_object_posn_vel( double pv[6], char *object, RefFrame refx,
			 double jd, double tdb )
{
    int i;
    for (i = 0; i < strlen(object); i++) {
	object[i] = toupper((int)object[i]);
    }
    if (!strcmp(object, "MOON")) {
	return get_moon_posn_vel (pv, refx, jd, tdb, 1);
    } else if (!strcmp(object, "SUN")) {
	return get_sun_posn_vel (pv, refx, jd, tdb, 1);
    } else {
	Object planet;
	if (!strcmp(object, "MERCURY")) {
	    planet = Mercury;
	} else if (!strcmp(object, "VENUS")) {
	    planet = Venus;
	} else if (!strcmp(object, "MARS")) {
	    planet = Mars;
	} else if (!strcmp(object, "JUPITER")) {
	    planet = Jupiter;
	} else if (!strcmp(object, "SATURN")) {
	    planet = Saturn;
	} else if (!strcmp(object, "URANUS")) {
	    planet = Uranus;
	} else if (!strcmp(object, "NEPTUNE")) {
	    planet = Neptune;
	} else if (!strcmp(object, "PLUTO")) {
	    planet = Pluto;
	} else {
	    return 0;
	}
	return get_planet_posn_vel (pv, planet, refx, jd, tdb, 1);
    }
}

// SS_DOPPLER_FRACTION
//      Returns the observer's radial velocity component with respect
// to the solar system barycenter for the specified J2000 coordinates and
// time(s).  The delay is in the sense that
//      barycentric velocity = observed velocity + "V / C" * C
// 
// 	Inputs:	RA(J2000), Dec(J2000) in fractional hours and degrees
// 		MJD
// 		UTC in fraction of a day for the first delay
// 		Number of values wanted
// 		Interval between values in UTC seconds
// 
// 	Output:	V / C for each time requested

int ss_doppler_fraction ( double ra2000, double dec2000, long mjd, double utc,
			  int number_of_values, double interval,
			  double *doppler_set )
{
    int number_done = 0;

    // get direction cosines
    double dcos[3];
    double ra_rad = ra2000 * TWOPI / 24.0;
    double dec_rad = dec2000 * TWOPI / 360.0;
    iauS2c ( ra_rad, dec_rad, dcos );

    double pv1[6], pv[6];
    long tt_mjd, tdb_mjd;
    double tt, tdb;
    long mjd_temp = mjd;
    double utc_temp = utc;
    int i;
    for (i = 0; i < number_of_values; i++) {
	utc2tt(mjd_temp, utc_temp, tt_mjd, tt);
	utc2tdb(mjd_temp, utc_temp, tdb_mjd, tdb);
	double jd_tt = (double)tt_mjd + 2400000.5;
	double jd = (double)tdb_mjd + 2400000.5;
	int status = get_earth_posn_vel(pv, SS_Barycentric, jd, tdb, 1);
	if (status == 0) {
	    return number_done;
	}
	get_observatory_posn_vel(pv1, jd_tt, tt, 1);
	double vel = 0.0;
	int j;
	for (j = 0; j < 3; j++) {
	    vel += (pv1[j + 3] + pv[j + 3]) * dcos[j] / 86400.0;
	}
	doppler_set[i] = vel / C_LIGHT;
	number_done++;
	utc_temp += interval / 86400.0;
	if (utc_temp >= 1.0) {
	    mjd_temp++;
	    utc_temp -= 1.0;
	}
    }
    return number_done;
}

// SS_OBSERVER_POSITION_VELOCITY
//      Returns the observer's position and velocity with respect
// to the solar system barycenter in J2000 rectangular coordinates.
// 
// 	Inputs:	MJD
// 		UTC in fraction of a day for the first delay
// 		Number of values wanted
// 		Interval between values in UTC seconds
// 
// 	Output:	x, y, z, dx/dt, dy/dt, dz/dt for each time requested
//		The units are kilometers and km/s.

int ss_observer_position_velocity ( long mjd, double utc,
				    int number_of_values, double interval,
				    double *x, double *y, double *z,
				    double *dx, double *dy, double *dz)
{
    int number_done = 0;

    double pv1[6], pv[6];
    long tt_mjd, tdb_mjd;
    double tt, tdb;
    long mjd_temp = mjd;
    double utc_temp = utc;
    int i;
    for (i = 0; i < number_of_values; i++) {
	utc2tt(mjd_temp, utc_temp, tt_mjd, tt);
	utc2tdb(mjd_temp, utc_temp, tdb_mjd, tdb);
	double jd_tt = (double)tt_mjd + 2400000.5;
	double jd = (double)tdb_mjd + 2400000.5;
	int status = get_earth_posn_vel(pv, SS_Barycentric, jd, tdb, 1);
	if (status == 0) {
	    return number_done;
	}
	get_observatory_posn_vel(pv1, jd_tt, tt, 1);

	x[i] = pv1[0] + pv[0];
	y[i] = pv1[1] + pv[1];
	z[i] = pv1[2] + pv[2];
	dx[i] = (pv1[3] + pv[3]) / 86400.0;
	dy[i] = (pv1[4] + pv[4]) / 86400.0;
	dz[i] = (pv1[5] + pv[5]) / 86400.0;

	number_done++;
	utc_temp += interval / 86400.0;
	if (utc_temp >= 1.0) {
	    mjd_temp++;
	    utc_temp -= 1.0;
	}
    }
    return number_done;
}

// SS_OBJECT_DOPPLER
//      Returns topocentric radial velocity component for the specified
// object and time(s).
// 
// 	Inputs:	Object name (we'll need a list of recognized names)
// 		MJD
// 		UTC in fraction of a day for the first velocity
// 		Number of velocities wanted
// 		Interval between velocities in UTC seconds
// 
// 	Output:	V / C for each time requested

int ss_object_doppler ( char *object, long mjd, double utc,
		       int number_of_velocities, double interval,
		       double *doppler_set )
{
    int number_done = 0;

    double pv[6], pv1[6];
    long tt_mjd, tdb_mjd;
    double tt, tdb;
    long mjd_temp = mjd;
    double utc_temp = utc;
    int i;
    for (i = 0; i < number_of_velocities; i++) {
	utc2tt(mjd_temp, utc_temp, tt_mjd, tt);
	utc2tdb(mjd_temp, utc_temp, tdb_mjd, tdb);
	double jd_tt = (double)tt_mjd + 2400000.5;
	double jd = (double)tdb_mjd + 2400000.5;
	// get the light travel distance
	int status = get_object_posn_vel(pv, object, Topocentric, jd, tdb);
	double d = sqrt(pv[0] * pv[0] + pv[1] * pv[1] +
			pv[2] * pv[2]);

	// get the observer's position and velocity for the time of
	// observation
	status = get_earth_posn_vel(pv, SS_Barycentric, jd, tdb, 1);
	get_observatory_posn_vel(pv1, jd_tt, tt, 1);
	int j;
	for (j = 0; j < 6; j++) {
	    pv[j] += pv1[j];
	}
	// get the object's position and velocity for the time that the
	// observed light began its journey
	double tdb_advanced = tdb - d / (C_LIGHT * 86400.0);
	status = get_object_posn_vel(pv1, object, SS_Barycentric, jd,
					 tdb_advanced);
	if (status == 0) {
	    return number_done;
	}

	// then use the line between the endpoints of the light's travel in
	// the barycentric frame for the direction, and use the object's
	// velocity at the beginning and the earth's velocity at the end for a
	// relative velocity vector.
	for (j = 0; j < 6; j++) {
	    pv[j] -= pv1[j];
	}
	d = sqrt(pv[0] * pv[0] + pv[1] * pv[1] + pv[2] * pv[2]);
	double vel = 0.0;
	for (j = 0; j < 3; j++) {
	    vel += pv[j + 3] * pv[j] / (d * 86400.0);
	}
	doppler_set[i] = vel / C_LIGHT;
	number_done++;
	utc_temp += interval / 86400.0;
	if (utc_temp >= 1.0) {
	    mjd_temp++;
	    utc_temp -= 1.0;
	}
    }
    return number_done;
}

// PULSE_DELAY
//      Returns the topocentric pulse delay and delay derivative with respect
// to the solar system barycenter for the specified J2000 coordinates and
// time(s).  Eventually, we will want another client service to return
// polynomial coefficients for the pulse frequency over a specified interval
// in a form useful to the spectral processor, but this can wait for a more
// detailed specification.  The delay is in the sense that
//      barycentric TOA = observed TOA + delay
// 
// 	Inputs:	RA(J2000), Dec(J2000) in fractional hours and degrees
// 		MJD
// 		UTC in fraction of a day for the first delay
// 		Number of delays wanted
// 		Interval between delays in UTC seconds
// 
// 	Output:	Delay in UTC seconds

int ss_pulse_delay ( double ra2000, double dec2000, long mjd, double utc,
		     int number_of_delays, double interval, double *delay_set )
{
    int number_done = 0;

    long mjd_temp = mjd;
    double utc_temp = utc;
    double ra_rad = ra2000 * TWOPI / 24.0;
    double dec_rad = dec2000 * TWOPI / 360.0;
    int i;
    for (i = 0; i < number_of_delays; i++) {
	delay_set[i] = -get_geometric_delay (SS_Barycentric, Topocentric,
					     mjd_temp, utc_temp,
					     ra_rad, dec_rad);
	number_done++;
	utc_temp += interval / 86400.0;
	if (utc_temp >= 1.0) {
	    mjd_temp++;
	    utc_temp -= 1.0;
	}
    }
    return number_done;
}

void ss_utc_to_tdb ( long mjd, double utc, long &tdb_mjd, double &tdb )
{
    utc2tdb (mjd, utc, tdb_mjd, tdb);
}

// SS_UTC_TO_LAST
//      Returns the local apparent sidereal time.
// 
// 	Inputs:	MJD
// 		UTC in fraction of a day
// 
// 	Output: LAST in fraction of a day

double ss_utc_to_last ( long mjd, double utc )
{
    return utc2last ( mjd, utc );
}

// SS_LAST_TO_UTC
//      Returns the UTC for the specified local apparent sidereal time.
// 
// 	Inputs: MJD
// 		LAST in fraction of a day
// 
// 	Output: UTC in fraction of a day.  Actually, a pointer to an array
//		of two double precisions numbers is returned. Since two UTC's
//		can satisfy the MJD/LAST combination for 4 minutes out of each
//		day, the second number is non-negative and contains the second
//		possibility if the ambiguity exists.

double *ss_last_to_utc ( long mjd, double last )
{
    static double utc[2];

    double gast = last - get_observatory_longitude() / TWOPI;
    if (gast < 0.0) {
        gast += 1.0;
    }
    if (gast >= 1.0) {
        gast -= 1.0;
    }
    double gast_zero = utc2gast (mjd, 0.0);
    double gast_diff = gast - gast_zero;
    if (gast_diff < 0.0) gast_diff += 1.0;

    utc[0] = gast_diff * 0.997269625;
    utc[1] = -1.0;
    if (gast_diff < 0.002737851) {
	utc[1] = (gast_diff + 1.0) * 0.997269625;
    }
    return utc;
}

// SS_J2000_TO_EPOCH
//      Returns the RA and Dec for the specified epoch and given J2000
// coordinates.  The correction is for precession and nutation.
// 
// 	Inputs: Modified Julian date
//		UTC in fraction of a day
// 		RA (J2000) in fractional hours
// 		Dec (J2000) in fractional degrees
// 
// 	Output: RA, Dec of epoch in fractional hours and degrees

PositionT ss_j2000_to_epoch ( long mjd, double utc,
			   double ra2000, double dec2000 )
{
    PositionT pos;

    long tt_mjd;
    double tt;
    utc2tt(mjd, utc, tt_mjd, tt);

    double date = (double)mjd + tt;

    double rmatpn[3][3], vect[3], new_vec[3], ra, dec;
    double ra_rad = ra2000 * TWOPI / 24.0;
    double dec_rad = dec2000 * TWOPI / 360.00;
    iauS2c (ra_rad, dec_rad, vect);
    //slaPrenut (2000.0, date, rmatpn);
    iauPnm00a (2400000.5, date, rmatpn);
    //slaDmxv (rmatpn, vect, new_vec);
    iauRxp (rmatpn, vect, new_vec);
    //slaDcc2s (new_vec, &ra, &dec);
    iauC2s (new_vec, &ra, &dec);

    pos.ra = ra * 24.0 / TWOPI;
    if (pos.ra < 0.0) pos.ra += 24.0;
    pos.dec = dec * 360.0 / TWOPI;
    pos.epoch = 2000.0 + (date - 51545.0) / 365.25;

    return pos;
}

// SS_EPOCH_TO_J2000
//      Returns the J2000 RA and Dec given the coordinates for the specified
// epoch.  The correction is for precession and nutation.
// 
// 	Inputs: Modified Julian date
//		UTC in fraction of a day
// 		RA (epoch) in fractional hours
// 		Dec (epoch) in fractional degrees
// 
// 	Output: RA(J2000), Dec(J2000) in fractional hours and degrees

PositionT ss_epoch_to_j2000 ( long mjd, double utc, double ra, double dec )
{
    PositionT pos;

    long tt_mjd;
    double tt;
    utc2tt(mjd, utc, tt_mjd, tt);

    double date = (double)mjd + tt;

    double rmatpn[3][3], vect[3], new_vec[3], ra2000, dec2000;
    double ra_rad = ra * TWOPI / 24.0;
    double dec_rad = dec * TWOPI / 360.00;
    iauS2c (ra_rad, dec_rad, vect);
    //slaPrenut (2000.0, date, rmatpn);
    iauPnm00a (2400000.5, date, rmatpn);
    //slaDimxv (rmatpn, vect, new_vec);
    iauTrxp (rmatpn, vect, new_vec);
    //slaDcc2s (new_vec, &ra2000, &dec2000);
    iauC2s (new_vec, &ra2000, &dec2000);

    pos.ra = ra2000 * 24.0 / TWOPI;
    if (pos.ra < 0.0) pos.ra += 24.0;
    pos.dec = dec2000 * 360.0 / TWOPI;
    pos.epoch = 2000.0;

    return pos;
}

// SS_ADD_ABERRATION
//      Returns the given geometric RA and Dec corrected to apparent RA and Dec
//      corrected for diurinal and annual aberration.
// 
// 	Inputs: MJD
// 		UTC in fraction of a day
// 		RA (geometric) in fractional hours
// 		Dec (geometric) in fractional degrees
// 
// 	Output: Apparent RA in fractional hours
// 		Apparent Dec in fractional degrees

PositionT ss_add_aberration ( long mjd, double utc, double ra, double dec )
{
    PositionT pos;

    double delta_ra, delta_dec;
    ss_aberration ( mjd, utc, ra, dec, delta_ra, delta_dec );

    pos.ra = ra + delta_ra;
    if (pos.ra >= 24.0) pos.ra -= 24.0;
    if (pos.ra < 0.0) pos.ra += 24.0;
    pos.dec = dec + delta_dec;
    pos.epoch = 2000.0;

    return pos;
}

// SS_REMOVE_ABERRATION
//      Returns the given geometric RA and Dec corrected to apparent RA and Dec
//      corrected for diurinal and annual aberration.
// 
// 	Inputs: MJD
// 		UTC in fraction of a day
//		Apparent RA in fractional hours
// 		Apparent Dec in fractional degree

// 	Output:	RA (geometric) in fractional hours
// 		Dec (geometric) in fractional degrees

PositionT ss_remove_aberration ( long mjd, double utc, double ra, double dec )
{
    PositionT pos;

    double delta_ra, delta_dec;
    ss_aberration ( mjd, utc, ra, dec, delta_ra, delta_dec );

    pos.ra = ra - delta_ra;
    if (pos.ra >= 24.0) pos.ra -= 24.0;
    if (pos.dec < 0.0) pos.ra += 24.0;
    pos.dec = dec - delta_dec;
    pos.epoch = 2000.0;

    return pos;
}

// SS_ABERRATION
//      Returns the aberration corrections given geometric RA and Dec.
// 
// 	Inputs: MJD
// 		UTC in fraction of a day
// 		RA (geometric) in fractional hours
// 		Dec (geometric) in fractional degrees
// 
// 	Output: RA correction in fractional hours
// 		Dec correction in fractional degrees

void ss_aberration ( long mjd, double utc, double ra, double dec,
		     double &delta_ra, double &delta_dec )
{
    delta_ra = delta_dec = 0.0;

    long tt_mjd, tdb_mjd;
    double tt, tdb;
    utc2tt(mjd, utc, tt_mjd, tt);
    utc2tdb(mjd, utc, tdb_mjd, tdb);

    double tpv[6];
    double jd_tt = (double)tt_mjd + 2400000.5;
    double jd = (double)tdb_mjd + 2400000.5;
    if (!get_observatory_posn_vel(tpv, jd_tt, tt, 1)) {
	printf ("aberration: error getting observatory's position\n");
	return;
    }
    double epv[6];
    if (!get_earth_posn_vel(epv, SS_Barycentric, jd, tdb, 1)) {
	printf ("aberration: error getting earth's position\n");
	return;
    }
    double v_o_c[3];
    int i;
    for (i = 0; i < 3; i++) {
	v_o_c[i] = (tpv[i + 3] + epv[i + 3]) / (86400.0 * C_LIGHT);
    }
    double ra_rad = TWOPI * ra / 24.0;
    double dec_rad = TWOPI * dec / 360.0;
    delta_ra = (-v_o_c[0] * sin(ra_rad) + v_o_c[1] * cos(ra_rad)) /
			cos(dec_rad);
    delta_dec = -v_o_c[0] * cos(ra_rad) * sin(dec_rad) -
			v_o_c[1] * sin(ra_rad) * sin(dec_rad) +
			 v_o_c[2] * cos(dec_rad);
    delta_ra *= 24.0 / TWOPI;
    delta_dec *= 360.0 / TWOPI;
}

// SS_SET_OBSERVER_COORDINATES
//	Sets the observer's (x, y, z) coordinates in meters with respect to
// the IERS earth coordinate system.  The default is the GBT.  These
// coordinates are used whenever the observer's location or velocity is
// required for a calculation.  The x axis is toward the equator/Greenwich
// meridian intersection, y is toward the equator 90 degrees east, and z is
// toward the north pole.

void ss_set_observer_coordinates ( double x, double y, double z )
{
     set_observatory_xyz ( x, y, z );
}

// GEOCENTRIC_OBSERVER_TRACK
//      Returns the observer's rectangular J2000 coordinates and velocities,
// in km and km/sec with respect to the center of the earth.
// 
// 	Inputs:	MJD
// 		UTC in fraction of a day for the first position
// 		Number of values wanted
// 		Interval between values in UTC seconds
// 
// 	Output:	x, y, z, x_rate, y_rate, and z_rate for each time requested

int ss_geocentric_observer_track ( long mjd, double utc,
				   int number_of_values, double interval,
				   double *obs_x, double *obs_y, double *obs_z,
				   double *obs_x_rate, double *obs_y_rate,
				   double *obs_z_rate)
{
    int number_done = 0;

    double pv[6];
    long tt_mjd;
    double tt;
    long mjd_temp = mjd;
    double utc_temp = utc;
    int i;
    for (i = 0; i < number_of_values; i++) {
	utc2tt(mjd_temp, utc_temp, tt_mjd, tt);
	double jd_tt = (double)tt_mjd + 2400000.5;
	int status = get_observatory_posn_vel(pv, jd_tt, tt, 1);
	if (status == 0) {
	    return number_done;
	}
	obs_x[i] = pv[0];
	obs_y[i] = pv[1];
	obs_z[i] = pv[2];
	obs_x_rate[i] = pv[3] / 86400.0;
	obs_y_rate[i] = pv[4] / 86400.0;
	obs_z_rate[i] = pv[5] / 86400.0;

	number_done++;
	utc_temp += interval / 86400.0;
	if (utc_temp >= 1.0) {
	    mjd_temp++;
	    utc_temp -= 1.0;
	}
    }
    return number_done;
}

// BARYCENTRIC_OBSERVER_TRACK
//      Returns the observer's rectangular J2000 coordinates and velocities,
// in km and km/sec with respect to the solar system barycenter.
// 
// 	Inputs:	MJD
// 		UTC in fraction of a day for the first position
// 		Number of values wanted
// 		Interval between values in UTC seconds
// 
// 	Output:	x, y, z, x_rate, y_rate, and z_rate for each time requested

int ss_barycentric_observer_track ( long mjd, double utc,
				    int number_of_values, double interval,
				    double *obs_x, double *obs_y,
				    double *obs_z,
				    double *obs_x_rate, double *obs_y_rate,
				    double *obs_z_rate)
{
    int number_done = 0;

    double pv[6], pv1[6];
    long tt_mjd, tdb_mjd;
    double tt, tdb;
    long mjd_temp = mjd;
    double utc_temp = utc;
    int i;
    for (i = 0; i < number_of_values; i++) {
	utc2tt(mjd_temp, utc_temp, tt_mjd, tt);
	utc2tdb(mjd_temp, utc_temp, tdb_mjd, tdb);
	double jd_tt = (double)tt_mjd + 2400000.5;
	double jd = (double)tdb_mjd + 2400000.5;
	int status = get_earth_posn_vel(pv1, SS_Barycentric, jd, tdb, 1);
	if (status == 0) {
	    return number_done;
	}
	get_observatory_posn_vel(pv, jd_tt, tt, 1);
	obs_x[i] = pv1[0] + pv[0];
	obs_y[i] = pv1[1] + pv[1];
	obs_z[i] = pv1[2] + pv[2];
	obs_x_rate[i] = (pv1[3] + pv[3]) / 86400.0;
	obs_y_rate[i] = (pv1[4] + pv[4]) / 86400.0;
	obs_z_rate[i] = (pv1[5] + pv[5]) / 86400.0;

	number_done++;
	utc_temp += interval / 86400.0;
	if (utc_temp >= 1.0) {
	    mjd_temp++;
	    utc_temp -= 1.0;
	}
    }
    return number_done;
}

// BARYCENTRIC_EARTH_TRACK
//      Returns the Earth's rectangular J2000 coordinates and velocities,
// in km and km/sec with respect to the solar system barycenter.
// 
// 	Inputs:	MJD
// 		UTC in fraction of a day for the first position
// 		Number of values wanted
// 		Interval between values in UTC seconds
// 
// 	Output:	x, y, z, x_rate, y_rate, and z_rate for each time requested

int ss_barycentric_earth_track ( long mjd, double utc,
				 int number_of_values, double interval,
				 double *earth_x, double *earth_y,
				 double *earth_z,
				 double *earth_x_rate, double *earth_y_rate,
				 double *earth_z_rate)
{
    int number_done = 0;

    double pv[6];
    long tdb_mjd;
    double tdb;
    long mjd_temp = mjd;
    double utc_temp = utc;
    int i;
    for (i = 0; i < number_of_values; i++) {
	utc2tdb(mjd_temp, utc_temp, tdb_mjd, tdb);
	double jd = (double)tdb_mjd + 2400000.5;
	int status = get_earth_posn_vel(pv, SS_Barycentric, jd, tdb, 1);
	if (status == 0) {
	    return number_done;
	}
	earth_x[i] = pv[0];
	earth_y[i] = pv[1];
	earth_z[i] = pv[2];
	earth_x_rate[i] = pv[3] / 86400.0;
	earth_y_rate[i] = pv[4] / 86400.0;
	earth_z_rate[i] = pv[5] / 86400.0;

	number_done++;
	utc_temp += interval / 86400.0;
	if (utc_temp >= 1.0) {
	    mjd_temp++;
	    utc_temp -= 1.0;
	}
    }
    return number_done;
}

// BARYCENTRIC_OBJECT_TRACK
//      Returns the specified object's rectangular J2000 coordinates and
// velocities,in km and km/sec with respect to the solar system barycenter.
// 
// 	Inputs:	Object name
//		MJD
// 		UTC in fraction of a day for the first position
// 		Number of values wanted
// 		Interval between values in UTC seconds
// 
// 	Output:	x, y, z, x_rate, y_rate, and z_rate for each time requested

int ss_barycentric_object_track ( char *object, long mjd, double utc,
				  int number_of_values, double interval,
				  double *x, double *y, double *z,
				  double *x_rate, double *y_rate,
				  double *z_rate)
{
    int number_done = 0;

    double pv[6];
    long tdb_mjd;
    double tdb;
    long mjd_temp = mjd;
    double utc_temp = utc;
    int i;
    for (i = 0; i < number_of_values; i++) {
	utc2tdb(mjd_temp, utc_temp, tdb_mjd, tdb);
	double jd = (double)tdb_mjd + 2400000.5;
	int status = get_object_posn_vel(pv, object, SS_Barycentric,
					 jd, tdb);
	if (status == 0) {
	    return number_done;
	}
	x[i] = pv[0];
	y[i] = pv[1];
	z[i] = pv[2];
	x_rate[i] = pv[3] / 86400.0;
	y_rate[i] = pv[4] / 86400.0;
	z_rate[i] = pv[5] / 86400.0;

	number_done++;
	utc_temp += interval / 86400.0;
	if (utc_temp >= 1.0) {
	    mjd_temp++;
	    utc_temp -= 1.0;
	}
    }
    return number_done;
}
