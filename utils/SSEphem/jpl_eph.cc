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
#include <string.h>
#include <math.h>
#include "jpl_eph.h"
#include "astrtime.h"

// These subroutines handle the interpretation of the JPL DE/LE planetary
// ephemerides and compute planetary positions, earth nutation, and lunar
// libration values.

de_headerT hdr;                           // file header structure
char       _file_name[120] = "";
FILE       *stream = NULL;                // ephemeris file stream

void set_ephemeris_file_path ( char *file_path )
{
    close_ephemeris();
    open_ephemeris(file_path);
}

int open_ephemeris (char *file_nm = 0)
{
    if (stream != NULL) {
        printf ("open_ephemeris: Ephemeris already open\n");
        return 0;
    }
    if (file_nm != 0) {
        strcpy(_file_name, file_nm);
    } else if (strlen(_file_name) == 0) {
        strcpy(_file_name, getenv("PYTHONPATH"));
	strcat(_file_name, "/JPLEphem/DEc403");
    }
    stream = fopen(_file_name, "rb");
    if (stream == NULL) {
        printf ("open_ephemeris: Cannot open file '%s'\n", _file_name);
        return 0;
    }
    fread(&hdr, sizeof(hdr), 1, stream);
    return 1;
}

void close_ephemeris()
{
    if (stream != NULL) {
        fclose(stream);
	stream = NULL;
    }
}

// This routine returns the geocentric position and velocity with respect to
// the selected reference frame for the specified time.  The return value is a
// six-element vector, three rectangular coordinates, in km, and three
// velocity components, in km/day, (x,y,z,xv,yv,zv), where the x axis is toward
// the J2000 equinox, and the z axis is normal to the J2000 equator.
// The requested time is carried in two double words for more precision.
// Ordinarily, the time origin will be an integer Julian date and the offset
// will contain the fractional day.  'time_offset' need not be between
// zero and one.  The velocity flag can be set to 0, if velocity information
// is not needed.

int get_earth_posn_vel (
		      double pv[6],	  // vector for return values
		      RefFrame ref_frame, // position-velocity reference frame
                      double time_org,    // Julian date time origin
                      double time_offset, // time offset from origin in days
                                          // epoch = time_org + time_offset
		      int    vel_flag)    // 0 if no velocity data needed
{
    if (ref_frame == Topocentric) {
	if (!get_observatory_posn_vel(pv, time_org, time_offset, vel_flag)) {
	    printf (
		"get_earth_posn_vel: error getting observatory's position\n");
	    return 0;
	}
	int i;
	for (i = 0; i < 6; i++) {
	    pv[i] = -pv[i];
	}
	return 1;
    }
    if (ref_frame == Geocentric) {
	printf (
         "get_earth_posn_vel: Geocentric earth position is a null value\n");
	int i;
	for (i = 0; i < 6; i++) {
	    pv[i] = 0.0;
	}
	return 1;
    }
    double geocentric_moon[6];
    if (!get_posn_vel (geocentric_moon, GeoCMoon, time_org, time_offset,
		       vel_flag)) {
	return 0;
    }
    if (ref_frame == EM_Barycentric) {
	int i;
        for (i = 0; i < 6; i++) {
            pv[i] = -geocentric_moon[i] / (1.0 + hdr.emrat);
        }
        return 1;
    }
    // get the barycentric earth from the earth-moon barycenter and the
    // geocentric moon
    if (!get_posn_vel (pv, EM_Bary, time_org, time_offset, vel_flag)) {
	return 0;
    }
    int i;
    for (i = 0; i < 6; i++) {
        pv[i] -= geocentric_moon[i] / (1.0 + hdr.emrat);
    }
    if (ref_frame == SS_Barycentric) {
        return 1;
    }
    // finally subtract sun's barycentric position, if requested.
    if (ref_frame == Heliocentric) {
	double sun_pv[6];
        if (!get_posn_vel (sun_pv, Sun, time_org, time_offset, vel_flag)) {
	    return 0;
	}
        for (i = 0; i < 6; i++) {
            pv[i] -= sun_pv[i];
        }
        return 1;
    }
    printf ("get_earth_posn_vel: Reference frame for earth unrecognized\n");
    return 0;     // an error if this point reached
}

// This routine returns the moon's position and velocity with respect to
// the selected reference frame for the specified time.  The return value is a
// six-element vector, three rectangular coordinates, in km, and three
// velocity components, in km/day, (x,y,z,xv,yv,zv), where the x axis is toward
// the J2000 equinox, and the z axis is normal to the J2000 equator.
// The requested time is carried in two double words for more precision.
// Ordinarily, the time origin will be an integer Julian date and the offset
// will contain the fractional day.  'time_offset' need not be between
// zero and one.  The velocity flag can be set to 0, if velocity information
// is not needed.

int get_moon_posn_vel (
		      double pv[6],	  // vector for return values
		      RefFrame ref_frame, // position-velocity reference frame
                      double time_org,    // Julian date time origin
                      double time_offset, // time offset from origin in days
                                          // epoch = time_org + time_offset
		      int    vel_flag)    // 0 if no velocity data needed
{
    if (!get_posn_vel (pv, GeoCMoon, time_org, time_offset, vel_flag)) {
	return 0;
    }
    if (ref_frame == Geocentric) {
        return 1;
    }
    double obs_pv[6];
    if (ref_frame == Topocentric) {
	if (!get_observatory_posn_vel(obs_pv, time_org, time_offset,
				      vel_flag)) {
	    printf (
		"get_moon_posn_vel: error getting observatory's position\n");
	    return 0;
	}
	int i;
	for (i = 0; i < 6; i++) {
	    pv[i] -= obs_pv[i];
	}
	return 1;
    }
    // compute earth-moon-barycentric moon position
    int i;
    for (i = 0; i < 6; i++) {
        pv[i] = hdr.emrat * pv[i] / (1.0 + hdr.emrat);
    }
    if (ref_frame == EM_Barycentric) {
        return 1;
    }
    // get the barycentric moon from the barycentric earth-moon barycenter
    // and the earth-moon-barycentric moon
    double emb_pv[6];
    if (!get_posn_vel (emb_pv, EM_Bary, time_org, time_offset, vel_flag)) {
	return 0;
    }
    for (i = 0; i < 6; i++) {
        pv[i] += emb_pv[i];
    }
    if (ref_frame == SS_Barycentric) {
        return 1;
    }
    // finally subtract the sun's position, if requested.
    if (ref_frame == Heliocentric) {
	double sun_pv[6];
        if (!get_posn_vel (sun_pv, Sun, time_org, time_offset, vel_flag)) {
	    return 0;
	}
        for (i = 0; i < 6; i++) {
            pv[i] -= sun_pv[i];
        }
        return 1;
    }
    printf ("get_moon_posn_vel: Reference frame for moon unrecognized\n");
    return 0;     // an error if this point reached
}

// This routine returns the sun's position and velocity with respect to
// the selected reference frame for the specified time.  The return value is a
// six-element vector, three rectangular coordinates, in km, and three
// velocity components, in km/day, (x,y,z,xv,yv,zv), where the x axis is toward
// the J2000 equinox, and the z axis is normal to the J2000 equator.
// The requested time is carried in two double words for more precision.
// Ordinarily, the time origin will be an integer Julian date and the offset
// will contain the fractional day.  'time_offset' need not be between
// zero and one.  The velocity flag can be set to 0, if velocity information
// is not needed.

int get_sun_posn_vel (
		      double pv[6],	  // vector for return values
		      RefFrame ref_frame, // position-velocity reference frame
                      double time_org,    // Julian date time origin
                      double time_offset, // time offset from origin in days
                                          // epoch = time_org + time_offset
		      int    vel_flag)    // 0 if no velocity data needed
{
    if (ref_frame == Heliocentric) {
	printf (
         "get_sun_posn_vel: Heliocentric sun position is a null value\n");
	int i;
	for (i = 0; i < 6; i++) {
	    pv[i] = 0.0;
	}
	return 1;
    }
    if (!get_posn_vel (pv, Sun, time_org, time_offset, vel_flag)) {
	return 0;
    }
    if (ref_frame == SS_Barycentric) {
        return 1;
    }
    // get the heliocentric earth-moon barycenter and use the reverse vector
    double emb_pv[6];
    if (!get_posn_vel (emb_pv, EM_Bary, time_org, time_offset, vel_flag)) {
	return 0;
    }
    int i;
    for (i = 0; i < 6; i++) {
        pv[i] -= emb_pv[i];
    }
    if (ref_frame == EM_Barycentric) {
        return 1;
    }
    // now move from earth-moon barycenter to geocenter
    double gc_moon_pv[6];
    if (!get_posn_vel (gc_moon_pv, GeoCMoon, time_org, time_offset, vel_flag)) {
	return 0;
    }
    for(i = 0; i < 6; i++) {
        pv[i] += gc_moon_pv[i] / (1.0 + hdr.emrat);
    }
    if (ref_frame == Geocentric) {
        return 1;
    }
    double obs_pv[6];
    if (ref_frame == Topocentric) {
	if (!get_observatory_posn_vel(obs_pv, time_org, time_offset,
				      vel_flag)) {
	    printf (
		"get_moon_posn_vel: error getting observatory's position\n");
	    return 0;
	}
	int i;
	for (i = 0; i < 6; i++) {
	    pv[i] -= obs_pv[i];
	}
	return 1;
    }
    printf ("get_sun_posn_vel: Reference frame for moon unrecognized\n");
    return 0;     // an error if this point reached
}

// This routine returns a planet's position and velocity with respect to
// the selected reference frame for the specified time.  The return value is a
// six-element vector, three rectangular coordinates, in km, and three
// velocity components, in km/day, (x,y,z,xv,yv,zv), where the x axis is toward
// the J2000 equinox, and the z axis is normal to the J2000 equator.
// The requested time is carried in two double words for more precision.
// Ordinarily, the time origin will be an integer Julian date and the offset
// will contain the fractional day.  'time_offset' need not be between
// zero and one.  The velocity flag can be set to 0, if velocity information
// is not needed.

int get_planet_posn_vel (
		      double pv[6],	  // vector for return values
                      Object planet,      // planet name
		      RefFrame ref_frame, // position-velocity reference frame
                      double time_org,    // Julian date time origin
                      double time_offset, // time offset from origin in days
                                          // epoch = time_org + time_offset
		      int    vel_flag)    // 0 if no velocity data needed
{
    if ((planet != Mercury) && (planet != Venus) && (planet != Mars) &&
        (planet != Jupiter) && (planet != Saturn) && (planet != Uranus) &&
        (planet != Neptune) && (planet != Pluto) && (planet != EM_Bary)) {
        printf ("get_planet_posn_vel: unrecognized planet name\n");
        return 0;
    }
    if (!get_posn_vel (pv, planet, time_org, time_offset, vel_flag)) {
	return 0;
    }
    if (ref_frame == SS_Barycentric) {
        return 1;
    }
    if (ref_frame == Heliocentric) {
	double sun_pv[6];
        if (!get_posn_vel (sun_pv, Sun, time_org, time_offset, vel_flag)) {
	    return 0;
	}
	int i;
        for (i = 0; i < 6; i++) {
            pv[i] -= sun_pv[i];
        }
        return 1;
    }
    // move to earth-moon barycenter
    double emb_pv[6];
    if (!get_posn_vel (emb_pv, EM_Bary, time_org, time_offset, vel_flag)) {
	return 0;
    }
    int i;
    for (i = 0; i < 6; i++) {
        pv[i] -= emb_pv[i];
    }
    if (ref_frame == EM_Barycentric) {
        return 1;
    }
    // move to geocenter
    double gc_moon_pv[6];
    if (!get_posn_vel (gc_moon_pv, GeoCMoon, time_org, time_offset, vel_flag)) {
	return 0;
    }
    for (i = 0; i < 6; i++) {
        pv[i] += gc_moon_pv[i] / (1.0 + hdr.emrat);
    }
    if (ref_frame == Geocentric) {
        return 1;
    }
    double obs_pv[6];
    if (ref_frame == Topocentric) {
	if (!get_observatory_posn_vel(obs_pv, time_org, time_offset,
				      vel_flag)) {
	    printf (
		"get_moon_posn_vel: error getting observatory's position\n");
	    return 0;
	}
	int i;
	for (i = 0; i < 6; i++) {
	    pv[i] -= obs_pv[i];
	}
	return 1;
    }
    printf ("get_planet_posn_vel: Reference frame for planet unrecognized\n");
    return 0;     // an error if this point reached
}

// This routine figures out where in the ephemeris file are the coefficients
// for the specified object and time, passes the coefficients and time to
// the interpolation routine, and returns the interpolated three-dimensional
// position and velocity.  The return value is a six-element vector, three
// barycentric rectangular coordinates, in km, and three velocity components,
// in km/day, (x,y,z,xv,yv,zv), where the x axis is toward the J2000 equinox,
// and the z axis is normal to the J2000 equator.  The file must
// have been opened and header loaded before calling this routine
// The requested time is carried in two double words for more precision.
// Ordinarily, the time origin will be an integer Julian date and the offset
// will contain the fractional day.  'time_offset' need not be between
// zero and one.  The velocity flag can be set to 0, if velocity information
// is not needed.

int get_posn_vel (double pv[6],		// vector for return values
		  Object object,	// defined by Object enum
		  double time_org,	// Julian date time origin
		  double time_offset,	// time offset from origin in days
                                        // epoch = time_org + time_offset
		  int    vel_flag)	// 0 if no velocity data needed
{
    if ((object == Nutation) || (object == Libration)) {
	printf ("get_posn_vel: Nutation and libration not handled by this");
	printf ("    routine.  Use get_nutation() or get_libration.\n");
	return 0;
    }
    // get the full interval record buffer from the coefficients file for this
    // time for all objects.
    double *coeff_buffer = get_coeff_buffer(time_org, time_offset);
    if (coeff_buffer == 0) {
	return 0;
    }
    // get the coefficients for the specified object and subinterval range,
    // the fractional subinterval time fraction, the number of coefficients,
    // and the subinterval duration offset for use by the interpolation
    // routine.  Note the the coefficients list contains three sets of
    // coefficients, one for each (x,y,z) component.
    double subint_frac;
    int num_coeff;
    double sub_interval;
    double *coeff = get_coeff_subinterval(object, coeff_buffer, time_org,
					  time_offset, subint_frac, num_coeff,
					  sub_interval);
    if (coeff == 0) return 0;

    int i;
    for (i = 0; i < 3; i++) {
        double *pvx = interp (coeff, num_coeff, vel_flag, sub_interval,
                                                               subint_frac);
        if (pvx == 0) return 0;
        pv[i]     = pvx[0];
        pv[i + 3] = pvx[1];
        coeff += num_coeff;
    }
    return 1;
}

// This routine figures out where in the ephemeris file are the nutation
// coefficients for the specified time, passes the coefficients and time to
// the interpolation routine, and returns the interpolated two-dimensional
// offsets and velocity.  The return value is a four-element vector, the
// longitude and obliquity in radians of the nutated pole and a velocity
// component for each in radians / day.  The file must have been opened and
// header loaded before calling this routine.
// The requested time is carried in two double words for more precision.
// Ordinarily, the time origin will be an integer Julian date and the offset
// will contain the fractional day.  'time_offset' need not be between
// zero and one.  The velocity flag can be set to 0, if velocity information
// is not needed.

int get_nutation (double pv[4],		// vector for return values
		  double time_org,	// Julian date time origin
		  double time_offset,	// time offset from origin in days
                                        // epoch = time_org + time_offset
		  int    vel_flag)	// 0 if no velocity data needed
{
    // get the full interval record buffer from the coefficients file for this
    // time for all objects.
    double *coeff_buffer = get_coeff_buffer(time_org, time_offset);
    if (coeff_buffer == 0) return 0;

    // get the coefficients for nutation in the subinterval range,
    // the fractional subinterval time fraction, the number of coefficients,
    // and the subinterval duration offset for use by the interpolation
    // routine.  Note the the coefficients list contains three sets of
    // coefficients, one for each (x,y,z) component.
    double subint_frac;
    int num_coeff;
    double sub_interval;
    double *coeff = get_coeff_subinterval(Nutation, coeff_buffer, time_org,
					  time_offset, subint_frac, num_coeff,
					  sub_interval);
    if (coeff == 0) return 0;

    int i;
    for (i = 0; i < 2; i++) {
        double *pvx = interp (coeff, num_coeff, vel_flag, sub_interval,
                                                               subint_frac);
        if (pvx == 0) return 0;
        pv[i]     = pvx[0];
        pv[i + 2] = pvx[1];
        coeff += num_coeff;
    }
    return 1;
}

int get_libration(double pv[6])		  // vector for return values
{
    printf ("Libration retrieval not implemented\n");
    int i;
    for (i = 0; i < 6; i++) {
	pv[i] = 0.0;
    }
    return 0;
}

double *get_coeff_subinterval (Object object, double *coeff_buffer,
			       double time_org, double time_offset,
			       double &subint_frac, int &num_coeff,
			       double &sub_interval)
{
    int index = (int)object;

    // first point to the beginning of this object's area
    static double *sub_coeff;
    sub_coeff = coeff_buffer + hdr.obj[index].offset;

    // figure out which subinterval the time is in
    double start_offset = time_org - coeff_buffer[0];
    start_offset += time_offset;
    double time_frac = start_offset / hdr.interval;
    double dnum_sub = (double)hdr.obj[index].num_subinterv;
    // truncate to catch time_frac = 1.0
    double dt1 = (double)((int)time_frac);
    // compute subinterval index
    int si = (int)(dnum_sub * time_frac - dt1);

    // move the pointer to the subinterval coefficients
    num_coeff = hdr.obj[index].num_coeff;
    sub_coeff += num_coeff * hdr.obj[index].num_components * si;

    // compute the subinterval fractional time
    sub_interval = hdr.interval / dnum_sub;
    subint_frac = (start_offset - (double)si * sub_interval) / sub_interval;
    if ((subint_frac < 0.0) || (subint_frac > 1.0)) {
	printf (
       "get_coeff_subinterval: error computing subinterval fraction %18.16f\n",
	    subint_frac);
	return 0;
    }
    return sub_coeff;
}

// This routine creates memory for the coefficients list, keeps track of the
// time range of the list currently in memory, and reads a new record, if the
// requested time isn't in the current range.  The returned pointer is to the
// cofficients list.  This list covers a time range given by hdr.interval.  The
// list contains coefficients for all objects and may contain coefficient sets
// for subintervals for each object.  Objects can have a different numbers
// of subintervals and numbers of coefficients in a subinterval, but the number
// of coefficients is the same for all subintervals for one object.  The file
// layout is
//  File record:
//     Object_1, Object_2, Object_3,.... Object_N
//                                        where N = numObj in enum Object
//       The beginning of the coefficients for Objectn in the full set is
//       given by hdr.obj[n].offset as the number of double words.
//  Objectn_:
//     Subinterval_1, Subinterval_2,.... Subinterval_M
//                                        where M = hdr.obj[n].num_subinterv
//  Subinterval_m:
//     CoeffList_1, CoeffList_2, [CoeffList_3]
//       Each Subinterval contains either 2 or 3 sets of double precision
//       coefficients, one for each rectangular coordinate component (3) or
//       each polar offset component (2).  The number of coefficients in each
//       CoeffList is hdr.obj[n].num_coeff.  The time interval for each
//       subinterval is hdr.interval / hdr.obj[n].num_subinterv.

double *get_coeff_buffer(double time_org, double time_offset)
{
    if (stream == 0) {
        int stat = open_ephemeris();
	if (stat == 0) return 0;
    }
    // save the time range currently in memory to reduce repeated file
    // reads of the same record
    static double rec_epoch1 = 1.0e10;
    static double rec_epoch2 = -1.0e10;
    static double *coeff_buffer = 0;
    if (coeff_buffer == 0) {
        coeff_buffer = new double[hdr.num_coeff];
    }
    // compute the time offset in days from the beginning of the record in
    // memory
    double start_offset = time_org - rec_epoch1;
    start_offset += time_offset;
    // read a different record only if the current record doesn't cover the
    // requested time
    if ((start_offset < 0.0) || ((rec_epoch1 + start_offset) > rec_epoch2)) {
	// hdr.epoch1 and hdr.epoch2 are the start and end epochs of the full
	// ephemeris file
        start_offset = (time_org - hdr.epoch1) * 1.0038;
        start_offset += time_offset;
        if ((start_offset < 0.0) ||
                      ((hdr.epoch1 + start_offset) > hdr.epoch2)) {
            printf ("get_posn_vel: requested time outside ephemeris range\n");
            return 0;
        }
	// start by assuming that the coefficient records lie end to end in the
	// file with no overlap
        long record_number = (long)(start_offset / hdr.interval);
	long rec_len = hdr.num_coeff * sizeof(double);
        long foffset = (long)sizeof(hdr) + record_number * rec_len;
        if (fseek (stream, foffset, SEEK_SET)) {
            printf ("get_posn_vel: file seek error\n");
            return 0;
        }
        fread(coeff_buffer, sizeof(double), hdr.num_coeff, stream);
	// make sure we have the right record in case our no-overlap assumption
	// was wrong
        start_offset = time_org - coeff_buffer[0];
        start_offset += time_offset;
	int safety = 20;
	while ((start_offset < 0.0) && (foffset >= (sizeof(hdr) + rec_len)) &&
					(--safety > 0)) {
	    foffset -= rec_len;
	    if (fseek (stream, foffset, SEEK_SET)) {
		printf ("get_posn_vel: file seek error (backward search)\n");
		return 0;
	    }
	    fread(coeff_buffer, sizeof(double), hdr.num_coeff, stream);
	    start_offset = time_org - coeff_buffer[0];
	    start_offset += time_offset;
	}
	// or maybe there was a gap between records in the file that made us
	// miss the right record
	safety = 20;
	double end_offset = time_org - coeff_buffer[1];
	end_offset += time_offset;
	while((end_offset > 0.0) && (--safety > 0)) {
	    foffset += rec_len;
	    if (fseek (stream, foffset, SEEK_SET)) {
		printf ("get_posn_vel: file seek error (forward search)\n");
		return 0;
	    }
	    fread(coeff_buffer, sizeof(double), hdr.num_coeff, stream);
	    end_offset = time_org - coeff_buffer[1];
	    end_offset += time_offset;
	}
        start_offset = time_org - coeff_buffer[0];
        start_offset += time_offset;
        rec_epoch1 = coeff_buffer[0];
        rec_epoch2 = coeff_buffer[1];
	//printf("start_offset: %f\n", start_offset);
	//printf("rec_epoch1: %f\n", rec_epoch1);
	//printf("rec_epoch2: %f\n", rec_epoch2);
	// alas, might we have failed after all?
	if ((start_offset < 0.0) ||
	    ((rec_epoch1 + start_offset) > rec_epoch2)) {
	    printf ("get_posn_vel: time range not found in ephemeris\n");
	    return 0;
	}
    }
    return coeff_buffer;
    
}

// This routine computes a position vector component and vector component
// derivative from a given set of coefficients and the relative time within
// the interval covered by those coefficients.  The returned value is a
// two-element, double precision number (position and velocity) interpolated
// from the coefficients of the chebyshev polynomial equation.  The units are
// kilometers and kilometers / day, except for nutation and libration, which
// are in radians and radians / day.

// The Chebyshev polynomial is the truncated series
//           y = c[n] * T[n](x)      n = 0 -> N-1
// where
//           T[0] = 1,  T[1] = x,  T[2] = 2 * x^2 - 1,
// and then
//           T[n] = 2 * x * T[n-1] - T[n-2]
// The values c[n] are the coefficients stored in the ephemerides for
// each object, subinterval, and position component.  See "Numerical Recipes"
// in FORTRAN, page 147.

double *interp (double *coeff,      // coefficient list
                int    num_coeff,   // number of coefficients in the list
                int    vel_flag,    // 0 = no velocity, 1 = velocity returned
                double interval,    // interval covered by the coefficients
                                    // (full interval / num subintervals)
                double time)        // fractional time in coefficients
                                    // interval (0.0 <= time <= 1.0)
{
    if (time < 0.0 || time > 1.0) {
        printf("interp: time out of range\n");
        return 0;
    }
    // tc is the normalized chebyshev time (-1 <= tc <= 1)
    double tc = 2.0 * time - 1.0;

    // evaluate 'num_coef' polynomial values
    double pc[18];
    pc[0] = 1.0;
    double twot = tc + tc;
    pc[1] = tc;
    int i;
    for (i = 2; i < num_coeff; i++) {
        pc[i] = twot * pc[i - 1] - pc[i - 2];
    }
    // sum the polynomial to get position for each component
    // the sum is done in reverse order to preserve accuracy of smaller terms
    static double val[2];
    val[0] = val[1] = 0.0;
    for (i = num_coeff - 1; i >= 0; i--) {
        val[0] = val[0] + pc[i] * coeff[i];
    }
    if (vel_flag < 1) {
        return val;
    }
    // if velocity interpolation is wanted, generate derivative polynomials
    double vc[18];
    vc[0] = 0.0;
    vc[1] = 1.0;
    vc[2] = twot + twot;
    for (i = 3; i < num_coeff; i++) {
        vc[i] = twot * vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2];
    }
    // sum the first derivative of the series to get the velocity
    for (i = num_coeff - 1; i >= 1; i--) {
        val[1] = val[1] + vc[i] * coeff[i];
    }
    val[1] = 2.0 * val[1] / interval;
    return val;
}

double get_au()
{
    return hdr.au;
}
