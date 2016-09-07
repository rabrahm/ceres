#include "Python.h"
#include "arrayobject.h"
#include "jplephem.h"
#include "jpl_eph.h"
#include "delay.h"
#include "astrtime.h"
#include "python_client.h"
#include <stdio.h>

extern void set_ephemeris_dir ( char *path, char *file_name )
{
    char leap_sec_file_path[200];
    char iers_file_path[200];
    char ephem_file_path[200];
    strcpy(leap_sec_file_path, path);
    strcat(leap_sec_file_path, "/leapsec.tab");
    strcpy(iers_file_path, path);
    strcat(iers_file_path, "/iers.tab");
    strcpy(ephem_file_path, path);
    strcat(ephem_file_path, "/");
    strcat(ephem_file_path, file_name);
    set_leap_sec_file_name (leap_sec_file_path);
    set_iers_file_name (iers_file_path);
    set_ephemeris_file_path(ephem_file_path);
}

extern PyObject *object_track ( char *object, long mjd, double utc,
				int number_of_positions,
				double interval )
{
    static int np = 0;
    // return values
    // RA,dec (radians) with e-terms included
    static double *ra;
    static double *dec;
    static double *ra_rate;
    static double *dec_rate;
    if (np != number_of_positions) {
        delete ra;
	delete dec;
	delete ra_rate;
	delete dec_rate;
	np = number_of_positions;
    }
    ra = new double[np];
    dec = new double[np];
    ra_rate = new double[np];
    dec_rate = new double[np];

    int stat = ss_object_direction (object, mjd ,utc, np,
				    interval,
				    ra, dec, ra_rate, dec_rate);
    PyObject *ra_list = PyList_New(np);
    PyObject *dec_list = PyList_New(np);
    PyObject *ra_rate_list = PyList_New(np);
    PyObject *dec_rate_list = PyList_New(np);

    int i;
    for (i = 0; i < np; i++) {
        PyList_SetItem(ra_list, i, PyFloat_FromDouble(ra[i]));
        PyList_SetItem(dec_list, i, PyFloat_FromDouble(dec[i]));
        PyList_SetItem(ra_rate_list, i, PyFloat_FromDouble(ra_rate[i]));
        PyList_SetItem(dec_rate_list, i, PyFloat_FromDouble(dec_rate[i]));
    }
    return Py_BuildValue("{s:i, s:O, s:O, s:O, s:O}", "status", stat,
			 "ra", ra_list, "dec", dec_list,
			 "ra_rate", ra_rate_list, "dec_rate", dec_rate_list);
}

extern PyObject *object_doppler ( char *object, long mjd, double utc,
				  int number_of_velocities, double interval )
{
    static int np = 0;
    // return value
    static double *frac;
    if (np != number_of_velocities) {
        delete frac;
	np = number_of_velocities;
	frac = new double[np];
    }

    PyObject *frac_list = PyList_New(np);

    int stat = ss_object_doppler (object, mjd ,utc, number_of_velocities,
				  interval, frac);

    int i;
    for (i = 0; i < np; i++) {
        PyList_SetItem(frac_list, i, PyFloat_FromDouble(frac[i]));
    }
    return Py_BuildValue("{s:i, s:O}", "status", stat, "frac", frac_list);
}

extern PyObject *doppler_fraction ( double ra2000, double dec2000,
				 long mjd, double utc, int
				 number_of_values, double interval )
{
    static int np = 0;
    // return value
    static double *frac;
    if (np != number_of_values) {
        delete frac;
	np = number_of_values;
	frac = new double[np];
    }

    PyObject *frac_list = PyList_New(np);

    int stat = ss_doppler_fraction (ra2000, dec2000, mjd, utc,
				    number_of_values, interval, frac);

    int i;
    for (i = 0; i < np; i++) {
        PyList_SetItem(frac_list, i, PyFloat_FromDouble(frac[i]));
    }
    return Py_BuildValue("{s:i, s:O}", "status", stat, "frac", frac_list);
}

extern PyObject *observer_position_velocity ( long mjd, double utc,
					      int number_of_values,
					      double interval )
{
    static int np = 0;
    // return values
    static double *x ;
    static double *y ;
    static double *z ;
    static double *dx ;
    static double *dy ;
    static double *dz ;
    if (np != number_of_values) {
        if (np != 0) {
	  delete x;
	  delete y;
	  delete z;
	  delete dx;
	  delete dy;
	  delete dz;
	}
	np = number_of_values;
	x = new double[np];
	y = new double[np];
	z = new double[np];
	dx = new double[np];
	dy = new double[np];
	dz = new double[np];
    }

    PyObject *x_list = PyList_New(np);
    PyObject *y_list = PyList_New(np);
    PyObject *z_list = PyList_New(np);
    PyObject *dx_list = PyList_New(np);
    PyObject *dy_list = PyList_New(np);
    PyObject *dz_list = PyList_New(np);

    int stat = ss_observer_position_velocity (mjd, utc, np, interval,
					      x, y, z, dx, dy, dz);

    int i;
    for (i = 0; i < np; i++) {
        PyList_SetItem(x_list, i, PyFloat_FromDouble(x[i]));
        PyList_SetItem(y_list, i, PyFloat_FromDouble(y[i]));
        PyList_SetItem(z_list, i, PyFloat_FromDouble(z[i]));
        PyList_SetItem(dx_list, i, PyFloat_FromDouble(dx[i]));
        PyList_SetItem(dy_list, i, PyFloat_FromDouble(dy[i]));
        PyList_SetItem(dz_list, i, PyFloat_FromDouble(dz[i]));
    }
    return Py_BuildValue("{s:i, s:O, s:O, s:O, s:O, s:O, s:O}",
			 "status", stat,
			 "x", x_list, "y", y_list, "z", z_list,
			 "dx", dx_list, "dy", dy_list, "dz", dz_list);
}

extern PyObject *pulse_delay ( double ra2000, double dec2000, long mjd,
			    double utc, int number_of_delays,
			    double interval )
{
    static int np = 0;
    // return value
    static double *delay;
    if (np != number_of_delays) {
        delete delay;
	np = number_of_delays;
	delay = new double[np];
    }

    PyObject *delay_list = PyList_New(np);

    int stat = ss_pulse_delay (ra2000, dec2000, mjd, utc, np, interval, delay);

    int i;
    for (i = 0; i < np; i++) {
        PyList_SetItem(delay_list, i, PyFloat_FromDouble(delay[i]));
    }
    return Py_BuildValue("{s:i, s:O}", "status", stat, "delay", delay_list);
}

extern PyObject *utc_to_tdb ( long mjd, double utc )
{
    // return values
    long  tdb_mjd;
    double tdb;

    ss_utc_to_tdb (mjd, utc, tdb_mjd, tdb);

    return Py_BuildValue("{s:l, s:d}", "tdb_mjd", tdb_mjd, "tdb", tdb);
}

extern double utc_to_last ( long mjd, double utc )
{
    return ss_utc_to_last (mjd, utc);
}

extern PyObject *last_to_utc ( long mjd, double last )
{
    // return values
    double *utc = ss_last_to_utc (mjd, last);

    return Py_BuildValue("{s:d, s:d}", "utc", utc[0], "utc_alt", utc[1]);
}

extern PyObject *epoch_to_j2000 ( long mjd, double utc, double ra,
				  double dec )
{
    PositionT pos = ss_epoch_to_j2000 (mjd, utc, ra, dec);

    return Py_BuildValue("{s:d, s:d, s:d}", "epoch", pos.epoch,
			 "ra", pos.ra, "dec", pos.dec);
}

extern PyObject *j2000_to_epoch ( long mjd, double utc,
				  double ra2000, double dec2000 )
{
    PositionT pos = ss_j2000_to_epoch (mjd, utc, ra2000, dec2000);

    return Py_BuildValue("{s:d, s:d, s:d}", "epoch", pos.epoch,
			 "ra", pos.ra, "dec", pos.dec);
}

extern PyObject *add_aberration ( long mjd, double utc, double ra, double dec )
{
    PositionT pos = ss_add_aberration (mjd, utc, ra, dec);

    return Py_BuildValue("{s:d, s:d}", "ra", pos.ra, "dec", pos.dec);
}

extern PyObject *remove_aberration ( long mjd, double utc, double ra,
				     double dec )
{
    PositionT pos = ss_remove_aberration (mjd, utc, ra, dec);

    return Py_BuildValue("{s:d, s:d}", "ra", pos.ra, "dec", pos.dec);
}

extern int set_observer_coordinates ( double x, double y, double z )
{
    ss_set_observer_coordinates (x, y, z);
    return 0;
}

extern PyObject *geocentric_observer_track ( long mjd, double utc,
					     int number_of_values,
					     double interval )
{
    static int np = 0;
    // return values
    static double *obs_x;
    static double *obs_y;
    static double *obs_z;
    static double *obs_x_rate;
    static double *obs_y_rate;
    static double *obs_z_rate;
    if (np != number_of_values) {
        delete obs_x;
        delete obs_y;
        delete obs_z;
        delete obs_x_rate;
        delete obs_y_rate;
        delete obs_z_rate;
	np = number_of_values;
	obs_x = new double[np];
	obs_y = new double[np];
	obs_z = new double[np];
	obs_x_rate = new double[np];
	obs_y_rate = new double[np];
	obs_z_rate = new double[np];
    }

    PyObject *x_list = PyList_New(np);
    PyObject *y_list = PyList_New(np);
    PyObject *z_list = PyList_New(np);
    PyObject *x_rate_list = PyList_New(np);
    PyObject *y_rate_list = PyList_New(np);
    PyObject *z_rate_list = PyList_New(np);

    int stat = ss_geocentric_observer_track (mjd, utc, np, interval,
					     obs_x, obs_y, obs_z,
					     obs_x_rate, obs_y_rate,
					     obs_z_rate);
    int i;
    for (i = 0; i < np; i++) {
        PyList_SetItem(x_list, i, PyFloat_FromDouble(obs_x[i]));
        PyList_SetItem(y_list, i, PyFloat_FromDouble(obs_y[i]));
        PyList_SetItem(z_list, i, PyFloat_FromDouble(obs_z[i]));
        PyList_SetItem(x_rate_list, i, PyFloat_FromDouble(obs_x_rate[i]));
        PyList_SetItem(y_rate_list, i, PyFloat_FromDouble(obs_y_rate[i]));
        PyList_SetItem(z_rate_list, i, PyFloat_FromDouble(obs_z_rate[i]));
    }
    return Py_BuildValue("{s:i, s:O, s:O, s:O, s:O, s:O, s:O}",
			 "status", stat,
			 "x", x_list, "y", y_list, "z", z_list,
			 "x_rate", x_rate_list, "y_rate", y_rate_list,
			 "z_rate", z_rate_list);
}


extern PyObject *barycentric_observer_track ( long mjd, double utc,
					      int number_of_values,
					      double interval )
{
    static int np = 0;
    // return values
    static double *obs_x;
    static double *obs_y;
    static double *obs_z;
    static double *obs_x_rate;
    static double *obs_y_rate;
    static double *obs_z_rate;
    if (np != number_of_values) {
        if (np != 0) {
	  delete obs_x;
	  delete obs_y;
	  delete obs_z;
	  delete obs_x_rate;
	  delete obs_y_rate;
	  delete obs_z_rate;
	}
	np = number_of_values;
	obs_x = new double[np];
	obs_y = new double[np];
	obs_z = new double[np];
	obs_x_rate = new double[np];
	obs_y_rate = new double[np];
	obs_z_rate = new double[np];
    }

    PyObject *x_list = PyList_New(np);
    PyObject *y_list = PyList_New(np);
    PyObject *z_list = PyList_New(np);
    PyObject *x_rate_list = PyList_New(np);
    PyObject *y_rate_list = PyList_New(np);
    PyObject *z_rate_list = PyList_New(np);

    int stat = ss_barycentric_observer_track (mjd, utc, np, interval,
					      obs_x, obs_y, obs_z,
					      obs_x_rate, obs_y_rate,
					      obs_z_rate);
    int i;
    for (i = 0; i < np; i++) {
        PyList_SetItem(x_list, i, PyFloat_FromDouble(obs_x[i]));
        PyList_SetItem(y_list, i, PyFloat_FromDouble(obs_y[i]));
        PyList_SetItem(z_list, i, PyFloat_FromDouble(obs_z[i]));
        PyList_SetItem(x_rate_list, i, PyFloat_FromDouble(obs_x_rate[i]));
        PyList_SetItem(y_rate_list, i, PyFloat_FromDouble(obs_y_rate[i]));
        PyList_SetItem(z_rate_list, i, PyFloat_FromDouble(obs_z_rate[i]));
    }
    return Py_BuildValue("{s:i, s:O, s:O, s:O, s:O, s:O, s:O}",
			 "status", stat,
			 "x", x_list, "y", y_list, "z", z_list,
			 "x_rate", x_rate_list, "y_rate", y_rate_list,
			 "z_rate", z_rate_list);
}

extern PyObject *barycentric_earth_track ( long mjd, double utc,
					   int number_of_values,
					   double interval )
{
    static int np = 0;
    // return values
    static double *earth_x;
    static double *earth_y;
    static double *earth_z;
    static double *earth_x_rate;
    static double *earth_y_rate;
    static double *earth_z_rate;
    if (np != number_of_values) {
        delete earth_x;
        delete earth_y;
        delete earth_z;
        delete earth_x_rate;
        delete earth_y_rate;
        delete earth_z_rate;
	np = number_of_values;
	earth_x = new double[np];
	earth_y = new double[np];
	earth_z = new double[np];
	earth_x_rate = new double[np];
	earth_y_rate = new double[np];
	earth_z_rate = new double[np];
    }

    PyObject *x_list = PyList_New(np);
    PyObject *y_list = PyList_New(np);
    PyObject *z_list = PyList_New(np);
    PyObject *x_rate_list = PyList_New(np);
    PyObject *y_rate_list = PyList_New(np);
    PyObject *z_rate_list = PyList_New(np);

    int stat = ss_barycentric_earth_track (mjd, utc, np, interval,
					   earth_x, earth_y, earth_z,
					   earth_x_rate, earth_y_rate,
					   earth_z_rate);
    int i;
    for (i = 0; i < np; i++) {
        PyList_SetItem(x_list, i, PyFloat_FromDouble(earth_x[i]));
        PyList_SetItem(y_list, i, PyFloat_FromDouble(earth_y[i]));
        PyList_SetItem(z_list, i, PyFloat_FromDouble(earth_z[i]));
        PyList_SetItem(x_rate_list, i, PyFloat_FromDouble(earth_x_rate[i]));
        PyList_SetItem(y_rate_list, i, PyFloat_FromDouble(earth_y_rate[i]));
        PyList_SetItem(z_rate_list, i, PyFloat_FromDouble(earth_z_rate[i]));
    }
    return Py_BuildValue("{s:i, s:O, s:O, s:O, s:O, s:O, s:O}",
			 "status", stat,
			 "x", x_list, "y", y_list, "z", z_list,
			 "x_rate", x_rate_list, "y_rate", y_rate_list,
			 "z_rate", z_rate_list);
}

extern PyObject *barycentric_object_track ( char *object, long mjd, double utc,
					    int number_of_values,
					    double interval)
{
    static int np = 0;
    // return values
    static double *x;
    static double *y;
    static double *z;
    static double *x_rate;
    static double *y_rate;
    static double *z_rate;
    if (np != number_of_values) {
        delete x;
        delete y;
        delete z;
        delete x_rate;
        delete y_rate;
        delete z_rate;
	np = number_of_values;
	x = new double[np];
	y = new double[np];
	z = new double[np];
	x_rate = new double[np];
	y_rate = new double[np];
	z_rate = new double[np];
    }

    PyObject *x_list = PyList_New(np);
    PyObject *y_list = PyList_New(np);
    PyObject *z_list = PyList_New(np);
    PyObject *x_rate_list = PyList_New(np);
    PyObject *y_rate_list = PyList_New(np);
    PyObject *z_rate_list = PyList_New(np);

    int stat = ss_barycentric_object_track (object, mjd, utc, np, interval,
					    x, y, z,
					    x_rate, y_rate, z_rate);
    int i;
    for (i = 0; i < np; i++) {
        PyList_SetItem(x_list, i, PyFloat_FromDouble(x[i]));
        PyList_SetItem(y_list, i, PyFloat_FromDouble(y[i]));
        PyList_SetItem(z_list, i, PyFloat_FromDouble(z[i]));
        PyList_SetItem(x_rate_list, i, PyFloat_FromDouble(x_rate[i]));
        PyList_SetItem(y_rate_list, i, PyFloat_FromDouble(y_rate[i]));
        PyList_SetItem(z_rate_list, i, PyFloat_FromDouble(z_rate[i]));
    }
    return Py_BuildValue("{s:i, s:O, s:O, s:O, s:O, s:O, s:O}",
			 "status", stat,
			 "x", x_list, "y", y_list, "z", z_list,
			 "x_rate", x_rate_list, "y_rate", y_rate_list,
			 "z_rate", z_rate_list);
}
