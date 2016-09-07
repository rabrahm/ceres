// Astronomical time conversion routines header file

#ifndef ASTRTIME_H
#define ASTRTIME_H

#define MAX_NUM_LEAP_SEC     100
#define MAX_NUM_IERS_OFFSETS 2000
#ifndef PI
#define PI                   3.1415926535897932384626433
#endif
#define TWOPI                (2.0 * PI)
#define DEG2RAD              (PI / 180.0)
#define AU2METERS	     0.14959787066e+12

typedef struct {
    double utc_offset;
    long   mjd;
    char   yr[10];
    char   month[20];
    char   day[10];
} leap_secT;

typedef struct {
    double ut1_offset;           // UT1 - UTC in seconds
    double x;                    // pole offset X in arcsec
    double y;                    // pole offset X in arcsec
    long   mjd;                  // MJD of above offsets
} iersT;

typedef struct {
    double x;                    // maters in Greenwich meridian plane
    double y;                    // meters east from Greenwich meridian plane
    double z;                    // meters north from equatorial plane
} observatoryT;


int utc2tai (long mjd, double utc, long &tai_mjd, double &tai);

int tai2utc (long tai_mjd, double tai, long &mjd, double &utc);

int utc2tt (long mjd, double utc, long &tt_mjd, double &tt);

int tt2utc (long tt_mjd, double tt, long &mjd, double &utc);

int utc2tdb (long mjd, double utc, long &tdb_mjd, double &tdb);

int tdb2utc (long tdb_mjd, double tdb, long &mjd, double &utc);

void utc2ut1 (long mjd, double utc, long &ut1_mjd, double &ut1);

void utc2ut0 (long mjd, double utc, long &ut0_mjd, double &ut0);

void utc2ut2 (long mjd, double utc, long &ut2_mjd, double &ut2);

double utc2gmst (long mjd, double utc);

double utc2gast (long mjd, double utc);

double utc2lmst (long mjd, double utc);

double utc2last (long mjd, double utc);

int set_leap_sec_file_name (char *file_name);

int open_leap_sec_file();

void close_leap_sec_file();

int load_leap_sec_data();

double get_leap_sec(long mjd, double utc);

int set_iers_file_name (char *file_name);

int open_iers_file();

void close_iers_file();

int load_iers_data();

double get_ut1_offset(long mjd, double utc);

int get_pole_offsets(long mjd, double utc, double &x, double &y);

void set_latlong2xyz (double geodetic_long, double geodetic_lat,
                      double geoid_height);

void set_observatory_xyz (double x, double y, double z);

double get_observatory_longitude();

double get_observatory_geoc_latitude();

double get_observatory_radius();

double get_observatory_x();

double get_observatory_y();

double get_observatory_z();

int get_observatory_posn_vel(double pv[6], double jd, double tt,
			     int vel_flag=1);
#endif  // ASTRTIME_H
