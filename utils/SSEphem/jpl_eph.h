// This header defines the data structure for the binary files, written
// and readable in C or C++, which contain the JPL solar system epemeris
// data.  The structure is similar to the DExxx or unx+XXXX.YYY binary
// file exported by JPL, but the C and C++ files do not contain any of
// the information that is unique to FORTRAN.  The JPL FORTRAN code
// equivalent of each variable is given.  The structure data packing
// details are machine and compiler dependent, so any binary files using
// these structures should be written and read by the same processor type
// and compiler executables.

// Objects or other quantities in the DE200 coefficients list.  The
// coefficients describe the motion of the objects with respect to the sun
// and the mean equator and equinox of date, except for the geocentric moon,
// whose origin is the center of the earth.  *** These must be in the order in
// which they appear in the coefficients file record, since this is the only
// place where this order definition is encoded. ***

#ifndef JPL_EPH_H
#define JPL_EPH_H

#ifndef PI
#define PI 3.1415926535897932384626433
#endif
#define RAD2DEG (360.0 / (2.0 * PI))
#define C_LIGHT 0.299792458e+06

enum Object { Mercury, Venus, EM_Bary, Mars, Jupiter, Saturn,
              Uranus, Neptune, Pluto,
              GeoCMoon,    // geocentric moon
              Sun,         // sun with respect to solar system barycenter
              Nutation,    // nutations in longitude and obliquity
              Libration,   // lunar librations
              numObj };

enum RefFrame { Heliocentric, SS_Barycentric, EM_Barycentric, Geocentric,
	        Topocentric };

typedef struct {
    int    offset;              // coefficients offset from beginning of table
                                // for this object in # of double words. These
                                // are 1 less than the JPL values to fit C
                                // index conventions.
    int    num_coeff;           // number of coefficients in the subinterval
    int    num_subinterv;       // number of subintervals in the coefficient
                                // block for this object
    int    num_components;      // number of position vector components
} obj_coeffT;

typedef struct {
    char title1[85];		// Title lines
    char title2[85];            // CHARACTER*6  TTL(14,3)
    char title3[85];
    char constant_name[400][7]; // CHARACTER*6 CNAM (400)
                                // Only NCON places filled
    double epoch1;              // DOUBLE PRECISION  SS(3)
    double epoch2;
    double interval;
    int    num_constants;       // number of active values in constant_name
                                // and constant_value  INTEGER NCON
    double au;                  // astron. unit     DOUBLE PRECISION  AU
    double emrat;               // earth/moon ratio DOUBLE PRECISION  EMRAT
    obj_coeffT obj[numObj];     // offset, # subintervals, and # of
                                // coefficients INTEGER  IPT (3,12), LPT(3)
    int    DEnumber;            // JPL DE serial number INTEGER  NUMDE
    double constant_value[400]; // DOUBLE PRECISION CVAL (400)
                                // Only NCON places filled
    int    num_coeff;           // total number of coefficients in a full
                                // interval all objects  INTEGER NCOEFF
} de_headerT;

void set_ephemeris_file_path ( char *file_path );

int open_ephemeris (char *file_nm);

void close_ephemeris();

int get_earth_posn_vel (double pv[6], RefFrame ref_frame, double time_org,
			double time_offset, int vel_flag=1);

int get_moon_posn_vel (double pv[6], RefFrame ref_frame, double time_org,
		       double time_offset, int vel_flag=1);

int get_sun_posn_vel (double pv[6], RefFrame ref_frame, double time_org,
		      double time_offset, int vel_flag=1);

int get_planet_posn_vel (double pv[6], Object planet, RefFrame ref_frame,
			 double time_org, double time_offset, int vel_flag=1);

int get_posn_vel (double pv[6], Object object, double time_org,
		  double time_offset, int vel_flag=1);

int get_nutation (double pv[4], double time_org, double time_offset,
		  int vel_flag=1);

int get_libration(double pv[6]);

double *get_coeff_subinterval (Object object, double *coeff_buffer,
			       double time_org, double time_offset,
			       double &subint_frac, int &num_coeff,
			       double &sub_interval);

double *get_coeff_buffer(double time_org, double time_offset);

double *interp (double *coeff, int num_coeff, int vel_flag, double interval,
                double time);

double get_au();

#endif  // JPL_EPH_H
