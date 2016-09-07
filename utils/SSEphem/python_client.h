
typedef struct {
    double ra;
    double dec;
    double epoch;
} PositionT;

int ss_object_direction ( char *object, long mjd, double utc,
			  int number_of_positions, double interval,
			  double *ra, double *dec,
			  double *ra_rate, double *dec_rate );

int get_object_posn_vel( double pv[6], char *object, RefFrame refx,
			 double jd, double tdb );

int ss_object_doppler ( char *object, long mjd, double utc,
		        int number_of_velocities, double interval,
		        double *frac );

int ss_doppler_fraction ( double ra2000, double dec2000, long mjd, double utc,
			  int number_of_values, double interval,
			  double *frac );

int ss_observer_position_velocity ( long mjd, double utc,
				    int number_of_values, double interval,
				    double *x, double *y, double *z,
				    double *dx, double *dy, double *dz);

int ss_pulse_delay ( double ra2000, double dec2000, long mjd, double utc,
		     int number_of_delays, double interval, double *delay );

void ss_utc_to_tdb ( long mjd, double utc, long &tdb_mjd, double &tdb );

double ss_utc_to_last ( long mjd, double utc );

double *ss_last_to_utc ( long mjd, double last );

PositionT ss_epoch_to_j2000 ( long mjd, double UTC, double ra, double dec );

PositionT ss_j2000_to_epoch ( long mjd, double utc,
			      double ra2000, double dec2000 );

PositionT ss_add_aberration ( long mjd, double utc, double ra, double dec );

PositionT ss_remove_aberration ( long mjd, double utc, double ra, double dec );

void ss_aberration ( long mjd, double utc, double ra, double dec,
		     double &delta_ra, double &delta_dec );

void ss_set_observer_coordinates ( double x, double y, double z );

int ss_geocentric_observer_track(long mjd, double utc,
				  int number_of_values, double interval,
				  double *obs_x, double *obs_y, double *obs_z,
				  double *obs_x_rate, double *obs_y_rate,
				  double *obs_z_rate);

int ss_barycentric_observer_track(long mjd, double utc,
				  int number_of_values, double interval,
				  double *obs_x, double *obs_y, double *obs_z,
				  double *obs_x_rate, double *obs_y_rate,
				  double *obs_z_rate);

int ss_barycentric_earth_track(long mjd, double utc,
			       int number_of_values, double interval,
			       double *earth_x, double *earth_y,
			       double *earth_z,
			       double *earth_x_rate, double *earth_y_rate,
			       double *earth_z_rate);

int ss_barycentric_object_track(char *object, long mjd, double utc,
				int number_of_values, double interval,
				double *x, double *y, double *z,
				double *x_rate, double *y_rate,
				double *z_rate);
