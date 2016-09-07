extern void set_ephemeris_dir ( char *path, char *file_name );

extern PyObject *object_track ( char *object, long mjd, double utc,
				    int number_of_positions,
				    double interval );

extern PyObject *object_doppler ( char *object, long mjd, double utc,
				  int number_of_velocities, double interval );

extern PyObject *doppler_fraction ( double ra2000, double dec2000,
				 long mjd, double utc, int
				 number_of_values, double interval );

extern PyObject *observer_position_velocity ( long mjd, double utc,
					    int number_of_values,
					    double interval );

extern PyObject *pulse_delay ( double ra2000, double dec2000, long mjd,
			    double utc, int number_of_delays,
			    double interval );

extern PyObject *utc_to_tdb ( long mjd, double utc );

extern double utc_to_last ( long mjd, double utc );

extern PyObject *last_to_utc ( long mjd, double last );

extern PyObject *epoch_to_j2000 ( long mjd, double UTC, double ra,
				  double dec );

extern PyObject *j2000_to_epoch ( long mjd, double utc,
				  double ra2000, double dec2000 );

extern PyObject *add_aberration ( long mjd, double utc, double ra,
				  double dec );

extern PyObject *remove_aberration ( long mjd, double utc, double ra,
				     double dec );

extern int set_observer_coordinates ( double x, double y, double z );

extern PyObject *geocentric_observer_track ( long mjd, double utc,
					     int number_of_values,
					     double interval );

extern PyObject *barycentric_observer_track ( long mjd, double utc,
					      int number_of_values,
					      double interval );

extern PyObject *barycentric_earth_track ( long mjd, double utc,
					   int number_of_values,
					   double interval );

extern PyObject *barycentric_object_track ( char *object, long mjd, double utc,
					    int number_of_values,
					    double interval);
