// Header file for delay.cc, the calculation of delay with respect to solar
// system locations.

#ifndef DELAY_H
#define DELAY_H

double get_geometric_delay (RefFrame ref1, RefFrame ref2,
			    long mjd, double utc,
			    double ra2000, double dec2000);

#endif	// DELAY_H
