// This file is part of SSEphem.
//
//  SSEphem is a collection of C/C++ source files that provide an interface
//  to the JPL Solar System Ephemeris.
//  Copyright (C) 2009 Associated Universities, Inc. Washington DC, USA
//
//  SSEphem is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  SSEphem is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with SSEphem.  If not, see <http://www.gnu.org/licenses/>.

// This program tests the geometric delay routine that uses the JPL
// solar system ephemeris and topocentric location routines.  The test
// data are times-of-arrival for PSR 1929+10 and the period and period
// derivative obtained from the 'tempo' program.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "jpl_eph.h"
#include "astrtime.h"
#include "delay.h"
#include "sofa.h"
#include "ptw_aux.h"

int	   num_toas = 0;

struct ToaHdr {
    char   name[100];
    double ra2000;
    double dec2000;
    double period;
    double p_dot;
    double p_dot_dot;
    double epoch;
    double dm;
} toa_hdr;

struct Toa {
    double freq;
    long   mjd;
    double utc;
    double weight;
} toa[100];

int get_toas (char *file_name)
{
    FILE *toa_stream = fopen(file_name, "rt");
    if (toa_stream == NULL) {
        printf ("open_toa_file: Cannot open file '%s'\n", file_name);
        return 0;
    }
    int flag = 0;
    char label[50], eq[20];
    double d1, d2, d3, ra1950, dec1950;
    if (fscanf (toa_stream, "%s %s %s", label, eq, toa_hdr.name) != 3) flag++;
    if (flag || strcmp(label, "pulsar_name") || strcmp(eq, "=")) {
	printf ("Error reading pulsar name\n");
	return 0;
    }
    if (fscanf (toa_stream, "%s %s %lf %lf %lf",
                             label, eq, &d1, &d2, &d3) != 5)
	flag++;
    if (flag || strcmp(label, "ra1950") || strcmp(eq, "=")) {
	printf ("Error reading right ascension\n");
	return 0;
    } else {
	ra1950 = TWOPI * ((((d3 / 60.0) + d2) / 60.0 + d1) / 24.0);
    }
    if (fscanf (toa_stream, "%s %s %lf %lf %lf",
                                    label, eq, &d1, &d2, &d3) != 5)
	flag++;
    if (flag || strcmp(label, "dec1950") || strcmp(eq, "=")) {
	printf ("Error reading declination\n");
	return 0;
    } else {
	double sign = 1.0;
	if (d1 < 0.0) sign = -1.0;
	dec1950 = TWOPI * ((((d3 / 60.0) + d2) / 60.0 + fabs(d1)) / 360.0);
	dec1950 *= sign;
    }
    if (fscanf (toa_stream, "%s %s %lf", label, eq, &toa_hdr.period) != 3)
	flag++;
    if (flag || strcmp(label, "period") || strcmp(eq, "=")) {
	printf ("Error reading period\n");
	return 0;
    }
    if (fscanf (toa_stream, "%s %s %lf", label, eq, &toa_hdr.p_dot) != 3)
	flag++;
    if (flag || strcmp(label, "p_dot") || strcmp(eq, "=")) {
	printf ("Error reading p_dot\n");
	return 0;
    }
    if (fscanf (toa_stream, "%s %s %lf", label, eq, &toa_hdr.p_dot_dot) != 3)
	flag++;
    if (flag || strcmp(label, "p_dot_dot") || strcmp(eq, "=")) {
	printf ("Error reading p_dot_dot\n");
	return 0;
    }
    if (fscanf (toa_stream, "%s %s %lf", label, eq, &toa_hdr.epoch) != 3)
	flag++;
    if (flag || strcmp(label, "epoch") || strcmp(eq, "=")) {
	printf ("Error reading epoch\n");
	return 0;
    }
    if (fscanf (toa_stream, "%s %s %lf", label, eq, &toa_hdr.dm) != 3)
	flag++;
    if (flag || strcmp(label, "dispersion_measure") || strcmp(eq, "=")) {
	printf ("Error reading dispersion measure\n");
	return 0;
    }
    double b_epoch = iauEpb(2400000.5, toa_hdr.epoch);
    ptwFk45z(ra1950, dec1950, b_epoch, &toa_hdr.ra2000, &toa_hdr.dec2000);

    int ti = 0;
    int nf;
    while ((nf = fscanf (toa_stream, "%lf %ld %lf %lf", &toa[ti].freq,
                      &toa[ti].mjd, &toa[ti].utc, &toa[ti].weight)) != EOF) {
        if (nf != 4) {
            printf ("Error reading TOA entry %d\n", ti+1);
            break;
        }
        ti++;
    }
    num_toas = ti;
    fclose(toa_stream);
    return 1;
}

void print_toas()
{
    printf ("Name:      %s\n", toa_hdr.name);
    printf ("RA(2000):  %lf hours\n", 24.0 * toa_hdr.ra2000 / TWOPI);
    printf ("DEC(2000): %lf degrees\n", 360.0 * toa_hdr.dec2000 / TWOPI);
    printf ("Period:    %lf seconds\n", toa_hdr.period);
    printf ("P-dot:     %lg sec/sec\n", toa_hdr.p_dot);
    printf ("P-dot-dot: %lg sec/sec^2\n", toa_hdr.p_dot_dot);
    printf ("Epoch:     %8.2lf\n", toa_hdr.epoch);
    printf ("DM:        %8.5lf\n", toa_hdr.dm);
    printf(" Freq MHz   MJD      UTC TOA      Weight\n");
    int i;
    for (i = 0; i < num_toas; i++) {
        printf ("%9.3lf %6ld %15.12lf %7.2lf\n", toa[i].freq,
                         toa[i].mjd, toa[i].utc, toa[i].weight);
    }
}

double get_num_pulses(double epoch, double p_freq, double p_f_dot, long mjd,
                      double utc, double &pulse_phase)
{
    long tdb_mjd;
    double tdb;
    utc2tdb(mjd, utc, tdb_mjd, tdb);
    double t = (double)tdb_mjd - epoch;
    t = (t + tdb) * 86400.0;
    pulse_phase = p_freq * t + p_f_dot * t * t / 2.0;

    return p_freq + p_f_dot * t;
}

int main (int argc, char *argv[])
{
    if (argc < 2) {
        printf ("Usage: delayT ephemeris_file_name TOA_file_name\n");
        return 1;
    } else {
        printf ("\nPulsar delay test using PSR 1929+10 TOA's\n");
    }
    if (!open_ephemeris(argv[1])) {
        return 1;
    }
    if (!get_toas(argv[2])) {
	return 1;
    }
//    print_toas();

    long mjd = 50094;     // July 14, 1996
    double utc = 0.0;
    double delay = get_geometric_delay (SS_Barycentric, Topocentric,
			    mjd, utc, toa_hdr.ra2000, toa_hdr.dec2000);
    printf ("Delay = %lf sec\n", delay);

printf ("     2960597188  1 49465.38768 R 23-APR-94  0.000    0.54\n");
    double frac[50];
    int ti;
    for (ti = 0; ti < num_toas; ti++) {
        double pulse_phase;
        double p_freq = get_num_pulses(toa_hdr.epoch,
                                       4.414676683554012548, -2.25160031e-14,
                                       toa[ti].mjd, toa[ti].utc, pulse_phase);
        delay = get_geometric_delay (SS_Barycentric, Topocentric,
			             toa[ti].mjd, toa[ti].utc,
                                     toa_hdr.ra2000, toa_hdr.dec2000);
//        delay = 0.0;
        pulse_phase -= delay * p_freq;
        printf ("Ph = %lf\n", pulse_phase);
        frac[ti] = fmod(pulse_phase, 1.0);
    }
printf ("     3184684270 17 50052.88851 G  1-DEC-95   0.001    0.54\n");
    double avg = 0.0;
    for (ti = 0; ti < num_toas; ti++) {
        avg += frac[ti];
    }
    avg /= (double)num_toas;
    double ssq = 0.0;
    for (ti = 0; ti < num_toas; ti++) {
        ssq += (frac[ti] - avg) * (frac[ti] - avg);
        printf ("P_frac = %8.6lf, Resid. = %7.1lf us\n", frac[ti], (frac[ti] - avg) * 1.0e6 / 4.41467);
    }
    printf ("RMS = %5.1lf\n", sqrt(ssq / ((double)(num_toas - 1))) * 1.0e6 / 4.41467);
    close_ephemeris();
    return 0;
}
