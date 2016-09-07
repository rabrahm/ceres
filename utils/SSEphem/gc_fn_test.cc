
// This is the test program for the functions in python_client_fn.cc.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "jpl_eph.h"
#include "astrtime.h"
#include "python_client.h"

int main (int argc, char *argv[])
{
    if (argc < 2) {
        printf ("Usage: gc_fn_test ephemeris_file_name\n");
        return 1;
    }
    if (!open_ephemeris(argv[1])) {
        return 1;
    }
    long mjd = 49719;
    double utc = 0.001 / 24.0;
    long tdb_mjd;
    double tdb;
    ss_utc_to_tdb ( mjd, utc, tdb_mjd, tdb );
    printf ("tdb_mjd = %d, tdb = %f\n", tdb_mjd, tdb * 24.0);

    double last  = ss_utc_to_last ( mjd, utc );
    printf ("last = %f\n", last * 24.0);

    double *utx = ss_last_to_utc ( mjd, last );
    printf ("utc(0) from last = %f\n", utx[0] * 24.0);
    printf ("utc(1) from last = %f\n", utx[1] * 24.0);

    double rax = 15.0;
    double decx = 60.0;
    PositionT pos = ss_j2000_to_epoch (mjd, utc, rax, decx);
    printf ("RA = %f, Dec = %f, Epoch = %f\n", pos.ra, pos.dec, pos.epoch);
    rax = pos.ra;
    decx = pos.dec;
    pos = ss_epoch_to_j2000 (mjd, utc, rax, decx);
    printf ("RA = %f, Dec = %f, Epoch = %f\n", pos.ra, pos.dec, pos.epoch);

/*
    long tt_mjd = 49719;
    double tt = 1.0 / 24.0;
    tt2utc (tt_mjd, tt, mjd, utc);
    double *df = ss_doppler_fraction ( rax, decx, mjd, utc, 1, 20.0 );
    printf ("Doppler velocity = %6.3f\n", df[0] * 299792.0);

    double *pd = ss_pulse_delay ( rax, decx, mjd, utc, 1, 20.0 );
    printf ("Pulse delay = %f sec.\n", pd[0]);

    double *rap, *decp, *rrp, *drp;
//    set_observer_coordinates (.8828800208, -4.9244824385, 3.9441306438);
    int np = ss_object_direction ( "Pluto", mjd, utc, 1, 20.0,
				      rap, decp, rrp, drp);
    printf ("Num pos = %d\n", np);
    printf ("RA = %f, Dec = %f, RA rate = %f, Dec rate = %f\n",
	    rap[0], decp[0], rrp[0], drp[0]);
    
    PositionT pos = ss_j2000_to_epoch (mjd, utc, rap[0], decp[0]);
    printf ("RA = %f, Dec = %f, Epoch = %f\n", pos.ra, pos.dec, pos.epoch);

    PositionT pa = ss_add_aberration (mjd, utc, pos.ra, pos.dec );
    printf ("RA = %f, Dec = %f\n", pa.ra, pa.dec);

    printf ("LAST = %f\n", 24.0 * utc_to_last (mjd, utc));

    int i;
    for (i = 0; i < 12; i++) {
	PositionT p = ss_add_aberration(50090, 0.20, ra, dec);
	double ra_diff = 15.0 * 3600.0 * (ra - p.ra);
	if (ra_diff < -500.0) ra_diff += 360.0 * 3600.0;
	if (ra_diff > 500.0) ra_diff -= 360.0 * 3600.0;
	double dec_diff = 3600.0 * (dec - p.dec);

	printf ("RA = %f, Dec = %f\n", p.ra, p.dec);
	printf ("RA(aber) = %5.2f, Dec(aber) = %5.2f\n", ra_diff, dec_diff);
	ra += 2.0;
    }
*/
    close_ephemeris();
    return 0;
}
