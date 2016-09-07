
void ptwGeoc ( double p, double h, double *r, double *z );

void ptwEqecl ( double dr, double dd, double date,
                double *dl, double *db );

void ptwFk45z ( double r1950, double d1950, double bepoch,
                double *r2000, double *d2000 );

void ptwFk425 ( double r1950, double d1950, double dr1950,
                double dd1950, double p1950, double v1950,
                double *r2000, double *d2000, double *dr2000,
                double *dd2000, double *p2000, double *v2000 );

void ptwObs ( int n, char *c, char *name, double *w, double *p, double *h );

void ptwCldj ( int iy, int im, int id, double *djm, int *j );

void ptwCaldj ( int iy, int im, int id, double *djm, int *j );
