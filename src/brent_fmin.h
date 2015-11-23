// 1-d optimization by Brent's method
#ifndef BRENT_FMIN_H
#define BRENT_FMIN_H

double qtl2_Brent_fmin(double ax, double bx, double (*f)(double, void *),
                       void *info, double tol);

#endif // BRENT_FMIN_H
