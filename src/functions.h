#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "cell.h"

long double expm1(long double x);

long double F_c(long double theta);
long double F_irr(long double star_luminocity, long double r);
long double t_surf(Cell cell, long double t_midplane); 
long double planck(long double lambda, long double t);

long double t_eff(long double star_luminocity, long double star_radius);

#endif