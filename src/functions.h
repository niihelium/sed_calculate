#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "cell.h"

long double expm1(long double x);

long double F_c(long double theta);
long double F_irr(long double star_luminocity, long double r);
long double t_surf(Cell cell, long double t_midplane); 
long double planck(long double lambda, long double t);

long double t_eff(long double star_luminocity, long double star_radius);

long double unlog10(double long value);
long double interpolate_log(long double first_x, long double first_y, long double second_x, long double second_y,  long double x);
bool double_equals(long double a, long double b, long double epsilon);

#endif