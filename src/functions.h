#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "cell.h"

double F_c(double theta);
double F_irr(double star_luminocity, double r);
double t_surf(Cell cell, double t_midplane); 
double planck(double nu, double t_surf);

#endif