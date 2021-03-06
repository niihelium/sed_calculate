/*
 * constants.h
 *
 *  Created on: Dec 23, 2014
 *      Author: nii
 */
 #include "math.h"

//physical constants
static const long double SIGMA = 5.67040040e-05l; //Stefan Boltzmann constant, erg*cm-2*s-1*K-4
static const long double K_BOLTZ = 1.38064881e-16l; //Boltzmann constant, erg*K
static const long double h = 6.62606957e-27l; //Plank constant, erg*s
static const long double t_background = 1.0e+01l; //10 k
static const long double C = 2.99792458e+10l; //light speed
static const long double C_sq = C*C; //light speed^2
static const long double sun_radius = 6.96e+10l; //sun radius in cm
static const long double mkm_to_hz = 2.99792458e+14l;// 1 mkm = mkm_to_hz Hz
static const long double G = 6.6725985e-8; //Gravitational constant cm3 g-1 s-2

//functions
static const long double PIl3 = 3.0l*M_PIl;
static const long double PIl4 = 4.0l*M_PIl;
static const long double h2 = 2.0l*h;

static const long double M_star = 0.75l * 1.99e+33l; //star mass in g



//sed_calculate
static const long double d = 3.08568000e+20l; //100pc distance to object cm
static const long double d_sq = d*d; //100pc distance to object cm
static const long double phi0 = 0.0l;





