/*
 * functions.cpp
 *
 *  Created on: Dec 25, 2014
 *      Author: nii
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//boost
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
//classes
#include "cell.h"
#include "constants.cpp"

double F_c() {
	double F_c;
	F_c = Lambda / SIGMA * t_midplane * (theta / 1 + pow(theta, 2)); //5 formula
	return F_c;
}

double F_irr() {
	double F_irr;
	F_irr = (L / 4 * M_PI * pow(r, 2)) * cos(gamma_irr); //eq 8 cos(gamma_irr) = 0.05 L = Arate 3+4
	return F_irr;
}

double t_surf(double t_background, double t_midplane) {
	double t_surf;
	t_surf = 1 / 2 * (F_c() * t_midplane) * (1 / 1 + t) + t_background // 15 formula
			+ F_irr() / SIGMA;
	return t_surf;
}

/*
double planck(double nu, double t_surf) {
	double B;
	B = (2 * h * pow(nu, 3) / pow(c, 2))
			* (1 / pow(M_E, h * nu / K_BOLtZ * t_surf) - 1); //Planck law
	return B;
}*/

