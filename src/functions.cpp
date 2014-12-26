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
//classes
#include "functions.h"
#include "cell.h"
#include "constants.cpp"

double F_c(double theta) {
	double F_c;
	F_c =2 + 20*pow(tan(theta), -1)/(3*M_PIl);
	return F_c;
}

double F_irr(double star_luminocity, double r) {
	double F_irr;
	F_irr = (star_luminocity / 4 * M_PIl * pow(r, 2)) * 0.5; //cos(gamma_irr); //eq 8 cos(gamma_irr) = 0.05 L/star_luminocity = Arate 3+4
	return F_irr;
}

double t_surf(Cell cell, double star_luminocity) {
	double t_surf;
	//t_surf = cell.T_midplane_+t_background;
	t_surf = 1 / 2 * (F_c(cell.theta_) * cell.T_midplane_) * (1 / 1 + cell.theta_) + t_background // 15 formula
			+ F_irr(star_luminocity, cell.r_) / SIGMA;
	return t_surf;
}


double planck(double nu, double t_surf) {
	double B;
	B = (2 * h * pow(nu, 3) / pow(C, 2))
			* (1 / pow(M_E, h * nu / K_BOLTZ * t_surf) - 1); //Planck law
	return B;
}

