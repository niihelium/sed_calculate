/*
 * functions.cpp
 *
 *  Created on: Dec 25, 2014
 *      Author: nii
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
//classes
#include "functions.h"
#include "cell.h"
#include "constants.cpp"

 using namespace std;

 long double expm1(long double x){
        if (fabs(x) < 1e-5)
                return x + 0.5*x*x;
        else
                return exp(x) - 1.0;
}


//  .d88888b.           888                              888 d8b                   
//d88P" "Y88b          888                              888 Y8P                   
//888     888          888                              888                       
//888     888 888  888 888888 .d88b.  888d888       .d88888 888 .d8888b   .d8888b 
//888     888 888  888 888   d8P  Y8b 888P"        d88" 888 888 88K      d88P"    
//888     888 888  888 888   88888888 888          888  888 888 "Y8888b. 888      
//Y88b. .d88P Y88b 888 Y88b. Y8b.     888          Y88b 888 888      X88 Y88b.    
// "Y88888P"   "Y88888  "Y888 "Y8888  888           "Y88888 888  88888P'  "Y8888P


long double F_c(long double tau) {
	long double F_c;
	F_c =(2.0 + 20.0*(atan(tau)))/(PIl3);
	return F_c;
}

long double F_irr(long double star_luminocity, long double r) {
	long double F_irr;
	F_irr = (star_luminocity / (PIl4 * pow(r, 2.0))) * 0.5; //cos(gamma_irr); //eq 8 cos(gamma_irr) = 0.05 L/star_luminocity = Arate 3+4
	return F_irr;
}

long double t_surf(Cell cell, long double star_luminocity) {
	//std::cout << star_luminocity;
	//std::cin.ignore();
	long double v1 = 0.5 * F_c(cell.tau_) * pow(cell.T_midplane_, 4) * (1.0 / 1.0 + cell.tau_);
	//cout << "t_surf() v1=" << v1 << endl;
	long double v2 = pow(t_background, 4) // 15 formula
			+ F_irr(star_luminocity, cell.r_) / SIGMA;
	//cout << "t_surf() v2=" << v2 << endl << endl;
	long double t_surf = v1 + v2;
	return t_surf;
}

/*long double planck(long double nu, long double t) { // t = 1 ???
	//cout <<"in planck() - nu: " << nu << " t: " << t << endl;
	long double B;
	long double v1 = 2.0*h*C_sq/(pow(nu, 5.0));
	long double v2 = h * C / (nu*K_BOLTZ * t);
	//cout <<"in planck() - v1: " << v1 << " v2: " << v2 << endl; 
	long double exponent_minus_one = expm1(v2);
	//cout <<"in planck() - exponent: " << exponent_minus_one << endl;
	B = v1 / exponent_minus_one; //Planck law
	//cout <<"in planck() - B: " << B << endl;
	//cin.ignore();
	return B;	
}*/
long double planck(long double nu, long double t){
	long double B;
	long double v1 = (h2*pow(nu, 3.0))/C_sq;
	long double v2 = h*nu/(K_BOLTZ*t);
	long double exponent_minus_one = expm1(v2);
	B = v1/exponent_minus_one;
	return B;
}

/*long double planck_nu(long double nu, long double t_surf) { // t_surf = 1 ???
	//std::cout <<"in planck() - nu: " << nu << " t_surf: " << t_surf; 
	//std::cin.ignore();
	long double B;
	B = (2 * h * pow(nu, 3) / pow(C, 2))*(1 / pow(M_E, h * nu / K_BOLTZ * t_surf) - 1); //Planck law
	return B;
}*/



//8888888                                         8888888b.  d8b                   
//  888                                           888  "Y88b Y8P                   
//  888                                           888    888                       
//  888   88888b.  88888b.   .d88b.  888d888      888    888 888 .d8888b   .d8888b 
//  888   888 "88b 888 "88b d8P  Y8b 888P"        888    888 888 88K      d88P"    
//  888   888  888 888  888 88888888 888          888    888 888 "Y8888b. 888      
//  888   888  888 888  888 Y8b.     888          888  .d88P 888      X88 Y88b.    
//8888888 888  888 888  888  "Y8888  888          8888888P"  888  88888P'  "Y8888P 








 //.d8888b.                    888                    888       .d8888b.  888                     
//d88P  Y88b                   888                    888      d88P  Y88b 888                     
//888    888                   888                    888      Y88b.      888                     
//888         .d88b.  88888b.  888888 888d888 8888b.  888       "Y888b.   888888  8888b.  888d888 
//888        d8P  Y8b 888 "88b 888    888P"      "88b 888          "Y88b. 888        "88b 888P"   
//888    888 88888888 888  888 888    888    .d888888 888            "888 888    .d888888 888     
//Y88b  d88P Y8b.     888  888 Y88b.  888    888  888 888      Y88b  d88P Y88b.  888  888 888     
// "Y8888P"   "Y8888  888  888  "Y888 888    "Y888888 888       "Y8888P"   "Y888 "Y888888 888     
                                                                                                
long double t_eff(long double star_luminocity, long double star_radius){
	std::cout << "t_eff(): star_luminocity=" << star_luminocity << " star_radius=" << star_radius << std::endl;
	long double t_sub = star_luminocity/(PIl4*pow(star_radius, 2.0)*SIGMA);
	std::cout << "t_sub " << t_sub << std::endl;
	long double t = pow(t_sub, 0.25);
	return t;
}