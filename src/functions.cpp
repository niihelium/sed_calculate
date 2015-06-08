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
        if (fabs(x) < 1e-5l)
                return x + 0.5l*x*x;
        else
                return exp(x) - 1.0l;
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
	long double F_c =2.0l + (20.0l*atan(tau))/(PIl3);
	return F_c;
}

long double F_irr(long double star_luminocity, long double r) {
	long double F_irr = (star_luminocity / (PIl4 *r*r)) * 0.5l; //cos(gamma_irr); //eq 8 cos(gamma_irr) = 0.05 L/star_luminocity = Arate 3+4
	return F_irr;
}

long double t_surf(Cell cell, long double star_luminocity) {
	//std::cout << star_luminocity;
	//std::cin.ignore();
	long double Fc05 = 0.5l *F_c(cell.tau_);//ok
	long double tmp4 = pow(cell.T_midplane_, 4.0l); //ok //ORLY???????????????????????????????//
	long double d1tau = (1.0l / (1.0l + cell.tau_)); //ok
	
	long double v1 = Fc05 * tmp4 * d1tau;
	/*cout << "v1=" << v1 << endl;
	cin.ignore();*/
	long double v2 = pow(t_background, 4.0l); // 15 formula
	long double v3 = F_irr(star_luminocity, cell.r_) / SIGMA;
	//cout << "t_surf() v2=" << v2 << endl << endl;
	long double t_surf = v1 + v2 + v3;

	/*if(t_surf >= pow(500.0l, 4)){
		cout << " in t_surf(): t_surf=" << pow(t_surf, 0.25l) << endl;
		cout << " in t_surf(): t_surf=" << t_surf << " v1=" << v1 << " v2=" << v2 << endl;
		cin.ignore();
	}*/
	/*cout << "t_surf=" << pow(t_surf, 0.25l) << endl;
	cin.ignore();*/
	return t_surf;
}

long double planck(long double nu, long double t){
	long double B;
	long double v1 = (h2*pow(nu, 3.0l))/C_sq;
	long double v2 = h*nu/(K_BOLTZ*t);
	long double exponent_minus_one = expm1(v2);
	B = v1/exponent_minus_one;
	return B;
}



//8888888                                         8888888b.  d8b                   
//  888                                           888  "Y88b Y8P                   
//  888                                           888    888                       
//  888   88888b.  88888b.   .d88b.  888d888      888    888 888 .d8888b   .d8888b 
//  888   888 "88b 888 "88b d8P  Y8b 888P"        888    888 888 88K      d88P"    
//  888   888  888 888  888 88888888 888          888    888 888 "Y8888b. 888      
//  888   888  888 888  888 Y8b.     888          888  .d88P 888      X88 Y88b.    
//8888888 888  888 888  888  "Y8888  888          8888888P"  888  88888P'  "Y8888P 

long double omega(long double r){
	long double t = (G*M_star)/pow(r, 3.0l);
	long double omega = pow(t, 0.25l);
	return omega;
}

long double nu_sigma(long double a_rate, long double r, long double star_radius){
	long double nu_sigma = (a_rate/PIl3)*(1.0l - pow((star_radius/r), 0.25l));
	//cout << "nu_sigma=" << nu_sigma << endl;
	//cin.ignore();
	return nu_sigma;
}

long double t_surf_in(Cell cell, long double star_luminocity, long double star_radius, long double a_rate) {
	long double v1 = (9.0l/SIGMA*8.0l)*pow(omega(cell.r_), 2.0l)*nu_sigma(a_rate, cell.r_, star_radius);
	long double t_surf_in = v1 + pow(t_background, 4.0l) + F_irr(star_luminocity, cell.r_)/SIGMA;
	return t_surf_in;
}

 //.d8888b.                    888                    888       .d8888b.  888                     
//d88P  Y88b                   888                    888      d88P  Y88b 888                     
//888    888                   888                    888      Y88b.      888                     
//888         .d88b.  88888b.  888888 888d888 8888b.  888       "Y888b.   888888  8888b.  888d888 
//888        d8P  Y8b 888 "88b 888    888P"      "88b 888          "Y88b. 888        "88b 888P"   
//888    888 88888888 888  888 888    888    .d888888 888            "888 888    .d888888 888     
//Y88b  d88P Y8b.     888  888 Y88b.  888    888  888 888      Y88b  d88P Y88b.  888  888 888     
// "Y8888P"   "Y8888  888  888  "Y888 888    "Y888888 888       "Y8888P"   "Y888 "Y888888 888     
//FIXME ok                       
long double t_eff(long double star_luminocity, long double star_radius){
	long double t_sub = star_luminocity/(PIl4*(star_radius*star_radius)*SIGMA);
	long double t = pow(t_sub, 0.25l);
	return t;
}

//MISC

long double unlog10(double long value)
{
	return pow(10.0l, value);
}

long double interpolate_log(long double first_x, long double first_y, long double second_x, long double second_y,  long double x)
{
	first_y = log10(first_y);
	second_y = log10(second_y);
	first_x = log10(first_x);
	second_x = log10(second_x);
	x = log10(x);
	long double y = ((x-first_x)*(second_y-first_y))/(second_x-first_x)+first_y;
	//cout << "In interpolate_my: y=" << y << endl;
	return unlog10(y);
}

bool double_equals(long double a, long double b, long double epsilon)
{
    return fabs(a - b) < epsilon;
} 