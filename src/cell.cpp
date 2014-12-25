/*
* cell.cpp
*
*  Created on: Nov 8, 2014
*      Author: nii
*/

#include <iostream>
#include "cell.h"
using namespace std;


Cell::Cell(int line_index, double r, double phi, double s, double sigma, double T_midplane, double theta, double t_surface){
     line_index_ = line_index;
     sigma_ = sigma;
     T_midplane_ = T_midplane;
     theta_ = theta;
     r_ = r;
     phi_ = phi;
     s_ = s;
     t_surface_ = t_surface; 
};


/*ostream& operator<<(ostream& os, const Cell& cell) {
os << cell.line_index_ << cell.r_ << cell.phi_ << cell.T_midplane_;
return os;
}*/

/*Cell::Cell(int line_index, double r, double phi, double s, double sigma, double T_midplane, double theta) {
line_index_ = line_index;

sigma_ = sigma;
T_midplane_ = T_midplane;
theta_ = theta;

r_ = r;
phi_ = phi;
s_ = s;
}*/
