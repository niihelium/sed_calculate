/*
* cell.cpp
*
*  Created on: Nov 8, 2014
*      Author: nii
*/

#include <iostream>
#include "cell.h"
using namespace std;


Cell::Cell(int line_index, long double r, long double phi, long double s, long double sigma, long double T_midplane, long double tau){
     line_index_ = line_index;
     sigma_ = sigma;
     T_midplane_ = T_midplane;
     tau_ = tau;
     r_ = r;
     phi_ = phi;
     s_ = s;
};


/*ostream& operator<<(ostream& os, const Cell& cell) {
os << cell.line_index_ << cell.r_ << cell.phi_ << cell.T_midplane_;
return os;
}*/

/*Cell::Cell(int line_index, long double r, long double phi, long double s, long double sigma, long double T_midplane, long double theta) {
line_index_ = line_index;

sigma_ = sigma;
T_midplane_ = T_midplane;
theta_ = theta;

r_ = r;
phi_ = phi;
s_ = s;
}*/
