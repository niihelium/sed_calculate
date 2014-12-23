// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <list>
//boost
#include <lexical_cast.hpp>
#include <foreach.hpp>
//classes
#include "cell.cpp"
#include "constants.h"

using namespace std;

static const char *grid = "../data/grid.dat";
static const char *var030 = "../data/var030.dat";

static const double gamma_incl = 0;
static const double t_background = 1.00000000e+01; //10 k
static const double d = 3.08568000e+20; //100pc distance to object cm
static const double phi0 = 0;



int lines_count;

//Calculation variables
//varXXX.dat
double system_age;

list<Cell> cells_outer;
list<int> line_index;

void readOuterData() {
	//grid.dat
	double rt, phit, St;
	//varXXX.dat
	double sigmat, t_midplane, thetat;
	int curr_line_index = 0;

	ifstream file_grid(grid);
	ifstream file_var(var030);

	file_var >> system_age;
	cout.precision(22);
	cout << fixed << system_age << endl;

	cout << system_age;
	cin.ignore();

	while (true) {
		++curr_line_index;

		cout << "Current line:" <<curr_line_index << endl;

		file_grid >> rt >> phit >> St;
		file_var >> sigmat >> t_midplane >> thetat;

		//cout << "rt:"<< rt << "phit:" << phit << "St:" << St << endl;
		//cout << "rt:"<< sigmat << "t_midplane:" << t_midplane<< "thetat:" << thetat << endl;

		if (St != 0 && rt >= 0.970332E-04) { //Reading line if S != 0 and r >= 20 AU
			cells_outer.insert(cells_outer.end(),
					Cell(curr_line_index, rt, phit, St, sigmat, t_midplane,	thetat));
		}
		if(file_grid.eof() || file_var.eof()) break;
	}
	file_grid.close();
	file_var.close();
}

/*void readVarXXXData() {

 int curr_line_index = 0;
 list<int>::iterator it = line_index.begin();


 double sigmat, t_midplane, St;


 while () { //while not eof
 ++curr_line_index;              //increase counter current line
 if(curr_line_index == *it){
 sigma.insert(sigma.end(), sigmat);
 t_midplane.insert(t_midplane.end(), phit);
 SOuter.insert(SOuter.end(), St);
 }
 ++it;
 }

 }*/

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
 t_surf = 1 / 2 * (F_c * t_midplane) * (1 / 1 + t) + t_background // 15 formula
 + F_irr / SIGMA;
 return t_surf;
 }

 double planck(double nu, double t_surf) {
 double B;
 B = (2 * h * pow(nu, 3) / pow(c, 2))
 * (1 / pow(M_E, h * nu / K_BOLtZ * t_surf) - 1); //Planck law
 return B;
 }

 double calcOuterDiscCell(Cell cell) {
 double B = planck(nu, t_surf(t_background, cell.T_midplane_));
 double F;
 F = (cell.s_ / d ^ 2) * B
 * (1 - pow(M_E, -epsilon * Knu * (1 / cos(gamma_incl)))); // 14 formula (planck+ t_surf)
 return F;

 }

int main(int argc, char **argv) {
	/*	int number = countLines();
	 cout << "Lines in file: " << number << endl;*/
	cout << "Reading data" << endl;
	readOuterData();
	cout << "Data reading complete" << endl;
	cout << system_age << endl;
	for (list<Cell>::iterator it = cells_outer.begin(); it != cells_outer.end(); it++){
	    cout << *it;
	    cout << endl;
	}

	for(Cell & cell : cells_outer){
		calcOuterDiscCell(cell);
	}

	cin.ignore();
	return 0;
}
