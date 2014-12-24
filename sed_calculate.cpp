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
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
//classes
#include "cell.cpp"
#include "constants.h"
#include "functions.cpp"

using namespace std;

static const char *grid = "../data/grid.dat";
static const char *var030 = "../data/var030.dat";
static const char *opacities = "../data/opacities_oh2.tbl";

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

list<double> k_nu_temp;
list<double> k_nu;

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

		cout << "Current line:" << curr_line_index << endl;

		file_grid >> rt >> phit >> St;
		file_var >> sigmat >> t_midplane >> thetat;

		//cout << "rt:"<< rt << "phit:" << phit << "St:" << St << endl;
		//cout << "rt:"<< sigmat << "t_midplane:" << t_midplane<< "thetat:" << thetat << endl;

		if (St != 0 && rt >= 0.970332E-04) { //Reading line if S != 0 and r >= 20 AU
			cells_outer.insert(cells_outer.end(),
					Cell(curr_line_index, rt, phit, St, sigmat, t_midplane,
							thetat, t_surf(t_background, t_midplane)));
		}
		if (file_grid.eof() || file_var.eof())
			break;
	}
	file_grid.close();
	file_var.close();
}

void readKnu() {
	double k_nu_temp_t, k_nu_t;

	ifstream file_opacities(opacities);
	while (true) {
		file_opacities >> k_nu_temp_t >> k_nu_t;

		k_nu_temp.insert(k_nu_temp.end(), k_nu_temp_t);
		k_nu.insert(k_nu.end(), k_nu_t);

		if (file_opacities.eof())
			break;
	}

}

double getKnu(double t) {
	int index;
	double Knu, sum, min;
	min = t;
	for (list<double>::iterator it = k_nu_temp.begin(); it != k_nu_temp.end();
			it++) {
		sum = abs(*it - t);
		if (sum < min) {
			min = sum;
			index = it - k_nu_temp.begin();
		}
	}

	if (k_nu.size() > index) {
		list<double>::iterator it = k_nu.begin();
		advance(it, index);
		// 'it' points to the element at index 'N'
		return *it;
	}

	return 0;
}

double calcOuterDiscCell(Cell cell, int nu) {
	double B = planck(nu, cell.t_surface_);
	double F;
	F = (cell.s_ / d ^ 2) * B
			* (1
					- pow(M_E,
							-epsilon * getKnu(cell.t_surface_)
									* (1 / cos(gamma_incl)))); // 14 formula (planck+ t_surf)
	return F;
}

int main(int argc, char **argv) {
	/*	int number = countLines();
	 cout << "Lines in file: " << number << endl;*/
	cout << "Reading data" << endl;
	readOuterData();
	cout << "Data reading complete" << endl;
	cout << system_age << endl;
	for (list<Cell>::iterator it = cells_outer.begin(); it != cells_outer.end();
			it++) {
		cout << *it;
		cout << endl;
	}

	for (Cell & cell : cells_outer) {
		for (int nu = pow(10, -1); nu < pow(10, 3); nu = nu + 0.1)
			calcOuterDiscCell(cell, nu);
	}

	cin.ignore();
	return 0;
}
