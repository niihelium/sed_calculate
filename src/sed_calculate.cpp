// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <list>
#include <iterator>
//boost
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
//classes
#include "functions.cpp"
#include "data_reader.cpp"

using namespace std;

static const double gamma_incl = 0;
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



double getKnu(double t) {
	int index;
	double Knu, sum, min;
	min = t;
	for (list<double>::iterator it = k_nu_temp.begin(); it != k_nu_temp.end();
			it++) {
		sum = abs(*it - t);
		if (sum < min) {
			min = sum;
			index = distance(k_nu_temp.begin(), it);
			cout << "Temp is: " << *it << endl;
		}
	}
	cout << "Index is: " << index << endl;

	if (k_nu.size() > index) {
		list<double>::iterator it = k_nu.begin();
		for (int i = 0; i < index; ++i)
		{			
			it++;
		}
		Knu = *it;
	}

	return Knu;
}

/*double calcOuterDiscCell(Cell cell, int nu) {
	double B = planck(nu, cell.t_surface_);
	double F;
	F = (cell.s_ / d ^ 2) * B
			* (1
					- pow(M_E,
							-epsilon * getKnu(cell.t_surface_)
									* (1 / cos(gamma_incl)))); // 14 formula (planck+ t_surf)
	return F;
}*/

int main(int argc, char **argv) {
	/*	int number = countLines();
	 cout << "Lines in file: " << number << endl;*/
	cout << "Reading data" << endl;
	DataReader::readOuterData();
	DataReader::readKnu();
	cout << "Data reading complete" << endl;
	cout.precision(22);
	/*for (list<double>::iterator it = k_nu_temp.begin(); it != k_nu_temp.end();
				it++) {
		cout << *it << endl;
	}*/

	cout << "K_nu is:" << getKnu(4);

	/*for (list<double>::iterator it = k_nu.begin(); it != k_nu.end();
					it++) {
			cout << *it << endl;
		}*/
	//cout << system_age << endl;
	/*for (list<Cell>::iterator it = cells_outer.begin(); it != cells_outer.end();
			it++) {
		cout << *it;
		cout << endl;
	}*/

	/*for (Cell & cell : cells_outer) {
		for (int nu = pow(10, -1); nu < pow(10, 3); nu = nu + 0.1)
			calcOuterDiscCell(cell, nu);
	}*/

	cin.ignore();
	return 0;
}
