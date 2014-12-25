// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cmath>
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
#include "./data_reader.h"
#include "./functions.h"

using namespace std;

static const double gamma_incl = 0;
static const double d = 3.08568000e+20; //100pc distance to object cm
static const double phi0 = 0;

int lines_count;

//Calculation variables
//varXXX.dat
double system_age;
double star_luminocity;

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
			//cout << "Temp is: " << *it << endl;
		}
	}
	//cout << "Index is: " << index << endl;

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

double calcOuterDiscCell(Cell cell, double nu) {
	double B = planck(nu, t_surf(cell, star_luminocity));
	double F;
	F = (cell.s_ / pow(d, 2)) * B* (1 - pow(M_E,-cell.sigma_ * getKnu(t_surf(cell, star_luminocity)) * (1 / cos(gamma_incl)))); // 14 formula (planck+ t_surf)
	return F;
}

int main(int argc, char **argv) {
	/*	int number = countLines();
	 cout << "Lines in file: " << number << endl;*/
	cout << "Reading data" << endl;
	DataReader::readOuterData(system_age, cells_outer);
	DataReader::readKnu(k_nu_temp, k_nu);
	DataReader::readArate(system_age, star_luminocity);
	cout << "Data reading complete" << endl;
	//cout << system_age << endl;
	/*for (list<Cell>::iterator it = cells_outer.begin(); it != cells_outer.end();
			it++) {
		cout << *it;
		cout << endl;
	}*/
	double spectre[100];
	int spectre_pos;
	cout << "First cell: "  << cells_outer.front().line_index_ << " Last cell: " << cells_outer.back().line_index_ << endl;
	for (Cell & cell : cells_outer) {
		cout << cell.line_index_ << endl;
		spectre_pos = 0;
		for (double nu = 1; nu < pow(10, 3); nu = nu + 10){
			double counting_result =calcOuterDiscCell(cell, nu);
			//cout << counting_result << endl;
			if (counting_result == counting_result && !std::isinf(counting_result))
				spectre[spectre_pos] = spectre[spectre_pos] + counting_result;
			spectre_pos++;
		}
		
	}

	ofstream fout("result.dat");
	if(fout.is_open()){
    //file opened successfully so we are here
    cout << "File Opened successfully!!!. Writing data from array to file" << endl;

		for(int i = 0; i < 100; i++)
		{
      		fout << spectre[i] << endl; //writing ith character of array in the file
		}
    cout << "Array data successfully saved into the file test.txt" << endl;
	}
	cout << "done" << endl;
	cin.ignore();
	return 0;
}
