#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <list>
#include <limits>
#include <cmath>

#include "./data_reader.h"
#include "./functions.h"
#include "./constants.cpp"

using namespace std;

const string DataReader::grid = "../data/grid.dat";
const string DataReader::var030 = "../data/var030.dat";
const string DataReader::opacities = "../data/opacities_oh2.tbl";
const string DataReader::Arate = "../data/ARate.dat";

void DataReader::readOuterData(long double &system_age, list<Cell> &cells_outer){
			//grid.dat
			long double rt, phit, St;
			//varXXX.dat
			long double sigmat, t_midplane, tau;
			int curr_line_index = 0;

			ifstream file_grid(grid.c_str());
			ifstream file_var(var030.c_str());

			file_var >> system_age;
			long double dumb;
			cout.precision(6);
			//cout << scientific;
			//cout << scientific << system_age << endl;

			while (true) {
				++curr_line_index;				
				file_grid >> rt >> phit >> St;
				file_var >> sigmat >> t_midplane >> tau >> dumb >> dumb >> dumb >> dumb;
				//Convert to SGS 
				rt = rt*3.08567758128e+18;
				St = St*3.08567758128e+18;
				if (St != 0 && rt >= 0.970332E-04) { //Reading line if S != 0 and r >= 20 AU
					cells_outer.insert(cells_outer.end(),
							Cell(curr_line_index, rt, phit, St, sigmat, t_midplane, tau));
				}
				if (file_grid.eof() || file_var.eof())
					break;
			}
			file_grid.close();
			file_var.close();
			cout << "Reading outer data done" << endl;
}

void DataReader::readKlambda(list<long double> &k_nu_temp, list<long double> &k_nu){
			long double k_nu_temp_t, k_nu_t;
			ifstream file_opacities(opacities.c_str());
			//First additional point for extrapolation
			k_nu_temp.insert(k_nu_temp.end(), (long double)0.1);
			k_nu.insert(k_nu.end(), (long double)5e+5);
			while (true) {
				file_opacities >> k_nu_temp_t >> k_nu_t;
				k_nu_temp.insert(k_nu_temp.end(), k_nu_temp_t);
				k_nu.insert(k_nu.end(), k_nu_t);
				if (file_opacities.eof())
					break;
			}
}

long double DataReader::readArate(long double system_age, long double &star_luminocity, long double &star_radius){ //working good
	cout << "system_age " << system_age << endl;
	int curr_line_index;
	long double age_t, L_accr, L_photo, L_full, diff, min, dumb;
	min = std::numeric_limits<long double>::max();
	ifstream file_Arate(Arate.c_str());

	string line;
	while (getline(file_Arate, line)){
		++curr_line_index;
	    istringstream iss(line);
	    if (!(iss >> age_t >> dumb >> L_accr >> L_photo >> dumb >> dumb >> dumb >> dumb >> star_radius)) break; // error

	    star_radius = star_radius * sun_radius;
	    diff = abs(age_t - system_age);
		if (diff<min){
			min = diff;
	    	L_full = L_accr + L_photo;
	    	L_full = L_full*3.827e+33;
	    	star_luminocity = L_full;
		} else {
			break;
		}

	}
/*	cout << "readArate(): star_luminocity=" << L_full << endl;
	cin.ignore();*/
	return L_full;
}
