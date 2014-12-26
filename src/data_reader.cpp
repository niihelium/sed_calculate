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

using namespace std;

const string DataReader::grid = "../data/grid.dat";
const string DataReader::var030 = "../data/var030.dat";
const string DataReader::opacities = "../data/opacities_oh2.tbl";
const string DataReader::Arate = "../data/ARate.dat";

void DataReader::readOuterData(double &system_age, list<Cell> &cells_outer){
			//grid.dat
			double rt, phit, St;
			//varXXX.dat
			double sigmat, t_midplane, thetat;
			int curr_line_index = 0;

			ifstream file_grid(grid.c_str());
			ifstream file_var(var030.c_str());

			file_var >> system_age;
			double dumb;
			cout.precision(6);
			//cout << scientific;
			//cout << scientific << system_age << endl;

			while (true) {
				++curr_line_index;				

				file_grid >> rt >> phit >> St;
				file_var >> sigmat >> t_midplane >> thetat >> dumb >> dumb >> dumb >> dumb;
				

				if (St != 0 && rt >= 0.970332E-04) { //Reading line if S != 0 and r >= 20 AU
					//cout << "Current line:" << curr_line_index << endl;
					//cout << "rt: "<< rt << " phit: " << phit << " St: " << St << endl;
					//cout << "sigmat: "<< sigmat << " t_midplane: " << t_midplane<< " thetat: " << thetat << endl;
					cells_outer.insert(cells_outer.end(),
							Cell(curr_line_index, rt, phit, St, sigmat, t_midplane,
									thetat));
				}
				if (file_grid.eof() || file_var.eof())
					break;
			}
			file_grid.close();
			file_var.close();
			cout << "Reading outer data done" << endl;
}

void DataReader::readKnu(list<double> &k_nu_temp, list<double> &k_nu){
			double k_nu_temp_t, k_nu_t;
			ifstream file_opacities(opacities.c_str());
			while (true) {
				file_opacities >> k_nu_temp_t >> k_nu_t;

				k_nu_temp.insert(k_nu_temp.end(), k_nu_temp_t);
				k_nu.insert(k_nu.end(), k_nu_t);

				if (file_opacities.eof())
					break;
			}
}

double DataReader::readArate(double system_age, double &star_luminocity){ //working good
	//cout << "readArate" << endl;
	int curr_line_index;
	double age_t, L_accr, L_photo, L_full, diff, min, dumb;
	min = std::numeric_limits<double>::max();
	//cout << "Min: " << min << endl;
	//cout << "system_age: " << system_age << endl;
	ifstream file_Arate(Arate.c_str());

	string line;
	while (getline(file_Arate, line)){
		++curr_line_index;
	    istringstream iss(line);
	    if (!(iss >> age_t >> dumb >> L_accr >> L_photo)) break; // error
	    
	    diff = abs(age_t - system_age);
		if (diff<min){
			min = diff;
	    	L_full = L_accr + L_photo;
		} else {
			break;
		}
	}
	//cout << "age_t: " << age_t << " L_accr: " << L_accr << " L_photo: " << L_photo << endl;
	//cout << curr_line_index << "diff: " << diff << endl;
	return L_full;
}
