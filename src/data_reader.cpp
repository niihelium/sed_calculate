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

void DataReader::readInnerData(long double &system_age, list<Cell> &cells_inner){
	long double rt, phit, St;
	long double sigmat, t_midplane, tau;
	int curr_line_index = 0;

	ifstream file_grid(grid.c_str());
	ifstream file_var(var030.c_str());
	file_var >> system_age;
	++curr_line_index; //1
	long double dumb;
	cout.precision(6);

	while (true){
		++curr_line_index;	//2
		file_grid >> rt >> phit >> St; //All well
		file_var >> sigmat >> t_midplane >> tau >> dumb >> dumb >> dumb >> dumb;
		if (St != 0 && rt < 0.970332E-04){
			rt = rt*3.08567758128e+18l;
			St = St*9.52140614e+36l;
			cells_inner.insert(cells_inner.end(),
					Cell(curr_line_index, rt, phit, St, sigmat, t_midplane, tau));
		}else if(rt >= 0.970332E-04 || file_grid.eof() || file_var.eof()){
			break;
		}
		
	}
	file_grid.close();
	file_var.close();
	cout << "Reading inner data done" << endl;

}

void DataReader::readOuterData(long double &system_age, list<Cell> &cells_outer){
	//grid.dat
	long double rt = 0;	
	long double phit, St;
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

	while (!file_grid.eof() || !file_var.eof()) {
		++curr_line_index;				
		file_grid >> rt >> phit >> St; //All well
		file_var >> sigmat >> t_midplane >> tau >> dumb >> dumb >> dumb >> dumb; //All well
		if (St != 0 && rt >= 0.969627362e-4l && rt < 0.242406841e-2l) { //Reading line if S != 0 and r >= 20 AU
			
		//cout << curr_line_index << " rt=" << rt << " phit=" << phit << " St=" << St << endl;
		//cin.ignore();
			//Convert to SGS 
			rt = rt*3.08567758128e+18l;
			St = St*9.52140614e+36l;
			
			cells_outer.insert(cells_outer.end(),
					Cell(curr_line_index, rt, phit, St, sigmat, t_midplane, tau));
		} 
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

long double DataReader::readArate(long double system_age, long double &star_luminocity, long double &star_radius, long double &a_rate){ //working good
	cout << "system_age " << system_age << endl;
	long double age_t, L_accr, L_photo, L_full, diff, min, dumb;
	min = std::numeric_limits<long double>::max();
	ifstream file_Arate(Arate.c_str());

	string line;
	while (getline(file_Arate, line)){

	    istringstream iss(line);
	    if (!(iss >> age_t >> a_rate >> L_accr >> L_photo >> dumb >> dumb >> dumb >> dumb >> star_radius)) break; // Looks good

	    
	    //diff = abs(age_t - system_age);
		//if (diff<min){
	    if(double_equals(age_t, 0.290000e+0l, 0.00001)){
			//min = diff;
			//star_radius = star_radius * sun_radius;
			star_radius = star_radius*sun_radius;

			a_rate = a_rate*1.99e+33l;// in year
			a_rate = a_rate/3.15569e+7l; // in second
	    	L_full = L_accr + L_photo;
	    	L_full = L_full*3.827e+33;
	    	star_luminocity = L_full;
	    	break;
		} /*else {
			break;
		}*/
		
	}
/*	cout << "age_t=" << age_t << " a_rate=" << a_rate << " L_accr=" << L_accr << " L_photo=" << L_photo << " star_radius=" << star_radius << cout;
	cin.ignore();*/
	return L_full;
}
