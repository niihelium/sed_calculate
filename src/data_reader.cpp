#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <list>

#include "data_reader.h"
#include "cell.h"
#include "functions.cpp"

using namespace std;

const string DataReader::grid = "../data/grid.dat";
const string DataReader::var030 = "../data/var030.dat";
const string DataReader::opacities = "../data/opacities_oh2.tbl";

void DataReader::readOuterData(double &system_age, list<Cell> &cells_outer){
			//grid.dat
			double rt, phit, St;
			//varXXX.dat
			double sigmat, t_midplane, thetat;
			int curr_line_index = 0;

			ifstream file_grid(grid.c_str());
			ifstream file_var(var030.c_str());

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

