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
#include <time.h>
//classes
#include "./data_reader.h"
#include "./functions.h"
#include "./constants.cpp"

#define MAX_DATE 12

using namespace std;

int lines_count;

//Calculation variables
//varXXX.dat
long double system_age;
long double star_luminocity;
long double star_radius;

list<Cell> cells_outer;
list<int> line_index;

list<long double> k_lambda_temp;
list<long double> k_lambda;

long double lambda_min = 1e-1;
long double lambda_max = 1e+3;
long double precision = 100.0;
long double log_lambda_min = log10(lambda_min);
long double log_lambda_max = log10(lambda_max);
long double delta_log_lambda = (log_lambda_max - log_lambda_min)/precision;


long double getKlambda(long double lambda) {
	//cout << "getKlambda is: " << endl;
	int index;
	long double Klambda, sum, min;
	min = lambda;
	for (list<long double>::iterator it = k_lambda_temp.begin(); it != k_lambda_temp.end();
			it++) {
		sum = abs(*it - lambda);
		if (sum < min) {
			min = sum;
			index = distance(k_lambda_temp.begin(), it);
			//cout << "Temp is: " << *it << endl;
		}
	}
	//cout << "Index is: " << index << endl;

	if (k_lambda.size() > index) {
		list<long double>::iterator it = k_lambda.begin();
		for (int i = 0; i < index; ++i)
		{			
			it++;
		}
		Klambda = *it;
	}
	cout << "in getKlambda Klambda: " << Klambda << " lambda: " << lambda << endl;
	cin.ignore();
	return Klambda;
}

long double calcOuterDiscCell(Cell cell, long double nu) {/*This function must be calculatet by micrones, not frequency*/
	long double surf_t = pow(t_surf(cell, star_luminocity), 0.25); //veryvery high tsurf
	//cout << "in calcOuterDiscCell() surf_t=" << surf_t << endl;
	long double B = planck(nu, surf_t);
	//cout << "in calcOuterDiscCell() B=" << B << endl;
	long double lambda = C/nu;
	long double F = (cell.s_ / d_sq) * B* (1 - pow(M_E,-cell.sigma_ * getKlambda(lambda) * (1 / cos(gamma_incl)))); // 14 formula (planck+ t_surf)
	return F;
}

long double calcStar(long double nu, long double star_radius){
	long double t = t_eff(star_luminocity, star_radius);
	//cout << "Teff in star: " << t << endl;
	cout<< "calcStar - t: " << t << endl;
	long double B = planck(nu, t);
	//cout << "Planck in star: " << B << endl;
	long double F = ((M_PIl*pow(star_radius, 2.0))/d_sq)*B;
	//cout << "F in star: " << F << endl;
	return F;
}

void readData(){
	cout << "Reading data" << endl;
	DataReader::readOuterData(system_age, cells_outer);
	DataReader::readKlambda(k_lambda_temp, k_lambda);
	DataReader::readArate(system_age, star_luminocity, star_radius);
	cout << "Data reading complete" << endl;
}

void starRoutine(){
	ofstream star_out("star_result_1.dat");
	long double lambda = lambda_min;
	long double  i = 0.0;
	while (lambda <= lambda_max){
			lambda = pow(10, log_lambda_min+i*delta_log_lambda)*0.0001;
			long double frequency = C/lambda; // cm/s /cm = 1/s

			long double counting_result = calcStar(frequency, star_radius);
			counting_result = counting_result*frequency;
			cout << "counting_result in star: "<< counting_result << endl;
			//	cin.ignore();
			if (counting_result == counting_result && !std::isinf(counting_result))
				star_out << lambda << " " << frequency << " " << counting_result << endl;		
			++i;		
	}
	star_out.close();
}

void discOutputRoutine(long double spectre[]){
	ofstream fout("result.dat");
	if(fout.is_open()){
    	cout << "File Opened successfully. Writing data from array to file" << endl;

		for(int i = 0; i < precision; i++)
		{
      		fout << spectre[i] << endl; //writing ith character of array in the file      		
		}
    	cout << "Array data successfully saved into the file result.dat" << endl;
	}
}

void discRoutine(long double spectre[] ){
	cout << "First cell: "  << cells_outer.front().line_index_ << " Last cell: " << cells_outer.back().line_index_ << endl;
	int spectre_pos = 0.0;
	for (Cell & cell : cells_outer) {
		cout << cell.line_index_ << "\r" << flush;
		//printf("\r %c", cell.line_index_);
		spectre_pos = 0.0;
		long double lambda = lambda_min;
		while (lambda < lambda_max){
			lambda = pow(10, log_lambda_min+spectre_pos*delta_log_lambda)*0.0001;
			long double frequency = C/lambda; // cm/s /cm = 1/s

			long double counting_result =calcOuterDiscCell(cell, frequency);
			//cout << "discRoutine(): counting_result=" << counting_result << endl;
			//if (counting_result == counting_result && !std::isinf(counting_result)){
				spectre[spectre_pos] = spectre[spectre_pos] + counting_result;
				//cout << "discRoutine(): spectre[spectre_pos]=" << spectre[spectre_pos] << " spectre_pos=" << spectre_pos <<endl;
			//}
			++spectre_pos;
		}		
	}
	discOutputRoutine(spectre);
}



int main(int argc, char **argv) {
	readData();

	long double spectre[(int)precision];
	for (int i = 0; i < precision; ++i)
	{
		spectre[i] = 0;
	}

	//starRoutine();
	discRoutine(spectre);
	
	
	cout << "done" << endl;

	return 0;
}
