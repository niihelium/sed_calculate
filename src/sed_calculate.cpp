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

//Lets calculate wavelengths before main program #################################
long double lambda_min = 1e-1l;
long double lambda_max = 1e+3l;
long double precision = 100.0l;

long double * wavelengths = new long double [(int)precision];

long double log_lambda_min = log10(lambda_min);
long double log_lambda_max = log10(lambda_max);
long double delta_log_lambda = (log_lambda_max - log_lambda_min)/precision;

void calculate_logscale_waves(){
	for (int i = 0; i < precision; ++i){	
		long double lambda = pow(10, log10(lambda_min) + delta_log_lambda*i);
		wavelengths[i] = lambda;
	}
}
//#####################################################################################

//Lets calculate opacities before main program ########################################
list<long double> k_lambda_temp;
list<long double> k_lambda;

long double * k_lambda_precalculated = new long double [(int)precision];

long double getKlambda(long double lambda){
	//cout << "getKlambda is: " << endl;
	int index;
	long double Klambda;
	bool accurate = false;
	for (list<long double>::iterator it = k_lambda_temp.begin(); it != k_lambda_temp.end();
			it++) {

		if (double_equals(lambda, *it, 0.00001l)){
			accurate = true;
			index = distance(k_lambda_temp.begin(), it);
			break;
		}
		else if(lambda < (*it)){
			accurate = false;
			index = distance(k_lambda_temp.begin(), it);
			break;
		}			
	}
	
	list<long double>::iterator it_ktemp = k_lambda_temp.begin();
	list<long double>::iterator it_k = k_lambda.begin();

	if (k_lambda.size() > index) {
		for (int i = 0; i < index; ++i){			
			it_ktemp++;
			it_k++;
		}
		if (accurate){
			Klambda = *it_k;
		}else{
			Klambda = interpolate_log(*it_ktemp, *it_k, *prev(it_ktemp, 1), *prev(it_k, 1), lambda);
		}
	}

	return Klambda;
}

void calculate_Klambda(){
	ofstream fout("opacities_extra.dat");
	for (int i = 0; i < precision; ++i){
		k_lambda_precalculated[i] = getKlambda(wavelengths[i]);
		fout << wavelengths[i] << " " << k_lambda_precalculated[i] << endl;
	}
	fout.close();
}
//#####################################################################################

void readData(){
	cout << "Reading data" << endl;
	DataReader::readOuterData(system_age, cells_outer);
	DataReader::readKlambda(k_lambda_temp, k_lambda);
	DataReader::readArate(system_age, star_luminocity, star_radius);
	cout << "Data reading complete" << endl;
}

//Star calculation adapted to new sequence ############################################
long double calcStar(long double nu, long double star_radius, long double t){
	//Calculate star effective temperature by 17 equation
	//Planck result
	long double B = planck(nu, t);
	//Function 16 result
	long double F = ((M_PIl*(star_radius*star_radius))/d_sq)*B;
	return F;
}

void starRoutine(){
	//TODO I counting all in LOGSCALE so, I shouldnt draw graph in logscale should I rewrite without logscale?
	//Assign output file for star spectre
	ofstream star_out("star_result.dat");
	//Assign starting spectre position at minal lambda
	long double lambda = lambda_min;
	long double t = t_eff(star_luminocity, star_radius);
	for (int i = 0; i < precision; ++i){
			//WTF?? Something related to logariphmic scale
			lambda = wavelengths[i]; //Lambda in micrones = OK
			//Convert length to freq
			long double frequency = C/(lambda*0.0001l); // cm/s /cm = 1/s Freq in Hz = OK

			//cin.ignore();
			//Calculation result
			long double counting_result = calcStar(frequency, star_radius, t);
			//Multiply by freq for graph value relation
			counting_result = counting_result*frequency;
			//	cin.ignore();
			// Check is result really exist not random and not inf
			if (counting_result < 1e-11l) //Looks like it's right
				counting_result = 0;
			if (counting_result == counting_result && !std::isinf(counting_result))
				star_out << lambda << " " << frequency << " " << counting_result << endl;			
	}
	star_out.close();
}
//#####################################################################################

//Looks like outer disk calculation adapted to new sequence ###########################

ofstream tfout("temps_result.dat");

long double calcOuterDiscCell(Cell cell, long double nu, int n) {

	long double surf_t = pow(t_surf(cell, star_luminocity), 0.25l);
	long double B = planck(nu, surf_t);
	long double one_div_cos_gamma = 1.0l; // =1.0l /cos(gamma_incl);
	long double F = (cell.s_ / d_sq) * B* (1.0l - exp(-cell.sigma_ * k_lambda_precalculated[n] * (one_div_cos_gamma))); // 14 formula (planck+ t_surf)
	return F;
}

void discOutputRoutine(long double spectre[]){
	ofstream fout("outer_result.dat");
	if(fout.is_open()){
    	cout << "File Opened successfully. Writing data from array to file" << endl;

		for(int i = 0; i < precision; ++i){
      		fout << wavelengths[i] << " " << spectre[i] << endl; //writing ith character of array in the file      		
      		//fout << wavelengths[i] << " " << ((spectre[i] < 1e-11l) ? (0) : (spectre[i])) << endl;
		}
    	cout << "Array data successfully saved into the file result.dat" << endl;
	}
}

void discRoutine(long double spectre[] ){
	//parallel_for()
	for (Cell & cell : cells_outer){
		cout << cell.line_index_ << "\r" << flush;
		long double lambda = lambda_min;
		for (int i = 0; i < precision; ++i){
			lambda = wavelengths[i];
			long double frequency = C/(lambda*0.0001); // cm/s /cm = 1/s

			long double counting_result = calcOuterDiscCell(cell, frequency, i);
			//cout << "discRoutine(): counting_result=" << counting_result << endl;
			if (counting_result == counting_result && !std::isinf(counting_result)){
				spectre[i] = spectre[i] + counting_result*frequency;
			}
		}		
	}
	discOutputRoutine(spectre);
}
//#####################################################################################
void preparation(){
	//Reading all nesessary data from .dat files
	readData();
	calculate_logscale_waves();
	calculate_Klambda();
}

int main(int argc, char **argv){
	preparation();

	//Initialising output spectre and nullify it
	long double spectre[(int)precision];
	for (int i = 0; i < precision; ++i){
		spectre[i] = 0;
	}

	//All actions related to star
	//starRoutine();
	//All actions related to outer disc
	discRoutine(spectre);
	
	
	cout << "done" << endl;

	return 0;
}
