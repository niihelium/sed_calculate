#ifndef DATA_READER_H
#define DATA_READER_H

#include "cell.h"

class DataReader
{
	private:
		static const std::string grid;
		static const std::string var030;
		static const std::string opacities;
		static const std::string Arate;
	public:
		static void readOuterData(double &system_age, std::list<Cell> &cells_outer); 

		static void readKnu(std::list<double> &k_nu_temp, std::list<double> &k_nu);

		static double readArate(double system_age, double &star_luminocity);


};
#endif