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
		static void readOuterData(long double &system_age, std::list<Cell> &cells_outer); 

		static void readKlambda(std::list<long double> &k_nu_temp, std::list<long double> &k_nu);

		static long double readArate(long double system_age, long double &star_luminocity, long double &star_radius);


};
#endif