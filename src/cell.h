#ifndef CELL_H
#define CELL_H
 
class Cell
{
private:
    Cell() { }// private default constructor
 
public:
	Cell(int line_index, long double r, long double phi, long double s, long double sigma, long double T_midplane, long double tau);

	int line_index_;

	//varXXX.dat
	long double sigma_;
	long double T_midplane_;
	long double tau_;

	//grid.dat
	long double r_;
	long double phi_;
	long double s_;
};
 
#endif