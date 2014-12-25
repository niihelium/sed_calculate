#ifndef CELL_H
#define CELL_H
 
class Cell
{
private:
    Cell() { }// private default constructor
 
public:
	Cell(int line_index, double r, double phi, double s, double sigma, double T_midplane, double theta, double t_surface);

	int line_index_;

	//varXXX.dat
	double sigma_;
	double T_midplane_;
	double theta_;

	//grid.dat
	double r_;
	double phi_;
	double s_;

	//calculated data
	double t_surface_;
};
 
#endif