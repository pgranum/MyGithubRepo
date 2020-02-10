#include <omp.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <vector>
#include <sstream>
#include "TVector.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

//using namespace std;
const double PI=3.1415926535897932384626433832795;
const double mu0=4.*PI*1.0e-7;
double dx;

double sq(double x); // define a function in scope. Actual def further down

class Vector3D {
	public:
		// properties
		double x,y,z;
		// methods
		double getLength(void);
		void setVector(double a, double b, double c);
		
};
double Vector3D::getLength(void){
	return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
};
void Vector3D::setVector(double a, double b, double c){
	x = a; y = b; z = c;
	return;
};


class Magnet {
	public:
		double current; // in Amps
};

class Loop : public Magnet {
	public:
		// properties:
		Vector3D centre;
		double radius;
		double width;

		// methods:
		double getCircumference(void);

};
double Loop::getCircumference(void){
	// this is a function of a class, indicated by the ::
	// it has access to the class variables
	return PI*2*radius;
};

class StraightWire : public Magnet {
	public:
		// properties
		Vector3D startPoint;
		Vector3D endPoint;

		// methods
		std::vector<Vector3D> wirePoints(void);
};
std::vector<Vector3D> StraightWire::wirePoints (void){
	Vector3D wireVec = startPoint; // endPoint - startPoint;
	int n = floor(wireVec.getLength() + dx); // change to /
	std::cout << n << std::endl;
	std::vector<Vector3D> wirePoints;
	wirePoints.push_back(wireVec);
	return wirePoints;
};

int	main(int argc, char *argv[]) { // .exe dx
	// some people use char **argv[] ...
	dx = atof(argv[1]); // string to double
	double x = 3;
	Vector3D v;
	v.x = 2;
	Loop MA;
	MA.centre.x = 4;
	MA.radius = 1;
	//~ std::vector<double> myVec = {0.5, 0.7, 1}; // making a vector with elements of type double.
	//~ myVec.push_back(0.5) // the way to add elements to a vector.
	double circ = MA.getCircumference();
	double test = MA.getCircumference();
	StraightWire myWire;
	myWire.startPoint.x = 3;
	myWire.startPoint.y = 4;
	double a = 5;
	myWire.startPoint.setVector(a,a,a);
	double array[3]{1.0,2.0,3.0};
	std::cout << "array 0"<< array[0] << std::endl;
	
	std::cout << "x times x = " << myWire.startPoint.x << std::endl;
	std::cout << "x times x = " << myWire.startPoint.getLength() << std::endl;
	return 0;
};


double sq(double x) {
	return x*x;
};
