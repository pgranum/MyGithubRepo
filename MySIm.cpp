#include <omp.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <vector>
#include <sstream>
#include <ctime>
#include <iomanip>

#include "TVector3.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h" // Joe enabled my makefile to link to ROOT
#include "Rtypes.h" // Root data types
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TPolyMarker3D.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h> // some matrix and vector operations


//using namespace std; // not allowed acc to Joe
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

void printGslMatrix(gsl_matrix * A){
    int n = A->size1;
    int m = A->size2;
    std::cout << std::endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            double A_ij = gsl_matrix_get(A,i,j);
            if (abs(A_ij) < 1e-6){
                std::cout << 0 << "\t";
            } else{
                std::cout << std::setprecision(3) << A_ij << "\t";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
void printGslVector(const gsl_vector * a){
    int n = a->size;
    std::cout << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << std::setprecision(5) << gsl_vector_get(a,i) << std::endl;
    }
}
void setGslVector(gsl_vector * a, double x, double y, double z){
    gsl_vector_set(a,0,x);
    gsl_vector_set(a,1,y);
    gsl_vector_set(a,2,z);
}
void setVecInMatrix(gsl_matrix * A, int i, gsl_vector * a){
    int n = A->size1;
    for (int j = 0; j < n; ++j) {
        gsl_matrix_set(A,j,i,gsl_vector_get(a,j));
    }
}
void plotGslMatrix(gsl_matrix * A, TPolyMarker3D *pm3d){
    int n = A->size1;
    int m = A->size2;
    if (n != 3) {
        std::cout << "error in plotGslMatrix. n != 3"
    }
    else {
        for (int i = 0; i < m; ++i) {
            pm3d->SetPoint( i, gsl_matrix_get(wire1.wirePoints()), y, z )
        }
    }
}

double fRand(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

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
		gsl_vector * startPoint = gsl_vector_alloc(3);
		gsl_vector * endPoint = gsl_vector_alloc(3);

		// methods
        gsl_matrix * wirePoints(void);
};
gsl_matrix * StraightWire::wirePoints (void){
	gsl_vector * wireVector = gsl_vector_alloc(3);
    gsl_vector_memcpy(wireVector,endPoint);
    gsl_vector_sub(wireVector,startPoint); // endPoint - startPoint
    printGslVector(wireVector);

    double norm = gsl_blas_dnrm2(wireVector);
    int n = ceil(norm/dx);
    const double scale = 1.0/n;

    gsl_vector_scale(wireVector,scale);
    printGslVector(wireVector);

    gsl_matrix * wirePoints = gsl_matrix_alloc(3,n+1);
    gsl_vector * temp = gsl_vector_alloc(3);
    gsl_vector_memcpy(temp,startPoint);
    setVecInMatrix(wirePoints,0,startPoint);

    for (int i = 0; i < n; ++i) {
        gsl_vector_add(temp,wireVector);
        setVecInMatrix(wirePoints,i+1,temp);
    }
    std::cout << "wirePoint called and finished" << std::endl; // debugging
    return wirePoints;
};

int	main(int argc, char *argv[]) { // .exe dx
	// some people use char **argv[] ...
	dx = atof(argv[1]); // string to double
	srand((unsigned)time(NULL)); // for the random number generator
	TApplication *app = new TApplication("app", 0, 0); // for showing TCanvas
	
	{ // test area
//	double x = 3;
//	Vector3D v;
//	v.x = 2;
//	Loop MA;
//	MA.centre.x = 4;
//	MA.radius = 1;
//	//~ std::vector<double> myVec = {0.5, 0.7, 1}; // making a vector with elements of type double.
//	//~ myVec.push_back(0.5) // the way to add elements to a vector.
//	double circ = MA.getCircumference();
//
//	StraightWire myWire;
//	myWire.startPoint.x = 3;
//	myWire.startPoint.y = 4;
//	double a = 5;
//	myWire.startPoint.setVector(a,a,a);
//
//	double array[3]{1.0,2.0,3.0};
//	std::cout << "array 0"<< array[0] << std::endl;
//
//	TVector3 vec;
//	vec(0) = 10;
//	vec(2) = 11;
//	vec(1) = 12;
//	{
//	int n = 5;
//    gsl_matrix * A = gsl_matrix_alloc(n,n); // generate a symmetric matrix
//    //~ delete a;
//    gsl_vector * a = gsl_vector_alloc(n);
//    double rand;
//    for (int i = 0; i < n; ++i) {
//        for (int j = i; j < n; ++j) {
//            rand = fRand(0.1,10);
//            gsl_matrix_set(A,i,j,rand);
//            gsl_matrix_set(A,j,i,rand);
//        }
//        gsl_vector_set(a,i,rand);
//    }
//
//    std::cout << "Generating a random symmetric matrix A" << std::endl;
//    printGslMatrix(A);
//    printGslVector(a);
//	}
//	vec.Print();
//	std::cout << "x times x = " << myWire.startPoint.x << std::endl;
//	std::cout << "x times x = " << myWire.startPoint.getLength() << std::endl;
//
//	{
//	TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
//	c1->SetGrid();
//	c1->GetFrame()->SetBorderSize(12);
//	const Int_t n = 10;
//	Float_t x[n]  = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
//	Float_t y[n]  = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
//	Float_t ex[n] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
//	Float_t ey[n] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
//    TGraphErrors *gr = new TGraphErrors(n,x,y,ex,ey);
//    gr->SetTitle("TGraphErrors Example");
//    gr->SetMarkerColor(4);
//    gr->SetMarkerStyle(21);
//    gr->Draw("ALP");
//    c1->Update();
//    c1->SaveAs("graph.png");
//	}

	} // test area

	StraightWire wire1;
    setGslVector(wire1.endPoint,0,0,5);
    setGslVector(wire1.startPoint,0,0,-5);
    printGslMatrix(wire1.wirePoints());

    TCanvas *c1 = new TCanvas("c1","wire graph",800,600);
    c1->SetGrid();
	c1->GetFrame()->SetBorderSize(12);
    TPolyMarker3D *pm3d = new TPolyMarker3D(11);

    for (int i = 0; i < 11; ++i) {
        pm3d->SetPoint( i, gsl_matrix_get(wire1.wirePoints()), y, z )
    }

//	app->Run(); // Running the TApplication
	return 0;
};


double sq(double x) {
	return x*x;
};
