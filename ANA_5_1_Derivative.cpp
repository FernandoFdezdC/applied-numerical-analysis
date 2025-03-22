// NUMERICAL DERIVATIVE

// Date: 28/11/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
using namespace std;

double f(double x){				// Function to integrate
	return sin(x)*exp(x);
}
double dfdx(double x){
	return exp(x)*(sin(x)+cos(x));
}

int main(){
	cout << setprecision(9);	// Number of digits in the output
	
	double x0 = 1.9;			// Point to evaluate f'(x) and f''(x)
	double h = 0.05;			// step
	
	// FIRST DERIVATIVE:
	double FD = (f(x0+h)-f(x0))/h;					// + O(h)
	double BD = (f(x0)-f(x0-h))/h;					// + O(h)
	double CF = (f(x0+h)-f(x0-h))/(2*h);			// + O(h^2)
	double AFD = (-3*f(x0)+4*f(x0+h)-f(x0+2*h))/(2*h);	// + O(h^2)
	double ABD = (3*f(x0)-4*f(x0-h)+f(x0-2*h))/(2*h);	// + O(h^2)
	cout << "f'(x)=\n";
	cout << "FORWARD DIFFERENCE FORMULA: " << FD << " + O(h)\n";
	cout << "BACKWARD DIFFERENCE FORMULA: " << BD << " + O(h)\n";
	cout << "CENTRAL DIFFERENCE FORMULA: " << CF << " + O(h^2)\n";
	cout << "ACCURATE FORWARD DIFFERENCE FORMULA: " << AFD << " + O(h^2)\n";
	cout << "Real value: "<< dfdx(x0);
	
	// SECOND DERIVATIVE:
	double FD2 = (f(x0+2*h)-2*f(x0+h)+f(x0))/(pow(h,2));	// + O(h)
	double BD2 = (f(x0)-2*f(x0-h)+f(x0-2*h))/(pow(h,2));	// + O(h)
	double CD2 = (f(x0+h)+f(x0-h)-2*f(x0))/(pow(h,2));		// + O(h^2)
	cout << "\nf''(x)=\n";
	cout << "FORWARD DIFFERENCE FORMULA: " << FD2 << " + O(h)\n";
	cout << "BACKWARD DIFFERENCE FORMULA: " << BD2 << " + O(h)\n";
	cout << "CENTRAL DIFFERENCE FORMULA: " << CD2 <<  " + O(h^2)";
	
	return 0;
}
