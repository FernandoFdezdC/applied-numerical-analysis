// FIXED-POINT METHOD FOR NON-LINEAR EQUATIONS

// Date: 5/11/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;

double f(double x){			// function f(x) whose roots are to be found
	return cos(x);
}

int main(){
	cout << setprecision(9);
	// Equation to be solved: f(x) = 0
	std::cout << "f(x) = cos(x)" << endl;
	
	int n = 0;			// counter of iterations
	
	// first try:
	double a = 0.7;
	
	a = f(a);
	double eps; 		// precision of the (i+1)-th approximation
	double TOL = 1.0e-8;		// tolerance of the solution's precision
		
	do {
		n++;
		eps = fabs(f(a)-a);
		a = f(a);
		cout << a <<"\n";
	} while (eps > TOL);
	
	// Print the solution
	std::cout << "x=" << a << " is a solution to the equation f(x) = x";
	std::cout << "\nNecessary iterations: " << n;
	
	return 0;
}
