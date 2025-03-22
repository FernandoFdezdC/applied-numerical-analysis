// SECANT METHOD FOR NON-LINEAR EQUATIONS

// Date: 5/11/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;

double f(double x){			// function f(x) whose roots are to be found
	double y = log(x) - exp(-x) + pow(x,2) - 5;
  	return y;
}

int main(){
	cout << setprecision(9);
	
	// Equation to be solved: f(x) = 0
	std::cout << "f(x) = ln(x) - exp(-x) + x^2 - 5" << endl;
	
	int i = 2; 			// index of the iteration
	
	double TOL = 1.0e-5;			// tolerance of the solution's precision
	
	// first tries:
	double x0 = 2.0;
	double x1 = 2.2;
	double x2;		// next iteration
	
	// It's better if |f(x_{0})| > |f(x_{1})|
	if (fabs(f(x0)) < fabs(f(x1))){
		double temp;
		temp = x0;
    	x0 = x1;
    	x1 = temp;
	}
	do {
		x2 = x1 - (x0-x1)/(f(x0)-f(x1))*f(x1);	// Secant method approximation
		x0=x1;
		x1=x2;
		i++;
	} while (fabs(f(x2)) > TOL);

	// Print the solution
	std::cout << x2 << " is a solution to the equation f(x) = 0.";
	std::cout << "\nNumber of iterations made: n=" << i;
	std::cout << "\nNote: If f(x) is not continuous, the method may fail.";
	
	return 0;
}
