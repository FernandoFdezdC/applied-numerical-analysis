// METHOD OF FALSE POSITION (REGULA FALSI) FOR NON-LINEAR EQUATIONS

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
	
	// first tries:
	double x0 = 2.0;
	double x1 = 2.1;
	double x2;		// next iteration
	
	if (f(x0)*f(x1) > 0){
		std::cout << "ERROR: Regula Falsi method cannot be applied if f(a)*f(b) > 0.";
		exit(0);		// Ends the program
	}
	else if (f(x0)*f(x1) == 0){
		if (f(x0) == 0){
			std::cout << x0 << " is a solution to the equation f(x) = 0";
			exit(0);		// Ends the program
		}
		else if (f(x1) == 0){
			std::cout << x1 << " is a solution to the equation f(x) = 0";
			exit(0);		// Ends the program
		}
	}
	
	double TOL = 1.0e-5;		// precision of the solution
	
	do {
		x2 = (x0*f(x1)-x1*f(x0))/(f(x1)-f(x0));	 // Regula Falsi approximation
		if (f(x0)*f(x2) < 0)
			x1 = x2;
		else if (f(x0)*f(x2) > 0)
			x0 = x2;
	} while (fabs(f(x2)) > TOL);
	
	// Print the solution
	std::cout << "\nx = " << x2 << " is a solution to the equation f(x) = 0";
	std::cout << "\nNote: If f(x) is not continuous, the method may fail.";
	
	return 0;
}
