// BISECTION METHOD FOR NON-LINEAR EQUATIONS

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
	double a = 2;
	double b = 3;
	
	if (f(a)*f(b) > 0){
		std::cout << "ERROR: Bisection method cannot be applied if f(a)*f(b) > 0.";
		exit(0);		// ends the program
	}
	else if (f(a)*f(b) == 0){
		if (f(a) == 0){
			std::cout << a << " is a solution to the equation f(x) = 0";
			exit(0);		// ends the program
		}
		else if (f(b) == 0){
			std::cout << b << " is a solution to the equation f(x) = 0";
			exit(0);		// ends the program
		}
	}
	
	double TOL = 1.0e-9;			// tolerance of the solution's precision
	double c = (a+b)/2;				// Bisection approximation
	double eps = fabs(b-c); 			// precision of the 3rd approximation
		
	while (eps > TOL){
		if (f(a)*f(c) < 0)
			b = c;
		else if (f(a)*f(c) > 0)
			a = c;
		if (f(c) == 0){
			std::cout << c << " is a solution to the equation f(x) = 0";
			break;		// ends the while loop
		}
		c = (a+b)/2;
		eps = fabs(b-c); 			// precision of the (i+1)-th approximation
	}
	
	// Print the solution:
	std::cout << "x = " << c << " is a solution to the equation f(x) = 0.";
	std::cout << "\nNote: The method may produce a false root if f(x) is discontinuous in [x_{1},x_{2}].";
	
	return 0;
}
