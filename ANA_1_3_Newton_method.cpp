// NEWTON'S METHOD FOR NON-LINEAR EQUATIONS

// Date: 28/2/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;

double f(double x){			// function f(x) whose roots are to be found
	return 3*x + sin(x)-exp(x);
}

double dfdx(double x){		// derivative of function f(x)
	return 3+cos(x)-exp(x);
}

int main(){
	cout << setprecision(9);
	
	// Equation to be solved: f(x) = 0
	std::cout << "f(x) = 3x + sin(x) - e^x" << endl;

	int n = 1;			// counter of iterations
	
	double x0=2.0,   x1;	// first try, x0;
	double TOL=1e-5;		// required tolerance

	double diff;			// difference between successive values of x
	do {
		if (f(x0) == 0 || dfdx(x0) == 0){
			std::cout << "ERROR: Newton's method cannot be applied if f(a)*f(b) > 0.";
			exit(0);		// Ends the program
		}
		x1 = x0 - f(x0)/dfdx(x0);
		diff = x1-x0;
		x0=x1;
		n++;
	} while(fabs(diff)>TOL || fabs(f(x1))>TOL);
	
	// Print the solution
	std::cout << x1 << " is a solution to the equation f(x) = 0.";
	std::cout << "\nNecessary iterations: " << n;
	std::cout << "\nNote: The method may converge to a root different from the expected one or diverge if the starting value is not close enough to the root.";
	
	return 0;
}
