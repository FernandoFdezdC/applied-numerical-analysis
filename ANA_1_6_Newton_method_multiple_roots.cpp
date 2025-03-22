// NEWTON'S METHOD FOR NON-LINEAR EQUATIONS
// WHEN THE FUNCTION HAS A MULTIPLE ROOT

// Date: 11/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;

double f(double x){			// function f(x) whose roots are to be found
	return (x-1)*(exp(x-1)-1);
}

double dfdx(double x){		// derivative of function f(x)
	return x*exp(x-1)-1;
}

int main(){
	cout << setprecision(9);
	
	// Equation to be solved: f(x) = 0
	std::cout << "f(x) = 3x + sin(x) - e^x" << endl;

	int n = 1;			// counter of iterations
	int k=2; 		// multiplicity of the root
	
	double x0=2.0,   x1;	// first try, x0
	double TOL=1e-4;		// required tolerance

	double diff;			// difference between successive values of x
	do {
		if (f(x0) == 0 || dfdx(x0) == 0){
			std::cout << "ERROR: Newton's method cannot be applied if f(a)*f(b) > 0.";
			exit(0);		// Ends the program
		}
		x1 = x0 - k*f(x0)/dfdx(x0);
		diff = x1-x0;
		cout<<x0<<"\n";
		x0=x1;
		n++;
	} while(fabs(diff)>TOL || fabs(f(x1))>TOL);
	
	// Print the solution
	std::cout << x1 << " is a solution to the equation f(x) = 0.";
	std::cout << "\nNecessary iterations: " << n;
	std::cout << "\nNote: The method may converge to a root different from the expected one or diverge if the starting value is not close enough to the root.";
	
	return 0;
}
