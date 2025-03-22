// HORNER'S METHOD FOR POLYNOMIAL EQUATIONS

// Date: 10/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;

double f(double x){			
	return 3*x + sin(x)-exp(x);
}

int main(){
	cout << setprecision(9);
	
	// Equation to be solved: f(x) = 0
	std::cout << "f(x) = 2x^3+x^2-3x-3" << endl;
	
	int n = 3;			// degree of polynomial f(x)
	
	double f[n+1] = {-3.0,-3.0,1.0,2.0};	// Coefficients of the polynomial f(x)
		// the leftmost coefficient is the constant coefficient and
		// the rightmost one if the coefficient of x^n
	
	int i=n;
	double x0=2.0;
	double r[n+1];		// vector with residues
	r[n] = f[n];
	
	// Value of the polynomial f(x) at x0:
	while (i>=1){
		r[i-1] = r[i]*x0+f[i-1];
		// cout << r[i]<<"\n";
		i--;
	}
	
	// Print the solution
	std::cout << "Polynomial f(x_{0}="<<x0<<") = "<< r[0];

	return 0;
}
