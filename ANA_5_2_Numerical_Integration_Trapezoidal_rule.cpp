// NUMERICAL INTEGRATION
// NEWTON-COTES FORMULAE
// TRAPEZOIDAL RULE

// Date: 28/11/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;


double f(double x){				// Function to integrate
	double f = sin(x);
	return f;
}

int main(){
	cout << setprecision(9);	// Number of digits in the output
	
	double I = 0;		// Integral of f(x) from a to b
	double a = 0;		// Lower value
	double b = M_PI;	// Upper value
	int n = 4; 			// Number of divisions for the integration
	
	int i = 1;			// x_i
	while (i<=n) {
		I = I + (b-a)*(f(a+(b-a)*(i-1)/n)+f(a+(b-a)*i/n))/(2*n);
		i++;
	}
	
	cout << "TRAPEZOIDAL RULE:\n";
	cout << "Integral of f(x) from " << a << " to " << b << " = " << I;
   	
	return 0;
}
