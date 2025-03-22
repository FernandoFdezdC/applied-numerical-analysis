// NUMERICAL INTEGRATION
// NEWTON-COTES FORMULAE
// SIMPSON'S 3/8 RULE

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
	int n = 6; 			// number of divisions for the integration. Must be multiple of 3
	
	int i = 3;			// x_i
	
	while (i<=n) {
		I = I + 3*(b-a)*(f(a+(b-a)*(i-3)/n)+3*f(a+(b-a)*(i-2)/n)+3*f(a+(b-a)*(i-1)/n)+f(a+(b-a)*i/n))/(8*n);
		
		i=i+3;
	}
	
	cout << "SIMPSON'S' RULE:\n";
	cout << "Integral of f(x) from " << a << " to " << b << " = " << I;
   	
	return 0;
}
