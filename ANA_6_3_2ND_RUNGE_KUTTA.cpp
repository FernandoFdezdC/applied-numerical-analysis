// SOLUTION OF ODE'S WITH 2ND ORDER RUNGE-KUTTA METHOD

// DATE: 9/12/2020
// AUTHOR: Fernando Fern�ndez del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;


double f(double x, double y){				// ODE y' = f(x,y)
	return -2*x-y;
}

int main(){
	cout << setprecision(9);
	
	double h = 0.1;				// Lenght of the interval of integration
	double x0 = 0.0, xf = 5.0;	// Interval of integration
	double y0 = -1.0;			// Initial value of y(x)
	
	int n = round((xf-x0)/h);	// Number of intervals
	double x[n+1];				// Vector with values of x at the end of each interval
	double y[n+1];				// Vector with solution y
	y[0] = y0;
	// std::cout << "Number of intervals: " << n <<"\n\n";
	for (int i=0 ; i<=n ; i++){
		x[i] = x0 + i*h;
	}
	
	// a+b=1, \alpha b = \frac{1}{2}, \beta b = \frac{1}{2} must hold
	double b = pow(3.0,-1); 		// Optimal value for b
	// double b = 0.5; 				// Arbitrary value for b
	double a = 1-b, alpha = 1/(2*b), beta = 1/(2*b);
	cout << a << "\n" << alpha << "\n" << beta << "\n";
	
	std::cout << "2ND ORDER RUNGE-KUTTA METHOD:\n";
	double k1, k2;			// Auxiliar values
	
	for (int i=0 ; i<=n-1 ; i++){
		k1 = h*f(x[i],y[i]);
		k2 = h*f(x[i]+alpha*h,y[i]+beta*k1);
		y[i+1] = y[i] + a*k1 + b*k2;
	}
	
	std::cout << "\nx_n=[";
	for (int i=0 ; i<n ; i++){
		std::cout << x[i] << ", ";
	}
	cout << x[n] << "];\ny_n=[";
	for (int i=0 ; i<n ; i++){
		std::cout << y[i] << ", ";
	}
	cout << y[n] << "];";
	
	// Write the solution in a file
	// ofstream file;
	// file.open ("P14_2nd_runge_kutta.txt");
	// file << "2ND ORDER RUNGE-KUTTA METHOD:\n";
	
	// file.close();
	
	return 0;
}
