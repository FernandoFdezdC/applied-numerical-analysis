// SOLUTION OF ODE'S WITH 3RD ORDER RUNGE-KUTTA METHOD

// DATE: 9/12/2020
// AUTHOR: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;


double f(double x, double y){				// ODE y' = f(x,y)
	double f = y - x;
	return f;
}

int main(){
	cout << setprecision(9);
	
	double h = 0.1;				// Lenght of the intervals to integrate
	double x0 = 0.0, xf = 0.8;	// Interval of integration
	double y0 = 0.5;			// Initial value of y(x)
	
	int n = round((xf-x0)/h);	// Number of intervals
	double y[n];			// Solution y
	y[0] = y0;
	
	cout << "3RD ORDER RUNGE-KUTTA METHOD:\n";
	double k1, k2, k3;			// Auxiliar values
	
	for (int i=0 ; i<n ; i++){
		k1 = f(x0+i*h,y[i]);
		k2 = f(x0+i*h+0.5*h,y[i]+k1*0.5*h);
		k3 = f(x0+i*h+h,y[i]-k1*h+2*k2*h);
		y[i+1] = y[i] + (k1+4*k2+k3)*h/6;
		// cout << "y[" << i << "] = " << y[i] << endl;		// Checking
	}
	
	double yf = y[n];
	cout << "y(" << xf << ") = " << yf;
	
	// Write the solution in a file
	// ofstream file;
	// file.open ("P14_3rd_runge_kutta.txt");
	// file << "3RD ORDER RUNGE-KUTTA METHOD:\n";
	
	// file.close();
	
	return 0;
}
