// SOLUTION OF ODE'S WITH IMPROVED EULER'S METHOD

// DATE: 24/6/2021
// AUTHOR: Fernando Fernández del Cerro

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
	
	std::cout << "IMPROVED EULER'S METHOD:\n";
	
	for (int i=0 ; i<=n-1 ; i++){
		y[i+1] = y[i] + f(x[i],y[i])*h;
		y[i+1] = y[i] + (f(x[i],y[i])+f(x[i+1],y[i+1]))*h/2;
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
	// file.open ("P14_euler.txt");
	// file << "EULER'S METHOD:\n";
	// file.close();
	
	return 0;
}
