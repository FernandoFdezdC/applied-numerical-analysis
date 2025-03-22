// SOLUTION OF ODE BEGGINING WITH 4TH ORDER RUNGE-KUTTA-FEHLBERG
// METHOD AND CONTINUING WITH ADAMS-MOULTON METHOD

// DATE: 25/6/2021
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
	
	double h = 0.2;				// Lenght of the interval of integration
	double x0 = 0.0, xf = 0.6;	// Interval of integration
	double y0 = -1.0;			// Initial value of y(x)
	
	int n = round((xf-x0)/h);	// Number of intervals
	double x[n+1];				// Vector with values of x at the end of each interval
	double y[n+1];				// Vector with solution y with global error O(h^4)
	y[0] = y0;
	// std::cout << "Number of intervals: " << n <<"\n\n";
	for (int i=0 ; i<=n ; i++){
		x[i] = x0 + i*h;
	}
	
	// 4TH ORDER RUNGE-KUTTA-FEHLBERG METHOD:
	double k1, k2, k3, k4, k5, k6;			// Auxiliar values
	double E;			// Error
	for (int i=0 ; i<=2 ; i++){
		k1 = h*f(x[i],y[i]);
		k2 = h*f(x[i] + 0.25*h, y[i] + 0.25*k1);
		k3 = h*f(x[i] + 3.0*h/8.0, y[i] + 3.0*k1/32.0 + 9.0*k2/32.0);
		k4 = h*f(x[i] + 12.0*h/13.0, y[i] + 1932.0*k1/2197.0 - 7200.0*k2/2197.0 + 7296.0*k3/2197.0);
		k5 = h*f(x[i] + h, y[i] + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 - 845.0*k4/4104.0);
		k6 = h*f(x[i] + 0.5*h, y[i] -8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0);
		y[i+1] = y[i] + 16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0 - 9.0*k5/50.0 + 2.0*k6/55.0;
		E = k1/360.0 - 128.0*k3/4275.0 - 2197.0*k4/75240.0+k5/50.0+2.0*k6/55.0;
	}
	
	// ADAMS METHOD:
	for (int i=3 ; i<=n-1 ; i++){
		y[i+1] = y[i] + h*(55.0*f(x[i],y[i]) - 59.0*f(x[i-1],y[i-1]) + 37.0*f(x[i-2],y[i-2]) - 9.0*f(x[i-3],y[i-3]))/24.0;  // Predictor
		y[i+1] = y[i] + h*(9.0*f(x[i+1],y[i+1]) + 19.0*f(x[i],y[i]) - 5.0*f(x[i-1],y[i-1]) + f(x[i-2],y[i-2]))/24.0;  // Corrector
	}
	
	cout << "ADAMS-MOULTON METHOD:\n";
	
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
	// file.open ("P14_4th_runge_kutta.txt");
	// file << "4TH ORDER RUNGE-KUTTA METHOD:\n";
	
	// file.close();
	
	return 0;
}
