// SOLUTION OF ODE WITH 4TH ORDER RUNGE-KUTTA-MERSON METHOD

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
	
	double h = 0.1;				// Lenght of the interval of integration
	double x0 = 0.0, xf = 5.0;	// Interval of integration
	double y0 = -1.0;			// Initial value of y(x)
	
	int n = round((xf-x0)/h);	// Number of intervals
	double x[n+1];				// Vector with values of x at the end of each interval
	double y[n+1];				// Vector with solution y with global error O(h^4)
	y[0] = y0;
	// std::cout << "Number of intervals: " << n <<"\n\n";
	for (int i=0 ; i<=n ; i++){
		x[i] = x0 + i*h;
	}
	
	cout << "4TH ORDER RUNGE-KUTTA-FEHLBERG METHOD:\n";
	double k1, k2, k3, k4, k5;			// Auxiliar values
	double E;			// Error
	for (int i=0 ; i<=n-1 ; i++){
		k1 = h*f(x[i],y[i]);
		k2 = h*f(x[i] + h/3.0, y[i] + k1/3.0);
		k3 = h*f(x[i] + h/3.0, y[i] + k1/6.0 + k2/6.0);
		k4 = h*f(x[i] + h/2.0, y[i] + k1/8.0 + 3.0*k3/8.0);
		k5 = h*f(x[i] + h, y[i] + k1/2.0 - 3.0*k3/2.0 + 2.0*k4);
		y[i+1] = y[i] + (k1 + 4.0*k4 + k5)/6.0;
		E = (2.0*k1 - 9.0*k3 + 8.0*k4 - k5)/30.0;
		// cout << k1<<"\n"<< k2<<"\n"<< k3<<"\n"<< k4<<"\n"<< k5<<"\n"<< k6<<"\n"<< y4[1]<<"\n"<< y5[1]<<"\n"<<E<<"\n";
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
	// file.open ("P14_4th_runge_kutta.txt");
	// file << "4TH ORDER RUNGE-KUTTA METHOD:\n";
	
	// file.close();
	
	return 0;
}
