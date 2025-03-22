// SOLUTION OF SYSTEM OF ODE'S WITH 4TH ORDER RUNGE-KUTTA-FEHLBERG METHOD

// DATE: 25/6/2021
// AUTHOR: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;


double f(double x, double y){				// ODE: y' = f(x,y)
	return -2*x-y;
}

int main(){
	cout << setprecision(9);
	
	double h = 0.1;				// Lenght of the interval of integration
	double x0 = 0.0, xf = 5.0;	// Interval of integration
	double y0 = -1.0;			// Initial value of y(x)
	
	int n = round((xf-x0)/h);	// Number of intervals
	double x[n+1];				// Vector with values of x at the end of each interval
	double y4[n+1];				// Vector with solution y with global error O(h^4)
	double y5[n+1];				// Vector with solution y with global error O(h^5)
	y5[0] = y0;
	// std::cout << "Number of intervals: " << n <<"\n\n";
	for (int i=0 ; i<=n ; i++){
		x[i] = x0 + i*h;
	}
	
	cout << "4TH ORDER RUNGE-KUTTA-FEHLBERG METHOD:\n";
	double k1, k2, k3, k4, k5, k6;			// Auxiliar values
	double E;			// Error
	for (int i=0 ; i<=n-1 ; i++){
		k1 = h*f(x[i],y5[i]);
		k2 = h*f(x[i] + 0.25*h, y5[i] + 0.25*k1);
		k3 = h*f(x[i] + 3.0*h/8.0, y5[i] + 3.0*k1/32.0 + 9.0*k2/32.0);
		k4 = h*f(x[i] + 12.0*h/13.0, y5[i] + 1932.0*k1/2197.0 - 7200.0*k2/2197.0 + 7296.0*k3/2197.0);
		k5 = h*f(x[i] + h, y5[i] + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 - 845.0*k4/4104.0);
		k6 = h*f(x[i] + 0.5*h, y5[i] -8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0);
		y4[i+1] = y5[i] + 25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0;
		y5[i+1] = y5[i] + 16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0 - 9.0*k5/50.0 + 2.0*k6/55.0;
		E = k1/360.0 - 128.0*k3/4275.0 - 2197.0*k4/75240.0+k5/50.0+2.0*k6/55.0;
		// cout << k1<<"\n"<< k2<<"\n"<< k3<<"\n"<< k4<<"\n"<< k5<<"\n"<< k6<<"\n"<< y4[1]<<"\n"<< y5[1]<<"\n"<<E<<"\n";
	}
	
	std::cout << "\nx_n=[";
	for (int i=0 ; i<n ; i++){
		std::cout << x[i] << ", ";
	}
	cout << x[n] << "];\ny_n=[";
	for (int i=0 ; i<n ; i++){
		std::cout << y5[i] << ", ";
	}
	cout << y5[n] << "];";
	
	// Write the solution in a file
	// ofstream file;
	// file.open ("P14_4th_runge_kutta.txt");
	// file << "4TH ORDER RUNGE-KUTTA METHOD:\n";
	
	// file.close();
	
	return 0;
}
