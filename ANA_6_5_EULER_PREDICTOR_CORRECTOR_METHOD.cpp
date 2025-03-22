// SOLUTION OF ODE USING EULER PREDICTOR-CORRECTOR METHOD

// DATE: 25/6/2021
// AUTHOR: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;


double f(double t, double x, double y){				// ODE x' = f(t,x,y)
	return x*y+t;
}
double g(double t, double x, double y){				// ODE y' = g(t,x,y)
	return t*y+x;
}

int main(){
	cout << setprecision(9);
	
	double h = 0.1;				// Lenght of the interval of integration
	double t0 = 0.0, tf = 0.1;	// Interval of integration
	double x0 = 1.0;			// Initial value of x(t)
	double y0 = -1.0;			// Initial value of y(t)
	
	int n = round((tf-t0)/h);	// Number of intervals
	double t[n+1];				// Vector with values of t at the end of each interval
	double x[n+1],  			// Vectors with solution x(t), y(t)
	y[n+1];
	x[0] = x0;
	y[0] = y0;
	// std::cout << "Number of intervals: " << n <<"\n\n";
	for (int i=0 ; i<=n ; i++){
		t[i] = t0 + i*h;
	}
	
	// EULER PREDICTOR-CORRECTOR METHOD:
	for (int i=0 ; i<=n-1 ; i++){
		x[i+1] = x[i] + f(t[i],x[i],y[i])*h;
		cout << x[i+1]<<endl;
		y[i+1] = y[i] + g(t[i],x[i],y[i])*h;
		cout << y[i+1]<<endl;
		x[i+1] = x[i] + (f(t[i],x[i],y[i]) + f(t[i+1],x[i+1],y[i+1]))*h/2;
		cout << x[i+1]<<endl;
		y[i+1] = y[i] + (g(t[i],x[i],y[i]) + g(t[i+1],x[i+1],y[i+1]))*h/2;
		cout << y[i+1]<<endl;
	}
	
	cout << "EULER PREDICTOR-CORRECTOR METHOD:\n";
	
	std::cout << "\nt_n=[";
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << t[i] << ", ";
	}
	cout << t[n] << "];\nx_n=[";
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << x[i] << ", ";
	}
	cout << x[n] << "];\ny_n=[";
	for (int i=0 ; i<=n-1 ; i++){
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
