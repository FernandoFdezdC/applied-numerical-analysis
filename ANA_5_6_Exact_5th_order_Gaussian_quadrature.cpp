// NUMERICAL INTEGRATION
// NEWTON-COTES FORMULAE
// GAUSSIAN QUADRATURE

// Date: 28/11/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int n = 5;			// Number of points for the gaussian quadrature

double f(double x){								// Function to integrate
	double f = pow(1-pow(sin(5*M_PI/360)*sin(x),2),-0.5);
	return f;
}

void x_w(double x[n], double w[n], int n){		// Set points and weights
	if (n == 1){
		x[0] = 0.0;
		w[0] = 2.0;
	}
	else if (n == 2){
		x[0] = -pow(3.0,-0.5);
		x[1] = pow(3.0,-0.5);
		w[0] = 1;
		w[1] = 1;
	}
	else if (n == 3){
		x[0] = -pow(3.0/5.0,0.5);
		x[1] = 0;
		x[2] = pow(3.0/5.0,0.5);
		w[0] = 5.0/9.0;
		w[1] = 8.0/9.0;
		w[2] = 5.0/9.0;
	}
	else if (n == 4){
		x[0] = -pow(3.0/7.0+2.0*pow(6.0/5.0,0.5)/7.0,0.5);
		x[1] = -pow(3.0/7.0-2.0*pow(6.0/5.0,0.5)/7.0,0.5);
		x[2] = pow(3.0/7.0-2.0*pow(6.0/5.0,0.5)/7.0,0.5);
		x[3] = pow(3.0/7.0+2.0*pow(6.0/5.0,0.5)/7.0,0.5);
		w[0] = (18.0-pow(30,0.5))/36.0;
		w[1] = (18.0+pow(30,0.5))/36.0;
		w[2] = (18.0+pow(30,0.5))/36.0;
		w[3] = (18.0-pow(30,0.5))/36.0;
	}
	else if (n == 5){
		x[0] = -pow(5+2*pow(10.0/7.0,0.5),0.5)/3;
		x[1] = -pow(5-2*pow(10.0/7.0,0.5),0.5)/3;
		x[2] = 0;
		x[3] = pow(5-2*pow(10.0/7.0,0.5),0.5)/3;
		x[4] = pow(5+2*pow(10.0/7.0,0.5),0.5)/3;
		w[0] = (322-13*pow(70,0.5))/900;
		w[1] = (322+13*pow(70,0.5))/900;
		w[2] = 128.0/225;
		w[3] = (322+13*pow(70,0.5))/900;
		w[4] = (322-13*pow(70,0.5))/900;
	}
}

int main(){
	cout << setprecision(9);	// Number of digits in the output
	
	double I = 0;				// Integral of f(x) from a to b
	double a = 0;				// Lower value
	double b = M_PI/2;			// Upper value
	cout << setprecision(9);
	cout << "GAUSSIAN QUADRATURE WITH N = " << n << " DIVISIONS:\n";
	cout << "f(x) = 1/sqrt(1-(sin(5*PI/360)*sin(x))^2)\n";
	// Set points and weights
	double x[n], w[n];
	x_w(x,w,n);
	
	cout << "Points -- Weights:\n";
	
	// Calculate integral I:
	for(int i=0 ; i<=n-1 ; i++){
		I = I + w[i]*f((b-a)*x[i]/2+(a+b)/2);
		cout << x[i] << " -- " << w[i] << "\n";
	}
	
	I = (b-a)*I/2;
	
	cout << "GAUSSIAN QUADRATURE:\n";
	cout << "Integral of f(x) from " << a << " to " << b << " = " << I;
   	
	return 0;
}
