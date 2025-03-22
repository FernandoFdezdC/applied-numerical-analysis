// Quadratic back and forth method for finding function minima

// Date: 27/2/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;


double func(double x){			// function f(x) whose minimum is to be found
	return exp(x) + 2 - cos(x);
}

void lagrange(double *a, double *x, double *f){		// Lagrange interpolating polynomial of degree 2
	// P_2(x) = ax^2+bx+c
	a[2]=(x[1]*f[0]+x[0]*f[2]+x[2]*f[1]-x[1]*f[2]-x[2]*f[0]-x[0]*f[1])/((x[0]-x[1])*(x[0]-x[2])*(x[1]-x[2]));		// a
	a[1]=(pow(x[0],2)*f[1]+pow(x[2],2)*f[0]+pow(x[1],2)*f[2]-pow(x[2],2)*f[1]-pow(x[0],2)*f[2]-pow(x[1],2)*f[0])/((x[0]-x[1])*(x[0]-x[2])*(x[1]-x[2]));		// b
	a[0]=(pow(x[0],2)*x[1]*f[2]+x[0]*pow(x[2],2)*f[1]+pow(x[1],2)*x[2]*f[0]-x[1]*pow(x[2],0)*f[0]-x[2]*pow(x[0],2)*f[1]-x[0]*pow(x[1],2)*f[2])/((x[0]-x[1])*(x[0]-x[2])*(x[1]-x[2]));		// c
}



int main(){
	// Function to minimize, f(x)
	std::cout << "f(x) = e^x + 2 - \cos{x}\nMinimum is at:\n";
	cout << setprecision(9);
	
	double x0=-2, x1=-1, x2=0;	// First approximations
	double TOL=1e-4;			// tolerance
	double x_min;
	
	while (fabs(x2-x0)>TOL){
		double x[3]={x0,x1,x2};
		double f[3]={func(x0),func(x1),func(x2)};
		double a[3];		// vector with Lagrange interpolating polynomial coefficients
		lagrange(a, x, f);		// Lagrange polynomial
		
		x_min = -a[1]/(2*a[2]);	// minimum of the interpolating polynomial
		cout<<x0<<" "<<x1<<" "<<x2<<" "<<x_min<<endl;
		
		if (x_min>=x2){ x0=x1 ; x1=x2 ; x2=x_min; }
		else if (x_min>=x1 && x_min<x2){
			if((x2-x_min)>(x_min-x0)){ x2=x_min; }
			else{ x0=x1; x1=x_min; }
		}
		else if (x_min>=x0 && x_min<x1){
			if((x2-x_min)>(x_min-x0)){ x2=x1; x1=x_min; }
			else{ x0=x_min; }
		}
		else if (x_min<x0){ x2=x1 ; x1=x0 ; x0=x_min; }
	}
	
	// Print the solution:
	std::cout << x_min;
	
	return 0;
}
