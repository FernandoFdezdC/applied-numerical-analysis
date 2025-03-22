// FINITE DIFFERENCE METHOD FOR 2ND ORDER BOUNDARY VALUE PROBLEMS

// DATE: 6/7/2021
// AUTHOR: Fernando Fernández del Cerro

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <math.h>
using namespace std;

const int n=4;				// Number of intervals

void Tridiag(double a[], double b[], double c[], double f[], double *u, int n){
	// Solve system of n-1 linear equations with tridiagonal system matrix.
	// u[n+1]: vector with solution
	
	double alpha[n+1], beta[n+1], gamma[n+1], z[n+1];
	
	for (int i=1 ; i<=n-2 ; i++) {
	    gamma[i] = c[i];
	}
	beta[1] = b[1];
	for (int i=2 ; i<=n-1 ; i++) {
	    alpha[i] = a[i]/beta[i-1];
		beta[i] = b[i] -alpha[i]*c[i-1];
	}
	z[1] = f[1];
	for (int i=2 ; i<=n-1 ; i++){
		z[i] = f[i] - alpha[i]*z[i-1];
	}
	// Find solution u:
	u[n-1] = z[n-1]/beta[n-1];
	for (int i=n-2 ; i>=1 ; i--){
		u[i]=(z[i]-c[i]*u[i+1])/beta[i];
	}
}

int main(){
	cout << setprecision(9);
	
	double x0=1.0, xf=3.0;			// Interval of integration
	double u_x0=2.0, u_xf=-1.0;		// initial values for dependent variable u
	double h=(xf-x0)/n;			// Integration step
	
	double* u = new double[n+1];		// Vector with solution
	double x[n+1];			// Vector with nodes
	for (int i=1 ; i<=n-1 ; i++) {
	    x[i]=x[0] + i*h;}
	
	double a[n+1], b[n+1], c[n+1], f[n+1];		// System coefficients
	
	x[0]=x0;
	u[0]=u_x0; 		// boundary condition
	x[n]=xf;
	u[n]=u_xf; 		// boundary condition
	
	// Fill vectors a, b, c & f:
	for (int i=1 ; i<=n-1 ; i++) {
	    a[i]= 1.0;
	    b[i]= -2.0-h*h*(1.0-x[i]/5.0);
	    c[i]= 1.0;
	    f[i]= h*h*x[i];
//		cout << a[i] << ", " << b[i] << ", "<<c[i] <<endl;
	}
	// Correct first and last equations:
	f[1] = f[1] -a[1]*u[0];
	f[n-1] = f[n-1] -c[n-1]*u[n];
	
	for (int i=1 ; i<=n-1 ; i++) {
		cout << f[i] <<endl;
	}
	
	Tridiag(a, b, c, f, u, n);		// Solve tridiagonal system of equations
	
	// DISPLAY SOLUTION:
	std::cout << "SOLUTION IS\nx_n=[";
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << x[i] << ", ";
	}
	cout << x[n] << "];\nu=[";
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << u[i] << ", ";
	}
	cout << u[n] << "];";
	
	// system("PAUSE");
	
	return 0;
}
