// LU DOOLITTLE FACTORISATION FOR TRIGIADONAL SYSTEMS OF EQUATIONS

// Date: 5/11/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

void Tridiag(double a[], double b[], double c[], double f[], double x[], int n){
	
	// x: solution vector of the tridiagonal system
	// n: number of equations
	
	double alpha[n], beta[n], z[n];
	beta[0] = b[0];
	alpha[0] = 0;
	z[0] = f[0];
	
	for (int i = 1 ; i <= n-1 ; i++){
		alpha[i] = a[i]/beta[i-1];
		beta[i] = b[i] - alpha[i]*c[i-1];
		z[i] = f[i] - alpha[i]*z[i-1];
	}
	
	x[n-1] = z[n-1]/beta[n-1];
	
	for (int j = n-2 ; j >= 0 ; j--)
		x[j] = (z[j] - c[j]*x[j+1])/beta[j];
	
}

int main(){
	cout << setprecision(9);
	
	int n = 1000; 		// Number of equations of the tridiagonal system
	
	double a[n], b[n], c[n], f[n], x[n];
	
	// Define a_n, b_n and c_n:
	for (int i = 0 ; i <= n-1 ; i++){
		a[i] = -1.0;
		b[i] = 4.0;
		c[i] = -1.0;
	}
	
	// Define f:
	f[0] = 100.0;
	f[n-1] = 100.0;
	for (int i = 1 ; i <= n-2 ; i++)
		f[i] = 200.0;
	
	// Solve the tridiagonal system:
	Tridiag(a, b, c, f, x, n);
	
	// Test the solution:
	std::cout << b[0]*x[0] + c[0]*x[1] << "\n";
	for (int i = 1; i <= n-2; i++)
		std::cout << a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] << "\n";
	std::cout << a[n-1]*x[n-2] + b[n-1]*x[n-1] << "\n\n";
	
	// Print the solution x of the tridiagonal system:
	cout << "SOLUTION IS:\n";
	for (int i=0 ; i<=n-1 ; i++)
		std::cout << x[i] << "\n";
	
	return 0;
}
