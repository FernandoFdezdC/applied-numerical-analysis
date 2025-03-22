// DIVIDED DIFFERENCES

// Date: 9/3/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n=5; 			// amount of data to be interpolated

double int_pol(double x0, double *x, double *a, int n){
	// Calculate at x0 n-th degree polynomial
	// that best fits the given data 
	double P = a[0];
	double fact=1.0;
	for (int j=1 ; j<=n ; j++){
   		fact = fact*(x0-x[j-1]);
   		P += a[j]*fact;
   		cout << P;
	}
	return P;
}

int main(){
	cout << setprecision(9);
	
	// Consider the n-th degree polynomial written as
	// P_n(x) = a_0+(x-x_0)a_1+(x-x_0)(x-x_1)a_2+\ldots+(x-x_0)(x-x_1)\ldots(x-x_{n-1})a_n
	
	double x[n]={3.2,2.7,1.0,4.8,5.6},		// data to be interpolated
	y[n] = {22.0,17.8,14.2,38.3,51.7};		// data to be interpolated
	
	double a[n]; 		// coefficients of the changed polynomial
	
	// Divided differences table:
	double f[n][n];		// matrix with divided differences
	for (int i=0 ; i<=n-1 ; i++){
		f[i][0]=y[i];		// first column of divided differences table
	}
	// Neville's method:
	for (int j=1 ; j<=n-1 ; j++){
		for (int i=0 ; i<=n-j-1 ; i++){
			f[i][j]=(f[i+1][j-1]-f[i][j-1])/(x[i+j]-x[i]);		// set the other columns equal to 0
		}
		for (int i=n-j ; i<=n-1 ; i++){
			f[i][j]=0;		// set the other elements equal to 0
		}
	}
	cout << "Divided differences table:\n";
	for (int i=0 ; i<=n-1 ; i++){
		cout <<x[i]<<" ";
		for (int j=0 ; j<=n-1 ; j++){
			std::cout << f[i][j] << " ";
		}
		std::cout << "\n";
	}
	// Divided differences table. END.
	
	for (int j=0 ; j<=n-1 ; j++){
   		a[j] = f[0][j];
	}
	
	// Display the solution:
	cout << "\nThe divided differences coefficients are:";
	for (int j=0 ; j<=n-1 ; j++){
   		cout << "\na_" << j << " = " << a[j];
	}
	
	// Interpolated value at x0=3.0:
	double x0=3.0;
	double val = int_pol(x0,x,a,3);
	cout <<"\n\nValue at "<<x0<<" of the interpolating polynomial: "<<val;
	
	return 0;
}
