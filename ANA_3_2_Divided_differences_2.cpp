// DIVIDED DIFFERENCES

// Date: 9/3/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n=5; 			// amount of data to be interpolated

void* int_pol(double x0, double *x, double *a, double *ret, int n, int N){
	// Calculate at x0 the N-th degree polynomial
	// that best fits the given data 
	double P = a[0];
	double fact=1.0;
	for (int j=1 ; j<=N ; j++){
   		fact = fact*(x0-x[j-1]);
   		P += a[j]*fact;
	}
	double Err;
	if (N<n){
		fact = fact*(x0-x[N]);
		Err = a[N+1]*fact;} 	// next-term rule estimated error
	else { Err = 0.0;
	}
	ret[0]=P;		// Value of the polynomial together
	ret[1]=Err;		// with its estimated error
}

int main(){
	cout << setprecision(9);
	
	// Consider the n-th degree polynomial written as
	// P_n(x) = a_0+(x-x_0)a_1+(x-x_0)(x-x_1)a_2+\ldots+(x-x_0)(x-x_1)\ldots(x-x_{n-1})a_n
	
	double x[n]={1.10,2.00,3.50,5.00,7.10},		// data to be interpolated
	y[n] = {0.6981,1.4715,2.1287,2.0521,1.4480};		// data to be interpolated
	
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
	
	// Interpolated value at x0=1.75:
	double x0=1.75;
	double ret[2];		// returned vector with estimated value and error
	int_pol(x0,x,a,ret,n,3);
	cout <<"\n\nValue at "<<x0<<" of the interpolating polynomial: "<<ret[0];
	cout <<"\nEstimated error: "<<ret[1];
	
	return 0;
}
