// CALCULATE CHEBYSHEV POLYNOMIALS COEFFICIENTS

// Date: 1/3/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int n = 10;			// Degree of the last polynomial calculated

void print(double M[n+1][n+1], int n){			// Displays nxn matrix M
	for (int i=0 ; i<=n ; i++){
		for (int j=0 ; j<=n ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

int main(){
	cout << setprecision(9);
	// T_0(x)=1, 		T_1(x)=x
	// T_n(x) = a[n][0]+a[n][1]*x+...+a[n][n]*x^n
	// T_{n+1}(x) = 2xT_n(x)-T_{n-1}(x)
	// PROPERTY: T_n(\cos{\theta})=\cos{n\theta}
	
	double a[n+1][n+1];		// Coefficients of the Legendre polynomials
							// First entry: degree of polynomial.
							// Second entry: coefficients of the polynomial
	a[0][0] = 1.0;
	a[1][0] = 0.0;
	a[1][1] = 1.0;
	
	// Recurrence relation:
	// For all i>=2 we have:
	// a_{i,0}=-a_{i-2,0}
	// a_{i,j}=2a_{i-1,j-1}-a_{i-2,j}, with j>=1
	for(int i=0 ; i<=1 ; i++){
		for(int j=i+1 ; j<=n+1 ; j++){
			a[i][j] = 0;
		}
	}
	for(int i=2 ; i<=n ; i++){
		a[i][0] = -a[i-2][0];
		for(int j=1 ; j<=n ; j++){
			a[i][j] = 2*a[i-1][j-1]-a[i-2][j];
		}
	}
	cout << "CHEBYSHEV COEFFICIENTS:\n";
	print(a,n);
	
	return 0;
}
