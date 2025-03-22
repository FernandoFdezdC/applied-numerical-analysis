// B-SPLINE CURVES

// Date: 23/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n=3; 			// number of points-1
const int N = 100; 		// resolution of the parameter u

int fact(int x){		// Calculates factorial of number x
	int fact = 1;
	for (int j=2 ; j<=x ; j++){
		fact = fact*j;
	}
	return fact;
}

int main(){
	cout << setprecision(9);
	
	// P(u) = \sum_{k=-1}^{2}b_k p_{i+k}
	
	double x[n+1]={0.0,2.0,4.5,5.2},		// data to be interpolated
	y[n+1] = {0.0,1.0,2.0,0.4};		// data to be interpolated
	
	
	// CONSTRUCT BEZIER CURVE:
	double u[N+1]; 		// Parameter u
	for (int j=0 ; j<=N ; j++){
		u[j] = j*pow(N,-1);
		// cout << u[j] << endl;
	}
	// cout << fact(0)<<endl;
	// cout << fact(1)<<endl;
	double X[n-1][N+1], Y[n-1][N+1];
	for (int i=1 ; i<=n-2 ; i++){
		for (int j=0 ; j<=N ; j++){
			X[i][j] = pow(1-u[j],3.0)*x[i-1]/6.0  + (pow(u[j],3.0)/2.0 - u[j]*u[j] + 2.0/3.0)*x[i]  + (-pow(u[j],3.0)/2.0 + u[j]*u[j]/2.0 + u[j]/2.0 + 1.0/6.0)*x[i+1]  + pow(u[j],3.0)*x[i+2]/6.0 ;
			Y[i][j] = pow(1-u[j],3.0)*y[i-1]/6.0  + (pow(u[j],3.0)/2.0 - u[j]*u[j] + 2.0/3.0)*y[i]  + (-pow(u[j],3.0)/2.0 + u[j]*u[j]/2.0 + u[j]/2.0 + 1.0/6.0)*y[i+1]  + pow(u[j],3.0)*y[i+2]/6.0 ;
		}
	}
	// cout << u[N]<<endl;
	// Construction of Bezier curve. END.

	// Print Bezier curve data:
	std::cout << "X = [";
	for (int i=1 ; i<=n-2 ; i++){
		for (int j=0 ; j<=N-1 ; j++){
			std::cout << X[i][j] << ", ";
		}
	}
	std::cout << X[n-2][N] << "];\nY = [";
	for (int i=1 ; i<=n-2 ; i++){
		for (int j=0 ; j<=N-1 ; j++){
			std::cout << Y[i][j] << ", ";
		}
	}
	std::cout << Y[n-2][N] << "];";
	
	return 0;
}
