// BEZIER CURVES

// Date: 23/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n=3; 			// degree of the Bezier polynomial
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
	
	// P(u) = \sum_{i=1}^{n}choose(n,i)(1-u)^{n-i} u^i p_i
	
	double x[n+1]={0.1,0.2,0.3,0.2},		// data to be interpolated
	y[n+1] = {0.0,0.2,0.1,0.0};		// data to be interpolated
	
	
	// CONSTRUCT BEZIER CURVE:
	double u[N+1]; 		// Parameter u
	for (int j=0 ; j<=N ; j++){
		u[j] = j*pow(N,-1);
		// cout << u[j] << endl;
	}
	// cout << fact(0)<<endl;
	// cout << fact(1)<<endl;
	double X[N+1], Y[N+1];
	for (int i=0 ; i<=N ; i++){
		X[i] = pow(1-u[i],n)*x[0];
		Y[i] = pow(1-u[i],n)*y[0];
		for (int j=1 ; j<=n ; j++){
			X[i] = X[i] + fact(n)*pow(fact(j)*fact(n-j),-1)*pow(1-u[i],n-j)*pow(u[i],j)*x[j];
			Y[i] = Y[i] + fact(n)*pow(fact(j)*fact(n-j),-1)*pow(1-u[i],n-j)*pow(u[i],j)*y[j];
		}
	}
	// Construction of Bezier curve. END.

	// Print Bezier curve data:
	cout<<"X = [";
	for (int i=0 ; i<=N-1 ; i++){
		std::cout << X[i];
		cout << ", ";
	}
	cout << X[N] << "];\nY = [";
	for (int i=0 ; i<=N-1 ; i++){
		std::cout << Y[i];
		cout << ", ";
	}
	cout << Y[N] << "];";
	
	return 0;
}
