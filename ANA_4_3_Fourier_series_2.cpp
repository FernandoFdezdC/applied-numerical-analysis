// FOURIER SERIES

// Date: 24/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int N = 200;			// Resolution of the interval of approximation
const int M = 32;			// Last harmonic taken to approximate


double f(double x){
	return exp(x);
}

int main(){
	cout << setprecision(9);
	
	// CONSTRUCTION OF THE FOURIER SERIES:
	double x[N+1]; 		// interval of representation
	double f[N+1];		// Fourier series
	for (int i=0 ; i<=N ; i++){
		x[i] = 2*M_PI*i*pow(N,-1)-M_PI;
		f[i] = M_PI/2;
		for (int n=1 ; n<=M ; n++){
			f[i] = f[i] - 4*cos((2*n-1)*x[i])/(M_PI*pow(2*n-1,2));
		}
	}
	
	// PRINT THE FOURIER APPROXIMATION:
	std::cout << "\nx=[";
	for (int i=0 ; i<=N-1 ; i++){
		cout << x[i] << ", ";
	}
	cout << x[N] << "];\n";
	std::cout << "\nf(x)=[";
	for (int i=0 ; i<=N-1 ; i++){
		std::cout << f[i] << ", ";
	}
	cout << f[N] << "];";
	
	return 0;
}
