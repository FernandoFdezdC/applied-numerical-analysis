// Calculate order n determinant with various methods

// Date: 27/2/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
#include<string>
using namespace std;
const int n = 4; 			// order of the square matrix
const int v = 4;

void print(int M[][n+1], int m, int p){			// Displays mxp matrix M
	for (int i=0 ; i<=m-1 ; i++){
		for (int j=0 ; j<=p-1 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

int fact(int m){				// Calculate factorial of number m
	int f = 1;			// Factorial
	for (int i=2 ; i<=m ; i++){
		f = f*i;
	}
	return f;			// Return factorial of m
}


int main()
{
	cout << setprecision(9);
	
	double M[n][n];				// Enter order n matrix M
	ifstream file;
    file.open("Matrix.txt", ios::in);
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			file >> M[i][j];
		}
	}
	file.close();
	// cout << "M = \n";
	// print(M,n);			// Display matrix M
	
	// BY DEFINITION:
	// Permutation matrix:
	int Perm[fact(n)][n+1];		// Matrix with every permutation of n together with its sign
	for (int j=0 ; j<=n-1 ; j++){ Perm[0][j]=j+1; }
	Perm[0][n]=1;
	for (int i=1 ; i<=fact(n)-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){ Perm[i][j]=Perm[i-1][j]; }
		Perm[i][n]=-Perm[i-1][n];  // Sign of the permutation
	}
	print(Perm,fact(n),n+1);
	cout << "\n";
	
	// Determinant:
	double det=0;
	double prod;
	for (int i=0 ; i<=fact(n)-1 ; i++){
		prod=1;
		for (int j=0 ; j<=n-1 ; j++){
			 prod=prod*M[Perm[i][j]-1][j];
		}
		prod = prod*Perm[i][n];		// Multiply by the sign of the permutation
		// cout << prod << endl;
		det = det + prod;
	}
	
	std::cout << "\nDeterminant of matrix M = " << det;
	
	// COFACTOR EXPANSION:
	
	int j=0;
	double det;
	for (int i=0 ; i<=n-1 ; i++){
		for (int k=0 ; k<=n-1 ; k++){
			double Mjk[n-1][n-1];
			for (int j=0 ; j<=n-1 ; j++){
				= ;
			}
			
			det = det + M[j][k]*pow(-1,j+k)*cofact(Mjk, n-1);
		}
	}
	return det;
	
	// Vandermonde matrix's determinant:
	double x[v] = {3.2,2.7,1.0,4.8};		// Vandermonde values
	double V=1;					// determinant
	for (int i=0 ; i<=v-2 ; i++){
		for(int j=i+1 ; j<=v-1 ; j++){
			V = V*(x[j]-x[i]);
		}
	}
	
	// Print the solution
	//cout << "Difference: " << x[n]-x[n-1] <<endl;
}
