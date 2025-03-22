// CROUT FACTORISATION FOR ILL-CONDITIONED SYSTEMS

// Date: 5/11/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int n=3;			// number of equations in the system

void print(double M[n][n], int n){			// Displays matrix M
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

int main(){
	double M[n][n];				// Matrix M such that Mx=b
	ifstream file;
    file.open("Matrix.txt", ios::in);
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			file >> M[i][j];
		}
	}
	file.close();
	std::cout << "M=\n";
	print(M,n);			// Display matrix M
	
	double b[n] = {-1.61, 7.22, -3.38};		// vector b such that Mx=b
	
	// LU factorisation:
	double L[n][n];			// Lower matrix
	double U[n][n];			// Upper matrix
	// M = LU
	
	double S;
  	for (int k=0 ; k<=n-1 ; k++){
  		for (int j=0 ; j<=k-1 ; j++){
			L[j][k] = 0;
			U[k][j] = 0;
		}
	}			// Upper and lower matrices
	
  	for (int k=0 ; k<=n-1 ; k++){
  		
  		U[k][k] = 1;		// Crout method
  		
  		S = 0;
  		for (int s=0 ; s<=k-1 ; s++){
  			S += L[k][s]*U[s][k];
		}
  		L[k][k] = M[k][k] - S;
  		
  		for (int j=k+1 ; j<=n-1 ; j++){
  			S = 0;
  			for (int s=0 ; s<=k-1 ; s++){
  				S += L[k][s]*U[s][j];
			}
  			U[k][j] = (M[k][j] - S)/(L[k][k]);
		}
  		
  		for (int i=k+1 ; i<=n-1 ; i++){
  			S = 0;
  			for (int s=0 ; s<=k-1 ; s++){
  				S += L[i][s]*U[s][k];
			}
  			L[i][k] = M[i][k] - S;
		}
	}
	
	// Display matrices L and U:
	cout << "L = \n";
	print(L,n);
	cout << "U = \n";
	print(U,n);
	
	// Solve the system LUx=b:
	double Z[n];		// Lz=b;
	double X[n];		// Ux=z;
	
	for (int k=0 ; k<=n-1 ; k++){
		S = 0;
  		for (int s = 0; s <= k-1 ; s++){
  			S += L[k][s]*Z[s];
  		}
		Z[k] = (b[k] - S)/(L[k][k]);
	}
	
	for (int k=n-1 ; k>=0 ; k--){
		S = 0;
  		for (int s=k+1 ; s<=n-1 ; s++){
  			S += U[k][s]*X[s];
  		}
		X[k] = (Z[k] - S)/(U[k][k]);
	}
	

	// Display the solution to the equation Mx=b
	for (int j=0; j<=n-1 ; j++){
		std::cout << "x_" << j+1 << " = " << X[j] << endl;
	}
	
	return 0;
}
