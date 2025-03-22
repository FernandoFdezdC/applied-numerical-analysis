// GAUSS-SEIDEL METHOD FOR SYSTEMS OF LINEAR EQUATIONS

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

double norm(double v[n], int nor){				// norm of vector v[n]
	double norm = 0;
	if (nor==0){			// Supremum norm
		for (int j=0 ; j<=n-1 ; j++){
			if (fabs(v[j])>norm){
				norm = fabs(v[j]);
			}
		}
	}
	else if (nor==1){
		for (int j=0 ; j<=n-1 ; j++){
			norm += fabs(v[j]);
		}
	}
	else if (nor>=2){
		for (int j=0 ; j<=n-1 ; j++){
			norm += pow(v[j],nor);}
		norm = pow(norm,pow(nor,-1));	
	}
	return norm;	
}

int main(){
	cout << setprecision(9);
	
	double M[n][n];				// Matrix M such that Mx=b
	ifstream file;
    file.open("Matrix.txt", ios::in);
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			file >> M[i][j];
		}
	}
	file.close();
	cout << "M = \n";
	print(M,n);			// Display matrix M
	
   	double b[n] = {-1, 2, 3};		// vector b such that Mx=b
   	
   	// Check if M is strictly diagonally dominant
   	for (int i=0 ; i<=n-1 ; i++){
   		double S=0;
   		for (int j=0 ; j<=n-1 ; j++){
   			if (i != j){
   				S += abs(M[i][j]);}
		}
		if (abs(M[i][i]) <= S){
			std::cout << "WARNING: Coefficient matrix M is not strictly diagonally dominant" << endl;
			break;
		}
	}
	
	double TOL = 1.0e-6;			// Tolerance
	double c = 0;					// counter of iterations
   	float oldx[n], newx[n];			// vector with the iterated solutions
	double diff[n]; 				// difference vector between two estimations
   	for (int i=0 ; i<=n-1 ; i++){
	   oldx[i] = 0;
	   newx[i] = oldx[i];
	} // canonical initial try
   	
   	do {
   		for (int i=0 ; i<=n-1 ; i++){
   			oldx[i] = newx[i];
		}
   		for (int i=0 ; i<=n-1 ; i++){
   			newx[i] = b[i]/(M[i][i]);
   			for (int j=0 ; j<=n-1 ; j++){
			   	if (j != i){
				   newx[i]=newx[i]-(M[i][j]/M[i][i])*newx[j];
				}
			}
		}
		
		c = c+1;			// counter adds a iteration
		for (int k=0 ; k<=n-1 ; k++){
			diff[k] = newx[k] - oldx[k];
		} // difference vector
	} while (norm(diff,0) > TOL);
	
	// Display matrix M
	cout << "M = \n";
	print(M,n);
	
	// Display the solution to the equation Mx=b
	cout << "Solution:\n";
   	for (int j=0 ; j<=n-1 ; j++){
   		cout << "x_" << j+1 << " = " << newx[j] << "\n";
	}
	
	return 0;
}
