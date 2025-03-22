// NEWTON-RAPHSON GENERALISED METHOD FOR SYSTEMS OF NON-LINEAR EQUATIONS

// Date: 5/11/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;
const int n = 2;			// number of equations


double* LU_method(double M[n][n], double b[n], double *X, int n){
	
	// n: number of equations
	// X[n]: solution of the system of linear equations M*X=b
	
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
		
	// Solve the system LUx=b:
	double Z[n];		// L*Z=b and U*X=Z
	
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
	return X;
}

double* increment(double v[n], double *arr,  int n){
	
	double x = v[0];
	double y = v[1];
	
	double J[n][2];				// Jacobian matrix
	
	// Evaluate in Jacobian matrix:
	J[0][0] = -0.5+y*cos(x*y)/2;
	J[0][1] = -1/(4*M_PI)+x*cos(x*y)/2;
	J[1][0] = -2*exp(1)+(2-1/(2*M_PI))*exp(2*x);
	J[1][1] = exp(1)/M_PI;
	
	// Functions evaluated in x_1, x_2, ..., x_n
	double f = sin(x*y)/2-y/(4*M_PI)-x/2;
	double g = (1-1/(4*M_PI))*(exp(2*x)-exp(1))+exp(1)*y/M_PI - 2*exp(1)*x;
	
	double b[n] = {-f,-g};
	
	double X[n];
	// Solve the system J*delta = b using LU method:
	arr = LU_method(J,b,X,n);
	return arr;
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
	
   	// System of non-linear equations:
	// f(x,y) = sin(x*y)/2-y/(4*M_PI)-x/2 = 0
	// g(x,y) = (1-1/(4*M_PI))*(exp(2*x)-exp(1))+exp(1)*y/M_PI - 2*exp(1)*x = 0
   	
   	double TOL = 1.0e-6;			// Required tolerance
   	double c = 0;					// counter of iterations
   	double oldx[n], newx[n];		// vector with the iterated solutions
   	double diff[n]; 				// difference vector between two estimations
   	double arr[n];					// Newton difference
   	
   	// Initialization:
   	oldx[0] = 0.6;
   	oldx[1] = 3;
   	for (int i=0 ; i<=n-1 ; i++){
	   newx[i] = oldx[i];
	}
   	
   	do {
   		for (int i=0 ; i<=n-1 ; i++){
   			oldx[i] = newx[i];
		}
   		
   		double* delta = increment(oldx,arr,n);

   		for (int i=0 ; i<=n-1 ; i++){		// i: index of the unknown variable
			newx[i] = oldx[i] + delta[i];
		}

		c = c+1;			// counter adds an iteration
		for (int k=0 ; k<=n-1 ; k++){
			diff[k] = delta[k];
		} // difference vector
	} while (norm(diff,0) > TOL);
   	
   	// Display the solution:
   	cout << "Solution:\n";
   	for (int i=0 ; i<=n-1 ; i++){
   		cout << "x_" << i+1 << " = " << newx[i] << "\n";
	}
   	
	return 0;
}
