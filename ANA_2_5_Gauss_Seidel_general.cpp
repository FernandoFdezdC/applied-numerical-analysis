// GAUSS-SEIDEL GENERALISED METHOD FOR SYSTEMS OF NON-LINEAR EQUATIONS

// Date: 5/11/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int n = 3;			// number of equations

double F(double v[n], int i){
	double F;
	double x = v[0];
	double y = v[1];
	double z = v[2];
	if (i == 1)
		F = (1-y*sin(z))/4;
	else if (i == 2)
		F = exp(-cos(z)/2-pow(x,0.5));
	else if (i == 3)
		F = (1+y*log(4*x)-exp(z))/2;
	return F;
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
   	// x_i = F(x_1, x_2, ..., x_n,i)	for i = 0, 1, ..., n-1
   	
   	double TOL = 1.0e-6;			// Tolerance
   	double c = 0;					// counter of iterations
   	double oldx[n], newx[n];		// vector with the iterated solutions
   	double diff[n]; 				// difference vector between two estimations
   	double v[n];					// auxiliary vector
   	
   	// Initialization:
   	oldx[0] = 1/3;
   	oldx[1] = 1/3;
   	oldx[2] = M_PI/6;
   	for (int i=0 ; i<=n-1 ; i++){
	   newx[i] = oldx[i];
	}
   	
   	do {
   		for (int i=0 ; i<=n-1 ; i++){
   			oldx[i] = newx[i];
		}
		
		for (int i=0 ; i<=n-1 ; i++){
	   		for (int j=0 ; j<=i-1 ; j++){
				v[j] = newx[j];
			}
			for (int j=i ; j<=n-1 ; j++){
				v[j] = oldx[j];
			}
			newx[i] = F(v,i+1);
		}
		
		c = c+1;		// counter adds an iteration
		for (int k=0 ; k<=n-1 ; k++){
			diff[k] = newx[k] - oldx[k];
		} // difference vector
	} while (norm(diff,0) > TOL);
	
   	// Display the solution:
   	cout << "Solution:\n";
   	for (int i=0 ; i<=n-1 ; i++){
   		cout << "x_" << i+1 << " = " << newx[i] << "\n";
	}
   	
	return 0;
}
