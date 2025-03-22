// SOLVING A SYSTEM BY ITERATION

// Date: 5/11/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;
const int n = 2;			// number of equations


double norm(double v[n], int nor){		// norm of vector v[n]
	double norm = 0;
	if (nor==0){			// Supremum norm
		for (int j=0 ; j<=n-1 ; j++){
			if (fabs(v[j])>norm){
				norm = fabs(v[j]);}}}
	else if (nor==1){
		for (int j=0 ; j<=n-1 ; j++){
			norm += fabs(v[j]);}}
	else if (nor>=2){
		for (int j=0 ; j<=n-1 ; j++){
			norm += pow(v[j],nor);}
		norm = pow(norm,pow(nor,-1.0));	}
	return norm;	
}



int main(){
	cout << setprecision(9);
	
   	// System of non-linear equations:
	// f(x,y) = e^x - y = 0
	// g(x,y) = xy - e^x = 0
   	
   	double TOL = 1.0e-8;		// Required tolerance
   	double c = 0;				// counter of iterations
   	double x[n];		// vector with the iterated solutions
   	double diff[n]; 		// difference vector between two estimations
   	double arr[n];				// Newton difference
   	
   	// Initialization:
   	x[1] = 2.0;			// Start with y=2;
   	x[0] = log(x[1]);
   	
   	do {
   		cout<<x[0]<<", "<<x[1]<<"\n";
   		diff[1] = fabs(exp(x[0])/x[0]-x[1]);
   		x[1] = exp(x[0])/x[0];
   		diff[0] = fabs(x[0]-log(x[1]));
		x[0] = log(x[1]);
		c++;			// counter adds an iteration
	} while (norm(diff,0) > TOL);
   	
   	// Print the solution
   	cout << "Solution:\n";
   	for (int i=0 ; i<=n-1 ; i++){
   		cout << "x_" << i+1 << " = " << x[i] << "\n";
	}
   	std::cout << "\nNecessary iterations: " << c;
   	
	return 0;
}
