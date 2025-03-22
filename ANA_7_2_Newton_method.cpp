// NEWTON'S METHOD FOR FINDING VARIOUS VARIABLES FUNCTIONS MINIMA

// Date: 28/2/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

double f(double x, double y){	// function f(x,y) whose minimum is to be found
	return pow(x,2)+2*pow(y,2)+x*y+3*x;
}

void grad(double *gradient, double x, double y){
	gradient[0]=2*x+y+3;		// x component of gradient
	gradient[1]=4*y+x;			// y component of gradient
}

void hess_1(double hessian_1[2][2], double x, double y){	// Hessian matrix
	hessian_1[0][0]=4.0/7.0;
	hessian_1[0][1]=-1.0/7.0;
	hessian_1[1][0]=-1.0/7.0;
	hessian_1[1][1]=2.0/7.0;
}

void* prod(double hessian_1[2][2], double *gradient, double *diff){	// Hessian matrix
	diff[0]=-hessian_1[0][0]*gradient[0]-hessian_1[0][1]*gradient[1];
	diff[1]=-hessian_1[1][0]*gradient[0]-hessian_1[1][1]*gradient[1];
}

double norm(double v[2], int p){				// norm of vector v[n]
	double norm = 0;
	if (p==0){			// Supremum norm
		for (int j=0 ; j<=1 ; j++){
			if (fabs(v[j])>norm){
				norm = fabs(v[j]);
			}
		}
	}
	else if (p==1){
		for (int j=0 ; j<=1 ; j++){
			norm += fabs(v[j]);
		}
	}
	else if (p>=2){
		for (int j=0 ; j<=1 ; j++){
			norm += pow(v[j],p);}
		norm = pow(norm,pow(p,-1));	
	}
	return norm;	
}

int main(){
	// Function to find its minimum, z=f(x,y)
	std::cout << "f(x) = x^2+2y^2+xy+3x\nMinimum is at:\n\n";
	cout << setprecision(9);
	
	double x=0, y=0;		// starting point
	double diff[2], h, k;	// difference between two successive estimations
	double TOL=1e-3;		// tolerance of search
	
	double gradient[2];
	double hessian_1[2][2];
	grad(gradient,x,y);
	hess_1(hessian_1,x,y);			// Inverse of Hessian matrix
	
	do {
		grad(gradient, x, y);
		hess_1(hessian_1, x, y);
		prod(hessian_1,gradient,diff);		// difference vector
		h=diff[0];      k=diff[1];
		// cout << h<<" "<<k<<" "<<endl;
		x=x+h;   y=y+k;					// next iteration
		// cout << x<<" "<<y<<" "<<endl;
	} while (norm(diff,0) > TOL);
	
	// Print the solution:
	std::cout <<endl<<x+h/2<<" "<<y+k/2;
	
	return 0;
}
