// Method of steepest descent for finding function minima

// Date: 28/2/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

double f(double x, double y){			// function f(x,y) whose minimum is to be found
	return pow(pow(x,2)-2*y,2)+pow(x-y,2)+x+5;
}

void grad(double *gradient, double x, double y){
	gradient[0]=4*pow(x,3)-8*x*y+2*x-2*y+1;		// x component of gradient
	gradient[1]=-4*pow(x,2)+10*y-2*x;			// y component of gradient
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
	std::cout << "f(x) = (x^2-2y)^2+(x-y)^2+x+5\nMinimum is at:\n\n";
	cout << setprecision(9);
	
	double x=-1, y=-1;		// starting point
	
	double h=0.2, k;		// step size for search in x and y direction
	double TOL = 1e-7;		// tolerance of search
	
	double gradient[2];
	
	double diff[2];		// difference between two successive estimations
	do {
		grad(gradient, x, y);
		if(gradient[0]>0){ h=-fabs(h); }
		else{ h=fabs(h); }
		k=gradient[1]*h/gradient[0];
		// cout<<h<<" "<<k<<endl;
		if((f(x+h,y+k)-f(x,y))<0){
			x=x+h;   y=y+k;
			// cout << x<<" "<<y<<" "<<f(x,y)<<endl;
		}
		else{
			x=x+h;   y=y+k;		h=h/2;
			// cout << x<<" "<<y<<" "<<f(x,y)<<endl;
		}
		diff[0]=h;		diff[1]=k;
	} while (norm(diff,0) > TOL);
	
	// Print the solution:
	std::cout <<endl<< x+h/2 <<" "<<y+k/2;
	
	return 0;
}
