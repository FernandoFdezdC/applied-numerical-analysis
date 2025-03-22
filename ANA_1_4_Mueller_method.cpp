// MUELLER'S METHOD FOR NON-LINEAR EQUATIONS

// Date: 11/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;

double f(double x){			// function f(x) whose roots are to be found
	return 3*x + sin(x)-exp(x);
}

double dfdx(double x){		// derivative of function f(x)
	return 3+cos(x)-exp(x);
}

int main(){
	cout << setprecision(9);
	
	// Equation to be solved: f(x) = 0
	std::cout << "f(x) = 3x + sin(x) - e^x" << endl;

	int n = 1;			// counter of iterations
	
	// first tries:
	double x2=-2.0,   x0=0.0,	x1=1.0;
	double TOL=1e-8;		// required tolerance
	
	
	double f2, f0, f1;
	double h1, h2, gamma, a, b, c, root;
	do {
		f2=f(x2),	f0=f(x0),	f1=f(x1);
		h1 = x1-x0; h2=x0-x2; gamma=h2/h1;
		// Coefficients of the parabola:
		c=f0;
		a=(gamma*f1-f0*(1+gamma)+f2)/(gamma*h1*h1*(1+gamma));
		b=(f1-f0-a*h1*h1)/h1;
		// Roots of the polynomial:
		if (b<0){
			root = x0-2*c/(b-pow(b*b-4*a*c,1.0/2.0));
		}
		else if (b==0){
			root = x0+pow(-a/c,1.0/2.0);
		}
		else{
			root = x0-2*c/(b+pow(b*b-4*a*c,1.0/2.0));
		}
		if (root>x0){
			x2=x0;
			x0=x1;
			x1=root;
		}
		else{
			x1=x0;
			x0=root;
		}
		cout<<root<<"\n";
		n++;
	} while(fabs(f(root)) > TOL);
	
	// Print the solution
	std::cout << root << " is a solution to the equation f(x) = 0.";
	std::cout << "\nNecessary iterations: " << n;
	
	return 0;
}
