// NEWTON'S METHOD APPLIED TO FIND SQUARE ROOT OF NUMBER N

// Date: 11/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;


int main()
{
	cout << setprecision(9);
	
	// Find the square root of N
	double N = 101.0;
	std::cout << "x^2=" << N << endl;
	
	int n=1;			// counter of iterations
	
	double x0=N,   x1;	// first try, x0;
	double TOL=1e-7;		// required tolerance
	
	if (x0 == 0){
		std::cout << "ERROR: Newton's method cannot be applied if x_0=0.";
		exit(0);		// Ends the program
	}
	
	double diff;	// difference between successive values of x
	do {
		x1 = (x0 + N/x0)/2;
		diff = x1-x0;
		x0=x1;
		n++;
	} while(fabs(diff)>TOL || fabs(pow(x1,2)-N)>TOL);
	
	// Print the solution
	std::cout << x1 << " is the square root of "<<N;
	std::cout << ".\nNecessary iterations: " << n;
	std::cout << "\nNote: The method may converge to a root different from the expected one or diverge if the starting value is not close enough to the root.";
	
	return 0;
}
