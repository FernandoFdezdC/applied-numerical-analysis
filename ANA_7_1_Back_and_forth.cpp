// Back and forth method for finding function minima

// Date: 27/2/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;


double f(double x){			// function f(x) whose minimum is to be found
	return exp(x) + 2 - cos(x);
}

int main(){
	// Function to find its minimum, f(x)
	std::cout << "f(x) = e^x + 2 - \cos{x}\nMinimum is at:\n\n";
	cout << setprecision(9);
	
	double a=-3, b=1;			// interval of minimum searching
	
	double x=a;
	double h=(b-a)/4;
	double TOL =1e-4;			// tolerance
	
	while (fabs(h) > TOL){
		while ((f(x+h)-f(x))<0){
			cout << x <<", "<<f(x)<<"\n";
			cout << x+h <<", "<<f(x+h)<<"\n";
			x=x+h;
		}
		cout << x <<", "<<f(x)<<"\n";
		cout << x+h <<", "<<f(x+h)<<"\n";
		x=x+h;
		h=-h/4;
		cout << x << endl;
		cout << h << endl;
	}
	
	// Print the solution:
	std::cout <<endl<< x+h/2;
	
	return 0;
}
