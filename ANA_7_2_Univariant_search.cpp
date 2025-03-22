// Univariant method for finding function minima

// Date: 27/2/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

double f(double x, double y){			// function f(x,y) whose minimum is to be found
	return pow(pow(x,2)-2*y,2)+pow(x-y,2)+x+5;
}

int main(){
	// Function to find its minimum, z=f(x,y)
	std::cout << "f(x) = (x^2-2y)^2+(x-y)^2+x+5\nMinimum is at:\n\n";
	cout << setprecision(9);
	
	double x=-1, y=-1;			// starting point
	
	double h=0.1;			// step size for search
	double TOL = 1e-7;			// tolerance
	
	while (fabs(h) > TOL){
		while ((f(x+h,y)-f(x,y))<0){
			cout << x<<" "<<y<<" "<<f(x,y)<<endl;
			x=x+h;
		}
		x=x+h;
		while ((f(x,y+h)-f(x,y))<0){
			cout << x<<" "<<y<<" "<<f(x,y)<<endl;
			y=y+h;
		}
		y=y+h;
		h=-h/2;
		// cout<<h<<endl;
	}
	
	// Print the solution:
	std::cout <<endl<< x+h/2 <<" "<<y+h/2;
	
	return 0;
}
