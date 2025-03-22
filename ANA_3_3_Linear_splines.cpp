// LINEAR SPLINES

// Date: 9/3/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n=5; 			// amount of data to be interpolated

double int_pol(double x0, double *x, double a[2][n-1]){
	// Calculate at x0 the linear spline
	// that best fits the given data 
	double P;
	if (x0<x[0]){ P = a[0][0]+a[1][0]*x0; }
	else if (x0>x[n-1]){ P = a[0][n-2]+a[1][n-2]*x0; }
	else{
		for (int i=0 ; i<=n-2 ; i++){
	   		if ((x0-x[i])*(x0-x[i+1])<0){
	   				P = a[0][i]+a[1][i]*x0;
			   	}
		}
	}
	return P;		// returned value
}

int main(){
	cout << setprecision(9);
	
	double x[n]={1.10,2.00,3.50,5.00,7.10},		// data to be interpolated
	y[n] = {0.6981,1.4715,2.1287,2.0521,1.4480};		// data to be interpolated
	
	double a[2][n-1]; 		// matrix with y-intercept and slope of each spline

	// Construct linear spline:
	for (int i=0 ; i<=n-2 ; i++){
   		a[0][i] = -x[i]*(y[i+1]-y[i])/(x[i+1]-x[i])+y[i];	// y-intercept
   		a[1][i] = (y[i+1]-y[i])/(x[i+1]-x[i]);		// slope
	}
	
	// Interpolated values in interval [0,10]:
	double x0=0.0;	
	double val;
	cout<<"[";
	do{
	val = int_pol(x0,x,a); // returned value
	// cout <<"\n\nValue at "<<x0<<" of the interpolating polynomial: "<<val;
	cout <<val<<", ";
	x0=x0+0.01;
	}while (x0<=9.99);
//	cout<<endl<<x0<<endl;
	val = int_pol(x0,x,a); // returned value
	cout <<val<<"]";
	
	return 0;
}
