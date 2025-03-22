// SOLUTION OF SYSTEM OF ODE'S WITH 4TH ORDER RUNGE-KUTTA METHOD

// DATE: 9/12/2020
// AUTHOR: Fernando Fernández del Cerro

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <math.h>
using namespace std;

const double m1 = 2.0, m2 = 3.5, k1 = 2.5, k2 = 3.5;	// Problem parameters
const int m=4; 					// Number of 1st order equations

double v1(double t, double y[m]){			// Equations of velocity
	return y[2];
}
double v2(double t, double y[m]){			// Equations of velocity
	return y[3];
}
double a1(double t, double y[m]){				// Equations of acceleration
	return -k1*y[0]/m1-k2*(y[0]-y[1])/m1;
}
double a2(double t, double y[m]){
	return k2*(y[0]-y[1])/m2;
}

double (*Af[m]) (double t, double y[m])={v1, v2, a1, a2};

int main(){
	cout << setprecision(15);
	
	double h = 1.0;								// Interval of time for integration
	double t0 = 0.0, tf = 100.0;				// Interval of integration
	
	int n = round((tf-t0)/h);							// Number of divisions
	double y[m];										// Vector with solution for each step
	y[0] = 3.0, y[1] = 4.0, y[2]= 0.0, y[3] = 0.0;		// Initial conditions
	double t[n+1]; 										// Time vector
	for (int i=0 ; i<=n-1 ; i++){
		t[i] = t0 + i*h;}
	double y1[n+1], y2[n+1], v1[n+1], v2[n+1];			// Solution vectors y1, y2, v1, v2
	
	
	cout << "4TH ORDER RUNGE-KUTTA METHOD:\n";
	
	double k1[m], k2[m], k3[m], k4[m], yy[m];		// Auxiliar variables
	
	
	// RUNGE-KUTTA METHOD:
	for (int i=0 ; i<=n ; i++){
		y1[i] = y[0];
		y2[i] = y[1];
		v1[i] = y[2];
		v2[i] = y[3];
		
		for (int j=0 ; j<=m-1 ; j++)    { k1[j]=h*Af[j](t[i],y);}
	    for (int j=0 ; j<=m-1 ; j++)    { yy[j]=y[j]+0.5*k1[j];}
	    for (int j=0 ; j<=m-1 ; j++)    { k2[j]=h*Af[j](t[i]+0.5*h,yy);}
	    for (int j=0 ; j<=m-1 ; j++)    { yy[j]=y[j]+0.5*k2[j];}
	    for (int j=0 ; j<=m-1 ; j++)    { k3[j]=h*Af[j](t[i]+0.5*h,yy); }
	    for (int j=0 ; j<=m-1 ; j++)    { yy[j]=y[j]+k3[j];}
	    for (int j=0 ; j<=m-1 ; j++)    { k4[j]=h*Af[j](t[i]+h,yy);}
	    
		for (int j=0 ; j<=m-1 ; j++)
    	{ y[j] = y[j] + (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j])/6.0;}
	}
	
	cout << "EULER PREDICTOR-CORRECTOR METHOD:\n";
	
	std::cout << "\nt=[";
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << t[i] << ", ";
	}
	cout << t[n] << "];\ny1=[";
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << y1[i] << ", ";
	}
	cout << y1[n] << "];\ny2=[";
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << y2[i] << ", ";
	}
	cout << y2[n] << "];\nv1=[";
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << v1[i] << ", ";
	}
	cout << v1[n] << "];\nv2=[";
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << v2[i] << ", ";
	}
	cout << v2[n] << "];";
	
	// system("PAUSE");
	
	return 0;
}
