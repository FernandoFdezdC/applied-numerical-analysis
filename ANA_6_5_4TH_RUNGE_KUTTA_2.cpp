// SOLUTION OF SYSTEM OF ODE'S WITH 4TH ORDER RUNGE-KUTTA METHOD

// DATE: 20/12/2020
// AUTHOR: Fernando Fernández del Cerro

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <math.h>
using namespace std;

const int m = 4; 					// Number of 1st order equations

double v1(double t, double y[m]){				// Equations of velocity
	return y[2];
}
double v2(double t, double y[m]){
	double v2 = 0;
	return v2;
}
double a1(double t, double y[m]){				// Equations of acceleration
	double a1 = -2/t + y[2] - 2;
	return a1;
}
double a2(double t, double y[m]){
	double a2 = 0;
	return a2;
}

double (*Af[m]) (double t, double y[m])={v1, v2, a1, a2};

int main(){
	cout << setprecision(15);
	
	double h = 0.01;											// Interval of time for integration
	double t0 = 0.5, tf = 1.0;							// Interval of integration
	
	int n = round((tf-t0)/h);								// Number of divisions
	double y[m];											// Vector with solution for each step
	y[0] = 3, y[1] = 0, y[2]= 0, y[3] = 0;		// Initial conditions
	double t[n]; 											// Time vector
	double y1[n], y2[n], v1[n], v2[n];					// Solution vectors y1, y2, v1, v2
	
	cout << "4TH ORDER RUNGE-KUTTA METHOD:\n";
	
	int j;
	double k1[m], k2[m], k3[m], k4[m], yy[m];			// Auxiliar variables
	
	// Write the solution in a file
	ofstream file1, file2;
	if (!file1 || !file2){
		cout << "No se pueden abrir los archivos";
    	system("PAUSE");
	}
	file1.open ("P15_y1.txt");
	file2.open ("P15_y2.txt");
	file1 << "4TH ORDER RUNGE-KUTTA METHOD. y1(t) with h=" << h << ":\n";
	file2 << "4TH ORDER RUNGE-KUTTA METHOD. y2(t) with h=" << h << ":\n";
	file1 << "[";
	file2 << "[";

	// Runge-Kutta iterations:
	for (int i=0 ; i<n ; i++){
		t[i] = t0+i*h;
		file1 << y[0] << ", ";
		file2 << y[1] << ", ";
		y1[i] = y[0];
		y2[i] = y[1];
		v1[i] = y[2];
		v2[i] = y[3];
		
		for (j=0; j<m; j++)    { k1[j]=h*Af[j](t[i],y);}
	    for (j=0; j<m; j++)    { yy[j]=y[j]+0.5*k1[j];}
	    for (j=0; j<m; j++)    { k2[j]=h*Af[j](t[i]+0.5*h,yy);}
	    for (j=0; j<m; j++)    { yy[j]=y[j]+0.5*k2[j];}
	    for (j=0; j<m; j++)    { k3[j]=h*Af[j](t[i]+0.5*h,yy); }
	    for (j=0; j<m; j++)    { yy[j]=y[j]+k3[j];}
	    for (j=0; j<m; j++)    { k4[j]=h*Af[j](t[i]+h,yy);}
	    
		for (j=0; j<m; j++)
    	{ y[j]=y[j] + (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j])/6.0;}
	}
	
	
	file1 << y[0] << "]";
	file2 << y[1] << "]";
	file1.close();
	file2.close();
	
	// system("PAUSE");
	return 0;
}
