// SHOOTING METHOD FOR 2ND ORDER BOUNDARY VALUE PROBLEMS

// DATE: 6/7/2021
// AUTHOR: Fernando Fernández del Cerro

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <math.h>
using namespace std;

const int m = 2; 				// Number of 1st order equations
const double h = 0.2;			// Independent variable step for integration

double dudx(double x, double y[m]){			// Equation for the derivative u'
	return y[1];
}
double d2udx2(double x, double y[m]){			// Equation for the second derivative u''
	return x+(1-x/5.0)*y[0]*y[1];
}

double (*Af[m]) (double x, double y[m])={dudx, d2udx2};

void r_k_y(double a, double u_a, double p0, int n, double *x, double *u, double *dudx){
	// RUNGE-KUTTA METHOD:
	double y[m];					// Vector with solution for each step
	// x[n+1]: independent variable vector
	// y[0]: dependent variable.
	// y[1]: derivative of dependent variable.
	// u[n+1], dudx[n+1]: vectors with solution
	y[0] = u_a, y[1]= p0;			// Initial conditions
	
	double k1[m], k2[m], k3[m], k4[m], yy[m];	// Auxiliar variables
	
	for (int i=0 ; i<=n ; i++){
		u[i] = y[0];
		dudx[i] = y[1];
		
		for (int j=0 ; j<=m-1  ; j++)    { k1[j]=h*Af[j](x[i],y);}
	    for (int j=0 ; j<=m-1 ; j++)    { yy[j]=y[j]+0.5*k1[j];}
	    for (int j=0 ; j<=m-1 ; j++)    { k2[j]=h*Af[j](x[i]+0.5*h,yy);}
	    for (int j=0 ; j<=m-1 ; j++)    { yy[j]=y[j]+0.5*k2[j];}
	    for (int j=0 ; j<=m-1 ; j++)    { k3[j]=h*Af[j](x[i]+0.5*h,yy); }
	    for (int j=0 ; j<=m-1 ; j++)    { yy[j]=y[j]+k3[j];}
	    for (int j=0 ; j<=m-1 ; j++)    { k4[j]=h*Af[j](x[i]+h,yy);}
	    
		for (int j=0 ; j<=m-1 ; j++)
    	{ y[j] = y[j] + (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j])/6.0;}
	}
}

int main(){
	cout << setprecision(9);
	
	double TOL = 1e-6; 							// REQUIRED TOLERANCE
	
	double a=1.0, b=3.0;			// Interval of integration
	double u_a=2.0, u_b=-1.0;		// initial values for dependent variable u
	int n = round((b-a)/h);			// Number of divisions
	double p0 = (u_b-u_a)/(b-a);	// first initial slope guess
	double p1;						// second initial slope guess
	
	double* u = new double[n+1];		// Vector with solution
	double* dudx = new double[n+1];		// Vector with solution
	
	double* x = new double[n+1];	// Vector with values of x at the end of each interval
	for (int i=0 ; i<=n ; i++){
		x[i] = a + i*h;
	}
	
	r_k_y(a, u_a, p0, n, x, u, dudx);		// APPLY RUNGE-KUTTA METHOD with initial slope p0
	double E0 = u[n] -u_b;			// Error at the end of the interval for p0
	
	if (p0 == 0.0){
		if (E0 > 0.0){
			p1 = -1.0;}		// second initial slope guess
		else if (E0 < 0.0){
			p1 = 1.0;}
	}
	else {
		if (E0*p0 > 0.0){
			p1 = p0/2.0;}
		else if (E0*p0 < 0.0){
			p1 = 2.0*p0;}
	}
	
	r_k_y(a, u_a, p1, n, x, u, dudx);		// APPLY RUNGE-KUTTA METHOD
	double E1 = u[n]-u_b;		// Error at the end of the interval for p1
	
	double p2 = p1 - E1*(p1-p0)/(E1-E0);	// third initial slope guess
	r_k_y(a, u_a, p2, n, x, u, dudx);		// APPLY RUNGE-KUTTA METHOD
	double E2 = u[n]-u_b;		// Error at the end of the interval for p2
	
	while (fabs(E2)>=TOL){
		p0 = p1;
		E0 = E1;
		p1 = p2;
		E1 = E2;
		p2 = p1 - E1*(p1-p0)/(E1-E0);
		r_k_y(a, u_a, p2, n, x, u, dudx);		// RUNGE-KUTTA METHOD
		E2 = u[n]-u_b;
		cout << "Error= " << E2 << "\n";
	}
	cout << "E0 = " << E0 << "\n" << "E1 = " << E1 << "\n" << "E2 = " << E2 << "\n";
	
	std::cout << "p2 = " << p2 << "\nu_p2 = " << u[n] << "\n";
	
	
	// DISPLAY SOLUTION:
	std::cout << "\nSOLUTION IS\nx_n=[";
	for (int i=0 ; i<n ; i++){
		std::cout << x[i] << ", ";
	}
	cout << x[n] << "];\nu(x)=[";
	for (int i=0 ; i<n ; i++){
		std::cout << u[i] << ", ";
	}
	cout << u[n] << "];\nu'=[";
	for (int i=0 ; i<n ; i++){
		std::cout << dudx[i] << ", ";
	}
	cout << dudx[n] << "];";
	
	// system("PAUSE");
	
	return 0;
}
