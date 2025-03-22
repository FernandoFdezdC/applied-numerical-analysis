// NEWTON'S METHOD FOR NON-LINEAR EQUATIONS
// SEARCH FOR COMPLEX ROOTS

// Date: 28/2/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
// #include <complex>
using namespace std;

void* div(double *division, double *z1, double *z2){			// division of two complex numbers
	// z_1/z_2 = (a+ib)/(c+id)=(ac+bd)/(c^2+d^2)+i(bc-ad)/(c^2+d^2)
	division[0]=(z1[0]*z2[0]+z1[1]*z2[1])/(z2[0]*z2[0]+z2[1]*z2[1]);		// Real part of the division
	division[1]=(z1[1]*z2[0]-z1[0]*z2[1])/(z2[0]*z2[0]+z2[1]*z2[1]);		// Imaginary part of the division
}

void* prod(double *product, double *z1, double *z2){			// product of two complex numbers
	// z_1 z_2 = (a+ib)*(c+id)=(ac-bd)+i(bc+ad)
	product[0]=z1[0]*z2[0]-z1[1]*z2[1];		// Real part of the product
	product[1]=z1[1]*z2[0]+z1[0]*z2[1];		// Imaginary part of the product
}

void* add(double *sum, double *z1, double *z2){			// addition of two complex numbers
	// z_1+z_2 = (a+ib)+(c+id)=(a+c)+i(b+d)
	sum[0]=z1[0]+z2[0];		// Real part of the addition
	sum[1]=z1[1]+z2[1];		// Imaginary part of the addition
}

void* function(double *f, double* z0){			// function f(x) whose roots are to be found
	double z2[2];		prod(z2,z0,z0);
	double z3[2];		prod(z3,z2,z0);
	double two[2]={2.0,0.0};
	double twoz2[2];	prod(twoz2,z2,two);
	double _one[2]={-1.0,0.0};
	double _z[2];		prod(_z,z0,_one);
	double five[2]={5.0,0.0};
	double sum1[2];		add(sum1,_z,five);
	double sum2[2];		add(sum2,sum1,twoz2);
	double sum3[2];		add(sum3,sum2,z3);
	f[0]=sum3[0]; 		f[1]=sum3[1];
}

void* dfdz(double *f, double *z0){		// derivative of function f(x)
	double z2[2];		prod(z2,z0,z0);
	double three[2]={3.0,0.0};
	double four[2]={4.0,0.0};
	double threez2[2];	prod(threez2,three,z2);
	double fourz[2];	prod(fourz,four,z0);
	double _one[2]={-1.0,0.0};
	double sum1[2];		add(sum1,fourz,_one);
	double sum2[2];		add(sum2,sum1,threez2);
	f[0]=sum2[0]; 		f[1]=sum2[1];
}

double norm(double z[2]){		// modulus of complex number z=x+iy
	double norm = pow(pow(z[0],2)+pow(z[1],2),1.0/2.0);
	return norm;
}

int main(){
	cout << setprecision(9);
	
	// Equation to be solved: f(x) = 0
	std::cout << "f(z) = z^3+2z^2-z+5" << endl;
	
	int n = 1;			// counter of iterations
	
	// z=x+iy
	double z0[2]={1.0,-1.0};	// first try, z0
	double z1[2];				// next approximation, z1			
	double TOL=1e-5;			// required tolerance
	
	double diff[2];		// difference between successive values of z
	double fz0[2], dfdz0[2];
	double division[2]; 		// division between f(z0) and dfdz(z0)
	
	do {
		function(fz0,z0);
		dfdz(dfdz0,z0);
		div(division,fz0,dfdz0);
		z1[0] = z0[0]-division[0];	z1[1] = z0[1]-division[1];
		diff[0] = z1[0]-z0[0]; 	diff[1]=z1[1]-z0[1];
		z0[0]=z1[0];		z0[1]=z1[1];
		std::cout <<"("<<z1[0]<<")+i("<<z1[1]<<")\n";
		n++;
	} while(norm(diff)>TOL || norm(fz0)>TOL);
	
	// Print the solution
	std::cout << "\n(" << z1[0] << ")+i(" << z1[1] <<
	") is a complex solution to the equation f(z) = 0.\nNecessary iterations: " << n;
	std::cout << "\nError in the last iteration: " << norm(diff) << endl;
	
	return 0;
}
