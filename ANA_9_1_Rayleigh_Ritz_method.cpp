// RAYLEIGH-RITZ METHOD

// Date: 8/3/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n=4;			// number of constants to find

void u(double *c, double x, double a, double b){			// Approximated function u
	double u = c[0]+c[1]*(x-a)+c[2]*(x-a)*(x-b)+c[3]*(x-a)*(x-a)*(x-b);
	// u must be the sum of linearly independent functions and
	// the boundary conditions must be met for u
}

double Q(double x){		// Function Q
	return 1.0;
}

double F(double x){		// Function F
	return 3.0*x*x;
}

double integralA00(double a, double b, double *x, double *w){
	// integral from a to b of function Q(x-a)^2(x-b)^2
	double I = 0.0;		// Integral of function from a to b
	// Calculate integral I:
	double zeta;
	for(int i=0 ; i<=4 ; i++){
		zeta = (b-a)*x[i]/2.0 + (a+b)/2.0;
		I = I + w[i]*Q(zeta)*pow((zeta-a),2.0)*pow((zeta-b),2.0);
	}
	I = (b-a)*I/2.0;
	return I;
}
double integralA01(double a, double b, double *x, double *w){
	// integral from a to b of function Q(x-a)^3(x-b)^2
	double I = 0.0;		// Integral of function from a to b
	// Calculate integral I:
	double zeta;
	for(int i=0 ; i<=4 ; i++){
		zeta = (b-a)*x[i]/2.0 + (a+b)/2.0;
		I = I + w[i]*Q(zeta)*pow((zeta-a),3.0)*pow((zeta-b),2.0);
	}
	I = (b-a)*I/2.0;
	return I;
}
double integralA11(double a, double b, double *x, double *w){
	// integral from a to b of function Q(x-a)^4(x-b)^2
	double I = 0.0;		// Integral of function from a to b
	// Calculate integral I:
	double zeta;
	for(int i=0 ; i<=4 ; i++){
		zeta = (b-a)*x[i]/2.0 + (a+b)/2.0;
		I = I + w[i]*Q(zeta)*pow((zeta-a),4.0)*pow((zeta-b),2.0);
	}
	I = (b-a)*I/2.0;
	return I;
}
double integral01(double a, double b, double *x, double *w){
	// integral from a to b of function Q(x-a)(x-b)
	double I = 0.0;		// Integral of function from a to b
	// Calculate integral I:
	double zeta;
	for(int i=0 ; i<=4 ; i++){
		zeta = (b-a)*x[i]/2.0 + (a+b)/2.0;
		I = I + w[i]*Q(zeta)*(zeta-a)*(zeta-b);
	}
	I = (b-a)*I/2.0;
	return I;
}
double integral02(double a, double b, double *x, double *w){
	// integral from a to b of function Q(x-a)^2(x-b)
	double I = 0.0;		// Integral of function from a to b
	// Calculate integral I:
	double zeta;
	for(int i=0 ; i<=4 ; i++){
		zeta = (b-a)*x[i]/2.0 + (a+b)/2.0;
		I = I + w[i]*Q(zeta)*pow((zeta-a),2.0)*(zeta-b);
	}
	I = (b-a)*I/2.0;
	return I;
}
double integral03(double a, double b, double *x, double *w){
	// integral from a to b of function F(x-a)(x-b)
	double I = 0.0;		// Integral of function from a to b
	// Calculate integral I:
	double zeta;
	for(int i=0 ; i<=4 ; i++){
		zeta = (b-a)*x[i]/2.0 + (a+b)/2.0;
		I = I + w[i]*F(zeta)*(zeta-a)*(zeta-b);
	}
	I = (b-a)*I/2.0;
	return I;
}
double integral12(double a, double b, double *x, double *w){
	// integral from a to b of function Q(x-a)^3(x-b)
	double I = 0.0;		// Integral of function from a to b
	// Calculate integral I:
	double zeta;
	for(int i=0 ; i<=4 ; i++){
		zeta = (b-a)*x[i]/2.0 + (a+b)/2.0;
		I = I + w[i]*Q(zeta)*pow((zeta-a),3.0)*(zeta-b);
	}
	I = (b-a)*I/2.0;
	return I;
}
double integral13(double a, double b, double *x, double *w){
	// integral from a to b of function F(x-a)^2(x-b)
	double I = 0.0;		// Integral of function from a to b
	// Calculate integral I:
	double zeta;
	for(int i=0 ; i<=4 ; i++){
		zeta = (b-a)*x[i]/2.0 + (a+b)/2.0;
		I = I + w[i]*F(zeta)*pow((zeta-a),2.0)*(zeta-b);
	}
	I = (b-a)*I/2.0;
	return I;
}

double* solve(double M[n-2][n-2], double *b, double *x){	
	// Solve system of equations
	// M[2][2]: Matrix M such that M*sol=b
    // x[2]: solution

   	double pivote, temp;
   	int filapivote;
   	int n=2; 	// number of equations
   	
   	for(int j=0 ; j<=n-2 ; j++) {
		// bucle desde 0 hasta n-2 para recorrer todas las columnas excepto la última.
		pivote = fabs(M[j][j]);
		filapivote = j;
		for(int i=j+1 ; i<=n-1 ; i++) { // encuentra la fila pivote dentro de cada columna. pivoteo
			if (fabs(M[i][j]) > pivote) {
				pivote=fabs(M[i][j]);
				filapivote=i;
			}
		} // final bucle en i
		if (filapivote != j ) { // intercambia filas en caso de ser necesario
			for(int k=0 ; k<=n-1 ; k++) { // intercambia filas dadas por j y filapivote
				temp = M[j][k];
				M[j][k] = M[filapivote][k];
				M[filapivote][k] = temp;
			} // final bucle en k
			temp = b[j];
			b[j] = b[filapivote];
			b[filapivote] = temp;
		}
		for (int i=j+1 ; i<=n-1 ; i++) { // Calcula y almacena las razones de coeficientes. Matriz L.
			M[i][j] = M[i][j]/M[j][j];
			for (int k=j+1 ; k<=n-1 ; k++) { // Calcula los otros terminos, resultantes de hacer la resta
				M[i][k]=M[i][k]-M[i][j]*M[j][k];
			} // final bucle en k
			b[i]=b[i]-M[i][j]*b[j];
		} // final bucle en i
	} // final bucle en j (el del comienzo)
   	
   	// Back substitution --> Solution
	x[n-1] = b[n-1]/M[n-1][n-1];
	for (int j=n-2 ; j>=0 ; j--) {
		x[j]=b[j];
		for (int k=j+1 ; k<=n-1 ; k++) {
			x[j]=x[j]-x[k]*M[j][k];
		} // final bucle en k
		x[j]=x[j]/M[j][j];
	} // final bucle en j
	
	return x;
}

int main(){
	cout << setprecision(9);
	
	// ODE to solve:
	// y''+Q(x)y=F(x)
	// Functional of the equation:
	// I[u]=\int_{a}^{b} \left[\left(\frac{du}{dx}\right)^2-Qu^2+2Fu\right] dx
	// Supose u(x) = c_0v_0(x)+c_1v_1(x)+...+c_nv_n(x)
	
	double a=0.0, b=2.0;			// boundary points
	double u_a = 0.0, u_b = 3.5;	// Dirichlet boundary conditions
	
	double c[n];		// optimal constants to find
	// Supose u(x)=c_0+c_1(x-a)+c_2(x-a)(x-b)+c_3(x-a)^2(x-b)
	
	c[0] = u_a ; c[1] = (u_b-u_a)/(b-a);
	
	double A[n-2][n-2];		// coefficient matrix of the system that gives c_2 and c_3
	double f[n-2];
	
	// Set points and weights of 5 point Gaussian quadrature:
	double x[5]={-pow(5.0+2.0*pow(10.0/7.0,0.5),0.5)/3.0, -pow(5.0-2.0*pow(10.0/7.0,0.5),0.5)/3.0, 0, pow(5.0-2.0*pow(10.0/7.0,0.5),0.5)/3.0, pow(5.0+2.0*pow(10.0/7.0,0.5),0.5)/3.0};
	double w[5]={(322.0-13.0*pow(70.0,0.5))/900.0, (322.0+13.0*pow(70.0,0.5))/900.0, 128.0/225.0, (322.0+13.0*pow(70.0,0.5))/900.0, (322.0-13.0*pow(70.0,0.5))/900.0};
	
	
	A[0][0] = integralA00(a,b,x,w)-pow(b-a,3.0)/3.0;
	A[0][1] = integralA01(a,b,x,w)-pow(b-a,4.0)/6.0;
	A[1][0] = integralA01(a,b,x,w)-pow(b-a,4.0)/6.0;
	A[1][1] = integralA11(a,b,x,w)-2.0*pow(b-a,5.0)/15;
	f[0] = -c[0]*integral01(a,b,x,w)-c[1]*integral02(a,b,x,w)+integral03(a,b,x,w);
	f[1] = -c[0]*integral02(a,b,x,w)-c[1]*integral12(a,b,x,w)+integral13(a,b,x,w);
	
	// System of equations to solve:
	cout << A[0][0]<<"c_2+"<<A[0][1]<<"c_3="<<f[0]<<"\n";
	cout << A[1][0]<<"c_2+"<<A[1][1]<<"c_3="<<f[1]<<"\n";
	
	// Obtain solution of the system:
	double sol[n-2];
	double* cons= solve(A, f, sol);
	c[2]=cons[0];
	c[3]=cons[1];
	cout <<"\n"<<c[2]<<"\n"<<c[3]<<"\n";
	
	// Coefficients of the polynomial that best fits the solution:
	double p[4]={c[0]-a*c[1]+a*b*c[2]-a*a*b*c[3], c[1]-a*c[2]-b*c[2]+(a*a+2.0*a*b)*c[3], c[2]-(2.0*a+b)*c[3], c[3]};
	cout<<"\nCoefficients:\n"<<p[0]<<"\n"<<p[1]<<"\n"<<p[2]<<"\n"<<p[3];
	
	return 0;
}
