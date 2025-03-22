// FINITE ELEMENTS FOR ORDINARY-DIFFERENTIAL EQUATIONS
// WITH NEUMANN CONDITIONS
// (if Q=0, then, the solution is not unique up to a constant)

// Date: 1/7/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n = 4; 			// number of subintervals (elements)

void solve(double M[n+1][n+1], double b[n+1], double *x, int n){
	// M[n+1][n+1], Matrix M such that Mx=b
    // x[n+1]: solution

   	double pivote, temp;
   	int filapivote;
   	
   	for(int j=0 ; j<=n-1 ; j++) {
		// bucle desde 0 hasta n-1 para recorrer todas las columnas excepto la última.
		pivote = fabs(M[j][j]);
		filapivote = j;
		for(int i=j+1 ; i<=n ; i++) { // encuentra la fila pivote dentro de cada columna. pivoteo
			if (fabs(M[i][j]) > pivote) {
				pivote=fabs(M[i][j]);
				filapivote=i;
			}
		} // final bucle en i
		if (filapivote != j ) { // intercambia filas en caso de ser necesario
			for(int k=0 ; k<=n ; k++) { // intercambia filas dadas por j y filapivote
				temp = M[j][k];
				M[j][k] = M[filapivote][k];
				M[filapivote][k] = temp;
			} // final bucle en k
			temp = b[j];
			b[j] = b[filapivote];
			b[filapivote] = temp;
		}
		for (int i=j+1 ; i<=n ; i++) { // Calcula y almacena las razones de coeficientes. Matriz L.
			M[i][j] = M[i][j]/M[j][j];
			for (int k=j+1 ; k<=n ; k++) { // Calcula los otros terminos, resultantes de hacer la resta
				M[i][k]=M[i][k]-M[i][j]*M[j][k];
			} // final bucle en k
			b[i]=b[i]-M[i][j]*b[j];
		} // final bucle en i
	} // final bucle en j (el del comienzo)
   	
   	// Back substitution --> Solution
	x[n] = b[n]/M[n][n];
	for (int j=n-1 ; j>=0 ; j--) {
		x[j]=b[j];
		for (int k=j+1 ; k<=n ; k++) {
			x[j]=x[j]-x[k]*M[j][k];
		} // final bucle en k
		x[j]=x[j]/M[j][j];
	} // final bucle en j
}

void print(double M[n+1][n+1], int n){		// Displays matrix M
	for (int i=0 ; i<=n ; i++){
		for (int j=0 ; j<=n ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void print_vec(double b[n+1], int n){		// Displays vector b
	for (int i=0 ; i<=n ; i++){
		std::cout << b[i] << "\n";
	}
}

double Q(double x){		// Function Q
	return -(x+1.0);
}

double F(double x){		// Function F
	return -exp(-x)*(x*x-x+2.0);
}

int main(){
	cout << setprecision(9);
	
	// ODE to solve:
	// y''+Q(x)y=F(x)
	
	double x_0=2.0, x_n=4.0;		// Boundary values
	double dydx_a=0.0, dydx_b=-2.0*exp(-4.0);		// Neumann conditions
	
	double x[n+1];			// nodes of the intervals
	x[0]=x_0;		x[n]=x_n;	// boundary values
	for (int k=1 ; k<=n-1 ; k++){
		x[k] = x_0 + k*(x_n-x_0)/n;		// evenly spaced nodes
				// nodes need not be evenly spaced
	}
	
	double h[n];			// size of elements
	for (int k=0 ; k<=n-1 ; k++){
		h[k] = x[k+1] -x[k];
	}
	
	std::cout << "The nodes are at:\n";
	print_vec(x,n);
	
	double Qav[n], Fav[n];		// Average values
	
	for (int k=0 ; k<=n-1 ; k++){
		Qav[k] = Q(x[k] + h[k]/2.0);	// Average value of Q \approx value of
				// of Q at the midpoint of the interval
		Fav[k] = F(x[k] + h[k]/2.0);	// Average value of F \approx value of
				// F at the midpoint of the interval
	}
	
	double K[n+1][n+1];		// System matrix that solves for all c_i's
	for (int i=0 ; i<=n ; i++){
		for (int j=0 ; j<=n ; j++){
			K[i][j] = 0.0;			// set matrix equal to null
		}
	}
	
	K[0][0] = 1.0/h[0] -Qav[0]*h[0]/3.0;
	K[0][1] = -1.0/h[0] -Qav[0]*h[0]/6.0;
	for (int k=1 ; k<=n-1 ; k++){
		K[k][k] = 1.0/h[k] + 1.0/h[k-1] - (Qav[k-1]*h[k-1] + Qav[k]*h[k])/3.0;
		K[k][k-1] = -1.0/h[k-1] -Qav[k-1]*h[k-1]/6.0;
		K[k][k+1] = -1.0/h[k] -Qav[k]*h[k]/6.0;
	}
	K[n][n] = 1.0/h[n-1] -Qav[n-1]*h[n-1]/3.0;
	K[n][n-1] = -1.0/h[n-1] -Qav[n-1]*h[n-1]/6.0;
	
	// Calculate vector b:
	double b[n+1];
	
	b[0] = -Fav[0]*h[0]/2.0;
	for (int k=1 ; k<=n-1 ; k++){
		b[k] = -Fav[k-1]*h[k-1]/2.0 -Fav[k]*h[k]/2.0;
	}
	b[n] = -Fav[n-1]*h[n-1]/2.0;
	
	cout << "\nSystem Matrix K=\n";
	print(K,n);
	cout<<"\nVector b=\n";
	print_vec(b,n);
	
	// APPLY NEUMANN BOUNDARY CONDITIONS:
	b[0] = b[0] - dydx_a;
	b[n] = b[n] + dydx_b;
	
	cout << "\nNew System Matrix K=\n";
	print(K,n);
	cout<<"\nNew vector b=\n";
	print_vec(b,n);
	
	// Obtain solution of the (tridiagonal) system of equations:
	double* u = new double[n+1];
	solve(K,b,u,n);
	
	cout << "\nTHE SOLUTION IS:\nu = [";
	cout << u[0] << ", ";
	for (int i=1 ; i<=n-1 ; i++){
	 	cout << u[i] <<", ";}
	cout << u[n] <<"];";
	
	return 0;
}
