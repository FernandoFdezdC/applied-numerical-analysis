// A BETTER RATIONAL FUNCTION FOR exp(x)
// NOT COMPLETE

// Date: 24/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int n = 2;			// Degree of the numerator
const int m = 1;			// Degree of the denominator
const int N = n+m;

void print(double M[N][N], int N){			// Displays NxN matrix M
	for (int i=0 ; i<=N-1 ; i++){
		for (int j=0 ; j<=N-1 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void print_vecN(double d[N], int N){		// Displays vector d
	for (int i=0 ; i<=N-1 ; i++){
		std::cout << d[i] << "\n";
	}
}

void print_vecn(double a[n+1], int n){		// Displays vector a
	for (int i=0 ; i<=n ; i++){
		std::cout << a[i] << "\n";
	}
}

void print_vecm(double b[m+1], int m){		// Displays vector b
	for (int i=0 ; i<=m ; i++){
		std::cout << b[i] << "\n";
	}
}

double* solve(double M[N][N], double *d, double *x, int N){	
	// Solve system of equations
	// M[N][N]: Matrix M such that M*sol=d
    // x[N]: solution

   	double pivote, temp;
   	int filapivote;
   	
   	for(int j=0 ; j<=N-2 ; j++) {
		// bucle desde 0 hasta N-2 para recorrer todas las columnas excepto la última.
		pivote = fabs(M[j][j]);
		filapivote = j;
		for(int i=j+1 ; i<=N-1 ; i++) { // encuentra la fila pivote dentro de cada columna. pivoteo
			if (fabs(M[i][j]) > pivote) {
				pivote=fabs(M[i][j]);
				filapivote=i;
			}
		} // final bucle en i
		if (filapivote != j ) { // intercambia filas en caso de ser necesario
			for(int k=0 ; k<=N-1 ; k++) { // intercambia filas dadas por j y filapivote
				temp = M[j][k];
				M[j][k] = M[filapivote][k];
				M[filapivote][k] = temp;
			} // final bucle en k
			temp = d[j];
			d[j] = d[filapivote];
			d[filapivote] = temp;
		}
		for (int i=j+1 ; i<=N-1 ; i++) { // Calcula y almacena las razones de coeficientes. Matriz L.
			M[i][j] = M[i][j]/M[j][j];
			for (int k=j+1 ; k<=N-1 ; k++) { // Calcula los otros terminos, resultantes de hacer la resta
				M[i][k]=M[i][k]-M[i][j]*M[j][k];
			} // final bucle en k
			d[i]=d[i]-M[i][j]*d[j];
		} // final bucle en i
	} // final bucle en j (el del comienzo)
   	
   	// Back substitution --> Solution
	x[N-1] = d[N-1]/M[N-1][N-1];
	for (int j=N-2 ; j>=0 ; j--) {
		x[j]=d[j];
		for (int k=j+1 ; k<=N-1 ; k++) {
			x[j]=x[j]-x[k]*M[j][k];
		} // final bucle en k
		x[j]=x[j]/M[j][j];
	} // final bucle en j
	return x;
}

// double f(double x){
// 	return exp(x);
// }

int main(){
	cout << setprecision(9);
	
	double a[n+1], b[m+1];		// Coefficients of the Padé polynomials
	
	b[0] = 1.0;
	double c[N+1] = {1.2661,1.1302,0.2715,0.0443}; // Chebyshev first N+1 coefficients of the series of f(x)
	a[0] = c[0];
	
	// State system of equations to solve
	double d[N];
	for (int i=0 ; i<=N-1 ; i++){
		d[i]=-c[i+1];
	}
	double M[N][N];
	for (int i=0 ; i<=N-1 ; i++){
		for (int j=0 ; j<=N-1 ; j++){
			M[i][j]=0;		// Set matrix equal to null
		}
	}
	for (int i=0 ; i<=n-1 ; i++){
		M[i][i]=-1;		// State certain coefficients of the matrix
	}
	// cout <<b[0]<<endl;
	for (int i=0 ; i<=N-1 ; i++){
		for (int j=n ; j<=N-1 && j<=n+i ; j++){
			M[i][j] = c[i-j+n];
		}
	}
	// cout << endl<<b[0]<<endl;
	print(M,N); cout<<"\n";
	print_vecN(d,N);
	// solution of the system of equations:
	double sol[N];
   	double* solution=solve(M,d,sol,N);
   	for(int i=1 ; i<=n ; i++){
		a[i]=solution[i-1];}
	for(int i=1 ; i<=m ; i++){
		b[i]=solution[n+i-1];}
	
	cout<<"\na_{i}=\n"; print_vecn(a,n);
	cout<<"\nb_{i}=\n"; print_vecm(b,m);
	
	return 0;
}
