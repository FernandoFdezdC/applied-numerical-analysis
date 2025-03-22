// Lagrange polynomials of order n
// with Gauss elimination method

// Date: 26/2/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n = 11; 			// amount of data to be interpolated

void print(double M[n][n], int n){			// Displays matrix M
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void print_vec(double b[n], int n){			// Displays vector b
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << b[i] << "\n";
	}
}

double func(double x){		// Function that fits the data
	return pow(cos(x),10);
}
double pn(double a[n], double x){	// Evaluate the value of the Lagrange interpolating polynomial in x
	double pn = a[0];
	for (int i=1 ; i<=n-1 ; i++){
		pn = pn + a[i]*pow(x,i);
	}
  	return pn;
}
double error(double a[n], double x){
	double error = func(x)-pn(a,x);	// error function E(x)=f(x)-P_n(x)
	return error;
}


int main(){
	cout << setprecision(9);
	
	double x[n], f[n];		// data to be interpolated
	double x0=-2, xn=2;		// interval for interpolation
	
	for (int i=0 ; i<=n-1 ; i++){
		x[i]=(xn-x0)*i/(n-1)+x0;
		f[i]=func(x[i]);
	}
	
	double V[n][n];				// Vandermonde matrix V such that V\vec{a}=\vec{f}
	
	// Vandermonde matrix:
	for (int i=0 ; i<=n-1 ; i++){
		V[i][0]=1;
		for(int j=1 ; j<=n-1 ; j++){
			V[i][j] = pow(x[i],j);
		}
	}
	cout << "V =\n";
	print(V,n);				// Display matrix V
	
	double a[n];			// Solution of the system of equations
	
	
	// Solve V\vec{a}=\vec{f} with Gauss elimination method:
	double pivote, temp;
   	int filapivote;
   	for(int j=0 ; j<=n-2 ; j++) {
		// bucle desde 0 hasta n-2 para recorrer todas las columnas excepto la última.
		pivote = fabs(V[j][j]);
		filapivote = j;
		for(int i=j+1 ; i<=n-1 ; i++) { // encuentra la fila pivote dentro de cada columna. pivoteo
			if (fabs(V[i][j]) > pivote) {
				pivote=fabs(V[i][j]);
				filapivote=i;
			}
		} // final bucle en i
		if (filapivote != j ) { // intercambia filas en caso de ser necesario
			for(int k=0 ; k<=n-1 ; k++) { // intercambia filas dadas por j y filapivote
				temp = V[j][k];
				V[j][k] = V[filapivote][k];
				V[filapivote][k] = temp;
			} // final bucle en k
			temp = f[j];
			f[j] = f[filapivote];
			f[filapivote] = temp;
		}
		for (int i=j+1 ; i<=n-1 ; i++) { // Calcula y almacena las razones de coeficientes. Matriz L.
			V[i][j] = V[i][j]/V[j][j];
			for (int k=j+1 ; k<=n-1 ; k++) { // Calcula los otros terminos, resultantes de hacer la resta
				V[i][k]=V[i][k]-V[i][j]*V[j][k];
			} // final bucle en k
			f[i]=f[i]-V[i][j]*f[j];
		} // final bucle en i
	} // final bucle en j (el del comienzo)
   	// Back substitution --> Solution
	a[n-1] = f[n-1]/V[n-1][n-1];
	for (int j=n-2 ; j>=0 ; j--) {
		a[j]=f[j];
		for (int k=j+1 ; k<=n-1 ; k++) {
			a[j]=a[j]-a[k]*V[j][k];
		} // final bucle en k
		a[j]=a[j]/V[j][j];
	} // final bucle en j
	
	cout << "Polynomial coefficients a_0 + a_1 x + ... + a_n x^n:\n";
	
	for (int j=0 ; j<=n-1 ; j++){		// Display the solution to the equation Ax=b
   		std::cout << "a_" << j << " = " << a[j] << endl;
	}
	
	// Open a file:
	ofstream file;
	file.open ("ANA_Lagrange_polynomials.txt");
	// Write the solution in a file:
	file << "When f(x)=cos(x)^10, the Lagrange coefficients with "<< n<< " points between " << x0 << " and " << xn <<" are:\n";
	for (int j=0 ; j<=n-1 ; j++){
   		file << "a_" << j << " = " << a[j] << endl;
	}
	// For matlab:
	file << "[";
	cout << "[";
	for (int j=0 ; j<=n-2 ; j++){
   		file << a[j] << ", ";
   		cout << a[j] << ", ";
	}
	file << a[n-1] << "];";
	cout << a[n-1] << "];";
	file.close();
	
	return 0;
}
