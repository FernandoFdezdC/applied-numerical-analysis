// GAUSS ELIMINATION METHOD TO SOLVE SYSTEMS OF EQUATIONS

// Date: 5/11/2020
// Author: Fernando Fernández del Cerro

#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int n=3;			// number of equations in the system

void print(double M[n][n], int n){			// Displays matrix M
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void print_vec(double b[n], int n){		// Displays vector b
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << b[i] << "\n";
	}
}

int main(){
	double M[n][n];				// Matrix M such that Mx=b
	ifstream file;
    file.open("Matrix.txt", ios::in);
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			file >> M[i][j];
		}
	}
	file.close();
	std::cout << "M=\n";
	print(M,n);			// Displays matrix M
	
   	double b[n] = {-1, 0, 3};		// vector b such that Mx=b
   	double x[n];					// Solution
   	
   	double pivote, temp;
   	int filapivote;
   	
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
	
	cout << "diag(M) = \n";
	print(M,n);			// Display diagonalized matrix M
	cout << "b = \n";
	print_vec(b,n);		// Display vector b when diagonalized
	cout << endl;
   	for (int j=0 ; j<=n-1 ; j++){		// Display the solution to the equation Ax=b
   		std::cout << "x_" << j+1 << " = " << x[j] << endl;
	}
	
	return 0;
}
