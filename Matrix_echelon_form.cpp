// FIND ECHELON FORM OF AN NXM MATRIX

// Date: 1/3/2020
// Author: Fernando Fernández del Cerro

#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int n=3;			// dimension of the matrix
const int m=4;

void print(double M[n][m], int n, int m){	// Displays matrix nxm M
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=m-1 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

int main(){
	double M[n][m];				// Matrix M such that Mx=b
	ifstream file;
    file.open("Matrix.txt", ios::in);
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=m-1 ; j++){
			file >> M[i][j];
		}
	}
	file.close();
	cout << "M=\n";
	print(M,n,m);			// Display matrix M
	
	// Escalonar la matrix M:
	double pivote, temp;
   	int filapivote=0;
   	int columnapivote=0;
   	
   	while(columnapivote<=m-2) {
		// bucle desde 0 hasta n-2 para recorrer todas las columnas excepto la última.
		pivote=fabs(M[filapivote][columnapivote]);
		// cout << filapivote<<" "<<columnapivote<<endl;
		if(pivote<1e-50){
			for(int i=filapivote+1 ; i<=n-1 ; i++) { // encuentra la fila pivote dentro de cada columna. pivoteo
				if (fabs(M[i][columnapivote]) > 1e-50) {
					// intercambia filas en caso de ser necesario
					for(int k=0 ; k<=m-1 ; k++) { // intercambia fila i y j
						temp = M[i][k];
						M[i][k] = M[filapivote][k];
						M[filapivote][k] = temp;
					}
					pivote = fabs(M[filapivote][columnapivote]);
					cout << endl;
					print(M,n,m);		// Muestra el intercambio
					break;
				}
			} // final bucle en i
		}
		if(pivote<1e-50){
			columnapivote++;
			continue;
		}
		for (int i=filapivote+1 ; i<=n-1 ; i++) { // Calcula y almacena las razones de coeficientes. Matriz L.
			if(fabs(M[i][columnapivote])>1e-50){
				for (int k=columnapivote+1 ; k<=m-1 ; k++) { // Hace resta de filas
				// cout<<M[i][k]<<" "<<M[filapivote][k]<<" "<<M[i][columnapivote]<<" "<<M[filapivote][columnapivote]<<" ";
					M[i][k]=M[i][k]-M[filapivote][k]*M[i][columnapivote]/M[filapivote][columnapivote];
				// cout<<M[i][k]<<endl;
				} // final bucle en k
				// cout<<endl;
				M[i][columnapivote]=0;
			}
		} // final bucle en i
		cout<<endl;		// Muestra los pasos
		print(M,n,m);
		filapivote++;
   	} // final bucle en columnapivote
   	
	cout << "\ndiagonalized(M) = \n";
	print(M,n,m);			// Display diagonalized matrix M
	
	return 0;
}
