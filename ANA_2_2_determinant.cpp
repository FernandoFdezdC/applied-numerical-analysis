// GAUSS-JORDAN ELIMINATION METHOD TO FIND DETERMINANT
// OF SQUARE MATRIX

// Date: 28/2/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int n=3;			// dimension of the matrix

void print(double M[n][n], int n){			// Displays matrix M
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

int main(){
	cout << setprecision(9);
	
	double M[n][n];				// Matrix M such that Mx=b
	ifstream file;
    file.open("Matrix.txt", ios::in);
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			file >> M[i][j];
		}
	}
	file.close();
	
	cout << "M=\n";
	print(M,n);			// Display matrix M
	
	double det=1;		// at first, this carries the sign of the
						// determinant of matrix M
	
	// Escalonar la matriz:
   	double pivote, temp;
   	int filapivote=0;
   	int columnapivote=0;
   	
   	while(columnapivote<=n-2) {
		// bucle desde 0 hasta n-2 para recorrer todas las columnas excepto la última.
		pivote=fabs(M[filapivote][columnapivote]);
		// cout << filapivote<<" "<<columnapivote<<endl;
		if(pivote<1e-50){
			for(int i=filapivote+1 ; i<=n-1 ; i++) { // encuentra la fila pivote dentro de cada columna. pivoteo
				if (fabs(M[i][columnapivote]) > 1e-50) {
					// intercambia filas en caso de ser necesario
					for(int k=0 ; k<=n-1 ; k++) { // intercambia fila i y j
						temp = M[i][k];
						M[i][k] = M[filapivote][k];
						M[filapivote][k] = temp;
					}
					pivote = fabs(M[filapivote][columnapivote]);
					cout << endl;
					print(M,n);		// Muestra el intercambio
					det=-det; 		// el determinante cambia de signo
									// al intercambiar filas
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
				for (int k=columnapivote+1 ; k<=n-1 ; k++) { // Hace resta de filas
				// cout<<M[i][k]<<" "<<M[filapivote][k]<<" "<<M[i][columnapivote]<<" "<<M[filapivote][columnapivote]<<" ";
					M[i][k]=M[i][k]-M[filapivote][k]*M[i][columnapivote]/M[filapivote][columnapivote];
				// cout<<M[i][k]<<endl;
				} // final bucle en k
				// cout<<endl;
				M[i][columnapivote]=0;
			}
		} // final bucle en i
		cout<<endl;		// Muestra los pasos
		print(M,n);
		filapivote++;
   	} // final bucle en columnapivote
   	
	for(int i=0 ; i<=n-1 ; i++) {
		// bucle desde la primera hasta la última fila para recorrerlas todas y multiplicar
		det=det*M[i][i];		// Producto de los términos diagonales de
								// la matrix diagonalizada
	} // final bucle en i
   	
	cout << "\ndiagonalized(M) = \n";
	print(M,n);			// Display diagonalized matrix M
	cout << "det(M) = "<<det; // Display determinant of M
	
	return 0;
}
