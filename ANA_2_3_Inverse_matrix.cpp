// GAUSS-JORDAN ELIMINATION METHOD APPLIED
// TO FIND THE INVERSE OF A SQUARE MATRIX

// Date: 28/2/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int n=3;			// dimension of the matrix
const int m=2*n;

void print(double M[n][m], int n, int m){	// Displays nxm matrix M
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=m-1 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void final_print(double M[n][m], int n, int m){	// Displays right side of augmented matrix M
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=n ; j<=m-1 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

int main(){
	cout << setprecision(9);
	
	double M[n][2*n];				// Matrix M such that Mx=b
	ifstream file;
    file.open("Matrix.txt", ios::in);
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			file >> M[i][j];
		}
	}
	file.close();
	// Augment matrix M with identity:
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			if (i==j){
				M[i][j+n]=1.0;
			}
			else{
				M[i][j+n]=0.0;
			}
		}
	}
	cout << "M=\n";
	print(M,n,2*n);			// Display matrix M
	
   	// Escalonar la matrix M:
	double pivote, temp;
   	int filapivote=0.0;
   	int columnapivote=0.0;
   	
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
   	filapivote=0;	columnapivote=0;
   	cout<<filapivote<<endl<<columnapivote<<endl;
   	// FORMA REDUCIDA:
   	while(filapivote<=n-1) {
		// bucle desde 0 hasta n-1 para recorrer todas las filas excepto la última.
		pivote=fabs(M[filapivote][columnapivote]);
		// cout << filapivote<<" "<<columnapivote<<endl;
		for(int j=columnapivote ; j<=n-1 ; j++) { // encuentra la fila pivote dentro de cada columna. pivoteo
			if (fabs(M[filapivote][j]) > 1e-50) {
				columnapivote=j;
				pivote = fabs(M[filapivote][columnapivote]);
				break;
			}
		} // final bucle en j
		for (int i=filapivote-1 ; i>=0 ; i--) { // Calcula y almacena las razones de coeficientes. Matriz L.
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
   	filapivote=0;	columnapivote=0;
   	cout<<filapivote<<endl<<columnapivote<<endl;
   	// NORMALIZATION:
   	while(filapivote<=n-1) {
		// bucle desde 0 hasta n-1 para recorrer todas las filas excepto la última.
		// cout << filapivote<<" "<<columnapivote<<endl;
		for(int j=columnapivote ; j<=n-1 ; j++) { // encuentra la fila pivote dentro de cada columna. pivoteo
			if (fabs(M[filapivote][j]) > 1e-50) {
				columnapivote=j;
				pivote = M[filapivote][columnapivote];
				// Normaliza la fila:
				for (int k=columnapivote ; k<=m-1 ; k++) { // Hace resta de filas
					M[filapivote][k]=M[filapivote][k]/pivote;
				} // final bucle en k
				break;
			}
		} // final bucle en j
		cout<<endl;		// Muestra el paso
		print(M,n,m);
		filapivote++;
   	} // final bucle en columnapivote
   	
	cout << "\ninverse(M) = \n";
	final_print(M,n,2*n);			// Display inverse of matrix M
	
	return 0;
}
