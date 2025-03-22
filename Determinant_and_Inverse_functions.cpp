// FINITE ELEMENTS FOR PARTIAL-DIFFERENTIAL EQUATIONS

// Date: 26/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n = 3;				// dimension of matrix
const int m = 2*n;
//const int n = 10; 			// number of subintervals (elements)




void inv(double A[n][n], double **invA, int n){
	// A[n][n], Matrix A
    // invA: Inverse matrix
	
	double M[n][2*n];		// augmented matrix
	// Augment matrix M with identity:
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			M[i][j] = A[i][j] ; }}
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			if (i==j){
				M[i][j+n]=1.0;}
			else{
				M[i][j+n]=0.0;}
		}
	}
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
		filapivote++;
   	} // final bucle en columnapivote
   	filapivote=0;	columnapivote=0;
   	// cout<<filapivote<<endl<<columnapivote<<endl;
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
		filapivote++;
   	} // final bucle en columnapivote
   	filapivote=0;	columnapivote=0;
   	// cout<<filapivote<<endl<<columnapivote<<endl;
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
		filapivote++;
   	} // final bucle en columnapivote
	
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			invA[i][j] = M[i][j+n];
		}
	}
}

double det(double M[n][n], int n){
	// M[n][n]: Matrix M
    // det: determinant of matrix M
	double det=1.0;		// at first, this carries the sign of the
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
		filapivote++;
   	} // final bucle en columnapivote
   	
	for(int i=0 ; i<=n-1 ; i++) {
		// bucle desde la primera hasta la última fila para recorrerlas todas y multiplicar
		det=det*M[i][i];		// Producto de los términos diagonales de
								// la matrix diagonalizada
	} // final bucle en i
	
	return det;
}

double* product(double **M, double c[n], double *x, int n){
	// x = M*c
	for (int i=0 ; i<=n-1 ; i++){
		x[i] = 0.0;
		for (int j=0 ; j<=n-1 ; j++){
			x[i] = x[i] + M[i][j]*c[j];
		}
	}
   	return x;
}

void print33(double **M){		// Displays dynamic 3x3 matrix M
	for (int i=0 ; i<=2 ; i++){
		for (int j=0 ; j<=2 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void print(double M[n][n], int n){		// Displays nxn matrix M
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void* print_vec3(double *b){		// Displays dynamic vector b of 3 elements
	for (int i=0 ; i<=2 ; i++){
		std::cout << b[i] << "\n";
	}
}

void print_vec(double b[n], int n){		// Displays vector b
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << b[i] << "\n";
	}
}

double Q(double x, double y){		// Function Q(x,y)
	return x*y/2.0;
}

double F(double x, double y){		// Function F(x,y)
	return x+y;
}

int main(){
	cout << setprecision(9);
	
	// Elliptic PDE to solve:
	// $u_{xx}+u_{yy}+Q(x,y)u=F(x,y)$
	// on region $R$ with boundary conditions $u(x,y)=u_0$ on $L_1$,
	// $\frac{\partial u}{\partial n}=\alpha u + \beta$ on $L_2$
	// Approximate solution: v(x,y)
	
	double c[3] = {100.0, 200.0, 300.0};
	double r[2]={0.0,0.0}, t[2]={0.0,1.0}, s[2]={2.0,0.0}; // Coordinates of the triangle vertices 
	
	double M[3][3]; 			// Matrix M
	for (int i=0 ; i<=2 ; i++){
		M[i][0] = 1.0;
	}
	for (int i=0 ; i<=1 ; i++){
		M[0][i+1] = r[i];
		M[1][i+1] = s[i];
		M[2][i+1] = t[i];
	}
	cout << "M=\n";
	print(M,n);
	
	double **invM;
	invM = new double *[3];
	for (int i=0 ; i<3 ; i++){
		invM[i] = new double[3];
	}
	inv(M,invM,n);
	
	cout<<"\n";
	print33(invM);
	
	double x[n];
	double* a = product(invM,c,x,n);
	cout << "\n";
	print_vec3(a);
	
	double A[n], B[n], C[n];
	for (int i=0 ; i<=2 ; i++){
		A[i] = invM[0][i];
		B[i] = invM[1][i];
		C[i] = invM[2][i];
	}
	
	double detM = det(M,n);		// Calculate determinant of matrix M
	cout << "\ndet(M)= "<<detM << "\n";
	
	double xc=(r[0]+s[0]+t[0])/3.0,
	yc=(r[1]+s[1]+t[1])/3.0;			// Position of centroid
	
	double Qav = Q(xc,yc),
	Fav = F(xc,yc);				// Average values of F and Q
	
	cout << "\nQ_{av}= " << Qav << "\nF_{av}= " << Fav<<"\n\n";
	
	// CONSTRUCT MATRIX K AND VECTOR b:
	double K[n][n], b[n];
	for (int j=0 ; j<=2 ; j++){
		// Area = 0.5*det(M)
		b[j] = -0.5*detM*Fav/3.0;
		for (int k=0 ; k<=2 ; k++){
			if (j==k){
				K[j][j] = 0.5*detM*(B[j]*B[j] + C[j]*C[j] - Qav/6.0);
			}
			else {
				K[j][k] = 0.5*detM*(B[j]*B[k] + C[j]*C[k] - Qav/12.0);
			}
		}
	}
	
	print(K,n);
	cout << "\n";
	print_vec(b,n);
	
	// ASSEMBLE THE SYSTEM:
	
	
	
	
	
	
	
	
	
	
	
	return 0;
}

