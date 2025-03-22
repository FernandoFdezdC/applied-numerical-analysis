// CUBIC SPLINES. 3RD CONDITION.

// Date: 9/3/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n=5; 			// amount of data to be interpolated

double int_pol(double x0, double *x, double a[4][n-1]){
	// Calculate at x0 the linear spline
	// that best fits the given data 
	double P;
	if (x0<x[0]){ P=a[0][0]+a[1][0]*(x0-x[0])+a[2][0]*pow((x0-x[0]),2.0)+a[3][0]*pow((x0-x[0]),3.0); }
	else if (x0>x[n-2]){ P=a[0][n-2]+a[1][n-2]*(x0-x[n-2])+a[2][n-2]*pow((x0-x[n-2]),2.0)+a[3][n-2]*pow((x0-x[n-2]),3.0); }
	else{
		for (int i=0 ; i<=n-2 ; i++){
	   		if ((x0-x[i])*(x0-x[i+1])<0){
	   				P = a[0][i]+a[1][i]*(x0-x[i])+a[2][i]*pow((x0-x[i]),2.0)+a[3][i]*pow((x0-x[i]),3.0);
			   	}
			else if((x0-x[i])==0){ P=a[0][i]; }
		}
	}
	return P;		// returned value
}

void print(double M[n-2][n-2], int n){			// Displays matrix M
	for (int i=0 ; i<=n-3 ; i++){
		for (int j=0 ; j<=n-3 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}
void print_4xn_1(double M[4][n-1],int n){
	for (int i=0 ; i<=3 ; i++){
		for (int j=0 ; j<=n-2 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}
void print_vec(double b[n-2], int n){		// Displays vector b
	for (int i=0 ; i<=n-3 ; i++){
		std::cout << b[i] << "\n";
	}
}

double* solve(double M[n-2][n-2], double *b, double *x, int n){	
	// Solve system of equations
	// M[n-2][n-2]: Matrix M such that M*sol=b
    // x[n-2]: solution

   	double pivote, temp;
   	int filapivote;
   	n=n-2; 	// number of equations
   	
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
	n=n+2;
	return x;
}

int main(){
	cout << setprecision(9);
	
	// g_i(x)=a_i(x-x_i)^3+b_i(x-x_i)^2+c_i(x-x_i)+d_i
	
	double x[n]={0.0,1.0,1.25,1.5,2.25},		// data to be interpolated
	y[n] = {2.0000,4.4366,20.0,6.7134,13.9130};		// data to be interpolated
	
	double a[4][n-1]; 		// matrix with coefficients of the cubic spline
	double h[n-1];			// vector with width of the i-th interval in i-th position
	for (int i=0 ; i<=n-2 ; i++){ h[i] = x[i+1]-x[i]; }
	
	// CONSTRUCT CUBIC SPLINE:
	double f[n-2];
	for (int i=0 ; i<=n-3 ; i++){
		f[i]=6.0*((y[i+2]-y[i+1])/(x[i+2]-x[i+1])-(y[i+1]-y[i])/(x[i+1]-x[i]));
	}
	double K[n-2][n-2];
	for (int i=0 ; i<=n-3 ; i++){
		for (int j=0 ; j<=n-3 ; j++){
			K[i][j]=0;		// Set matrix equal to null
		}
	}
	K[0][1]=h[1];
	K[0][0]=3*h[0]+2*h[1];
	for (int i=1 ; i<=n-4 ; i++){
		K[i][i-1]=h[i];
		K[i][i]=2*(h[i]+h[i+1]);		// Set tridiagonal matrix
		K[i][i+1]=h[i+1];
	}
	K[n-3][n-4]=h[n-3];			// Set tridiagonal matrix
	K[n-3][n-3]=2*h[n-3]+3*h[n-2];
	print(K,n);
	print_vec(f,n);
	// solution of the tridiagonal system:
	double sol[n-2];
   	double* solution=solve(K,f,sol,n);
   	double S[n];
   	for(int i=1 ; i<=n-2 ; i++){
		S[i]=solution[i-1];}
	S[0]=S[1];
	S[n-1]=S[n-2];
   	
   	for (int i=0 ; i<=n-2 ; i++){
   		a[0][i] = y[i];						// d_i
	   	a[1][i] = (y[i+1]-y[i])/(h[i])-(2*h[i]*S[i]+h[i]*S[i+1])/6.0;		// c_i
   		a[2][i] = S[i]/2.0;					// b_i
   		a[3][i] = (S[i+1]-S[i])/(6*h[i]);	// a_i
	}
	// Construction of cubic spline. END.
	cout<<endl;
	print_4xn_1(a,n);
	cout<<endl;
	// Interpolated values in interval [0,2.5]:
	double x0=0.0;	
	double val;
	cout<<"[";
	do{
	val = int_pol(x0,x,a); // returned value
	// cout <<"\n\nValue at "<<x0<<" of the interpolating polynomial: "<<val;
	cout <<val<<", ";
	x0=x0+0.01;
	}while (x0<=2.49);
//	cout<<endl<<x0<<endl;
	val = int_pol(x0,x,a); // returned value
	cout <<val<<"]";
	
	return 0;
}
