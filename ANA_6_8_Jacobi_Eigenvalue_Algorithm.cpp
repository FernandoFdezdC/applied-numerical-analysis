// JACOBI EIGENVALUE ALGORITHM

// Date: 21/11/2020
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int n=3;
const double TOL=1e-12;

void print(double M[n][n], int n){			// Displays matrix M
	for (int i=0 ; i<n ; i++){
   		for (int j=0 ; j<n ; j++){
   			cout << M[i][j] << " ";
   		}
   		cout << "\n";
	}
}

int main(){
	cout << setprecision(9);
	
	double A[n][n];
	ifstream fentrada("Matrix.txt", ios::in);
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			fentrada >> A[i][j];
		}
	}
	fentrada.close();
 	
   	double U[n][n];		// Identity matrix
   	
   	for (int i=0 ; i<n ; i++){			// U = I
   		for (int j=i ; j<n ; j++){
   			if (i == j)
   				U[i][j] = 1;
			else{
   				U[i][j] = 0;
   				U[j][i] = 0;
   			}
   		
   		}
	}
   	
   	// Displays matrix U and A
   	cout << "U = \n";
   	print(U,n);
   	cout << "A = \n";
   	print(A,n);
   	cout << endl;
   	
   	int count = 0;
   	
   	int i,j;
   	
	double theta, aij, uki, ukj, aii;
	
	do{
	
	count = count + 1;
	
	aij = 0; 			// maximum value of matrix A above the principal diagonal
   	for (int x=0 ; x<n ; x++){
   		for (int y=x+1 ; y<n ; y++){
   			if (fabs(A[x][y]) > aij){
   				aij = fabs(A[x][y]);
   				i = x;
   				j = y;
			   }
   		}
	}
   	
   	if (fabs(A[i][i]-A[j][j]) < TOL*fabs(A[i][i]))
   		theta = atan(1.0);
	else
	   	theta = 0.5*atan(2*A[i][j]/(A[i][i]-A[j][j]));
	
	// Check the solution:
//	cout << "For count = " << count <<": \n";
//	cout << "theta = " << theta << "\n" << i << "\n" << j << "\n";
//	print(A);
//	cout << "Maximum value of A above principal diagonal: " << aij << ". Continue? "<< (aij>=TOL) << "\n";
   	
   	for (int k=0 ; k<n ; k++){
   		uki=U[k][i];
   		ukj=U[k][j];
   		U[k][i] = uki*cos(theta) + ukj*sin(theta);
   		U[k][j] = -uki*sin(theta) + ukj*cos(theta);
   		if (k!=i && k!=j){
	   		A[k][i] = A[i][k]*cos(theta) + A[j][k]*sin(theta);
	   		A[k][j] = -A[i][k]*sin(theta) + A[j][k]*cos(theta);
	   		A[j][k] = A[k][j];
	   		A[i][k] = A[k][i];
   		}
	   }
	aii = A[i][i];
   	A[i][i] = aii*pow(cos(theta),2) + A[j][j]*pow(sin(theta),2) + 2*A[i][j]*sin(theta)*cos(theta);
	A[j][j] = aii*pow(sin(theta),2) + A[j][j]*pow(cos(theta),2) - 2*A[i][j]*sin(theta)*cos(theta);
   	A[i][j] = 0;
   	A[j][i] = 0;
	
	} while (aij > TOL); 
   	
   	cout << "\nNumber of iterarions: " << count << endl;
   	cout << "The eigenvectors are the columns of the following matrix:\n";
   	print(U,n);
   	cout<<endl;
   	cout << "\nThe eigenvalues are:\n";
   	for (int j=0 ; j<=n-1 ; j++){
   		cout << A[j][j] << "\n";
	}
	cout << endl;
	print(A,n);
   	
	return 0;
}
