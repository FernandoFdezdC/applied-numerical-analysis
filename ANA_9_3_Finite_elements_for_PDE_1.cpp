// FINITE ELEMENTS FOR PARTIAL-DIFFERENTIAL EQUATIONS

// Date: 27/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n = 5; 			// number of nodes


void inv(double M[3][3], double **invM){
	// M[3][3]: 3x3 matrix M
    // invM: Inverse matrix
	
	double detM = M[0][0]*M[1][1]*M[2][2] + M[0][1]*M[1][2]*M[2][0] + 
	M[2][1]*M[1][0]*M[0][2] - M[2][0]*M[1][1]*M[0][2] - 
	M[2][1]*M[1][2]*M[0][0] - M[1][0]*M[0][1]*M[2][2];	// determinant of matrix M
	
	invM[0][0] = (M[1][1]*M[2][2] -M[1][2]*M[2][1])/detM;
	invM[0][1] = -(M[0][1]*M[2][2] -M[0][2]*M[2][1])/detM;
	invM[0][2] = (M[0][1]*M[1][2] -M[0][2]*M[1][1])/detM;
	invM[1][0] = -(M[1][0]*M[2][2] -M[1][2]*M[2][0])/detM;
	invM[1][1] = (M[0][0]*M[2][2] -M[0][2]*M[2][0])/detM;
	invM[1][2] = -(M[0][0]*M[1][2] -M[0][2]*M[1][0])/detM;
	invM[2][0] = (M[1][0]*M[2][1] -M[1][1]*M[2][0])/detM;
	invM[2][1] = -(M[0][0]*M[2][1] -M[0][1]*M[2][0])/detM;
	invM[2][2] = (M[0][0]*M[1][1] -M[0][1]*M[1][0])/detM;
}

double det(double M[3][3]){
	// M[3][3]: 3x3 matrix M
    // detM: determinant of 3x3 matrix M
    
	double detM = M[0][0]*M[1][1]*M[2][2] + M[0][1]*M[1][2]*M[2][0] + 
	M[2][1]*M[1][0]*M[0][2] - M[2][0]*M[1][1]*M[0][2] - 
	M[2][1]*M[1][2]*M[0][0] - M[1][0]*M[0][1]*M[2][2];
	
	return detM;
}

void product(double **M, double c[n], double *a){
	// a = M*c
	for (int i=0 ; i<=2 ; i++){
		a[i] = 0.0;
		for (int j=0 ; j<=2 ; j++){
			a[i] = a[i] + M[i][j]*c[j];
		}
	}
}

void print33(double **M){		// Displays dynamic 3x3 matrix M
	for (int i=0 ; i<=2 ; i++){
		for (int j=0 ; j<=2 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void print(double M[3][3]){		// Displays 3x3 matrix M
	for (int i=0 ; i<=2 ; i++){
		for (int j=0 ; j<=2 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void print_vec3(double *b){		// Displays dynamic vector b of 3 elements
	for (int i=0 ; i<=2 ; i++){
		std::cout << b[i] << "\n";
	}
}

void print_vec(double b[3]){		// Displays vector b of 3 elements
	for (int i=0 ; i<=2 ; i++){
		std::cout << b[i] << "\n";
	}
}

//double Q(double x){		// Function Q
//	return -(x+1);
//}
//
//double F(double x){		// Function F
//	return -exp(-x)*(x*x-x+2);
//}

int main(){
	cout << setprecision(9);
	
	// Elliptic PDE to solve:
	// $u_{xx}+u_{yy}+Q(x,y)u=F(x,y)$
	// on region $R$ with boundary conditions $u(x,y)=u_0$ on $L_1$,
	// $\frac{\partial u}{\partial n}=\alpha u + \beta$ on $L_2$
	// Approximate solution: v(x,y)
	
	double c[3] = {100.0,200.0,300.0};
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
	cout << "M=\n";		print(M);			// Display matrix M
	
	double **invM;
	invM = new double *[3];
	for (int i=0 ; i<3 ; i++){
		invM[i] = new double[3];
	}
	inv(M,invM);
	
	cout<<"\ninv(M)=\n";      print33(invM);
	
	double* a = new double[3];
	product(invM,c,a);
	
	cout << "\na=\n";	print_vec3(a);
	
	return 0;
}
