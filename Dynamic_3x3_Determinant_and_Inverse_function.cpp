// FINITE ELEMENTS FOR PARTIAL-DIFFERENTIAL EQUATIONS

// Date: 27/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n = 10; 			// number of subintervals (elements)


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


int main(){
	cout << setprecision(9);
	
	double M[3][3];				// Matrix M such that Mx=b
	ifstream file;
    file.open("Matrix.txt", ios::in);
	for (int i=0 ; i<=2 ; i++){
		for (int j=0 ; j<=2 ; j++){
			file >> M[i][j];
		}
	}
	file.close();
	
	cout << "M=\n";
	print(M);			// Display matrix M
	
	double detM = det(M); 		// Determinant of matrix M
	cout << "\ndet(M)= " << detM << "\n";
	
	double **invM;
	invM = new double *[3];
	for (int i=0 ; i<3 ; i++){
		invM[i] = new double[3];
	}
	inv(M,invM);
	
	cout<<"\ninv(M)=\n";      print33(invM);
	
	
	return 0;
}
