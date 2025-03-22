// FINITE ELEMENTS FOR PARTIAL-DIFFERENTIAL EQUATIONS

// Date: 27/6/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n = 5; 			// number of nodes
const int m = 3;			// number of elements


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

void printn(double M[n][n], int n){			// Displays nxn matrix M
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void print_vecn(double b[n], int n){		// Displays nxn vector b
	for (int i=0 ; i<=n-1 ; i++){
		std::cout << b[i] << "\n";
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
	
	
	double r[2]={0.0,0.0}, s[2]={2.0,0.0}, t[2]={0.0,1.0};		// Coordinates of the triangle vertices 
	
	double pr[m][2], ps[m][2], pt[m][2];
	for (int p=0 ; p<=m-1 ; p++){
		// The first [p] entry refers to the element index
		for (int q=0 ; q<=1 ; q++){
			pr[p][q] = r[q];
			ps[p][q] = s[q];
			pt[p][q] = t[q];
		}
	}
	
	int o[m][3]; 			// Order of the points in each element
	// Element 1:
	o[0][0] = 1; 	o[0][1] = 2; 	o[0][2] = 4;
	// Element 2:
	o[1][0] = 2; 	o[1][1] = 3; 	o[1][2] = 4;
	// Element 3:
	o[2][0] = 4; 	o[2][1] = 5; 	o[2][2] = 1;
	
	
	// CALCULATE K FOR EVERY ELEMENT:
	
	double pK[m][3][3], 		// Matrix K as a function of element p (first entry)
			pb[m][3];			// Vector b as a function of element p (first entry)
	double M[3][3];
	double detM;
	double **invM;
	double A[3], B[3], C[3];
	double xc, yc;
	double Qav, Fav;
	double K[3][3], b[3];
	
	for (int p=0 ; p<=m-1 ; p++){
		// The first [p] entry refers to the element index
		cout << "ELEMENT " << p+1 << ":\n";
		for (int i=0 ; i<=2 ; i++){
			M[i][0] = 1.0;
		}
		for (int i=0 ; i<=1 ; i++){
			M[0][i+1] = pr[p][i];
			M[1][i+1] = ps[p][i];
			M[2][i+1] = pt[p][i];
		}
		cout << "M=\n"; 	print(M);			// Display matrix M
		
		detM = det(M); 		// Determinant of matrix M
		cout << "\ndet(M)= " << detM << "\n";
		
//		double **invM;
		invM = new double *[3];
		for (int i=0 ; i<3 ; i++){
			invM[i] = new double[3];
		}
		inv(M,invM);		// Inverse matrix
		
		cout<<"\ninv(M)=\n";      print33(invM);
		
		for (int i=0 ; i<=2 ; i++){
			A[i] = invM[0][i];
			B[i] = invM[1][i];
			C[i] = invM[2][i];
		}
		
		xc=(r[0]+s[0]+t[0])/3.0;
		yc=(r[1]+s[1]+t[1])/3.0;			// Position of centroid
		
		Qav = Q(xc,yc);			// Average value of Q
		Fav = F(xc,yc);			// Average value of F
		
		cout << "\nQ_{av}= " << Qav << "\nF_{av}= " << Fav<<"\n";
		
		// CONSTRUCT MATRIX K AND VECTOR b:
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
		cout << "\nK=\n"; 		print(K);
		cout << "\nb=\n"; 		print_vec(b);
		cout << "\n";
		
		// Save the obtained K and b into the system with all K's and b's:
		for (int i=0 ; i<=2 ; i++){
			for (int j=0 ; j<=2 ; j++){
				pK[p][i][j] = K[i][j];
			}
			pb[p][i] = b[i];
		}
	}
	
	// ASSEMBLE THE SYSTEM:
	
	double sK[n][n], 		// System matrix
			sb[n]; 			// System vector
	// we always go counterclockwise around the element in selecting the nodes.
	
	for (int i=0 ; i<=n-1 ; i++){
		for (int j=0 ; j<=n-1 ; j++){
			sK[i][j] = 0.0;
			for (int p=0 ; p<=m-1 ; p++){
				// cout << "ELEMENT " << p+1 << ":\n";
				for (int q=0 ; q<=2 ; q++){
					// cout << o[p][q] << "\n";
					if (o[p][q] == i+1){
						for (int r=0 ; r<=2 ; r++){
							if (o[p][r] == j+1){
								// cout << "K[" << p+1 << "][" << q+1 << "][" << r+1 <<"]\n";
								sK[i][j] = sK[i][j] + pK[p][q][r];
							}
						}
					}
				}
			// cout <<"\n";
			}
		}
	}
	for (int i=0 ; i<=n-1 ; i++){
		sb[i] = 0.0;
		for (int p=0 ; p<=m-1 ; p++){
			for (int q=0 ; q<=2 ; q++){
				if (o[p][q] == i+1){
					sb[i] = sb[i] + pb[p][q];	
				}
			}
		}
	}
	
	cout <<"\nSystem Matrix = \n";
	printn(sK,n);	cout<<"\n";
	print_vecn(sb,n);
	
	return 0;
}
