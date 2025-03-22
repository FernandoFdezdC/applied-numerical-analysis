// FINITE ELEMENTS FOR PARTIAL-DIFFERENTIAL EQUATIONS

// Date: 5/7/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n = 4; 			// number of nodes
const int m = 2;			// number of elements
const int s = 4;			// number of sides on the edge


void solve(double M[n][n], double b[n], double *x, int n){
	// M[n][n], Matrix M such that Mx=b
    // x[n]: solution
	
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
}

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
void printn_2(double M[n-2][n-2], int n){			// Displays nxn matrix M
	for (int i=0 ; i<=n-3 ; i++){
		for (int j=0 ; j<=n-3 ; j++){
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
void print_vecn_2(double b[n-2], int n){		// Displays nxn vector b
	for (int i=0 ; i<=n-3 ; i++){
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
	return -y/10.0;
}

double F(double x, double y){		// Function F(x,y)
	return x/4.0 + y -12.0;
}

double alpha1(double x, double y){		// Function alpha(x,y)
	return 0.1;
} double beta1(double x, double y){		// Function beta(x,y)
	return 1.0;
} double alpha2(double x, double y){		// Function alpha(x,y)
	return 0.0;
} double beta2(double x, double y){		// Function beta(x,y)
	return 0.0;
} double alpha4(double x, double y){		// Function alpha(x,y)
	return 0.0;
} double beta4(double x, double y){		// Function beta(x,y)
	return 2.0;
}

void midpoint(double v[s][2], double **ms, int s){
	for (int i=0 ; i<=s-2 ; i++){
		ms[i][0] = (v[i+1][0] + v[i][0])/2.0;
		ms[i][1] = (v[i+1][1] + v[i][1])/2.0;
	}
	ms[s-1][0] = (v[0][0] + v[s-1][0])/2.0;
	ms[s-1][1] = (v[0][1] + v[s-1][1])/2.0;
}

void distance(double v[s][2], double *length, int s){		// distance between succesive vectors at the edge
	for (int i=0 ; i<=s-2 ; i++){
//		cout << v[i][1] - v[i][1] << "\n";
		length[i] = pow(pow(v[i+1][1]-v[i][1],2.0) + pow(v[i+1][0]-v[i][0],2.0),1.0/2.0);
//		cout << length[i]<<"\n";
	}
	length[s-1] = pow(pow(v[0][1]-v[s-1][1],2.0) + pow(v[0][0]-v[s-1][0],2.0),1.0/2.0);;
}

int main(){
	cout << setprecision(9);
	
	// Elliptic PDE to solve:
	// $u_{xx}+u_{yy}+Q(x,y)u=F(x,y)$
	// on region $R$ with boundary conditions $u(x,y)=u_0$ on $L_1$,
	// $\frac{\partial u}{\partial n}=\alpha u + \beta$ on $L_2$
	// Approximate solution: v(x,y)
	
	double pr[m][2], ps[m][2], pt[m][2];	// position of local points for each element
	// Element 1:
	pr[0][0] = 0.0;		pr[0][1] = 0.0;		// Global point 2
	ps[0][0] = 5.0;		ps[0][1] = 3.0;		// Global point 4
	pt[0][0] = -1.0;	pt[0][1] = 4.0;		// Global point 1
	// Element 2:
	pr[1][0] = 0.0;		pr[1][1] = 0.0;		// Global point 2
	ps[1][0] = 4.0;		ps[1][1] = 0.0;		// Global point 3
	pt[1][0] = 5.0;		pt[1][1] = 3.0;		// Global point 4
	
	
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
		
		xc=(pr[p][0]+ps[p][0]+pt[p][0])/3.0;
		yc=(pr[p][1]+ps[p][1]+pt[p][1])/3.0;			// Position of centroid
		
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
	
	int o[m][3]; 			// Order of the points in each element
							// (connectivity relations)
	// Element 1:
	o[0][0] = 2; 	o[0][1] = 4; 	o[0][2] = 1;
	// Element 2:
	o[1][0] = 2; 	o[1][1] = 3; 	o[1][2] = 4;
	
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
	printn(sK,n);	cout<<"b=\n";
	print_vecn(sb,n);
	
	
	// APPLY NON-DIRICHLET BOUNDARY CONDITIONS:
	
	// There are four different sides with vertices:
	double v[s][2];			// position of each point on the edge
	v[0][0]=-1.0; v[0][1]=4.0; 		// Position of point 1
	v[1][0]=0.0; v[1][1]=0.0; 		// Position of point 2
	v[2][0]=4.0; v[2][1]=0.0; 		// Position of point 3
	v[3][0]=5.0; v[3][1]=3.0; 		// Position of point 4
	
	double **ms;
	ms = new double *[s];
	for (int i=0 ; i<=s-1 ; i++){
		ms[i] = new double[2];
	}
	midpoint(v, ms, s);
	cout << "\nMidpoint of side 1 = (" << ms[0][0]<<", "<<ms[0][1]<<")\n";
	cout << "Midpoint of side 2 = (" << ms[1][0]<<", "<<ms[1][1]<<")\n";
	cout << "Midpoint of side 4 = (" << ms[3][0]<<", "<<ms[3][1]<<")\n";
	
	// Length of sides with non-Dirichlet boundary conditions:
	double *L = new double[s];
	distance(v,L,s);
	
	// L[0]: side 1. The other sides are named after it in counterclockwise order.
	
	cout << "\nl41= " << L[3] << "\nl12= " << L[0] << "\nl23= " << L[1];
	
	double alphaav[s], betaav[s];
	alphaav[0] = alpha1(ms[0][0], ms[0][1]);
	betaav[0] = beta1(ms[0][0], ms[0][1]);		// Side 1
	alphaav[1] = alpha2(ms[1][0], ms[1][1]);
	betaav[1] = beta2(ms[1][0], ms[1][1]);		// Side 2
	alphaav[3] = alpha4(ms[3][0], ms[3][1]);
	betaav[3] = beta4(ms[3][0], ms[3][1]);		// Side 4
	
	// ROW 1:
	sK[0][0] = sK[0][0] - alphaav[0]*L[0]/3.0 - alphaav[3]*L[3]/3.0;
	sK[0][1] = sK[0][1] - alphaav[0]*L[0]/6.0;
	sK[0][3] = sK[0][3] - alphaav[3]*L[3]/6.0;
	sb[0] = sb[0] + betaav[0]*L[0]/2.0 + betaav[3]*L[3]/2.0;
	// ROW 2:
	sK[1][0] = sK[1][0] - alphaav[0]*L[0]/6.0;
	sK[1][1] = sK[1][1] - alphaav[0]*L[0]/3.0 - alphaav[1]*L[1]/3.0;
	sK[1][2] = sK[1][2] - alphaav[1]*L[1]/6.0;
	sb[1] = sb[1] + betaav[0]*L[0]/2.0 + betaav[1]*L[1]/2.0;
	// ROW 3:
	sK[2][1] = sK[2][1] - alphaav[1]*L[1]/6.0;
	sK[2][2] = sK[2][2] - alphaav[1]*L[1]/3.0;
	sb[2] = sb[2] + betaav[1]*L[1]/2.0;
	// ROW 4:
	sK[3][0] = sK[3][0] - alphaav[3]*L[3]/6.0;
	sK[3][3] = sK[3][3] - alphaav[3]*L[3]/3.0;
	sb[3] = sb[3] + betaav[3]*L[3]/2.0;
	
	cout <<"\n\nNew system Matrix = \n";
	printn(sK,n);	cout<<"b=\n";
	print_vecn(sb,n);
	
	// APPLY DIRICHLET BOUNDARY CONDITIONS:
	double u_3 = 100.0, 		u_4 = 100.0; 		// Known values
	
	sb[2] = u_3, sb[3] = u_4;
	sK[2][0] = 0.0, sK[2][1] = 0.0, sK[2][2] = 1.0, sK[2][3] = 0.0;
	sK[3][0] = 0.0, sK[3][1] = 0.0, sK[3][2] = 0.0, sK[3][3] = 1.0;
	sb[0] = sb[0] - sK[0][2]*u_3 - sK[0][3]*u_4;
	sb[1] = sb[1] - sK[1][2]*u_3 - sK[1][3]*u_4;
	sK[0][2] = 0.0, sK[0][3] = 0.0;
	sK[1][2] = 0.0, sK[1][3] = 0.0;
	
	cout <<"\n\nNew system Matrix = \n";
	printn(sK,n);	cout<<"b=\n";
	print_vecn(sb,n);
	
	// Obtain solution of the system of equations:
	double* u = new double[n];
	solve(sK,sb,u,n);
	
	cout << "\nTHE SOLUTION IS:\nu = [";
	for (int i=0 ; i<=n-2 ; i++){
	 	cout << u[i] <<", ";}
	cout << u[n-1] <<"];";
	
	return 0;
}
