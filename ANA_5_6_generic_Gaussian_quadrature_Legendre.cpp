// GAUSS-LEGENDRE QUADRATURE WITH LEGENDRE POLYNOMIALS
// APPLICATION TO THE SIMPLE PENDULUM PROBLEM

// Date: 7/1/2021
// Author: Fernando Fernández del Cerro

// #include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int n = 5;			 // Degree of the last polynomial calculated
const int N = 20;			 // Intervals for application of Newton method
const double theta_m=70;			// Angular amplitude of the pendulum (º)

void print(double M[n+1][n+1], int n){			// Displays matrix M
	for (int i=0 ; i<=n ; i++){
		for (int j=0 ; j<=n ; j++){
			std::cout << M[i][j] << " ";
		}
		std::cout << "\n";
	}
}

double f(double x){					// Function to integrate
	double f = pow(1-pow(sin(theta_m*M_PI/360)*sin(x),2),-0.5);
	return f;
}

double L(double a[n+1][n+1], double x, int m){	// Value of the degree m Legendre polynomial in x
	double L = a[m][0];
	for (int i=1 ; i<=m ; i++){
		L = L + a[m][i]*pow(x,i);
	}
  	return L;
}

double dLdx(double a[n+1][n+1], double x, int m){	// Derivative of the degree m Legendre polynomial in x
	double dL = a[m][1];
	for (int i=2 ; i<=m ; i++){
		dL = dL + i*a[m][i]*pow(x,i-1);
	}
  	return dL;
}

double* newton(double *newx, double oldx[N+1], double p[n+1][n+1], double TOL, int N, int m){		// Newton-Raphson process
	// oldx: starting approximation to the root of pn(x)
	// TOL: precision needed < 1

	double diff;	// difference between one approximation and the following one
	double c;		// count of iterations

	for (int j=0 ; j<=N ; j++){
		newx[j] = oldx[j];
		cout << "First value: " << oldx[j] << ". ";
		c = 0, diff = 1;
		do {
			newx[j] = oldx[j] - L(p,oldx[j],m)/dLdx(p,oldx[j],m);
			diff = fabs(newx[j] - oldx[j]);
			oldx[j] = newx[j];
			c++;
		} while (fabs(diff)>TOL||fabs(L(p,newx[j],m))>TOL);
		cout << "Root = " << newx[j] << ". Necessary iterations: " << c << "\n";
	}
	return newx;					// Best approximations to the roots
}

void* sort(double *list, int n){	// Sorts list of n numbers
	double temp;
	for(int i=0 ; i<=n-1 ; i++){
		for(int j=i+1 ; j<=n-1 ; j++)
		{
			if(list[i] < list[j])
			{		// Interchange numbers that are not sorted:
				temp = list[i];
				list[i] = list[j];
				list[j] = temp;
			}
		}
	}
}

int main(){
	cout << setprecision(11);

	// LEGENDRE POLYNOMIALS:
	// L_0(x)=1, 		L_1(x)=x
	double p[n+1][n+1];		// Coefficients of the Legendre polynomials
							// First entry: degree of polynomial.
							// Second entry: coefficients of the polynomial
	p[0][0] = 1.0;
	p[1][0] = 0.0;
	p[1][1] = 1.0;
	// Recurrence relation:
	for(int i=0 ; i<=1 ; i++){
		for(int j=i+1 ; j<=n+1 ; j++){
			p[i][j] = 0;
		}
	}
	for(int i=2 ; i<=n ; i++){
		p[i][0] = -(i-1)*(p[i-2][0])/i;
		for(int j=1 ; j<=n ; j++){
			p[i][j] = ((2*i-1)*p[i-1][j-1]-(i-1)*p[i-2][j])/i;
		}
	}
	cout << "Coefficients of Legendre polyonomial up to "<<n<<" degree:\n";
	print(p,n);
	// LEGENDRE POLYNOMIALS. END.

	// Equation to solve: L_m(x) = 0
	// L_m(x) = p[m][0]+p[m][1]*x+...+p[m][m]*x^m
	std::cout << "L_m(x) = p[m][0]+p[m][1]*x+...+p[m][m]*x^m" << endl;

	// Open a file:
	// ofstream file;
	// file.open ("solution.txt");

	// ROOTS OF THE LEGENDRE POLYNOMIALS:
	double TOL = 1.0e-11;		 // precision of the result
	double a = -1.0, b = 1.0;		 // interval where roots are searched


	double h = (b-a)/N;	 		 // difference between starting values for the approximate root x
	double oldx[N+1], newx[N+1]; // values of the approximate roots

	double weight[N+1];

	for (int m=n ; m<=n ; m++){
		oldx[0] = a;
		for (int j=1 ; j<=N ; j++){
			oldx[j] = oldx[j-1] + h;
			// 	cout << oldx[j] << "\n";
		}
		cout << "\nRoots of " << m << " degree Legendre polynomial:\n";
		// Apply Newton-Raphson process to oldx to the m degree Legendre polynomial:
		double* root = newton(newx, oldx, p, TOL, N, m);
		sort(root, N+1);		// Sort the found roots from higher to lower
		cout <<"\nSorted list of roots:\n";
		for (int j=0 ; j<=N ; j++){
			cout << root[j]<< "\n";
		}
		cout <<"\n";

		// GAUSSIAN QUADRATURE:
		double I = 0;				// Integral of f(x) from a to b
		double x0 = 0;				// Lower value of the integral
		double x1 = M_PI/2;			// Upper value of the integral

		cout << "GAUSSIAN QUADRATURE WITH N = " << n << " DIVISIONS:\n";
		cout << "f(x) = \frac{1}{\sqrt{1-\left(sin{\frac{\theta_m}{2}}sin{x}\right)^2}}\n";
		// Points and weights:
		double x[m], w[m];
		cout <<"Every root without repetition:\n";
		x[0] = root[0];		cout<<x[0]<<"\n";
		int s=0;		// entry of the vector with non-repeated roots
		for (int j=0 ; j<=N ; j++){
			for(int i=0 ; i<=s ; i++){
				if(fabs(root[j]-x[s])>10*TOL){
					x[s+1] = root[j];
					s++;
					cout << x[s]<<"\n";
				}
			}
		}
		cout << "\nNumber of different roots="<<s+1<<"\n\n";

		for (int j=0 ; j<=s ; j++){
			w[j] = 2/((1-pow(x[j],2))*pow(dLdx(p,x[j],m),2));	// gaussian weights
		}
		cout << "Points -- Weights:\n";

		// Calculate integral I:
		for(int i=0 ; i<=m-1 ; i++){
			I = I + w[i]*f((x1-x0)*x[i]/2+(x0+x1)/2);
			cout << x[i] << " -- " << w[i] << "\n";
		}

		I = (x1-x0)*I/2;

		cout << "\nGAUSSIAN QUADRATURE:\n";
		cout << "Integral of f(x) from " << x0 << " to " << x1 << " = " << I;
		cout << "\n\frac{T}{T_0}=" << 2*I/M_PI;
		// GAUSSIAN QUADRATURE. END.
	}

	// Write the solution in a file:
	// file << x[j][1] << " is a solution to the equation f(x) = 0 if x_0 = " << x[j][0] << ". Necessary steps: " << x[j][2] << endl;
	// For LaTeX:
	// file << x[j][0] << " & " << x[j][1] << " & " << x[j][2] << " \\\\" << " \\hline"<<endl;
	// file.close();

	return 0;
}
