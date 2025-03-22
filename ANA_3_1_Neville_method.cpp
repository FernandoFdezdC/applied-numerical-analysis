// NEVILLE'S METHOD FOR INTERPOLATION

// Date: 8/3/2021
// Author: Fernando Fernández del Cerro

#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

const int n = 5; 			// amount of data to be interpolated

void* sort(double *list, double *a, double *b, int n){	// Sorts list of n numbers
	// but lists a and b are changed to fit the correct order of list
	double temp, temp_a, temp_b;
	for(int i=0 ; i<=n-1 ; i++){		
		for(int j=i+1 ; j<=n-1 ; j++)
		{
			if(list[i] > list[j])
			{		// Interchange numbers that are not sorted:
				temp = list[i];
				list[i] = list[j];
				list[j] = temp;
				temp_a = a[i];
				a[i] = a[j];
				a[j] = temp_a;
				temp_b = b[i];
				b[i] = b[j];
				b[j] = temp_b;
			}
		}
	}
}

int main(){
	cout << setprecision(9);
	
	double x[n]={10.1, 22.2, 32.0, 41.6, 50.5},		// data to be interpolated
	f[n] = {0.17537, 0.37784, 0.52992, 0.66393, 0.63608};		// data to be interpolated
	double x0 = 27.5; 			// point to calculate best f(x_0)
	
	double diff[n];			// List of numbers |x_0-x_i|
	for (int i=0 ; i<=n-1 ; i++){
		diff[i]=fabs(x0-x[i]);
	}
	sort(diff,x,f,n);
	
	double P[n][n];
	for (int i=0 ; i<=n-1 ; i++){
		P[i][0]=f[i];		// first column of Neville table
	}
	// Neville's method:
	for (int j=1 ; j<=n-1 ; j++){
		for (int i=0 ; i<=n-j-1 ; i++){
			P[i][j]=((x0-x[i])*P[i+1][j-1]+(x[i+j]-x0)*P[i][j-1])/(x[i+j]-x[i]);		// set the other columns equal to 0
		}
		for (int i=n-j ; i<=n-1 ; i++){
			P[i][j]=0;		// set the other elements equal to 0
		}
	}
	
	cout << "Neville table:\n";
	for (int i=0 ; i<=n-1 ; i++){
		cout << i <<" "<<x[i]<<" ";
		for (int j=0 ; j<=n-1 ; j++){
			std::cout << P[i][j] << " ";
		}
		std::cout << "\n";
	}
	
	return 0;
}
