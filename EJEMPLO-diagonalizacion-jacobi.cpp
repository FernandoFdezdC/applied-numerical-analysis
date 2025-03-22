#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
using namespace std;
const int n=2;
const double tolerancia=1e-7;
void jacobi (double A[n][n], double U[n][n], double L[n]);
void muestramatriz(double M[n][n], int f, int c);
void muestravector(double V[n], int d);

int main()
{
//double A[n][n]={8,-1,3,-1,-1,6,2,0,3,2,9,1,-1,0,1,7};
double U[n][n], L[n];
char archivo_entrada[30];
double A[n][n];
    // char archivo_salida[30];
cout << "Introduzca nombre del archivo de entrada: ";
cin >> archivo_entrada;
// Crear un flujo de entrada y llamarlo fentrada
ifstream fentrada(archivo_entrada, ios::in);
for (int i=0; i<n; i++){
  for (int j=0; j<n; j++){  
fentrada >> A[i][j];}
}
fentrada.close(); // cerrar el archivo

jacobi(A,U,L);
// En L estan los autovalores y en U los autovectores     
cout << "Autovalores" << endl;
muestravector(L,n);
//cout << "Matriz A" << endl;
//muestramatriz(A,n,n);
//cout << endl << "Matriz U" << endl;
cout << "Autovectores" << endl;
muestramatriz(U,n,n);
system("PAUSE");
return 0;
}

void jacobi (double A[n][n], double U[n][n], double L[n])
{
double theta, si, co, uki, aii, max;
int i, j, iter=0;
// matriz U -> I
for (int i=0; i<n ;i++) {
    for (int j=0; j<n; j++)
    { if (i==j) { U[i][j]=1.0;}
      else { U[i][j]=0.0;}
    }
    }
do 
{
// Buscamos el elemento maximo por encima de la diagonal
max=0;//fabs(A[0][1]);i=0, j=1;
for (int ib=0; ib < n ;ib++) {
    for (int jb=ib+1; jb < n; jb++)
    { if (fabs(A[ib][jb]) >= max) 
      { max=fabs(A[ib][jb]), i=ib, j=jb;}
      }
     }
// llevamos a cabo la rotacion Rij
if (A[i][i]==A[j][j])
   { theta = M_PI/4;}   // pi/4
else
   { theta = 0.5*atan(2.0*A[i][j]/(A[i][i]-A[j][j]));}
si=sin(theta);
co=cos(theta);
for (int k=0; k<n; k++)
{ 
uki=U[k][i];
U[k][i]=uki*co + U[k][j]*si;
U[k][j]=-uki*si + U[k][j]*co;
if (k!=i && k!=j)
{
 A[k][i]=A[i][k]*co + A[j][k]*si;
 A[k][j]=-A[i][k]*si + A[j][k]*co;
 A[j][k]=A[k][j];
 A[i][k]=A[k][i];
}
}  // cerramos bucle for en k
aii=A[i][i];
A[i][i]=aii*pow(co,2.0) + A[j][j]*pow(si,2.0) + 2.0*A[i][j]*si*co;    
A[j][j]=aii*pow(si,2.0) + A[j][j]*pow(co,2.0) - 2.0*A[i][j]*si*co;    
A[i][j]=0.0; // Se impone que sea cero despues de la rotacion
A[j][i]=0.0; // Se impone que sea cero despues de la rotacion
iter=iter+1;
} while (max > tolerancia);     
/* al final U es la matriz producto de todos las rotaciones y sus columnas 
seran los autovectores, y los elementos diagonales de A son los autovalores que
pasamos al vector L */
cout << "Numero de iteraciones=" << iter << endl << endl;
for (int k=0; k < n; k++) { L[k]=A[k][k];}
}

void muestramatriz(double M[n][n], int f, int c)
{
for (int i=0; i<f; i++) {
    for (int j=0; j<c; j++) {
    cout << M[i][j] << " ";
     }
     cout << endl;
 }
}

void muestravector(double V[n], int d)
{
for (int k=0; k < d; k++) 
{ cout << V[k] << endl; }
  cout << endl;
}

