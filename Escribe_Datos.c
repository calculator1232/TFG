#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <time.h> //semilla cambia con el tiempo

#define hbar 1
#define me 1
#define eps 1e-06

int main(){
    int Niter=10000;
    int N,M;
    double E,V,theta,deltaL,L;
    FILE *datos;
    datos=fopen("Datos.txt","w");
    N=1;
  //  M=70000;//Para log
    M=50000;
   // M=500;
    E=160;
    V=50;//EN meV
 // E=50.01;
 //  V=50;//EN meV
   // L=50;//en nm porque las k van en nm^-1
   // L=1500;lnT=-1.1 con alpha=0.5
   // L=26000; //G=0.15 con alpha=0.69
   // L=90;//T=0.46 en standar)
  //  L=725;//T=0.04
  L=2000;
    //L=2000 para -lnT=-1.1
    deltaL=0.02;
    theta=0.0;
    fprintf(datos, "%i\t %i\t %lf\t %lf\t %lf\t %lf\n",N,M,E,V,L,theta);
    printf( "%i\t %i\t %lf\t %lf\t %lf\t %lf\n",N,M,E,V,L,theta);

    fclose(datos);
return 0;
}
