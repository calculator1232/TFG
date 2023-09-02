#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <time.h> //semilla cambia con el tiempo

#define SACAR_PROBABILIDADES
#define hbar 1
#define me 1
#define eps 1e-06
#define numL 25
#define numE 70
#define CAMBIAR_L
//#define CAMBIAR_E
//****HISTOGRAMA******
#define N_Intervalos 50
#define N_Datos_Max 1000000

double PI=acos(-1.0);
double alpha_Levy=0.5;
typedef struct{
    double complex r1, t, r2;
} S_matrix;//Lo usamos para definir las matrices S, con 3 componentes complejas r1,t y r2

void readata (int *,int *,double *, double *, double *, double *, char *);//Lee los datos del fichero char *
void creaplt(char *);//Crea a partir de un fichero, un plt
void creaplt_proms(char *,double );//Crea a partir de un fichero, un plt de los promedios de T
void creaplt_proms_log(char *,double );//Crea a partir de un fichero, un plt de los promedios de log T
void creaplt_fluctuaciones(char *,double );//Crea a partir de un fichero, un plt de los promedios de log T
void creaplt_probabilidades_log(char *,double,double ,char *,double ,double);//Para ver la distribución de probabilidad para 2 valores dados de <-ln T>
void creaplt_probabilidades_T(char *,double ,double, char *,double, double);//Para ver la distribución de probabilidad para 2 valores dados de <T> y con los dos <-ln T> para esa L
void creaplt_probabilidades_T_paraesoslog(char *,double ,char *,double );//Para ver la distribución de probabilidad para los 2 valores  de <T> dados para los dos <-ln T> anteriores
void creaplt_probabilidades_log_paraesosT(char *,double ,double,char *,double,double );//Para ver la distribución de probabilidad   para los 2 valores  de <-ln T> dados para los dos <T> anteriores
S_matrix compose(S_matrix, S_matrix);//Compone 2 matrices
S_matrix translate (S_matrix, double, double);//Traslada una matriz S una distancia d, se usa para considerar la propagación de la onda
S_matrix barrier(double , double , double );//Devuelve la matriz S_matrix correspondiente a una barrera de anchura a y altura V
int check(S_matrix);//Comprueba que S*S^T =I
void Histograma(double *, double *,int ,int ,double *,double *, double *);

///Generador de números reales aleatorios de Parisi-Rapuano
#define NormRANu (2.3283063671E-10F)
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;
extern float random(void); //al final del código se desarrollan
extern void ini_ran(int);

double random_Levy_general(double);

int main (){
   // char argv[]="Datos.txt";
    ini_ran(time(NULL));//La semilla cambia con la hora
    int i,j,N,M;//Número total de barreras
    double escalamiento,
              k,    //Componente x del vector de ondas en el vacío
              k_y,  //Componente y del vector de ondas
              E,    //Energía de la onda
              E_0,  //Energía inicial
              V,    //Altura de las barreras
              a,    // anchura aleatoria de las barreras
              d,    // distancia entre barreras aleatoria
              L,    //Longitud total del sistema (va variando)
              L_0,  //Longitud inicial del sistema
              deltaL, //Incremento en la longitud del sistema para ver la dependencia con esta
              deltaE,  //Incremento en la energía incidente para ver la dependencia con esta
              theta_i,//ÁNGULO DE INCIDENCIA EN EL SISTEMA
              theta_j,//ÁNGULO INCIDENTE EN LA BARRERA J, NO LO CONSIDERAMOS DE MOMENTO
              theta_e,//ÁNGULO DE SALIDA
              T,    //TRANSMITANCIA
              fluctuacion, //Dispersión de T
              contador_log_transm,contador_transm,contador_transmcuadrada,promediolog,logsmall,logbig,promedioT,promedioTcuadrada,Tsmall,Tbig,Tparalogsmall,TparalogBig,logparaTBig,logparaTSmall,
              minimo,maximo,delta,//Valores para el histograma
              LparalogSmallbuscado,LparalogBigbuscado,LparaTSmallBuscado,LparaTbigbuscado,//Valores de longitud del sistema para las T y log T en los que se quiere incidir
              LogparaTBig,LogparaTSmall,
              promlogTsmall,promlogTbig,promTbig,promTsmall;//VALORES ENTORNO A LOS QUE BUSCAR PARA SACAR LOS GRÁFICOS DE P(-log T) y P(T)
    S_matrix S;//Sistema total
    FILE *datos;//archivo de salida
    FILE *graf;//PLT DE SALIDA
    FILE *promedios,*promedios_log,*FLUCTUACIONES;//promedios para misma L
    FILE *histlogsmall,*histTparalogsmall,*histlogbig,*histTparalogBig,*histlogparaTSmall,*histlogparaTBig;//2 columnas de tamaño M con los logaritmos guardados
     FILE *histTsmall,*histTbig;//2 columnas de tamaño M con los logaritmos guardados
     FILE *transmisiones_TBig,*transmisiones_TSmall,*logs_para_TSmall,*logs_para_TBig;//Anota las transmisiones y logs  para cuando el <T> es el buscado
    FILE *logs_Big,*logs_Small,*T_para_logSmall,*T_para_logBig;//Anota las transmisiones y logs  para cuando el <log T> es el buscado
    FILE *numBars;
//FIN DE DECLARACIONES INICIALES
   readata(&N, &M, &E, &V, &L,&theta_i, "Datos.txt");
    printf("Semilla=%d \n eps=%lf \n Datos:M=%i\t E=%lf\t V=%lf\t L=%lf\n",time(NULL),eps,M,E,V,L);
         double logaritmos[M],transmisiones[M], Hist[N_Intervalos],barreras[M];//Para almacenar log T y T para un solo L y los histogramas que se plotearán
    k_y=sqrt(2*me*E)*sin(theta_i)/hbar;
    k=sqrt(2*me*E-hbar*hbar*k_y*k_y)/(hbar);
     printf("k=%lf\n",k);
   // deltaL=100;
    L_0=L;
    E_0=E;
   escalamiento=1;
   //deltaL=0.2;
   deltaL=10000;
   deltaE=1;
    promlogTsmall=-1.1;
     promlogTbig=-6.37;
     promTsmall=0.15;
      promTbig=0.46;
    datos=fopen("datos_SCH.dat","w");
    promedios=fopen("promedios_SCH.dat","w");
    promedios_log=fopen("promedios_log_SCH.dat","w");
    FLUCTUACIONES=fopen("fluctuaciones_SCH.dat","w");
       printf("Entra al bucle \n");
#ifdef CAMBIAR_L
    for(int iter=0;iter<numL;iter++){
                    contador_transm=contador_log_transm=0.0;
        for(j=0;j<M;j++)
            {
            d=0;
           a=escalamiento*random_Levy_general(alpha_Levy);//da la anchura de la primera barrera
              // printf("\n a vale= %lf, L*log(random())=%lf \n",a , L*log(random()));
            S=barrier(a,V,k);//Construye la matriz S
            N=1;
       //   printf("S inicial r=%lf + i*%lf \t t=%lf + i*%lf\t r2=%lf+ i*%lf\n",creal(S.r1),cimag(S.r1), creal(S.t),cimag(S.t),creal(S.r2),cimag(S.r2));
            //  printf("\n S construida \n");
           while(1>0)
            {
                d+=a+escalamiento*random_Levy_general(alpha_Levy);
                a=escalamiento*random_Levy_general(alpha_Levy);//Barrera nueva
              //  printf("a=%lf y log random=%lf\n",a, -a*2*N/L);
            //  printf("r=%lf \t t=%lf \t r2=%lf \n",S.r1, S.t, S.r2);
              if((d+a)>L){
                    /*
                    Si se supera la longitud del sistema, en vez de componer S con la nueva barrera, se rellena con vacío el espacio entre d y L
                    y salimos del bucle para ver si es unitaria
                    */
                    S=compose(S,translate(barrier(L-d,0.0,k),L,k));
                    N++;
                    break;}
                S=compose(S,translate(barrier(a,V,k),d,k));//Se compone la matriz S con la nueva barrera desplazada a la distancia d
                N++;
                //printf("T=%lf\n",T);
            }

            if(check(S)==0)//Comprueba que sea unitaria S
            {
           // printf("S unitaria, se escribe\n");
            T=cabs(S.t*conj(S.t));
            ///Calculamos <T> y <-log T>
            /** Para resolver problemas con T=0->log T divergente, ponemos que si T<10^-6 se toma -log (10^-6)*/
            if (T<eps) fprintf(datos,"%lf \t %lf\n",T,k*L,-log(eps));
            else fprintf(datos,"%lf \t %lf\n",T,k*L,-log(T));
             contador_transm+=T;
            if (T<eps) contador_log_transm-=log(eps);
            else contador_log_transm-=log(T);
            #ifdef SACAR_PROBABILIDADES //Apuntamos los valores de los logaritmos y transmisiones para cada longitud L
             if (T<eps)logaritmos[j]=log(eps);
             else logaritmos[j]=log(T);
            transmisiones[j]=T;
            #endif
            barreras[j]=N;
            }
            else {printf("S no es unitaria");exit(1);}
        }

       // printf(" %i \t L vale %lf\n",iter,L);
       promediolog=-contador_log_transm/((double)M);
       promedioT=contador_transm/((double)M);
        fprintf(promedios_log,"%lf \t %lf\n",k*L,promediolog);
        fprintf(promedios,"%lf \t %lf\n",k*L,promedioT);
    #ifdef SACAR_PROBABILIDADES
            //********************HISTOGRAMAS********************//

                if(promediolog>=(promlogTsmall-0.0010) && promediolog<=(promlogTsmall+0.0010)){
                        histlogsmall=fopen("LogsSmall.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(logaritmos,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                           fprintf(histlogsmall,"%lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                     //  fprintf(histlogsmall,"%lf \t %lf\n",i*delta+minimo,Hist[i]);
                                  // for(int x=0;x<M;x++)fprintf(histlogsmall,"%lf \n", logaritmos[x]);
                                    logsmall=promediolog;
                                       // printf(" Histograma logaritmo dibujado\n");
                                        fclose(histlogsmall);
                        logs_Small=fopen("logs_para_logSmall.txt","w");
                            for (i=0;i<M;i++) fprintf(logs_Small," %lf\n",logaritmos[i]);
                            fclose(logs_Small);
                    //*************HISTOGRAMA DE P( T) RESPECTO A  T para el <-log T> buscado********************
                        histTparalogsmall=fopen("TparalogSmall.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(transmisiones,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histTparalogsmall," %lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    Tparalogsmall=promedioT;
                                //     printf(" Histograma promedio dibujado\n",iter,L);
                                fclose(histTparalogsmall);
                            T_para_logSmall=fopen("Transmisiones_para_logSmall.txt","w");
                            for (i=0;i<M;i++) fprintf(T_para_logSmall," %lf\n",transmisiones[i]);
                            fclose(T_para_logSmall);
                            numBars=fopen("Numero_barreras.txt", "w");
                             for (i=0;i<M;i++) fprintf(numBars," %lf\n",barreras[i]);
                            fclose(numBars);
                                LparalogSmallbuscado=L;
                   }
            if(promediolog>=(promlogTbig-0.10) && promediolog<=(promlogTbig+0.10)){
                            histlogbig=fopen("LogsBig.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(logaritmos,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histlogbig,"%lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    logbig=promediolog;
                                        //printf(" Histograma logaritmo dibujado\n");
                                          fclose(histlogbig);
                            logs_Big=fopen("logs_para_logBig.txt","w");
                            for (i=0;i<M;i++) fprintf(logs_Big," %d \t %lf\n",i,logaritmos[i]);
                            fclose(logs_Big);
                        //*************HISTOGRAMA DE P(T) RESPECTO A T para el <-log T> buscado********************
                            histTparalogBig=fopen("TparalogBig.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(transmisiones,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histTparalogBig,"%lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    TparalogBig=promedioT;
                                //     printf(" Histograma promedio dibujado\n",iter,L);
                                fclose(histTparalogBig);
                                    T_para_logBig=fopen("Transmisiones_para_logBig.txt","w");
                            for (i=0;i<M;i++) fprintf(T_para_logBig," %d \t %lf\n",i,transmisiones[i]);
                            fclose(T_para_logBig);
                                LparalogBigbuscado=L;               }
            if(promedioT>=(promTsmall-0.0001) && promedioT<=(promTsmall+0.0001)){
                            histTsmall=fopen("TSmall.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(transmisiones,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histTsmall,"%lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    Tsmall=promedioT;
                                //     printf(" Histograma promedio dibujado\n",iter,L);
                                fclose(histTsmall);
                             transmisiones_TSmall=fopen("Transmisiones_TSmall.txt","w");
                            for (i=0;i<M;i++) fprintf(transmisiones_TSmall," %d \t %lf\n",i,transmisiones[i]);
                            fclose(transmisiones_TSmall);
                        //*************HISTOGRAMA DE P(-log T) RESPECTO A -log T para el <T> buscado********************
                            histlogparaTSmall=fopen("LogparaTSmall.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(logaritmos,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histlogparaTSmall,"%lf \t %lf \n",i*delta+minimo+delta/2,Hist[i]);
                                    logparaTSmall=promediolog;
                                      logs_para_TSmall=fopen("Logs_para_TSmall.txt","w");
                            for (i=0;i<M;i++) fprintf(logs_para_TSmall," %d \t %lf\n",i,-logaritmos[i]);
                            fclose(logs_para_TSmall);
                                //     printf(" Histograma promedio dibujado\n",iter,L);
                                fclose(histlogparaTSmall);
                                LparaTSmallBuscado=L;
                                LogparaTSmall;
                   }
            if(promedioT>=(promTbig-0.01) && promedioT<=(promTbig+0.01)){
                        histTbig=fopen("TBig.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(transmisiones,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histTbig," %lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    Tbig=promedioT;
                                 //   printf(" Histograma promedio dibujado\n",iter,L);
                                  fclose(histTbig);
                            transmisiones_TBig=fopen("Transmisiones_TBig.txt","w");
                            for (i=0;i<M;i++) fprintf(transmisiones_TBig," %d \t %lf\n",i,transmisiones[i]);
                            fclose(transmisiones_TBig);
                         //*************HISTOGRAMA DE P(-log T) RESPECTO A -log T para el <T> buscado********************
                        histlogparaTBig=fopen("logparaTBig.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(logaritmos,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histlogparaTBig," %lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    logparaTBig=promediolog;
                                //     printf(" Histograma promedio dibujado\n",iter,L);
                                fclose(histlogparaTBig);
                                LparaTbigbuscado=L;
                                LogparaTBig=promediolog;
           }
         #endif // SACAR_PROBABILIDADES
        L+=deltaL;//CAMBIAMOS LA LONGITUD DEL SISTEMA
        if ((int)((L-L_0)/deltaL)%(numL/10)==0) printf("Completado al %lf porciento\n",((double)(iter+1)/((double)(numL)))*100);
    }
    printf(" Sale de bucle de L\n");
        fclose(datos);
        fclose(promedios);
        fclose(promedios_log);
        printf(" cierra ficheros\n");
      //  printf("\n t=%lf + i*%lf \n r=%lf + i*%lf", creal(S.t), cimag(S.t),creal(S.r1), cimag(S.r1));
        creaplt("datos_SCH.dat");
        creaplt_proms("promedios_SCH.dat",E);
        creaplt_proms_log("promedios_log_SCH.dat",E);
    #ifdef SACAR_PROBABILIDADES
        creaplt_probabilidades_log("LogsBig.txt", logbig,LparalogBigbuscado,"LogsSmall.txt", logsmall,LparalogSmallbuscado);
        creaplt_probabilidades_log_paraesosT("LogparaTBig.txt", logparaTBig,LparaTbigbuscado,"LogparaTSmall.txt", logparaTSmall,LparaTSmallBuscado);
        creaplt_probabilidades_T("TBig.txt", Tbig,LogparaTBig,"TSmall.txt", Tsmall,LogparaTSmall);
         creaplt_probabilidades_T_paraesoslog("TparalogBig.txt", TparalogBig,"TparalogSmall.txt", Tparalogsmall);
    #endif // SACAR_PROBABILIDADES
        printf("\n done dependencia con L");
#endif // CAMBIAR_L
/**DEPENDENCIA CON LA ENERGÍA INCIDENTE**/
#ifdef CAMBIAR_E
L=L_0;
    for(int iter=0;iter<numE;iter++){
                    contador_transm=contador_log_transm=contador_transmcuadrada=0.0;
        for(j=0;j<M;j++)
            {
            d=0;
           a=escalamiento*random_Levy_general(alpha_Levy);//da la anchura de la primera barrera
              // printf("\n a vale= %lf, L*log(random())=%lf \n",a , L*log(random()));
            S=barrier(a,V,k);//Construye la matriz S
       //   printf("S inicial r=%lf + i*%lf \t t=%lf + i*%lf\t r2=%lf+ i*%lf\n",creal(S.r1),cimag(S.r1), creal(S.t),cimag(S.t),creal(S.r2),cimag(S.r2));
            //  printf("\n S construida \n");
           while(1>0)
            {
                d+=a+escalamiento*random_Levy_general(alpha_Levy);
                a=escalamiento*random_Levy_general(alpha_Levy);//Barrera nueva
              //  printf("a=%lf y log random=%lf\n",a, -a*2*N/L);
            //  printf("r=%lf \t t=%lf \t r2=%lf \n",S.r1, S.t, S.r2);
              if((d+a)>L){
                    /*
                    Si se supera la longitud del sistema, en vez de componer S con la nueva barrera, se rellena con vacío el espacio entre d y L
                    y salimos del bucle para ver si es unitaria
                    */
                    S=compose(S,translate(barrier(L-d,0.0,k),L,k));
                    break;}
                S=compose(S,translate(barrier(a,V,k),d,k));//Se compone la matriz S con la nueva barrera desplazada a la distancia d
                //printf("T=%lf\n",T);
            }

            if(check(S)==0)//Comprueba que sea unitaria S
            {
           // printf("S unitaria, se escribe\n");
            T=cabs(S.t*conj(S.t));
            ///Calculamos <T> y <-log T>
            /** Para resolver problemas con T=0->log T divergente, ponemos que si T<10^-6 se toma -log (10^-6)*/
            if (T<eps) fprintf(datos,"%lf \t %lf\n",T,L,-log(eps));
            else fprintf(datos,"%lf \t %lf\n",T,L,-log(T));
             contador_transm+=T;
             contador_transmcuadrada+=T*T;
            if (T<eps) contador_log_transm-=log(eps);
            else contador_log_transm-=log(T);
            #ifdef SACAR_PROBABILIDADES //Apuntamos los valores de los logaritmos y transmisiones para cada longitud L
             if (T<eps)logaritmos[j]=log(eps);
             else logaritmos[j]=log(T);
            transmisiones[j]=T;
            #endif
            }
            else {printf("S no es unitaria");exit(1);}
        }

       // printf(" %i \t L vale %lf\n",iter,L);
       promediolog=-contador_log_transm/((double)M);
       promedioT=contador_transm/((double)M);
       promedioTcuadrada=contador_transmcuadrada/((double)M);
       fluctuacion=sqrt(promedioTcuadrada-promedioT*promedioT);
        fprintf(promedios_log,"%lf \t %lf\n",E,promediolog);
        fprintf(promedios,"%lf \t %lf\n",E,promedioT);
         fprintf(FLUCTUACIONES,"%lf \t %lf\n",E,fluctuacion);
    #ifdef SACAR_PROBABILIDADES
            //********************HISTOGRAMAS********************//

                if(promediolog>=(promlogTsmall-0.0010) && promediolog<=(promlogTsmall+0.0010)){
                        histlogsmall=fopen("LogsSmall.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(logaritmos,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                           fprintf(histlogsmall,"%lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                     //  fprintf(histlogsmall,"%lf \t %lf\n",i*delta+minimo,Hist[i]);
                                  // for(int x=0;x<M;x++)fprintf(histlogsmall,"%lf \n", logaritmos[x]);
                                    logsmall=promediolog;
                                       // printf(" Histograma logaritmo dibujado\n");
                                        fclose(histlogsmall);
                    //*************HISTOGRAMA DE P( T) RESPECTO A  T para el <-log T> buscado********************
                        histTparalogsmall=fopen("TparalogSmall.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(transmisiones,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histTparalogsmall," %lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    Tparalogsmall=promedioT;
                                //     printf(" Histograma promedio dibujado\n",iter,L);
                                fclose(histTparalogsmall);
                                LparalogSmallbuscado=L;
                   }
            if(promediolog>=(promlogTbig-0.10) && promediolog<=(promlogTbig+0.10)){
                            histlogbig=fopen("LogsBig.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(logaritmos,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histlogbig,"%lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    logbig=promediolog;
                                        //printf(" Histograma logaritmo dibujado\n");
                                          fclose(histlogbig);
                        //*************HISTOGRAMA DE P(T) RESPECTO A T para el <-log T> buscado********************
                            histTparalogBig=fopen("TparalogBig.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(transmisiones,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histTparalogBig,"%lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    TparalogBig=promedioT;
                                //     printf(" Histograma promedio dibujado\n",iter,L);
                                fclose(histTparalogBig);
                                LparalogBigbuscado=L;               }
            if(promedioT>=(promTsmall-0.0001) && promedioT<=(promTsmall+0.0001)){
                            histTsmall=fopen("TSmall.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(transmisiones,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histTsmall,"%lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    Tsmall=promedioT;
                                //     printf(" Histograma promedio dibujado\n",iter,L);
                                fclose(histTsmall);
                             transmisiones_TSmall=fopen("Transmisiones_TSmall.txt","w");
                            for (i=0;i<M;i++) fprintf(transmisiones_TSmall," %d \t %lf\n",i,transmisiones[i]);
                            fclose(transmisiones_TSmall);
                        //*************HISTOGRAMA DE P(-log T) RESPECTO A -log T para el <T> buscado********************
                            histlogparaTSmall=fopen("LogparaTSmall.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(logaritmos,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histlogparaTSmall,"%lf \t %lf \n",i*delta+minimo+delta/2,Hist[i]);
                                    logparaTSmall=promediolog;
                                      logs_para_TSmall=fopen("Logs_para_TSmall.txt","w");
                            for (i=0;i<M;i++) fprintf(logs_para_TSmall," %d \t %lf\n",i,-logaritmos[i]);
                            fclose(logs_para_TSmall);
                                //     printf(" Histograma promedio dibujado\n",iter,L);
                                fclose(histlogparaTSmall);
                                LparaTSmallBuscado=L;
                                LogparaTSmall;
                   }
            if(promedioT>=(promTbig-0.01) && promedioT<=(promTbig+0.01)){
                        histTbig=fopen("TBig.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(transmisiones,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histTbig," %lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    Tbig=promedioT;
                                 //   printf(" Histograma promedio dibujado\n",iter,L);
                                  fclose(histTbig);
                            transmisiones_TBig=fopen("Transmisiones_TBig.txt","w");
                            for (i=0;i<M;i++) fprintf(transmisiones_TBig," %d \t %lf\n",i,transmisiones[i]);
                            fclose(transmisiones_TBig);
                         //*************HISTOGRAMA DE P(-log T) RESPECTO A -log T para el <T> buscado********************
                        histlogparaTBig=fopen("logparaTBig.txt","w");
                                    for(i=0;i<N_Intervalos;i++) Hist[i]=0.0;
                                    Histograma(logaritmos,Hist,M,N_Intervalos,&delta,&minimo,&maximo); //Calculamos el histograma
                                     for(i=0;i<N_Intervalos;i++) //Escribimos en un archivo los datos
                                            fprintf(histlogparaTBig," %lf \t %lf\n",i*delta+minimo+delta/2,Hist[i]);
                                    logparaTBig=promediolog;
                                //     printf(" Histograma promedio dibujado\n",iter,L);
                                fclose(histlogparaTBig);
                                LparaTbigbuscado=L;
                                LogparaTBig=promediolog;
           }
         #endif // SACAR_PROBABILIDADES
        E+=deltaE;//CAMBIAMOS LA LONGITUD DEL SISTEMA
        if ((int)((E-E_0)/deltaE)%(numE/10)==0) printf("Completado al %lf porciento\n",((double)(iter+1)/((double)(numE)))*100);
    }
    printf(" Sale de bucle de E o L\n");
        fclose(datos);
        fclose(promedios);
        fclose(promedios_log);
        printf(" cierra ficheros\n");
      //  printf("\n t=%lf + i*%lf \n r=%lf + i*%lf", creal(S.t), cimag(S.t),creal(S.r1), cimag(S.r1));
        creaplt("datos_SCH.dat");
        creaplt_proms("promedios_SCH.dat",k*L);
        creaplt_proms_log("promedios_log_SCH.dat",k*L);
       creaplt_fluctuaciones("fluctuaciones_SCH.dat",k*L);
    #ifdef SACAR_PROBABILIDADES
        creaplt_probabilidades_log("LogsBig.txt", logbig,LparalogBigbuscado,"LogsSmall.txt", logsmall,LparalogSmallbuscado);
        creaplt_probabilidades_log_paraesosT("LogparaTBig.txt", logparaTBig,LparaTbigbuscado,"LogparaTSmall.txt", logparaTSmall,LparaTSmallBuscado);
        creaplt_probabilidades_T("TBig.txt", Tbig,LogparaTBig,"TSmall.txt", Tsmall,LogparaTSmall);
         creaplt_probabilidades_T_paraesoslog("TparalogBig.txt", TparalogBig,"TparalogSmall.txt", Tparalogsmall);
    #endif // SACAR_PROBABILIDADES
        printf("\n done");

#endif // CAMBIAR_E
    return 0;
}

void ini_ran(int SEMILLA){ //Parisi-Rapuano, NO TOCAR
    int INI,FACTOR,SUM,i;
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;
    for(i=0;i<256;i++) {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
}

float random(void){ //Nº aleatorio entre [0,1)
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;
    return r;
}


void readata(int *n, int *m, double *E, double *V, double *L, double *theta_i, char*n_arch){
    FILE *arch;
    arch=fopen(n_arch,"r");
    fscanf(arch, "%i\t %i\t %lf\t %lf\t %lf\t %lf",n,m,E,V,L,theta_i);
    //printf("%i\t%i\t%lf\t%lf\t%lf\n",*n,*m,*E,*V,*L);
}


void creaplt(char *n_arch){
    FILE *graf;
    double T,L;
    graf=fopen("grafico_SCH.plt","w");
    fprintf(graf, "set title 'Transmitividad respecto a k_x *L' \n");
    fprintf(graf, "set xlabel 'Vector de onda en componente x * Longitud del sistema' \n");
    fprintf(graf, "set ylabel 'T' \n");
    fprintf(graf, "plot '%s' u 1:2\n", n_arch);
fclose(graf);
}
#ifdef CAMBIAR_L
void creaplt_proms(char *n_arch,double E){
    FILE *graf;
    double T,L;
    graf=fopen("grafico_proms_SCH.plt","w");
    fprintf(graf, "set title 'Promedio de transmitividad respecto a k_x *L para E=%lf' \n",E);
    fprintf(graf, "set xlabel 'Vector de onda en componente x * Longitud del sistema' \n");
    fprintf(graf, "set ylabel '<T>' \n");
    fprintf(graf, "unset key \n");
    fprintf(graf, "plot '%s' u 1:2 w lines\n", n_arch);
fclose(graf);
}
void creaplt_proms_log(char *n_arch,double E){
    FILE *graf;
    double T,L;
    graf=fopen("grafico_proms_log_SCH.plt","w");
    fprintf(graf, "set title 'Promedios de logaritmo de la transmitividad respecto a k_x*L para E=%lf' \n",E);
    fprintf(graf, "set xlabel 'Vector de onda en componente x * Longitud del sistema' \n");
    fprintf(graf, "set ylabel '<ln T>' \n");
    fprintf(graf, "plot '%s' u 1:2 w lines\n", n_arch);
fclose(graf);
}

void creaplt_probabilidades_log(char *n_arch1,double logBig,double LBig,char *n_arch2,double logSmall, double Lsmall){
    FILE *graf1;
    graf1=fopen("distribucion_logs_Big_SCH.plt","w");
    fprintf(graf1, "set title 'Distribución de probabilidad del logaritmo de la transmitividad para L=%lf y para <ln T>=%lf' \n",LBig,logBig);
    fprintf(graf1, "set xlabel 'ln T' \n");
    fprintf(graf1, "set ylabel 'P(ln T)' \n");
    fprintf(graf1, "plot '%s' u 1:2 w boxes\n", n_arch1);
//Cambio a la distribución de logs small
 graf1=fopen("distribucion_logs_small_SCH.plt","w");
    fprintf(graf1, "set title 'Distribución de probabilidad del logaritmo de la transmitividad para L=%lf para <ln T>=%lf' \n",Lsmall,logSmall);
    fprintf(graf1, "set xlabel 'ln T' \n");
    fprintf(graf1, "set ylabel 'P(ln T)' \n");
    fprintf(graf1, "plot '%s' u 1:2 w boxes\n", n_arch2);
fclose(graf1);
}

void creaplt_probabilidades_log_paraesosT(char *n_arch1,double logBig,double LBig,char *n_arch2,double logSmall, double Lsmall ){
    FILE *graf1;
    graf1=fopen("distribucion_logs_para_TBig_SCH.plt","w");
    fprintf(graf1, "set title 'Distribución de probabilidad del logaritmo de la transmitividad para L=%lf para <-ln T>=%lf' \n",LBig,logBig);
    fprintf(graf1, "set xlabel 'ln T' \n");
    fprintf(graf1, "set ylabel 'P(ln T)' \n");
    fprintf(graf1, "plot '%s' u 1:2 w boxes\n", n_arch1);
 graf1=fopen("distribucion_logs_para_Tsmall_SCH.plt","w");
    fprintf(graf1, "set title 'Distribución de probabilidad del logaritmo de la transmitividad respecto a L=%lf para <-ln T>=%lf' \n",Lsmall,logSmall);
    fprintf(graf1, "set xlabel 'ln T' \n");
    fprintf(graf1, "set ylabel 'P(ln T)' \n");
    fprintf(graf1, "plot '%s' u 1:2 w boxes\n", n_arch2);
fclose(graf1);
}
void creaplt_probabilidades_T(char *n_arch1,double TBig,double logparaTBig, char *n_arch2,double TSmall,double logparaTSmall){
    FILE *graf1;
    double s=logparaTBig;
    graf1=fopen("distribucion_T_Big_SCH.plt","w");
    fprintf(graf1, "set title 'Distribución de probabilidad de la transmitividad respecto a L para <T>=%lf' \n",TBig);
    fprintf(graf1, "set xlabel 'T' \n");
    fprintf(graf1, "set ylabel 'P(T)' \n");
    fprintf(graf1, "C=0.3\n");//PARA EL FIT
    fprintf(graf1, "s=%lf\n",s);
 fprintf(graf1, "f(x)=C*(sqrt(acosh(1/sqrt(x)))*exp(-((acosh(1/sqrt(x)))*(acosh(1/sqrt(x))))/s))/(x**(3/2)*sqrt(sqrt((1-x)))) \n");
    fprintf(graf1, "fit f(x) '%s' u 1:2 via C,s\n", n_arch1);
    fprintf(graf1, "plot '%s' u 1:2 w boxes\n", n_arch1);
     fprintf(graf1, "replot f(x)");
s=logparaTSmall;
 graf1=fopen("distribucion_T_small_SCH.plt","w");
    fprintf(graf1, "set title 'Distribución de probabilidad de la transmitividad respecto a L para <T>=%lf' \n",TSmall);
    fprintf(graf1, "set xlabel 'T' \n");
    fprintf(graf1, "set ylabel 'P( T)' \n");
      fprintf(graf1, "C=0.3\n");//PARA EL FIT
    fprintf(graf1, "s=%lf\n",s);
    fprintf(graf1, "f(x)=C*(sqrt(acosh(1/sqrt(x)))*exp(-((acosh(1/sqrt(x)))*(acosh(1/sqrt(x))))/s))/(x**(3/2)*sqrt(sqrt((1-x)))) \n");
    fprintf(graf1, "fit f(x) '%s' u 1:2 via C,s\n", n_arch2);
    fprintf(graf1, "plot '%s' u 1:2 w boxes\n", n_arch2);
      fprintf(graf1, "replot f(x)");
fclose(graf1);
}

void creaplt_probabilidades_T_paraesoslog(char *n_arch1,double TBig,char *n_arch2,double TSmall){
    FILE *graf1;
    double s=-log(TBig);
    graf1=fopen("distribucion_T_paralogBig_SCH.plt","w");
    fprintf(graf1, "set title 'Distribución de probabilidad de la transmitividad respecto a L para <T>=%lf' \n",TBig);
    fprintf(graf1, "set xlabel 'T' \n");
    fprintf(graf1, "set ylabel 'P(T)' \n");
    fprintf(graf1, "C=0.3\n");//PARA EL FIT
    fprintf(graf1, "s=%lf\n",s);
 fprintf(graf1, "f(x)=C*(sqrt(acosh(1/sqrt(x)))*exp(-((acosh(1/sqrt(x)))*(acosh(1/sqrt(x))))/s))/(x**(3/2)*sqrt(sqrt((1-x)))) \n");
    fprintf(graf1, "fit f(x) '%s' u 1:2 via C,s\n", n_arch1);
    fprintf(graf1, "plot '%s' u 1:2 w boxes\n", n_arch1);
     fprintf(graf1, "replot f(x)");
s=-log(TSmall);
 graf1=fopen("distribucion_T_paralogSmall_SCH.plt","w");
    fprintf(graf1, "set title 'Distribución de probabilidad de la transmitividad respecto a L para <T>=%lf' \n",TSmall);
    fprintf(graf1, "set xlabel 'T' \n");
    fprintf(graf1, "set ylabel 'P( T)' \n");
      fprintf(graf1, "C=0.3\n");//PARA EL FIT
    fprintf(graf1, "s=%lf\n",s);
    fprintf(graf1, "f(x)=C*(sqrt(acosh(1/sqrt(x)))*exp(-((acosh(1/sqrt(x)))*(acosh(1/sqrt(x))))/s))/(x**(3/2)*sqrt(sqrt((1-x)))) \n");
    fprintf(graf1, "fit f(x) '%s' u 1:2 via C,s\n", n_arch2);
    fprintf(graf1, "plot '%s' u 1:2 w boxes\n", n_arch2);
      fprintf(graf1, "replot f(x)");
fclose(graf1);
}
#endif // CAMBIAR_L

#ifdef CAMBIAR_E
void creaplt_proms(char *n_arch,double L){
    FILE *graf;
    double T;
    graf=fopen("grafico_proms_SCH.plt","w");
    fprintf(graf, "set title 'Promedio de transmitividad respecto a E para L=%lf' \n",L);
    fprintf(graf, "set xlabel 'Energía incidente' \n");
    fprintf(graf, "set ylabel '<T>' \n");
    fprintf(graf, "unset key \n");
    fprintf(graf, "plot '%s' u 1:2 w lines\n", n_arch);
fclose(graf);
}
void creaplt_proms_log(char *n_arch,double L){
    FILE *graf;
    double T;
    graf=fopen("grafico_proms_log_SCH.plt","w");
    fprintf(graf, "set title 'Promedios de logaritmo de la transmitividad respecto a E para L=%lf' \n",L);
    fprintf(graf, "set xlabel 'Energía incidente' \n");
    fprintf(graf, "set ylabel '<ln T>' \n");
    fprintf(graf, "plot '%s' u 1:2 w lines\n", n_arch);
fclose(graf);
}
void creaplt_fluctuaciones(char *n_arch,double L){
    FILE *graf;
    double T;
    graf=fopen("grafico_fluctuaciones_SCH.plt","w");
    fprintf(graf, "set title 'Fluctuación de la transmitividad respecto a E para L=%lf' \n",L);
    fprintf(graf, "set xlabel 'Energía incidente' \n");
    fprintf(graf, "set ylabel 'Fluctuación' \n");
    fprintf(graf, "plot '%s' u 1:2 w lines\n", n_arch);
fclose(graf);
}
#endif // CAMBIAR_E

S_matrix compose(S_matrix S1, S_matrix S2){
    S_matrix composed;
    composed.r1=S1.r1+S1.t*S2.r1*S1.t/(1-S1.r2*S2.r1);
    composed.r2=S2.r2+S2.t*S1.r2*S2.t/(1-S1.r2*S2.r1);
    composed.t=S2.t*S1.t/(1-S1.r2*S2.r1);
    return composed;
}

S_matrix translate(S_matrix S, double d, double k){
    S_matrix translated;
    translated.r1=S.r1*cexp(2*k*d*I);
    translated.r2=S.r2*cexp(-2*k*d*I);
    translated.t=S.t;
    return translated;
}

S_matrix barrier(double a, double V, double k){
    /*
    kanorma corresponde al valor de k*a normalizado para entrar en el intervalo [0,2PI]
     y kbaranorma es lo mismo para el vector de onda en el interior de una barrera de potencial V *a
    */
    double K,kbar;
    double complex den,kanorma,kbaranorma;
    int kaentrepientero,kbaraentrepientero;
    S_matrix S;
    K=sqrt(2*me*V/hbar/hbar);
    kbar=sqrt(k*k-K*K);
    kaentrepientero=((int)k*a)/((int)(2*PI));
    kanorma=k*a-((double)kaentrepientero*2*PI);
    kbaraentrepientero=((int)kbar*a)/((int)(2*PI));
    kbaranorma=k*a-((double)kbaraentrepientero*2*PI);
   // printf("kbar=%lf\n", kbar);
 //  printf("cos(kbar*a)=%lf, \t k*k+kbar*kbar=%lf \t (2*k*kbar)*sin(kbar*a)=%lf y la parte imaginaria es %lf\n", cos(kbar*a), k*k+kbar*kbar,(2*k*kbar)*sin(kbar*a), (k*k+kbar*kbar)/(2*k*kbar)*sin(kbar*a));
    den=(cos(kbaranorma)-I*(k*k+kbar*kbar)/(2*k*kbar)*sin(kbaranorma));
   // printf("den=%lf \t cos(kbar*a)=%lf\n", den,cos(kbar*a));
    S.r1=-I*K*K/(2*k*kbar)*sin(kbaranorma)/den;
    S.r2=S.r1*cexp(-2*I*kanorma);
    S.t=cexp(-I*kanorma)/den;
    return S;
}

int check(S_matrix S){
    if(cabs(S.r1*conj(S.r1)+S.t*conj(S.t)-1)<eps||cabs(S.r1*conj(S.t)+S.t*conj(S.r2))<eps||cabs(conj(S.r1)*(S.t)+S.r2*conj(S.t))<eps||cabs(S.r2*conj(S.r2)+S.t*conj(S.t)-1)<eps)
        return 0;
    else return 1;
}


void Histograma(double *data, double *Hist,int NDatos,int n_Intervalos,double *d,double *m, double *M){
/* Genera un histograma. Calcula, en funcion de los datos el minimo y maximo
de los mismos para ajustar mejor los intervalos.
*data-> input Datos sobre los que se genera el histograma
*Hist-> output Histograma calculado
N_datos-> input Numero de datos
N_Intervalos-> input Numero de intervalos del histograma
*d -> output Medida de cada intervalo del histograma
*m -> output Valor minimo de los datos
*M -> output Valor maximo de los datos
*/
int i,Indice;
double del,min,max,Norm;
for(i=0;i<N_Intervalos;i++)//Inicializo
    Hist[i]=0;
//min=10000000;
//max=-10000000;
min=max=0;
for(i=0;i<NDatos;i++) //Calculo el minimo y el maximo
{
if(data[i]<min)min=data[i];
if(data[i]>max)max=data[i];
}
del=(max-min)/N_Intervalos;
if(del==0)
{
printf("Error: No se pueden calcular los intervalos; Max=%lf,Min=%lf\n",max,min);
exit(1);
}
for(i=0;i<NDatos;i++) //*********** Calculo el histograma: Nucleo del Programa**********
{
Indice=(data[i]-min)/del;
Hist[Indice]++;

}
// ************** Fin del nucleo del programa ***************
*d=del;
*m=min;
*M=max;
/* Ahora normalizo */
// Recordar: 1=A*sum(h_i * delta_i)=A*delta*sum(h_i)=A*delta*N => A=1/(delta*N)
Norm=1.0/(NDatos*del);
for(i=0;i<N_Intervalos;i++)
Hist[i]*=Norm;
}

double random_Levy_general(double alpha_Levy){
    double x,Z,W,beta,theta,theta_0,gamma,delta,randomW,thetanorma;
    int thetaentrepientero;
  //  ini_ran(time (NULL));
    // printf("atrapado\n");
     theta=random()*PI-PI/2.0;
     thetaentrepientero=((int)theta)/((int)(2*PI));
    thetanorma=theta-((double)thetaentrepientero*2*PI);
    gamma=1.0;//DIST. estandarizada
    delta=0.0;
    randomW=random();
    W=-log(randomW);
    //printf("theta=%lf \t, W=%lf\t, randomW=%lf, beta=%lf\n",theta,W,randomW, beta);
  //  beta=-log(1-random())+1;//En general
   if (alpha_Levy<1) {
        beta=1.0;//Queremos que esté cargado a un lado
        theta_0=PI/2.0;
   }
    else {
        beta=0.0;
        theta_0=0.0;
     }
  //theta_0=(atan(beta*tan(PI*(alpha_Levy)/2.0)))/(alpha_Levy);//SI BETA ES DISTINTA DE 1 o 0
//theta_0=PI/2.0;//Si beta=1
    if(alpha_Levy==1)//Z=2.0/PI*((PI/2.0+beta*theta)*tan(theta)-beta*log(PI/2.0*W*cos(theta)/(PI/2.0+beta*theta)));
            Z=tan(thetanorma)-beta*log(PI/2.0*W*cos(thetanorma)/(PI/2.0+beta*theta));
        else Z=sin((alpha_Levy)*(theta_0+theta))*pow(((cos((alpha_Levy)*theta_0+((alpha_Levy)-1.0)*theta))/W),((1.0-(alpha_Levy))/(alpha_Levy)))/pow((cos((alpha_Levy)*theta_0)*cos(theta)),(1.0/(alpha_Levy)));
    if (Z<0.0) Z=-Z;
return Z;
}
