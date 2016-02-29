#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double *initArray(int n_puntos);

double dp1(double q1);
double dp3(double q1,double q3);
void solve1(double *q1, double *p1, double dt, double Nsteps);
void solve3(double *q1, double *q3, double *p3, double dt, double Nsteps);
void solve1RK(double *q1, double *p1, double dt, double Nsteps);
void solve3RK(double *q1, double *q3, double *p1, double *p3, double dt, double Nsteps);
void print3(double *p1, double *q3, double *p3, int nPuntos, FILE *f);
void printE(double dt, double *q1,double *p1, double Nsteps, FILE *f);
double epsilon;



int main(){
  double *q1,*q3,*p1,*p3,*q1RK,*p1RK,*q3RK,*p3RK,*q1A,*q3A,*p1A,*p3A;
  double dt=0.006;
  int nSteps=466667;
  int nInitCond=83;
  //Inicializa el RNG
  srand48(time(NULL));

  //Inicializar Arrays
  epsilon =0.6;//Este valor dio los mapas que se visualizan mejor
  q1=initArray(nSteps);
  q3=initArray(nSteps);
  p1=initArray(nSteps);
  p3=initArray(nSteps);

  q1RK=initArray(nSteps);
  p1RK=initArray(nSteps);
  q3RK=initArray(nSteps);
  p3RK=initArray(nSteps);

  q1A=initArray(nSteps);
  p1A=initArray(nSteps);
  q3A=initArray(nSteps);
  p3A=initArray(nSteps);
  
  //Resolver el movimiento de las masas grandes
  q1[0]=0.35355339;
  solve1(q1,p1,dt,nSteps);

  q1RK[0]=q1[0];
  solve1RK(q1RK,p1RK,dt,nSteps);

  q1A[0]=0.4325;
  solve1(q1A,p1A,dt,nSteps);
  
  //Crea los archivos donde se guardara la informacion  
  FILE *fileSim = fopen("tresCuerpos.dat", "w");
  FILE *fileRK = fopen("tresCuerposRK.dat", "w");
  FILE *fileA = fopen("tresCuerposA.dat", "w");
  
  FILE *fileE = fopen("energia.dat", "w");
  FILE *fileERK = fopen("energiaRK.dat", "w");

  //Guarda los resultados de la energia para cada metodo de integracion
  printE(dt,q1,p1,nSteps,fileE);
  printE(dt,q1RK,p1RK,nSteps,fileERK);
	 
  int i;
  for (i=0;i<nInitCond;i++){

    //Pone condiciones iniciales aleatorias
    q3[0]=drand48()*6-3;
    p3[0]=drand48()*2-1;
    q3RK[0]=q3[0];
    p3RK[0]=p3[0];
    q3A[0]=drand48()*2-1;
    p3A[0]=drand48()-0.5;
    
    solve3(q1,q3,p3,dt,nSteps);
    solve3RK(q1RK,q3RK,p1RK,p3RK,dt,nSteps);
    solve3(q1A,q3A,p3A,dt,nSteps);
    
    print3(p1,q3,p3,nSteps,fileSim);
    print3(p1RK,q3RK,p3RK,nSteps,fileRK);
    print3(p1A,q3A,p3A,nSteps,fileA);
  }
    
  fclose(fileSim);
  fclose(fileRK);
  fclose(fileA);  
  fclose(fileE);
  fclose(fileERK);
  return 0;
}

//Derivada de p1
double dp1(double  q1){
  return -2*q1/pow(4*q1*q1+epsilon*epsilon,3/2);
}

//Derivada de p3
double dp3(double q1,double q3){
  return (q1-q3)/pow((q1-q3)*(q1-q3)+epsilon*epsilon/4,3/2) - (q1+q3)/pow((q1+q3)*(q1+q3)+epsilon*epsilon/4,3/2); 
}  

//Integrador simplectico
void solve1(double *q1, double *p1, double dt, double Nsteps){
  int i;
  double a0=-pow(2,1/3)/(2-pow(2,1/3));
  double a1=1/(2-pow(2,1/3));
  
  for (i=0;i<Nsteps-1;i++){
    p1[i+1]=p1[i] +0.5*a0*dt*dp1(q1[i]);
    q1[i+1]=q1[i] +    a0*dt*p1[i+1];
    p1[i+1]=p1[i+1] +0.5*a0*dt*dp1(q1[i+1]);
    
    p1[i+1]=p1[i+1] +0.5*a1*dt*dp1(q1[i+1]);
    q1[i+1]=q1[i+1] +    a1*dt*p1[i+1];
    p1[i+1]=p1[i+1] +0.5*a1*dt*dp1(q1[i+1]);
    
    p1[i+1]=p1[i+1] +0.5*a0*dt*dp1(q1[i+1]);
    q1[i+1]=q1[i+1] +    a0*dt*p1[i+1];
    p1[i+1]=p1[i+1] +0.5*a0*dt*dp1(q1[i+1]);
    
      }
}
//Integrador simplectico
void solve3(double *q1, double *q3, double *p3, double dt, double Nsteps){
  int i;
  double a0=-pow(2,1/3)/(2-pow(2,1/3));
  double a1=1/(2-pow(2,1/3));
  
  for (i=0;i<Nsteps-1;i++){
    p3[i+1]=p3[i] +0.5*a0*dt*dp3(q1[i], q3[i]);
    q3[i+1]=q3[i] +    a0*dt*p3[i+1];
    p3[i+1]=p3[i+1] +0.5*a0*dt*dp3(q1[i+1], q3[i+1]);
    
    p3[i+1]=p3[i+1] +0.5*a1*dt*dp3(q1[i+1], q3[i+1]);
    q3[i+1]=q3[i+1] +    a1*dt*p3[i+1];
    p3[i+1]=p3[i+1] +0.5*a1*dt*dp3(q1[i+1], q3[i+1]);
    
    p3[i+1]=p3[i+1] +0.5*a0*dt*dp3(q1[i+1], q3[i+1]);
    q3[i+1]=q3[i+1] +    a0*dt*p3[i+1];
    p3[i+1]=p3[i+1] +0.5*a0*dt*dp3(q1[i+1], q3[i+1]);
  }
  
}
//Integrador Runge Kutta
void solve1RK(double *q1, double *p1, double dt, double Nsteps){
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  
  int i;
  for (i=0;i<Nsteps-1;i++){
    k1 = p1[i];  
    l1 = dp1(q1[i]);  

    k2 = p1[i] + 0.5*l1*dt;
    l2 = dp1(q1[i] + 0.5*p1[i]*dt);

    k3 = p1[i] + 0.5*l2*dt;
    l3 = dp1(q1[i] + 0.5*p1[i]*dt);

    k4 = p1[i] + l3*dt;
    l4 = dp1(q1[i] + p1[i]*dt);

    q1[i+1] = q1[i] +  (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*dt;
    p1[i+1] = p1[i] +  (l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0)*dt;
  }
}
//Integrador Runge Kutta
void solve3RK(double *q1, double *q3, double *p1, double *p3, double dt, double Nsteps){
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  
  int i;
  for (i=0;i<Nsteps-1;i++){
    k1 = p3[i];  
    l1 = dp3(q1[i],q3[i]);  

    k2 = p3[i] + 0.5*l1*dt;
    l2 = dp3(q1[i] + 0.5*p1[i]*dt, q3[i] + 0.5*k1*dt);

    k3 = p3[i] + 0.5*l2*dt;
    l3 = dp3(q1[i] + 0.5*p1[i]*dt, q3[i] + 0.5*k2*dt);

    k4 = p3[i] + l3*dt;
    l4 = dp3(q1[i] + p1[i]*dt, q3[i] + k3*dt);

    q3[i+1] = q3[i] +  (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*dt;
    p3[i+1] = p3[i] +  (l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0)*dt;
  }
}

double *initArray(int n_puntos){
  double *array;
  int i;
  if(!(array = malloc(n_puntos * sizeof(double)))){
    printf("Problema en initArray\n");
    exit(1);
  }
  for(i=0;i<n_puntos;i++){
    array[i] = 0.0;
  }
  return array;
}


void print3(double *p1, double *q3, double *p3, int nPuntos, FILE *f){
   int i;
   for (i=0;i<nPuntos-1;i++){
     if(p1[i]*p1[i+1]<=0.0){
       fprintf(f,"%f %f\n",q3[i],p3[i]);
     }
   }
 }


void printE(double dt, double *q1, double *p1, double Nsteps, FILE *f){
  int i;
  double E;
  for (i=0;i<Nsteps;i++){
    E=0.5*p1[i]*p1[i]-1.0/(2*pow(4*q1[i]*q1[i]+epsilon*epsilon,0.5));
    fprintf(f,"%f %f\n",i*dt,E);
  }
}
