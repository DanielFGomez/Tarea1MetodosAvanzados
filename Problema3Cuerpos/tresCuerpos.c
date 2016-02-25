#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double *initArray(int n_puntos);
void printArray(double * array1, double * array2, double * array3, double * array4, int n_puntos);
void copy(double *origen, double *destino, int n_puntos);
void *resetArray(double *array,int n_puntos);

double fun(double x);
double dp1(double q1);
double dp3(double q1,double q3);
void solve1(double *q1, double *p1, double dt, double Nsteps);
void solve3(double *q1, double *q3, double *p3, double dt, double Nsteps);
void print3(double *q1, double *q3, double *p3, int nPuntos);
double epsilon;



int main(){
	double *q1,*q3,*p1,*p3;
	double dt=0.006;
	//paper son 466667
	int nSteps=500000;
	int nInitCond=83;
	srand48(time(NULL));
	
	epsilon = 1;
	q1=initArray(nSteps);
	q3=initArray(nSteps);
	p1=initArray(nSteps);
	p3=initArray(nSteps);

	//Tomado del paper
	q1[0]=0.35355339059;
	solve1(q1,p1,dt,nSteps);
	

	int i;
	for (i=0;i<nInitCond;i++){
	  resetArray(q3,nSteps);
	  resetArray(p3,nSteps);
	  q3[0]=drand48()*1;
	  solve3(q1,q3,p3,dt,nSteps);
	  print3(p1,q3,p3,nSteps);
	}
	
	return 0;
}

double dp1(double  q1){
  return -2*q1/pow(4*q1*q1+epsilon*epsilon,3/2);
}

double dp3(double q1,double q3){
  return (q1-q3)/pow((q1-q3)*(q1-q3)+epsilon*epsilon/4,3/2) - (q1+q3)/pow((q1+q3)*(q1+q3)+epsilon*epsilon/4,3/2); 
}  

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
/*
void solve3RK(double *q1, double *q3, double *p3, double dt, double Nsteps){
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double x_in;
  double v_in;
  x_in = *x;
  v_in = *v;
  int i;
  for (i=0;i<Nsteps-1
  k1 = deriv_x(t,x_in,v_in);  
  l1 = deriv_v(t,x_in,v_in);  

  k2 = deriv_x(t + delta_t*0.5, x_in + k1*delta_t*0.5, v_in + l1*delta_t*0.5);
  l2 = deriv_v(t + delta_t*0.5, x_in + k1*delta_t*0.5, v_in + l1*delta_t*0.5) ;

  k3 = deriv_x(t + delta_t*0.5, x_in + k2*delta_t*0.5, v_in + l2*delta_t*0.5);
  l3 = deriv_v(t + delta_t*0.5, x_in + k2*delta_t*0.5, v_in + l2*delta_t*0.5);

  k4 = deriv_x(t + delta_t, x_in + k3*delta_t, v_in + l3*delta_t);
  l4 = deriv_v(t + delta_t, x_in + k3*delta_t, v_in + l3*delta_t);

  x_in += (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*delta_t;
  v_in += (l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0)*delta_t;
}
*/

void copy(double *origen, double *destino, int n_puntos){
	int i;
	for(i=0;i<n_puntos;i++){
	destino[i] = origen[i];
	}
}

void printArray(double * array1, double * array2, double * array3, double * array4, int n_puntos){
	int i;
	for(i=0;i<n_puntos;i++){
	  printf("%f %f %f %f\n", array1[i],array2[i],array3[i],array4[i]);
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

void *resetArray(double *array,int n_puntos){
	int i;
	for(i=0;i<n_puntos;i++){
	array[i] = 0.0;
	}
}

 void print3(double *p1, double *q3, double *p3, int nPuntos){
   int i;
   for (i=0;i<nPuntos-1;i++){
     if(p1[i]*p1[i+1]<=0.0){
       printf("%f %f\n",q3[i],p3[i]); //Asumo que si pasa de positivo a negativo o viceversa es pq paso por el 0.
     }
   }
 }
