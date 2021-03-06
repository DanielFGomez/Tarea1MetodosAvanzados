#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MIN(x,y) x<y?x:y

/*
C adaptation from C++ code written by Richard J. Gonsalves.
Found here: http://www.physics.buffalo.edu/phy411-506/topic7/index.html
*/

typedef void (*solver)(void);
void solve(solver stepAlgorithm, double tMax, char *filename);

double L = 4.0;                   // length of shock tube
double gama = 1.4;                // ratio of specific heats
int N = 5000;                     // number of grid points

double CFL = 0.4;                 // Courant-Friedrichs-Lewy number

double **U = NULL;                // solution with 3 components
double **newU = NULL;             // new solution
double **F = NULL;                // flux with 3 components

double h;                         // lattice spacing
double tau;                       // time step
double c;                         // speed

int step;

void allocate();
double cMax();
void initialize();
void boundaryConditions(double **U);
void upwindGodunovStep();
void upwindRoeStep();
void LaxFriedrichsStep();
void Riemann(double *U4, double *U1, double *F);
void Roe(double *UL, double *UR, double *F);

int main()
{
   solve(upwindGodunovStep, 1.0, "UpwindGodunov");
//   solve(upwindRoeStep, 1.0, "UpwindRoe");
//   solve(LaxFriedrichsStep, 1.0, "LaxFriedrichs");

   return 0;
}

void solve(solver stepAlgorithm, double tMax, char *filename)
{
  initialize();

  double t = 0.0;
  int j;
  double X, rho, u, e, P;

  tau = CFL*h/cMax();

  while (t < tMax) {
  	boundaryConditions(U);
  	tau = CFL*h/cMax();
	stepAlgorithm();
  	t += tau;
  }

  // write data file with final state
  FILE *out;
  char filename_tmp[1024];
  sprintf(filename_tmp, "%s_finalstate.dat", filename);
  if(!(out = fopen(filename_tmp, "w"))){
	fprintf(stderr, "problem opening file %s\n", filename);
	exit(1);
  }

  for (j = 0; j < N; j++) {

    rho = U[j][0];
    u = U[j][1]/rho;
    e = U[j][2];
    P = (gama - 1.0) * (e - 0.5 * rho * u * u);

    X = j*h;

    fprintf(out, "%f\t%f\t%f\t%f\t%f\n", X, rho, u, e, P);
  }
}

void initialize() {

  int j;
  double rho,p,u,e;
  allocate();

  h = 1.0 * L / (N - 1);

  for (j = 0; j < N; j++) {

    rho = 1;
    p = 1;
    u = 0;

    if (j > N / 2){
      rho = 0.125;
      p = 0.1;
    }

    e = p/(gama-1) + rho*u*u/2.0;

    U[j][0] = rho;
    U[j][1] = rho*u;
    U[j][2] = e;

    F[j][0] = rho*u;
    F[j][1] = rho*u*u+p;
    F[j][2] = u*(e+p);
  }

  tau = CFL*h/cMax();
  step = 0;
}

void allocate() {

  int j;

  U = malloc(N * sizeof(double *));
  newU = malloc(N * sizeof(double *));
  F = malloc(N * sizeof(double *));

  for (j = 0; j < N; j++) {
    U[j] = malloc(3 * sizeof(double));
    newU[j] = malloc(3 * sizeof(double));
    F[j] = malloc(3 * sizeof(double));
  }
}

double cMax() {

  double uMax = 0;
  double rho, u, p, c;
  int i;

  for (i = 0; i < N; i++) {
    if (U[i][0] == 0)
	   continue;

    rho = U[i][0];
    u = U[i][1]/rho;
    p = (U[i][2]-rho*u*u/2)*(gama-1);

    c = sqrt(gama*fabs(p)/rho);

    if (uMax < (c + fabs(u)))
	   uMax = c + fabs(u);
  }

  return uMax;
}


void boundaryConditions(double **U) {
  // reflection boundary conditions at the tube ends
  U[0][0] = U[1][0];
  U[0][1] = -U[1][1];
  U[0][2] = U[1][2];

  U[N-1][0] = U[N-2][0];
  U[N-1][1] = -U[N-2][1];
  U[N-1][2] = U[N-2][2];
}


void upwindGodunovStep() {

  int i, j;

  // find fluxes using Riemann solver
  for (j = 0; j < N - 1; j++){
    Riemann(U[j], U[j + 1], F[j]);
  }

  // update U
  for (j = 1; j < N - 1; j++){
    for (i = 0; i < 3; i++){
      U[j][i] -= tau/h*(F[j][i] - F[j-1][i]);
    }
  }
}


void upwindRoeStep() {

  int i, j;
  double rho_u, u_u, e_u, p_u, a_u;
  double rho_d, u_d, e_d, p_d, a_d;

  // find fluxes using Roe approximation of Riemann solver
  for (j = 0; j < N - 1; j++){
    Roe(U[j], U[j + 1], F[j]);
  }

  // update U
  for (j = 1; j < N - 1; j++){
    for (i = 0; i < 3; i++){
      rho_u = U[j][0];
      u_u = U[j][1]/U[j][0];
      e_u = U[j][2];
      p_u = (gama - 1.0) * (e_u - rho_u * u_u * u_u / 2.0);
      a_u = gama * p_u / rho_u;

      rho_d = U[j-1][0];
      u_d = U[j-1][1]/U[j-1][0];
      e_d = U[j-1][2];
      p_d = (gama - 1.0) * (e_d - rho_d * u_d * u_d / 2.0);
      a_d = gama * p_d / rho_d;
      if(a_u > 0 &&  a_d > 0){
	U[j][i] -= tau/h*(F[j][i] - F[j-1][i]);
      } else if(a_u < 0 && a_d < 0){
	U[j][i] -= tau/h*(F[j+1][i] - F[j][i]);
      } else if(a_u < 0 && a_d > 0){
	U[j][i] -= tau/h*(F[j+1][i] - F[j-1][i]);
      }
    }
  }
}


void LaxFriedrichsStep() {
  int i, j;
  double rho, u, e, p;
    // compute flux F from U
    for (j = 0; j < N; j++) {
      rho = U[j][0];
      u = U[j][1]/U[j][0];
      e = U[j][2];
      p = (gama - 1.0) * (e - rho * u * u / 2.0);

      F[j][0] = rho * u;
      F[j][1] = rho * u * u + p;
      F[j][2] = (e + p) * u;
    }

    // Lax-Friedrichs step
    for (j = 1; j < N - 1; j++)
      for (i = 0; i < 3; i++){
	newU[j][i] = (U[j + 1][i] + U[j - 1][i]) / 2.0 - tau / h * (F[j + 1][i] - F[j - 1][i]);
      }
    boundaryConditions(newU);
    
    // update U from newU
    for (j = 1; j < N - 1; j++)
      for (i = 0; i < 3; i++){
	U[j][i] = newU[j][i];
      }
}


void Riemann(double *U4, double *U1, double *F) {

  const double gamma = 1.4;
  const double g1 = (2*gamma) / (gamma-1);
  const double g2 = (gamma+1) / (2*gamma);
  const double g3 = (gamma+1) / (gamma-1);
  const double tol = 1e-10;

  int i;

  // compute primitive variables
  double rho1, u1, p1, a1;
  rho1 = U1[0];
  u1 = U1[1] / rho1;
  p1 = (U1[2] - 0.5 * rho1 * u1 * u1) * (gamma - 1.0);
  a1 = sqrt(gamma * p1 / rho1);

  double rho4, u4, p4, a4;
  rho4 = U4[0];
  u4 = U4[1] / rho4;
  p4 = (U4[2] - 0.5 * rho4 * u4 * u4) * (gamma - 1.0);
  a4 = sqrt(gamma*p4/rho4);

  // apply the secant method
  double x, y, fx, fy, z, fz;
  x = 0.001 * p4 / p1;
  y = 1.0 * p4 / p1;
  fx = pow(p4/p1, g1) - pow(x, g1)*(1 + ((gamma - 1) / (2 * a4))*(u4 - u1 - (a1 / gamma)*(x - 1)/sqrt(g2 * (x - 1) + 1)));
  fy = pow(p4/p1, g1) - pow(y, g1)*(1 + ((gamma - 1) / (2 * a4))*(u4 - u1 - (a1 / gamma)*(y - 1)/sqrt(g2 * (y - 1) + 1)));

  while(fabs(x-y) > tol){

    z = y - fy * (y - x) / (fy - fx);

    fz = pow(p4/p1, g1) - pow(z, g1)*(1 + ((gamma - 1) / (2 * a4))*(u4 - u1 - (a1 / gamma)*(z - 1)/sqrt(g2 * (z - 1) + 1)));

    x = y;
    fx = fy;
    y = z;
    fy = fz;
  }

  // compute shock
  double rho2, u2, p2, a2;
  p2 = p1 * x;
  u2 = u1 + (a1 / gamma)*((x - 1) / sqrt(g2 * (x - 1) + 1));
  a2 = a1 * sqrt(x * (g3 + x) / (1 + g3 * x));
  rho2 = gamma * p2 / (a2 * a2);

  // compute contact
  double rho3, u3, p3, a3;
  p3 = p2;
  u3 = u2;
  a3 = (u4 + 2.0 * a4 / (gamma - 1) - u3)*(gamma - 1) / 2.0;
  rho3 = gamma*p3/(a3*a3);

  // compute limits
  double s1, s2, s3, s4;
  s1 = u1 + a1 * sqrt(g2 * (x - 1) + 1);
  s2 = u2;
  s3 = u3 - a3;
  s4 = u4 - a4;

  // compute fluxes
  double f1, f2, f3;
  double rho, a, u, p;
  if(s4 > 0) {
    f1 = rho4 * u4;
    f2 = rho4 * u4 * u4 + p4;
    f3 = 0.5 * rho4 * u4 * u4 * u4 + rho4 * a4 * a4 * u4 / (gamma - 1.0);
  } else if (s3 > 0) {
    u = ((gamma - 1.0) * u4 + 2.0 * a4) / (gamma + 1.0);
    a = u;
    p = p4*pow(a / a4, 2.0 * gamma / (gamma - 1.0));
    if (a < 0 || p < 0) {
      printf("Negative a or p in Riemann");
    }
    rho = gamma * p / (a * a);
    f1 = rho * u;
    f2 = rho * u * u + p;
    f3 = 0.5 * rho * u * u * u + rho * a * a * u / (gamma - 1.0);
  } else if (s2 > 0) {
    f1 = rho3 * u3;
    f2 = rho3 * u3 * u3 + p3;
    f3 = 0.5 * rho3 * u3 * u3 * u3 + rho3 * a3 * a3 * u3 / (gamma - 1.0);
  } else if (s1 > 0) {
    f1 = rho2 * u2;
    f2 = rho2 * u2 * u2 + p2;
    f3 = 0.5 * rho2 * u2 * u2 * u2 + rho2 * a2 * a2 * u2 / (gamma - 1.0);
  } else {
    f1 = rho1 * u1;
    f2 = rho1 * u1 * u1 + p1;
    f3 = 0.5 * rho1 * u1 * u1 * u1 + rho1 * a1 * a1 * u1 / (gamma - 1.0);
  }

  F[0] = f1;
  F[1] = f2;
  F[2] = f3;
}

void Roe(double *UR, double *UL, double *F) {
  int i;

  const double gamma = 1.4;
  
  // compute primitive variables for left and right sides
  double rho_r, u_r, p_r, h_r;

  rho_r = UR[0];
  u_r = UR[1] / rho_r;
  p_r = (UR[2] - rho_r * u_r * u_r / 2) * (gamma - 1);
  h_r = (UR[2] + p_r) / rho_r;
  
  double rho_l, u_l, p_l, h_l;
  rho_l = UL[0];
  u_l = UL[1] / rho_l;
  p_l = (UL[2] - rho_l * u_l * u_l / 2) * (gamma - 1);
  h_l = (UL[2] + p_l) / rho_l;

  // compute average properties
  double rho_rl, u_rl, p_rl, h_rl, a_rl;

  rho_rl = sqrt(rho_r * rho_l);
  u_rl = (sqrt(rho_r) * u_r + sqrt(rho_l) * u_l) / (sqrt(rho_r) + sqrt(rho_l));
  h_rl = (sqrt(rho_r) * h_r + sqrt(rho_l) * h_l) / (sqrt(rho_r) + sqrt(rho_l));
  a_rl = sqrt((gamma - 1.0) * (h_rl - 0.5 * u_rl * u_rl));

  // compute the wave speeds
  double l1, l2, l3;

  l1 = u_rl;
  l2 = u_rl + a_rl;
  l3 = u_rl - a_rl;

  // compute the wave strengths
  double d_rho, d_p, d_u;
  double dv1, dv2, dv3;

  d_rho = rho_r - rho_l;
  d_p = p_r - p_l;
  d_u = u_r - u_l;

  dv1 = d_rho - d_p / (a_rl * a_rl);
  dv2 = d_u + d_p / (rho_rl * a_rl);
  dv3 = d_u - d_p / (rho_rl * a_rl);

  // compute the characteristic vectors

  double r1_1, r1_2, r1_3;
  double r2_1, r2_2, r2_3;
  double r3_1, r3_2, r3_3;
  
  r1_1 = 1;
  r1_2 = u_rl;
  r1_3 = u_rl * u_rl * 0.5;

  r2_1 = (rho_rl / (2.0 * a_rl));
  r2_2 = (rho_rl / (2.0 * a_rl)) * (u_rl + a_rl);
  r2_3 = (rho_rl / (2.0 * a_rl)) * (h_rl + a_rl * u_rl);

  r3_1 = -(rho_rl / (2.0 * a_rl));
  r3_2 = -(rho_rl / (2.0 * a_rl)) * (u_rl - a_rl);
  r3_3 = -(rho_rl / (2.0 * a_rl)) * (h_rl - a_rl * u_rl);

  // compute fluxes
  double f1, f2, f3;
  double minl1, minl2, minl3;

  minl1=MIN(0.0, l1);
  minl2=MIN(0.0, l2);
  minl3=MIN(0.0, l3);

  f1 = rho_l * u_l + (r1_1 * minl1 * dv1 + r2_1 * minl2 * dv2 + r3_1 * minl3 * dv3);
  f2 = rho_l * u_l * u_l + p_l + (r1_2 * minl1 * dv1 + r2_2 * minl2 * dv2 + r3_2 * minl3 * dv3);
  f3 = rho_l * h_l * u_l + (r1_3 * minl1 * dv1 + r2_3 * minl2 * dv2 + r3_3 * minl3 * dv3);


//  double rho, a, u, p;
//  if(l3 > 0) {
//    f1 = 0;
//    f2 = 0;
//    f3 = 0;
//  } else if (l1 > 0 && l3 < 0) {
//    f1 = rho_l * u_l + r3_1 * l3 * dv3;
//    f2 = rho_l * u_l * u_l + p_l + r3_2 * l3 * dv3;
//    f3 = rho_l * h_l * u_l + r3_3 * l3 * dv3;
//  } else if (l1 < 0 && l2 > 0) {
//    f1 = rho_l * u_l + r2_1 * l2 * dv2;
//    f2 = rho_l * u_l * u_l + p_l + r2_2 * l2 * dv2;
//    f3 = rho_l * h_l * u_l + r2_3 * l2 * dv2;
//  } else if (l2 < 0) {
//    f1 = 0;
//    f2 = 0;
//    f3 = 0;
//  } else {
//    f1 = 0;
//    f2 = 0;
//    f3 = 0;
//  }

  F[0] = f1;
  F[1] = f2;
  F[2] = f3;
}
