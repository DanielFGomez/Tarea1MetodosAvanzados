#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "roe.h"
#define MIN(x,y) x<y?x:y

/*
  Roe's Riemann Solver    
 */ 


void Roe(double *UR, double *UL, double *F) {
  int i;

  const double gamma = 1.4;
  
    // compute primitive variables for left and right sides
  double rho_r = UR[0];
  double u_r = UR[1] / rho_r;
  double p_r = (UR[2] - rho_r * u_r * u_r / 2) * (gamma - 1);
  double h_r = (UR[2] + p_r) / rho_r;
  double rho_l = UL[0];
  double u_l = UL[1] / rho_l;
  double p_l = (UL[2] - rho_l * u_l * u_l / 2) * (gamma - 1);
  double h_l = (UL[2] + p_l) / rho_l;

  // switch states if necessary so high pressure is on left
  int revflag = FALSE;
  if (p_l < p_r) {
    double swap = p_r; p_r = p_l; p_l = swap;
    swap = u_r; u_r = -u_l; u_l = -swap;
    swap = rho_r; rho_r = rho_l; rho_l = swap;
    revflag = TRUE;
  }

  // compute average properties
  double rho_rl = sqrt(rho_r * rho_l);
  double u_rl = (sqrt(rho_r) * u_r + sqrt(rho_l) * u_l) / (sqrt(rho_r) + sqrt(rho_l));
  double h_rl = (sqrt(rho_r) * h_r + sqrt(rho_l) * h_l) / (sqrt(rho_r) + sqrt(rho_l));
  double a_rl = sqrt((gamma - 1) * (h_rl - 0.5 * u_rl * u_rl));

  // compute the wave velocities
  double l1 = u_rl;
  double l2 = u_rl + a_rl;
  double l3 = u_rl - a_rl;

  // compute the wave strengths
  double d_rho = rho_r - rho_l;
  double d_p = p_r - p_l;
  double d_u = u_r - u_l;

  double dv1 = d_rho - d_p / (a_rl * a_rl);
  double dv2 = d_u + d_p / (rho_rl * a_rl);
  double dv3 = d_u - d_p / (rho_rl * a_rl);

  // compute the characteristic vectors
  double r1_1 = 1;
  double r1_2 = u_rl;
  double r1_3 = u_rl * u_rl * 0.5;

  double r2_1 = (rho_rl / (2 * a_rl));
  double r2_2 = (rho_rl / (2 * a_rl)) * (u_rl + a_rl);
  double r2_3 = (rho_rl / (2 * a_rl)) * (h_rl + a_rl * u_rl);

  double r3_1 = -(rho_rl / (2 * a_rl));
  double r3_2 = -(rho_rl / (2 * a_rl)) * (u_rl - a_rl);
  double r3_3 = -(rho_rl / (2 * a_rl)) * (h_rl - a_rl * u_rl);

  // compute fluxes
  double f1, f2, f3, a, u, rho;
  if (revflag) {
            f1 = -rho_l * u_l - r1_1 * MIN(0,l1) * dv1 - r2_1 * MIN(0,l2) * dv2 - r3_1 * MIN(0,l3) * dv3;
            f2 = rho_l * u_l * u_l + p_l - r1_2 * MIN(0,l1) * dv1 - r2_2 * MIN(0,l2) * dv2 - r3_2 * MIN(0,l3) * dv3;
            f3 = -rho_l * h_l * u_l - r1_3 * MIN(0,l1) * dv1 - r2_3 * MIN(0,l2) * dv2 - r3_3 * MIN(0,l3) * dv3;
    } else {
            f1 = rho_l * u_l + r1_1 * MIN(0,l1) * dv1 + r2_1 * MIN(0,l2) * dv2 + r3_1 * MIN(0,l3) * dv3;
            f2 = rho_l * u_l * u_l + p_l + r1_2 * MIN(0,l1) * dv1 + r2_2 * MIN(0,l2) * dv2 + r3_2 * MIN(0,l3) * dv3;
            f3 = rho_l * h_l * u_l + r1_3 * MIN(0,l1) * dv1 + r2_3 * MIN(0,l2) * dv2 + r3_3 * MIN(0,l3) * dv3;
    }
		


  F[0] = f1;
  F[1] = f2;
  F[2] = f3;
  if(f1>0 || f2>0 || f3>0)
  	printf("%f %f %f\n",f1,f2,f3);
}
