#include <math.h>
#include <complex.h>
#include <stdio.h>
/* Illustrative constants and a CKM calculation. Not a complete physics solver.
   Equations in comments are not compiled into an executable theory.
   Values are historical illustrative inputs, not a current global fit. */

/* 27 illustrative values, including neutrino differences and gravity. */
double P[]={
 .3575,.6514,1.221,              /* g1,g2,g3 */
 125.25,246.22,                   /* mH,vH */
 .00216,1.27,172.69,              /* mu,mc,mt */
 .00467,.0934,4.18,               /* md,ms,mb */
 .000511,.10566,1.7768,           /* me,mmu,mtau */
 7.53e-5,2.453e-3,               /* dm21sq,dm32sq */
 .2257,.814,.117,.353,            /* CKM: lam,A,rho,eta */
 .5903,.8587,.1496,-1.601,        /* PMNS: t12,t23,t13,dCP */
 0.,                              /* theta_QCD */
 6.674e-11,1.1e-52               /* G_N, Lambda */
};

/* Gauge: SU(3)xSU(2)xU(1). Fermion reps per generation: */
int R[][3]={{3,2,1},{3,1,4},{3,1,-2},{1,2,-3},{1,1,-6}};
/* hypercharge Y = R[i][2]/6. 3 generations. Higgs: (1,2,3) -> Y=1/2 */

/* Illustrative unbarred Wolfenstein inputs define the three angles and phase.
   Standard complex CKM parametrization; not an electroweak simulation. */
void ckm(double complex V[9]) {
 double lam=P[16], A=P[17], rho=P[18], eta=P[19];
 double s12=lam, s23=A*lam*lam;
 double s13=s23*lam*hypot(rho,eta), delta=atan2(eta,rho);
 double c12=sqrt(1-s12*s12), c23=sqrt(1-s23*s23), c13=sqrt(1-s13*s13);
 double complex ep=cexp(I*delta), em=conj(ep);
 V[0]=c12*c13; V[1]=s12*c13; V[2]=s13*em;
 V[3]=-s12*c23-c12*s23*s13*ep;
 V[4]=c12*c23-s12*s23*s13*ep; V[5]=s23*c13;
 V[6]=s12*s23-c12*c23*s13*ep;
 V[7]=-c12*s23-s12*c23*s13*ep; V[8]=c23*c13;
}

/* Yukawa: y_i = sqrt(2)*m_i/vH. Higgs potential: V=-mH^2/2*|p|^2+mH^2/(2*vH^2)*|p|^4 */
/* Lagrangian: -F^2/4 + i*psi_bar*gamma*D*psi + |D*phi|^2 - V(phi) - y*psi_L*phi*psi_R */
/* GR: S=(R-2*Lambda)/(16*pi*G)*sqrt(-g)*d4x. G_uv+Lambda*g_uv=8*pi*G*T_uv */
/* QM: i*hbar*d|psi>/dt=H|psi>, P(a)=|<a|psi>|^2 */

int main(void) {
 double complex V[9]; ckm(V);
 printf("parameters %zu\n",sizeof(P)/sizeof(P[0]));
 for(int i=0;i<9;i++) printf("%.17g %.17g\n",creal(V[i]),cimag(V[i]));
 return 0;
}
