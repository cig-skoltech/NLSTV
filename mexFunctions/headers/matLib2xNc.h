#include <string.h> /* needed for memcpy() */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"
#include "epph.h" /* This is the head file that contains the implementation of the used functions for the projections*/
#include <blas.h>

// #if !defined(_WIN32)
// #define dgemm dgemm_
// #endif

#define MAXITER 100 //Maximum allowed number of iterations.
#define EPS_M 1e-15 //Machine double-point precision.

#ifndef max
#define max(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef min
#define min(a, b) ((a)<(b)?(a):(b))
#endif

#ifndef COND_ABS
#define COND_ABS(a, b) (b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a))
#endif

#ifndef sign
#define sign(a) (a >= 0 ? (a > 0 ? 1 : 0) : -1)
#endif



void eigen2x2Sym(double A[4], double U[2], double E[2]) {
  
  double k;
  double trace_A;
  double sqrt_delta_A;
  
  if (fabs(A[1]) < 1e-15){
    E[0]=A[0];
    E[1]=A[3];
    U[0]=1.0;
    U[1]=0.0;}
  else
  {
    trace_A=A[0]+A[3];
    sqrt_delta_A=sqrt((A[0]-A[3])*(A[0]-A[3])+4*A[1]*A[1]);
    E[0]=0.5*(trace_A+sqrt_delta_A);
    E[1]=0.5*(trace_A-sqrt_delta_A);
    k=sqrt((E[0]-A[0])*(E[0]-A[0])+A[1]*A[1]);
    U[0]=A[1]/k;
    U[1]=(E[0]-A[0])/k;
  }
}

void matmulAAT(double *A, double *C, mwSignedIndex m, mwSignedIndex n) {
  // A is a matrix of size m x n and C=A*A^T is of size m x m.
  char *chA = "N";
  char *chAT = "T";
  /* scalar values to use in dgemm */
  double one = 1.0, zero = 0.0;
  
  dgemm(chA, chAT, &m, &m, &n, &one, A, &m, A, &m, &zero, C, &m);  
  
}


void matmul(double *A, double *B, double *C, mwSignedIndex m, mwSignedIndex n, mwSignedIndex p) {
  /* A is a matrix of size m x p and B is a matrix of size p x n 
    C=A*B is a matrix of size m x n. */
  
  char *chn = "N";
  /* scalar values to use in dgemm */
  double one = 1.0, zero = 0.0;
  
  dgemm(chn, chn, &m, &n, &p, &one, A, &m, B, &p, &zero, C, &m);  
  
}

void matmulD(double *A, double *D, double *C, mwSignedIndex m, mwSignedIndex p) {
  // A and C are matrices of size m x p and D is a vector of size p x 1 
  // C=A*diag(D)
  
  int j;
  double *B=(double *)mxCalloc(p*p, sizeof(double));
  //Make a diagonal matrix of size p x p out of D
  
  for (j=0;j<p;j++)
    B[j*(p+1)]=D[j];
  
  /* form of op(A) & op(B) to use in matrix multiplication */
  char *chn = "N";
  /* scalar values to use in dgemm */
  double one = 1.0, zero = 0.0;
  
  dgemm(chn, chn, &m, &p, &p, &one, A, &m, B, &p, &zero, C, &m);  
  
  mxFree(B);
}

double pnorm(double *x, double k, double p) {
//x:vector of size k
//p: the order of the norm
  
  double res=0;
  int i;
  
  if (p < 1)
    mexErrMsgTxt("The order of the norm should be greater or equal to one");
  
  if (mxIsFinite(p)){
    for (i=0; i < k; i++)
      res+=pow(fabs(x[i]), p);
    res=pow(res, 1.0/p);}
  else{
    for (i=0; i < k; i++)
      res=max(res, fabs(x[i]));}
  
  return res;
}


double rootfind(double * x_, double *c_, int *iter_step, double * v_, int k_, double p_, double c0, double tau, double x1, double x2, double tol)
/*Using Brent's method, find the root of the (epp-tau) function known to lie between x1 and x2. The
 * root, returned as rootfind, will be refined until its accuracy is tol.*/
{
  int iter;
  double a=x1, b=x2, c=x2, d, e, min1, min2;
  double fa, fb;
  
  epp(x_, c_, iter_step, v_, k_, a, p_, c0);
  fa=pnorm(x_, k_, p_)-tau;
  
  
  epp(x_, c_, iter_step, v_, k_, b, p_, c0);
  fb=pnorm(x_, k_, p_)-tau;
  
  double fc, p, q, r, s, tol1, xm;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    mexErrMsgTxt("The specified interval doesn't contain a root");
  fc=fb;
  for (iter=1;iter<=MAXITER;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a; // Rename a, b, c and adjust bounding interval
      fc=fa; //d.
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS_M*fabs(b)+0.5*tol; //Convergence check.
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa; //Attempt inverse quadratic interpolation.
      if (a == c) {
        p=2.0*xm*s;
        q=1.0-s;
      } else {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q; //Check whether in bounds.
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d; //Accept interpolation.
        d=p/q;
      } else {
        d=xm; //Interpolation failed, use bisection.
        e=d;
      }
    } else { //Bounds decreasing too slowly, use bisection.
      d=xm;
      e=d;
    }
    a=b; //Move last best guess to a.
    fa=fb;
    if (fabs(d) > tol1) //Evaluate new trial root.
      b += d;
    else
      b += COND_ABS(tol1, xm);
    
    epp(x_, c_, iter_step, v_, k_, b, p_, c0);
    fb=pnorm(x_, k_, p_)-tau;
  }
  mexErrMsgTxt("Maximum number of iterations exceeded in rootfind");
  return 0.0; //Never get here.
}

void svd2_x_nc(double *X, mwSize nc, double *U, double *S){
  int j;
  double J[4];
  matmulAAT(X, J, 2, nc);//J=X*X^T
  eigen2x2Sym(J, U, S);// Compute the eigenvalues and the left eigenvector of J.
  
  for (j=0;j<2;j++){
    if ((S[j]<0) || (fabs(S[j])<=EPS_M))
      S[j]=0;
    else
      S[j]=sqrt(S[j]);}// The singular values of X are equal to sqrt(E) of J.
}

void svr2_x_nc(double *Xp, double *X, mwSize nc, double *U, double *S, double *Sp){
  int j;
  double R[4];
  //R=U*P(S)*S^(-1)*U^T where P(S) is the projection of S.
  double D[2];
  for (j=0; j < 2; j++){
    if (fabs(S[j])<=EPS_M)
      D[j]=0;
    else
      D[j]=Sp[j]/S[j];
  }
  
  R[0]=U[0]*U[0]*D[0]+U[1]*U[1]*D[1];
  R[1]=U[0]*U[1]*(D[0]-D[1]);
  R[2]=U[0]*U[1]*(D[0]-D[1]);
  R[3]=U[1]*U[1]*D[0]+U[0]*U[0]*D[1];
  
  matmul(R, X, Xp, 2, nc, 2);// Xp=R*X;
}



void projectS2(double *Xp, double *X, mwSize nc, double rho){
  double normF;
  unsigned int j;
  double J[4];
  
  matmulAAT(X, J, 2, nc);//J=X*X^T
  
  normF=sqrt(J[0]+J[3]);// Frobenius norm of X ||X||_F=sqrt(trace(X*X^T))
  
  if (normF <= rho){
    //matrix reconstruction
    memcpy(Xp, X, (2*nc)*sizeof(double));}
  else{
    //matrix projection
    for (j=0; j < 2*nc; j++)
      Xp[j]=(X[j]/normF)*rho;
  }
}


void projectS1(double *Xp, double *X, mwSize nc, double rho){
  
  double norm_l1;
  double S[2]; // Singular Values
  double Sp[2];// Projected Singular Values
  double U[2]; // Left Singular vector
  
  svd2_x_nc(X, nc, U, S);
  
  norm_l1=S[0]+S[1];//Nuclear norm.
  
  if (norm_l1 <= rho){
    //matrix reconstruction
    memcpy(Xp, X, (2*nc)*sizeof(double));}
  else{
    //matrix projection
    //Projection of the singular values
    double gamma;
    if (rho < fabs(S[0]-S[1]))
      gamma=max(S[0],S[1])-rho;
    else
      gamma=(norm_l1-rho)/2;
    
    Sp[0]=max(S[0]-gamma, 0.0);
    Sp[1]=max(S[1]-gamma, 0.0);
    
    //matrix reconstruction
    svr2_x_nc(Xp, X, nc, U, S, Sp);
  }
}

void projectSinf(double *Xp, double *X, mwSize nc, double rho){
  
  double norm_linf;
  double S[2]; // Singular Values
  double Sp[2];// Projected Singular Values
  double U[2]; // Left Singular vector
  
  svd2_x_nc(X, nc, U, S);
  
  norm_linf=max(S[0],S[1]);//Spectral norm.
  
  if (norm_linf <= rho){
    //matrix reconstruction
    memcpy(Xp, X, (2*nc)*sizeof(double));}
  else{
    //matrix projection
    //Projection of the singular values
    Sp[0]=min(S[0], rho);
    Sp[1]=min(S[1], rho);
    
    //matrix reconstruction
    svr2_x_nc(Xp, X, nc, U, S, Sp);
  }
}

void projectSp(double *Xp, double *X, mwSize nc, double rho, double p, double c0){
  
  double norm_lp;
  double S[2]; // Singular Values
  double Sp[2];// Projected Singular Values
  double U[2]; // Left Singular vector
  
  svd2_x_nc(X, nc, U, S);
  
  norm_lp=pnorm(S,2,p);//lp norm.
  
  if (norm_lp <= rho){
    //matrix reconstruction
    memcpy(Xp, X, (2*nc)*sizeof(double));}
  else{
    int steps[2];
    double c[1];
    double q=p/(p-1);
    double ub=pnorm(S,2,q);
    //matrix projection
    //Projection of the singular values
    double rho_opt=rootfind(Sp, c, steps, S, 2, p, c0, rho, 0, ub, 1e-8);
    epp(Sp, c, steps, S, 2, rho_opt, p, c0);
    
    //matrix reconstruction
    svr2_x_nc(Xp, X, nc, U, S, Sp);
  }
}


void proxS2(double *Xp, double *X, mwSize nc, double rho){
  unsigned int j;
  double normF=0;
  
  for(j=0;j<2*nc;j++)
    normF=normF+X[j]*X[j];
           
  normF=sqrt(normF);

  if (normF==0){
    //matrix reconstruction
    memcpy(Xp, X, 2*nc*sizeof(double));}
  else{
    double k=max((normF-rho),0.0)/normF;
    for (j=0; j < 2*nc; j++)
      Xp[j]=X[j]*k;}
}


void proxSp(double *Xp, double *X, mwSize nc, double rho, double p, double c0){
  
  double U[2]; //Left singular vectors of X.
  double S[2]; //Singular values of X.
  double Sp[2];//Projected singular values.

  svd2_x_nc(X, nc, U, S);
     
  int steps[2];
  double c[1];
  epp(Sp, c, steps, S, 2, rho, p, c0);
  
  //matrix reconstruction
  svr2_x_nc(Xp, X, nc, U, S, Sp);
}


void projectlp(double *xp, double *x, double n, double rho, double p, double * c, double c0 ){
  // n=numel(x)
  double ub=pnorm(x, n, p);
  if (ub <= rho){
    unsigned int j;
    //matrix reconstruction
    for (j=0; j < n; j++)
      xp[j]=x[j];}
  else{
    int steps[2];
    double q=p/(p-1);
    ub=pnorm(x, n, q);
    double rho_opt=rootfind(xp, c, steps, x, n, p, c0, rho, 0, ub, 1e-8);
    epp(xp, c, steps, x, n, rho_opt, p, c0);
  }
}

void projectlinf(double *xp, double *x, double n, double rho, double p){
  // n=numel(x)
  unsigned int j;
  double ub=pnorm(x, n, p);
  if (ub <= rho){
    //vector stays as it is
    for (j=0; j < n; j++)
      xp[j]=x[j];}
  else{
    for (j=0; j < n; j++)
      xp[j]=sign(x[j])*min(fabs(x[j]), rho);}
}





