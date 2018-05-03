#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <math.h>
#include <omp.h>
#include "matrix.h"



#define TRUE 1
#define FALSE 0


/*mex -v Nsum_mex.cpp
 * CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims
 * -outdir  source/ */

/* In a mxArray to access the element X[i][j][z] you can do it by referring
 * to the element X[i+j*dims[0]+z*dims[0]*dims[1]]
 */

double Nsum(double *X, double *V, double *idx, int n, int I, int i)
{
  int k;
  int ctr=0;
  int idx2;
  
  for (k=0; k < n; k++)
  {
    idx2=idx[I+ctr];
    V[i]=V[i]+X[idx2];
    ctr+=1;
  }
  
}

double NsumC(double *Xr, double *Xi, double *Vr, double *Vi, double *idx, int n, int I, int i)
{
  int k;
  int ctr=0;
  int idx2;
  
  for (k=0; k < n; k++)
  {
    idx2=idx[I+ctr];
    Vr[i]=Vr[i]+Xr[idx2];
    Vi[i]=Vi[i]+Xi[idx2];
    ctr+=1;
  }
  
}


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  
  if (nrhs!=4)
    mexErrMsgTxt("Four input arguments are expected.\n");
  /*set up input arguments */
  
  const mxClassID cid = mxGetClassID(prhs[0]);
  
  if (cid!=mxDOUBLE_CLASS)
    mexErrMsgTxt("Unsupported data type. Accepted type of data is double.");
  
  unsigned int j;
  int i;
  double *Xr, *Xi, *Vr, *Vi;
  
  Xr = mxGetPr(prhs[0]); // matrix input of size N x Nw where
  // N = N1 x N2 x N3 ...
  
  if (mxIsComplex(prhs[0])==TRUE){
    Xi = mxGetPi(prhs[0]);
  }
  
  
  double *idx = mxGetPr(prhs[1]); // vector of size (N x Nw) x 1.
  double *n = mxGetPr(prhs[2]); // vector of size N x 1.
  double *I = mxGetPr(prhs[3]); // vector of size N x 1.
  
  int  number_of_dims=mxGetNumberOfDimensions(prhs[0]);
  const mwSize *dims=mxGetDimensions(prhs[0]);
  size_t num_of_X=mxGetNumberOfElements(prhs[0]);
  
  
  int Nw=dims[number_of_dims-1];
  int N=num_of_X/Nw; // N= N1 x N2 ...
  
  mwSize dims_V[number_of_dims-1]; // dimensions of the output matrix V
  for (j=0; j < number_of_dims-1; j++)
    dims_V[j]=dims[j];
  
  if (mxIsComplex(prhs[0])==TRUE)
    plhs[0]=mxCreateNumericArray(number_of_dims-1, dims_V, mxDOUBLE_CLASS, mxCOMPLEX);
  else
    plhs[0]=mxCreateNumericArray(number_of_dims-1, dims_V, mxDOUBLE_CLASS, mxREAL);
  
  if (plhs[0] == NULL)
    mexErrMsgTxt("Could not create mxArray.\n");
  
  if (mxIsComplex(prhs[0])==TRUE){
    Vr=(double *)mxGetPr(plhs[0]);
    Vi=(double *)mxGetPi(plhs[0]);
  }
  else{
    Vr=(double *)mxGetPr(plhs[0]);
  }
  
 
  //int id1,id2;
  
  if (mxIsComplex(prhs[0])==TRUE){
    #pragma omp parallel for shared(Xr,Xi,Vr,Vi) private(i)
    for (i=0; i < N; i++)
    {
      //id1=n[i];
      //id2=I[i];
      NsumC(Xr,Xi,Vr,Vi,idx,(int) n[i],(int) I[i],i);
    }
  }
  else{
    #pragma omp parallel for shared(Xr,Vr) private(i)
    for (i=0; i < N; i++)
    {
      //id1=n[i];
      //id2=I[i];
      Nsum(Xr,Vr,idx,(int) n[i],(int) I[i],i);
      
    }
  }
  
}
