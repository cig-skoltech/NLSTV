#include <stdlib.h>
#include <mex.h>
#include <math.h>
#include <omp.h>
#include "matrix.h"
#include "matLib2xNc.h"

/*mex -v mexFunctions/source/proxSpMat2xNc.c -lmwblas
 * CFLAGS="\$CFLAGS -fopenmp -std=c99" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims
 * -ImexFunctions/headers/ -outdir  mexFunctions/source/ */

/* In a mxArray to access the element X[i][j][z] you can do it by referring
 * to the element X[i+j*dims[0]+z*dims[0]*dims[1]]
 */

/* =========================================================================
 * %
 * %  Author: s.lefkimmiatis@skoltech.ru
 * %
 * % =========================================================================*/


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  
  /*set up input arguments */
  if (nrhs < 3 || nrhs > 4)
    mexErrMsgTxt("At least 3 and at most 4 input arguments are expected.\n");
  
  int i;
  
  const mxClassID cid = mxGetClassID(prhs[0]);
  
  if (cid!=mxDOUBLE_CLASS)
    mexErrMsgTxt("Unsupported data type. Accepted type of data is double.");
  
  double* X= mxGetPr(prhs[0]); // matrix input of size nx x ny x ... x 2 x nc
  int  number_of_dims=mxGetNumberOfDimensions(prhs[0]);
  const mwSize *dims=mxGetDimensions(prhs[0]);
  
  mwSize nc;
  if (number_of_dims < 4)
    nc=1;
  else
    nc=dims[number_of_dims-1];
  
  size_t numel_X=mxGetNumberOfElements(prhs[0]);
  int num_of_mat=numel_X/(2*nc);//number of input matrices.
  
  if ((number_of_dims >= 4) && (dims[number_of_dims-2]!=2))
    mexErrMsgTxt("The (n-1)th dimension of the input matrix should be equal to 2.\n");
  
  if ((number_of_dims < 4) && (dims[number_of_dims-1]!=2))
    mexErrMsgTxt("The last dimension of the input matrix should be equal to 2.\n");
    
  double p = mxGetScalar(prhs[1]);
  if (p < 1)
    mexErrMsgTxt("The order of the norm should be greater or equal to 1.\n");
    
  double* rho =(double *)mxGetPr(prhs[2]);
  size_t num_of_rho=mxGetNumberOfElements(prhs[2]);
  if (!(num_of_rho==num_of_mat || num_of_rho==1))
    mexErrMsgTxt("rho should be either a scalar or equal to the number of the input matrices.\n");
  
  double c0;
  if (nrhs < 4)
    c0=0.0;
  else
    c0=mxGetScalar(prhs[3]);
  
  plhs[0]= mxCreateNumericArray(number_of_dims, dims, mxDOUBLE_CLASS, mxREAL); //Regularized matrix solution.
  if (plhs[0] == NULL)
    mexErrMsgTxt("Could not create mxArray.\n");
  
  double *Xp=mxGetPr(plhs[0]);
  
  int k;
  double tmp[2*nc];
  double tmp2[2*nc];
  
  if(p==2){
    #pragma omp parallel for shared(X, Xp) private(i,k,tmp,tmp2)
    for(i=0; i < num_of_mat; i++){
      for (k=0;k<2*nc;k++)
        tmp[k]=X[i+num_of_mat*k];
      
      if (num_of_rho==1)
        proxS2(tmp2, tmp, nc, rho[0]);
      else
        proxS2(tmp2, tmp, nc, rho[i]);
      
      for (k=0;k<2*nc;k++)
        Xp[i+num_of_mat*k]=tmp2[k];
    }
  }
  else{
 #pragma omp parallel for shared(X, Xp) private(i,k,tmp,tmp2)
    for(i=0; i < num_of_mat; i++){
      for (k=0;k<2*nc;k++)
        tmp[k]=X[i+num_of_mat*k];
      
      if (num_of_rho==1)
        proxSp(tmp2, tmp, nc, rho[0], p, c0);
      else
        proxSp(tmp2, tmp, nc, rho[i], p, c0);
      
      for (k=0;k<2*nc;k++)
        Xp[i+num_of_mat*k]=tmp2[k];
    }
  }
}

