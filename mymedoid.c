#include "mex.h"

void mymedoid(double *outmse, double *ymtr, size_t M, size_t N)
{
    mwSize m;
    mwSize mm;
    mwSize n;
    double armean;
    for(m=0;m<M;m++){
        for(mm=0;mm<M;mm++){
            armean=0.0;
            for(n=0;n<N;n++){
                armean+=(ymtr[m+M*n]-ymtr[mm+M*n])*(ymtr[m+M*n]-ymtr[mm+M*n]);
            }
            outmse[m+M*mm]=sqrt(armean);
        }
    }
    
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{ double *ymtr,*outmse;
 size_t M;
 size_t N;
 //size_t numel;
  M=mxGetM(prhs[0]);
  N=mxGetN(prhs[0]);
  ymtr=mxGetPr(prhs[0]);
  //numel=mxGetN(prhs[0]);
  plhs[0]=mxCreateDoubleMatrix(M, M, mxREAL);
  outmse = mxGetPr(plhs[0]);
  /*call the routine*/
 mymedoid(outmse,ymtr,M,N);
  
  
  
  /* code here */
}