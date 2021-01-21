#include "mex.h"

void findfrontier(double *outmatrix, size_t nobj, double *EC, size_t M)
{ 
    
   mwSize sum,sum2,sum3,sum4,m,n,k;
    for(n=0;n<M;n++){
        outmatrix[n]=0.0;
    }
    for(m=0;m<M;m++){
		for(n=m+1;n<M;n++ ){
			sum=0;sum2=0;sum3=0;sum4=0;
			for(k=0;k<nobj;k++){
				 if(EC[m+k*M]>=EC[n+k*M]){
					 sum=sum+1;}
				 if(EC[m+k*M]>EC[n+k*M]){
					 sum2=sum2+1;}
				 if(EC[n+k*M]>=EC[m+k*M]){
					 sum3=sum3+1;}
				 if(EC[n+k*M]>EC[m+k*M]){
					 sum4=sum4+1;}
			}
			 if(sum==0 && sum2==0){
				 outmatrix[n]=outmatrix[n]+1;}
     
             if(sum3==0 && sum4==0){
				 outmatrix[m]=outmatrix[m]+1;}
        
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{ double *EC,*outmatrix;
  mwSize M;
  mwSize nobj;
  M=mxGetM(prhs[0]);
  EC=mxGetPr(prhs[0]);
  nobj=mxGetN(prhs[0]);
  plhs[0]=mxCreateDoubleMatrix(1, M, mxREAL);
  outmatrix = mxGetPr(plhs[0]);
  /*call the routine*/
  findfrontier(outmatrix,nobj,EC,M);

  

/* code here */
}