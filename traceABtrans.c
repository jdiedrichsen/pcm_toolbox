/*==========================================================
%  Mex version to compute the trace of the multiplication of A and B' 
%  
%  The function uses the trace trick to speed up estimation times:
%  trace(A*B')=sum(sum(A.*B))
%
% INPUT:
%   A: NxK matrix
%   B: NxK matrix
%
% OUTPUT:
%   tr : the trace of A times B'
%
% WARING: for speed reasons the mex version of this function skips all
% checks involving the size of the matrices, these checks should be
% implemented in the calling function
% (c) Naveed Ejaz 2013 n.ejaz@ucl.ac.uk
 *========================================================*/

#include "mex.h" 

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    // (1) check for proper number of input and output arguments 
    if(nrhs!=2)
        mexErrMsgIdAndTxt("MyToolbox:traceAB_faster:nrhs","2 inputs required.");
    
    int N,K;
    int i,j, idx1,idx2;
    double *ptrA, *ptrB;
    double tr=0;
    
    // getting dimensions
    N = mxGetM(prhs[0]);
    K = mxGetN(prhs[0]);   
    
    // array pointers
    ptrA = mxGetPr(prhs[0]);      
    ptrB = mxGetPr(prhs[1]);      
    
    // simple matrix multiplication
    for (j=0; j<K; j++) { 
        idx1=j*N;
        for(i=0; i<N; i++) {
            idx2=i+idx1;
            tr += ptrA[idx2]*ptrB[idx2];
        } 
    }

    // (5) assigning outputs and returning to matlab
    plhs[0] = mxCreateDoubleScalar(tr);
    return;
}

