/*==========================================================
%  traceAB.c
%  Mex version to compute sum(sum(A.*B')
%
% function tr=traceAB(A,B);;
% INPUT:
%   A: NxK matrix
%   B : KxN matrix
%
% OUTPUT:
%   tr : sum(sum(A.*B')
%
% WARING: for speed reasons the mex version of this function skips all
% checks, they should be implemented in the calling function (mvpattern_covcomp_diag1)
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
    int i,j, idx;
    double *ptrA, *ptrB;
    double tr=0;
    
    // getting dimensions
    N = mxGetM(prhs[0]);
    K = mxGetN(prhs[0]);   
    
    // array pointers
    ptrA = mxGetPr(prhs[0]);      
    ptrB = mxGetPr(prhs[1]);      
    
    // simple matrix multiplication
    for(i=0; i<N; i++)
    {
        idx = i*K;
        for (j=0; j<K; j++)
                tr += ptrA[i+j*N]*ptrB[j+idx];
    }

    // (5) assigning outputs and returning to matlab
    plhs[0] = mxCreateDoubleScalar(tr);
    return;
}

