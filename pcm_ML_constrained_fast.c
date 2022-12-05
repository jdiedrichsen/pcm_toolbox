/*==========================================================
%  ML_constrained_fast.c
%  Mex version of ML_constrained2  
% constrained maximum likelihood estimation of regression parameters
% function [COV,VA]=ML_constrained_fast(rv,vv,Cc,CcCc);
%   Solves:
%   y=Cx;
%   where C=sum(c(i)*Cp{i});
% C is a M x N
% INPUT:
%   YX: Sufficient statistics <y*x'>  MxN Matrix
%   XX: Sufficient statistics <x*x'>: NxN Matrix
%   Cc: Cell{h,1} array of Components of C matrix (sparse, all have to be MxN)
%   CcCc: Cell array {i,j} arrary of Cc{i}*Cc{j}'
%          To save time only the elements with j>=i are needed
% OUTPUT:
%   COV: Hx1 matrix of covariances
%   VA:  HxH matrix of variances
% Solution is given by VA\COV
%
% WARING: for speed reasons the mex version of this function skips all
% checks
% (c) Joern Diedrichsen 2012 j.diedrichsen@ucl.ac.uk
 *========================================================*/

#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int m,n;            // Sizes of arrays 
    int H;              // Number of parameters 
    int i,j;            // Column and row indices into YX data 
    int h,g;            // Column and row indices intp parameters 
    double *YX;         // YX statistics 
    double *XX;         // XX statistics 
    double *COVd;       // Hx1 martix 
    double *VAd;        // HxH matrix 
    
    mxArray *Cc;        // Pointers to sparse structure matrix 
    mwIndex *Cc_Ir;     // Row indicator of sparse matrix 
    mwIndex *Cc_Jc;     // Jc[j] is the index in Cc_data and Cc_Ir where the jth column starts 
    double *Cc_data;    // Data of sparse matrix 
    
    // check for proper number of arguments 
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:ML_constrained_fast:nrhs","Four inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:ML_constrained_fast:nlhs","Two output required.");
    }
  
    // Get sizes: Be careful: all checks are omitted: will simply crash if not called properly
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    H = mxGetNumberOfElements(prhs[2]); // Number of parameters 
     
    // Get Pointer to data 
    YX = mxGetPr(prhs[0]);
    XX = mxGetPr(prhs[1]);
    
    // Allocate the output data 
    plhs[0] = mxCreateDoubleMatrix(H,1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(H,H, mxREAL); 
    COVd = mxGetPr(plhs[0]);
    VAd = mxGetPr(plhs[1]);

    
    for (h=0;h<H;h++) { 
        Cc=mxGetCell(prhs[2],h);
        Cc_Ir=mxGetIr(Cc); 
        Cc_Jc=mxGetJc(Cc); 
        Cc_data=mxGetPr(Cc); 

        // Get Covariances 
        for (j=0;j<n;j++) { 
            for (i=Cc_Jc[j];i<Cc_Jc[j+1];i++) {
                COVd[h]+=YX[Cc_Ir[i]+j*m]*Cc_data[i];
            } 
        } 

        // Get Variance sum(sum(CcCc{i,i}.*XX))
        Cc=mxGetCell(prhs[3],h+h*H);  // Diagonal element
        Cc_Ir=mxGetIr(Cc); 
        Cc_Jc=mxGetJc(Cc); 
        Cc_data=mxGetPr(Cc); 
        for (j=0;j<n;j++) { 
            for (i=Cc_Jc[j];i<Cc_Jc[j+1];i++) {
                VAd[h+h*H]+=XX[Cc_Ir[i]+j*n]*Cc_data[i];
            } 
        } 

        // Get Variance sum(sum(CcCc{i,j}.*XX))
        // Assume that only CcCc{i,j} j>i are non-empty 
        for (g=h+1;g<H;g++) { 
            Cc=mxGetCell(prhs[3],h+g*H);  // Diagonal element
            Cc_Ir=mxGetIr(Cc); 
            Cc_Jc=mxGetJc(Cc); 
            Cc_data=mxGetPr(Cc); 
            for (j=0;j<n;j++) { 
                for (i=Cc_Jc[j];i<Cc_Jc[j+1];i++) {
                    VAd[h+g*H]+=XX[Cc_Ir[i]+j*n]*Cc_data[i];
                } 
            } 
            VAd[g+h*H]=VAd[h+g*H];
        } 
    }    
}
