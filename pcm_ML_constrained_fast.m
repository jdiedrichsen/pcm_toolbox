function [COV,VA]=ML_constrained_fast(YX,XX,Cc,CcCc);
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
H=length(Cc);
COV = zeros(H,1);
VA = zeros(H, H);
for i = 1:H
    COV(i,1)=sum(sum(Cc{i}.*YX));           % trace-trick
    VA(i,i)=sum(sum(CcCc{i,i}.*XX));        % trace-trick
    for j = i+1:H                           % off-diagonal terms
        VA(i,j)=sum(sum(CcCc{i,j}.*XX));    % trace-trick
        VA(j,i)=VA(i,j);                    % asumme symmetry
    end;
end;
