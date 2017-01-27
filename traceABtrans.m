%trace of A times B'
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
%
% See also traceAB