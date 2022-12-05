function [param] = ML_constrained(s2,s3,Cp)
% function [param] = ML_constrained(s2,s3,Cp)
% constrained maximum likelihood estimation of regression parameters
% y=Cx;
% where C=sum(c(i)*Cp{i});
% INPUT:
%   s2: Sufficient statistics <x*y'>
%   s3: Sufficient statistics <x*x'>
%   Cp: Cell array of Components of C matrix
% OUTPUT: 
%   param: new parameter estimates 
num_bparams = length(Cp);
COV = zeros(num_bparams,1);
V = zeros(num_bparams, num_bparams);
for i = 1:num_bparams
    COV(i,1)=sum(sum(Cp{i}.*s2'));  % use trace-trick: trace(A*B)=sum(sum(A.*B'))
    for j = 1:num_bparams
        V(i,j)=sum(sum((Cp{j}'*Cp{i}).*s3));  % again: trace(A*B)=sum(sum(A.*B'))
    end;
end;
param = V\COV;

