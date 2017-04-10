function [G,dGdtheta] = ra_modelpred_scale(theta)
% function [G,dGdtheta] = ra_modelpred_scale(theta)
% This is an example of a nonlinear pattern component model 
% For an experiment with 5 fingers, measured at 4 different tapping
% speeds.
% The first 14 paramters determine the stucture of the Finger patterns
% encoded in A. With OM = A*A'. 
% Then there are 3 different scaling parameters for the 1-3 speed. The
% scaling parameter for the 4th speed is fixed to one. 
% So the prediction of the ith finger for the jth speed is 
% y_i,j = s_j * f_i 
% Where f_i is a full pattern and s_j a scalar.
% It then calculates the derivative of the G matrix in respect to the
% parameters. 
fingerParams   = theta(1:14);
scaleParams    = exp(theta(15:17));
indx1          = double(tril(true(5),0));
indx1(indx1>0) = [-1 1:14];
indx2          = indx1';
A              = zeros(5);
A(indx1==-1)   = 1; 
M              = kron([scaleParams;1],eye(5));
for i = 1:14
    A(indx1==i | indx2==i) = fingerParams(i); 
end; 
OM = A*A';          
G  = M*OM*M';  % Second moment matrix

% Build dertivatives 
for i=1:14
    dA  = zeros(5);
    dA(indx1==i | indx2==i) = 1; 
    dOM = dA*A'+A*dA;
    dGdtheta(:,:,i) = M*dOM*M';   
end; 

for i = 1:3
    sc    = zeros(4,1);
    sc(i) = 1;
    dM    = kron(sc,eye(5));
    dGdtheta(:,:,14+i) = dM*OM*M'+M*OM*dM'; % derivative for tapping speed i
    dGdtheta(:,:,14+i) = dGdtheta(:,:,14+i)*scaleParams(i); % scaled derivative 
end;
