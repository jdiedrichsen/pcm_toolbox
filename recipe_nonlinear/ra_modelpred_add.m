function [G,dGdtheta] = ra_modelpred_add(theta)
% function [G,dGdtheta] = ra_modelpred_add(theta)
% This is an example of a nonlinear pattern component model 
% For an experiment with 5 fingers, measured at 4 different tapping
% speeds.
% The first 14 paramters determine the stucture of the Finger patterns
% encoded in A. With OM = A*A'. 
% Then there are 3 different additive parameters for the 1-3 speed. The
% additive parameter for the 4th speed is fixed to one. 
% So the prediction of the ith finger for the jth speed is 
% y_i,j = f_i + s_q
% Where f_i is a full pattern and s_q a background additive pattern.
% It then calculates the derivative of the G matrix in respect to the
% parameters. 
fingerParams   = [1;theta(1:14)];
addParams      = exp(theta(15:17));
indx1          = double(tril(true(5),0));
indx1(indx1>0) = [1:15];
indx2          = indx1';
M              = [kron(ones(4,1),eye(5)) kron([0;addParams],ones(5,1))];

A = zeros(6);
for i = 1:15
    d                      = zeros(5);
    d(indx1==i | indx2==i) = 1;
    d(6,6)                 = 0; 
    dA(:,:,i)              = d; 
    A                      = A+dA(:,:,i)*fingerParams(i); 
end;
A(6,6) = 1; 
OM     = A*A'; 
G      = M*OM*M';  % Second moment matrix

% % Build dertivatives 
% finger params
for i = 1:15 
    if (i>1)
        dOM  = dA(:,:,i)*A'+A*dA(:,:,i)'; 
        dGdtheta(:,:,i-1) = M*dOM*M';
    end; 
end; 
% additive params
for i=1:3
    add                = zeros(4,1);
    add(i+1)           = 1;
    dM                 = [zeros(20,5) kron(add,ones(5,1))];
    dGdtheta(:,:,14+i) = dM*OM*M'+M*OM*dM'; % derivative for tapping speed i
    dGdtheta(:,:,14+i) = dGdtheta(:,:,14+i)*addParams(i); % scaled derivative 
end;
