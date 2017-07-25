function [G,dGdtheta] = ra_modelpred_addsc(theta)
% function [G,dGdtheta] = ra_modelpred_addsc(theta)
% This is an example of a nonlinear pattern component model 
% For an experiment with 5 fingers, measured at 4 different tapping
% speeds.
% The first 14 paramters determine the stucture of the Finger patterns
% encoded in A. With OM = A*A'. 
% Then there are 3 different scaling parameters for the 1-3 speed. The
% scaling parameter for the 4th speed is fixed to one. 
% Finally, there are 3 different additive parameters, also one per speed.
% Again, the additive parameter for the 4th speed is fixed to one.
% So the prediction of the ith finger for the jth speed is 
% y_i,j = s_j * f_i + s_q
% Where f_i is a full pattern and s_j a scalar and s_q an additive
% background pattern.
% It then calculates the derivative of the G matrix in respect to the
% parameters. 

fingerParams   = [1;theta(1:14)];
scaleParams    = exp(theta(15:17));
addParams      = exp(theta(18:20));
indx1          = double(tril(true(5),0));
indx1(indx1>0) = [1:15];
indx2          = indx1';
A              = zeros(5);
M              = [kron([scaleParams;1],eye(5)) kron([0;addParams],ones(5,1))];
% Finger parameters
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

% % Build derivatives
for i = 1:15 
    if (i>1)
        dOM               = dA(:,:,i)*A'+A*dA(:,:,i)'; 
        dGdtheta(:,:,i-1) = M*dOM*M';
    end; 
end; 
% Scale parameters 
for i = 1:3
    sc                 = zeros(4,1);
    sc(i)              = 1;
    dM                 = [kron(sc,eye(5)) zeros(20,1)];
    dGdtheta(:,:,14+i) = dM*OM*M'+M*OM*dM';
    dGdtheta(:,:,14+i) = dGdtheta(:,:,14+i)*scaleParams(i);
end;
% Additive parameters 
for i = 1:3
    add                = zeros(4,1);
    add(i+1)           = 1;
    dM                 = [zeros(20,5) kron(add,ones(5,1))];
    dGdtheta(:,:,17+i) = dM*OM*M'+M*OM*dM';
    dGdtheta(:,:,17+i) = dGdtheta(:,:,17+i)*addParams(i);
end;