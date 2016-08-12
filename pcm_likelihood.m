function [negLogLike,dnldtheta] = pcm_likelihood(theta,YY,Model,Z,X,P,varargin);
% function [negLogLike,dnldtheta,L,dLdtheta] = pcm_likelihood(theta,YY,Gc,Z,X,P,varargin);
% Returns negative log likelihood for the model, and the derivatives in
% respect to the model parameters. 
%   Y = Z U + X B + E 
%       U : Random effect ~ N(0,G) 
%       B : Fixed effects 
%       E : Random noise ~ N(0,S*sigma2)  
%
% INPUT:
%      theta:   Vector of (log-)model parameters
%      YY:      NxN Matrix of data
%      Model:   Model specification. Model can be either given by 
%               a. KxKxH Tensor of component matrices: G = sum (log(theta_i) * Model(:,:,i)); 
%               b. Model Structure containing the fields 
%                  Model.numGparams:  Number of model parameters (without the noise scaling) 
%                  Model.modelpred:   Function handle, returning [G,dGdtheta]= modelpred(theta); 
%      Z:       NxK Design matrix - relating the trials (N) to the random effects (K) 
%      X:       Fixed effects design matrix - will be accounted for by ReML
%      P:       Number of voxels
% VARARGIN:
%      'S':    Explicit noise covariance matrix structure matrix. The For speed, 
%              this is a cell array that contains
%              S.S:     Structure of noise 
%              S.invS:  inverse of the noise covariance matrix 
% OUTPUT:
%      negLogLike:  Negative Log likelihood of all subject's data
%                   We use the negative log liklihood to be able to plug the function into
%                   minimize or other optimisation routines.
%      dnldtheta:   Derivative of the negtive log-likelihood
%      L:           Log likelihood (not inverted) for all the subject
%      dLdtheta:    Derivate of Log likelihood for each subject
%
%   Joern Diedrichsen & Atsushi Yokoi, 6/2016, joern.diedrichsen@googlemail.com
%

N = size(YY,1);
K = size(Z,2); 
S = []; 

vararginoptions(varargin,{'S'}); 


if (isstruct(Model))
    H      = Model.numGparams+1; 
    [G,Gc] = Model.fcn(theta(1:M.numParams)); 
else
    Gc    =  Model; 
    H     =  size(Gc,3)+1;       % Number of Hyperparameters (+ 1 for noise term)
    G     =  zeros(K);
    for i = 1:H-1;
        G = G + Gc(:,:,i)*exp(theta(i));
    end
end; 

% Find the inverse of V - while dropping the zero dimensions in G
[u,s] = eig(G);
dS    = diag(s);
idx   = dS>eps;
Zu     = Z*u(:,idx);
% Apply the matrix inversion lemma. The following statement is the same as 
% V   = (Z*G*Z' + S.S*exp(theta(H))); % As S is not identity, matrix inversion lemma does not have big advantage here (ay)?
% iV  = pinv(V);   
if (isempty(S))
    iV    = (eye(N)-Zu/(diag(1./dS(idx))*exp(theta(H))+Zu'*Zu)*Zu')./exp(theta(H)); % Matrix inversion lemma
else
    iV    = (S.invS-S.invS*Zu/(diag(1./dS(idx))*exp(theta(H))+Zu'*S.invS*Zu)*Zu'*S.invS)./exp(theta(H)); % Matrix inversion lemma
end; 

% For ReML, compute the modified inverse iVr
if (~isempty(X))
    iVX   = iV * X;
    iVr   = iV - iVX*((X'*iVX)\iVX');  
else 
    iVr   = iV; 
end;

% Computation of (restricted) likelihood 
ldet  = -2* sum(log(diag(chol(iV))));        % Safe computation of the log determinant (V) Thanks to code from D. lu
l     = -P/2*(ldet)-0.5*traceABtrans(iVr,YY);
if (~isempty(X)) % Correct for ReML estimates
    l = l - P*sum(log(diag(chol(X'*iV*X))));  % - P/2 log(det(X'V^-1*X));
end;

% Precompute some matrices
A     = iVr*Z;
B     = YY*A/P;

% Get the derivatives for all the parameters
for h = 1:H-1
    C  = (A*Gc(:,:,i));
    dLdtheta(h,1) = -P/2*(traceABtrans(C,Z)-traceABtrans(C,B));
    dLdtheta(h,1) = dLdtheta(h,1)*exp(theta(h));
end

% Get the derivatives for the Noise parameters
if (isempty(S))
    dLdtheta(H,1)     = -P/2*traceABtrans(iVr,(speye(N)-YY*iVr/P));
else 
    dLdtheta(H,1)     = -P/2*traceABtrans(iVr*S.S,(speye(N)-YY*iVr/P));
end; 
dLdtheta(H,1) = dLdtheta(H,1)*exp(theta(H));

% invert sign 
negLogLike = -l;
dnldtheta   = -dLdtheta;
