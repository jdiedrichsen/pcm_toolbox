function [negLogLike,dnldtheta] = pcm_likelihood(theta,YY,M,Z,X,P,varargin);
% function [negLogLike,dnldtheta,L,dLdtheta] = pcm_likelihood(theta,YY,M,Z,X,P,varargin);
% Returns negative log likelihood for the model, and the derivatives in
% respect to the model parameters for an individual subject / dataset 
%
% INPUT:
%      theta:   Vector of (log-)model parameters: These include model
%               parameters, noise parameter, and (option) run parameter
%      YY:      NxN Matrix of inner product of the data
%      M:       Model specification. Structures with fields
%                  Model.type:        fixed, component, feature, nonlinear
%                  Model.numGparams:  Number of model parameters (without the noise or run parameter)
%                  ... 
%      Z:       NxK Design matrix - relating the trials (N) to the random effects (K)
%      X:       Fixed effects design matrix - will be accounted for by ReML
%      P:       Number of voxels
% VARARGIN:
%      'S':    Explicit noise covariance matrix structure matrix. The For speed,
%              this is a cell array that contains
%              S.S:     Structure of noise
%              S.invS:  inverse of the noise covariance matrix
% 
% OUTPUT:
%      negLogLike:  Negative Log likelihood of all subject's data
%                   We use the negative log liklihood to be able to plug the function into
%                   minimize or other optimisation routines.
%      dnldtheta:   Derivative of the negative log-likelihood
%      L:           Log likelihood (not inverted) for all the subject
%      dLdtheta:    Derivate of Log likelihood for each subject
%
%   Joern Diedrichsen & Atsushi Yokoi, 6/2016, joern.diedrichsen@googlemail.com
%

N = size(YY,1);
K = size(Z,2);
S = [];

pcm_vararginoptions(varargin,{'S'});


% Get G-matrix and derivative of G-matrix in respect to parameters
if (isstruct(M))
    [G,dGdtheta] = pcm_calculateG(M,theta(1:M.numGparams));
else
    G=M;
    M=[];
    M.numGparams=0;
end;

% If Run effect is to ne modelled as a random effect - add to G and
% design matrix
noiseParam = theta(M.numGparams+1);

% Find the inverse of V - while dropping the zero dimensions in G
[u,s] = eig(G);
dS    = diag(s);
idx   = dS>eps;
Zu     = Z*u(:,idx);
% Apply the matrix inversion lemma. The following statement is the same as
% V   = (Z*G*Z' + S.S*exp(theta(H))); % As S is not identity, matrix inversion lemma does not have big advantage here (ay)?
% iV  = pinv(V);
if (isempty(S))
    iV    = (eye(N)-Zu/(diag(1./dS(idx))*exp(noiseParam)+Zu'*Zu)*Zu')./exp(noiseParam); % Matrix inversion lemma
else
    iV    = (S.invS-S.invS*Zu/(diag(1./dS(idx))*exp(noiseParam)+Zu'*S.invS*Zu)*Zu'*S.invS)./exp(noiseParam); % Matrix inversion lemma
end;
iV  = real(iV); % sometimes iV gets complex

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
for i = 1:M.numGparams
    C  = (A*dGdtheta(:,:,i));
    dLdtheta(i,1) = -P/2*(traceABtrans(C,Z)-traceABtrans(C,B));
end

% Get the derivatives for the Noise parameters
indx             = M.numGparams+1;  % Which number parameter is it?
if (isempty(S))
    dLdtheta(indx,1)     = -P/2*traceABtrans(iVr,(speye(N)-YY*iVr/P))*exp(noiseParam);
else
    dLdtheta(indx,1)     = -P/2*traceABtrans(iVr*S.S,(speye(N)-YY*iVr/P))*exp(noiseParam);
end;

% invert sign
negLogLike = -l;
dnldtheta   = -dLdtheta;
