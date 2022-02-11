function [U,varU] = pcm_estimateU(M,theta,Y,Z,X,varargin);
% function [U,G] = pcm_estimateU(M,theta,Y,M,Z,X,varargin);
% Returns the random effects estimate of a PCM model, using the current set
% of parameters of the model 
%
% INPUT:
%       M:      Model specification. Model can be either given by
%               b. Model Structure containing the fields
%                  Model.type:        fixed, component, feature, nonlinear
%                  Model.... 
%      theta:   Vector of (log-)model parameters: These include model
%               parameters, noise parameter, and (option) run parameter
%      Y:       NxP Matrix of data
%      Z:       NxK Design matrix - relating the trials (N) to the random effects (K)
%      X:       Fixed effects design matrix - subtracted out before 
% VARARGIN:
%      'runEffect',B:  design matrice for the run effect,
%                   which is modelled as a individual subject-specific random effect.
%      'S':    Explicit noise covariance matrix structure matrix. The For speed,
%              this is a cell array that contains
%              S.S:     Structure of noise
%              S.invS:  inverse of the noise covariance matrix
%              if empty, this defaults to the indentity matrix 
% OUTPUT:
%      U:     KxP matrix   BLUP estimates of the activity patterns for the K 
%                          experimental conditions.
%      varU:  KxK matrix   Variance-covariance of the estimates for columns U 
%   Joern Diedrichsen 3/2017, joern.diedrichsen@googlemail.com

[N,P] = size(Y);
K = size(Z,2);
S = [];
runEffect =[];
pcm_vararginoptions(varargin,{'S','runEffect'});

% Get G-matrix and derivative of G-matrix in respect to parameters
if (isstruct(M))
    G = pcm_calculateG(M,theta(1:M.numGparams));
else
    G=M;
    M=[];
    M.numGparams=0;
end

% If Run effect is to ne modelled as a random effect - add to G and
% design matrix
noiseParam = theta(M.numGparams+1);

if (~isempty(runEffect))
    numRuns = size(runEffect,2);
    runParam = theta(M.numGparams+2);    % Subject run effect parameter
    G = pcm_blockdiag(G,eye(numRuns)*exp(runParam));  % Include run effect in G
    Z = [Z runEffect];                 % Include run effect in design matrix
else
    numRuns = 0;                                % No run effects modelled
end

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
end

% For ReML, compute the modified inverse iVr
if (~isempty(X))
    iVX   = iV * X;
    iVr   = iV - iVX*((X'*iVX)\iVX');
else
    iVr   = iV;
end

% Compute the random effects 
U=G*Z'*iVr*Y;
varU = G*Z'*iVr*Z*G; 
if (~isempty(runEffect)) 
    U=U(1:K,:); 
    varU=varU(1:K,1:K); 
end
