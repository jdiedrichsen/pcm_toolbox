function [U, varU] = pcm_estimateRegression(Z,Y,comp,X,theta,varargin)
% function [U, varU] = pcm_estimateRegression(Z,Y,comp,X,theta,varargin)
% [theta_hat,G_pred,INFO]=pcm_fitModelRegression(Z,Y,comp,varargin);
% Gets estimates (U) for a regression model with hyperparameters estimated. 
% INPUT:
%   Z: [Matrix #Observations x #Regressors]
%       Designmatrix for the random effects (the ones that are regualized)
%   Y: [Matrix #Observations x #Voxels]
%       Observed/estimated beta regressors from one subject.
%       Preferably multivariate noise-normalized beta regressors.
%   comp: [Vector #Regressor] 
%       Indicator (1..K) of which regressor belongs to which regressor 
%       group 
%   X: [Matrix #Observations x #fixedRegressors]
%       Design matrix of fixed effects (not regularized). Usually this
%       is the intercept or intercept for each run. 
%   theta_hat [max(comp)+1  x 1 Vector]: 
%       Vector of estimated hyperparamters. 
% OPTION (VARARGIN):
%   'S',[#Obs x #Obs] matrix: 
%       Covariance structure of noise (default identity) 
%  RESULTS: 
%   U: [Matrix #Regressors x # voxels] 
%       Estimates for Random effects (BLUPs or regularized regressor
%       coefficients)
%   varU: [Matrix #Regressors x # Regressors]
%       Variance - covariance matrix of the regressors 

S               = [];
pcm_vararginoptions(varargin,{'S'});

% Split the paraemters in noise and ridge parameter
nComp = max(comp);
nParam = length(theta);
modelParam = theta(1:nComp);
noiseParam = theta(nComp+1:end);
[N,P] = size(Y); 

% Generate the G-matrix
G  = exp(theta(comp));
iG = 1./ G;

% WAY 1: Find the inverse of V
% Apply the matrix inversion lemma. The following statement is the same as
% V   = (Z*diag(G)*Z' + eye(N)*exp(noiseParam)); % As S is not identity, matrix inversion lemma does not have big advantage here (ay)?
% iV  = pinv(V);
if (isempty(S))  
    if min(theta(comp))<-20
        iV = pinv((Z*diag(G)*Z' + eye(N)*exp(noiseParam)));
    else
        iV    = (eye(N)-Z/(diag(iG)*exp(noiseParam)+Z'*Z)*Z')./exp(noiseParam); % Matrix inversion lemma
    end
else
    if min(theta(comp))<-20
        iV = pinv((Z*diag(G)*Z' + S.S*exp(noiseParam)))
    else
        iV    = (S.invS-S.invS*Z/(diag(diag(iG)*exp(noiseParam)+Z'*S.invS*Z)*Z'*S.invS))./exp(noiseParam); % Matrix inversion lemma
    end
end

% For ReML, compute the modified inverse iVr
if (~isempty(X))
    % Ensure that X is full rank
    [U,SX,V] = svd(X,0); 
    X = U(:,(diag(SX)>eps)); 
    iVX   = iV * X;
    iVr   = iV - iVX*((X'*iVX)\iVX');
else
    iVr   = iV;
end

% Compute the random effects
GZt = bsxfun(@times,Z',G); % G*Z' 
U=GZt*iVr*Y;
if (nargout>1)
    varU = GZt*iVr*GZt'; 
end

% WAY 2: Do it over the ridge regression way 
% This is identical 
% (By matrix inversion lemma) 
% if (~isempty(X))
%     R = eye(N)-X*pinv(X); 
%     Z = R * Z; 
%     Y = R * Y; 
% end
% U = (Z' * Z + exp(noiseParam).*diag(iG))\(Z'*Y);
