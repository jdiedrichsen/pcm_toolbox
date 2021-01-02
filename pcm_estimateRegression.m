function [U, varU] = pcm_estimateRegression(Z,Y,comp,X,theta,varargin);
% [theta_hat,G_pred,INFO]=pcm_fitModelRegression(Z,Y,comp,varargin);
% Gets estimates for a regression model with estimated theta 

S               = [];

pcm_vararginoptions(varargin,{'S'});

% Split the paraemters in noise and ridge parameter
nComp = max(comp);
nParam = length(theta);
modelParam = theta(1:nComp);
noiseParam = theta(nComp+1:end);
[N,P] = size(Y); 

% Find the inverse of V
% Apply the matrix inversion lemma. The following statement is the same as
% V   = (Z*G*Z' + S.S*exp(theta(H))); % As S is not identity, matrix inversion lemma does not have big advantage here (ay)?
% iV  = pinv(V);
G  = exp(modelParam(comp));
iG = 1./ G;
if (isempty(S))
    iV    = (eye(N)-Z/(diag(iG)*exp(noiseParam)+Z'*Z)*Z')./exp(noiseParam); % Matrix inversion lemma
else
    iV    = (S.invS-S.invS*Z/(diag(diag(iG)*exp(noiseParam)+Z'*S.invS*Z)*Z'*S.invS))./exp(noiseParam); % Matrix inversion lemma
end

% For ReML, compute the modified inverse iVr
if (~isempty(X))
    % Ensure that X is full rank
    [U,SX,V] = svd(X,0); 
    X = U(:,find(diag(SX)>eps)); 
    iVX   = iV * X;
    iVr   = iV - iVX*((X'*iVX)\iVX');
else
    iVr   = iV;
end

% Compute the random effects 
U=G*Z'*iVr*Y;
if (nargout>1)
    varU = G*Z'*iVr*Z*G; 
end