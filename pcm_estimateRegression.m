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

% Generate the G-matrix
G  = exp(modelParam(comp));
iG = 1./ G;

% WAY 1: Find the inverse of V
% Apply the matrix inversion lemma. The following statement is the same as
% V   = (Z*G*Z' + S.S*exp(theta(H))); % As S is not identity, matrix inversion lemma does not have big advantage here (ay)?
% iV  = pinv(V);
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
GZt = bsxfun(@times,Z',G); % G*Z' 
U=GZt*iVr*Y;
if (nargout>1)
    varU = GZt*iVr*GZt'; 
end

% WAY 2: Do it over the ridge regression way 
% (By matrix inversion lemma) 
% if (~isempty(X))
%     R = eye(N)-X*pinv(X); 
%     Z = R * Z; 
%     Y = R * Y; 
% end
% U = (Z' * Z + exp(noiseParam).*diag(iG))\(Z'*Y);
