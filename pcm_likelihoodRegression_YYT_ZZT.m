function [negLogLike,dnl,d2nl] = pcm_likelihoodRegression_YYT_ZZT(theta,Z,YY,P,comp,X,varargin);
% function [negLogLike,dnl,d2nl] = pcm_likelihoodRegression_YYT_ZTZ(theta,Z,YY,P,comp,X,varargin);
% Returns negative log likelihood for the model, and the derivatives in
% respect to the model parameters. This is version of pcm_likelihood, which
% is designed specifically for determining the ridge coefficient for groups
% of features. 
% This version of the likelihood function is fastest if P>N and Q > N
%
% INPUT:
%      theta:   Vector of (log-)model parameters: These include model
%               parameters, noise parameter, and (optional) run parameter
%      Y:       NxP Matrix of data 
%      Z:       NxK Design matrix - relating the trials (N) to the random effects (K)
%      comp:    1xK vector of groups that determine the groups of features
%               in the design matrix that correspond to a theta coefficient
%      X:       Fixed effects design matrix - will be accounted for by ReML
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
%      dnl      :   Derivative of the negative log-likelihood in respect to
%                   the parameters
%      d2nl     :   Expected second derivative of the negative
%                   log-likelihood
%
%   Joern Diedrichsen & Atsushi Yokoi, 6/2016, joern.diedrichsen@googlemail.com
%

N = size(YY,1);
K = size(Z,2);
S = [];
nComp = max(comp); 
nParam = length(theta); 
pcm_vararginoptions(varargin,{'S'});

% Split the paraemters in noise and ridge parameter 
modelParam = theta(1:nComp); 
noiseParam = theta(nComp+1:end);

% Find the inverse of V 
% Apply the matrix inversion lemma. The following statement is the same as
% V   = (Z*G*Z' + S.S*exp(theta(H))); % As S is not identity, matrix inversion lemma does not have big advantage here (ay)?
% iV  = pinv(V);
G  = exp(modelParam(comp)); 
iG = 1./ G; 
if (isempty(S))
    iV    = (eye(N)-Z/(diag(iG)*exp(noiseParam)+Z'*Z)*Z')./exp(noiseParam); % Matrix inversion lemma
    Zw = bsxfun(@times,Z,sqrt(G)');
    lam = eig(Zw' * Zw); 
    ldet  = sum(log(lam+exp(noiseParam))) + (N - K)*noiseParam;
else
    iV    = (S.invS-S.invS*Z/(diag(diag(iG)*exp(noiseParam)+Z'*S.invS*Z)*Z'*S.invS))./exp(noiseParam); % Matrix inversion lemma
    ldet  = -2* sum(log(diag(chol(iV))));        % Safe computation of the log determinant (V) Thanks to code from D. lu
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
B     = YY * iVr;
l     = -P/2*(ldet)-0.5*trace(B);
if (~isempty(X)) % Correct for ReML estimates
    l = l - P*sum(log(diag(chol(X'*iV*X))));  % - P/2 log(det(X'V^-1*X));
end;
negLogLike = -l; % Invert sign

% Calculate the first derivative
if (nargout>1)
    % Get the derivatives for all the parameters
    for i = 1:length(modelParam)
        iVdV{i} = iVr*Z(:,comp==i)*Z(:,comp==i)'*exp(modelParam(i));
        dLdtheta(i,1) = -P/2*trace(iVdV{i})+1/2*traceABtrans(iVdV{i},B);
    end
    
    % Get the derivatives for the Noise parameters
    i = length(modelParam)+1;  % Which number parameter is it?
    if (isempty(S))
        dVdtheta{i}          = eye(N)*exp(noiseParam);
    else
        dVdtheta{i}          = S.S*exp(noiseParam);
    end;
    iVdV{i}     = iVr*dVdtheta{i};
    dLdtheta(i,1) = -P/2*trace(iVdV{i})+1/2*traceABtrans(iVdV{i},B);
    
    % invert sign
    dnl   = -dLdtheta;
end;

% Calculate expected second derivative?
if (nargout>2)
    for i=1:nParam
        for j=i:nParam
            d2nl(i,j)=-P/2*traceAB(iVdV{i},iVdV{j});
            d2nl(j,i)=d2nl(i,j);
        end;
    end;
    d2nl=-d2nl;
end;
