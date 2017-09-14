function [negLogLike,dnl,d2nl] = pcm_likelihoodIndivid(theta,YY,M,Z,X,P,varargin);
% function [negLogLike,dnl,d2nl] = pcm_likelihoodIndivid(theta,YY,M,Z,X,P,varargin);
% Returns negative log likelihood for the model, and the derivatives in
% respect to the model parameters for an individual subject / dataset
% It works very similar to pcm_likelihood, only that a run-effect can be
% specified to be modelled as an extra random effect
%
% INPUT:
%      theta:   Vector of (log-)model parameters: These include model
%               parameters, noise parameter, and (optional) run parameter
%      YY:      NxN Matrix of outer product of the activity data (Y*Y')
%      M:       Model structure with fields
%                  Model.type:        fixed, component, feature, nonlinear
%                  Model.numGparams:  Number of model parameters (without the noise or run parameter)
%                  ...
%      Z:       NxK Design matrix - relating the trials (N) to the random effects (K)
%      X:       Fixed effects design matrix - will be accounted for by ReML
%      P:       Number of voxels
% VARARGIN:
%      'runEffect',B:  design matrice for the run effect,
%                   which is modelled as a individual subject-specific random effect.
%      'S':    Explicit noise covariance matrix structure matrix. The For speed,
%              this is a cell array that contains
%              S.S:     Structure of noise
%              S.invS:  inverse of the noise covariance matrix
%       'calcLikelihood': If set to zero, will skip calculating the
%              log-liklihood
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
OPT.S = [];
OPT.runEffect =[];

OPT=pcm_vararginoptions(varargin,OPT,{'S','runEffect'});


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

if (~isempty(OPT.runEffect))
    numRuns = size(OPT.runEffect,2);
    runParam = theta(M.numGparams+2);    % Subject run effect parameter
    G = pcm_blockdiag(G,eye(numRuns)*exp(runParam));  % Include run effect in G
    Z = [Z OPT.runEffect];                 % Include run effect in design matrix
else
    numRuns = 0;                                % No run effects modelled
end;


% Find the inverse of V - while dropping the zero dimensions in G
[u,s] = eig(G);
dS    = diag(s);
idx   = dS>eps;
Zu     = Z*u(:,idx);
% Apply the matrix inversion lemma. The following statement is the same as
% V   = (Z*G*Z' + S.S*exp(theta(H)));
% iV  = pinv(V);
if (isempty(OPT.S))
    iV    = (eye(N)-Zu/(diag(1./dS(idx))*exp(noiseParam)+Zu'*Zu)*Zu')./exp(noiseParam); % Matrix inversion lemma
else
    iV    = (OPT.S.invS-OPT.S.invS*Zu/(diag(1./dS(idx))*exp(noiseParam)+Zu'*OPT.S.invS*Zu)*Zu'*OPT.S.invS)./exp(noiseParam); % Matrix inversion lemma
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
    negLogLike = -l; % Invert sign


% Calculate the first derivative
if (nargout>1)
    A     = iVr*Z;      % Precompute some matrices
    B     = YY*iVr;
    % Get the derivatives for all the parameters
    for i = 1:M.numGparams
        iVdV{i} = A*pcm_blockdiag(dGdtheta(:,:,i),zeros(numRuns))*Z';
        dLdtheta(i,1) = -P/2*trace(iVdV{i})+1/2*traceABtrans(iVdV{i},B);
    end
    
    % Get the derivatives for the Noise parameters
    i             = M.numGparams+1;  % Which number parameter is it?
    if (isempty(OPT.S))
        dVdtheta{i}          = eye(N)*exp(noiseParam);
    else
        dVdtheta{i}          = OPT.S.S*exp(noiseParam);
    end;
    iVdV{i}     = iVr*dVdtheta{i};
    dLdtheta(i,1) = -P/2*trace(iVdV{i})+1/2*traceABtrans(iVdV{i},B);
    
    % Get the derivatives for the block parameter
    if (~isempty(OPT.runEffect) && ~isempty(OPT.runEffect))
        i = M.numGparams+2;  % Which number parameter is it?
        %C          = A*pcm_blockdiag(zeros(size(G,1)),eye(numRuns));
        iVdV{i}     = A*pcm_blockdiag(zeros(K),eye(numRuns))*Z'*exp(runParam);
        dLdtheta(i,1) = -P/2*trace(iVdV{i})+1/2*traceABtrans(iVdV{i},B);
    end;
    
    % invert sign
    dnl   = -dLdtheta;
    numTheta=i;
end;

% Calculate expected second derivative?
if (nargout>2)
    for i=1:numTheta
        for j=i:numTheta;
            d2nl(i,j)=-P/2*traceABtrans(iVdV{i},iVdV{j});
            d2nl(j,i)=d2nl(i,j);
        end;
    end;
    d2nl=-d2nl;
end;

