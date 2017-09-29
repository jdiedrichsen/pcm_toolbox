function [W,A] = pcm_estimateW(M,theta,Y,Z,X,varargin);
% function [W,A] = pcm_estimateW(M,theta,Y,Z,X,varargin);
% Returns the voxel-feature weights for a PCM feature model. It calculates the BLUP 
% estimate of of W, given the best current model fit. 
%
% INPUT:
%       M:      Model specification. Model must be of type='feature'; 
%               The function also returns the best feature design matrix under the current model  
%      theta:   Vector of (log-)model parameters: These include model
%               parameters, noise parameter, and (option) run parameter
%      Y:       NxP Matrix of data
%      Z:       NxK Design matrix - relating the trials (N) to conditions (K)
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
%      W:     QxP matrix. Voxel-feature weights for each of the a flexible
%             feature model. 
%      A:     Feature matrix for the optimal model fit. This is
%             sum(Ac_i*theta_i); 
%   Joern Diedrichsen 3/2017, joern.diedrichsen@googlemail.com
%

[N,P] = size(Y);
K = size(Z,2);
S = [];
runEffect =[];
pcm_vararginoptions(varargin,{'S','runEffect'});

% Get G-matrix and derivative of G-matrix in respect to parameters
if (isstruct(M))
    if (~strcmp(M.type,'feature'))
        error('voxel-feature weights can only be estiamted for feature models'); 
    end; 
    A = bsxfun(@times,M.Ac,permute(theta(1:M.numGparams),[3 2 1]));
    A = sum(A,3); 
else
    A=M;
    M=[];
    M.numGparams=0;
end;
G = A*A'; 

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
    iV    = (eye(N)-Zu/(diag(1./dS(idx))*exp(noiseParam)+Zu'*Zu)*Zu')./exp(noiseParam); % Matrix inversion lemma
else
    iV    = (S.invS-S.invS*Zu/(diag(1./dS(idx))*exp(noiseParam)+Zu'*S.invS*Zu)*Zu'*S.invS)./exp(noiseParam); % Matrix inversion lemma
end;

% For ReML, compute the modified inverse iVr
if (~isempty(X))
    iVX   = iV * X;
    iVr   = iV - iVX*((X'*iVX)\iVX');
else
    iVr   = iV;
end;

% Compute the random effects 
W=A'*Z'*iVr*Y;
