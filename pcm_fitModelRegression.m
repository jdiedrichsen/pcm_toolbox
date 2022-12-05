function [theta_hat,INFO]=pcm_fitModelRegression(Z,Y,comp,X,varargin)
% [theta_hat,INFO]=pcm_fitModelRegression(Z,Y,comp,X,varargin);
% Estimates the hyperparameters for a Tikhonov-regression model
% Will fit max(comp) hyperparameters 
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
% OPTION (VARARGIN):
%   'maxInteration',int : 
%       Maximum number of iterations 
%   'verbose',int: 
%       Verbosity level: 0: no printout 
%   'S',[#Obs x #Obs] matrix: 
%       Covariance structure of noise (default identity) 
%   'theta0', vec: 
%       Vector of starting values for theta 
% OUTPUT: 
%   theta_hat [max(comp)+1  x 1 Vector]: 
%       Vector of estimated hyperparamters. All are log(var). Theta_hat(1..K-1) 
%       are the log(var(regressor group)). The last one log(var(noise))
%       The optimal regualrization coeffecient for the group i is 
%       exp(theta_hat(end))/exp(theta_hat(i))
%   INFO: Structure
%       Information structure of the fit. 
%       INFO.like: log-Likelihood
%       INFO.iter: Number of iterations 
%       INFO.thetaH: History of parameters across iterations 
%       INFO.likeH: History of likelihood across iterations 
maxIteration    = [];      % MaxIteration - otherwise algorithm-specific defaults 
likeThres       = [];    % if log-likelihood decreases less that this value, it stops (for pcm_NR)
verbose         = 1;   % 1: Indicating the subject 2: detailed feedback
S               = [];
theta0          = [];

pcm_vararginoptions(varargin,{'','isCheckDeriv','maxIteration',...
    'verbose','S','fitAlgorithm','theta0','fitScale','likeThres'});

% Set options for fitting algorithm 
fitOPT.verbose = verbose; 
if (~isempty(maxIteration))
    fitOPT.maxIteration = maxIteration;
end 
if (~isempty(likeThres))
    fitOPT.likeThres = likeThres; 
end


% Ensure that X is full rank
[U,SX,V] = svd(X,0); 
X = U(:,find(diag(SX)>eps)); 

% Get starting values (very imortant for noise variance) 
nParam = max(comp)+1;
if (isempty(theta0)) 
    theta0 = zeros(nParam,1);
    % Get noise variance approximately into the right size 
    [N,P] = size(Y); 
    R=Y-X*pinv(X)*Y;
    noiseVar = sum(sum(R.*R,1))/(N*P); 
    theta0(nParam,1) = log(noiseVar);
end


% Optimize the hyperparameters
fcn = @(x) pcm_likelihoodRegression_YTY_ZTZ(x,Z,Y,comp,X,'S',S);
[theta_hat,INFO.like,INFO.iter,INFO.regFinal,INFO.regH,INFO.thetaH, INFO.likeH] = ...
      pcm_NR(theta0,fcn,fitOPT); 