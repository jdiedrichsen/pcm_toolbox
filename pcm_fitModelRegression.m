function [theta_hat,INFO]=pcm_fitModelRegression(Z,Y,comp,X,varargin);
% [theta_hat,G_pred,INFO]=pcm_fitModelRegression(Z,Y,comp,varargin);
% Fits a Tikhonov-regression model to the data Y 

runEffect       = 'fixed';
isCheckDeriv    = 0;
maxIteration    = [];      % MaxIteration - otherwise algorithm-specific defaults 
likeThres       = [];    % if log-likelihood decreases less that this value, it stops (for pcm_NR)
verbose         = 1;   % 1: Indicating the subject 2: detailed feedback
S               = [];
fitAlgorithm    = [];
theta0          = [];
fitScale        = 0; 
pcm_vararginoptions(varargin,{'runEffect','isCheckDeriv','maxIteration',...
    'verbose','S','fitAlgorithm','theta0','fitScale','likeThres'});

% Set options for fitting algorithm 
fitOPT.verbose = verbose; 
if (~isempty(maxIteration))
    fitOPT.maxIteration = maxIteration;
end 
if (~isempty(likeThres))
    fitOPT.likeThres = likeThres; 
end

nParam = max(comp)+1;
theta0 = zeros(nParam,1); 
fcn = @(x) pcm_likelihoodRegression(x,Z,Y,comp,X,'S',S);
[theta_hat,INFO.like,INFO.iter,INFO.regFinal,INFO.regH,INFO.thetaH]= pcm_NR(theta0,fcn,fitOPT); 