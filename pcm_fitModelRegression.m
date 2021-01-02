function [theta_hat,INFO]=pcm_fitModelRegression(Z,Y,comp,X,varargin);
% [theta_hat,G_pred,INFO]=pcm_fitModelRegression(Z,Y,comp,varargin);
% Fits a Tikhonov-regression model to the data Y 

likeFcn         = 'auto';
maxIteration    = [];      % MaxIteration - otherwise algorithm-specific defaults 
likeThres       = [];    % if log-likelihood decreases less that this value, it stops (for pcm_NR)
verbose         = 1;   % 1: Indicating the subject 2: detailed feedback
S               = [];
fitAlgorithm    = 'newton';
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
tic; 
fcn = @(x) pcm_likelihoodRegression_YTY_ZTZ(x,Z,Y,comp,X,'S',S);
[theta_hat,INFO.like,INFO.iter,INFO.regFinal,INFO.regH,INFO.thetaH, INFO.likeH] = ...
      pcm_NR(theta0,fcn,fitOPT); 
INFO.time=toc; 