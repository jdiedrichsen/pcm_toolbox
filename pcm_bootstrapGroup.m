function [theta_boot,theta_g]=pcm_bootstrapGroup(Y,M,partitionVec,conditionVec,varargin);
% function [T,theta_hat,G_pred]=pcm_bootstrapGroup(Y,M,partitionVec,conditionVec,varargin);
% Bootstraps the fit of a single group model across subjects. For details
% in inout parameters see pcm_fitModelGroup. 
%==========================================================================
% INPUT:
%        Y: {#Subjects}.[#Conditions x #Voxels]
%            Observed/estimated beta regressors from each subject.
%            Preferably multivariate noise-normalized beta regressors.
%
%        M: a single model(s) to be fitted 
%            define multiple models if desired. Contains subfields:
%              .type:        Type of the model to be fitted
%                             'fixed': Fixed structure without parameter
%                                (except scale and noise for each subject)
%                             'component': G is a sum of linear components, specified by Gc
%                             'feature': G=A*A', with A a linear sum of weighted components Ac
%                             'nonlinear': Nonlinear model with own function
%                                 to return G-matrix and derivatives
%                             'noiseceiling': Uses the mean estimated
%                                 G-matrix from crossvalidation 
%                             'freechol':  Free model in Cholesky form 
%              .numGparams:  Scalar that defines the number of parameters
%                             included in model.
%              .theta0:      Vector of starting values for theta. If not given,
%                              the function attempts to estimate these from a
%                              crossvalidated version Values usually estimated from
%                             observed second-moment matrix. Can estimate
%                             these parameters using 'pcm_modelpred_free_startingval'
%         for more fields see the manual for model specification.                   

%
%   partitionVec: {#Subjects} Cell array with partition assignment vector
%                   for each subject. Rows of partitionVec{subj} define
%                   partition assignment of rows of Y{subj}.
%                   Commonly these are the scanning run #s for beta
%                   regressors. If 1 vector is given, then partitions are 
%                   assumed to be the same across subjects 
%
%   conditionVec: {#Subjects} Cell array with condition assignment vector
%                   for each subject. Rows of conditionVec{subj} define
%                   condition assignment of rows of Y{subj}. If 1 vector is given, 
%                   then conditions are  assumed to be the same across subjects
%                   If the (elements of) conditionVec are matrices, it is
%                   assumed to be the design matrix Z, allowing the
%                   specification individualized models. 
%--------------------------------------------------------------------------
% OPTION:
%   'bootstrapIter': Number of bootstrap interations (default 1000). 
%   'runEffect': How to deal with effects that may be specific to different
%                imaging runs:
%                  'random': Models variance of the run effect for each subject
%                            as a seperate random effects parameter.
%                  'fixed': Consider run effect a fixed effect, will be removed
%                            implicitly using ReML (default) 
%   'fitScale':    Fit additional scaling parameter for each of the
%                  subjects? Defaults to true. This makes sense, as signal/ noise is not 
%                  the same across subjects. 
%                  However, when you want to test strongly against a null model that
%                  does not predict any difference between different conditions, 
%                  then the scaling parameter makes the alternative model
%                  more flexible  
%   'S',S         : Structure of the NxN noise covariance matrix -
%                   otherwise independence is assumed
%   'fitAlgorithm': Either 'NR' or 'minimize' - provides over-write for
%                   model specific algorithms 
%--------------------------------------------------------------------------
% OUTPUT:
%    theta_boot:  Nparams x boostrapIter array of paramater bootstrap samples. 
%    theta_g:     Nparams vector of group parameters 

runEffect       = 'fixed';
bootstrapIter   = 1000; 
verbose         = 1;    
fitScale        = 1;   % Fit an additional scaling parameter for each subject? 
S               = [];  % Structure of noise matrix 
fitAlgorithm    = [];  % Over-write on model specific fit Algorithm 
theta0          = [];
scalePrior      = 10; 
pcm_vararginoptions(varargin,{'bootstrapIter','runEffect',...
                      'verbose','fitScale','S','fitAlgorithm','theta0','scalePrior'});
numSubj     = numel(Y);

% Determine number of models 
if (iscell(M)) 
    error('Can only handle a single model for boostrap') 
end
M = {M};
numModels = 1; 

% Determine optimal algorithm for each of the models 
if (~isempty(fitAlgorithm)) 
    M{1}.fitAlgorithm = fitAlgorithm; 
end; 
M = pcm_optimalAlgorithm(M); 

[T,theta_g]=pcm_fitModelGroup(Y,M,partitionVec,conditionVec,'runEffect',runEffect,...
    'verbose',verbose,'fitScale',fitScale,'S',S);
theta_boot = zeros(length(theta_g{1}),bootstrapIter); 
for n=1:bootstrapIter
    i = randsample(numSubj,numSubj,true);
    [T,theta_b]=pcm_fitModelGroup(Y(i),M,partitionVec(i),conditionVec(i),'runEffect',runEffect,...
    'verbose',verbose,'fitScale',fitScale,'S',S);
    theta_boot(:,n)=theta_b{1};
end
theta_g = theta_g{1};


