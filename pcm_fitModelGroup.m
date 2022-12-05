function [T,theta_hat,G_pred,theta0]=pcm_fitModelGroup(Y,M,partitionVec,conditionVec,varargin);
% function [T,theta_hat,G_pred]=pcm_fitModelCrossval(Y,M,partitionVec,conditionVec,varargin);
% Fits pattern component model(s) specified by M to data from a number of
% subjects.
% The model parameters are shared - the noise parameters are not.
% If provided with a G_hat and Sig_hat, it also automatically estimates close
% startingvalues.
%==========================================================================
% INPUT:
%        Y: {#Subjects}.[#Conditions x #Voxels]
%            Observed/estimated beta regressors from each subject.
%            Preferably multivariate noise-normalized beta regressors.
%
%        M: cell arrays {#Models} of model(s) to be fitted 
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
%   'isCheckDeriv: Check the derivative accuracy of theta params. Done using
%                  'pcm_checkderiv'. This function compares input to finite
%                  differences approximations. See function documentation.
%
%   'MaxIteration': Number of max minimization iterations. Default is 1000.
%
%   'verbose':      Optional flag to show display message in the command
%                   line (e.g., elapsed time). Default is 1.
% 
%   'S',S         : Structure of the NxN noise covariance matrix -
%                   otherwise independence is assumed
%   
%   'fitAlgorithm': Either 'NR' or 'minimize' - provides over-write for
%                   model specific algorithms 
%--------------------------------------------------------------------------
% OUTPUT:
%   T:      Structure with following subfields:
%       SN:             Subject number
%       likelihood:     Group-fit likelihood (with subj included)
%       noise:          Noise Variance 
%       scale:          Scale parameter for each subject (if fitScale set to 1)  
%       run:            Run variance (if runEffect = 'random'); 
%
%    theta_hat:  Estimated parameters at the overall fitting (including
%                noise, scale, and run parameters). A mx1 cell array   
%    Gpred:      Predicted second moment matrix for the model from group
%                fit for each model. A mx1 cell array 

runEffect       = 'fixed';
isCheckDeriv    = 0;
MaxIteration    = 1000;
verbose         = 1;    
fitScale        = 1;   % Fit an additional scaling parameter for each subject? 
S               = [];  % Structure of noise matrix 
fitAlgorithm    = [];  % Over-write on model specific fit Algorithm 
theta0          = [];
scalePrior      = 10; 
pcm_vararginoptions(varargin,{'runEffect','isCheckDeriv','MaxIteration',...
                      'verbose','fitScale','S','fitAlgorithm','theta0','scalePrior'});
numSubj     = numel(Y);

% Determine number of models 
if (~iscell(M)) 
    M={M}; 
end; 
numModels = numel(M); 


% Preallocate output structure
T.SN = [1:numSubj]';
T.iterations = zeros(numSubj,1); 
T.time = zeros(numSubj,1); 

% Determine optimal algorithm for each of the models 
if (~isempty(fitAlgorithm)) 
    for m=1:numModels
        M{m}.fitAlgorithm = fitAlgorithm; 
    end; 
end; 
M = pcm_optimalAlgorithm(M); 

% Set up all parameters for the upcoming fit 
[Z,B,X,YY,S,N,P,G_hat,noise0,run0]=...
    pcm_setUpFit(Y,partitionVec,conditionVec,'runEffect',runEffect,'S',S);

% -----------------------------------------------------
% Perform the overall group fit across all subjects
% -----------------------------------------------------
for m = 1:numModels
    if (verbose)
        if isfield(M{m},'name');
            fprintf('Overall fitting model:%s\n',M{m}.name);
        else
            fprintf('Overall fitting model:%d\n',m);
        end;
    end;
    tic; 
    
    % Get starting guess for theta if not provided
    if (~isfield(M{m},'theta0'))
        M{m}.theta0 = pcm_getStartingval(M{m},mean(G_hat,3));   
    end; 
    
    % Use normal linear regression to get scaling parameter for the
    % subject
    switch (M{m}.type)
        case 'freedirect'
            G0 = pcm_makePD(mean(G_hat,3)); 
            M{m}.numGparams=0; 
            M{m}.Gc = G0;
        otherwise 
            G0 = pcm_calculateG(M{m},M{m}.theta0);
    end; 
    g0 = G0(:);
    if (fitScale) 
        for s = 1:numSubj
            g_hat         = G_hat(:,:,s);
            g_hat         = g_hat(:); 
            scaling       = (g0'*g_hat)/(g0'*g0);
            if ((scaling<10e-5)||~isfinite(scaling)); scaling = 10e-5; end;      % Enforce positive scaling
            scale0(s,m)   = log(scaling);
        end;
    end; 
    
    % Put together the vector of starting value 
    if (numel(theta0)<m) 
        theta0{m} = [M{m}.theta0;noise0]; 
        if(fitScale) 
            theta0{m} = [theta0{m};scale0(:,m)]; 
        end; 
        if (strcmp(runEffect,'random'))
            theta0{m} = [theta0{m};run0];
        end; 
    end; 
    
    % Set options 
    OPT.runEffect=B; 
    OPT.S = S; 
    OPT.fitScale = fitScale; 
    OPT.scalePrior = scalePrior; 
    
    % Now do the fitting
    switch (M{m}.fitAlgorithm)
        case 'minimize'  % Use minimize to find maximum liklhood estimate runEffect',B{s});
            fcn = @(x) pcm_likelihoodGroup(x,YY,M{m},Z,X,P,OPT);
            [theta_hat{m},~,T.iterations(:,m)] = minimize(theta0{m}, fcn, MaxIteration);
        case 'NR' 
            fcn = @(x) pcm_likelihoodGroup(x,YY,M{m},Z,X,P,OPT);
            [theta_hat{m},~,T.iterations(:,m),T.reg(:,m)] = pcm_NR(theta0{m}, fcn);
        otherwise 
            error('unknown fitting Algorith: %s',M{m}.fitAlgorithm);
    end; 

    % retrieve parameters 
    T.noise(:,m)      = exp(theta_hat{m}(M{m}.numGparams+1:M{m}.numGparams+numSubj));
    if (fitScale) 
        T.scale(:,m)      = exp(theta_hat{m}(M{m}.numGparams+numSubj+1:M{m}.numGparams+2*numSubj));
    end; 
    if (strcmp(runEffect,'random'))
        T.run(:,m)  = exp(theta_hat{m}(M{m}.numGparams+(1+fitScale)*numSubj+1:M{m}.numGparams+(2+fitScale)*numSubj));
    end;
    G_pred{m}        = pcm_calculateG(M{m},theta_hat{m}(1:M{m}.numGparams));
    
    % Compute log-likelihood under the estimated parameters
    [~,~,~,T.likelihood(:,m)] = fcn(theta_hat{m});
    
    T.time(:,m) = toc; 
    % This is an optional check if the dervivate calculation is correct
    if (isCheckDeriv)
        d = pcm_checkderiv(fcn,theta_hat{m}-0.01,0.0000001);
    end;
end;
