function [T,theta_hat,G_pred,INFO]=pcm_fitModelIndivid(Y,M,partitionVec,conditionVec,varargin);
% function [T,theta_hat,G_pred]=pcm_fitModelIndivid(Y,M,partitionVec,conditionVec,varargin);
% Fits pattern component model(s) specified by M to data from a number of
% subjects.
% The model parameters are all individually fit.
%==========================================================================
% INPUT:
%        Y: {#Subjects}.[#Conditions x #Voxels]
%            Observed/estimated beta regressors from each subject.
%            Preferably multivariate noise-normalized beta regressors.
%
%        M: {#Models} Cell array with structure that defines model(s). Each
%        may contain the fields
%              .type:        Type of the model to be fitted
%                             'fixed':     Fixed structure without parameter (except scale for each subject)
%                             'component': G is a sum of linear components
%                             'feature':   G=A*A', with A a linear sum of weighted feature components
%                             'nonlinear': Nonlinear model with own function to return derivatives
%                             'freechol':  Free model in Cholesky form
%              .numGparams:  Scalar that defines the number of parameters
%                             included in model.
%              .theta0:      Vector of starting theta parameters to calculate predicted
%                             model G.
%         for more fields see the manual for model specification.
%
%   partitionVec: {#Subjects} Cell array with partition assignment vector
%                   for each subject. Rows of partitionVec{subj} define
%                   partition assignment of rows of Y{subj}.
%                   Commonly these are the scanning run #s for beta
%                   regressors.
%                   If a single vector is provided, it is assumed to be the
%                   same for all subjects
%
%   conditionVec: {#Subjects} Cell array with condition assignment vector
%                   for each subject. Rows of conditionVec{subj} define
%                   condition assignment of rows of Y{subj}.
%                   If a single vector is provided, it is assumed to me the
%                   same for all subjects.
%                   If the (elements of) conditionVec are matrices, it is
%                   assumed to be the design matrix Z, allowing the
%                   specification individualized models.
%
%--------------------------------------------------------------------------
% OPTION:
%   'runEffect': How to deal with effects that may be specific to different
%                imaging runs:
%                  'random': Models variance of the run effect for each subject
%                            as a seperate random effects parameter.
%                  'fixed': Consider run effect a fixed effect, will be removed
%                            implicitly using ReML.
%                  'none': No modeling of the run Effect (not recommended for fMRI data) 
%
%   'isCheckDeriv: Check the derivative accuracy of theta params. Done using
%                  'checkderiv'. This function compares input to finite
%                  differences approximations. See function documentation.
%
%   'fitAlgorithm': Either 'NR' or 'minimize' - provides over-write for
%                   model specific algorithms
%   'maxIteration': Number of max minimization iterations. Default is 1000.
%
%   'likeThres':    Threshold for minimal change in the log-likelihood to
%                   stop iterations 
%   'verbose':      Optional flag to show display message in the command
%                   line. Default is 1. Setting to 2 gives more detailed
%                   feedback on pcm_NR
%
%   'S':            Optional specific covariance structure of the noise
%                   {#Subjects} = cell array of N_s x N_s matrices 
%                   or 
%                   S(#Subjects).S and .invS: Structure of the N_s x N_s 
%                   normal and inverse covariances matrices
%
%   'theta0':       Cell array of starting values (same format as theta{m})
%   'fitScale':     Fit a additional scale parameter for each subject?
%                   (0/1). Default is set to 0.
%--------------------------------------------------------------------------
% OUTPUT:
%   T:      Structure with following subfields:
%       SN:                 Subject number
%       likelihood:         log-likelihood
%       scale:              Scale parameter (if fitscale = 1)-exp(theta_s)
%       noise:              Noise parameter- exp(theta_eps)
%       run:                Run parameter (if run = 'random')
%       iterations:         Number of interations for model fit
%       time:               Elapsed time in sec
%
%   theta{m}     Cell array of estimated model parameters, each a
%                 #params x #numSubj matrix
%   G_pred{m}     Cell array of estimated G-matrices under the model

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

numSubj     = numel(Y);

% Set options for fitting algorithm 
fitOPT.verbose = verbose; 
if (~isempty(maxIteration))
    fitOPT.maxIteration = maxIteration;
end 
if (~isempty(likeThres))
    fitOPT.likeThres = likeThres; 
end 

% Determine number of models
if (~iscell(M))
    M={M};
end
numModels = numel(M);

% Preallocate output structure
T.SN = [1:numSubj]';

% Determine optimal algorithm for each of the models
if (~isempty(fitAlgorithm))
    for m=1:numModels
        M{m}.fitAlgorithm = fitAlgorithm;
    end
end
M = pcm_optimalAlgorithm(M);

% Set up all parameters for the upcoming fit
[Z,B,X,YY,S,N,P,G_hat,noise0,run0]=...
    pcm_setUpFit(Y,partitionVec,conditionVec,'runEffect',runEffect,'S',S);

% Loop over subject and provide inidivdual fits
for s = 1:numSubj
    
    % Now loop over models
    for m = 1:length(M)
        if (verbose)
            if isfield(M{m},'name');
                fprintf('Fitting Subj: %d model:%s\n',s,M{m}.name);
            else
                fprintf('Fitting Subj: %d model:%d\n',s,m);
            end
        end
        tic;
        
        % Get starting guess for theta0 is not provided
        if (~isfield(M{m},'theta0'))
            M{m}.theta0=pcm_getStartingval(M{m},G_hat(:,:,s));
        end

        % If naive noise ceiling model, use G_hat
        if strcmp(M{m}.type,'freedirect')
            G0 = pcm_makePD(mean(G_hat,3)); 
            M{m}.numGparams=0; 
            M{m}.Gc = G0;
        else 
            G0 = pcm_calculateG(M{m},M{m}.theta0);
        end
        g0 = G0(:);
        
        % If fitting scale parameterm get initial estimate using OLS 
        if (fitScale)
            g_hat         = G_hat(:,:,s);
            g_hat         = g_hat(:);
            scaling       = (g0'*g_hat)/(g0'*g0);
            if ((scaling<10e-5)||~isfinite(scaling)); scaling = 10e-5; end      % Enforce positive scaling
            scale0(s,m)   = log(scaling);
        end
                
        if (numel(theta0)<m || size(theta0{m},2)<s)
            th0m = [M{m}.theta0;noise0(s)];
            if(fitScale)
                th0m = [th0m;scale0(s,m)];
            end
            if (strcmp(runEffect,'random'))
                th0m = [th0m;run0(s)];
            end
            theta0{m}(:,s)=th0m;
        end; 
        
        % Set options for likelihood function 
        OPT.fitScale = fitScale;
        OPT.runEffect=B{s};
        if (~isempty(S))
            OPT.S = S(s);
        end;
                
        
        % Now do the fitting, using the preferred optimization routine
        switch (M{m}.fitAlgorithm)
            case 'minimize'  % Use minimize to find maximum liklhood estimate runEffect',B{s});
                fcn = @(x) pcm_likelihoodIndivid(x,YY{s},M{m},Z{s},X{s},P(s),OPT);
                if (isempty(maxIteration)) 
                    maxIteration=1000;
                end 
                [theta_hat{m}(:,s),fX,T.iterations(s,m)]      =  ...
                    minimize(theta0{m}(:,s), fcn, maxIteration);
                T.likelihood(s,m) =  -fX(end);  %invert the sign
            case 'NR'
                fcn = @(x) pcm_likelihoodIndivid(x,YY{s},M{m},Z{s},X{s},P(s),OPT);
                [theta_hat{m}(:,s),T.likelihood(s,m),T.iterations(s,m),T.reg(s,m),INFO.regH{s,m},INFO.thetaH{s,m}]= ... 
                    pcm_NR(theta0{m}(:,s),fcn,fitOPT); 
        end
        
        G_pred{m}(:,:,s)  =  pcm_calculateG(M{m},theta_hat{m}(1:M{m}.numGparams,s));
        T.noise(s,m)      =  exp(theta_hat{m}(M{m}.numGparams+1,s));
        if (fitScale)
            T.scale(s,m)      = exp(theta_hat{m}(M{m}.numGparams+2,s));
        end
        
        if strcmp(runEffect,'random')
            T.run(s,m)      =  exp(theta_hat{m}(M{m}.numGparams+2+fitScale,s));
        end
        T.time(s,m)       = toc;
        
        % This is an optional check if the dervivate calculation is correct
        if (isCheckDeriv)
            d = pcm_checkderiv(fcn,theta_hat{m}(:,s)-0.001,0.00001);
            fprintf('discrepency :%d\n',d);
            keyboard; 
        end
    end % for each model
end % for each subject
INFO.theta0=theta0;