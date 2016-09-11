function [T,M,Iter,G_hat]=pcm_fitModelIndivid(Y,M,partitionVec,conditionVec,varargin);
% function [T,M,Iter,Ghat]=pcm_fitModelIndivid(Y,M,partitionVec,conditionVec,varargin);
% Fits pattern component model(s) specified by M to data from a number of
% subjects .
% The model parameters are all individually fit.
%==========================================================================
% INPUT:
%        Y: {#Subjects}.[#Conditions x #Voxels]
%            Observed/estimated beta regressors from each subject.
%            Preferably multivariate noise-normalized beta regressors.
%
%        M: (#Models).subfields: Structure that defines model(s). Can
%            define multiple models if desired. Contains subfields:
%              .type:        Type of the model to be fitted
%                             'fixed': Fixed structure without parameter (except sale for each subject)
%                             'component': G is a sum of linear components
%                             'squareroot': G=A*A', with A a linear sum of weighted components
%                             'nonlinear': Nonlinear model with own function to return derivatives
%              .numGparams:  Scalar that defines the number of parameters
%                             included in model.
%              .theta0:      Vector of starting theta parameters to calculate predicted
%                             model G. Can estimate
%                             these parameters using 'pcm_modelpred_free_startingval'
%              .modelpred':  Modelling func. Must take theta values as vector
%                             and return predicated second moment matrix and
%                             derivatives in respect to parameters (for nonlinear models).
%              .Gc:          Linear component matrices (for type 'component')
%              .Ac           Linear component matrices (for type 'squareroot')
%
%   conditionVec: {#Subjects} Cell array with condition assignment vector
%                   for each subject. Rows of conditionVec{subj} define
%                   condition assignment of rows of Y{subj}.
%
%   partitionVec: {#Subjects} Cell array with partition assignment vector
%                   for each subject. Rows of partitionVec{subj} define
%                   partition assignment of rows of Y{subj}.
%                   Commonly these are the scanning run #s for beta
%                   regressors.
%--------------------------------------------------------------------------
% OPTION:
%   'runEffect': How to deal with effects that may be specific to different
%                imaging runs:
%                  'random': Models variance of the run effect for each subject
%                            as a seperate random effects parameter.
%                  'fixed': Consider run effect a fixed effect, will be removed
%                            implicitly using ReML.
%                  'remove': Forced removal of the run effect before
%                            random effects modelling - Simply adjusts the
%                            error covariance matrix to reflect he removal
%
%   'isCheckDeriv: Check the derivative accuracy of theta params. Done using
%                  'checkderiv'. This function compares input to finite
%                  differences approximations. See function documentation.
%
%   'MaxIteration': Number of max minimization iterations. Default is 1000.
%
%    'isCheckIter': Optional flag to display the summary of function
%                   iteration at the end of fitting. Default is 0.
%
%    'isCheckTime': Optional flag to display the time took for each model.
%                   Default is 1.
%--------------------------------------------------------------------------
% OUTPUT:
%   T:      Structure with following subfields:
%       SN:                 Subject number
%       likelihood:         Crossvalidated likelihood
%       likelihood_all:     Group-fit likelihood (with subj included)
%       scale0:             Starting value for scale
%       noise0:             Starting value for noise
%       run0:               Starting value for run
%       scale_all:
%       noise_all:
%       run_all:
%       iterations:         Number of interations for model fit
%       time:               Elapsed time in sec 
%
%   M(m):    Structure array of models - with appended fields
%       theta_all:  Estimated parameters at the overall fitting.
%       G_pred:     Predicted second moment matrix of model(s)
%       theta:      Estimated parameters (model + scaling/noise parameters)
%                   across the crossvalidation rounds .

runEffect       = 'random';
isCheckDeriv    = 0;
MaxIteration    = 1000;
Iter            = [];
verbose         = 1; 
vararginoptions(varargin,{'runEffect','isCheckDeriv','MaxIteration','verbose'});

numSubj     = numel(Y);
numModels   = numel(M);

% Preallocate output structure
T.SN = [1:numSubj]';

% Now loop over subject and provide inidivdual fits 
for s = 1:numSubj
    
    % Prepare matrices and data depnding on how to deal with run effect 
    [N(s,1),P(s,1)] = size(Y{s});
    Z{s}   = indicatorMatrix('identity_p',conditionVec{s});
    numCond= size(Z{s},2);
    switch (runEffect)
        case 'random'
            YY{s}  = (Y{s} * Y{s}');
            B{s}   = indicatorMatrix('identity_p',partitionVec{s});
            X{s}   = [];
            S      = [];   % Use indentity for covariance
        case 'fixed'
            YY{s}  = (Y{s} * Y{s}');
            B{s}  =  [];
            X{s}  =  indicatorMatrix('identity_p',partitionVec{s});
            S     =  [];    % Use indentity for covariance
        case 'remove'  % This currently doesn't work as intended - dimensionality reduction needed 
            Run    =  indicatorMatrix('identity_p',partitionVec{s});
            R      =  eye(N(s))-Run*((Run'*Run)\Run');
            YY{s}  = (R*Y{s} * Y{s}'*R');
            S(s).S = R*R'
            S(s).invS = pinv(S(s).S);
            B{s}  = [];
    end;
    
    % Estimate starting value run and the noise from a crossvalidated estimate of the second moment matrix 
    [G_hat(:,:,s),Sig_hat(:,:,s)] = pcm_estGCrossval(Y{s},partitionVec{s},conditionVec{s});
    sh = Sig_hat(:,:,s);
    run0(s)     = log((sum(sum(sh))-trace(sh))/(numCond*(numCond-1)));
    noise0(s)   = log(trace(sh)/numCond-exp(run0(s)));
    
    % Print where we are in the fitting process 
    if (verbose)
        fprintf('Subj %d\n',s);
    end; 
    
    % Now loop over models 
    for m = 1:length(M)
        if (verbose) 
            if (isfield(M,'name')) && (~isempty(M(m).name));
                fprintf('fitting model:%s\n',M(m).name);
            else
                fprintf('fitting model:%d\n',m);
            end;
            tic;
        end; 
        
        % Get starting guess for theta if not provided
        if (isfield(M(m),'theta0'))
            theta0 = M(m).theta0;
        else
            theta0 = pcm_getStartingval(M(m),mean(G_hat,3));   
        end;

        % Use normal linear regression to get scaling parameter for the
        % subject
        switch (M(m).type)
            case 'noiseceiling'
                G0 = mean(G_hat,3); 
                M(m).numGparams=0; 
                M(m).Gc = G0;
            otherwise 
                G0 = pcm_calculateG(M(m),theta0);
        end; 
        g0 = vec(G0);

        % Estimate starting scaling value for each subject
        g_hat         = vec(G_hat(:,:,s));
        scaling       = (g0'*g_hat)/(g0'*g0);
        if ((scaling<10e-6)||~isfinite(scaling)); scaling = 10e-6; end;      % Enforce positive scaling
        scale0(s,m)   = log(scaling);
        

        % Now set up the function that returns likelihood and derivative 
        if (isempty(S))
            fcn = @(x) pcm_groupLikelihood(x,{YY{s}},M(m),{Z{s}},{X{s}},P(s),'runEffect',{B{s}});
        else
            fcn = @(x) pcm_groupLikelihood(x,{YY{s}},M(m),{Z{s}},{},P(s),'runEffect',{B{s}},'S',S(s));
        end;
        
        % Set up overall starting values 
        switch (runEffect) 
            case {'fixed','remove'}
                x0  = [theta0;scale0(s,m);noise0(s)];
            case {'random'}
                x0  = [theta0;scale0(s,m);noise0(s);run0(s)];
        end; 
        
        % Use minimize to fine maximum liklhood estimate 
        [theta,fX,i]      =  minimize(x0, fcn, MaxIteration);
        M(m).theta(:,s)   =  theta(1:M(m).numGparams);
        M(m).G_pred       =  pcm_calculateG(M(m),M(m).theta(:,s));
        T.scale(s,m)      =  exp(theta(M(m).numGparams+1)); 
        T.noise(s,m)      =  exp(theta(M(m).numGparams+2)); 
        if strcmp(runEffect,'random')
            T.run(s,m)      =  exp(theta(M(m).numGparams+3)); 
        end; 
        T.likelihood(s,m) =  -fX(end);  %invert the sign 
        T.iterations(s,m) = i;
        T.time(s,m)       = toc; 
        if verbose
            fprintf('... done!');
            fprintf('\t Iterations %d, Elapsed time: %3.3f\n',T.iterations(s,m),T.time(s,m));
        end; 
    end; % for each model
end; % for each subject
