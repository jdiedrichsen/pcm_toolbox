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
%        M: {#Models} Cell array with structure that defines model(s). Each may contain the following 
%                       Subfield 
%              .type:        Type of the model to be fitted
%                             'fixed':     Fixed structure without parameter (except scale for each subject)
%                             'component': G is a sum of linear components
%                             'feature':   G=A*A', with A a linear sum of weighted feature components
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
%   partitionVec: {#Subjects} Cell array with partition assignment vector
%                   for each subject. Rows of partitionVec{subj} define
%                   partition assignment of rows of Y{subj}.
%                   Commonly these are the scanning run #s for beta
%                   regressors.
%                   If a single vector is provided, it is assumed to me the
%                   same for all subjects 
%
%   conditionVec: {#Subjects} Cell array with condition assignment vector
%                   for each subject. Rows of conditionVec{subj} define
%                   condition assignment of rows of Y{subj}.
%                   If a single vector is provided, it is assumed to me the
%                   same for all subjects 
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
%       likelihood:         likelihood
%       noise:              Noise parameter 
%       run:                Run parameter (if run = 'random') 
%       iterations:         Number of interations for model fit
%       time:               Elapsed time in sec 
%
%   M{m}:    Structure array of models - with appended fields
%       G_pred:     Predicted second moment matrix of model: 3-D array  
%                   with 1 slice per subject 
%       theta:      Estimated parameters (model + scaling/noise parameters)
%                   for individual subjects

runEffect       = 'random';
isCheckDeriv    = 0;
MaxIteration    = 1000;
Iter            = [];
verbose         = 1; 
pcm_vararginoptions(varargin,{'runEffect','isCheckDeriv','MaxIteration','verbose'});

numSubj     = numel(Y);
numModels   = numel(M);

% Preallocate output structure
T.SN = [1:numSubj]';

% Now loop over subject and provide inidivdual fits 
for s = 1:numSubj
    
    % If condition and partition Vectors are not cells, assume they are the
    % same 
    if (iscell(conditionVec)) 
        cV = conditionVec{s}; 
        pV = partitionVec{s}; 
    else 
        cV = conditionVec; 
        pV = partitionVec; 
    end; 
    
    % Prepare matrices and data depnding on how to deal with run effect 
    [N(s,1),P(s,1)] = size(Y{s});
    Z{s}   = pcm_indicatorMatrix('identity_p',cV);
    numCond= size(Z{s},2);
    switch (runEffect)
        case 'random'
            YY{s}  = (Y{s} * Y{s}');
            B{s}   = pcm_indicatorMatrix('identity_p',pV);
            X{s}   = [];
            S      = [];   % Use indentity for covariance
        case 'fixed'
            YY{s}  = (Y{s} * Y{s}');
            B{s}  =  [];
            X{s}  =  pcm_indicatorMatrix('identity_p',pV);
            S     =  [];    % Use indentity for covariance
        case 'remove'  % This currently doesn't work as intended - dimensionality reduction needed 
            Run    =  indicatorMatrix('identity_p',pV);
            R      =  eye(N(s))-Run*((Run'*Run)\Run');
            YY{s}  = (R*Y{s} * Y{s}'*R');
            S(s).S = R*R'
            S(s).invS = pinv(S(s).S);
            B{s}  = [];
    end;
    
    % Estimate starting value run and the noise from a crossvalidated estimate of the second moment matrix 
    [G_hat(:,:,s),Sig_hat(:,:,s)] = pcm_estGCrossval(Y{s},pV,cV);
    sh = Sig_hat(:,:,s);
    run0(s)     = real(log((sum(sum(sh))-trace(sh))/(numCond*(numCond-1))));
    noise0(s)   = real(log(trace(sh)/numCond-exp(run0(s))));
    
    % Now loop over models 
    for m = 1:length(M)
        if (verbose) 
            if isfield(M,'name');
                fprintf('Fitting Subj: %d model:%s\n',s,M{m}.name);
            else
                fprintf('Fitting Subj: %d model:%ds\n',s,m);
            end;
        end; 
        tic; 
        
        % Get starting guess for theta if not provided
        if (isfield(M{m},'theta0'))
            theta0 = M{m}.theta0;
        else
            theta0 = pcm_getStartingval(M{m},G_hat(:,:,s));   
        end;
        
        % Now set up the function that returns likelihood and derivative 
        if (isempty(S))
            fcn = @(x) pcm_likelihoodIndivid(x,YY{s},M{m},Z{s},X{s},P(s),'runEffect',B{s});
        else
            fcn = @(x) pcm_likelihoodIndivid(x,YY{s},M{m},Z{s},[],P(s),'runEffect',B{s},'S',S(s));
        end;
        
        % Set up overall starting values 
        switch (runEffect) 
            case {'fixed','remove'}
                x0  = [theta0;noise0(s)];
            case {'random'}
                x0  = [theta0;noise0(s);run0(s)];
        end; 
        
        % Use minimize to fine maximum liklhood estimate 
        [theta,fX,i]      =  minimize(x0, fcn, MaxIteration);
        M{m}.thetaIndiv(:,s)   =  theta(1:M{m}.numGparams);
        M{m}.G_pred       =  pcm_calculateG(M{m},M{m}.thetaIndiv(:,s));
        T.noise(s,m)      =  exp(theta(M{m}.numGparams+1)); 
        if strcmp(runEffect,'random')
            T.run(s,m)      =  exp(theta(M{m}.numGparams+2)); 
        end; 
        T.likelihood(s,m) =  -fX(end);  %invert the sign 
        T.iterations(s,m) = i;
        T.time(s,m)       = toc; 
    end; % for each model
end; % for each subject
