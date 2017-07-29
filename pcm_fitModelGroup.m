function [T,theta_hat,G_pred]=pcm_fitModelGroup(Y,M,partitionVec,conditionVec,varargin);
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
%                            implicitly using ReML.
%                  'remove': Forced removal of the run effect before
%                            random effects modelling - Simply adjusts the
%                            error covariance matrix to reflect he removal
%   'fitScale':    Fit additional scaling parameter for each of the
%                  subjects? Defaults to 1. This makes a lot of sense, if the scaling of the 
%                  data is not of the same intensity across subjects. 
%                  However, when you want to test strongly against a null model that
%                  does not predict any difference between different conditions, 
%                  then the scaling parameter makes the alternative model
%                  more flexible and it will more often. 
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
%--------------------------------------------------------------------------
% OUTPUT:
%   T:      Structure with following subfields:
%       SN:                 Subject number
%       likelihood:     Group-fit likelihood (with subj included)
%       noise:          Noise Variance 
%       scale:          Scale parameter for each subject (if fitScale set to 1)  
%       run:            Run variance (if runEffect = 'random'); 
%
%    theta_hat:  Estimated parameters at the overall fitting (including
%                noise and scale parameters). A mx1 cell array   
%    Gpred:      Predicted second moment matrix for the model from group
%                fit for each model. A mx1 cell array 

runEffect       = 'random';
isCheckDeriv    = 0;
MaxIteration    = 1000;
verbose         = 1;    
fitScale        = 1;   % Fit an additional scaling parameter for each subject? 
S               = [];  % Structure of noise matrix 
pcm_vararginoptions(varargin,{'runEffect','isCheckDeriv','MaxIteration',...
                      'verbose','fitScale','S'});

numSubj     = numel(Y);
numModels   = numel(M);

% Preallocate output structure
T.SN = [1:numSubj]';

% --------------------------------------------------------
% Figure out a starting values for the noise parameters
% --------------------------------------------------------
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
    
    % Check if conditionVec is condition or design matrix
    if size(cV,2)==1;
        Z{s}   = pcm_indicatorMatrix('identity_p',cV);
    else
        Z{s} = cV;
    end;
    
    % Set up the main matrices
    [N(s,1),P(s,1)] = size(Y{s});   
    numCond= size(Z{s},2);
    YY{s}  = (Y{s} * Y{s}');
    
    % Depending on the way of dealing with the run effect, set up data
    switch (runEffect)
        case 'random'
            B{s}   = pcm_indicatorMatrix('identity_p',pV);
            X{s}   = [];
        case 'fixed'
            B{s}  =  [];
            X{s}  =  pcm_indicatorMatrix('identity_p',pV);
    end;
    
    % Estimate crossvalidated second moment matrix to get noise and run
    [G_hat(:,:,s),Sig_hat(:,:,s)] = pcm_estGCrossval(Y{s},pV,cV);
    sh              = Sig_hat(:,:,s);
    run0(s,1)       = real(log((sum(sum(sh))-trace(sh))/(numCond*(numCond-1))));
    noise0(s,1)     = real(log(trace(sh)/numCond-exp(run0(s))));
end;

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
    
    % Get starting guess for theta if not provided
    if (isfield(M{m},'theta0'))
        theta0 = M{m}.theta0;
    else
        theta0 = pcm_getStartingval(M{m},mean(G_hat,3));   
    end;
    
    % Use normal linear regression to get scaling parameter for the
    % subject
    switch (M{m}.type)
        case 'noiseceiling'
            G0 = mean(G_hat,3); 
            M{m}.numGparams=0; 
            M{m}.Gc = G0;
        otherwise 
            G0 = pcm_calculateG(M{m},theta0);
    end; 
    g0 = G0(:);
    
    % Estimate starting scaling value for each subject
    if (fitScale) 
        for s = 1:numSubj
            g_hat         = G_hat(:,:,s);
            g_hat         = g_hat(:); 
            scaling       = (g0'*g_hat)/(g0'*g0);
            if ((scaling<10e-6)||~isfinite(scaling)); scaling = 10e-6; end;      % Enforce positive scaling
            scale0(s,m)   = log(scaling);
        end;
    end; 
    
    % Put together the vector of starting value 
    x0 = theta0; 
    x0=[x0;noise0]; 
    if(fitScale) 
        x0 = [x0;scale0(:,m)]; 
    end; 
    if (strcmp(runEffect,'random'))
        x0 = [x0;run0];
    end; 
    
    % Now do the fitting
    if (M{m}.numGparams==0) 
        fcn = @(x) pcm_likelihoodGroup(x,YY,G0,Z,X,P,'runEffect',B,'S',S,'fitScale',fitScale);
    else 
        fcn = @(x) pcm_likelihoodGroup(x,YY,M{m},Z,X,P,'runEffect',B,'S',S,'fitScale',fitScale);
    end;
    [theta_hat{m},fx]       = minimize(x0, fcn, MaxIteration);
    
    T.noise(:,m)      = exp(theta_hat{m}(M{m}.numGparams+1:M{m}.numGparams+numSubj));
    if (fitScale) 
        T.scale(:,m)      = exp(theta_hat{m}(M{m}.numGparams+numSubj+1:M{m}.numGparams+2*numSubj));
    end; 
    if (strcmp(runEffect,'random'))
        T.run(:,m)  = exp(theta_hat{m}(M{m}.numGparams+(1+fitScale)*numSubj+1:M{m}.numGparams+(2+fitScale)*numSubj));
    end;
    
    G_pred{m}        = pcm_calculateG(M{m},theta_hat{m}(1:M{m}.numGparams));
    
    % Compute log-likelihood under the estimated parameters
    [~,~,T.likelihood(:,m)] = fcn(theta_hat{m});
    
    % This is an optional check if the dervivate calculation is correct
    if (isCheckDeriv)
        d = pcm_checkderiv(fcn,thetaAll-0.01,0.0000001);
    end;
end;
