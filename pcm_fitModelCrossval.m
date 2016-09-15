function [T,M]=pcm_fitModelCrossval(Y,M,partitionVec,conditionVec,varargin);
% function [T,theta_all,G_pred,theta]=pcm_fitModelCrossval(Y,M,partitionVec,conditionVec,varargin);
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
%        M: (#Models).subfields: Structure that defines model(s). Can
%            define multiple models if desired. Contains subfields:
%              .type:        Type of the model to be fitted
%                             'fixed': Fixed structure without parameter
%                                (except scale and noise for each subject)
%                             'component': G is a sum of linear components,
%                                specified by Gc
%                             'squareroot': G=A*A', with A a linear sum of
%                                 weighted components Ac
%                             'nonlinear': Nonlinear model with own function
%                                 to return G-matrix and derivatives
%                             'noiseceiling': Uses the mean estimated
%                                 G-matrix from crossvalidation to get an
%                                 estimate of the best achievavle fit
%              .numGparams:  Scalar that defines the number of parameters
%                             included in model.
%              .theta0:      Vector of starting values for theta. If not given,
%                              the function attempts to estimate these from a
%                              crossvalidated version Values usually estimated from
%                             observed second-moment matrix. Can estimate
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
%
%   conditionVec: {#Subjects} Cell array with condition assignment vector
%                   for each subject. Rows of conditionVec{subj} define
%                   condition assignment of rows of Y{subj}.
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
%                  'pcm_checkderiv'. This function compares input to finite
%                  differences approximations. See function documentation.
%
%   'MaxIteration': Number of max minimization iterations. Default is 1000.
%
%   'verbose':      Optional flag to show display message in the command
%                   line (e.g., elapsed time). Default is 1.
% 
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
%
%   M(m):    Structure array of models - with appended fields
%       theta_all:  Estimated parameters at the overall fitting.
%       G_pred:     Predicted second moment matrix of model(s)
%       theta:      Estimated parameters (model + scaling/noise parameters)
%                   across the crossvalidation rounds .

runEffect       = 'random';
isCheckDeriv    = 0;
MaxIteration    = 1000;
verbose         = 1;
vararginoptions(varargin,{'runEffect','isCheckDeriv','MaxIteration','verbose'});

numSubj     = numel(Y);
numModels   = numel(M);

% Preallocate output structure
T.SN = [1:numSubj]';

% --------------------------------------------------------
% Figure out a starting values for the noise parameters
% --------------------------------------------------------
for s = 1:numSubj
    
    % Set up the main matrices
    [N(s,1),P(s,1)] = size(Y{s});
    Z{s}   = indicatorMatrix('identity_p',conditionVec{s});
    numCond= size(Z{s},2);
    
    % Depending on the way of dealing with the run effect, set up data
    switch (runEffect)
        case 'random'
            YY{s}  = (Y{s} * Y{s}');
            B{s}   = indicatorMatrix('identity_p',partitionVec{s});
            X{s}   = [];
            S      = [];   % Use identity for covariance
        case 'fixed'
            YY{s}  = (Y{s} * Y{s}');
            B{s}  =  [];
            X{s}  =  indicatorMatrix('identity_p',partitionVec{s});
            S     =  [];    % Use identity for covariance
        case 'remove'
            Run         =  indicatorMatrix('identity_p',partitionVec{s});
            R           =  eye(N(s))-Run*((Run'*Run)\Run');
            
            % Remove redundant dimention from S and Z
            for c=1:size(Run,2)
                idx = find(Run(:,c)==1);
                idxrem(c,1) = idx(1); % can be other row
            end
            R(idxrem,:) = [];
            Run(idxrem,:) = [];
            % partitionVec{s}(idxrem) = []; % this also affect estG
            
            % Orthogonalize Y and Z
            Z{s} = R*Z{s};
            %Y{s} = R*Y{s}; % changing Y{s} here affects following
            %                 calculation of G                        
            Yrem        = R*Y{s}; % should this be kept for further use?
            
            YY{s}       = (Yrem * Yrem');
            S(s).S      = R*R';             % Rotated noise covariance
            S(s).invS   = pinv(S(s).S);     % Pre-calculated inverse of noise covariance
            B{s}        = [];
            X{s}        = [];
    end;
    
    % Estimate crossvalidated second moment matrix to get noise and run
    [G_hat(:,:,s),Sig_hat(:,:,s)] = pcm_estGCrossval(Y{s},partitionVec{s},conditionVec{s});
    sh              = Sig_hat(:,:,s);
    run0(s,1)       = log((sum(sum(sh))-trace(sh))/(numCond*(numCond-1)));
    noise0(s,1)     = log(trace(sh)/numCond-exp(run0(s)));
end;

% -----------------------------------------------------
% Perform the overall group fit across all subjects
% -----------------------------------------------------
for m = 1:numModels
    if (verbose)
        if isfield(M,'name');
            fprintf('Overall fitting model:%s\n',M(m).name);
        else
            fprintf('Overall fitting model:%d\n',m);
        end;
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
    for s = 1:numSubj
        g_hat         = vec(G_hat(:,:,s));
        scaling       = (g0'*g_hat)/(g0'*g0);
        if ((scaling<10e-6)||~isfinite(scaling)); scaling = 10e-6; end;      % Enforce positive scaling
        scale0(s,m)   = log(scaling);
    end;
    
    % Put together the vector of starting value 
    if (strcmp(runEffect,'random'))
        x0      = [theta0;scale0(:,m);noise0;run0];
    elseif (strcmp(runEffect,'remove'))
        % Force initial parameters to be good
        noise0          = real(noise0);
        run0            = real(run0);
        idxr            = (scale0(:,m)<eps)|(~isfinite(scale0(:,m)));
        scale0(idxr,m)  = repmat(log(1),sum(idxr),1);
        scale0(:,m)     = real(scale0(:,m));
        x0              = [theta0;scale0(:,m);noise0];
    else
        x0      = [theta0;scale0(:,m);noise0];
    end;

    % Now do the fitting
    if (M(m).numGparams==0) 
        fcn = @(x) pcm_groupLikelihood(x,YY,G0,Z,X,P,'runEffect',B,'S',S);
    else 
        fcn = @(x) pcm_groupLikelihood(x,YY,M(m),Z,X,P,'runEffect',B,'S',S);
    end;
    theta0              = [theta0;scale0(:,m);noise0;run0];
    [thetaAll,fx]       = minimize(x0, fcn, MaxIteration);
    
    M(m).theta_all      = thetaAll(1:M(m).numGparams);      % Main model parameters
    T.scale_all(:,m)    = (thetaAll(M(m).numGparams+1:M(m).numGparams+numSubj));%exp(thetaAll(M(m).numGparams+1:M(m).numGparams+numSubj));
    T.noise_all(:,m)    = (thetaAll(M(m).numGparams+numSubj+1:M(m).numGparams+2*numSubj));%exp(thetaAll(M(m).numGparams+numSubj+1:M(m).numGparams+2*numSubj));
    if (strcmp(runEffect,'random'))
        T.run_all(:,m)  = (thetaAll(M(m).numGparams+2*numSubj+1:M(m).numGparams+3*numSubj));%exp(thetaAll(M(m).numGparams+2*numSubj+1:M(m).numGparams+3*numSubj));
    end;
    M(m).G_pred         = pcm_calculateG(M(m),M(m).theta_all);
    
    % Compute log-likelihood under the estimated parameters
    [~,~,T.likelihood_all(:,m)] = fcn(thetaAll);
    
    % This is an optional check if the dervivate calculation is correct
    if (isCheckDeriv)
        d = pcm_checkderiv(fcn,thetaAll-0.01,0.0000001);
    end;
end;

% -----------------------------------------------------
% Loop over subjects and obtain crossvalidated likelihoods
% -----------------------------------------------------
for s = 1:numSubj
    % Inform user
    if (verbose)
        fprintf('Subj %d\n',s);
    end;
    
    % determine training set
    notS    = [1:numSubj];
    notS(s) = [];
    
    % Now loop over models
    for m = 1:numModels
        if (verbose)
            if isfield(M,'name');
                fprintf('fitting model:%s\n',M(m).name);
            else
                fprintf('fitting model:%ds\n',m);
            end;
        end;
        tic;
        
        % Now fit to all the subject but the left-out one
        switch (M(m).type)
            case 'fixed'
                G = M(m).Gc;
                i = 0; 
            case 'noiseceiling'
                G = mean(G_hat(:,:,notS),3);    % uses the mean of all other subjects
                G = pcm_makePD(G);
                i = 0; 
            otherwise
                if (isempty(S))
                    fcn = @(x) pcm_groupLikelihood(x,{YY{notS}},M(m),{Z{notS}},{X{notS}},P(notS(:)),'runEffect',{B{notS}});
                else
                    fcn = @(x) pcm_groupLikelihood(x,{YY{notS}},M(m),{Z{notS}},{X{notS}},P(notS(:)),'runEffect',{B{notS}},'S',S(notS));
                end;
                if (strcmp(runEffect,'random'))
                    x0  = [M(m).theta_all;T.scale_all(notS,m);T.noise_all(notS,m);T.run_all(notS,m)];       % Start with G-params from group fit
                else
                    x0  = [M(m).theta_all;T.scale_all(notS,m);T.noise_all(notS,m)];       % Start with theta from group fit
                end;
                [theta,fX,i] =  minimize(x0, fcn, MaxIteration);
                M(m).theta(:,s)=theta(1:M(m).numGparams);
                G   = pcm_calculateG(M(m),M(m).theta(1:M(m).numGparams,s));
        end;
        
        % Now get the fit the left-out subject
        % (maximizing scale and noise coefficients)
        if (strcmp(runEffect,'random'))
            x0      = [T.scale_all(s,m);T.noise_all(s,m);T.run_all(s,m)];
        else
            x0      = [T.scale_all(s,m);T.noise_all(s,m)];
        end;
        if (isempty(S))
            fcn     = @(x) pcm_groupLikelihood(x,{YY{s}},G,{Z{s}},{X{s}},P(s),'runEffect',{B{s}});   % Minize scaling params only
        else
            fcn     = @(x) pcm_groupLikelihood(x,{YY{s}},G,{Z{s}},{X{s}},P(s),'runEffect',{B{s}},'S',S(s));   % Minize scaling params only
        end;
        
        [th{m}(:,s),fX] =  minimize(x0, fcn, MaxIteration);
        
        T.likelihood(s,m) = -fX(end);
        T.iterations(s,m) = i;
        T.time(s,m)       = toc;
        if verbose
            fprintf('... done!');
            fprintf('\t Iterations %d, Elapsed time: %3.3f\n',T.iterations(s,m),T.time(s,m));
        end;
        T.scale_cv(s,m) = exp(th{m}(1,s));
        T.noise_cv(s,m) = exp(th{m}(2,s));
        if (strcmp(runEffect,'random'))
            T.run_cv(s,m)   = exp(th{m}(3,s));
        end;
    end;
end;
