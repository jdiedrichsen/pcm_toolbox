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
%        M: Cell array (number of models) of Models Each is a structure with subfields.
%              .type:        Type of the model to be fitted
%                             'fixed': Fixed structure without parameter
%                                (except scale and noise for each subject)
%                             'component': G is a sum of linear components,
%                                specified by Gc
%                             'feature': G=A*A', with A a linear sum of
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
%              .Ac           Linear component matrices (for type 'feature')
%
%   partitionVec: {#Subjects} Cell array with partition assignment vector
%                   for each subject. Rows of partitionVec{subj} define
%                   partition assignment of rows of Y{subj}.
%                   Commonly these are the scanning run #s for beta
%                   regressors.If 1 vector is given, then partitions are
%                   assumed to be the same across subjects
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
%   'fitScale'      Introduce additional scaling parameter for each
%                   participant? - default is true
%   'isCheckDeriv:  Check the derivative accuracy of theta params. Done using
%                   'pcm_checkderiv'. This function compares input to finite
%                   differences approximations. See function documentation.
%
%   'MaxIteration': Number of max minimization iterations. Default is 1000.
%
%   'verbose':      Optional flag to show display message in the command
%                   line (e.g., elapsed time). Default is 1.
%   'groupFit',T:   Structure T from the group fit: This provides better starting
%                   values and can speed up the computation
% 
%   'Z':            User specified design matrix Z for flexible modelling.
%                   
%--------------------------------------------------------------------------
% OUTPUT:
%   T:      Structure with following subfields:
%       SN:                 Subject number
%       likelihood:         Crossvalidated likelihood
%       scale:              Fitted scaling parameter - exp(theta_scale)
%       noise:             	Fitted noise parameter - exp(theta_noise)
%       run:               	Fitted run parameter - exp(theta_run)
%
%   M{m}:    Cell array of models - with appended fields
%       theta:      numGParams x numSubj Matrix: Estimated parameters (model + scaling/noise parameters)
%                   across the crossvalidation rounds.

runEffect       = 'random';
isCheckDeriv    = 0;
MaxIteration    = 1000;
verbose         = 1;
groupFit        = [];
fitScale        = 1;
Z               = [];
pcm_vararginoptions(varargin,{'runEffect','isCheckDeriv','MaxIteration',...
    'verbose','groupFit','fitScale','Z'});

numSubj     = numel(Y);
numModels   = numel(M);

% Preallocate output structure
T.SN = [1:numSubj]';

% Check size of Z
if ~isempty(Z)&&numel(Z)~=numSubj;
    error('Invalid size of user-specified "Z" matrix!');
end

% --------------------------------------------------------
% Figure out a starting values for the noise parameters
% --------------------------------------------------------
for s = 1:numSubj
    
    % Get Condition and partition vector
    if (iscell(conditionVec))
        cV = conditionVec{s};
        pV = partitionVec{s};
    else
        cV = conditionVec;
        pV = partitionVec;
    end;
    
    % Set up the main matrices
    [N(s,1),P(s,1)] = size(Y{s});
    if isempty(Z)||numel(Z)<numSubj
        Z{s}   = pcm_indicatorMatrix('identity_p',cV);
    end
    numCond= size(Z{s},2);
    
    % Depending on the way of dealing with the run effect, set up data
    switch (runEffect)
        case 'random'
            YY{s}  = (Y{s} * Y{s}');
            B{s}   = pcm_indicatorMatrix('identity_p',pV);
            X{s}   = [];
            S      = [];   % Use identity for covariance
        case 'fixed'
            YY{s}  = (Y{s} * Y{s}');
            B{s}  =  [];
            X{s}  =  pcm_indicatorMatrix('identity_p',pV);
            S     =  [];    % Use identity for covariance
        case 'remove'
            Run         =  indicatorMatrix('identity_p',pV);
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
    [G_hat(:,:,s),Sig_hat(:,:,s)] = pcm_estGCrossval(Y{s},pV,cV);
    sh              = Sig_hat(:,:,s);
    run0(s,1)       = log((sum(sum(sh))-trace(sh))/(numCond*(numCond-1)));
    noise0(s,1)     = log(trace(sh)/numCond-exp(run0(s)));
    
    if (fitScale)
        G0            = mean(G_hat,3);
        g0            = G0(:);
        g_hat         = G_hat(:,:,s);
        g_hat         = g_hat(:);
        scaling       = (g0'*g_hat)/(g0'*g0);
        if ((scaling<10e-6)||~isfinite(scaling));
            scaling = 10e-6;
        end;      % Enforce positive scaling
        scale0(s,1)   = log(scaling);
    end;
    
end;

% -----------------------------------------------------
% Loop over subjects and obtain crossvalidated likelihoods
% -----------------------------------------------------
for s = 1:numSubj
    % determine training set
    notS    = [1:numSubj];
    notS(s) = [];
    
    % Now loop over models
    for m = 1:numModels
        if (verbose)
            if isfield(M{m},'name');
                fprintf('Crossval Subj: %d model:%s',s,M{m}.name);
            else
                fprintf('Crossval Subj: %d model:%d',s,m);
            end;
        end;
        tic;
        
        % Determine starting values for fit: If group fit is given, start
        % with those values
        if (~isempty(groupFit))
            noise0 = log(groupFit.noise(:,m));
            if (fitScale)
                scale0 = log(groupFit.scale(:,m));
            end;
            if (strcmp(runEffect,'random'))
                run0 = log(groupFit.run(:,m));
            end;
            theta0 = M{m}.thetaGroup;
        else % No group fit: determine starting values
            if (isfield(M{m},'theta0'))
                theta0 = M{m}.theta0;
            else
                theta0 = pcm_getStartingval(M{m},mean(G_hat,3));
            end;
        end;
        
        % Now fit to all the subject but the left-out one
        switch (M{m}.type)
            case 'fixed'
                if size(M{m}.Gc,3)>1
                    G = mean(M{m}.Gc(:,:,notS),3);
                else
                    G = M{m}.Gc;
                end
                i = 0;
            case 'noiseceiling'
                G = mean(G_hat(:,:,notS),3);    % uses the mean of all other subjects
                G = pcm_makePD(G);
                i = 0;
            otherwise
                if (isempty(S))
                    fcn = @(x) pcm_likelihoodGroup(x,{YY{notS}},M{m},{Z{notS}},{X{notS}},P(notS(:)),...
                        'runEffect',{B{notS}},'fitScale',fitScale);
                else
                    fcn = @(x) pcm_likelihoodGroup(x,{YY{notS}},M{m},{Z{notS}},{X{notS}},P(notS(:)),...
                        'runEffect',{B{notS}},'S',S(notS),'fitScale',fitScale);
                end;
                
                % Generate the starting vector
                x0 = [theta0;noise0(notS)];
                if (fitScale)
                    x0 = [x0;scale0(notS)];
                end;
                if (strcmp(runEffect,'random'))
                    x0  = [x0;run0(notS)];       % Start with G-params from group fit
                end;
                [theta,fX,i] =  minimize(x0, fcn, MaxIteration);
                M{m}.thetaCross(:,s)=theta(1:M{m}.numGparams);
                G   = pcm_calculateG(M{m},M{m}.thetaCross(1:M{m}.numGparams,s));
        end;
        
        % Now get the fit the left-out subject
        % (maximizing scale and noise coefficients)
        x0 = noise0(s);
        if fitScale
            x0 = [x0;scale0(s)];
        end;
        if (strcmp(runEffect,'random'))
            x0      = [x0;run0(s)];
        end;
        
        if (isempty(S))
            fcn     = @(x) pcm_likelihoodGroup(x,{YY{s}},G,{Z{s}},{X{s}},P(s),...
                'runEffect',{B{s}},'fitScale',fitScale);   % Minize scaling params only
        else
            fcn     = @(x) pcm_likelihoodGroup(x,{YY{s}},G,{Z{s}},{X{s}},P(s),...
                'runEffect',{B{s}},'S',S(s),'fitScale',fitScale);   % Minize scaling params only
        end;
        
        [th{m}(:,s),fX] =  minimize(x0, fcn, MaxIteration);
        
        T.likelihood(s,m) = -fX(end);
        T.iterations(s,m) = i;
        T.time(s,m)       = toc;
        if verbose
            fprintf('\t Iterations %d, Elapsed time: %3.3f\n',T.iterations(s,m),T.time(s,m));
        end;
        T.noise(s,m) = exp(th{m}(1,s));
        if (fitScale)
            T.scale(s,m) = exp(th{m}(2,s));
        end;
        if (strcmp(runEffect,'random'))
            T.run(s,m)   = exp(th{m}(fitScale+1,s));
        end;
    end;
end;
