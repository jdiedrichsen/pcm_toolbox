function [T,theta]=pcm_fitModelIndividCrossval(Y,M,partitionVec,conditionVec,varargin);
% function [T,theta]=pcm_fitModelIndividCrossval(Y,M,partitionVec,conditionVec,varargin);
% Fits pattern component model(s) specified in M to data from one (or more)
% subjects individually, using leave-one out crossvalidation within each
% subject.
%==========================================================================
% INPUT:
%        Y: [#Conditions x #Voxels]
%            Observed/estimated beta regressors from each subject.
%            Preferably multivariate noise-normalized beta regressors.
%            If it's a cell array, it is assuming multiple subjects, each
%            cell containing the data from one subject
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
%   partitionVec: Partition assignment vector
%                   Could be a cell array of multiple vector if its for multiple
%                   subjects and the parition structure is different.
%                   Rows of partitionVec{subj} defin partition assignment of rows of Y{subj}.
%                   Commonly these are the scanning run #s for beta
%                   regressors.
%                   If a single vector is provided, it is assumed to me the
%                   same for all subjects
%                   The runs are assume to be labeled 1-numRuns
%
%   conditionVec:  Condition assignment vector
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
%   'MaxIteration': Number of max minimization iterations. Default is 1000.
%   'S',S         : Structure of the NxN noise covariance matrix -
%                   otherwise it assumes independence here
%
%--------------------------------------------------------------------------
% OUTPUT:
%   T:      Detailed crossvalidation results
%       SN:                 Subject number
%       likelihood:         crossvalidated likelihood
%       partition:          Partition serving as test set
%       noise:              Noise parameter
%       run:                Run parameter (if run = 'random')
%       iterations:         Number of interations for model fit
%       time:               Elapsed time in sec
%
%  theta{model}:      Estimated parameters (model + scaling/noise parameters)
%                       for individual subjects and partitions

runEffect       = 'random';
isCheckDeriv    = 0;
MaxIteration    = 1000;
Iter            = [];
verbose         = 1;
S               = [];
pcm_vararginoptions(varargin,{'runEffect','MaxIteration','verbose','S'});

if (~iscell(Y))
    Y={Y};
end;
numSubj     = numel(Y);
if (~iscell(M))
    M{M};
end;
numModels   = numel(M);

n=1; % Result index
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
    
    if (isempty(S))
        SS=[];
    else
        SS=S{s};
    end;
    
    % Prepare matrices and data depnding on how to deal with run effect
    [N,P] = size(Y{s});
    Z     = pcm_indicatorMatrix('identity_p',cV);
    numCond= size(Z,2);
    switch (runEffect)
        case 'random'
            YY  = (Y{s} * Y{s}');
            B   = pcm_indicatorMatrix('identity_p',pV);
            X   = zeros(N,0);  % Make an indexable, empty matrix
        case 'fixed'
            YY  = (Y{s} * Y{s}');
            B  =  zeros(N,0);  % Make an indexable, empty matrix
            X  =  pcm_indicatorMatrix('identity_p',pV);
    end;
    
    % Estimate starting value run and the noise from a crossvalidated estimate of the second moment matrix
    [G_hat,Sig_hat] = pcm_estGCrossval(Y{s},pV,cV);
    sh = Sig_hat;
    run0     = real(log((sum(sum(sh))-trace(sh))/(numCond*(numCond-1))));
    noise0   = real(log(trace(sh)/numCond-exp(run0)));
    
    % Get starting values over models
    for m = 1:length(M)
        % Get starting guess for theta if not provided
        if (isfield(M{m},'theta0'))
            theta0{m} = M{m}.theta0;
        else
            theta0{m} = pcm_getStartingval(M{m},G_hat);
        end;
    end;
    
    % Now perform the cross-validation across different partitions
    numPart = max(pV);
    for p=1:numPart
        
        testIdx = pV==p;
        traiIdx = pV~=p;
        
        % Now set up the function that returns likelihood and derivative
        Xtest = X(testIdx,:);
        Xtest = Xtest(:,sum(abs(Xtest))>0);  % remove unecessary columns
        Xtrai = X(traiIdx,:);
        Xtrai = Xtrai(:,sum(abs(Xtrai))>0);  % remove unecessary columns
        Btest = B(testIdx,:);
        Btest = Btest(:,sum(abs(Btest))>0);  % remove unecessary columns
        Btrai = B(traiIdx,:);
        Btrai = Xtest(:,sum(abs(Btrai))>0);  % remove unecessary columns
        
        % Now loop over models
        for m = 1:length(M)
            if (verbose)
                if isfield(M,'name');
                    fprintf('Fitting Subj: %d model:%s\n',s,M{m}.name);
                else
                    fprintf('Fitting Subj: %d model:%d\n',s,m);
                end;
            end;
            tic;
            
            
            if (isempty(SS))
                fcn = @(x) pcm_likelihoodIndivid(x,YY(traiIdx,traiIdx),M{m},Z(traiIdx,:),Xtrai,P,'runEffect',Btrai);
            else
                fcn = @(x) pcm_likelihoodIndivid(x,YY(traiIdx,traiIdx),M{m},Z(traiIdx,:),Xtrai,P,'runEffect',Btrai,'S',SS(traiIdx,traiIdx));
            end;
            
            % Set up overall starting values
            switch (runEffect)
                case {'fixed','remove'}
                    x0  = [theta0{m};noise0];
                case {'random'}
                    x0  = [theta0{m};noise0;run0T];
            end;
            
            % Use minimize to fine maximum liklhood estimate
            [th,fX,i]      =  minimize(x0, fcn, MaxIteration);
            
            % Record the stats from fitting
            T.SN(n,1)         = s;
            T.partition(n,1)  = p;
            T.noise(n,m)      =  exp(th(M{m}.numGparams+1));
            if strcmp(runEffect,'random')
                T.run(n,m)      =  exp(th(M{m}.numGparams+2));
            end;
            T.iterations(n,m) = i;
            T.time(n,m)       = toc;
            theta{m}(n,:)     = th(1:M{m}.numGparams)';
            
            % Evaluation criterion: Simple log-likelihood
            if (isempty(SS))
                T.likelihood(n,m) =  pcm_likelihoodIndivid(th,YY(testIdx,testIdx),M{m},...
                    Z(testIdx,:),Xtest,P,'runEffect',Btest);
            else
                T.likelihood(n,m) =  pcm_likelihoodIndivid(th,YY(testIdx,testIdx),M{m},...
                    Z(testIdx,:),Xtest,P,'runEffect',Btest,'S',SS(testIdx,testIdx));
            end;
            %             [U,G,iV]=pcm_estimateU(th,Y{s}(traiIdx,:),M{m},...
            %                 Z(traiIdx,:),Xtrai,'runEffect',Btrai);
            %             keyboard;
        end; % for each model
        n=n+1;
    end; % For each partition
end; % for each subject
