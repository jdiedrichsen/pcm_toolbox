function [D,T,theta]=pcm_fitModelIndividCrossval(Y,M,partitionVec,conditionVec,varargin);
% function [D,T,theta]=pcm_fitModelIndividCrossval(Y,M,partitionVec,conditionVec,varargin);
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
%
%   'MaxIteration': Number of max minimization iterations. Default is 1000.
%
%   'S',S         : (Cell array of) NxN noise covariance matrices -
%                   otherwise independence is assumed
%   'Evaluation'  : Evaluation criteria of model: Either string or cell
%                   array (for multiple criteria)
%                   'crossvalR2': crossvalidated coefficient of determination (1-PSS/TSS)
%                   'crossvalR': correlation between predicted and observed
%                   'likelihood':  log-likelihood of observed under predicted
%                   'likelihood_rf': log-likelihood of observed under
%                                   predicted, with maximized signal and noise parameter
%--------------------------------------------------------------------------
% OUTPUT:
%   D:      Summary crossvalidation results (1 per subject)
%   T:      Detailed crossvalidation results (1 per condition)
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
MaxIteration    = 1100;
Iter            = [];
verbose         = 1;
S               = [];
pcm_vararginoptions(varargin,{'runEffect','MaxIteration','verbose','S'});

% Preallocate output structure
T.SN = [1:numSubj]';


% Get the number of models
if (~iscell(M))
    M{M};
end;
numModels   = numel(M);

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

% Loop over subject and provide inidivdual fits
for s = 1:numSubj
    
    % Set up crossvalidation scheme
    part = unique(partitionVec)';
    numPart = numel(part);
    if ischar(crossvalScheme)
        partI={};
        switch (crossvalScheme)
            case 'leaveOneOut'
                for i=1:numPart
                    partI{i}=part(i);
                end;
            case 'leaveTwoOut'
                for i=1:floor(numPart/2)
                    partI{i}=part((i-1)*2+[1:2]);
                end;
            case 'oddEven'
                for i=1:2
                    partI{i}=part(mod(part+i,2)==0);
                end;
        end;
    else
        partI = crossvalScheme; % Direct specificiation
    end;
    numFolds = numel(partI);
    
    % Now loop over models
    for m = 1:length(M)
        if (verbose)
            if isfield(M{m},'name');
                fprintf('Fitting Subj: %d model:%s\n',s,M{m}.name);
            else
                fprintf('Fitting Subj: %d model:%d\n',s,m);
            end;
        end;
        if (verbose)
            if isfield(M,'name');
                fprintf('Fitting Subj: %d model:%s\n',s,M{m}.name);
            else
                fprintf('Fitting Subj: %d model:%d\n',s,m);
            end;
        end;
        tic;
        
        % Set options
        if (~isempty(S))
            OPT.S = S;
        end;
        
        % Get starting guess for theta if not provided
        if (numel(theta0)<m)
            if (isfield(M{m},'theta0'))
                th0m = M{m}.theta0(1:M{m}.numGparams);
            else
                th0m = pcm_getStartingval(M{m},G_hat);
            end;
            
            if (strcmp(runEffect,'random'))
                x0 = [th0m;noise0;run0];
            else
                x0 = [th0m;noise0];
            end;
            theta0{m}=x0;
        end;
        
        % Set starting values
        x0 = theta0{m};
        
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
            
            switch (M{m}.fitAlgorithm)
                case 'minimize'  % Use minimize to find maximum liklhood estimate runEffect',B{s});
                    fcn = @(x) pcm_likelihoodIndivid(x,y*y',M{m},Zt,Xt,1,OPT);
                    [th(:,i),~,iter(i)] =  minimize(x0, fcn, MaxIteration);
                case 'NR'
                    fcn = @(x) pcm_likelihoodIndivid(x,y*y',M{m},Zt,Xt,1,OPT);
                    lik(i) = fcn(x0);
                    [th(:,i),~,iter(i)]= pcm_NR(x0,fcn,'verbose',verbose==2);
            end;
            
            
            
            
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
            
            % Get test data from that voxel
            y=Y(testI,p);
            OPT.runEffect=B(testI,:);
            Xt=reduce(X,testI);
            Zt=Z(testI,:);
            
            % Evaluate log-likelihood on left-out voxel
            switch (evalType)
                case 'uncond' % evaluates p(Y1 | theta)
                    lik(i) = -pcm_likelihoodIndivid(th(:,i),y*y',M{m},Zt,Xt,1,OPT);
                case {'simple','full'} % evaluate the prediction
                    lik(i) = pcm_crossvalLikelihood(M{m},th(:,i),Y(:,p),Z,X,...
                        trainI,testI,'type',evalType);
            end;
            
            % Use last iterations as a parameter starting value
            x0 = th(:,i);
        end;
        n=n+1;
    end; % For each partition
    % Summarize results across partitions for each subject
    indx = (T.SN==s);
    D.SN(s,1) = s;
    D.noise(s,:) = mean(T.noise(indx,:));
    D.likelihood(s,:) = sum(T.likelihood(indx,:)); % Partitions are independent
end; % for each subject
