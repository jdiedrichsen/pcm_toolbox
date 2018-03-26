function [T,DD,theta_hat,theta0]=pcm_fitModelIndividCrossval(Y,M,partitionVec,conditionVec,varargin);
% function [T,D,theta_hat]=pcm_fitModelIndividCrossval(Y,M,partitionVec,conditionVec,varargin);
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
%   'fitAlgorithm': Either 'NR' or 'minimize' - provides over-write for
%                   model specific algorithms
%   'verbose':      Optional flag to show display message in the command
%                   line. Default is 1. Setting to 2 gives more detailed
%                   feedback on pcm_NR
%   'theta0':       Cell array of starting values (same format as theta{m})
%   'evaluation'  : Evaluation criteria of model: Either string or cell
%                   array (for multiple criteria)
%                   'crossvalR2': crossvalidated coefficient of determination (1-PSS/TSS)
%                   'crossvalR': correlation between predicted and observed
%                   'likelihood':  log-likelihood of observed under predicted
%                   'likelihood_rf': log-likelihood of observed under
%                                   predicted, with maximized signal and noise parameter
%--------------------------------------------------------------------------
% OUTPUT:
%   T:      Summary crossvalidation results (1 per subject)
%       SN:                 Subject number
%       likelihood:         crossvalidated likelihood
%       noise:              Noise parameter
%       run:                Run parameter (if run = 'random')
%       iterations:         Number of interations for model fit
%       time:               Elapsed time in sec
%   D:  Detailed crossvalidation results (1 per condition)    
%       fold:               Crossvalidation fold 
%
%  theta_hat{model}:  Estimated parameters (model + scaling/noise parameters)
%                       for individual subjects 

runEffect       = 'random';
isCheckDeriv    = 0;
MaxIteration    = 1000;
Iter            = [];
verbose         = 1;
S               = [];
fitAlgorithm    = [];
crossvalScheme  = 'leaveTwoOut'; 
evaluation      = {'R2','R'}; 
theta0          = {}; 
pcm_vararginoptions(varargin,{'crossvalScheme','fitAlgorithm','runEffect',...
            'MaxIteration','verbose','S','evaluation','crossvalScheme','theta0'});

DD=[]; 
% Number of subejcts 
if (~iscell(Y))
    Y={Y};
end; 
numSubj     = numel(Y);

% Number of evaluation criteria 
if (~iscell(evaluation)) 
    evaluation={evaluation}; 
end; 
numEval = numel(evaluation); 

% Preallocate output structure
T.SN = [1:numSubj]';

% Get the number of models
if (~iscell(M))
    M={M};
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
    if iscell(partitionVec) 
        pV=partitionVec{s}; 
    else 
        pV=partitionVec; 
    end; 
    part = unique(pV)';
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
        tic;
        
        % Set options
        if (~isempty(S))
            OPT.S = S;
        end;
        
        % Get starting guess for theta if not provided
        if (numel(theta0)<m || size(theta0{m},2)<s) 
            if (isfield(M{m},'theta0'))
                th0m = M{m}.theta0(1:M{m}.numGparams);
            else
                th0m = pcm_getStartingval(M{m},G_hat(:,:,s));   
            end;
        
            if (strcmp(runEffect,'random'))
                theta0{m}(:,s) = [th0m;noise0(s);run0(s)];
            else
                theta0{m}(:,s) = [th0m;noise0(s)]; 
            end;
        end; 

        % Set starting values
        x0 = theta0{m}(:,s);
        th = nan(size(x0,1),numFolds); 
        
        % Now perform the cross-validation across different partitions
        for p=1:numFolds
            trainIdx = ~ismember(pV,partI{p}); 
            testIdx = ismember(pV,partI{p}); 
            
            % Get the data and design matrices for training set 
            Ytrain=Y{s}(trainIdx,:); 
            Xtrain=reduce(X{s},trainIdx);
            Ztrain=Z{s}(trainIdx,:);
            OPT.runEffect = reduce(B{s},trainIdx); 

            % Get test data and design matrices for the test set 
            Ytest=Y{s}(testIdx,:);
            Xtest=reduce(X{s},testIdx);
            Ztest=Z{s}(testIdx,:);
            Btest=reduce(B{s},testIdx);
            
            % Perform the initial fit to the training data 
            switch (M{m}.fitAlgorithm)
                case 'minimize'  % Use minimize to find maximum liklhood estimate runEffect',B{s});
                    fcn = @(x) pcm_likelihoodIndivid(x,Ytrain*Ytrain',M{m},Ztrain,Xtrain,P(s),OPT);
                    l=fcn(x0); 
                    [th(:,p),fX,D.iterations(p,m)] =  minimize(x0, fcn, MaxIteration);
                    D.likelihood_fit(p,m)=-fX(end); 
                case 'NR'
                    fcn = @(x) pcm_likelihoodIndivid(x,Ytrain*Ytrain',M{m},Ztrain,Xtrain,P(s));
                    [th(:,p),D.likelihood_fit(p,m),D.iterations(p,m)]= pcm_NR(x0,fcn,'verbose',verbose);
            end;
            
            % Record the stats from fitting
            D.SN(p,1)         = s;
            D.fold(p,1)       = p;
            D.noise(p,m)      =  exp(th(M{m}.numGparams+1,p));
            if strcmp(runEffect,'random')
                D.run(p,m)      =  exp(th(M{m}.numGparams+2,p));
            end;
            D.time(p,m)       = toc;
            
            % calculate prediction on 
            estU = pcm_estimateU(M{m},th,Ytrain,Ztrain,Xtrain,'runEffect',OPT.runEffect); 
            Ypred  = Ztest*estU;
            Ytestx = Ytest-Xtest*pinv(Xtest)*Ytest;
            for c = 1:numEval 
                switch (evaluation{c})
                    case 'likelihood_uncond' % evaluates p(Y1 | theta)
                        lik(i) = -pcm_likelihoodIndivid(th(:,p),y*y',M{m},Zt,Xt,1,OPT);
                    case 'likelihood_cond' % evaluate the prediction
                        D.likelihood_cond = pcm_crossvalLikelihood(M{m},th(:,p),Y(:,p),Z,X,...
                            trainI,testI,'type',evalType);
                    case 'R2'              % Predictive R2 
                        D.TSS(p,m)= sum(sum(Ytestx.*Ytestx));
                        D.RSS(p,m)= sum(sum((Ytestx-Ypred).^2)); 
                    case 'R'               % Predictive correlation 
                        D.SS1(p,m) = sum(sum(Ytestx.*Ytestx));
                        D.SS2(p,m) = sum(sum(Ypred.*Ypred));
                        D.SSC(p,m) = sum(sum(Ypred.*Ytestx));
                end; 
            end;
            
            % Use last iterations as a parameter starting value
            x0 = th(:,p);
        end;                % For each partition 
        theta_hat{m}(:,s)=mean(th,2);
    end;                    % For each model 
    DD=addstruct(DD,D); 
    
    % Summarize results across partitions for each subject
    T.noise(s,:)=mean(D.noise); 
    if strcmp(runEffect,'random')
        T.run(s,:) =  mean(D.run);
    end;
    T.time(s,:)    =  sum(D.time);
    T.iterations(s,:)    =  sum(D.iterations);
    for c = 1:numEval 
        switch (evaluation{c})
            case {'likelihood_uncond','likelihood_cond'}
                T.(evaluation{c})(s,:)    =  sum(D.(evaluation{c}));
            case 'R2'
                TSS = sum(D.TSS); 
                RSS = sum(D.RSS); 
                T.R2(s,:)    = 1-RSS./TSS;  
            case 'R'
                SSC = sum(D.SSC); 
                SS1 = sum(D.SS1); 
                SS2 = sum(D.SS2); 
                T.R(s,:)    = SSC./sqrt(SS1.*SS2); 
        end;
    end; 
end; % for each subject


function Xt=reduce(X,index); 
    Xt=X(index,:);
    Xt=Xt(:,sum(abs(Xt))>0); 
