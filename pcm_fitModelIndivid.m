function [T,theta_hat,G_pred]=pcm_fitModelIndivid(Y,M,partitionVec,conditionVec,varargin);
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
%
%   'isCheckDeriv: Check the derivative accuracy of theta params. Done using
%                  'checkderiv'. This function compares input to finite
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
%       likelihood:         likelihood
%       noise:              Noise parameter 
%       run:                Run parameter (if run = 'random') 
%       iterations:         Number of interations for model fit
%       time:               Elapsed time in sec 
%
%   theta{m}     Cell array of estimated model parameters, each a 
%                 #params x #numSubj matrix 
%   G_pred{m}     Cell array of estimated G-matrices under the model 

runEffect       = 'random';
isCheckDeriv    = 0;
MaxIteration    = 1000;
Iter            = [];
verbose         = 1; 
S               = []; 
pcm_vararginoptions(varargin,{'runEffect','isCheckDeriv','MaxIteration','verbose','S'});

numSubj     = numel(Y);
numModels   = numel(M);

% Preallocate output structure
T.SN = [1:numSubj]';

% Determine optimal algorithm for each of the models 
M = pcm_optimalAlgorithm(M); 

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
    
    % Check if conditionVec is condition or design matrix
    if size(cV,2)==1;
        Z{s}   = pcm_indicatorMatrix('identity_p',cV);
    else
        Z{s} = cV;
    end;
    
    % Prepare matrices and data depending on how to deal with run effect 
    [N(s,1),P(s,1)] = size(Y{s});
    numCond= size(Z{s},2);
    YY{s}  = (Y{s} * Y{s}');
    switch (runEffect)
        case 'random'
            B{s}   = pcm_indicatorMatrix('identity_p',pV);
            X{s}   = [];
        case 'fixed'
            B{s}  =  [];
            X{s}  =  pcm_indicatorMatrix('identity_p',pV);
    end;
    
    % Estimate starting value run and the noise from a crossvalidated estimate of the second moment matrix 
    [G_hat(:,:,s),Sig_hat(:,:,s)] = pcm_estGCrossval(Y{s},pV,cV);
    sh = Sig_hat(:,:,s);
    run0(s)     = real(log((sum(sum(sh))-trace(sh))/(numCond*(numCond-1))));
    noise0(s)   = real(log(trace(sh)/numCond-exp(run0(s))));
    
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
        
        % Get starting guess for theta if not provided
        if (isfield(M{m},'theta0'))
            theta0 = M{m}.theta0(1:M{m}.numGparams);
        else
            theta0 = pcm_getStartingval(M{m},G_hat(:,:,s));   
        end;
        
        % if naive noise ceiling model, use crossvalidated G as component 
        if strcmp(M{m}.type,'noiseceiling')
            M{m}.Gc = pcm_makePD(G_hat(:,:,s)); 
        end; 
        
        % Set up overall starting values 
        switch (runEffect) 
            case {'fixed','remove'}
                x0  = [theta0;noise0(s)];
            case {'random'}
                x0  = [theta0;noise0(s);run0(s)];
        end; 
        
        % Now do the fitting, using the preferred optimization routine 
        switch (M{m}.fitAlgorithm)
            case 'minimize'  % Use minimize to find maximum liklhood estimate runEffect',B{s});
                if (isempty(S)) 
                    fcn = @(x) pcm_likelihoodIndivid(x,YY{s},M{m},Z{s},X{s},P(s),'runEffect',B{s});
                else 
                    fcn = @(x) pcm_likelihoodIndivid(x,YY{s},M{m},Z{s},X{s},P(s),'runEffect',B{s});
                end; 
                [theta_hat{m}(:,s),fX,T.iterations(s,m)]      =  ...
                        minimize(x0, fcn, MaxIteration);
                T.likelihood(s,m) =  -fX(end);  %invert the sign 
            case 'NR_diag'
                numDiag = size(M{m}.MM,2); 
                [~,theta_hat{m}(:,s),~,T.likelihood(s,m),T.iterations(s,m)] = ...
                    pcm_NR_diag(Y{s},Z{s}*M{m}.MM,'X',X{s},'Gd',ones(numDiag,1));
            case 'NR_comp'
                [~,theta_hat{m}(:,s),~,T.likelihood(s,m),T.iterations(s,m)] = ...
                    pcm_NR_comp(Y{s},Z{s},'X',X{s},'Gc',M{1}.Gc,'h0',x0);
            case 'NR' 
                if (isempty(S)) 
                    fcn = @(x) pcm_likelihoodIndivid(x,YY{s},M{m},Z{s},X{s},P(s),'runEffect',B{s});
                else 
                    fcn = @(x) pcm_likelihoodIndivid(x,YY{s},M{m},Z{s},X{s},P(s),'runEffect',B{s});
                end; 
                [theta_hat{m}(:,s),T.likelihood(s,m),T.iterations(s,m)]=pcm_NR(x0,fcn); 
        end; 
            
        G_pred{m}(:,:,s)  =  pcm_calculateG(M{m},theta_hat{m}(1:M{m}.numGparams,s));
        T.noise(s,m)      =  exp(theta_hat{m}(M{m}.numGparams+1,s)); 
        if strcmp(runEffect,'random')
            T.run(s,m)      =  exp(theta_hat{m}(M{m}.numGparams+2,s)); 
        end; 
        T.time(s,m)       = toc; 
        
        % This is an optional check if the dervivate calculation is correct
        if (isCheckDeriv)
            d = pcm_checkderiv(fcn,theta_hat{m}(:,s)-0.01,0.00001);
            fprintf('discrepency :%d\n',d);
        end;
    end; % for each model
end; % for each subject