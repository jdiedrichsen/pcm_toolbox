function out = pcm_searchlight_fcn_component(Y, partitionVec, conditionVec, MF, CompIdx, SPMs, varargin)
%% function out = pcm_searchlight_fcn_component(Y, partitionVec, conditionVec, MF, SPMs, varargin)
% A pcm-searchlight function for fitting component models using group-cross-validation.
%
% Inputs:
%	Y: A 1xNsubj cell array of timeseries data extracted from single
%		searchlight, where Nsubj is number of subjects. Each data is a TxP matrix.
%
%	partitionVec: A 1xNsubj cell array containing partition vectors.
%
%	conditionVec: A 1xNsubj cell array containing either condition vectors
%						 or design matrices (Z).
%
%	MF : Model family.
%
%   CompIdx: Component indices.
%
%	SPMs: A 1xNsubj cell array of SPM structure that specifies the 1st-level GLM analysis.
%
% Optional Inputs:
%	runEffect: The way how the run-effect is treated in fitting. Default is 'fixed'.
%
%	verbose: If set to 1, the program prints status in the command line. Default is 0.
%
%	fitAlgorithm: Algorithm for fitting. Default is 'NR'.
%
%	prior: Prior probability for calculating component posterior. Default is 0.5;
%
% a-yokoi (2018)

runEffect   = 'fixed'; % run effect
verbose     = 0;
fitAlgorithm= 'NR';
prior = 0.5;
calcCeiling = 1;
vararginoptions(varargin,{'runEffect', 'verbose','fitAlgorithm', 'prior','calcCeiling'});

%-----------------------------------------------------------------%
% Prewhiten Y (need partition and condition vectors for additional input)
Yprewh = cell(numel(Y),1);
for s=1:numel(Y)
    nCond       = length(partitionVec{s});
    %Yprewh{s}  = rsa_noiseNormalizeBeta(Y{s},SPMs(s));
    Yprewh{s}   = rsa.spm.noiseNormalizeBeta(Y{s},SPMs(s));
    Yprewh{s}   = Yprewh{s}(1:nCond,:);
end

% Calc noise ceiling if required
if calcCeiling
    if verbose; fprintf('Calculating noise ceiling...'); end;    
    M{1}.type = 'freedirect';
    M{1}.numGparams = 0;
    M{1}.theta0 = [];
    M{1}.name = 'noice_ceiling';
    % Run group-fit
    [Tupper, theta] = pcm_fitModelGroup(Yprewh, M, partitionVec, conditionVec,...
        'runEffect', runEffect, ...
        'fitAlgorithm', fitAlgorithm, ...
        'verbose', verbose);
    
    % Do the crossvalidated group-fit
    [Tlower] = pcm_fitModelGroupCrossval(Yprewh, M, partitionVec, conditionVec,...
        'runEffect', runEffect, ...
        'fitAlgorithm', fitAlgorithm, ...
        'groupFit', theta);
    if verbose; fprintf('...done.\n'); end;
end

% Run group-fit
[~, theta] = pcm_fitModelGroup(Yprewh, MF, partitionVec, conditionVec,...
    'runEffect', runEffect, ...
    'fitAlgorithm', fitAlgorithm, ...
    'verbose', verbose);

% Do the crossvalidated group-fit
[Tc, theta_cv] = pcm_fitModelGroupCrossval(Yprewh, MF, partitionVec, conditionVec,...
    'runEffect', runEffect, ...
    'fitAlgorithm', fitAlgorithm, ...
    'groupFit', theta);

% Get component posterior
[Tc.posterior, Tc.logBF] = pcm_componentPosterior(Tc.likelihood, CompIdx, 'prior', prior);

% Make output
out = [Tc.logBF, Tc.posterior, Tc.likelihood]; % [lower noise-ceiling,upper noise-ceiling,model likelihoods...]
if calcCeiling
    out = [out, Tupper.likelihood, Tlower.likelihood];
end

end