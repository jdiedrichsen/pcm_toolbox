function [T,M,C,history,Tall] = pcm_stepwise_group(Y,Models,conditionVec,partitionVec,varargin);
%% function varargout = pcm_stepwise_group(Data,Model,varargin);
% Performs stepwise regression based on Bayes factor.
%
% 0. get noise ceiling and full-model first (could be at last)
% 1. start from single models
% 2. choose the best model
% 3. add a model from unused ones
% 4. repeat 2 and 3 untill the 2*ln(Bayes factor) between previous
%    best model gets smaller than 2
% The model for step N-1 is the best, good-enough model
%
% 2*ln K	K           Strength of evidence
%----------------------------------------------------------
% 0 to 2    1 to 3      not worth more than a bare mention
% 2 to 6    3 to 20     positive
% 6 to 10   20 to 150   strong
% >10       >150        very strong
%                                  (Kass and Raftery (1995))
%
% Input:
%   Y: Prewhitened beta from one ROI for all subjects
%
%   Models: Candidate model structure array with following fields;
%       modelpred:  function handle to component model
%       numGparams: number of parameters to be optimized
%       params:     additional fixed parameters for the model function
%       theta0:     initial parameter values
%       name:       model name
%
%   Condition:  vector of condition
%
%   Partition:  vector of membership for session
%
% Optional inputs:
%
%
% Output:
%
%
%
%
% Atsushi Yokoi (2016)

runEffect   = 'remove';    % run effect
Z           = [];   % if specific shape of Z is needed
verbose     = 0;    % show command output
fig         = 0;    % show history of Bayes factor
selection   = 'forward'; % 'backward'
vararginoptions(varargin,{'verbose','runEffect','Z','fig'});

%-----------------------------------------------------------------%
%-- Noise ceiling model
Ceiling.modelpred   = @(x)(mean(x,3)); % model for noise ceiling
Ceiling.numGparams  = 0;
Ceiling.x0          = [];
Ceiling.name        = 'Noise Ceiling';
Ceiling.type        = 'noiseceiling';
%-----------------------------------------------------------------%


%-----------------------------------------------------------------%
% Get noise ceiling first
%-----------------------------------------------------------------%
switch (runEffect)
    case 'remove'
        [C] = pcm_fitModelCrossvalS(Y,Ceiling,Z,partitionVec,...
            'verbose',verbose,...
            'runEffect','remove',...
            'isCheckDeriv',verbose);
    otherwise
        [C] = pcm_fitModelCrossval(Y,Ceiling,partitionVec,conditionVec,...
            'runEffect',runEffect,...
            'isCheckDeriv',verbose,...
            'verbose',verbose,...
            'Z',[]);
end
%-----------------------------------------------------------------%

switch (selection)
    case 'forward'
        %-----------------------------------------------------------------%
        % Start forward selection
        %-----------------------------------------------------------------%
        % initial model
        models      = [1:length(Models)];
        iter        = 1;
        modeluse{1} = {[1]}; % starting from null-model
        Nmodel      = numel(modeluse{iter}); M = [];
        for i=1:Nmodel % model to be fitted
            M(i).modelpred  = @(x)(pcm_combineModels(x,Models(modeluse{iter}{i})));
            M(i).numGparams = sum(cat(2,Models(modeluse{iter}{i}).numGparams));
            M(i).x0         = cat(2,Models(modeluse{iter}{i}).theta0)';
            M(i).name       = sprintf('%s.',Models(modeluse{iter}{i}).name);
        end
        
        % loop start here
        maxIter = 2^(length(Models)-1);
        %-----------------------------------------------------------------%
        % Iterate pcm model fitting until converge
        for iter = 1:maxIter
            if fig>0;
                h3=figure('name','2*log(K)',...
                    'unit','centimeters','position',[20, 20, 12, 12]);
            end;
            
            if iter==1; isFinish = 0; end;
            switch (runEffect)
                case 'remove'
                    [T,~,~,~,M] = pcm_fitModelCrossvalS(Y,M,Z,partitionVec,...
                        'verbose',verbose,...
                        'runEffect','remove',...
                        'isCheckDeriv',verbose);
                otherwise
                    [T,M] = pcm_fitModelCrossval(Y,M,partitionVec,conditionVec,...
                        'runEffect',runEffect,...
                        'isCheckDeriv',verbose,...
                        'verbose',verbose,...
                        'Z',[]);
            end
            
            % Store result
            Tall{iter} = T;
            Mall{iter} = M;
            %Gall{iter} = Gpred;
            
            % Calc 2log Bayes factor (2*ln(K)) and select best model
            logBF{iter} = 2*nanmean(T.likelihood,1);
            [maxLogBF(iter),maxindx] = max(logBF{iter});
            
            % Compare with previous best model
            if iter>1;
                dMaxLogBF = maxLogBF(iter)-maxLogBF(iter-1);
            else
                dMaxLogBF = maxLogBF(iter);
            end
            bestmodel   = modeluse{iter}{maxindx}; % choose best model at current stage
            unusedmodel = setdiff(models,bestmodel); % find unused models
            
            if (fig>0)
                figure(h3);plot(maxLogBF,'k-o');
                xlabel('Iteration');ylabel('2*log(K) for best model');
                drawnow;
            end
            
            %-----------------------------------------------------------------%
            % Update models
            %
            % If 2*log-Baysfactor of current best model
            %  - vs. the previous best model
            % is smaller than 2, fit full-model and end selection.
            if (dMaxLogBF<2&&iter>1)
                if isFinish==0
                    bestiter = iter-1;
                end
                isFinish = 1;
                fprintf('Stop criteria satisfied (no further improvement in Bayes factor).\n');
                
                if isempty(unusedmodel)
                    fprintf('Finish estimation.\n');
                    break;
                else
                    fprintf('Going to full-model...\n');
                    
                    modeluse{iter+1}{1} = [1:numel(Models)];
                    Nmodel = numel(modeluse{iter+1}); M = [];
                    for i=1:Nmodel % model to be fitted
                        M(i).modelpred  = @(x)(sh2.pcm_modelpred(x,Models(modeluse{iter+1}{i})));
                        M(i).numGparams = sum(cat(2,Models(modeluse{iter+1}{i}).numGparams));
                        M(i).x0         = cat(2,Models(modeluse{iter+1}{i}).theta0)';
                        M(i).name       = sprintf('%s.',Models(modeluse{iter+1}{i}).name);
                    end
                end
                % If 2*log-Baysfactor of current best model over the
                % previous best model is larger than 2, add unused
                % model to the current best model and continue.
            else
                if isFinish==0
                    bestiter = iter;
                end
                if isempty(unusedmodel) % if full-model is already done, stop here
                    fprintf('Stop criteria satisfied (reached to full-model).\n');
                    break;
                end
                
                for m=1:numel(unusedmodel)
                    modeluse{iter+1}{m} = [bestmodel,unusedmodel(m)];
                end
                Nmodel = numel(modeluse{iter+1}); M = [];
                for i=1:Nmodel % model to be fitted
                    M(i).modelpred  = @(x)(sh2.pcm_modelpred(x,Models(modeluse{iter+1}{i})));
                    M(i).numGparams = sum(cat(2,Models(modeluse{iter+1}{i}).numGparams));
                    M(i).x0         = cat(2,Models(modeluse{iter+1}{i}).theta0)';
                    M(i).name       = sprintf('%s.',Models(modeluse{iter+1}{i}).name);
                end
            end
            
            % end of loop here
        end
    case 'backward'
        error('Not implemented.');
end

if fig>0; close(h3);end


%-----------------------------------------------------------------%
% Summarise result and history of selection
% (model: null, single models, best + competing models, full- model)
%-----------------------------------------------------------------%
clear T M; T = []; M = [];
iteruse = unique([1,2,bestiter,iter]);
itermodels = modeluse; modeluse = [];
for i=iteruse
    T = addstruct(T,Tall{i},'column');
    M = cat(2,M,Mall{i});
    for m=1:numel(itermodels{i})
        modeluse = cat(2,modeluse,{itermodels{i}{m}});
    end
end

% history
history.allmodels   = itermodels; % history of forward selection steps
history.logBF       = logBF;
history.bestiter    = bestiter;
history.modeluse    = modeluse;
history.maxIter     = maxIter;

