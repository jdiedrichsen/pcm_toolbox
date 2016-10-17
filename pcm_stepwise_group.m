function varargout = pcm_stepwise_group(Data,Models,conditionVec,partitionVec,varargin);
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
%   Data: Prewhitened beta from one ROI for all subjects
%
%   Models: Candidate model structure array with following fields;
%       Compmodels:
%       Compparams:
%
%   Condition:  vector of condition

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

runeffect   = 1;    % remove run-effect
Z           = [];   % if specific shape of Z is needed
verbose     = 0;    % show command output

vararginoptions(varargin,{'verbose','runeffect','Z'});

%-----------------------------------------------------------------%
%-- Noise ceiling model
Ceiling.modelpred   = @(x)(mean(x,3)); % model for noise ceiling
Ceiling.numGparams  = 0;
Ceiling.x0          = [];
Ceiling.name        = 'Noise Ceiling';
%-----------------------------------------------------------------%


%-----------------------------------------------------------------%
% Get noise ceiling
%-----------------------------------------------------------------%
[Tceiling,Mceiling] = pcm_fitModelCrossval(Y,Ceiling,partitionVec,conditionVec,...
                        'runEffect','remove',...
                        'isCheckDeriv',1,...
                        'verbose',1);
%-----------------------------------------------------------------%

%-----------------------------------------------------------------%
% Adjust common starting values for component models
%-----------------------------------------------------------------%

%-- Define model components
models      = [1:length(Models)];
modeluse{1} = {[1]}; % starting from model-1 (assuming this is null-model)
compModels  = Models.models; % function names for component modeling
compNparams = Models.CompNparams; % number of parameters
compParams  = Models.CompParams;
compTheta0  = Models.theta0;

for i=1:numel(Models); % define model components
    cM(i).modelpred = Models.modelfcn{i};
    cM(i).numGparams= Models.compNparams(i);
    cM(i).params    = Models.compParams{i};
    cM(i).name      = Models.modelfcn{i};
end

%-----------------------------------------------------------------%
% Start forward selection
%-----------------------------------------------------------------%
% initial model
iter = 1;
modeluse{1} = {[1]}; % starting from null-model
Nmodel = numel(modeluse{iter}); M = [];
for i=1:Nmodel % model to be fitted
    M(i).modelpred  = @(x)(pcm_combineModels(x,cM(modeluse{iter}{i})));
    M(i).numGparams = sum(compNparams(modeluse{iter}{i}));
    M(i).x0         = cat(2,compTheta0{modeluse{iter}{i}})';
    M(i).name       = sprintf('%s.',compModels{modeluse{iter}{i}});
end

% loop start here
maxIter = 2^(length(cM)-1);
%-----------------------------------------------------------------%
% Iteerate pcm model fitting until converge
for iter = 1:maxIter
    if fig>0;
        h3=figure('name','2*log(K)',...
            'unit','centimeters','position',[20, 20, 12, 12]);
    end;
    
    if iter==1; isFinish = 0; end;
    if (remove)
        [T,theta_all,Gpred,theta,M] = pcm_fitModelCrossvalS(Y,M,Z,partitionVec,...
            'verbose',1,...
            'runEffect','remove',...
            'isCheckDeriv',1);
    else
        [T,theta_all,Gpred,theta] = pcm_fitModelCrossval(Y,M,conditionVec,partitionVec,...
            'isCheckIter',0,...
            'isCheckTime',1);
    end
    
    % Store result
    Tall{iter} = T;
    Mall{iter} = M;
    Gall{iter} = Gpred;
    
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
            
            modeluse{iter+1}{1} = [1:numel(compModels)];
            Nmodel = numel(modeluse{iter+1}); M = [];
            for i=1:Nmodel % model to be fitted
                M(i).modelpred  = @(x)(sh2.pcm_modelpred(x,cM(modeluse{iter+1}{i})));
                M(i).numGparams = sum(compNparams(modeluse{iter+1}{i}));
                M(i).x0         = cat(2,compTheta0{modeluse{iter+1}{i}})';
                M(i).name       = sprintf('%s.',compModels{modeluse{iter+1}{i}});
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
            M(i).modelpred  = @(x)(sh2.pcm_modelpred(x,cM(modeluse{iter+1}{i})));
            M(i).numGparams = sum(compNparams(modeluse{iter+1}{i}));
            M(i).x0         = cat(2,compTheta0{modeluse{iter+1}{i}})';
            M(i).name       = sprintf('%s.',compModels{modeluse{iter+1}{i}});
            
            % check derivative of model
            if (checkderiv==1&&M(i).numGparams>0&&iter<3)
                d = sh2_simulate('check_model',...
                    'model','sh2.pcm_modelpred',...
                    'theta0',M(i).x0,...
                    'P',{cM(modeluse{iter+1}{i})});
                if any(d>1e-5)
                    warning('derivative of model%d is incorrect!',i);
                end
            end
        end
        
    end
    
    % end of loop here
end
if fig>0; close(h3);end

%-----------------------------------------------------------------%
% Save result
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
T.SN    = sn';
T.reg   = repmat(reg,size(T.SN));
T.hemis = repmat(h,size(T.SN));
history.allmodels = itermodels; % history of forward selection steps
history.logBF = logBF;
history.bestiter = bestiter;

fname = sprintf('%s.hem%d.reg%d.glm%d.AB%s%s.mat',what,h,reg,glm,removename{remove+1},option);
save(fullfile(regDir,fname),...
    'T','C','M','cM','modeluse','history','Tall','Mall','Gall','this');
clear modeluse Tall Mall Gall logBF


