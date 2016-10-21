function [compLike,ceilings,models] = pcm_stepwise_getCompLike(Tall,Ce,history);
%% function compLike = pcm_stepwise_getCompLike(T,C,history);
% Calculates component log-likelihood from the result of pcm_stepwise_group
% 
% 
% 
% 
% 
% 
% a-yokoi (2016)

%=========================================================%
% Choose result based on history of foward selection and
%  calculate relative log-likelihood
%=========================================================%
likelihood      = [];
likelihood_all  = [];
bestmodels      = [];
dModel          = [];
Niter           = min(history.maxIter,history.bestiter+1);
allmodels       = unique(history.allmodels{end}{1});
for i=1:Niter
    logBFs = history.logBF{i};
    models = history.allmodels{i};
    
    [maxlogBF(i),idx]   = max(logBFs);
    
    bestmodels{end+1}   = models{idx};
    likelihood          = [likelihood,Tall{i}.likelihood(:,idx)];
    likelihood_all      = [likelihood_all,Tall{i}.likelihood_all(:,idx)];
end
for i=1:Niter-1
    dModel(i) = setdiff(bestmodels{i+1},bestmodels{i});
end
bestmodel   = bestmodels{history.bestiter};
missing     = setdiff(allmodels(2:end),dModel);    % missing model (allmodel(1) is null model)
dLike       = [diff(likelihood,1,2),...
    NaN(size(likelihood,1),length(missing))];
dLike_all   = [diff(likelihood_all,1,2),...
    NaN(size(likelihood_all,1),length(missing))];
[dModel,id] = sort([dModel,missing]);   % add missing model and sort

% re-order data
dLike       = [dLike(:,id)];
dLike_all   = [dLike_all(:,id)];
dLike       = [dLike,nansum(dLike(:,ismember(dModel,bestmodel)),2)];
dLike_all   = [dLike_all,nansum(dLike_all(:,ismember(dModel,bestmodel)),2)];


%=========================================================%
% Noise ceiling
%=========================================================%
mixLike     = [Tall{history.bestiter}.likelihood,Ce.likelihood];
mixLike_all = [Tall{history.bestiter}.likelihood_all,Ce.likelihood_all];
[maxupper,idx] = max(mixLike_all,[],2);
Ce.upper = maxupper - Tall{1}.likelihood;
for i=1:numel(idx)
    Ce.lower(i,1) = mixLike(i,idx(i))- Tall{1}.likelihood(i,1);
end

compLike = dLike;
ceilings = [Ce.lower,Ce.upper];
models   = [dModel,-1]; % -1: winning model



end