function [postProb,logBayes]=pcm_componentPosterior(likelihood,compI,varargin); 
% function p=pcm_componentPosterior(likelihood,compI); 
% Calculates a posterior probability for the presence of each 
% component from a family of models 
% INPUT: 
%   likelihood: (numSubj x numModels) estimated marginal log-likelihood of each model. This can be 
%               crossvalidated Type-II maximum likelihood or Type-II maximum likelihood 
%               corrected with BIC. 
%   compI:      (numModels x numComp) Set of indicator variables showing the presence / absence
%               of the components for each fitted (pcm_constructModelFamily)
% VARARGIN: 
%  prior:       Prior probability of each model
%               - A numModels x 1 vector of prior probabilities for the models 
%               - A numComponent x 1 vector or prior probabilities for each of the model components                
%               - A scalar is interpreted as the base-probability for each component to be present  
%               - If empty, it defaults to a flat prior over models 
% OUTPUT:      
%  postProb:    (numSubj x numComp) Posterior probability of the presence of all components for all subjects 
%  logBayes:    (numSubj x numComp) Log Bayes factor of the models assuming
%               presence of factor over the models asssuming absence 
% 2017 by Joern Diedrichsen & Eva Berlot 

prior = []; 
pcm_vararginoptions(varargin,{'prior'}); 

[numModels,numComp]=size(compI); 
if (size(likelihood,2)~=numModels)
    error('log-likelihood needs to be a numSubj x numModels matrix'); 
end; 

% Generate a log-prior for each MODEL
if (isempty(prior))
    logPrior = log(ones(1,numModels)./numModels);
elseif isscalar(prior) 
    logPrior = sum(compI.*log(prior)+(1-compI).*log(1-prior),2)'; 
elseif (numel(prior)==numComp) 
    logPrior = zeros(1,numModels); 
    for i=1:numComp 
        logPrior=logPrior+compI(:,i)'.*log(prior(i))+(1-compI(:,i))'.*log(1-prior(i)); 
    end; 
elseif (numel(prior)==numModels) 
    logPrior = log(prior); 
else 
    error('Prior must either have the length numModels, numComponents, 1 or 0'); 
end; 

logJoint = bsxfun(@plus,likelihood,logPrior);
logJoint = bsxfun(@minus,logJoint,max(logJoint,[],2))+50;   % Scale the log-joint, so it doesn't give errors with exp
logLike  = bsxfun(@minus,likelihood,max(likelihood,[],2))+50;

for i=1:numComp 
    indx = compI(:,i)==1; 
   % check if any of components all set to 1s/0s
   % calculate posterior probability of component
   % calculate logBayes factor in favor of component 
   if sum(indx)==length(indx) || sum(indx)==0
       logPosterior(1:length(logJoint),i)=logPrior(i);
       logBayes(1:length(logLike),i)=NaN;
   else
       logPosterior(:,i)=log(sum(exp(logJoint(:,indx)),2))-log(sum(exp(logJoint),2));
       logBayes(:,i)=log(sum(exp(logLike(:,indx)),2))-log(sum(exp(logLike(:,~indx)),2));
   end
end; 

postProb = exp(logPosterior); 