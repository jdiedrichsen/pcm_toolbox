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
%               - A numModelsx1 vector  
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

if (isempty(prior))
    logPrior = log(ones(1,numModels)./numModels);
elseif isscalar(prior) 
    logPrior = sum(compI.*log(prior)+(1-compI).*log(1-prior),2)'; 
end; 

logJoint = bsxfun(@plus,likelihood,logPrior);
logJoint = bsxfun(@minus,logJoint,max(logJoint,[],2))+50;   % Scale the log-joint, so it doesn't give errors with exp

for i=1:numComp 
    indx = compI(:,i)==1; 
   % check if any of components all set to 1s/0s
   % calculate posterior probability of component
   % calculate logBayes factor in favor of component 
   if sum(indx)==length(indx) || sum(indx)==0
       logPosterior(1:length(logJoint),i)=logPrior(i);
       logBayes(1:length(logJoint),i)=NaN;
   else
       logPosterior(:,i)=log(sum(exp(logJoint(:,indx)),2))-log(sum(exp(logJoint),2));
       logBayes(:,i)=log(sum(exp(logJoint(:,indx)),2))-log(sum(exp(logJoint(:,~indx)),2));
   end
end; 

postProb = exp(logPosterior); 