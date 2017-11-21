function [postProp,logBayes]=pcm_componentPosterior(likelihood,compI,varargin); 
% function p=pcm_componentPosterior(likelihood,compI); 
% Calculates a posterior proability for the presence of each 
% component from a family of models 
% INPUT: 
%   likelihood: (numSubj x numModels) estimated marginal log-likelihood of each model. This can be 
%               crossvalidated Type-II maximum likelihood or Type-II maximum likelihood 
%               corrected with BIC. 
%   compI:      (numModels x numComp) Set of indicator variables showing the presence / absence
%               of the components for each fitted (pcm_constructModelFamily)
% VARARGIN: 
%  prior:       Prior probability of each model (assumed to be flat over
%               models) 
% OUTPUT:  
%               
prior = []; 
pcm_vararginoptions(varargin,{'prior'}); 

[numModels,numComp]=size(compI); 

if (isempty(prior))
    prior = ones(1,numModels)./numModels; 
end; 

nCondLogL=bsxfun(@minus,likelihood,mean(likelihood,2)); % Conditional log-probability up to a constant scaling factor 
condProp  = exp(nCondLogL); 
jointProp = bsxfun(@times,condProp,prior); 

for i=1:numComp 
    indx = compI(:,i)==1; 
    postProp(:,i)=sum(jointProp(:,indx),2)./sum(jointProp,2); % Posterior proability of component 
    logBayes(:,i)=log(sum(jointProp(:,indx),2))-log(sum(jointProp(:,~indx),2)); % Log-Bayes factor in favor of component 
end; 