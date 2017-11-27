function [knockIN,knockOUT]=pcm_knockModels(likelihood,compI,varargin); 
% function [knockIn,knockOut]=pcm_knockModels(likelihood,compI,varargin); 
% Calculates importance of each feature by computing two values
%   knockIN likelihood - improvement in model performance
%                        if the feature is added to the NULL model
%   knockOUT likelihood - decrement in model performance
%                        if the feature is taken away from the FULL model
% Used for assessing the importance of each model feature
% INPUT: 
%   likelihood: (numSubj x numModels) estimated marginal log-likelihood of each model. This can be 
%               crossvalidated Type-II maximum likelihood or Type-II maximum likelihood 
%               corrected with BIC. 
%   compI:      (numModels x numFeature) Set of indicator variables showing the presence / absence
%               of the components for each fitted (pcm_constructModelFamily)
% OUTPUT:  
%   knockIN:    (numSubj x numFeatures) - knockIn likelihood values 
%   knockOUT:   (numSubj x numFeatures)

 
numComp=0;
indComp=[];
for c=1:size(compI,2)
    % check if feature always / never included in models
    if sum(compI(:,c))==length(compI) || sum(compI(:,c))==0
    else
        numComp=numComp+1;
        indComp=[indComp c];
   end
end
% specify knock-in and knock-out indicator matrices
KnockIn = indicatorMatrix('identity_p',[1:numComp]); 
KnockOut = 1 - indicatorMatrix('identity_p',[1:numComp]); 
% determine which are the 'null' and 'full' models
nullInd = find(sum(compI(:,indComp),2)==0);
fullInd = find(sum(compI(:,indComp),2)==numComp);

for i=1:numComp
    indx_knockIN=find(ismember(compI(:,indComp),KnockIn(i,:),'rows'));
    indx_knockOUT=find(ismember(compI(:,indComp),KnockOut(i,:),'rows'));
    knockIN(:,indComp(i))=likelihood(:,indx_knockIN)-likelihood(:,nullInd);
    knockOUT(:,indComp(i))=likelihood(:,indx_knockOUT)-likelihood(:,fullInd);
    %knockIN(:,i)=likelihood(:,indx_knockIN)-likelihood(:,nullInd);
    %knockOUT(:,i)=likelihood(:,indx_knockOUT)-likelihood(:,fullInd);
end

% check if dimensions are correct 
% used in cases the last feature is always / never included (i.e. not knocked-IN/OUT)
if size(knockIN,2)==size(compI,2)
else
    knockIN(:,size(compI,2))=zeros(length(knockIN),1);
    knockOUT(:,size(compI,2))=zeros(length(knockOUT),1);
end

 