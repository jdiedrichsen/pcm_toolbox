function [M,Z] = pcm_buildModelFromFeatures(featuresets,varargin); 
% function [M,Z] = pcm_buildModelFromFeatures(features,varargin); 
% Builds a PCM model from a set of model features
% INPUT: 
%   featuresets:  cell array of featuresets, each containing a NxQ_i matrix 
% VARARGIN: 
%   'style': Specific style in formulating the PCM-model 
%          - 'encoding_style': Puts the representational structure into
%                              design matrix Z, components in M are
%                              diagnonal matrices (i.e. ridge regression)
%          - 'rsa_style':      Puts the representational structure into the
%                              component matrices  of M, the design matrix
%                              Z becomes a simple indicator matrix for the
%                              condition number 
%   'name': Specifies name of model 
%   'type': Specifies output type of model. Choice are 
%           'component'     : Model components are written as Gc/Gd 
%           'feature'       : Model components are written as Ac (Gc = Ac*Ac') 
%           'componentSq'   : Model components are written as Gd, but using theta^2 instead of exp(theta)  
% OUTPUT: 
%   M: Model structure of the desired type 
%   
name  = [];
style = 'encoding_style';
type  = 'component'; 

pcm_vararginoptions(varargin,{'style','name','type'}); 


if (~iscell(featuresets))
    featuresets={featuresets};
end; 

numFeatSets = numel(featuresets);    % Number of feature sets, indendepently weighted 
numCond = size(featuresets{1},1);   % number of conditions

% total number of columns across all condition vectors / matrices
numFeat=0;
for f=1:numFeatSets
    if (size(featuresets{f},1)~=numCond)
        error('all feature sets need to have the same number of rows'); 
    end; 
    if (size(featuresets{f},2)>numCond)
        warning('Feature set %d is overspecified- i.e. more features than observations\nFor faster convergence it is recommended to reduced feature set with SVD',f); 
    end; 
    numFeat = numFeat + size(featuresets{f},2); 
end

% Start buidling model structure 
M.type = type;
M.name = name; 

% Depending on the style, build the model into design matrix or second moment 
switch (style) 
    case 'encoding_style' 
        Z    = zeros(numCond,numFeat);      % Design matrix 
        Gd   = zeros(numFeatSets,numFeat); % Indicator matrix for columns in Design matrix  
        Gc    = zeros(numFeat,numFeat,numFeatSets); 
        colCount = 1; 
        for f=1:numFeatSets
            Q = size(featuresets{f},2); 
            Z(:,colCount:colCount+Q-1) = featuresets{f};
            Gd(f,colCount:colCount+Q-1) = 1;
            Gc(colCount:colCount+Q-1,colCount:colCount+Q-1,f)=eye(Q); 
            colCount = colCount+Q; 
        end; 
        switch(type) 
            case {'component','componentSq'}
                M.Gd  = Gd; 
                M.Gc  = Gc; 
            case 'feature' 
                M.Gd  = Gd; 
                M.Ac  = Gc;
        end; 
    case 'rsa_style' 
        X    = zeros(numCond,numFeat); % Model feature matrix 
        colCount = 1; 
        for f=1:numFeatSets
            Q = size(featuresets{f},2); 
            X(:,colCount:colCount+Q-1) = featuresets{f};
            colCount = colCount+Q; 
        end; 
        [X,~,condVec]=unique(X,'rows','stable'); % Reduce the model by finding unique 'conditions' 
        numCond=size(X,1); 
        switch(type) 
            case {'component','componentSq'}
                M.Gc   = zeros(numCond,numCond,numFeatSets);
                colCount = 1; 
                for f=1:numFeatSets 
                    Q = size(featuresets{f},2); 
                    M.Gc(:,:,f) = X(:,colCount:colCount+Q-1)*X(:,colCount:colCount+Q-1)';
                    colCount = colCount+Q; 
                end; 
            case 'feature' 
                M.Ac   = zeros(numCond,numFeat,numFeatSets);
                colCount = 1; 
                for f=1:numFeatSets 
                    Q = size(featuresets{f},2); 
                    M.Ac(:,:,f) = X(:,colCount:colCount+Q-1);
                    colCount = colCount+Q; 
                end; 
        end; 
        Z = condVec;       % Deisgn matrix is the condition vector 
end; 
