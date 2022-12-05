function [G,Sig]=pcm_estGCrossval(B,partition,conditionVec,varargin)
% function [G,Sigma]=pcm_estGCrossval(Y,partition,conditionVec,varargin);
% Estimates the second moment matrix (G) using crossvalidation across partitions 
% (imaging runs). This yields an unbiased estimate of G- meaning that the
% expected value for the variances and covariances will be equal to the
% true variance and co-variances. 
% The estimated G will correspond exactly to the crossvalidated Mahalanobis 
% distance (see Walther et al., 2015):   
% C = pcm_indicatorMatrix('allpairs',[1:numCond]); 
% D = squareform(diag(C*G*C')); 
% 
% If the optional input argument X is given, then the data will be combined
% across partitions in an optimal way, taking into account the different
% variabilities of the estimators. 
% 
% If removeMean is set to 1, the routine will first remove the voxel mean
% in each partition (across conditions). If all conditions are present in
% all partitions, this will simply double-center the G matrix (all rows and
% columns will be zero). 
% 
% INPUT:
%  B           : Noise-normalized activation patterns, a N x P matrix
%                If X is provided for optimal weighting of the regressors,
%                it is important also to submit all regressors of
%                no-interest (i.e. intercepts, etc), such that their value
%                can be taken into account 
%  partition   : N x 1 integer value that indicates the partition for crossvalidation (typically run number)
%                These need to be between 1...M, zeros are being ignored 
%                if the Partition vector is shorter than N, it is assumed that the
%                last regressors are intercept for the different runs-as is
%                usual in SPM.
%  conditionVec: N x 1 vector of conditions, zeros will be ignored as
%                regressors of no interest. If conditionVec is shorter than 
%                N, it is assumed that all remaining numbers are 0. 
% VARARGIN: 
%  'X',X         : T x N Design matrix that is used to estimate from the first
%                level. This is an optional input parameter to optimally combine the beta weights 
%                across partitions. If temporal filtering / prewhitening is applied 
%                to the data, this needs to be the filtered design matrix. 
%                - If X is given, B needs to also contain the 
%                   Coefficients of no interest. 
%                - If X is not given, it assumes that the Beta-weights are i.i.d 
%  'removeMean',0: Removes the mean activitation for each voxel in each
%                partition. This Double-centers the G-matrix (i.e. makes
%                the mean of the rows and the mean of the columns 0). 
% OUTPUT: 
%   G          : Estimated second moment matrix
%   Sig        : a KxK covariance matrix of the beta estimates across
%               different imaging runs. 
% Joern Diedrichsen 
% 2/2016 

[N,numVox]          = size(B); 
part                = unique(partition)';
part(part==0)       = []; % Ignore the zero partitions 
if isvector(conditionVec)
    cond            = unique(conditionVec);
    cond(cond==0)   = [];
else
    cond            = 1:size(conditionVec,2);
end
numPart = numel(part);
numCond = numel(cond);%max(conditionVec); 

X=[]; 
removeMean = 0;  

pcm_vararginoptions(varargin,{'X','removeMean'}); 

if (removeMean & ~isempty(X))
    error('partition-wise mean subtraction with design matrix not implemented yet'); 
end; 


% Check on design matrix 
if (~isempty(X))  
    numReg     = size(X,2);             % Number of regressors in the first-level design matrix 
    if (numReg ~=N) 
        error('For optimal integration of beta weights, all N regressors (including no-interest) need to be submitted in Y'); 
    end; 
end; 

% Check length of partition vector  
missing = N-length(partition); 
if missing > 0 
    partition  = [partition;[1:missing]']; % Asssume that these are run intercepts 
end; 

% Check if condition vector is vector or matrix,
% then make second-level design matrix, pulling through the regressors of no-interest 
missing = N-length(conditionVec);
if ~isvector(conditionVec) % when design matrix has passed
    Z = conditionVec;
else
    if missing>0
        conditionVec = [conditionVec;zeros(missing,1)];
    end    
    Z = pcm_indicatorMatrix('identity_p',conditionVec);
end
% Deal with number of no-interest regressors
numNonInterest = sum(all(Z==0,2));      
Z(all(Z==0,2),end+[1:numNonInterest]) = eye(numNonInterest);
% numNonInterest = sum(conditionVec==0);      
% Z(conditionVec==0,end+[1:numNonInterest])=eye(numNonInterest);

% REmove the partition-wise mean (warning: may not work correctly with
% coniditions of no-interest in here). 
if (removeMean) 
   Zp=pcm_indicatorMatrix('identity_p',partition); 
   B=B-Zp*pinv(Zp)*B; % Remove run mean  
end;
    
% Preallocate the memory 
A = zeros(numCond,numVox,numPart);           % Allocate memory 
nA = zeros(numCond,numPart);                 % Condition by partition indicator
Bp = zeros(numCond,numVox);
nPart = zeros(numCond,numCond); 
% Estimate condition means within each run and crossvalidate 
for i=1:numPart 
    % Left-out partition 
    indxA = partition==part(i);
    Za    = Z(indxA,:); 
    Za    = Za(:,any(Za,1));       % restrict to regressors that are not all 0
    Ba    = B(indxA,:);            % Get regression coefficients 

    % remainder of conditions 
    indxB = partition~=part(i);
    Zb    = Z(indxB,:); 
    Zb    = Zb(:,any(Zb,1));    % Restrict to regressors that are not all 0 
    Bb    = B(indxB,:);
    
    % valid conditions of interest
    interestA = find(any(Z(indxA,1:numCond),1));
    interestB = find(any(Z(indxB,1:numCond),1));
    
    % Use design matrix if present to get GLS estimate 
   if (~isempty(X))
        Xa      = X(:,indxA);
        Xb      = X(:,indxB);
        indxX   = any(Xa,1);    % Restrict to regressors that are used in this partition
        Za      = Xa*Za; 
        Za      = Za(:,indxX); 
        Zb      = Xb*Zb; 
        Ba      = Xa(:,indxX)*Ba(indxX,:);
        Bb      = Xb*Bb; 
   end; 
    a     = pinv(Za)*Ba;
    b     = pinv(Zb)*Bb;
    %A(:,:,i) = a(1:numCond,:); 
    %G(:,:,i)= A(:,:,i)*b(1:numCond,:)'/numVox;      % Note that this is normalised to the number of voxels 
    A(interestA,:,i)    = a(1:length(interestA),:); 
    nA(interestA,i)   = 1; 
    Bp(interestB,:)     = b(1:length(interestB),:);
    G(:,:,i)= A(:,:,i)*Bp'/numVox;      % Note that this is normalised to the number of voxels 
    nPart=nPart+(G(:,:,i)~=0);          % Keep track of how many test partitions had this combination of conditions 
end; 
G=sum(G,3)./nPart;                      % Divide by the number of cross-products to get mean.

% If requested, also calculate the estimated variance-covariance 
% matrix from the residual across folds. 
if (nargout>1) 
    meanA = sum(A,3); 
    meanA = bsxfun(@rdivide,meanA,sum(nA,2)); 
    R=bsxfun(@minus,A,meanA);
    R=bsxfun(@times,R,permute(nA,[1 3 2])); 
    for i=1:numPart
        Sig(:,:,i)=R(:,:,i)*R(:,:,i)'/numVox;
        nSig(:,:,i)=nA(:,i)*nA(:,i)'; 
    end;
    Sig=sum(Sig,3)./(sum(nSig,3)-1);
end; 