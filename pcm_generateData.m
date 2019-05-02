function [Y,partVec,condVec] = pcm_generateData(Model,theta,varargin);
% function  [Y,part,conditions] =pcm_generateData(Model,theta,D,numSim,signal,noise,varargin);
% pcm_generateData: Simulate multivariate data using the generative model specified in Model
% Noise and Signal strength can be specified for each simulation and voxel
% separately. 
% INPUT: 
%   Model:   Model to generate data from 
%   theta:   numParams x 1 vector of parameters for Model 
% VARARGIN: 
%   'numPart',number of partitions  (default 8)
%   'numVox', number of independent voxels (default 50)
%   'numSim', number of simulations,all returned in cell array Y (default 1)
%   'signal', Signal variance: scalar, <numSim x 1>, <1xnumVox>, or <numSim x numVox> (default 0.1)
%   'noise', Noise  variance: scalar, <numSim x 1>, <1xnumVox>, or <numSim x numVox> (default 1)
%   'signalDist',fcnhnd:    Functionhandle to distribution function for signal (default normal)
%   'noiseDist',fcnhnd:     Functionhandle to distribution function for noise
%   'design',X:             - Design matrix (for encoding-style models) 
%                           - Condition vector (for RSA-style models) 
%                           Design matrix and Condition vector are assumed
%                           to be for 1 partition only. 
%                           If not specified - the function assumes a
%                           RSA-style model with G being numCond x numCond
% OUTPUT: 
%    Y:          Cell array{numSim} of data 
%    partVec:    Vector indicating the independent partitions 
%    condVec:    Vector of conditions for RSA-style model 
%                Design matrix for Encoding-style models
%  Note that the function uses an "exact" generation of the signal. Thus,
%  the true (noiseless) pattern consistent across partitions has exactly
%  the representational structure specified by the model 
% 2017 joern.diedrichsen@googlemail.com 
% ToDo: 
% - Options for adding temporal co-variance structure of the noise 
% - Options for adding spatial co-variance structuer of the noise 

% Defaults:
numPart = 8; 
numVox = 50; 
numSim = 1; 
signal = 0.1; 
noise  = 1; 
noiseDist = @(x) norminv(x,0,1);   % Standard normal inverse for Noise generation 
signalDist = @(x) norminv(x,0,1);  % Standard normal inverse for Signal generation 
design = [];
pcm_vararginoptions(varargin,{,'numPart','numVox','numSim','signal','noise','signalDist','noiseDist','design'}); 

% Make the overall generative model 
if (numel(theta)~=Model.numGparams)
    error('theta needs to be a numParams x 1 vector'); 
end; 
G = pcm_calculateG(Model,theta); 

if (isempty(design)) 
    numCond = size(G,1); 
    N          = numPart*numCond;          % Number of trials
    partVec    = kron([1:numPart]',ones(numCond,1));            % Partitions
    condVec    = kron(ones(numPart,1),[1:numCond]');            % Conditions
    Za         = kron(ones(numPart,1),eye(numCond));
else 
    if (size(design,2)==1) % RSA-style condition vector 
        ZZ        = pcm_indicatorMatrix('identity',design); 
        if (size(G,1)~=size(ZZ,2))
            error('For RSA-style models, the design needs to contain as many conditions as G'); 
        end; 
        condVec      = kron(ones(numPart,1),[1:ones(size(ZZ,2))]');            % Conditions
        numCond = size(ZZ,2);
        N       = numPart * size(ZZ,1); 
        partVec   =  kron([1:numPart]',ones(size(ZZ,1),1));
        Za        =  kron(ones(numPart,1),ZZ);
    else                    % Encoding-style design matrix  
        ZZ = design; 
        numCond = size(ZZ,2);
        if (size(G,1)~=size(ZZ,2))
            error('For Encoding-style models, the size(G) needs to be equal to the number of columns in design (feature)'); 
        end; 
        N       = numPart * size(ZZ,1); 
        partVec   =  kron([1:numPart]',ones(size(ZZ,1),1));
        Za        =  kron(ones(numPart,1),ZZ);
        condVec   =  Za; 
    end; 
end; 

for n = 1:numSim
    % Determine noise and signal for this simulation 
    if (size(signal,1)>1)
        thisSig = signal(n,:); 
    else 
        thisSig = signal; 
    end; 
    if (size(noise,1)>1)
        thisNoi = noise(n,:); 
    else 
        thisNoi = noise; 
    end; 
    
    % Generate true pattern from specified second moment 
    K = size(G,1); 
    pSignal = unifrnd(0,1,numCond,numVox); 
    U       = signalDist(pSignal); 
    E       = (U*U'); 
    Z       = E^(-0.5)*U;   % Make random orthonormal vectors 
    A       = pcm_diagonalize(G); 
    if (size(A,2)>numVox)
        error('not enough voxels to represent G'); 
    end; 
    trueU = A*Z(1:size(A,2),:)*sqrt(numVox); 
    trueU = bsxfun(@times,trueU,sqrt(thisSig));   % Multiply by (voxel-specific) signal scaling factor 
    
    % Now add the random noise 
    pNoise = unifrnd(0,1,N,numVox); 
    Noise  = bsxfun(@times,noiseDist(pNoise),sqrt(thisNoi)); 
    Y{n}  = Za*trueU + Noise;
end; 
    
