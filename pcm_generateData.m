function [Y,partVec,condVec] = pcm_generateData(Model,theta,D,numSim,signal,noise,varargin);
% function  [Y,part,conditions] =pcm_generateData(Model,theta,D,numSim,signal,noise,varargin);
% pcm_generateData: Simulate multivariate data using the generative model specified in Model
% Noise and Signal strength can be specified for each simulation and voxel
% separately. 
% INPUT: 
%   theta:   numParams x 1 vector of parameters for Model 
%   D:       Experimental structure with fields 
%       D.numPart = number of partititions 
%       D.numVox  = number of independent voxels
%   numSim:  number of simulations,all returned in cell array Y
%   signal:  Signal variance: scalar, <numSim x 1>, <1xnumVox>, or <numSim x numVox>
%   noise:   Noise  variance: scalar, <numSim x 1>, <1xnumVox>, or <numSim x numVox>
% VARARGIN: 
%   'signalDist',fcnhnd:    Functionhandle to distribution function for signal
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

noiseDist = @(x) norminv(x,0,1);   % Standard normal inverse for Noise generation 
signalDist = @(x) norminv(x,0,1);  % Standard normal inverse for Signal generation 
design = []; 
pcm_vararginoptions(varargin,{'noiseDist','signalDist','design'}); 

% Make the overall generative model 
if (size(theta,1)~=Model.numGparams || size(theta,2)~=1)
    error('theta needs to be a numParams x 1 vector'); 
end; 
G = pcm_calculateG(Model,theta); 

if (isempty(design)) 
    D.numCond = size(G,1); 
    D.N          = D.numPart*D.numCond;          % Number of trials
    partVec      = kron([1:D.numPart]',ones(D.numCond,1));            % Partitions
    condVec      = kron(ones(D.numPart,1),[1:D.numCond]');            % Conditions
    Za            = kron(ones(D.numPart,1),eye(D.numCond));
else 
    if (size(design,2)==1) % RSA-style condition vector 
        ZZ        = pcm_indicatorMatrix('identity',design); 
        if (size(G,1)~=size(ZZ,2))
            error('For RSA-style models, the design needs to contain as many conditions as G'); 
        end; 
        condVec      = kron(ones(D.numPart,1),[1:ones(size(ZZ,2))]');            % Conditions
        D.numCond = size(ZZ,2);
        D.N       = D.numPart * size(ZZ,1); 
        partVec   =  kron([1:D.numPart]',ones(size(ZZ,1),1));
        Za        =  kron(ones(D.numPart,1),ZZ);
    else                    % Encoding-style design matrix  
        ZZ = design; 
        D.numCond = size(ZZ,2);
        if (size(G,1)~=size(ZZ,2))
            error('For Encoding-style models, the size(G) needs to be equal to the number of columns in design (feature)'); 
        end; 
        D.N       = D.numPart * size(ZZ,1); 
        partVec   =  kron([1:D.numPart]',ones(size(ZZ,1),1));
        Za        =  kron(ones(D.numPart,1),ZZ);
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
    pSignal = unifrnd(0,1,D.numCond,D.numVox); 
    U       = signalDist(pSignal); 
    E       = (U*U'); 
    Z       = E^(-0.5)*U;   % Make random orthonormal vectors 
    A       = pcm_diagonalize(G); 
    if (size(A,2)>D.numVox)
        error('not enough voxels to represent G'); 
    end; 
    trueU = A*Z(1:size(A,2),:)*sqrt(D.numVox); 
    trueU = bsxfun(@times,trueU,sqrt(thisSig));   % Multiply by (voxel-specific) signal scaling factor 
    
    % Now add the random noise 
    pNoise = unifrnd(0,1,D.N,D.numVox); 
    Noise  = bsxfun(@times,noiseDist(pNoise),sqrt(thisNoi)); 
    Y{n}  = Za*trueU + Noise;
end; 
    
