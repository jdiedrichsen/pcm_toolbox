function Y = pcm_makeDataset(Model,theta,varargin);
% function  Y =pcm_makeDataset(Model,theta,varargin);
% pcm_generateData: Simulate multivariate data using the generative model specified in Model
% Noise and Signal strength can be specified for each simulation and voxel
% separately. 
% INPUT: 
%   Model:   Model to generate data from 
%   theta:   numParams x 1 vector of parameters for Model 
% VARARGIN: 
%   'numVox', number of independent voxels (default 50)
%   'numSim', number of simulations,all returned in cell array Y (default 1)
%   'signal':               - Signal variance: scalar, <numSim x 1>      (default 0.1)
%   'signalSpatialCov':     - voxel x voxel covariance matrix of signal (default identity)
%   'noise':                - Noise  variance: scalar, <numSim x 1>     (default 1)
%   'noiseTrialCov':        - Trial x trial covariance matrix of noise  (default identity)
%   'noiseSpatialCov':      - voxel x voxel covariance matrix of noise  (default identity)
%   'signalDist',fcnhnd:    - Functionhandle to distribution function for signal (default normal)
%   'noiseDist',fcnhnd:     - Functionhandle to distribution function for noise (default normal)
%   'design',X:             - Design matrix (for encoding-style models) 
%                           - Condition vector (for RSA-style models) 
%   'samesignal',false: Should we use exactly the same pattern for the
%                       signal for each simulation? 
%   'exact',true: Make the signal with exact second moment matrix G. If set
%                 to false, the signal will be simply drawn from a
%                 multivariate normal with variance G
% OUTPUT: 
%    Y:          Cell array{numSim} of data 
%  Note that the function uses an "exact" generation of the signal. Thus,
%  the true (noiseless) pattern consistent across partitions has exactly
%  the representational structure specified by the model 
% 2017 joern.diedrichsen@googlemail.com 
% Defaults:
numVox = 50; 
numSim = 1; 
signal = 0.1; 
signalSpatialCov =1; % Instead of Identity, we can use a scalar 1 
noise  = 1; 
noiseTrialCov = 1;   % Instead of Identity, we can use a scalar 1 
noiseSpatialCov = 1; % Instead of Identity, we can use a scalar 1 
noiseDist = @(x) norminv(x,0,1);   % Standard normal inverse for Noise generation 
signalDist = @(x) norminv(x,0,1);  % Standard normal inverse for Signal generation 
noiseTrial = 1; % shortcut for identity if no trial noise structure given
design = [];
samesignal = false; 
exact = true;                       % Make signal to have exactly G covariance (for the given spatial structure)
pcm_vararginoptions(varargin,{'numVox','numSim','signal','signalSpatialCov','noise','noiseTrialCov','noiseSpatialCov',...
                        'signalDist','noiseDist','noiseTrial','design','samesignal','exact'}); 

% Make the overall generative model 
if (size(theta,1)~=Model.numGparams)
    error('theta needs to be a numParams x 1 vector'); 
end; 
G = pcm_calculateG(Model,theta); 


if (isempty(design)) 
    error('Need to specify design'); 
else 
    N  = size(design,1); 
    if (size(design,2)==1) % RSA-style condition vector 
        Za        = pcm_indicatorMatrix('identity',design); 
        if (size(G,1)~=size(Za,2))
            error('For RSA-style models, the design needs to contain as many conditions as G'); 
        end; 
    else                    % Encoding-style design matrix  
        Za = design; 
        if (size(G,1)~=size(Za,2))
            error('For Encoding-style models, the size(G) needs to be equal to the number of columns in design (feature)'); 
        end; 
    end; 
end; 

% precompute spatial and temporal covariance matrix decompositions 
signalSpatialChol = cholcov(signalSpatialCov);     
noiseSpatialChol = cholcov(noiseSpatialCov);     
noiseTrialChol = cholcov(noiseTrialCov);     

for n = 1:numSim
    % Determine signal for this simulation 
    if (numel(signal) == numSim)
        thisSig = signal(n); 
    elseif (numel(signal)==1) 
        thisSig = signal; 
    else 
        error('Signal needs to be either a scalar or vector of them numSim'); 
    end; 
        
    % Determine noise for this simulation     
    if (numel(noise)==numSim) 
        thisNoi = noise(n); 
    elseif (numel(noise) == 1) 
        thisNoi = noise; 
    else 
        error('Signal needs to be either a scalar or vector of them numSim'); 
    end; 
    
    % Generate true pattern from specified second moment 
    % samesignal = true: generate it on the first simulation and keep the same 
    % samesignal = false: generate new every time 
    if (n==1 || ~samesignal) 
        K = size(G,1); 
        pSignal = unifrnd(0,1,K,numVox); 
        U       = signalDist(pSignal)*signalSpatialChol; 
        
        % If exact = true - make a matrix of random numbers with exactly
        % the correct covariance matrix 
        A       = pcm_diagonalize(G); % A*A' = G 
        if (exact) 
            E       = (U*U')/numVox; 
            U       = E^(-0.5)*U;   % Make random orthonormal vectors 
            if (size(A,2)>numVox)
                error('not enough voxels to represent G'); 
            end; 
        end; 
        trueU = A*U(1:size(A,2),:)*sqrt(thisSig); 
    end; 
    
    % Now add the random noise 
    pNoise = unifrnd(0,1,N,numVox); 
    Noise  = noiseTrialChol*noiseDist(pNoise)*noiseSpatialChol*sqrt(thisNoi); 
    Y{n}  = Za*trueU + Noise;
end;