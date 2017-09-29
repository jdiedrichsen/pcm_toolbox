function [Y,part,conditions] = pcm_generateData(Model,theta,D,numSim,scale,noise,varargin);
% function  [Y,part,conditions] =pcm_generateData(Model,theta,D,numSim,scale,noise);
% pcm_generateData: Simulate multivariate data using the generative model
% specified in M, using the parameters theta and the experimental
% paramaters as specified in D.
%
% Experimental structure D: 
%   D.numPart = number of partititions 
%   D.numVox  = number of independent voxels 

noiseDist = @(x) norminv(x,0,1);   % Standard normal inverse for Noise generation 
signalDist = @(x) norminv(x,0,1);  % Standard normal inverse for Signal generation 

pcm_vararginoptions(varargin,{'noiseDist','signalDist'}); 

% Make the overall generative model 
G = pcm_calculateG(Model,theta); 
D.numCond = size(G,1); 

D.N             = D.numPart*D.numCond;          % Number of trials
part            = kron([1:D.numPart]',ones(D.numCond,1));            % Partitions
conditions      = kron(ones(D.numPart,1),[1:D.numCond]');            % Conditions

% Get design matrix
Za  = kron(ones(D.numPart,1),eye(D.numCond));
if (length(scale)==1); 
    scale=ones(numSim,1)*scale; 
end; 
if (length(noise)==1); 
    noise=ones(numSim,1)*noise; 
end; 
for n = 1:numSim
        
    % Generate true pattern from predicted distances
    K = size(G,1); 
    pSignal = unifrnd(0,1,D.numCond,D.numVox); 
    U       = signalDist(pSignal); 
    E       = (U*U'); 
    Z       = E^(-0.5)*U;   % Make random orthonormal vectors 
    A       = pcm_diagonalize(G*scale(n)); 
    if (size(A,2)>D.numVox)
        error('not enough voxels to represent G'); 
    end; 
    trueU = A*Z(1:size(A,2),:)*sqrt(D.numVox); 
    pNoise = unifrnd(0,1,D.N,D.numVox); 
    Noise  = sqrt(noise(n)) * noiseDist(pNoise); 
    Y{n}  = Za*trueU + Noise;
end; 
    
