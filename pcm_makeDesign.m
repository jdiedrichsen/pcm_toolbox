function [condVec,partVec] = pcm_makeDesign(varargin);
% function [condVec,partVec] = pcm_makeDesign(varargin);
% VARARGIN: 
%   'numPart',number of partitions  (default 8)
% OUTPUT: 
%    condVec:    Vector of conditions for RSA-style model 
%    partVec:    Vector indicating the independent partitions 
numPart = 8; 
numCond = 5; 
N          = numPart*numCond;          % Number of trials
partVec    = kron([1:numPart]',ones(numCond,1));            % Partitions
condVec    = kron(ones(numPart,1),[1:numCond]');            % Conditions
