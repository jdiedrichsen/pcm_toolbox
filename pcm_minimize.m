function [G,theta,u,l,k]=pcm_minimize(Y,Z,M,varargin)
% pcm_minimize: estimate random-effects variance component coefficients using
% minimsation of cost function through line search (minimize).
%
% y_n = X b_n + Z u_n + e,
%         u ~ (b, G)
%           G = sum (theta * Gc)
%
% y: N x P observations
% Z: N x K random effects matrix
% X: N x Q fixed effects matrix
%
% N: numbers of observations (trials)
% P: numbers of experimental units (voxels)
% Q: number of random effects
% 
%   'M'            : Model structure to fit  
%
% VARARGIN:
%   'MaxIteration' : Maximal number of iterations
%   'theta0'           : Starting values for the parameters (Hx1 vector)
%                    (otherwise uses starting guess of zero)
%   'meanS'        : Remove the mean for each pattern component (a)
%                    (Logical flag, true by default)
%   'X'            : Fixed effects matrix that will be removed from y
%                    In this case, ReML will be used for the estimation of
%                    G.
%                    Important: X will also be removed from Z to orthogonilise the
%                    random effects from the constant effects
%   'gradcheck'    : Optional checking of gradient after maximization 
%
% OUTPUT:
%   G     : variance-covariance matrix
%   theta : Variance coefficients (one column for each iteration)
%   u     : hidden patterns components
%   l     : log-likelihood p(y|theta) for maximal theta, i.e. the maximal 
%               liklihood type II  for the best estimates of theta, integrated over u
%
% Examples:
%
% v.1.1:
%
% Copyright 2017 Joern Diedrichsen, joern.diedrichsen@googlemail.com 

% Defaults
%--------------------------------------------------------------------------
X            = [];                 % By default fixed effects are empty 
maxIteration = 1000; 
theta0           = []; 
gradcheck    = 0;  
% Variable argument otions
%--------------------------------------------------------------------------
vararginoptions(varargin, ...
    {'Gc','X','runEffect','maxIteration','theta0','gradcheck'});

% check input size
%--------------------------------------------------------------------------
[N,P]  =  size(Y);
[N2,K] =  size(Z);
if N2  ~= N
    error('Mismatched numbers of rows in data (%d) and design (%d)', N, N2)
end


% Minimize 
%--------------------------------------------------------------------------
fcn = @(x) pcm_likelihood(x,Y*Y',M,Z,X,P);
if (isempty(h0))
    theta0                  = zeros(M.numGparams+1,1);
end; 
[theta,fx]          = minimize(h0, fcn, MaxIteration);

% Check gradient if necessary 
%--------------------------------------------------------------------------
if (gradcheck) 
    checkderiv(fcn,h+0.02,0.0001); 
end;

% Provide other output functions 
G = pcm_calculateG(M,theta(1:M.numGparams));
u=[];
l=-fx(end); 