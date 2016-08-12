function [G,theta,u,l,k]=pcm_minimize(Y,Z,varargin)
% pcm_minimize: estimate random-effects variance component coefficients using
% minimsation of cost function through line search (minimize).
%
% y_n = X b_n + Z u_n + e,
%         u ~ (b, G)
%           G = sum (theta * Gc)
%
% y: N x P observations
% Z: N x K random effects matrix
% Z: N x Q fixed effects matrix
%
% N: numbers of observations (trials)
% P: numbers of experimental units (voxels)
% Q: number of random effects
%
% VARARGIN:
%   'num_iter'     : Maximal number of iterations
%   'Gc'           : Cell array (Hx1) of components of variance components 
%                    matrix G = sum(h_m Gc{m}), with m=1...H 
%   'h0'           : Starting values for the parameters (Hx1 vector)
%                    (otherwise uses starting guess based on Laird & Lange)
%   'TolL'         : Tolerance of the likelihood (l-l'), where l' is
%                    projected likelihood
%   'meanS'        : Remove the mean for each pattern component (a)
%                    (Logical flag, true by default)
%   'X'            : Fixed effects matrix that will be removed from y
%                    In this case, ReML will be used for the estimation of
%                    G.
%                    Important: X will also be removed from Z to orthogonilise the
%                    random effects from the constant effects, making the
%                    estimation of b_n independent of G
%   'gradcheck'    : Optional checking of gradient after maximization 
%
% OUTPUT:
%   G     : variance-covariance matrix
%   theta : Variance coefficients (one column for each iteration)
%   u     : hidden patterns components
%   l     : 1x2 vector of likelihoods 
%           l(1) = likelihood p(y|theta) for maximal theta, i.e. the maximal 
%               liklihood type II  for the best estimates of theta, integrated over u
%           l(2) = marginal liklihood p(y) based on the Laplace (Normal) 
%           approximation around the posterior mode of log(theta) 
%
% Examples:
% See mva_component_examples
%
% See also: mva_component_examples, spm_reml, spm_reml_sc, spm_sp_reml
% Where spm_* are from the SPM software, http://www.fil.ion.ucl.ac.uk/spm
%
% v.3.0:
%
% Copyright 2014 Joern Diedrichsen, j.diedrichsen@ucl.ac.uk

% Defaults
%--------------------------------------------------------------------------
X            = [];                 % By default fixed effects are empty 
MaxIteration = 1000; 
remove       = 0; 
h0           = []; 
gradcheck    = 0; 
% Variable argument otions
%--------------------------------------------------------------------------
vararginoptions(varargin, ...
    {'Gc','X','runEffect','MaxIteration','remove','h0','gradcheck'});

% check input size
%--------------------------------------------------------------------------
[N,P]  =  size(Y);
[N2,K] =  size(Z);
if N2  ~= N
    error('Mismatched numbers of rows in data (%d) and design (%d)', N, N2)
end

% Intialize the Model structure
%--------------------------------------------------------------------------
H     =  length(Gc)+1;       % Number of Hyperparameters (+ 1 for noise term)
for i =  1:H-1
    Gc{i} = sparse(Gc{i});
end;

% If fixed effects are given, remove fixed effects from data and random
% effects matrix
%--------------------------------------------------------------------------
if (~isempty(X) & remove)
    pX = pinv(X);
    Z  = Z-X*pX*Z;
    Y  = Y-X*pX*Y;
end;

% Minimize 
%--------------------------------------------------------------------------
fcn = @(x) pcm_likelihood(x,Y*Y',Gc,Z,X,P);
if (isempty(h0))
    h0                  = [0;0];
end; 
[h,fx]          = minimize(h0, fcn, MaxIteration);

% Check gradient if necessary 
%--------------------------------------------------------------------------
if (gradcheck) 
    checkderiv(fcn,h+0.02,0.0001); 
end;

% Reassemble G matrix 
%--------------------------------------------------------------------------
G     = sparse(K,K);
for i = 1:H-1;
    G = G + Gc{i}*exp(h(i));
end
G=full(G);

% Provide other output functions 
u=[];
l=-fx(end); 
theta=exp(h); 