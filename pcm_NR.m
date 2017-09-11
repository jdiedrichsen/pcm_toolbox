function [theta,l,k]=pcm_NR(theta0,likefcn,varargin)
% function [theta,l,k]=pcm_NR(theta0,likefcn,varargin)
% Newton-Raphson algorithm. 
%
% likefcn:          Function handle that returns the 
%                   a) Negative log-likelihood 
%                   b) First derivative of the negative log-likelihood 
%                   c) Expected second derivative of the negative log-likelhood
%
% VARARGIN:
%   'numIter'     : Maximal number of iterations
%   'theta0'      : Starting values for the parameters (Hx1 vector)
%   'TolL'         : Tolerance of the likelihood (l-l'), where l' is
%                    projected likelihood
%   'HessReg':      Regulariser on the Hessian to increase the stability of
%                   the fit (set to 1/256)
%
% OUTPUT:
%   theta : Variance coefficients (one column for each iteration)
%   l     : Log-likelihood of p(y|theta) for maximal theta 
%           This is Type II liklihood type II  the best estimates of theta, integrated over u
%
% See also: pcm_NR_diag, pcm_NR_comp
% v.1:
%
% Copyright 2017 Joern Diedrichsen, joern.diedrichsen@googlemail.com

% Defaults
%--------------------------------------------------------------------------
ac      = [];                      % Which terms do I want to include?
OPT.numIter = 32;                 % Maximal number of iterations
OPT.low     = -16;                % Low value of parameters to cut them out  
OPT.thres   = 1e-2;                     
OPT.HessReg = 1/256;                  % Regularisation on the Hessian (Fisher) matrix 

% Variable argument otions
%--------------------------------------------------------------------------
rsa.getUserOptions(varargin, {'theta0','HessReg','thres','low','numIter'});


% Initialize Interations 
%--------------------------------------------------------------------------
dF    = Inf;
H     = length(theta0); % Number of parameters 
OPT.HessReg = OPT.HessReg*speye(H,H);          % Prior precision (1/variance) of h
theta=theta0; 
as = true(H,1); 
for k = 1:OPT.numIter
    [nl(k),dFdh,dFdhh]=likefcn(theta);
    
    % Safety check if likelihood decreased
    %----------------------------------------------------------------------
    if (k>1 & nl(k)>nl(k-1))
        warning('likelihood decreased. Exiting...\n'); 
        break; 
    end; 
    
    % Add slight regularisation to second derivative 
    %----------------------------------------------------------------------
    dFdhh = dFdhh + OPT.HessReg;
    
    % Fisher scoring: update dh = inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    dtheta   =  dFdhh(as,as)\dFdh(as);
    theta(as) = theta(as) - dtheta;
    dF    = dFdh(as)'*dtheta;
    
    % convergence
    %----------------------------------------------------------------------
    if dF < OPT.thres
        break;
    else
        % as  = find(theta > low);
        % h(h<low)=low;
        % as  = as(:)';
    end
    thetaH(:,k)=theta;
end
% Return likelihood 
if(nargout>1)
    l = -likefcn(theta); 
end; 