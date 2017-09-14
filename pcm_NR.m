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
OPT.numIter = 32;                 % Maximal number of iterations
OPT.thres   = 1e-5;                     
OPT.HessReg = 1/256;                  % Regularisation on the Hessian (Fisher) matrix 

% Variable argument otions
%--------------------------------------------------------------------------
OPT=rsa.getUserOptions(varargin,OPT,{'HessReg','thres','low','numIter'});


% Initialize Interations 
%--------------------------------------------------------------------------
dF    = Inf;
H     = length(theta0); % Number of parameters 
OPT.HessReg = OPT.HessReg*eye(H,H);          % Prior precision (1/variance) of h
theta=theta0; 
for k = 1:OPT.numIter
    thetaH(:,k)=theta;
    [nl(k),dFdh,dFdhh]=likefcn(theta);
    
    % Safety check if likelihood decreased
    %----------------------------------------------------------------------
    if (k>1 & (nl(k)-nl(k-1))>eps)
        OPT.HessReg = OPT.HessReg*10; 
        % warning('likelihood decreased. Regularisation %2.3f\n',OPT.HessReg(1)); 
        theta = thetaH(:,k-1); 
        thetaH(:,k)=theta;
        nl(k)=nl(k-1); 
        dFdh = dFdh_old; 
        dFdhh = dFdhh_old; 
    else 
        OPT.HessReg = OPT.HessReg/10; 
    end; 
    dFdh_old=dFdh; 
    dFdhh_old = dFdhh; 
    
    % Add slight regularisation to second derivative 
    %----------------------------------------------------------------------
    dFdhh = dFdhh + OPT.HessReg;
    
    % Fisher scoring: update dh = inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    dtheta   =  dFdhh\dFdh;
    % dtheta   =  -spm_dx(-dFdhh,-dFdh,{1});
    theta = theta - dtheta;
    dF    = dFdh'*dtheta;
    
    % convergence
    %----------------------------------------------------------------------
    if dF < OPT.thres
        break;
    else
        % as  = find(theta > low);
        % h(h<low)=low;
        % as  = as(:)';
    end
    
end
% Return likelihood 
if(nargout>1)
    l = -likefcn(theta); 
end; 