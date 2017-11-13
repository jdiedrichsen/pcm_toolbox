function [theta,l,k,reg]=pcm_NR(theta0,likefcn,varargin)
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
OPT.numIter = 80;                 % Maximal number of iterations
OPT.thres   = 1e-4;               % Tolerance on decrease of likelihood 
OPT.HessReg = 1;                  % Regularisation on the Hessian (Fisher) matrix 
OPT.verbose = 0; 
OPT.regularization = 'L';         % Type of regularisation 'L':Levenberg, 'LM': Levenberg-Marquardt 

% Variable argument otions
%--------------------------------------------------------------------------
try
OPT=rsa.getUserOptions(varargin,OPT,{'HessReg','thres','low','numIter','verbose','regularization'});
catch
    error('This might need the RSA toolbox to work (https://github.com/rsagroup/rsatoolbox). Do you have it in the matlab path?')
end
% Set warning to error, so it can be caught 
warning('error','MATLAB:nearlySingularMatrix');


% Initialize Interations 
%--------------------------------------------------------------------------
dF    = Inf;
H     = length(theta0); % Number of parameters 
% OPT.HessReg = OPT.HessReg*eye(H,H);          % Prior precision (1/variance) of h
theta=theta0; 
for k = 1:OPT.numIter
    thetaH(:,k)=theta;
    regH(k)=OPT.HessReg;
    [nl(k),dFdh,dFdhh]=likefcn(theta);
    
    % Safety check if negative likelihood decreased
    %----------------------------------------------------------------------
    if (k>1 & (nl(k)-nl(k-1))>eps)      % If not.... 
        OPT.HessReg = OPT.HessReg*10;   % Increase regularisation 
        if (OPT.verbose) 
            fprintf('likelihood decreased. Regularisation %2.3f\n',OPT.HessReg(1)); 
        end; 
        theta = thetaH(:,k-1); 
        thetaH(:,k)=theta;
        nl(k)=nl(k-1); 
        dFdh = dFdh_old; 
        dFdhh = dFdhh_old; 
        dL = inf;                       % Definitely try again 
    else
        if (OPT.HessReg>1e-8) 
            OPT.HessReg = OPT.HessReg/10;
        end; 
        dL = inf; 
        if (k>1) 
            dL = nl(k-1)-nl(k); 
        end;
    end; 
    dFdh_old=dFdh; 
    dFdhh_old = dFdhh; 
    
    % Add slight Levenberg-Marquardt-style regularisation to second derivative 
    %----------------------------------------------------------------------
    switch (OPT.regularization)
        case 'LM'       % Levenberg-Marquardt 
            dFdhh = dFdhh + diag(diag(dFdhh))*OPT.HessReg;
        case 'L'        % Levenberg 
            dFdhh = dFdhh + eye(size(dFdhh,1))*OPT.HessReg;
    end; 
    
    % Fisher scoring: update dh = inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    try 
        dtheta   =  dFdhh\dFdh;
        theta = theta - dtheta;
    catch % Likely matrix close to singluar 
        OPT.HessReg = OPT.HessReg*10;
        dFdhh = dFdhh_old + diag(diag(dFdhh_old))*OPT.HessReg;
        dtheta   =  dFdhh\dFdh;
        theta = theta - dtheta;
        keyboard; 
    end; 
    % convergence
    %----------------------------------------------------------------------
    if dL < OPT.thres
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
reg=OPT.HessReg; 
