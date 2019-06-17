function [theta,l,k,reg,regH,thetaH]=pcm_NR(theta0,likefcn,varargin)
% function [theta,l,k]=pcm_NR(theta0,likefcn,varargin)
% Newton-Raphson algorithm.
% INPUT:
%   theta0            Vector of parameter starting values
%   likefcn:          Function handle that returns the
%                       a) Negative log-likelihood
%                       b) First derivative of the negative log-likelihood
%                       c) Expected second derivative of the negative log-likelhood
%
% VARARGIN:
%   'numIter'     : Maximal number of iterations
%   'thres'       : Stopping criterion on likelihood change
%   'HessReg':      Starting regulariser on the Hessian matrix (default starts at 1
%                   then being adjusted in decibels)
%   'regularization': 'L': Levenberg 'LM': Levenberg-Marquardt  'sEig':smallest Eigenvalue 
%   'verbose':     0: No feedback, 1:Important warnings, 2:full feedback regularisation
%
% OUTPUT:
%   theta : Variance coefficients
%   l     : Log-likelihood of p(y|theta) for maximal theta
%           This is a Type II maximal likelihood - maximal likelhood of theta, integrated over u
%   k     : Number of iterations
%   reg   : Final regularisation value
% See also: pcm_NR_diag, pcm_NR_comp
% v.1:
%
% Copyright 2017 Joern Diedrichsen, joern.diedrichsen@googlemail.com

% Defaults
%--------------------------------------------------------------------------
OPT.maxIteration = 80;                 % Maximal number of iterations
OPT.likeThres    = 1e-4;               % Tolerance on decrease of likelihood
OPT.HessReg      = 1;                  % Regularisation on the Hessian (Fisher) matrix
OPT.verbose      = 0;
OPT.regularization = 'sEig';         % Type of regularisation 'L':Levenberg, 'LM': Levenberg-Marquardt

% Variable argument otions
%--------------------------------------------------------------------------
OPT=pcm_getUserOptions(varargin,OPT,{'HessReg','likeThres','maxIteration','verbose','regularization'});

% Set warning to error, so it can be caught
CATCHEXP = {'MATLAB:nearlySingularMatrix','MATLAB:singularMatrix',...
    'MATLAB:illConditionedMatrix','MATLAB:posdef',...
    'MATLAB:nearlySingularMatrix','MATLAB:eig:matrixWithNaNInf'};
for i=1:numel(CATCHEXP)
    warning('error',CATCHEXP{i});
end;

% Initialize Interations
%--------------------------------------------------------------------------
dF    = Inf;
H     = length(theta0); % Number of parameters
% OPT.HessReg = OPT.HessReg*eye(H,H);          % Prior precision (1/variance) of h
theta=theta0;
for k = 1:OPT.maxIteration
    
    % If more than the first interation: Try to update theta
    if k>1        
        % Fisher scoring: update dh = inv(ddF/dhh)*dF/dh
        % if it fails increase regularisation until dFdhh is invertible
        % Add regularisation to second derivative
        %----------------------------------------------------------------------
        while true
            try
                switch (OPT.regularization)
                    case 'LM'       % Levenberg-Marquardt
                        H = dFdhh + diag(diag(dFdhh))*OPT.HessReg;
                        dtheta   =  H\dFdh;
                    case 'L'        % Levenberg
                        H = dFdhh + eye(size(dFdhh,1))*OPT.HessReg;
                        dtheta   =  H\dFdh;
                    case 'sEig' 
                        [VH,lH]=eig(dFdhh);
                        lH = diag(lH); 
                        lH(lH<OPT.HessReg)=OPT.HessReg;                 % Increase smallest eigenvalue 
                        dtheta  =  bsxfun(@times,VH,1./lH')*VH'*dFdh; 
                    end;
                    break;    
            catch ME% Likely matrix close to singluar
                if any(strcmp(ME.identifier,CATCHEXP))
                    if (OPT.verbose==2) 
                        fprintf('Ill-conditioned Hessian. Regularisation %2.3f\n',OPT.HessReg);
                    end;
                    OPT.HessReg = OPT.HessReg*10;
                    if (OPT.HessReg>100000)
                        if (OPT.verbose)
                            warning('Cant regularise second derivative.. Giving up\n');
                        end;
                        exitflag=3;  % Regularisation increased too much 
                        break;  % Give up
                    end;
                else
                    ME.rethrow;
                end;
            end;
        end;
        theta = theta - dtheta; 
    end;
    
    thetaH(:,k)=theta;
    regH(k)=OPT.HessReg;
    ME=[];
    try
        [nl(k),dFdh,dFdhh]=likefcn(theta);
    catch ME  % Catch errors based on invalid parameter settings
        if any(strcmp(ME.identifier,CATCHEXP))
            if (k==1)
                error('bad initial values for theta');
            else
                nl(k)=inf;         % Set new likelihood to -inf: take a step back
            end;
        else
            ME.rethrow;
        end;
    end;
    % Safety check if negative likelihood decreased
    %----------------------------------------------------------------------
    if (k>1 & (nl(k)-nl(k-1))>eps)      % If not....
        OPT.HessReg = OPT.HessReg*10;   % Increase regularisation
        if (OPT.verbose==2)
            fprintf('Step back. Regularisation %2.3f\n',OPT.HessReg);
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
    
    % convergence
    %----------------------------------------------------------------------
    if dL < OPT.likeThres
        break;
    end;
end
% Return likelihood
if(nargout>1)
    l = -likefcn(theta);
end;
reg=OPT.HessReg;