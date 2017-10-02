function [G,theta,u,l,k]=pcm_NR_comp(Y,Z,varargin)
% function [G,theta,u,l,k]=pcm_NR_comp(Y,Z,varargin)
% Estimate random-effects variance component coefficients using
% Newton-Raphson gradient descent.
% Diedrichsen, Ridgway, Friston & Wiestler (2011).
%
% y_n = X b_n + Z u_n + e,
%         u ~ (b, G)
%           G = sum (theta * Gc)
%
% Y: N x P observations
% Z: N x Q random effects matrix
%
% N: numbers of observations (trials)
% P: numbers of experimental units (voxels)
% Q: number of random effects
%
% VARARGIN:
%   'numIter'     : Maximal number of iterations
%   'Gc'           : QxQxH tensor of variance components
%                    matrix G = sum(h_m Gc(:,:,m)), with m=1...H
%   'h0'           : Starting values for the parameters (Hx1 vector)
%                    (otherwise uses starting guess based on Laird & Lange)
%   'TolL'         : Tolerance of the likelihood (l-l'), where l' is
%                    projected likelihood
%   'meanS'        : Remove the mean for each pattern component (a)
%                    (Logical flag, true by default)
%   'X'            : Fixed effects matrix that will be removed from y
%                    In this case, ReML will be used for the estimation of
%                    G.
%   'S':            Explicit noise covariance matrix structure matrix. The For speed, 
%                   this is a cell array that contains
%                   S.S:     Structure of noise 
%                   S.invS:  inverse of the noise covariance matrix 
%   'HessReg':      Regulariser on the Hessian to increase the stability of
%                   the fit (set to 1/256)
%
%
%
% OUTPUT:
%   G     : variance-covariance matrix
%   theta : Variance coefficients (one column for each iteration)
%   u     : hidden patterns components
%   l     : Log-likelihood of p(y|theta) for maximal theta 
%           This is Type II liklihood type II  the best estimates of theta, integrated over u
%
%
% See also: pcm_NR_diag, pcm_NR_comp
% v.1.1:
%
% Copyright 2017 Joern Diedrichsen, joern.diedrichsen@googlemail.com

% Defaults
%--------------------------------------------------------------------------
meanS   = 0;                    % Mean subtract
ac      = [];                      % Which terms do I want to include?
numIter = 32;                 % Maximal number of iterations
Gc      = {};
h0      = [];
low     = -16;
thres   = 1e-2;
HessReg = 1/256;                  % Precision on hyper prior
X       = [];                 % By default fixed effects are empty
S       = []; 

% Variable argument otions
%--------------------------------------------------------------------------
vararginoptions(varargin, ...
    {'Gc','meanS','h0','ac','HessReg','thres','X','numIter','S'});

% check input size
%--------------------------------------------------------------------------
[N,P]  =  size(Y);
[N2,Q] =  size(Z);
if N2  ~= N
    error('Mismatched numbers of rows in data (%d) and design (%d)', N, N2)
end

% Intialize the Model structure
%--------------------------------------------------------------------------
H     =  size(Gc,3)+1;       % Number of Hyperparameters (+ 1 for noise term)

% Caluclate Sufficient Stats on Y 
YY   =  Y*Y';                        % This is Suffient stats 1 (S1)
trYY =  sum(diag(YY));


% Figure out which terms to include into the model (set others to -32)
%--------------------------------------------------------------------------
if (isempty(ac))
    as    = 1:H;
    nas   = [];
else
    as    = find(ac);
    nas   = find(ac==0);
end;

% initialise h
%--------------------------------------------------------------------------
if (isempty(h0))
    rs     = Y-X*pinv(X)*Y; 
    Gd     =  zeros(Q,H-1);
    for i=1:H-1
        Gd(:,i)  = diag(Gc(:,:,i));
    end;
    h      =  ones(H,1)*low;
    u      =  pinv(Z)*rs;
    h(H,1) =  (trYY-traceAB(u'*Z',rs))/(P*N);
    D      =  u*u'/P-h(H)*pinv(Z'*Z)/P;  % Crude approx for variance-covariance matrix
    hD     =  diag(D);                   % Use the size of the diagnal values to get starting
    % Starting values for constrained estimates
    h(1:H-1,1) = real(log(pinv(Gd)*hD(:)));  % (xC'*xC)\(xC'*hD(:))
    h(h<low) = low+1;
else
    h      = h0;
end;
h(nas) = low;

% Initialize Interations and hyperpriors
%--------------------------------------------------------------------------
dF    = Inf;
dFdh  = zeros(H,1);
dFdhh = zeros(H,H);
HessReg = HessReg*speye(H,H);          % Prior precision (1/variance) of h


for k = 1:numIter
    % compute current estimate of covariance
    %----------------------------------------------------------------------
    G     = zeros(Q,Q);
    for i = as(1:end-1);
        G = G + Gc(:,:,i)*exp(h(i));
    end
    
    % Find the inverse of V - while dropping the zero dimensions in G
    [u,s] = eig(full(G));
    dS    = diag(s);
    idx   = dS>eps;
    Zu     = Z*u(:,idx);
    if (isempty(S)) 
        iV    = (eye(N)-Zu/(diag(1./dS(idx))*exp(h(H))+Zu'*Zu)*Zu')./exp(h(H)); % Matrix inversion lemma
    else 
        iV    = (S.invS-S.invS*Zu/(diag(1./dS(idx))*exp(h(H))+Zu'*S.invS*Zu)*Zu'*S.invS)./exp(h(H)); % Matrix inversion lemma
    end; 
    
    if (~isempty(X))
        iVX   = iV * X;
        iVr   = iV - iVX*((X'*iVX)\iVX');  % Correction for the fixed effects
    else
        iVr   = iV;
    end;
    
    % ReML estimate of hyperparameters
    
    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    A     = iVr*Z;
    B     = YY*A/P;
    for i = as(1:end-1)
        
        % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
        %------------------------------------------------------------------
        C{i}  = (A*Gc(:,:,i));
        CZ{i} = C{i}*Z';
        dFdh(i) = -P/2*(traceABtrans(C{i},Z)-traceABtrans(C{i},B));
    end
    dFdh(H) = -P/2*traceABtrans(iVr,(speye(N)-YY*iVr/P));
    
    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    if (isempty(S))
        iVrS = iVr; 
    else 
        iVrS = iVr*S.S; 
    end; 
    for i = 1:length(as)-1
        for j = i:length(as)-1
            
            % dF/dhh = -trace{iV*Z*G*Z'*iV*Z*G*Z'}
            %--------------------------------------------------------------
            dFdhh(as(i),as(j)) = -P/2*traceAB(CZ{as(i)},CZ{as(j)});
            dFdhh(as(j),as(i)) =  dFdhh(as(i),as(j));
        end
        dFdhh(as(i),H) = -P/2*traceABtrans(CZ{as(i)},iVrS);
        dFdhh(H,as(i)) =  dFdhh(as(i),H);
    end
    dFdhh(H,H) = -P/2*traceABtrans(iVrS,iVrS);
    
    
    % modulate
    %----------------------------------------------------------------------
    dFdh  = dFdh.*exp(h);
    dFdhh = dFdhh.*(exp(h)*exp(h)');
    
    % Add slight regularisation to second derivative 
    %----------------------------------------------------------------------
    dFdhh = dFdhh - HessReg;
    
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    % dh    = spm_dx(dFdhh(as,as),dFdh(as),{t}); % This is actually much
    % less accurate it seems than
    dh   = -dFdhh(as,as)\dFdh(as);
    h(as) = h(as) + dh;
    
    
    % predicted change in F - increase regularisation if increasing
    %----------------------------------------------------------------------
    pF    = dFdh(as)'*dh;
    dF    = pF;
    
    % convergence
    %----------------------------------------------------------------------
    % fprintf('%s %-23d: %10s%e [%+3.2f]\n','  ReML Iteration',k,'...',full(dF),t);
    if dF < thres
        break;
    else
        % eliminate redundant components (automatic selection)
        %------------------------------------------------------------------
        as  = find(h > low);
        h(h<low)=low;
        as  = as(:)';
    end
    hh(:,k)=h;
end

% return exp(h) and rescale
%--------------------------------------------------------------------------
theta  = h;
G     = zeros(Q);
for i = 1:H-1;
    G = G + full(Gc(:,:,i))*exp(theta(i));
end

% Find the inverse of V - while dropping the zero dimensions in G
% V = Z*G*Z' + I sigma
[u,s] = eig(full(G));
dS    = diag(s);
idx   = dS>eps;
Zu     = Z*u(:,idx);
if (isempty(S)) 
    iV    = (eye(N)-Zu/(diag(1./dS(idx))*exp(theta(H))+Zu'*Zu)*Zu')./exp(theta(H)); % Matrix inversion lemma
else 
    iV    = (S.invS-S.invS*Zu/(diag(1./dS(idx))*exp(theta(H))+Zu'*S.invS*Zu)*Zu'*S.invS)./exp(theta(H)); % Matrix inversion lemma
end; 

if (~isempty(X))
    iVX   = iV * X;
    iVr   = iV - iVX*((X'*iVX)\iVX');  % Correction for the fixed effects
else
    iVr   = iV;
end;

if nargout > 2
    u=G*Z'*iVr*Y;
end;
if ~all(isreal(u(:)));
    warning('U-estimates are not all real: should not happen'); 
    keyboard; 
end;

if nargout > 3
    ldet  = -2*sum(log(diag(chol(iV))));        % Safe computation of the log determinant (V) Thanks to code from D. lu
    l     = -P/2*(ldet)-0.5*traceABtrans(iVr,YY);
    if (~isempty(X)) % Correct for ReML estimates
        l = l - P*sum(log(diag(chol(X'*iV*X))));  % - P/2 log(det(X'V^-1*X));
    end; 
end;
