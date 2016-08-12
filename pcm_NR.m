function [G,theta,u,l,k]=pcm_NR(y,Z,varargin)
% pcm_NR: estimate random-effects variance component coefficients using
% Newton-Raphson gradient descent.
%
% function [G,theta,u,l,k]=pcm_NR(y,Z,varargin)
% Estimates the variance coefficients of the model described in:
% Diedrichsen, Ridgway, Friston & Wiestler (2011).
%
% y_n = X b_n + Z u_n + e,
%         u ~ (b, G)
%           G = sum (theta * Gc)
%
% y: N x P observations
% Z: N x Q random effects matrix
%
% N: numbers of observations (trials)
% P: numbers of experimental units (voxels)
% Q: number of random effects
%
% VARARGIN:
%   'numIter'     : Maximal number of iterations
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
%                    ??? Important: X will also be removed from Z to orthogonilise the
%                    random effects from the constant effects, making the
%                    estimation of b_n independent of G
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
meanS   = 0;                    % Mean subtract
ac      = [];                      % Which terms do I want to include?
numIter = 32;                 % Maximal number of iterations
Gc      = {};
h0      = [];
low     = -16;
thres   = 1e-2;
hE      = -8;                 % Hyper prior: Set of the floor to collect evidence
hP      = 1/256;              % Precision on hyper prior
X       = [];                 % By default fixed effects are empty

% Variable argument otions
%--------------------------------------------------------------------------
vararginoptions(varargin, ...
    {'Gc','meanS','h0','ac','hE','hP','thres','X','numIter'});

% check input size
%--------------------------------------------------------------------------
[N,P]  =  size(y);
[N2,Q] =  size(Z);
if N2  ~= N
    error('Mismatched numbers of rows in data (%d) and design (%d)', N, N2)
end

% Intialize the Model structure
%--------------------------------------------------------------------------
H     =  length(Gc)+1;       % Number of Hyperparameters (+ 1 for noise term)
for i =  1:H-1
    Gc{i} = sparse(Gc{i});
end;

% If necessary, subtract the fixed effects estimates (a)
%--------------------------------------------------------------------------
if (meanS)
    a  =  pinv(Z)*sum(y,2)/P;
    r  =  bsxfun(@minus,y,Z*a);
else
    r  =  y;
end;

YY   =  r*r';                        % This is Suffient stats 1 (S1)
trYY =  sum(diag(YY));


% Scale YY
%--------------------------------------------------------------------------
sY = 1; % trace(YY)/(P*N);
sYY = YY/sY;
trYY = trYY/sY;
rs=r/sqrt(sY);

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
    xC     =  zeros(Q,H-1);
    for i=1:H-1
        xC(:,i)  = diag(Gc{i});
    end;
    h      =  ones(H,1)*low;
    u      =  pinv(Z)*rs;
    h(H,1) =  (trYY-traceAB(u'*Z',rs))/(P*N);
    D      =  u*u'/P-h(H)*pinv(Z'*Z)/P;  % Crude approx for variance-covariance matrix
    hD     =  diag(D);                   % Use the size of the diagnal values to get starting
    % Starting values for constrained estimates
    h(1:H-1,1) = real(log(pinv(xC)*hD(:)));  % (xC'*xC)\(xC'*hD(:))
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
hE = hE*ones(H,1);             % Prior mean of h
hP = hP*speye(H,H);          % Prior precision (1/variance) of h


for k = 1:numIter
    % compute current estimate of covariance
    %----------------------------------------------------------------------
    G     = sparse(Q,Q);
    for i = as(1:end-1);
        G = G + Gc{i}*exp(h(i));
    end
    
    % Find the inverse of V - while dropping the zero dimensions in G
    [u,s] = eig(full(G));
    dS    = diag(s);
    idx   = dS>eps;
    Zu     = Z*u(:,idx);
    iV    = (eye(N)-Zu/(diag(1./dS(idx))*exp(h(H))+Zu'*Zu)*Zu')./exp(h(H)); % Matrix inversion lemma
    
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
    B     = sYY*A/P;
    for i = as(1:end-1)
        
        % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
        %------------------------------------------------------------------
        C{i}  = (A*Gc{i});
        CZ{i} = C{i}*Z';
        dFdh(i) = -P/2*(traceABtrans(C{i},Z)-traceABtrans(C{i},B));
    end
    dFdh(H) = -P/2*traceABtrans(iVr,(speye(N)-sYY*iVr/P));
    
    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = 1:length(as)-1
        for j = i:length(as)-1
            
            % dF/dhh = -trace{iV*Z*G*Z'*iV*Z*G*Z'}
            %--------------------------------------------------------------
            dFdhh(as(i),as(j)) = -P/2*traceAB(CZ{as(i)},CZ{as(j)});
            dFdhh(as(j),as(i)) =  dFdhh(as(i),as(j));
        end
        dFdhh(as(i),H) = -P/2*traceABtrans(CZ{as(i)},iVr);
        dFdhh(H,as(i)) =  dFdhh(as(i),H);
    end
    dFdhh(H,H) = -P/2*traceABtrans(iVr,iVr);
    
    
    % modulate
    %----------------------------------------------------------------------
    dFdh  = dFdh.*exp(h);
    dFdhh = dFdhh.*(exp(h)*exp(h)');
    
    % add hyperpriors
    %----------------------------------------------------------------------
    e     = h     - hE;
    dFdh  = dFdh  - hP*e;
    dFdhh = dFdhh - hP;
    
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
theta  = sY*exp(h);
G     = zeros(Q);
for i = 1:H-1;
    G = G + full(Gc{i})*theta(i);
end

% Find the inverse of V - while dropping the zero dimensions in G
% V = Z*G*Z' + I sigma
[u,s] = eig(full(G));
dS    = diag(s);
idx   = dS>eps;
Zu     = Z*u(:,idx);
iV    = (eye(N)-Zu/(diag(1./dS(idx))*theta(H)+Zu'*Zu)*Zu')./theta(H); % Matrix inversion lemma

if (~isempty(X))
    iVX   = iV * X;
    iVr   = iV - iVX*((X'*iVX)\iVX');  % Correction for the fixed effects
else
    iVr   = iV;
end;

if nargout > 2
    u=G*Z'*(iV)*r;
end;
if ~all(isreal(u(:)));
    keyboard;
end;

if nargout > 3
    ldet  = -2*sum(log(diag(chol(iV))));        % Safe computation of the log determinant (V) Thanks to code from D. lu
    l     = -P/2*(ldet)-0.5*traceABtrans(iVr,YY);
    if (~isempty(X)) % Correct for ReML estimates
        l = l - P*sum(log(diag(chol(X'*iV*X))));  % - P/2 log(det(X'V^-1*X));
    end; 
end;
