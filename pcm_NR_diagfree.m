function [G,th,u,l,k]=pcm_NR_diagfree(y,Z,varargin)
% function [G,th,u,l,k]=pcm_NR_diagfree(y,Z,varargin)
% pcm_NR_diagfree: estimate random-effects variance component coefficients
% using Newton-Raphson algorithm, where G has diagnonal structure and all
% all elements of the diagonal are free to vary differently
%
% y_n = X b_n + Z u_n + e,
%         u ~ (b, G)
%           G = sum (h * Gc)
%
% y: N x P observations
% Z: N x Q random effects matrix
% X: N x L fixed effects matrix (optional)
%
% N: numbers of observations (trials)
% P: numbers of experimental units (voxels)
% Q: number of random effects
%
% VARARGIN:
%   'num_iter'     : Maximal number of iterations
%   'h0'           : Starting values for the parameters (Hx1 vector)
%                    (otherwise uses starting guess based on Laird & Lange)
%   'TolL'         : Tolerance of the likelihood (l-l'), where l' is
%                    projected likelihood
%   'accel_method' : Acceleration method used
%                    'none'   : unaccelerated EM algorithm
%                    'Aitken' : Aitken acceleration with projected jump
%                               based on observed convergence
%   'meanS'        : Remove the mean for each pattern component (a)
%                    (Logical flag, true by default)
%   'X'            : Fixed effects matrix that will be removed from y
%                    In this case, ReML will be used for the estimation of
%                    G.
%                    Important: X will also be removed from Z to orthogonilise the
%                    random effects from the constant effects, making the
%                    estimation of b_n independent of G
%
% OUTPUT:
%   G     : variance-covariance matrix
%   th    : coefficients (one column for each iteration)
%   u     : hidden patterns components
%   l     : likelihood
%
% Examples:
% See mva_component_examples
%
% See also: mva_component_examples, spm_reml, spm_reml_sc, spm_sp_reml
% Where spm_* are from the SPM software, http://www.fil.ion.ucl.ac.uk/spm
%
% v.1.0:
%
% Copyright 2014 Joern Diedrichsen, j.diedrichsen@ucl.ac.uk


%--------------------------------------------------------------------------
meanS   = 1;                    % Mean subtract
ac      = [];                      % Which terms do I want to include?
numIter = 32;                 % Maximal number of iterations
Gc      = {};
h0      = [];
low     = -16;
thres   = 1e-2;
hE      = -8;                 % Hyper prior: Set of the floor to collect evidence
hP      = 1/256;
X       = [];                   % Design matrix to be removed 

% Variable argument otions
vararginoptions(varargin, ...
    {'meanS','Gc','h0','ac','hE','thres','X'});

% check input size
%--------------------------------------------------------------------------
[N,P]  =  size(y);
[N2,Q] =  size(Z);
if N2  ~= N
    error('Mismatched numbers of rows in data (%d) and design (%d)', N, N2)
end

% Intialize the Model structure
%--------------------------------------------------------------------------
H     =  Q+1;       % Number of Hyperparameters (+ 1 for noise term)
dF    = Inf;
dFdh  = zeros(H,1);
dFdhh = zeros(H,H);
hE = hE*ones(H,1);             % Prior mean of h
hP = hP*speye(H,H);          % Prior precision (1/variance) of h


% If fixed effects are given, remove fixed effects from data and random
% effects matrix
%--------------------------------------------------------------------------
if (~isempty(X))
    pX = pinv(X);
    Z  = Z-X*pX*Z;
    y  = y-X*pX*y;
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


% scale Gc and YY
%--------------------------------------------------------------------------
sY = trYY/(P*N);
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
    xC     =  eye(Q);
    h      =  ones(H,1)*low;
    u      =  pinv(Z)*rs;
    %     h(H,1) =  (trYY-sum(sum((u'*Z').*rs')))/(P*N);
    h(H,1) =  (trYY-traceAB(u'*Z',rs))/(P*N);
    D      =  u*u'/P-h(H)*pinv(Z'*Z)/P;
    hD     =  diag(D);
    % Starting values for constrained estimates
    h(1:H-1,1) = real(log((xC'*xC)\(xC'*hD(:))));
    h(h<low) = low+1;
else
    h      = h0;
end;
h(nas) = low;

% Start iteration
%--------------------------------------------------------------------------
for k = 1:numIter
    % compute current estimate of covariance
    %----------------------------------------------------------------------
    G     = diag(exp(h(1:Q)));
    
    % Find the inverse of V - while dropping the zero dimensions in G
    [u,s] = eig(G);
    dS    = diag(s);
    idx   = dS>eps;
    Zu    = Z*u(:,idx);
    iV    = (eye(N)-Zu/(diag(1./dS(idx))*exp(h(H))+Zu'*Zu)*Zu')./exp(h(H)); % Matrix inversion lemma
    if (~isempty(X))
        iVX   = iV * X;
        iVr   = iV - iVX*((X'*iVX)\iVX');  % Correction for the fixed effects
    else
        iVr   = iV;
    end;
    
    % l(k)=P/2*(log(det(iV)))-0.5*sum(sum(iV.*sYY,2));
    
    % M-step: ReML estimate of hyperparameters
    %======================================================================
    
    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    A     = iVr*Z;
    B     = sYY*A/P;
    dFdh(1:Q) = -P/2*(sum(A.*Z)-sum(A.*B))';
    dFdh(H) = -P/2*traceABtrans(iVr,(speye(N)-sYY*iVr/P));
    
    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    dFdhh(1:Q,1:Q) = -P/2*((Z'*A).*(Z'*A)); % Short form
    dFdhh(H,H)     = -P/2*traceABtrans(iVr,iVr);
    dFdhh(1:Q,H)   = -P/2*diag(Z'*iVr*A);
    dFdhh(H,1:Q)   = dFdhh(1:Q,H)';
    
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
    
    % test for convergence
    %----------------------------------------------------------------------
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
end;

% return exp(h) and rescale
%--------------------------------------------------------------------------
th  = sY*exp(h);

% Return G and u if necessary
%--------------------------------------------------------------------------
G  = diag(th(1:Q));
if nargout > 2
    u=G*Z'*(iV/sY)*r;
end;
if nargout > 3
    iV  = (eye(N)-Z/(inv(G)*th(H)+Z'*Z)*Z')./th(H); % Matrix inversion lemma
    ldet = - 2* sum(log(diag(chol(iV))));        % Safe computation of the log determinant Thanks to code from D. lu
    l   = P/2*(ldet)-0.5*traceABtrans(iV,YY);
    if (~isempty(X)) % Correct for ReML estimates
        l = l - sum(log(diag(chol(X'*W*X))));  % - 0.5 * log(det());
    end;
end;

