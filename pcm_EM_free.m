function [G,h,u,l,n,jumpI,a]=pcm_EM_free(y,Z,varargin)
% estimate random-effects variance component coefficients
%
% Usage: [G,h,u,l,n,jumpI,a]=pcm_EM_free(y,Z,varargin);
%
% Estimates the variance coefficients of the model described in:
% Laird, Lange & Stram (1987).
% y_n = X b_n + Z u_n + e,
%         u ~ (b, G)
% G is any arbitrary p.d. matrix
%
% y: N x P observations
% Z: N x Q random effects matrix
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
%   h     : coefficients (one column for each iteration)
%   u     : hidden patterns components
%   l     : likelihood
%   n     : number of iterations until convergence
%   jumpI : record of jumps in convergence
%   a     : Fixed effect means for the pattern components
%
% Examples:
%
% See also: mvpattern_covcomp, mva_component_examples, spm_reml, spm_reml_sc, spm_sp_reml
% Where spm_* are from the SPM software, http://www.fil.ion.ucl.ac.uk/spm
%
% v.1.0:
%
% Copyright 2012 Joern Diedrichsen, j.diedrichsen@ucl.ac.uk

% check input size
%--------------------------------------------------------------------------
[N,P]   =  size(y);
[N2,Q]  =  size(Z);
if N2   ~= N
    error('Mismatched numbers of rows in data (%d) and design (%d)', N, N2)
end

% defaults
%--------------------------------------------------------------------------
num_iter     =  600;       % Maximal number of iterations
TolL         =  0.00001;            % Tolerance on Likelihood
accel_method = 'Aitken';
meanS        = 1;                    % Mean subtract (across voxels)
h0           = [];                   % Forced starting value
X            = [];
% Variable argument otions
vararginoptions(varargin, ...
    {'num_iter','TolL','accel_method','meanS','h0','X'});

% Number of hyper parameters are elements in G-matrix
%--------------------------------------------------------------------------
H=Q*Q;    % Full G-matrix

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
    a=pinv(Z'*Z*P)*Z'*sum(y,2);
    r=y-repmat(Z*a,1,P);
else
    a=zeros(Q,1);
    b=zeros(size(Z,2),1);
    r=y-repmat(Z*b,1,P);
end;
rr=r*r';                        % This is Suffient stats 1 (S1)
trRR=sum(diag(rr));

% preallocate arrays for speed
%--------------------------------------------------------------------------
l=zeros(1,num_iter);
h=zeros(H+1,num_iter); % h(H+1) = Sigma2
delta_h=zeros(H+1,num_iter);

% Provide initial guess from Laird, Lange and Dempster
%--------------------------------------------------------------------------
if (isempty(h0))
    u          = pinv(Z'*Z)*Z'*r;
    h(H+1,1)   = (sum(sum(y.*y,2))-a'*Z'*sum(y,2)-sum(sum((u'*Z').*r')))/(P*N);
    G          = u*u'/P-h(H+1,1)*pinv(Z'*Z)/P;
    % Starting values for constrained estimates
    h(1:H,1)   = G(:);
else
    h(1:H+1,1) = h0;
end;

% Initialize
%--------------------------------------------------------------------------
n        = 1;
jump     = 1;
jumpI    = nan(1, num_iter);
jumpI(1) = 1;
diffL    = inf;

% Iterate
%--------------------------------------------------------------------------
while (n<num_iter && diffL>TolL)
    
    
    % Estep
    %--------------------------------------------------------------------------
    G  = reshape(h(1:H,1),Q,Q);
    V  = h(H+1,n)*eye(N)+Z*G*Z';
    W  = pinv(V);
    Wr = W * r;
    u  = G*Z'*Wr;
    try
        ldet = 2* sum(log(diag(chol(V))));        % Safe computation of the log determinant Thanks to code from D. lu
        l(n) = -P/2*(ldet)-0.5*sum(sum(Wr.*r,2));
        if (~isempty(X))                          % Correction for ReML estimation
            l(n) = l(n) - sum(log(diag(chol(X'*W*X))));  % - 0.5 * log(det());
        end;
    catch
        l(n)=-inf;
    end;
    
    % Check if likelihood decreased on the last iteration
    %--------------------------------------------------------------------------
    if (n>1 && l(n)<l(n-1) || h(H+1,n)<0)
        
        % Check if last iteration was a jump
        if (~jumpI(n)==1) % It wasn't: that's bad it should not decrease
            % warning('mvpattern_covcomp:EMdecrease', ...
            %    'EM decreased by %g', l(n-1)-l(n));
            diffL=0;
            if (l(n-1)-l(n)>0.01)
                % If it is only a small decrease, it may be rouding error
                error('likelihood decreased!');
            end;
        else
            % Last step was a jump: So go undo the jump and redo the E-step
            n=n-1;
            
            % Estep
            V  = h(H+1,n)*eye(N)+Z*G*Z';
            W  = pinv(V);
            Wr = W * r;
            u  = G * Z' * Wr;
            try
                ldet = 2* sum(log(diag(chol(V))));        % Safe computation of the log determinant Thanks to code from D. lu
                l(n) = -P/2*(ldet)-0.5*sum(sum(Wr.*r,2));
                if (~isempty(X))                          % Correction for ReML estimation
                    l(n) = ln - sum(log(diag(chol(X'*W*X))));  % - 0.5 * log(det());
                end;
            catch
                l(n) = -inf;
            end;
            
        end;
    end;
    
    % Mstep: Update on D (see Laird, Lange & Stram): do ML or ReML
    %--------------------------------------------------------------------------
    if (~isempty(X))  % Correct W for the influence of X
        W = W - W*X*pinv(X'*W*X)*X'*W;
    end;
    G  = (u*u'+G*(eye(Q)-Z'*W*Z*G))/P;
    h(1:H,n+1) = G(:);
    % Estimate sigma_e
    % h(H+1,n+1)=1/(N*P)*trace(rr-C*rv');
    rr = r-Z*u;
    h(H+1,n+1) = sum(sum(rr.*rr))/(P*N)+h(H+1,n)*trace(eye(N)-h(H+1,n)*W)/N;
    
    % Track the change in parameters
    %--------------------------------------------------------------------------
    % Track the change in likelihood as a stopping criterion.
    % Do not abort for small steps in liklihood, but only when the
    % real likelihood is close to the estimated maximal likelihood.
    % This prevents stopping of the iteration when there is only slow
    % progress (see McLaughlan & Krishnan, 1997. The EM alogithm and
    % Extensions)
    delta_h(:,n+1)=h(:,n+1)-h(:,n);    
    if (n-jump>2)                           % Check convergence by
        Rdl=(l(n)-l(n-1))/(l(n-1)-l(n-2));  % Ratio of differences in l
        lA=l(n-1)+1./(1-Rdl)*(l(n)-l(n-1)); % Predicted maximal likelihood
        diffL=lA-l(n);                      % Estimated deviation
    end;
    
    % prepare next iteration
    %--------------------------------------------------------------------------
    n=n+1;
    jumpI(n)=0;
    
    % See if a jump can be done based on recent progress
    %--------------------------------------------------------------------------
    if (strcmp(accel_method,'Aitken'))
        if (n-jump>3)
            lambda_hat=mean(delta_h(:,n)./delta_h(:,n-1));
            if (lambda_hat<1)
                h(:,n+1)=h(:,n-1)+1./(1-lambda_hat)*(h(:,n)-h(:,n-1));
                l(n)=l(n-1);
                n=n+1;
                jumpI(n)=1;
                jump=n;
            end
        end
    end
end
jumpI(n+1:end)=[]; % trim down from num_iter length if diffL converged
h=h(:,1:n-1);
l=l(1:n-1);


function vararginoptions(options,allowed_vars,allowed_flags)
% function vararginoptions(options,allowed_vars,allowed_flags);
% Deals with variable argument in
% INPUTS
%   options: cell array of a argument list passed to a function
%   allowed_vars: Variables that can be set
%   allowed_flags: Flags that can be set
%  vararginoptions assigns the value of the option to a variable with the
%  name option (in called workspace).
%  Flags are set to one:
% EXAMPLE:
%   the option-string 'var1',4,'var2',10,'flag'
%   causes the var1 and var2 to be set to 4 and 10 and flag to 1
%   if allowedvars are not given, all variables are allowed
% Joern Diedrichsen
% v1.0 9/13/05
checkflags=1;
checkvars=1;
if nargin<2
    checkvars=0;
end;
if nargin<3
    checkflags=0;
end;

c=1;
while c<=length(options)
    a=[];
    if ~ischar(options{c})
        error('Options must be strings on argument %d',c);
    end;
    if checkflags
        a=find(strcmp(options{c},allowed_flags), 1);
    end;
    if ~isempty(a)
        assignin('caller',options{c},1);
        c=c+1;
    else
        if checkvars
            a=find(strcmp(options{c},allowed_vars), 1);
            if (isempty(a))
                error(['unknown option:' options{c}]);
            end;
        end;
        if (c==length(options))
            error('Option %s must be followed by a argument',options{c});
        end;
        assignin('caller',options{c},options{c+1});
        c=c+2;
    end;
end;
