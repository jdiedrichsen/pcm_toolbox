function [G,h,u,l,n,jumpI,a]=pcm_EM(y,Z,varargin)
% pcm_EM: estimate random-effects variance component coefficients
% Usage: [G,h,u,l,n,jumpI,a]=mvpattern_covcomp(y,Z,varargin);
% Estimates the variance coefficients of the model described in:
% Diedrichsen, Ridgway, Friston & Wiestler (2011).
% 
% y_n = X b_n + Z u_n + e,
%         u ~ (a, G)
%                 G = A A'
%
% y: N x P observations (doubles) 
% Z: N x Q random effects matrix (doubles)
%
% N: numbers of observations (trials)
% P: numbers of experimental units (voxels)
% Q: number of random effects
%
% VARARGIN:
%   'num_iter'     : Maximal number of iterations
%   'Ac'           : Cell array (Hx1) of components of factor loading 
%                    matrix A = sum(h_m A_m), used to form G = A A'
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
% See mva_component_examples
%  
% See also: pcm_EM_free, pcm_NR, pcm_NR_diagfree
% This function was formly known as mvpattern_covcomp
% 
%
% v.2.0: speed up (30%) achieved by using mex-version of
%           ML_constrained_fast
%          
% Copyright 2011 Joern Diedrichsen, j.diedrichsen@ucl.ac.uk

% check input size
%--------------------------------------------------------------------------
[N,P]  =  size(y);
[N2,Q] =  size(Z);
if N2  ~= N
    error('Mismatched numbers of rows in data (%d) and design (%d)', N, N2)
end

% defaults 
%--------------------------------------------------------------------------
num_iter     = 600;         % Maximal number of iterations
Ac           = {};          % This is the variance-covariance structure
TolL         = 0.00001;     % Tolerance on Likelihood
accel_method = 'Aitken';    % Acceleration method         
meanS=1;                    % Mean subtract
h0=[];                      % Starting value 

% Variable argument otions
%--------------------------------------------------------------------------
pcm_vararginoptions(varargin, ...
    {'num_iter','Ac','TolL','accel_method','meanS','h0','X'});

% Intialize the Model structure
%--------------------------------------------------------------------------
H    = length(Ac);       % Number of Hyperparameters (without noise)          
Cc   = cell(H, 1);       % Z*Ac cell
sCc  = cell(H, 1);       % Sparse version of Z*Ac 

xC   = zeros(Q*Q,H);
for i=1:H
    Cc{i}   =  Z*Ac{i};           % Calculate regression coefficients for full data matrix 
    sCc{i}  =  sparse(Cc{i});     % Sparse version 
    xC(:,i) =  Ac{i}(:);          % Vector version for intial starting estimate 
end;

% Precompute cross-terms for M-step
%--------------------------------------------------------------------------
CcCc  = cell(H, H);
sCcCc = cell(H, H);
for i=1:H
    for j=i:H
        CcCc{i,j}   =  Cc{j}'*Cc{i}; % Calculate cross terms
        sCcCc{i,j}  =  sparse(CcCc{i,j}); 
    end;
end;

% If fixed effects are given, remove fixed effects from data and random
% effects matrix 
%--------------------------------------------------------------------------
if (~isempty(X))
    pX = pinv(X); 
    Z  = Z-X*pX*Z; 
    y  = y-X*pX*y; 
    error('ReML estimation not implemented yet!'); 
end; 


% If necessary, subtract the mean value of the random effects (a)
%--------------------------------------------------------------------------
if (meanS)
    a  =  pinv(Z'*Z*P)*Z'*sum(y,2);
    r  =  y-repmat(Z*a,1,P);
else
    a  =  zeros(Q,1);
    b  =  zeros(size(Z,2),1);
    r  =  y-repmat(Z*b,1,P);
end;
rr   =  r*r';                        % This is Suffient stats 1 (S1)
trRR =  sum(diag(rr)); 

% preallocate arrays for speed
%--------------------------------------------------------------------------
l       = zeros(1,num_iter);
h       = zeros(H+1,num_iter); % h(H+1) = Sigma2
delta_h = zeros(H+1,num_iter);

% Provide initial guess from Laird, Lange and Dempster
%--------------------------------------------------------------------------
if (isempty(h0))
    u        = pinv(Z'*Z)*Z'*r;
    h(H+1,1) = (sum(sum(y.*y,2))-a'*Z'*sum(y,2)-sum(sum((u'*Z').*r')))/(P*N);
    D        = u*u'/P-h(H+1,1)*pinv(Z'*Z)/P;
    hD       = diag(diag(D).^0.5);
    % Starting values for constrained estimates
    h(1:H,1) = real((xC'*xC)\(xC'*hD(:)));  % Remove any problems when one components is numerically 0
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

C=zeros(N,Q);
for i=1:H
    C  = C+Cc{i}*h(i,1);
end;

% Iterate
%--------------------------------------------------------------------------
while (n<num_iter && diffL>TolL)
    
    % Estep
    %--------------------------------------------------------------------------
    V  =  h(H+1,n)*eye(N)+C*C';
    Wr =  V\r;
    WC =  V\C;
    v  =  C'*Wr;
    rv =  r*v';                            % <rv'> suffienct stats (S2)
    vv = (v*v')+P*eye(Q)-P*C'*WC;          % <vv'> sufficient statistics (S3)
	try
        ldet = 2* sum(log(diag(chol(V))));        % Safe computation of the log determinant Thanks to code from D. lu
        l(n)=-P/2*(ldet)-0.5*sum(sum(Wr.*r,2));
    catch
        l(n)=-inf; 
    end; 
        
    % Check if likelihood decreased on the last iteration
    if (n>1 && (l(n)<l(n-1) || h(H+1,n)<0))
        
        % Check if last iteration was a jump
        if (~jumpI(n)==1) % It wasn't: that's bad it should not decrease
            % warning('mvpattern_covcomp:EMdecrease', ...
            %     'EM decreased by %g', l(n-1)-l(n));
            diffL=0;
            if (l(n-1)-l(n)>0.01)
                % If it is only a small decrease, it may be rouding error
                error('likelihood decreased!');
            end;
        else
            % Last step was a jump: So go undo the jump and redo the E-step
            n=n-1;
            
            % Build up current C-matrix
            C=zeros(N,Q);
            for i=1:H
                C=C+Cc{i}*h(i,n);
            end;
            
            % Estep
            V   =  h(H+1,n)*eye(N)+C*C';
            Wr  =  V\r;
            WC  =  V\C;
            v   =  C'*Wr;
            rv  =  r*v';
            vv  =  (v*v')+P*eye(Q)-P*C'*WC; % <vv'> sufficient statistics
            ldet = 2* sum(log(diag(chol(V))));
            l(n)= -P/2*(ldet)-0.5*sum(sum(Wr.*r,2));
        end;
    end;
    
    % Mstep: Constrained regression: here it uses new MEX function 
    %--------------------------------------------------------------------------
    [COV,VA]=pcm_ML_constrained_fast(rv,vv,sCc,sCcCc); 
    h(1:H,n+1) = VA\COV;
    
    % Based on the new h-parameters, build up C
    %--------------------------------------------------------------------------
    C=zeros(N,Q);
    for i=1:H
        C=C+Cc{i}*h(i,n+1);
    end;
    
    % Estimate sigma_e
    %--------------------------------------------------------------------------
    % h(H+1,n+1)=1/(N*P)*trace(rr-C*rv');
    h(H+1,n+1)=1/(N*P)*(trRR-sum(sum(C.*rv)));

    % Track the change in parameters
    %--------------------------------------------------------------------------
    delta_h(:,n+1)=h(:,n+1)-h(:,n);
    
    % Track the change in likelihood as a stopping criterion.
    % Do not abort for small steps in liklihood, but only when the
    % real likelihood is close to the estimated maximal likelihood.
    % This prevents stopping of the iteration when there is only slow
    % progress (see McLaughlan & Krishnan, 1997. The EM alogithm and
    % Extensions)
    %--------------------------------------------------------------------------
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
                % Build up new C-matrix 
                C=zeros(N,Q);
                for i=1:H
                    C=C+Cc{i}*h(i,n);
                end;
            end
        end
    end
end
jumpI(n+1:end)=[]; % trim down from num_iter length if diffL converged

% Now build the G-matrix of random effects
%--------------------------------------------------------------------------
% First build the full A matrix, then build G
% Note that: 
% sum(Ac_i*Ac_i'*h_i*h_i) != sum(Ac_i*h_i)*sum(Ac_i*h_i)'! 
h=h(:,1:n-1);
l=l(1:n-1);
A=zeros(Q,Q);
for i=1:H
    A=A+Ac{i}*h(i,end);
end;
G=A*A';

% Now estimate the regression coefficients 
%--------------------------------------------------------------------------
% u=G*Z*inv(Z*G*Z'+eye(N)*sigma2)*r
% This is equal by matrix inversion lemma to the ridge regression form: 
% u=inv(Z'*Z+inv(G)*sigma2)*Z'*r
u=G*Z'*Wr;


