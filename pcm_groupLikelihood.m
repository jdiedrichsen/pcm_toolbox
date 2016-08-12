function [negLogLike,dnldtheta,L,dLdtheta] = pcm_groupLikelihood(theta,YY,M,Z,X,P,varargin);
% function [negLogLike,dnldtheta,L,dLdtheta] = pcm_groupLikelihood(theta,YY,M,Z,X,P,varargin);
% Returns negative log likelihood of the data from a number of subjects (YY) under the
% representational model specified by M.
% The representational model specifies some group-level parameters, i.e the
% structure of the second-moment matrix is supposed to be stable across all
% subjects.
% To link the group prediction to individual second-moment matrices, a
% number of individual-level nuisance parameters are introduced:
%   s:        Scaling of the data
%   sigma2:   Noise variability
% Additionally, if 'runEffect' is set, the effect of each imaging run is
%   also modelled as a random effect, giving the last individual-subject
%   parameter
% INPUT:
%      theta:   Vector of model parameters + subject-specific nuisance parameters
%               1...numGParams           : model parameters
%               numGParams+1... +numSubj : log(scaling parameters)
%               numGParams+numSubj+1...  : log(noise variance)
%               numGParams+numSubj*2+1...: log(run effect variance)
%                                          - used if runEffect is provided
%      YY:      cell array {numSubj}: Y*Y': Outer product of the data (NxN)
%      M:       Model for G (same for the whole group):
%               a. this can either be a fixed KxK matrix (fixed model) or a
%               b. Structure with fields
%                   M.type      = 'fixed','component','squareroot','nonlinear'
%                   M.numGparam = number of free parameters of G
%                   M.modelpred = handle to a fcn [G, dGdtheta] = fcn(theta);
%                   M.Gc        = Component matrices
%                   M.Ac        = Component matrices for matrix square root
%      Z:       cell array {numSubj}: Design matrix for random effects,
%                    variables of interest - often condition to trials (N) matrix
%      X:       cell array {numSubj}: Design matrix for fixed effect,
%                     accounted for through ReML.
%      P:       numSubj x 1 vector: Number of voxels
% VARARGIN:
%      'runEffect',B:  cell array {numSubj} of design matrices for the run effect,
%                   which is modelled as a individual subject-specific random effect.
%      'S',S:       Optional structured noise matrix - default is the identity matrix
%                   S should be a structure-array (number of subjects) with two fields
%                   S.S: noise structure matrix
%                   S.invS: Inverse of noise structure matrix
%      'verbose',0/1: Set verbose mode -
% OUTPUT:
%      negLogLike:  Negative Log likelihood of all subject's data
%                   We use the negative log liklihood to be able to plug the function into
%                   minimize or other optimisation routines.
%      dnldtheta:   Derivative of the negative log-likelihood
%      L:           Log likelihood (not inverted) for all the subject
%      dLdtheta:    Derivate of Log likelihood for each subject
%
%   Joern Diedrichsen, 6/2016, joern.diedrichsen@googlemail.com
%

S         = [];
verbose   = 0;
runEffect = [];
vararginoptions(varargin,{'S','verbose','runEffect'});

% Get G-matrix and derivative of G-matrix in respect to parameters 
if (isstruct(M))
    [G,dGdtheta] = pcm_calculateG(M,theta(1:M.numGparams)); 
else 
    G=M;
    M=[]; 
    M.numGparams=0; 
end; 

% Loop over subjects
numSubj = length(YY);
for s=1:numSubj
    % Get parameter and sizes
    scaleParam = theta(M.numGparams+s);   % Subject Scaling Parameter
    noiseParam = theta(M.numGparams+numSubj+s);   % Subject Noise Parameter
    N  = size(Z{s},1);
    Gs = G*exp(scaleParam);         % Scale the subject G matrix up by individual scale parameter
    
    % If Run effect is to ne modelled as a random effect - add to G and
    % design matrix
    if (~isempty(runEffect{s}))
        numRuns = size(runEffect{s},2);
        runParam = theta(M.numGparams+numSubj*2+s);    % Subject run effect parameter
        Gs = blockdiag(Gs,eye(numRuns)*exp(runParam));  % Include run effect in G
        Z{s} = [Z{s} runEffect{s}];                 % Include run effect in design matrix
    else
        numRuns = 0;                                % No run effects modelled
    end;
    
    % Find the inverse of the G-matrix
    Gs = (Gs+Gs')/2;        % Symmetrize
    [u,lam_Gs] = eig(full(Gs));
    dS    = diag(lam_Gs);
    idx = dS>1e-5;
    Zu  = Z{s}*u(:,idx);
    
    % Give warning if G was not invertible
    if (verbose)
        if (any(dS<-1e-5))
            %error('negative eigenvalues in G');
            warning('negative eigenvalues in G: num=%d, min=%d',sum(dS<-1e-5),dS(dS<-1e-5)); % ay
            if isstruct(M)
                if isfield(M,'name') % ay
                    fprintf('\t Model: %s',M.name);
                else
                    fprintf('\t Model: %s',char(M.modelpred));
                end
            end
        end;
    end
    
    % Compute inv(V) over matrix inversion lemma - use particular noise
    % structure if given.
    if (~isempty(S));
        sS = S(s).S; % noise covariance
        iS = S(s).invS; % inverse noise covariance
        iV  = (iS-iS*Zu/(diag(1./dS(idx))*exp(noiseParam)+Zu'*iS*Zu)*Zu'*iS)./exp(noiseParam); % Matrix inversion lemma
    else
        iV  = (eye(N)-Zu/(diag(1./dS(idx))*exp(noiseParam)+Zu'*Zu)*Zu')./exp(noiseParam); % Matrix inversion lemma
    end;
    iV  = real(iV); % sometimes iV gets complex (ay)
    
    % Deal with the fixed effects, if they are present
    if (isempty(X) || isempty(X{s}))
        iVr   = iV;
    else
        iVX   = iV * X{s};
        iVr   = iV - iVX*((X{s}'*iVX)\iVX');
    end;
    
    % Precompute some matrices
    A     = iVr*Z{s};
    B     = YY{s}*A/P(s);
    % Get the derivatives for all the parameters
    for i = 1:M.numGparams
        C{i}  = (A*blockdiag(dGdtheta(:,:,i),zeros(numRuns)));
        dLdtheta(i,s) = -P(s)/2*(traceABtrans(C{i},Z{s})-traceABtrans(C{i},B))*exp(scaleParam);
    end
    
    % Get the derivatives for the scaling parameters
    indx             = M.numGparams+s;  % Which number parameter is it?
    C{indx}          = A*blockdiag(G,zeros(numRuns));
    dLdtheta(indx,s) = -P(s)/2*(traceABtrans(C{indx},Z{s})-traceABtrans(C{indx},B))*exp(scaleParam);
    
    % Get the derivatives for the Noise parameters
    indx             = M.numGparams+numSubj+s;  % Which number parameter is it?
    if (isempty(S))
        dLdtheta(indx,s)     = -P(s)/2*traceABtrans(iVr,(speye(N)-YY{s}*iVr/P(s)))*exp(noiseParam);
    else
        dLdtheta(indx,s)     = -P(s)/2*traceABtrans(iVr*S(s).S,(speye(N)-YY{s}*iVr/P(s)))*exp(noiseParam);
    end;
    
    % Get the derivatives for the block parameter
    if (~isempty(runEffect) && ~isempty(runEffect{s}))
        indx             = M.numGparams+numSubj*2+s;  % Which number parameter is it?
        C{indx}          = A*blockdiag(zeros(size(G,1)),eye(numRuns));
        dLdtheta(indx,s) = -P(s)/2*(traceABtrans(C{indx},Z{s})-traceABtrans(C{indx},B))*exp(runParam);
    end;
    
    % Calculate the overall Likelihood for this subject
    ldet  = -2* sum(log(diag(chol(iV))));        % Safe computation of the log determinant (V) Thanks to code from D. lu
    L(s)     = -P(s)/2*(ldet)-0.5*traceABtrans(iVr,YY{s});
    if (~isempty(X) && ~isempty(X{s})) % Correct for ReML estimates
        L(s) = L(s) - P(s)*sum(log(diag(chol(X{s}'*iV*X{s}))));  % - P/2 log(det(X'V^-1*X));
    end;
end

% Sum over participants
negLogLike = -sum(L);
dnldtheta   = -sum(dLdtheta,2);
