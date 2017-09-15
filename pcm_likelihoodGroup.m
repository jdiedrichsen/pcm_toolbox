function [negLogLike,dnl,d2nl,LogLike] = pcm_likelihoodGroup(theta,YY,M,Z,X,P,varargin);
% function [negLogLike,dnl,d2nl] = pcm_likelihoodGroup(theta,YY,M,Z,X,P,varargin);
% Returns negative log likelihood of the data from a group of subjects under the
% representational model specified by M. The function uses common model
% parameters (theta) for all subjects , i.e the structure of the second-moment
% matrix is supposed to be same across all subjects.
% To link the group prediction to individual second-moment matrices, a
% number of individual-level nuisance parameters are introduced:
%   s:        Scaling of the data
%   sigma2:   Noise variability
% Additionally, if 'runEffect' is set, the effect of each imaging run is
%   also modelled as a random effect, giving the third individual-subject
%   parameter
% INPUT:
%      theta:   Vector of model parameters + subject-specific nuisance parameters
%               1...numGParams           : model parameters
%               numGParams+1... +numSubj : log(noise variance)
%               numGParams+numSubj+1...  : log(scaling parameters) - used if fitScale
%               numGParams+numSubj*2+1...: log(run effect variance)- used if runEffect is provided
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
%      'fitScale'   Introduce additional scaling parameter for each
%                   participant? - default is true
%      'S',S:       Optional structured noise matrix - default is the identity matrix
%                   S should be a structure-array (number of subjects) with two fields
%                   S.S: noise structure matrix
%                   S.invS: Inverse of noise structure matrix
%      'verbose',0/1: Set verbose mode -
% OUTPUT:
%      negLogLike:  Negative Log likelihood of all subject's data
%                   We use the negative log liklihood to be able to plug the function into
%                   minimize or other optimisation routines.
%      dnl:         Derivative of the negative log-likelihood
%      d2nl:        expected second derivative of the negative log-likelihood
%      LogLike:     Log likelihood (not inverted) for individual subjects
%
%   Joern Diedrichsen, 6/2016, joern.diedrichsen@googlemail.com

OPT.S         = [];
OPT.verbose   = 0;
OPT.runEffect = [];
OPT.fitScale  = 1;
OPT = pcm_getUserOptions(varargin,OPT,{'S','verbose','runEffect','fitScale'});

% Get G-matrix and derivative of G-matrix in respect to parameters
if (isstruct(M))
    [G,dGdtheta] = pcm_calculateG(M,theta(1:M.numGparams));
else
    G=M;
    M=[];
    M.numGparams=0;
end;

% Preallocate arrays
numSubj = length(YY);
numParams = 1 + (~isempty(OPT.runEffect{1})) + (OPT.fitScale>0); % Number of parameters per subject 
numParams = numSubj * numParams + M.numGparams;                  % Total number of subjects 
dLdtheta=zeros(numParams,numSubj); 
d2L=zeros(numParams,numParams,numSubj); 

for s=1:numSubj    
    % Get parameter and sizes
    noiseParam = theta(M.numGparams+s);             % Subject Noise Parameter
    if (OPT.fitScale)
        scaleParam = theta(M.numGparams+numSubj+s);   % Subject Noise Parameter
    else
        scaleParam = 0;
    end;
    Gs = G*exp(scaleParam);         % Scale the subject G matrix up by individual scale parameter
    [N,K]  = size(Z{s});
    
    % If Run effect is to ne modelled as a random effect - add to G and
    % design matrix
    if (~isempty(OPT.runEffect{s}))
        numRuns = size(OPT.runEffect{s},2);
        runParam = theta(M.numGparams+numSubj*(1+fitScale)+s);    % Subject run effect parameter
        Gs = pcm_blockdiag(Gs,eye(numRuns)*exp(runParam));  % Include run effect in G
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
    if (OPT.verbose)
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
    if (~isempty(OPT.S));
        sS = OPT.S(s).S; % noise covariance
        iS = OPT.S(s).invS; % inverse noise covariance
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
    
    % Compute the (restricted) likelihood for this Subject
    ldet  = -2* sum(log(diag(chol(iV))));        % Safe computation of the log determinant (V) Thanks to code from D. lu
    LogLike(s)     = -P(s)/2*(ldet)-0.5*traceABtrans(iVr,YY{s});
    if (~isempty(X) && ~isempty(X{s})) % Correct for ReML estimates
        LogLike(s) = LogLike(s) - P(s)*sum(log(diag(chol(X{s}'*iV*X{s}))));  % - P/2 log(det(X'V^-1*X));
    end;
    
    % Calculate the first derivative
    if (nargout>1)
        A     = iVr*Z{s};
        B     = YY{s}*iVr;
        
        % Get the derivatives for all the parameters
        indx = [1:M.numGparams]; 
        for i = indx 
            iVdV{i} = A*pcm_blockdiag(dGdtheta(:,:,i),zeros(numRuns))*Z{s}'*exp(scaleParam);
            dLdtheta(i,s) = -P(s)/2*trace(iVdV{i})+1/2*traceABtrans(iVdV{i},B);
        end
        
        % Get the derivatives for the Noise parameters
        i = M.numGparams+1; 
        indx(i) = M.numGparams+s;  % Which number parameter is it?
        if (isempty(OPT.S))
            dVdtheta  = eye(N)*exp(noiseParam);
        else
            dVdtheta  = sS*exp(noiseParam);
        end;
        iVdV{indx(i)}     = iVr*dVdtheta;
        dLdtheta(indx(i),s) = -P(s)/2*trace(iVdV{indx(i)})+1/2*traceABtrans(iVdV{indx(i)},B);
        
        % Get the derivatives for the scaling parameters
        if (OPT.fitScale)
            i = i+1; 
            indx(i) = M.numGparams+numSubj+s;    % Which number parameter is it?
            iVdV{indx(i)}          = A*pcm_blockdiag(G,zeros(numRuns))*Z{s}'*exp(scaleParam);
            dLdtheta(indx(i),s)    = -P(s)/2*trace(iVdV{indx(i)})+1/2*traceABtrans(iVdV{indx(i)},B);
        end;
        
        % Get the derivatives for the block parameter
        if (~isempty(OPT.runEffect{s}))
            i = i+1; 
            indx(i) = M.numGparams+numSubj*(1+OPT.fitScale)+s;  % Which number parameter is it?
            iVdV{indx(i)}       = A*pcm_blockdiag(zeros(K),eye(numRuns))*Z{s}'*exp(runParam);
            dLdtheta(indx(i),s) = -P(s)/2*trace(iVdV{indx(i)})+1/2*traceABtrans(iVdV{indx(i)},B);
        end;
    end;
    
    % Determine second derivative for non-zero entries 
    if (nargout>2)
        for i=1:length(indx)
            for j=i:length(indx)
                d2L(indx(i),indx(j),s)=-P(s)/2*traceABtrans(iVdV{indx(i)},iVdV{indx(j)});
                d2L(indx(j),indx(i),s)=d2L(indx(i),indx(j),s); 
            end;
        end;
    end;
end

% Sum over participants
negLogLike = -sum(LogLike);
dnl        = -sum(dLdtheta,2);
d2nl       = -sum(d2L,3);
