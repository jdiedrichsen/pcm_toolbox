function [Z,B,X,YY,Ss,N,P,G_hat,noise0,run0]=pcm_setUpFit(Y,partitionVec,conditionVec,varargin);
% function [YY,Z,B,X,noise0,run0]=pcm_setUpFit(Y,conditionVec,partitionVec,runEffect);
%   Function to set up fitting of PCM model to a series of subjects
%   INPUT: 
%       Y: Cell array of data: Y{numSubj}(numObs x #numVoxels)
%       partitionVec: {#Subjects} Cell array with partition assignment vector
%                   for each subject. Rows of partitionVec{subj} define
%                   partition assignment of rows of Y{subj}.
%                   Commonly these are the scanning run #s for beta
%                   regressors.
%                   If a single vector is provided, it is assumed to be the
%                   same for all subjects 
%
%       conditionVec: {#Subjects} Cell array with condition assignment vector
%                   for each subject. Rows of conditionVec{subj} define
%                   condition assignment of rows of Y{subj}.
%                   If a single vector is provided, it is assumed to me the
%                   same for all subjects.
%                   If the (elements of) conditionVec are matrices, it is
%                   assumed to be the design matrix Z, allowing the
%                   specification individualized models. 
%  OPTIONS: 
%       'S':        Specific assumed noise structure - usually inv(XX'*XX),
%                   where XX is the first-level design matrix used to
%                   estimate the activation estimates. Either structure,
%                   matrix, or cell array 
%
%      'runEffect':   How to deal with effects that may be specific to different
%                   imaging runs:
%                  'random': Models variance of the run effect for each subject
%                            as a seperate random effects parameter.
%                  'fixed': Consider run effect a fixed effect, will be removed
%                            implicitly using ReML (default). 
%                  'none': No modeling of the run effect (not recommended
%                           for real fMRI data)
%  OUTPUT: 
%       Z:          Design matrix for random effect of interest 
%       B:          Design matrix for run-effects (random) 
%       X:          Design matrix for run-effect (fixed)
%       YY:         Sufficient stats on data 
%       Ss:          Structure with noise matrix
%       N:          Number of Observations per subject (numObs x 1) vector 
%       P:          Number of voxels per subject (numSubj x 1) vector 
%       G_hat:      Cross-validated estimate of second moment matrix for
%                   each subject 
%       noise0:     Starting value for the noise parameter for each subject       
%       run0:       For runEffect = 'random', starting value of run effect
%                   parameter 
% Joern Diedrichsen (joern.diedrichsen@googlemail.com) 

runEffect = 'fixed'; 
S         = []; 
pcm_vararginoptions(varargin,{'runEffect','S'});

% Determine number of subjects 
numSubj = numel(Y); 

% Loop over all subjects 
for s = 1:numSubj

    % Set up the main matrices
    [N(s,1),P(s,1)] = size(Y{s});   
    
    % Caluculate sufficient statistics 
    YY{s}  = (Y{s} * Y{s}');

    % If condition and partition Vectors are not cells, assume they are the same 
    if (iscell(conditionVec)) 
        cV = conditionVec{s}; 
        pV = partitionVec{s}; 
    else 
        cV = conditionVec; 
        pV = partitionVec; 
    end; 
    
    % Check if conditionVec is condition or design matrix
    if size(cV,2)==1;
        Z{s}   = pcm_indicatorMatrix('identity_p',cV);
    else
        Z{s} = cV;
    end;
    numReg= size(Z{s},2);
    
    % Depending on the way of dealing with the run effect, set up data
    % and determine initial values for run and noise effects 
    switch (runEffect)
        case 'random'
            B{s}   = pcm_indicatorMatrix('identity_p',pV);
            X{s}   = zeros(N(s),0);
            numPart=size(B{s},2);
            run0(s,1)=log(sum(sum((pinv(B{s})*Y{s}).^2))./(numPart*P(s))); 
            RX = eye(N(s))-B{s}*pinv(B{s}); 
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},pV,cV);
        case 'fixed'
            B{s}  =  zeros(N(s),0);
            run0  =  []; 
            X{s}  =  pcm_indicatorMatrix('identity_p',pV);
            numPart=size(X{s},2);
            RX = eye(N(s))-X{s}*pinv(X{s}); 
            G_hat(:,:,s) = pcm_estGCrossval(RX*Y{s},pV,cV);
        case 'none' 
            B{s}    =  zeros(N(s),0);
            run0    =  []; 
            X{s}    =  zeros(N(s),0); 
            numPart = 0;  % no mean accounted for  
            RX      = eye(N(s)); 
            G_hat(:,:,s) = pcm_estGCrossval(RX*Y{s},pV,cV);
    end;
    
    % Estimate noise variance 
    numReg   = size(Z{s},2); 
    RZ           = eye(N(s))-Z{s}*pinv(Z{s}); 
    noise0(s,1)  = sum(sum((RZ*RX*Y{s}).^2))/(P(s)*(N(s)-numReg-numPart));
    if (noise0(s)<=0) 
        error('Too many model factors to estimate noise variance. Consider removing terms or setting runEffect to ''none'''); 
    end; 
    noise0(s)=log(noise0(s)); 
    
    % Set up noise covariance structure 
    if (~isempty(S))
        if (~isstruct(S)) 
            if (iscell(S))
                Ss(s).S=S{s}; % Per Subject 
            else 
                Ss(s).S=S;   % Same for all subjects 
            end;
            Ss(s).invS = inv(Ss(s).S); 
        end; 
        noise0(s)=noise0(s)./mean(trace(Ss(s).S)); 
    else 
        Ss=[]; 
    end; 
end;