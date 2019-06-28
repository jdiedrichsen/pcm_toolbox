function M = pcm_buildCorrModel(varargin); 
% function M = pcm_buildCorrModel(varargin); 
% Builds different representational models to assess the correlation between
% two sets of patterns. The same items are measured under two different conditions, 
% and each item pattern is predicted to have a
% specific correlation with the one under the the other condition
% For fixed correlation models the
% 'type' should be set to 'nonlinear'. 
%  
% INPUT ARGUMENT (varargin) 
%   'r': either fixed correlation ([-1 - 1]) or 'flexible' 
%   'type': Either 'feature' for 'nonlinear' for specific model formulation 
%   'numItems': How many items measured for each condition 
%   'withinCov': How should the within Condition covariance structure be modelled? 
%        'iid': Same variance + independent for all 
%        'individual': Independnet, but with a individual scalable variances 
%        Matrix: Arbitrary, fixed covariance structure 
%   'condEffect': Model overall condition effect? (Default is yes) 
r     =  0.5;  % fixed correlation or 'flexible'
type  = 'nonlinear'; 
numItems = 6; 
numCond  = 2; 
condEffect = 1;                 
withinCov = 'individual';     % How to model the covariance structure within each condition  'iid' / 'individual'
name      = []; 

pcm_vararginoptions(varargin,{'r','type','numItems','numCond','itemVec',... 
    'condVec','withinCov','condEffect','name'});

M.r = r; 
M.condVec=kron([1:numCond]',ones(numItems,1)); 
M.itemVec=kron(ones(numCond,1),[1:numItems]');

switch (type) 
    case 'feature' 
        M.type         = 'feature';
        if (isempty(name))
            if (isnumeric(r))
                M.name         = sprintf('Corr_feature_%1.3f',r);
            else 
                M.name         = sprintf('Corr_feature_%s',r);
            end 
        else 
            M.name = name; 
        end 
        if (condEffect) 
            for i=1:numCond 
                M.Ac(:,i,i)    = double(M.condVec==i);
            end
        end; 
        if (strcmp(withinCov,'individual'))
            A=zeros(numItems); 
            for i=1:numItems
                indx=find(M.itemVec==i & M.condVec==1); 
                A(indx,indx,i)=1; 
            end
        elseif(strcmp(withinCov,'iid'))
            A = eye(numItems);
        elseif(isnumeric(withinCov))
            A = withinCov; 
        else 
            error('withinCov needs to be either ''iid'',''individual'', or a matrix'); 
        end
        n = size(A,3);      % How many parameters 
        i1=find(M.condVec==1);
        i2=find(M.condVec==2); 
        offset = numCond*condEffect; 
        if (isnumeric(r)) % Specified fixed correlation 
            th2=sqrt(1-r^2);
            th3=r;     % equation for determining thetas
            M.Ac(i1,offset+i1,offset+[1:n]) = A;   % Items for condition 1 
            M.Ac(i2,offset+i2,offset+n+[1:n]) = A*th2;  % Unique sess2 sequence pattterns
            M.Ac(i2,offset+i1,offset+n+[1:n]) = A*th3;  % Shared pattern from condition 1         
        elseif (strcmp(r,'flexible'))  
            M.Ac(i1,offset+i1,offset+[1:n]) = A;   % Items for condition 1 
            M.Ac(i2,offset+i2,offset+n+[1:n]) = A;  % Unique sess2 sequence pattterns
            M.Ac(i2,offset+i1,offset+2*n+[1:n]) = A;  % Shared pattern from condition 1         
        else
            error('Correlation must be either ''flexible'' or a fixed number'); 
        end
        M.numGparams   = size(M.Ac,3);

    case 'nonlinear' 
        M.type         = 'nonlinear';
        if (isempty(name))
            if (isnumeric(r))
                M.name         = sprintf('Corr_nonlin_%1.3f',r);
            else 
                M.name         = sprintf('Corr_nonlin_%s',r);
            end 
        else 
            M.name = name; 
        end 
        if (condEffect)
            for i=1:numCond 
                indx  = find(M.condVec==i); 
                M.Gc(indx,indx,i)    = 1;
            end 
        end
        if (strcmp(withinCov,'individual'))
            A=zeros(numItems); 
            for i=1:numItems
                indx=find(M.itemVec==i & M.condVec==1); 
                A(indx,indx,i)=1; 
            end
        elseif(strcmp(withinCov,'iid'))
            A = eye(numItems);
        elseif(isnumeric(withinCov))
            A = withinCov; 
            if (size(A,3))>1 
                error('arbitrary covariance structures within condition need to be fixed'); 
            end; 
        else 
             error('withinCov needs to be either ''iid'',''individual'', or a matrix'); 
        end; 
        n = size(A,3);      % How many parameters 
        i1=find(M.condVec==1);
        i2=find(M.condVec==2); 
        M.Gc(i1,i1,numCond*condEffect+[1:n]) = A;   % Covariance structure for condition 1 
        M.Gc(i2,i2,numCond*condEffect+n+[1:n]) = A;  % Covariance structure for condition 2 
        
        if (isnumeric(r)) % Specified fixed correlation 
            M.numGparams   = size(M.Gc,3);
        elseif (strcmp(r,'flexible'))  
            M.numGparams   = size(M.Gc,3)+1;
        else
            error('Correlation must be either ''flexible'' or a fixed number'); 
        end
        M.numCond=numCond; 
        M.numItems=numItems; 
        M.numWithinParam = n; 
        M.condEffect = condEffect;
        M.modelpred=@pcm_calculateGnonlinCorr; 
        M.theta0=zeros(M.numGparams,1); 
end; 
