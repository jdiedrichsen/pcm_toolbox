function [G,dGdtheta] = pcm_calculateG(M,theta)
% function [G,dGdtheta] = pcm_calculateG(M,theta)
% This function calculates the predicted second moment matrix (G) and the
% derivate of the second moment matrix in respect to the parameters theta. 
% INPUT: 
%       M:        Model structure 
%       theta:    Vector of parameters 
% OUTPUT: 
%       G:        Second moment matrix 
%       dGdtheta: Matrix derivatives in respect to parameters 
% Joern Diedrichsen, 2016 

if (~isstruct(M))
    G=M;   % Fixed, parameterless G
    dGdtheta = [];
else
    if length(theta)~=M.numGparams
        error('lenth of column-vector theta should be equal to number of G parameters'); 
    end; 
    if (~isfield(M,'type'))
        error('M should have a type of fixed / component / feature / nonlinear');
    end;
    switch (M.type)
        case {'fixed','freedirect'}
            G        = mean(M.Gc,3);
            dGdtheta =[]; 
        case 'component'
            dGdtheta=bsxfun(@times,M.Gc,permute(exp(theta),[3 2 1]));
            G = sum(dGdtheta,3); 
        case {'feature'}
            A = bsxfun(@times,M.Ac,permute(theta,[3 2 1]));
            A = sum(A,3); 
            G = A*A'; 
            for i=1:M.numGparams
                dA = M.Ac(:,:,i)*A';  
                dGdtheta(:,:,i) =  dA + dA';     
            end; 
        case 'freechol' 
            A         = zeros(M.numCond); 
            A(M.indx) = theta; 
            G = A*A';
            dGdtheta = zeros(M.numCond,M.numCond,M.numGparams); 
            for i=1:M.numGparams
                dGdtheta(M.row(i),:,i)  = A(:,M.col(i))'; 
                dGdtheta(:,:,i) = dGdtheta(:,:,i) + dGdtheta(:,:,i)';     
            end;             
        case 'nonlinear'
            [G,dGdtheta]=M.modelpred(theta(1:M.numGparams));
    end;
end;
