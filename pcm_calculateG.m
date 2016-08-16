function [G,dGdtheta] = pcm_calculateG(M,theta,G_hat)
% function [G,dGdtheta] = pcm_calculateG(M,theta,G_hat)
% This function calculates the predicted second moment matrix (G) and the
% derivate of the second moment matrix in respect to the parameters theta. 
% INPUT: 
%       M:        Model structure 
%       theta:    Vector of parameters 
%       G_hat:    Crossvalidated second moment matrix for subject
% OUTPUT: 
%       G:        Second moment matrix 
%       dGdtheta: Matrix derivatives in respect to parameters 
% Joern Diedrichsen, 2016 

if (~isstruct(M))
    G=M;   % Fixed, parameterless G
    dGdtheta = [];
else
    if length(theta)~=M.numGparams
        error('lenth of vector theta should be equal to number of G parameters'); 
    end; 
    if (~isfield(M,'type'))
        error('M should have a type of fixed / component / squareroot / nonlinear');
    end;
    switch (M.type)
        case 'fixed'
            G=M.Gc;
            dGdtheta =[]; 
        case 'component'
            dGdtheta=bsxfun(@times,M.Gc,permute(exp(theta),[3 2 1]));
            G = sum(dGdtheta,3); 
        case 'squareroot'
            A = bsxfun(@times,M.Ac,permute(theta,[3 2 1]));
            A = sum(A,3); 
            G = A*A'; 
            for i=1:M.numGparams
                dGdtheta(:,:,i) = M.Ac(:,:,1)*A' + A*M.Ac(:,:,1)';     
            end;                 
        case 'nonlinear'
            [G,dGdtheta]=M.modelpred(theta(1:M.numGparams));
        case 'noiseceiling'
            G = M.modelpred(G_hat);
            dGdtheta = [];
    end;
end;
