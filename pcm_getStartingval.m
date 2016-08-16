function theta0 = pcm_getStartingval(M,G_hat)
% function theta = pcm_getStartingvel(M,G_hat)
% Provides startingvalues for model fitting from the crossvalidated G-matrix
% estimate provided as an input to this function. 
% 
% INPUT: 
%       G_hat:      estimated KxK second momement matrix
% OUTPUT:  
%       theta0:     vector of model parameters 

switch (M.type) 
    case 'fixed'
        theta0=[]; 
    case 'component' 
        for i=1:size(M.Gc,3)
            g = M.Gc(:,:,i); 
            X(:,i)= g(:); 
        end; 
        h0 = pinv(X)*G_hat(:); % Use normal regression to estimate parameters 
        h0(h0<10e-6)=10e-6;        % Ensure positivity of the parameters 
        theta0 = log(h0); 
    case 'squareroot' 
        error('not implemented yet'); 
    case 'nonlinear' 
        error('cannot provide starting values for nonlinear models'); 
    case 'noiseceiling'
        if ~(isempty(M.theta))
            theta0 = M.theta0;
        else
            theta0 = [];
        end
end; 
        
