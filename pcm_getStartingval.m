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
    case {'fixed','noiseceiling'}
        theta0=[]; 
    case 'component' 
        for i=1:size(M.Gc,3)            % Use normal regression against crossvalidated G-matrix to estimate starting parameters 
            g = M.Gc(:,:,i); 
            X(:,i)= g(:); 
        end; 
        h0 = pinv(X)*G_hat(:); 
        h0(h0<10e-4)=10e-4;        % Ensure positivity of the parameters 
        theta0 = log(h0); 
        if (M.numGparams==0)
            theta0=[]; 
        end; 
    case {'squareroot','feature'} 
        for i=1:size(M.Ac,3); 
            g = M.Ac(:,:,i)*M.Ac(:,:,i)'; 
            X(:,i)= g(:); 
        end; 
        h0 = pinv(X)*G_hat(:); % Use normal regression to estimate parameters 
        h0(h0<10e-4)=10e-4;        % Ensure positivity of the parameters 
        theta0 = log(h0); 
    case 'freechol'      % Totally free model using cholesky decomposition 
        G_hat  = pcm_makePD(G_hat); 
        G = (G_hat+G_hat')/2;        % Symmetrize 
        [V,lam_G] = eig(full(G));
        dS    = diag(lam_G);
        dS(dS<eps) = eps; 
        Gpd     = V*diag(dS)*V'; 
        A      = chol(Gpd)'; 
        theta0 = A(M.indx); 
    case 'nonlinear' 
        error('cannot provide starting values for nonlinear models'); 
end; 
        