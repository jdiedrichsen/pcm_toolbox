function [theta,G,scaleParam] = pcm_free_startingval(G,varargin)
% function [theta,G,scaleParam] = pcm_modelpred_free_startingval(G,varargin)
% Provides startingvalues from a estimated G-matrix for the optimisation of 
% the liklihood function using pcm_modelLike. 
% The G matrix with a full set of parameters as A*A' 
% Therefore it uses K(K-1)/2+K free parameters. 
% As a default, it reduces the parameter set by one and replaces the first
% parameter by 1 to account for a free scaling parameter 
% INPUT: 
%       G_hat:      estimated KxK second momement matrix (does not have to
%                   be pd
% VARARGIN: 
%       'reduced',1: Reduce the parameters by setting the first parameter
%                    to 1? (default to reduce)
% OUTPUT:  theta:    vector of parameters to be entered in
%                    pcm_modelpred_free
%       G :          positive definite version of G   
%       scaleParam:  This is the scale parameter, such that A*A'*scaleParam=G

% Check if we need to add an additional theta 
reduced = 1; 
pcm_vararginoptions(varargin,{'reduced'}); 

% Make the G-estimate postive definite: 
G = (G+G')/2;        % Symmetrize 
[V,lam_G] = eig(full(G));
dS    = diag(lam_G);
idx   = dS>eps;
G     = V(:,idx)*lam_G(idx,idx)*V(:,idx)'; 

% Get matrix sqaure root and normalise to 1 if neccesary 
A = real(sqrtm(G)); 
if (reduced) 
    scaleParam=A(1,1).^2; 
    A=A./A(1,1); 
else 
    scaleParam=1; 
end; 

% Now vectorize 
n = size(A,1); 
indx=tril(true(n),0);
theta=A(indx);
if (reduced) 
    theta=theta(2:end); 
end; 
