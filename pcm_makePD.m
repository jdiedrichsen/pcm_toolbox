function Gpd = pcm_makePD(G)
% function Gpd = pcm_makePD(G)
% Enforces that G is semi-positive definite by removing possible negative
% eigenvectors from the square matrix 
% INPUT: 
%       G:      estimated KxK second momement matrix 
% OUTPUT:
%       Gpd :   semi-positive definite version of G   

% Make the G-estimate postive definite: 
G = (G+G')/2;        % Symmetrize 
[V,lam_G] = eig(full(G));
dS    = diag(lam_G);
idx   = dS>eps;
Gpd     = V(:,idx)*lam_G(idx,idx)*V(:,idx)'; 
