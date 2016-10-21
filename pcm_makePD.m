function Gpd = pcm_makePD(G,th)
% function Gpd = pcm_makePD(G,th)
% Enforces that G is semi-positive definite by removing possible negative
% eigenvectors from the square matrix 
% INPUT: 
%       G:      estimated KxK second momement matrix 
%       th:     threshold
% OUTPUT:
%       Gpd :   semi-positive definite version of G   

if nargin<2
    th = eps;
end
    
% Make the G-estimate postive definite: 
G = (G+G')/2;        % Symmetrize 
[V,lam_G] = eig(full(G));
dS    = diag(lam_G);
idx   = dS>th;
Gpd     = V(:,idx)*lam_G(idx,idx)*V(:,idx)'; 
