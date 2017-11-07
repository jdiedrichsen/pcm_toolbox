function A = pcm_diagonalize(G,th); 
% function A = pcm_diagonalize(G,th); 
% safe decomposition of G into A * A' 
% Removes eigenvectors below a value of th (default = eps)

if nargin<2
    th = eps;
end

% Diagnoalize G up to a threshold 
G = (G+G')/2;        % Symmetrize 
[V,lam_G] = eig(full(G));
dS        = diag(lam_G);
idx       = dS>th;
A       = V(:,idx)*sqrt(lam_G(idx,idx)); 