function [Y,l] =pcm_classicalMDS(G,varargin)
% function [Y,l] =pcm_classicalMDS(G,varargin)
% Classical MDS directly from the second moment matrix 
% This routine can deal appropriately with crossvalidateed 
% second moment matrices - i.e. matrices that are not positive definite. 
% 
% INPUT: 
%   G:  KxK matrix form of the second moment matrix estimate 
%   'contrast',Z 
%       Z is a contrast matrix that specifies what type of contrast is
%       optimized in the resultant representation. If the contrast involves K
%       independent conditions, the maximal dimensionality of the
%       representation is K-1. 
%  OUTPUT: 
%    Y:   (K,Q) coordinates of the K conditions on the Q dimensions of the
%               representational space.  
%    l:   (Qx1) eigenvalues for the importance of the different vectors. 

% NOTES: if you want to know the direction that a particular feature vector
% has in your representational space, you simply need to project the
% contrast vector onto the plotted space defined by Y. 
contrast = []; 
pcm_vararginoptions(varargin,{'contrast'});

% Get the eigenvalues from the full G-matrix and restruct the space, 
% retaining the negative eigenvalues. 
G       = (G+G')/2; 
[V,L]   = eig(G);                     % This makes a Eucledian distance to covariance matrix 
[l,i]   = sort(diag(L),1,'descend');           % Sort the eigenvalues
V       = V(:,i);
Y       = bsxfun(@times,V,sqrt(l'));

% Now rotate the projection to maximize discriminant for a certain contrast 
% See Notes/SVD.pages 
if (~isempty(contrast)) 
    numVec = rank(contrast); 
    H = contrast*pinv(contrast);  % Projection matrix 
    [V,L]=eig(conj(Y)'*H'*H*Y); 
    [l,i]   = sort(real(diag(L)),1,'descend');           % Sort the eigenvalues
    V       = V(:,i); 
    Y       = Y*V; 
else 
    V=eye(size(Y,2));
end; 

% 
indx = find(l<eps); 
Y(:,indx)=0;                % Kill superflous dimensions 
Y=real(Y); 
V(:,indx)=0;
V=real(V);
 