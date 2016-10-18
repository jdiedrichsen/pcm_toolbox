function G=pcm_normMean(G_hat); 
% Means estimated G-matrix after normalizing each 
factor = sqrt(sum(sum(G_hat.^2,1),2));
G = mean(bsxfun(@rdivide,G_hat,factor),3); 