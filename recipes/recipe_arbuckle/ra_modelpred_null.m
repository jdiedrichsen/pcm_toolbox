function [G,dGdtheta] = ra_modelpred_null(theta)
% Predicts distaces and G-matrix from the 18 parameters of the simple
% scaling model and returns the derivative in respect to the parameters
G        = ones(20)*exp(theta(1))+eye(20); 
dGdtheta = ones(20)*exp(theta(1));