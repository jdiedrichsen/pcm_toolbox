function [G, dGdtheta, Gcomponent] = pcm_combineModels(theta,Models,varargin)
% Linear component modeling of nonlinear G matrices; G = sum(Gc), dG = cat(dGc).
% 
% function [G, dGdtheta] = pcm_modelpred(theta,Models,varargin)
% 
% Inputs:
%   theta:  parameters for all models used
%   Models: model structures with fields;
%       modelpred:  function handle for a model which also gives G and dG
%       numGparams: number of parameters for the model
%       params:     additional parameters for the model function
% 
% 
% 
% Outputs:
%   G:          second moment matrix
%   dGdtheta:   derivative of G with respect to theta
%   Gcomponent: component G
% 
% 
% 
% a-yokoi (2016)


theta_ = theta;
Nmodel = length(Models);

G = [];
dGdtheta = [];
Gcomponent = [];

% Concatenate M and O then compute G
for m=1:Nmodel
    Nparam(m) = Models(m).numGparams;
    theta_m = theta_(1:Nparam(m));
    theta_(1:Nparam(m)) = [];
    
   [Gc, dG] = Models(m).modelpred(theta_m,Models(m).params{:}); 
   
   if isempty(G);
       G = zeros(size(Gc));
   end
   
   G = G+Gc;
   if ~isempty(dG)
       dGdtheta = cat(3,dGdtheta,dG);
   end
   if nargout>2
       Gcomponent = cat(3,Gcomponent,Gc);
   end
end

end