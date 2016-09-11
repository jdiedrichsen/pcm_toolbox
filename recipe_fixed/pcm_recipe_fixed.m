function varargout=pcm_recipe_fixed
% Example of a fit to two fixed pcm-models
% The example of the models and data is taken from Ejaz et al. (2015).
% Nature Neuroscience: "Hand usage predicts.." 
% We fit a null model (all finger patterns are independent) 
% two ficed models (Muscle and Natural stats) and the arbtrary combination
% of the two fixed models (component model). The last model is the noise
% ceiling. 
% The log-likelihoods of the models are presented scaled to the null model
% (0) and the upper noise ceiling (1) 

load data_recipe_fixed.mat
% Treat the run effect as random or fixed? 
% For run effect fixed, it does not matter if the
% whether its the centered or non-centered G matrix
% If the run effect is random, it does matter, as the model needs
% to also explain some part of the common activity across
runEffect  = 'fixed'; 

% Model 1: Null model for baseline: here we use a model which has all finger 
% Patterns be independent. 
M(1).type       = 'fixed';
M(1).numGparams = 0;
M(1).Gc         = eye(5);
M(1).name       = 'null'; 

% Model 2: Muscle model: derived from covariance structure of muscle
% activity during single finger movements 
M(2).type       = 'fixed';
M(2).numGparams = 0;
M(2).Gc         = Model(1).G_cent;
M(2).name       = 'muscle'; 

% Model 3: Natural statistics model: derived from covariance structure of
% natual movements 
M(3).type       = 'fixed';
M(3).numGparams = 0;
M(3).Gc         = Model(2).G_cent;
M(3).name       = 'usage'; 

% Model 4: Additive mixture between muscle and natural stats model
M(4).type       = 'component';
M(4).numGparams = 2;
M(4).Gc(:,:,1)  = Model(1).G_cent;
M(4).Gc(:,:,2)  = Model(2).G_cent;
M(4).name       = 'muscle + usage'; 

% Model 5: average of all other subjects (Noise ceiling)
M(5).type       = 'noiseceiling';
M(5).numGparams = 0;

% Fit the models 
T = pcm_fitModelCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'verbose',0);

% Provide a plot of the crossvalidated likelihoods 
T = pcm_plotModelLikelihood(T,M); 

varargout={T}; 
