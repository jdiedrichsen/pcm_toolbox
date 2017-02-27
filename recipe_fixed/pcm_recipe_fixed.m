function varargout=pcm_recipe_fixed
% Example of a fit to two fixed pcm-models
% The example of the models and data is taken from Ejaz et al. (2015).
% Nature Neuroscience: "Hand usage predicts.." 
% We fit a null model (all finger patterns are equally far away from each other) 
% and two fixed models (Muscle and Natural stats). 
% Finally, we fit an arbitrary combination
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
% Patterns be independent - i.e. all finger pairs are equally far away from
% each other 
M{1}.type       = 'component';
M{1}.numGparams = 1;
M{1}.Gc         = eye(5);
M{1}.name       = 'null'; 

% Model 2: Muscle model: derived from covariance structure of muscle
% activity during single finger movements 
M{2}.type       = 'component';
M{2}.numGparams = 1;
M{2}.Gc         = Model(1).G_cent;
M{2}.name       = 'muscle'; 

% Model 3: Natural statistics model: derived from covariance structure of
% natual movements 
M{3}.type       = 'component';
M{3}.numGparams = 1;
M{3}.Gc         = Model(2).G_cent;
M{3}.name       = 'usage'; 

% Model 4: Additive mixture between muscle and natural stats model
M{4}.type       = 'component';
M{4}.numGparams = 2;
M{4}.Gc(:,:,1)  = Model(1).G_cent;
M{4}.Gc(:,:,2)  = Model(2).G_cent;
M{4}.name       = 'muscle + usage'; 

% Model 5: Average of all other subjects (Noise ceiling)
M{5}.type       = 'noiseceiling';
M{5}.numGparams = 0;
M{5}.name       = 'noiseceiling'; 

%M{6}.type       = 'freechol'; 
%M{6}.numCond    = 5;
% M{6}.name       = 'noiseceiling'; 
% M{6}           = pcm_prepFreeModel(M{6}); 


% Fit the models on the group level 
[Tgroup,M] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);

% Fit the models through cross-subject crossvalidation
[Tcross,M] = pcm_fitModelCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',Tgroup,'fitScale',1);

% Provide a plot of the crossvalidated likelihoods 
T = pcm_plotModelLikelihood(Tcross,M,'upperceil',Tgroup.likelihood(:,5)); 

varargout={T,M}; 
