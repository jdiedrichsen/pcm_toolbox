function varargout=pcm_recipe_1(runEffect)
% Fit a fixed pcm-model
% For run Effect fixed, it does not matter if the
% whether its the centered or non-centered G matrix
% If the run effect is random, it does matter, as the model needs
% to also explain some part of the common activity across
% conditions
load recipe1_data.mat
% specify the model structure:
% Model 1: Muscle model
M(1).type       = 'fixed';
M(1).numGparams = 0;
M(1).Gc         = Model(1).G_cent;
% Model2: Muscle + Natural statistics model
M(2).type       = 'fixed';
M(2).numGparams = 0;
M(2).Gc         = Model(2).G_cent;
% Model 3: Additive mixture between muscle and natural stats model
M(3).type       = 'component';
M(3).numGparams = 2;
M(3).Gc(:,:,1)  = Model(1).G_cent;
M(3).Gc(:,:,2)  = Model(2).G_cent;
% "Model" 4: average of all other subjects (Noise ceiling)
M(4).type       = 'noiseceiling';
M(4).numGparams = 0;

T= pcm_fitModelCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'verbose',0);

% Plot results - maybe also afd the fitModelIndivid
