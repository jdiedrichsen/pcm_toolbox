function varargout=pcm_recipe_anova
% Example of a pcm-model to fit a multivariate anova example 
% The data is taken from Kornysheva et al. (2013). 
% Elife: "Hand usage predicts.." 
% Subjects were instructed to perform 9 temporal-spatial sequence of finger
% movements 
% There were 3 different temporal rythms and 3 different finger orders 
% The models are now whether the region encodes the temporal rythm, the
% finger order or the interaction of the two 
% Here we look at three different ROIs, each of which has a very different
% signature 

load data_recipe_anova.mat
% Here we treat the run effect as a fixed effect, 
% Effectively subtracting it out. Therefore we do not have to model 
% the average activity at all. 
runEffect  = 'fixed'; 

tem = [1;1;1;2;2;2;3;3;3]; % Temporal structure of the sequence  
ord = [1;2;3;1;2;3;1;2;3]; % Order of the sequence 
T = pcm_indicatorMatrix('identity',tem); 
O = pcm_indicatorMatrix('identity',ord); 

% Model 0: Null model for baseline: here we use a model the null of model
% of no encoding as a baseline 
M{1}.type       = 'fixed';
M{1}.numGparams = 0;
M{1}.Gc         = zeros(9);
M{1}.name       = 'null'; 

% Model 2: Rythm main effect 
M{2}.type       = 'component';
M{2}.numGparams = 1;
M{2}.Gc         = T*T';
M{2}.name       = 'T'; 

% Model 3: Order Main effect 
M{3}.type       = 'component';
M{3}.numGparams = 1;
M{3}.Gc         = O*O';
M{3}.name       = 'O'; 

% Model 4: Interaction Main effect 
M{4}.type       = 'component';
M{4}.numGparams = 1;
M{4}.Gc         = eye(9);
M{4}.name       = 'I'; 

% Model 5: Order + Movement, full model 
M{5}.type       = 'component';
M{5}.numGparams = 2;
M{5}.Gc(:,:,1)  = T*T';
M{5}.Gc(:,:,2)  = O*O';
M{5}.name       = 'T+O'; 

% Model 6: Order + Movement, full model 
M{6}.type       = 'component';
M{6}.numGparams = 3;
M{6}.Gc(:,:,1)  = T*T';
M{6}.Gc(:,:,2)  = O*O';
M{6}.Gc(:,:,3)  = eye(9);
M{6}.name       = 'T+O+I'; 

% Model 7: Average of all other subjects (approximate Noise ceiling)
M{7}.type       = 'noiseceiling';
M{7}.numGparams = 0;
M{7}.name       = 'noiseceiling'; 

% Fit the models on the group level 
[Tgroup,M] = pcm_fitModelGroup(Y{2},M,partVec,condVec,'runEffect',runEffect);

% Fit the models through cross-subject crossvalidation
% [Tcross,M] = pcm_fitModelCrossval(Y{1},M,partVec,condVec,'runEffect',runEffect,'groupFit',Tgroup);

% Provide a plot of the crossvalidated likelihoods 
T = pcm_plotModelLikelihood(Tgroup,M,'normalize',0); 

varargout={T,M}; 
