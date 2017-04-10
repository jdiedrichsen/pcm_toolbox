function varargout=pcm_recipe_finger
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

load data_recipe_finger7T.mat

% ----------------------------------------------------------------
% Data visualisation 
% 1. visualize the representational stucture: Get the crossvalidated
%    G-matrix for each data set and plot the average of this
for s=1:length(Y)
    G_hat(:,:,s)=pcm_estGCrossval(Y{s},partVec{s},condVec{s}); 
end; 
Gm = mean(G_hat,3); % Mean estimate  
subplot(2,3,1); 
H = eye(5)-ones(5)/5; 
imagesc(H*Gm*H'); 

% ----------------------------------------------------------------
% 2. visualize the the average G-matrix using an MDS plot 
C= pcm_indicatorMatrix('allpairs',[1:5]'); 
[COORD,l]=pcm_classicalMDS(Gm,'contrast',C); 
subplot(2,3,2); 
scatterplot(COORD(:,1),COORD(:,2),'label',[1:5]); 
axis equal; 

% ----------------------------------------------------------------
% 2. Visualize the the predicted G-matrix from the different models 
subplot(2,3,4); 
imagesc(Model(1).G_cent); 
title('Muscle'); 

subplot(2,3,5); 
imagesc(Model(2).G_cent); 
title('Naturalstats'); 

% ----------------------------------------------------------------
% Now build the models 
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

% Model 5: Free model as Noise ceiling
M{5}.type       = 'freechol'; 
M{5}.numCond    = 5;
M{5}.name       = 'noiseceiling'; 
M{5}           = pcm_prepFreeModel(M{5}); 

% Treat the run effect as random or fixed? 
% We are using a fixed run effect here, as we are not interested in the
% activity relative the the baseline (rest) - so as in RSA, we simply
% subtract out the mean patttern across all conditions. 
runEffect  = 'fixed'; 

% Fit the models on the group level 
[Tgroup,theta] = pcm_fitModelGroup(Y,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);

% Fit the models through cross-subject crossvalidation
[Tcross,thetaCr] = pcm_fitModelGroupCrossval(Y,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta,'fitScale',1);

% Provide a plot of the crossvalidated likelihoods 
subplot(2,3,[3 6]); 
T = pcm_plotModelLikelihood(Tcross,M,'upperceil',Tgroup.likelihood(:,5)); 

varargout={T,M}; 
