load data_recipe_finger7T.mat
M{1}.type       = 'component';
M{1}.numGparams = 1;
M{1}.Gc         = Model(1).G;
M{1}.name       = 'muscle';

YY = Y{1}*Y{1}';
n_channels = size(Y{1},2); 
Z = indicatorMatrix('identity',condVec{1});
X = indicatorMatrix('identity',partVec{1});
theta = [-0.5;0.5;0.1]; 

[lik,dL,dLL] = pcm_likelihoodIndivid(theta,YY,M{1},Z,[],n_channels,'fitScale',0,'runEffect',X);