% Testing script for the finger example 
load data_recipe_finger7T.mat
M{1}.type       = 'component';
M{1}.numGparams = 1;
M{1}.Gc         = Model(1).G;
M{1}.name       = 'muscle';

M{2}.type       = 'component';
M{2}.numGparams = 1;
M{2}.Gc         = Model(2).G;
M{2}.name       = 'muscle';

M{3}.type       = 'component';
M{3}.numGparams = 2;
M{3}.Gc(:,:,1)  = Model(1).G;
M{3}.Gc(:,:,2)  = Model(2).G;
M{3}.name       = 'muscle+nat';

[Z,B,X,YY,Ss,N,P,G_hat,noise0,run0]=pcm_setUpFit(Y,partVec,condVec,'runEffect','fixed');
theta0={zeros(15,1),zeros(15,1),zeros(16,1)};
[T1,theta1]=pcm_fitModelGroup(Y,M{3},partVec,condVec,'fitScale',1,'runEffect','fixed');
[T2,theta2]=pcm_fitModelGroupCrossval(Y,M{3},partVec,condVec,'fitScale',1,'runEffect','fixed','groupFit',theta1);
keyboard; 