function varargout=pcm_recipe_correlation
% Example for a simple model with two factors:
% Factor A: 1-5 finger
% Factor B: left or right hand
% Example data come from the left motor cortex (see Diedrichsen et al. 2013).
% Note that the strength of encoding for the right (contralateral) fingers is much larger
% than for the left (ipsilateral) fingers.
%
%
load data_recipe_correlation.mat
% Treat the run effect as random or fixed?
% For run effect fixed, it does not matter if the
% If the run effect is random, it does matter, as the model needs
% to also explain some part of the common activity across
runEffect  = 'random';

% Calculate empricial correlation 
for p=1:12 
    Z=pcm_indicatorMatrix('identity',condVec{p}); 
    b = pinv(Z)*Data{p}; 
    x1 = b(1:5,:); 
    x2 = b(6:10,:); 
    x1 = bsxfun(@minus,x1,mean(x1)); 
    x2 = bsxfun(@minus,x2,mean(x2));
    R=corr([x1' x2']); 
    meanR(p) = mean(diag(R,5)); 
end; 
fprintf('mean correlation = %2.2f (+-%2.2f)\n', mean(meanR),std(meanR)/sqrt(12)); 



% Model2: Reduced model like in a chol-decomposition
M{1}.type       = 'feature';
M{1}.numGparams = 5;
M{1}.Ac(:,1,1)  = [ones(5,1);zeros(5,1)]; % Inidividual component: right hand (b)
M{1}.Ac(:,2,2)  = [zeros(5,1);ones(5,1)]; % Common component: left hand (c)
M{1}.Ac(:,3:7,3)  = [eye(5);zeros(5)]; % Inidividual component: right hand (b)
M{1}.Ac(:,3:7,4)  = [zeros(5);eye(5)]; % Common component: left hand (c)
M{1}.Ac(:,8:12,5)  = [zeros(5);eye(5)]; % Common component: right hand (d)
M{1}.name       = 'correlation';
M{1}.theta0=[1 1 1 0 0.4]';

% Fit the models
% [T,M] = pcm_fitModelIndivid(Data,M,partVec,condVec,'runEffect',runEffect);
% T.theta = M{1}.thetaIndiv';
[T,theta_all] = pcm_fitModelIndivid(Data,M,partVec,condVec,'runEffect',runEffect);
T.theta = theta_all{1};
T.var1    = T.theta(:,3).^2;
T.var2    = T.theta(:,4).^2+T.theta(:,5).^2;
T.cov     = T.theta(:,3).*T.theta(:,4);
T.r       = T.cov./sqrt(T.var1.*T.var2);

% Provide a plot of the crossvalidated likelihoods
keyboard;
varargout={T};
