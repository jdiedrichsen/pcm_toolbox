function varargout=pcm_recipe_correlation
% This example uses a nonlinear correlation model to test for the correspondence of activity 
% patterns cross two conditions. 
% The data contains patterns for movements of 5 fingers for the left and
% right hand, respectively. 
% The data come from the 12 hemispheres studied for the finger movements
% of contralateral and ipsilateral hand (see Diedrichsen et al. 2013).
% Note that the strength of encoding for the contralateral fingers is much larger
% than for the ipsilateral fingers.
% The main variable of interest is the finger-to-finger correspondence between
% activity patterns for ipsi- and contralateral patterns. 
% We will test: 
%     a) 21 models, equally spaced from -1 to 1 correlation 
%     b) A flexible correlation model, that is functionally equivalent to
%     the feature model in pcm_recipe_feature
runEffect  = 'fixed';
nPart = 8; % Number of partitions 
nSubj = 10; % Number of subjects / simulations 
nVox  = 80; % Number of voxels

% --------------------------------------
% Generate data from a model with r =0.5 correlation 
Mtrue = pcm_buildCorrModel('type','nonlinear','withinCov','iid','numItems',2,'r','flexible','condEffect',0);
thetaTrue = [0 -0.5 1 ]' ; % true Theta's 
[G,dG] = pcm_calculateGnonlinCorr(thetaTrue,Mtrue)
keyboard;

Mtrue = pcm_buildCorrModel('type','nonlinear','withinCov','individual','numItems',5,'r',0.5);
thetaTrue = [-10 -10 -2 -2 -2 -2 -2 -3 -3 -3 -3 -3]' ; % true Theta's 
[Y,partVec,condVec] = pcm_generateData(Mtrue,thetaTrue,'numPart',nPart,'numVox',nVox,'numSim',nSubj,'signal',1,'noise',1);

% --------------------------------------
% Fixed models, equally spaced from -1 to 1 
nModel  = 21; 
r = linspace(-1,1,nModel); 
for i=1:nModel             
    M{i} = pcm_buildCorrModel('type','nonlinear','withinCov','individual','numItems',5,'r',r(i)); 
end

% --------------------------------------
% Flexible model
Mf{1} = pcm_buildCorrModel('type','nonlinear','withinCov','individual','numItems',5,'r','flexible'); 

% --------------------------------------
% Fit models from -1 to 1 to get average likelihood surface for the parameters 
[T,theta,G_hat] = pcm_fitModelIndivid(Y,M,partVec,condVec,'runEffect',runEffect);
            
% --------------------------------------
% Express all individual likelihood relative to the highest 
T.likelihood = bsxfun(@minus,T.likelihood,max(T.likelihood,[],2));
subplot(2,1,1); 
plot(r,T.likelihood); 
xlabel ('correlation'); 
ylabel ('log likelikhood'); 

% --------------------------------------
% Get full approximation of posterior probability distribution for each
% subject 
T.posterior = exp(T.likelihood); % Assuming uniform prior on [-1...1]
T.posterior = bsxfun(@rdivide,T.posterior,sum(T.posterior,2)); % Normalize to 1 
subplot(2,1,2); 
plot(r,T.posterior); 
xlabel ('correlation'); 
ylabel ('posterior'); 

% --------------------------------------
% Fit flexible model 
[T_flex,theta_flex,G_hat_flex] = pcm_fitModelIndivid(Y,Mf,partVec,condVec,'runEffect',runEffect);

% --------------------------------------
% Calculate the maximum a-posterior estimate of the correlation between pattern 
z = theta_flex{1}(Mf{1}.numGparams,:);  % Pick out correlation parameter 
r=(exp(2.*z)-1)./(exp(2.*z)+1);      % Take inverse Fisherz transform 
x=[r;r];
ylim = get(gca,'YLim'); 
y=repmat(ylim',1,numel(r)); 
line(x,y);                           
