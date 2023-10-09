function varargout=pcm_recipe_correlation
% This example uses pcm-correlation models to test for the correspondence of activity 
% patterns cross two conditions. 
% 
% The data contains patterns for movements of 5 fingers for the left and
% right hand, respectively. The data come from the 12 hemispheres studied for the finger movements
% of contralateral and ipsilateral hand (see Diedrichsen et al. 2013).
% Note that the strength of encoding for the contralateral fingers is much larger
% than for the ipsilateral fingers.
% The main variable of interest is the finger-to-finger correspondence between
% activity patterns for ipsi- and contralateral patterns. 
% We will test: 
%     a) 21 models, equally spaced from -1 to 1 correlation 
%     b) A flexible correlation model, that is functionally equivalent to
%     the feature model described in pcm_recipe_feature


load data_recipe_correlation.mat

runEffect  = 'fixed';

% --------------------------------------
% 1. Calculate empricial correlation
for p=1:12
    Z=pcm_indicatorMatrix('identity',condVec{p});
    b = pinv(Z)*Data{p};           % Estimate mean activities
    b(1:5,:)  = bsxfun(@minus,b(1:5,:) ,mean(b(1:5,:))); % Subtract mean per hand
    b(6:10,:) = bsxfun(@minus,b(6:10,:),mean(b(6:10,:)));
    G=cov(b');
    T.r_naive(p,1) = calcCorr(G);
end;

% --------------------------------------
% 2. Crossvalidated covariance martric and correlation
for p=1:12
    Z=pcm_indicatorMatrix('identity',condVec{p});
    % Subtract mean for each hand and run
    X = pcm_indicatorMatrix('identity',partVec{p}*2+(condVec{p}>5)-1);
    R=eye(size(X,1))-X*pinv(X);         % Residual forming matrix
    Gcv(:,:,p)=pcm_estGCrossval(R*Data{p},partVec{p},condVec{p});
    T.r_crossval(p,1)=calcCorr(pcm_makePD(Gcv(:,:,p)));
end;

% --------------------------------------
% Plot covarianc matrix and histograms 
rbins = linspace(0,1,21); 
subplot(4,2,1); 
imagesc(mean(Gcv,3)); 
subplot(4,2,3); 
histplot(T.r_naive,'catX',rbins);
subplot(4,2,4); 
histplot(T.r_crossval,'catX',rbins);

% --------------------------------------
% 3. Generate fixed models, equally spaced from 0 to 1 
nModel  = 51; 
r = linspace(0,1,nModel); 
for i=1:nModel             
    M{i} = pcm_buildCorrModel('type','nonlinear','withinCov','individual','numItems',5,'r',r(i)); 
end

% --------------------------------------
% 4. Generate flexible model
Mf{1} = pcm_buildCorrModel('type','nonlinear','withinCov','individual','numItems',5,'r','flexible'); 

% --------------------------------------
% 5. Fit models from -1 to 1 to get average likelihood surface for the parameters 
[T_fix,theta,G_hat] = pcm_fitModelIndivid(Data,M,partVec,condVec,'runEffect',runEffect);
            
% --------------------------------------
% Express all individual likelihood relative to the highest 
likeBaseline = mean(T_fix.likelihood,2); 
T_fix.likelihood = bsxfun(@minus,T_fix.likelihood,likeBaseline);
subplot(4,2,[5 6]); 
plot(r,T_fix.likelihood,'Color',[0.5 0.5 0.5]); 
xlabel ('correlation'); 
ylabel ('log likelikhood'); 

% --------------------------------------
% Get full approximation of posterior probability distribution for each
% subject 
T_fix.posterior = exp(T_fix.likelihood); % Assuming uniform prior on [-1...1]
T_fix.posterior = bsxfun(@rdivide,T_fix.posterior,sum(T_fix.posterior,2)); % Normalize to 1 
subplot(4,2,[7 8]); 
plot(r,T_fix.posterior,'Color',[0.5 0.5 0.5]); 
set(gca,'Box','off'); 
xlabel ('correlation'); 
ylabel ('posterior'); 

% --------------------------------------
% Fit flexible model 
[T_flex,theta_flex,G_hat_flex] = pcm_fitModelIndivid(Data,Mf,partVec,condVec,'runEffect',runEffect);

% --------------------------------------
% Calculate the maximum a-posterior estimate of the correlation between pattern 
z = theta_flex{1}(Mf{1}.numGparams,:);  % Pick out correlation parameter 
T.r_pcm=(exp(2.*z)-1)./(exp(2.*z)+1);      % Take inverse Fisherz transform 

subplot(4,2,[5 6]); 
hold on; 
% x=[r;r];
% ylim = get(gca,'YLim'); 
% y=repmat(ylim',1,numel(r)); 
% line(x,y);                  
plot(T.r_pcm,T_flex.likelihood-likeBaseline,'k.');
set(gca,'Box','off'); 
hold off; 
set(gcf,'PaperPosition',[2 2 5 10]); 
wysiwyg; 
keyboard; 

% --------------------------------------
% Get the correlation from a a covariance matrix
% By avagering across covariances and variances
function r=calcCorr(G);
d0=diag(G);
v1 = d0(1:5)';    % Variances contra
v2 = d0(6:10)';   % Variances ipsi
cv=diag(G,5);     % Covariance
r = mean(cv)/sqrt(mean(v1)*mean(v2));

