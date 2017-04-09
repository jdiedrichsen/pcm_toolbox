function varargout=pcm_recipe_correlation
% Example for a simple model with two factors:
% Factor A: 1-5 finger
% Factor B: left or right hand
% Example data come from the 12 hemispheres studied for the finger movements 
% of contralateral and ipsilateral hand (see Diedrichsen et al. 2013).
% Note that the strength of encoding for the right (contralateral) fingers is much larger
% than for the left (ipsilateral) fingers.
% The correlations are determined using a simple model that assumes the
% same variance for each 
%
load data_recipe_correlation.mat
% Treat the run effect as random or fixed?
% For run effect fixed, it does not matter if the
% If the run effect is random, it does matter, as the model needs
% to also explain some part of the common activity across
runEffect  = 'random';

% --------------------------------------
% Model1: Reduced model like in a chol-decomposition
M{1}.type       = 'feature';
M{1}.numGparams = 5;
M{1}.Ac(:,1,1)  = [ones(5,1);zeros(5,1)]; % Inidividual component: right hand (b)
M{1}.Ac(:,2,2)  = [zeros(5,1);ones(5,1)]; % Common component: left hand (c)
M{1}.Ac(:,3:7,3)  = [eye(5);zeros(5)]; % Inidividual component: right hand (b)
M{1}.Ac(:,3:7,4)  = [zeros(5);eye(5)]; % Common component: left hand (c)
M{1}.Ac(:,8:12,5)  = [zeros(5);eye(5)]; % Common component: right hand (d)
M{1}.name       = 'correlation';
M{1}.theta0=[1 1 1 0 0.4]';

% More complex feature model 
M{2}.type       = 'feature';
M{2}.numGparams = 17;
M{2}.Ac(:,1,1)  = [ones(5,1);zeros(5,1)]; % Inidividual component: right hand (b)
M{2}.Ac(:,2,2)  = [zeros(5,1);ones(5,1)]; % Common component: left hand (c)
for i=1:5
    A=zeros(5); 
    A(i,i)=1; 
    M{2}.Ac(:,3:7,2+i)  = [A;zeros(5)]; % Inidividual component: right hand (b)
    M{2}.Ac(:,3:7,7+i)  = [zeros(5);A]; % Common component: left hand (c)
    M{2}.Ac(:,8:12,12+i)  = [zeros(5);A]; % Common component: right hand (d)
    M{2}.name       = 'correlation2';
end; 
M{2}.theta0=ones(17,1);



% Fit the models
[D,theta,G_hat] = pcm_fitModelIndivid(Data,M,partVec,condVec,'runEffect',runEffect);

% Overall correlation 
T.Avar1    = theta{1}(:,3).^2;
T.Avar2    = theta{1}(:,4).^2+theta{1}(:,5).^2;
T.Acov     = theta{1}(:,3).*theta{1}(:,4);
T.Ar       = T.Acov./sqrt(T.Avar1.*T.Avar2);

T.Bvar1    = theta{2}(:,3:7).^2;
T.Bvar2    = theta{2}(:,8:12).^2+theta{2}(:,13:17).^2;
T.Bcov     = theta{2}(:,3:7).*theta{2}(:,8:12);
T.Br       = T.Bcov./sqrt(T.Bvar1.*T.Bvar2);
T.Bmr      = mean(T.Bcov,2)./sqrt(mean(T.Bvar1,2).*mean(T.Bvar2,2));


% --------------------------------------
% Calculate empricial correlation 
for p=1:12 
    Z=pcm_indicatorMatrix('identity',condVec{p}); 
    b = pinv(Z)*Data{p}; 
    x1 = b(1:5,:); 
    x2 = b(6:10,:); 
    x1 = bsxfun(@minus,x1,mean(x1)); 
    x2 = bsxfun(@minus,x2,mean(x2));
    R=corr([x1' x2']); 
    T.meanR(p,1) = mean(diag(R,5)); 
end; 
fprintf('mean correlation = %2.2f (+-%2.2f)\n', mean(T.meanR),std(T.meanR)/sqrt(12)); 

% --------------------------------------
% Crossvalidated correlation 
for p=1:12 
    Z=pcm_indicatorMatrix('identity',condVec{p}); 
    % Do the hard subtraction of each hand effect 
    X = indicatorMatrix('identity',partVec{p}*2+(condVec{p}>5)-1); 
    N=size(X,1); 
    R=eye(N)-X*pinv(X); 
    G_hat(:,:,p)=pcm_estGCrossval(R*Data{p},partVec{p},condVec{p}); 
    [T.Cv1(p,:),T.Cv2(p,:),T.Ccv(p,:),...
        T.Cr(p,:),T.Cmr(p,:)]=calcCorr(G_hat(:,:,p));
end; 
fprintf('mean correlation = %2.2f (+-%2.2f)\n', mean(T.Cr),std(T.Cr)/sqrt(12)); 


% Make Figure X of the paper 
subplot(1,4,1); 
imagesc_rectangle(mean(G_hat,3),'YDir','reverse');
drawline(5.5); 
drawline(5.5,'dir','horz');
subplot(1,4,2); 
histplot([T.Ar],'catX',[0.3:0.05:0.99]);
set(gca,'YLim',[0 12],'XLim',[0.4 1]);
subplot(1,4,3); 
histplot([T.Bmr],'catX',[0.3:0.05:0.99]);
set(gca,'YLim',[0 12],'XLim',[0.4 1]);
subplot(1,4,4); 
histplot([T.Cmr],'catX',[0.3:0.05:0.99]);
set(gca,'YLim',[0 12],'XLim',[0.4 1]); 

set(gcf,'PaperPosition',[2 2 9 1.8]);
wysiwyg; 

varargout={T}; 

function [v1,v2,cv,r,mr]=calcCorr(G); 
    d0=diag(G); 
    d5=diag(G,5); 
    v1 = d0(1:5)'; 
    v2 = d0(6:10)'; 
    cv = d5'; 
    r  = cv./sqrt(v1.*v2); 
    mr = mean(cv)/sqrt(mean(v1)*mean(v2)); 

