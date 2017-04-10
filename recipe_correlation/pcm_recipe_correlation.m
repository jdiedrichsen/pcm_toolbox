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
runEffect  = 'fixed';

% --------------------------------------
% Model1: Model with one feature component for all fingers
M{1}.type       = 'feature';
M{1}.numGparams = 5;
M{1}.Ac(:,1:5 ,1)  = [eye(5);zeros(5)];      % Contralateral finger patterns   (a)
M{1}.Ac(:,1:5 ,2)  = [zeros(5);eye(5)];      % Mirrored Contralateralpatterns  (b)
M{1}.Ac(:,6:10,3)  = [zeros(5);eye(5)];      % Unique Ipsilateral pattterns    (c)
M{1}.Ac(:,11  ,4)  = [ones(5,1);zeros(5,1)]; % Hand-specific component contra  (d)
M{1}.Ac(:,12  ,5)  = [zeros(5,1);ones(5,1)]; % Hand-specific component ipsi    (e)
M{1}.name       = 'correlation';
M{1}.theta0=[1 1 0.5 0.1 0.1 ]';		% Starting values

% --------------------------------------
% Model1: Model with a feature component for ech fingers
M{2}.type       = 'feature';
M{2}.numGparams = 17;
for i=1:5
    A=zeros(5);
    A(i,i)=1;
    M{2}.Ac(:,1:5 ,i)    = [A;zeros(5)]; % Contralateral finger patterns   (a)
    M{2}.Ac(:,1:5 ,5+i)  = [zeros(5);A]; % Mirrored Contralateralpatterns  (b)
    M{2}.Ac(:,6:10,10+i) = [zeros(5);A]; % Unique Ipsilateral pattterns    (c)
    M{2}.name       = 'correlation2';
end;
M{2}.Ac(:,11  ,16)  = [ones(5,1);zeros(5,1)]; % Hand-specific component contra  (d)
M{2}.Ac(:,12  ,17)  = [zeros(5,1);ones(5,1)]; % Hand-specific component ipsi    (e)
M{2}.theta0=ones(17,1);

% --------------------------------------
% Fit the models and calculate the correlation
[D,theta,G_hat] = pcm_fitModelIndivid(Data,M,partVec,condVec,'runEffect',runEffect);

% Get the correlations from the parameters for Model1 
var1        = theta{1}(1,:).^2;
var2        = theta{1}(2,:).^2+theta{1}(3,:).^2;
cov12       = theta{1}(1,:).*theta{1}(2,:);
T.r_model1  = (cov12./sqrt(var1.*var2))';
% Get the correlations from the parameters for Model2 
var1       = (theta{2}(1:5,:).^2)';
var2       = (theta{2}(6:10,:).^2+theta{2}(11:15,:).^2)';
cov12      = (theta{2}(1:5,:).*theta{2}(6:10,:))';
T.r_model2 =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
            
% --------------------------------------
% Calculate empricial correlation
for p=1:12
    Z=pcm_indicatorMatrix('identity',condVec{p});
    b = pinv(Z)*Data{p};           % Estimate mean activities
    b(1:5,:)  = bsxfun(@minus,b(1:5,:) ,mean(b(1:5,:))); % Subtract mean per hand
    b(6:10,:) = bsxfun(@minus,b(6:10,:),mean(b(6:10,:)));
    G=cov(b');
    T.r_naive(p,1) = calcCorr(G);
end;

% --------------------------------------
% Crossvalidated correlation
for p=1:12
    Z=pcm_indicatorMatrix('identity',condVec{p});
    % Subtract mean for each hand and run
    X = indicatorMatrix('identity',partVec{p}*2+(condVec{p}>5)-1);
    R=eye(size(X,1))-X*pinv(X);         % Residual forming matrix
    Gcv(:,:,p)=pcm_estGCrossval(R*Data{p},partVec{p},condVec{p});
    T.r_crossval(p,1)=calcCorr(pcm_makePD(Gcv(:,:,p)));
end;

% --------------------------------------
% Plot histgrams
edge=[0:0.1:1];
subplot(1,4,1);
hist(T.r_naive,edge);
set(gca,'XLim',[0 1],'YLim',[0 12]); 
title('Simple correlation');
subplot(1,4,2);
hist(T.r_model1,edge);
set(gca,'XLim',[0 1],'YLim',[0 12]); 
title('Model 1');
subplot(1,4,3);
hist(T.r_model2,edge);
set(gca,'XLim',[0 1],'YLim',[0 12]); 
title('Model 2');
subplot(1,4,4);
hist(T.r_crossval,edge);
set(gca,'XLim',[0 1],'YLim',[0 12]); 
title('Crossvalidated');
varargout={T};

% --------------------------------------
% Get the correlation from a a covariance matrix
% By avagering across covariances and variances
function r=calcCorr(G);
d0=diag(G);
v1 = d0(1:5)';    % Variances contra
v2 = d0(6:10)';   % Variances ipsi
cv=diag(G,5);     % Covariance
r = mean(cv)/sqrt(mean(v1)*mean(v2));
