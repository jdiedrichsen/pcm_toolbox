function varargout=pcm_recipe_feature
% This example uses a PCM feature model to determine correspondence of activity 
% patterns cross two conditions. 
% The data contains patterns for movements of 5 fingers for the left and
% right hand, respectively. 
% The data come from the 12 hemispheres studied for the finger movements
% of contralateral and ipsilateral hand (see Diedrichsen et al. 2013).
% Note that the strength of encoding for the contralateral fingers is much larger
% than for the ipsilateral fingers.
% The main variable of interest is the finger-to-finger correspondence between
% activity patterns for ipsi- and contralateral patterns. 
% Specifically, we want to 
%     a) estimate the correlation 
%     b) test the models of 
%           r=0 : no relationship 
%           r=1 : ipsi patterns are pure reactivations of contralateral
%                   patterns 
%           0<r<1: Ipsilateral patterns are related, but have distinct own
%                  component 
load ../recipe_correlation/data_recipe_correlation.mat
runEffect  = 'fixed';
alg = 'NR'; 

% --------------------------------------
% Model1: Model with independent contra and ipsilateral patterns (zero
% correlation)
M{1}.type       = 'feature';
M{1}.numGparams = 12;
for i=1:5
    A=zeros(5);
    A(i,i)=1;
    M{1}.Ac(:,1:5 ,i)    = [A;zeros(5)];       % Contralateral finger patterns   (theta_a)
    M{1}.Ac(:,6:10,5+i) = [zeros(5);A];        % Unique Ipsilateral pattterns    (theta_c)
    M{1}.name       = 'null';
end;
M{1}.Ac(:,11  ,11)  = [ones(5,1);zeros(5,1)];  % Hand-specific component contra  (theta_d)
M{1}.Ac(:,12  ,12)  = [zeros(5,1);ones(5,1)];  % Hand-specific component ipsi    (theta_e)
M{1}.theta0=ones(12,1);                        % Starting values: could be closer, but converges anyways 
M{1}.fitAlgorithm = alg; 

% --------------------------------------
% Model2: Model with a flexible correlation for each finger 
M{2}.type       = 'feature';
M{2}.numGparams = 17;
for i=1:5
    A=zeros(5);
    A(i,i)=1;
    M{2}.Ac(:,1:5 ,i)    = [A;zeros(5)];       % Contralateral finger patterns   (theta_a)
    M{2}.Ac(:,1:5 ,5+i)  = [zeros(5);A];       % Mirrored Contralateralpatterns  (theta_b)
    M{2}.Ac(:,6:10,10+i) = [zeros(5);A];       % Unique Ipsilateral pattterns    (theta_c)
    M{2}.name       = 'flex';
end;
M{2}.Ac(:,11  ,16)  = [ones(5,1);zeros(5,1)];  % Hand-specific component contra  (theta_d)
M{2}.Ac(:,12  ,17)  = [zeros(5,1);ones(5,1)];  % Hand-specific component ipsi    (theta_e)
M{2}.theta0=ones(17,1);
M{2}.fitAlgorithm = alg; 

% --------------------------------------
% Model3: Model with a fixed r=1 correlation (ipsilateral = scaled version of contralateral pattern) 
M{3}.type       = 'feature';
M{3}.numGparams = 12;
for i=1:5
    A=zeros(5);
    A(i,i)=1;
    M{3}.Ac(:,1:5 ,i)    = [A;zeros(5)]; % Contralateral finger patterns   (theta_a)
    M{3}.Ac(:,1:5 ,5+i)  = [zeros(5);A]; % Mirrored Contralateralpatterns  (theta_b)
    M{3}.name       = 'one';
end;
M{3}.Ac(:,6,11)  = [ones(5,1);zeros(5,1)]; % Hand-specific component contra  (theta_d)
M{3}.Ac(:,7,12)  = [zeros(5,1);ones(5,1)]; % Hand-specific component ipsi    (theta_e)
M{3}.theta0=ones(12,1);
M{3}.fitAlgorithm = alg; 
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
% 2. Crossvalidated correlation
for p=1:12
    Z=pcm_indicatorMatrix('identity',condVec{p});
    % Subtract mean for each hand and run
    X = pcm_indicatorMatrix('identity',partVec{p}*2+(condVec{p}>5)-1);
    R=eye(size(X,1))-X*pinv(X);         % Residual forming matrix
    Gcv(:,:,p)=pcm_estGCrossval(R*Data{p},partVec{p},condVec{p});
    T.r_crossval(p,1)=calcCorr(pcm_makePD(Gcv(:,:,p)));
end;

% --------------------------------------
% 3. Fit model 2  and infer correlations from the parameters 
[D,theta,G_hat] = pcm_fitModelIndivid(Data,M,partVec,condVec,'runEffect',runEffect);

% Get the correlations 
var1       = (theta{2}(1:5,:).^2)';
var2       = (theta{2}(6:10,:).^2+theta{2}(11:15,:).^2)';
cov12      = (theta{2}(1:5,:).*theta{2}(6:10,:))';
T.r_model2 =  mean(cov12,2)./sqrt(mean(var1,2).*mean(var2,2));
            
% --------------------------------------
% Plot histgrams
edge=[0:0.1:1];
subplot(1,4,1);
hist(T.r_naive,edge);
set(gca,'XLim',[0 1.05],'YLim',[0 12]); 
title('Simple correlation');
subplot(1,4,2);
hist(T.r_crossval,edge);
set(gca,'XLim',[0 1.05],'YLim',[0 12]); 
title('Cross-validated correlations');
subplot(1,4,3);
hist(T.r_model2,edge);
set(gca,'XLim',[0 1.05],'YLim',[0 12]); 
title('PCM-model');

% --------------------------------------
% Now do crossvalidated model comparision: 
[Tgroup,thetaGr] = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1);
[Tcross,thetaCr] = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',thetaGr,'fitScale',1);
subplot(1,4,4);
pcm_plotModelLikelihood(Tcross,M,'normalize',0,'Nnull',1); 
varargout={D,Tgroup,Tcross};

% --------------------------------------
% Get the correlation from a a covariance matrix
% By avagering across covariances and variances
function r=calcCorr(G);
d0=diag(G);
v1 = d0(1:5)';    % Variances contra
v2 = d0(6:10)';   % Variances ipsi
cv=diag(G,5);     % Covariance
r = mean(cv)/sqrt(mean(v1)*mean(v2));
