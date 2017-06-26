function [T,M,Ti,Mi] = pcm_recipe_nonlinear
% Fit series of non-linear & noise-ceiling models (see 3.3 in PCM_toolbox).
% Non-linear models are defined by function that returns an estimated G
% matrix and the derivatives of G with respect to theta parameters.
% Derivitives used for 'minimize' function with employs gradient descent to 
% find theta params that yeild G with highest likelihood for that model.
%
%
% Data for this recipe are from {REF}. Dataset contains M1 activity patterns 
% for 4 subjects. The data structure is as follows:
%   'betaW': Multivariate noise normalized beta patterns for each condition
%   'tt'   : Trial type of corresponding betaW (used for G crossvalidation)
%   'run'  : Scanning run for corresponding betaW (used for G crossvalidation)
% There are 160 regressors: 8 runs of 20 conditions. The 20 conditions are
% pressing one finger of the right hand either 2,4,8, or 16 times in 6
% seconds (5 fingers * 4 speeds = 20 conds). 
%
%
% The ratio between paired finger pattern distances remains stable as activity 
% increases (i.e. as # of presses increase), but the distances increase. In 
% this recipe, we fit one model of how the finger patterns may scale as the
% number of presses increases. Model:
%   'Scaling':  patterns multiplicatively scaled by constant dependent
%               on the number of presses (G_2 = s*G_1, where G is the second 
%               moment martix at one pressing speed, and s is # presses scaling 
%               constant)
%   'Additive': patterns are shifted from baseline by a constant background
%               patten dependent on the speed. Therefore, the covariances
%               will remain approximately consistent, but the variance of
%               each condition in the second moment will increase with
%               pressing speed.
%   'Combination': patterns are scaled multiplicatively and shifted from
%                  baseline by some additive background pattern.
% 
% In addition to the above three nonlinear models, we fit two additional
% models:
%   'Null':     model that returns data where all distances are equal (i.e.
%               a bad model). This is used as the zero point when scaling 
%               likelihoods to an interpretable value
%   'NoiseCeiling': naive averaging model that overfits data. Uses the
%                   observed G as a fixed model, meaning it is the best possible 
%                   model. Likelihood of the non-crossvalidated group fit
%                   of this model is set to be 1 when scaling likelihoods.
%

load data_recipe_nonlinear.mat % loads struct I
runEffect = 'random';   
% The run effect is considered a random effect with zero mean this is 
% important, as we want to preserve the information of where the baseline is 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (1) Estimate crossvalidated G from acivity patterns
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for s=1:length(I) % Each row of I is one subject
    N                = length(I(s).run);      % number of condition regressors
    Y{s}             = I(s).betaW(1:N,:);     % condition-specific activity pattern
    conditionVec{s}  = I(s).tt;       
    partitionVec{s}  = I(s).run;
    G_hat(:,:,s)     = pcm_estGCrossval(Y{s},I(s).run,I(s).tt);
end;
     

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (2) Guess starting theta values
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PCM toolbox cannot estimate starting values for nonlinear scaling.
% Starting vals can be arbitrary but can drastically increase computation
% time.
scale_vals = [log(0.30);log(0.62);log(0.85)];
add_vals   = [log(0.2);  log(0.62); log(1)];

% Get starting values for the finger structure (Omega). Because we're interested
% in how these patterns scale, we use the (co-)variances from the portion
% fo the second moment matrix (G) that correspond to one pressing speed.
% These values will be used to estimate the (co-)variances of the fingers
% at the other three pressing soeeds. Importantly, results will be the same 
% no matter what of the 4 different pressing conditions are used to
% determine the starting values. Here we take the 15 G values for 16
% presses. We can further reduce the number of parameters to minimize by
% scaling the parameters such that the first param is equal to 1.
G_mean    = mean(G_hat,3);
[Fx0,~,~] = pcm_free_startingval(G_mean([16:20],[16:20])); % scales params such that G(1,1) is 1.


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (3) Specify Models
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Null model- all distances equal 
M{1}.type       = 'fixed'; 
M{1}.name       = 'Null';
M{1}.Gc         = eye(20);
M{1}.numGparams = 0; % totally fixed model- no free params
%   Use likelihood fit of this model as 0 scaling point in each subject

% Scaling model- distances multiplied by constant scaler dependent on number of presses
M{2}.type       = 'nonlinear'; 
M{2}.name       = 'Scaling';
M{2}.modelpred  = @ra_modelpred_scale;
M{2}.numGparams = 17; % 14 free theta params in Fx0 and 3 free scaling params
M{2}.theta0     = [Fx0;scale_vals]; 

% Additive independent model- adds independent pattern (NOT mean
% pattern) that scales with pressing speed
M{3}.type       = 'nonlinear'; 
M{3}.modelpred  = @ra_modelpred_add;
M{3}.numGparams = 17;  % 14 free theta params in Fx0 and 3 additive params
M{3}.theta0     = [Fx0;add_vals];
M{3}.name       = 'Additive';

% Additive independent + Scaling model combo
M{4}.type       = 'nonlinear'; 
M{4}.modelpred  = @ra_modelpred_addsc;
M{4}.numGparams = 20; % 14 free theta params in Fx0 and 3 free scaling params + 3 additive params
M{4}.theta0     = [Fx0;scale_vals;add_vals];
M{4}.name       = 'Combination';

% Naive averaging model- noise ceiling
M{5}.type       = 'noiseceiling';   
M{5}.name       = 'noiseceiling';
M{5}.numGparams = 0; % totally fixed model- no free params
M{5}.theta0     = [];
%   Use likelihood fit of this model as 1 scaling point in each subject

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (4) Fit Models and plot group lvl results
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% First fit at the group level without crossvalidation, where likelihood of 
% one subject's data is assessed under the group level data (which contains
% their data)
[Tgroup,theta_all,~] = pcm_fitModelGroup(Y,M,partitionVec,conditionVec,...
    'runEffect',runEffect,'isCheckDeriv',0);
% Second fit in a crossvalidated manner, where likelihood of one subject's
% data is assessed under the group level data of all other subjects.
[Tcross,~,G_pred] = pcm_fitModelGroupCrossval(Y,M,partitionVec,conditionVec,...
    'runEffect',runEffect,'isCheckDeriv',0,'groupFit',theta_all);
 
% Can scale and plot group likelihoods of model fits.
T = pcm_plotModelLikelihood(Tcross,M,'upperceil',Tgroup.likelihood(:,5),'style','bar');
% Returns T with subfields for scaled likelihoods (relative to null model (M1)
% and noise ceiling (M5). 
        
% We can also plot and compare the real/observed and estimate (co-)variance
% matrices.
for s = 1:size(G_pred{2},3) % for each subject
    G_scaling(:,:,s)  = G_pred{2}(:,:,s).*Tcross.scale(s,2);
    G_additive(:,:,s) = G_pred{3}(:,:,s).*Tcross.scale(s,2);
    G_combo(:,:,s)    = G_pred{4}(:,:,s).*Tcross.scale(s,2);
end
G_scaling     = mean(G_scaling,3);
G_additive    = mean(G_additive,3);
G_combo       = mean(G_combo,3);
% plot group crossval fitted G_scaling against mean of G_hat
figure(2);
subplot(1,4,1);
imagesc(G_mean);
title('group G-hat')
subplot(1,4,2);
imagesc(G_scaling);
title('group scaling G')
subplot(1,4,3);
imagesc(G_additive);
title('group additive G')
subplot(1,4,4);
imagesc(G_combo);
title('group combination G')

% We can also evalute how effective the scaling parameter estimtes using
% simple line plots. For example, we can take the diagonal of G_mean and
% G_scaling:
figure(3); 
hold on; 
plot(diag(G_mean),'LineWidth',2); 
plot(diag(G_scaling),'LineWidth',2); 
hold off;
legend({'observed','scaling'});
legend boxoff
xlabel('condition number');
xlim([1 20]);
ylabel('variance (a.u.)');
box off        
 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (5) Fit Model to single subjects and plot fits for one subj
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
[Tindivid,~,G_pred_individ] = pcm_fitModelIndivid(Y,M,partitionVec,conditionVec,...
    'runEffect',runEffect,'isCheckDeriv',0);
figure(4); 
pcm_plotModelLikelihood(Tindivid,Mi,'subj',4,'normalize',0,'plotceil',0);
% We don't plot a lower noise ceiling because the noiseceiling model in 
% Tindivid is NOT crossvalidated, and so it is not appropriate for a lower 
% noise ceiling.
