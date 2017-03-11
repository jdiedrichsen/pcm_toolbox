function [T] = pcm_plotModelLikelihood(T,M,varargin)
% function [T] = pcm_plotModelLikelihood(T,M,varargin)
%
% Plots the scaled log likelihood of the pattern component
% model(s) specified in M as a bar or scatterplot. Errorbars reflect 'sem' 
% (standard error of the mean; default) or 'std' (standard deviation).
%
% Two options for the noise ceiling:
% - Plot the entire noise ceiling (upper and lower- assuming the 
%    'upperceil' option is used; see Options) as a gray patch.
% - Plot only the lower noise ceiling as a black line (if no 'upperceil' is
%    supplied).
%
% Likelihoods are scaled such that zero is the likelihood of the null model. 
% If 'normalize' is set to 1, then likelihoods are normalized to set upper 
% ceiling to 1. For normalization, you must supply 'upperceil'; group fit of 
% the noiseceiling model from pcm_fitModelGroup (no crossvalidation).
% 
% The upper noise ceiling is the group fit of the noiseceiling model.
% The lower noise ceiling is the fit of the crossvalidated noiseceiling model.
%
% Two options to plot individual subject fits from pcm_fitModelIndivid:
%   - submitted T structure contains info for only one subject
%   - OR, specify subject of T structure to plot using 'subj' option
%
% Default is to plot group average likelihoods for each model fit.
%
% Returns structure T with likelihood_norm field (the scaled likelihoods).
% 
% NOTE: If field T.likelihood_norm already exists, function still
% normalizes values and overwrites this field. Exercise caution.
%
% ----------------------------- Inputs ------------------------------------
%
%   M:  Model structure. Accepts 'Name' field to label bars in plot.
%   T:  Output structure of model fits from pcm_fitModelCrossval
%
% ----------------------------- Options -----------------------------------
%
%   'Nnull':        Null model #. Default is 1. 
%   'Nceil':        Noiseceiling model #. If set to NaN, upper noise ceiling is 
%                    not plotted as a line/patch. If no # is given nor is it
%                    set to NaN, the first model labeled 'noiseceiling' will 
%                    be used for scaling and plotting.
%   'upperceil':    Upper noise ceiling SNx1 vector (usually group fit from noise ceiling model)
%   'normalize':    Normalizes all likelihoods to upper noise ceiling (if given) 
%   'colors':       Cell of RGB cells for bar colors. Defaults included
%   'varfcn':       Error function used for errorbars. 'sem' (default) or 'std'
%   'mindx':        Models to plot as bars. Default is to plot all but those
%                    used for scaling likelihoods.
%   'subj':         Specify single subjcet to plot fit of. NOTE: T will still
%                    be returned with scaled likelihoods for all subjects.
%   'style':        Specificy whether the plot will be 'bar' (default) or 'scatter'. 
%
% ----------------------------- Outputs -----------------------------------
%
%       T: Input structure returned with scaled likelihood fit fields:
%           'likelihood_norm' : Scaled crossvalidated likelihoods 
% 

% - - - - - - -
% Defaults
% - - - - - - -
Nnull     = 1;
Nceil     = [];
varfcn    = 'sem';
mindx     = [];
subj      = [];
upperceil = [];
normalize = 1; 
style     = 'bar';
ceilColor = [0.8 0.8 0.8];
colors    = {[.7 0 0],...      % red
             [0 0 .7],...      % blue
             [.9 .6 0],...     % orange
             [0 0.6 0.6],...   % cyan
             [0.5 0 0.5],...   % purple
             [0.2 0.6 0.2]};   % green
pcm_vararginoptions(varargin,{'Nnull','Nceil','upperceil','colors',...
    'varfcn','mindx','Subj','normalize','style','ceilColor'});
% - - - - - - -
% Check inputs
% - - - - - - -
if ~(isstruct(T))
    error('T is not pcm result structure. Check inputs.')
end

numSubj = size(T.likelihood,1);
% Check subj (determine if plotting group average or individual's data).
if (~isempty(subj) || length(T.SN)==1)
    % User didn't specify specific subject, but T only has one subject.
    if isempty(subj)  
        subj = T.SN;
    end
    % Take only data for specified subject.
    sf = @(x) x(subj,:);
    T  = structfun(sf,T,'UniformOutput',false); 
elseif (isempty(subj))
    % Plotting group avg.
    subj = 1:numSubj; 
end; 

% Check to see if we are plotting makes sense and inform user (if necessary).
if (isfield(M{1},'thetaIndiv') && (length(subj)>1 || isempty(subj)))
    warning('You are plotting the group average of fits from pcm_fitModelIndivid. This is not recommended as likelihoods are not crossvalidated.');
end

% Locate models for likelihood scaling (first check type, then name).
if isnan(Nceil) 
    Nceil     = [];
    normalize = 0;
else
    if isempty(Nceil) % noise ceiling model
        for m=1:length(M)
            if (strcmp(M{m}.type,'noiseceiling'));
                Nceil = m; 
            end
        end
    end
end;

% Determine models being plotted.
if isempty(mindx)
    mindx = 1:length(M);
    mindx([Nnull,Nceil]) = [];  % don't plot null and noise ceiling (if appropriate)
end

% Check # of colors.
if length(mindx) > length(colors)
    error('More models to be plotted than colors available')
end;

% Set the error function.
switch varfcn
    case 'sem'
        vfcn = @(x)nanstd(x)/sqrt(length(x));
    case 'std'
        vfcn = @nanstd;
    otherwise
        error('Unkown variance function for errobars')
end

% Determine plotting and errorbar functions.
switch style
    case 'bar'
        plot_fcn = @(x,y)  bar(x,y,'FaceColor',colors{x},...
                                   'EdgeColor',colors{x});
        ebar_fcn = @(x,ub,lb) line([x,x-0.1; x,x+0.1],[ub,ub; (ub+lb)/2,ub],...
                                   'LineWidth',2,...
                                   'Color',[0 0 0]);
    case 'scatter'
        plot_fcn = @(x,y) plot(x,y,'MarkerFaceColor',colors{x},...
                                   'MarkerEdgeColor',[0 0 0],...
                                   'MarkerSize',12,...
                                   'Marker','o');
        ebar_fcn = @(x,ub,lb) line([x,x-0.1,x-0.1; x,x+0.1,x+0.1],[ub,ub,lb; lb,ub,lb],...
                                   'LineWidth',2,...
                                   'Color',[0 0 0]);
end

% - - - - - - -
% Scale Likelihoods
% - - - - - - -
% Scale likelihoods between null and ceiling fits.
% Make all data relative to the the null model.
T.likelihood_norm = bsxfun(@minus,T.likelihood,T.likelihood(:,Nnull)); % null fit is 0
% Correct the upper noise ceiling for the new zero point.
if (~isempty(upperceil))
    upperceil = upperceil - T.likelihood(:,Nnull); 
end; 
if (normalize) 
    if (isempty(upperceil))
        error(sprintf('The normalize option is currently set to 1 but no noise ceiling is provided. Please either: \n   - provide upper noise ceiling for normalization (see ''upperceil'' option) \n   - set ''normalize'' option to 0 (see ''normalize'' option)')); 
    end; 
    T.likelihood_norm = bsxfun(@rdivide,T.likelihood_norm,upperceil);     % ceiling fit is 1
    upperceil         = 1; % Set the upper noise ceiling to 1. 
end; 

% - - - - - - -
% Plot noise ceilings
% - - - - - - -
figure('Color',[1 1 1]);
i = length(mindx)+1;
% Find lower noise ceiling.
lowerceil = T.likelihood_norm(:,Nceil); 
if (~isempty(lowerceil) && ~isempty(upperceil)) 
    % Draw a patch encompassing area between upper and lower noise ceilings
    v = [0,mean(lowerceil); 0,mean(upperceil); i,mean(upperceil); i,mean(lowerceil)];
    f = [1:4];
    patch('Vertices',v,'Faces',f,'EdgeColor',ceilColor,...
        'FaceColor',ceilColor,'FaceAlpha',.75);
    line([0;i],[mean(lowerceil);mean(lowerceil)],'LineWidth',1.5,'Color','k','LineStyle','-.'); 
elseif (~isempty(lowerceil) && isempty(upperceil)) 
    % Only draw the lower noise ceiling (upper not given when fcn called).
    line([0;i],[mean(lowerceil);mean(lowerceil)],'LineWidth',1.5,'Color','k','LineStyle','-.'); 
end; 
hold on;

% - - - - - - -
% Plot scaled fits
% - - - - - - -
i = 1;
for m = mindx
    % Model fit across subjs (if appropriate)
    Y = mean(T.likelihood_norm(subj,m));
    % Error of fits
    U = vfcn(T.likelihood_norm(subj,m));
    % Set xtick label
    if isfield(M{m},'name') && (~isempty(M{m}.name))
        labels{i} = M{m}.name;
    else
        labels{i} = sprintf('Model %d',m);
    end
    % Plot errorbars (if group fits)
    if length(subj)>1; ebar_fcn(i,Y+U,Y-U); end
    % Plot this model's fit
    plot_fcn(i,Y);
    % update ticker
    i = i + 1;
end;

% Adjust figure settings
set(gca,'XTick',[1:i-1]);
set(gca,'XTickLabel',labels);
set(gca,'XLim',[0.5 i-0.5]);
ylabel('Relative Likelihood');