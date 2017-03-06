function [T] = pcm_plotModelLikelihood(T,M,varargin)
% % function [T] = pcm_plotModelLikelihood(T,M,varargin)
% Plots the scaled log likelihood of the pattern component
% model(s) specified in M as a barplot. Errorbars reflect standard error of the mean
% Noise ceiling model is plotted as a gray area 
%
% Likelihoods are scaled to the likelihood of the null model (0). 
% If 'normalize' is set to 1, then likelihoods are normalized to 
% set upper ceiling to 1. 
% 
% The upper bound of noise ceiling is the group fit of the noiseceiling model.
% The lower bound is the fit of the crossvalidated noiseceiling model.
%
% Two options to plot individual subject fits from pcm_fitModelIndivid:
%   - submitted T structure contains info for only one subject
%   - OR, specify subject of T structure to plot using 'Subj' option
%
% INPUT:
%       M:  Model structure. Accepts 'Name' field to label bars in plot.
%       T:  Output structure of model fits from pcm_fitModelCrossval
%
% OPTIONS:
%     'Nnull':    # Null model. Default is 1
%     'Nceil':    # Noiseceiling model. Default is last model. If set to NaN, no noise ceiling is plotted 
%     'upperceil':Upper noise ceiling SNx1 vector (usually group fit from noise ceiling model)
%     'normalize':Normalizes all likelihoods to upper noise ceiling 
%     'colors':   Cell of RGB cells for bar colors. Defaults included
%     'varfcn':   Error function used for errorbars. 'sem' or 'std' (default)
%     'mindx':    Models to plot as bars. Default is to plot all but those
%                 used for scaling likelihoods.
%     'Subj':     Specify single subjcet to plot fit of. NOTE: T will still
%                 be returned with scaled likelihoods for all subjects.
%
% OUTPUT:
%       T: Input structure returned with scaled likelihood fit fields:
%           'likelihood_norm' : Scaled crossvalidated likelihoods 
%
% SArbuckle, 2016

% defaults
Nnull  = 1;
Nceil  = [];
colors = {[.77 0 0],[0 0 .77],[.95 .6 0],[0 .8 0],[.8 0 .8],[.95 .95 0]};
varfcn = 'sem';
mindx  = [];
Subj   = [];
upperceil =[];
normalize = 1; 
pcm_vararginoptions(varargin,{'Nnull','Nceil','upperceil','colors',...
    'varfcn','mindx','Subj','normalize'});
% - - - - - - -
% Check inputs
% - - - - - - -
if ~(isstruct(T))
    error('T is not pcm result structure. Check inputs.')
end
numSubj = size(T.likelihood,1); 
if (isempty(Subj)) 
    Subj = 1:numSubj; 
end; 
switch varfcn
    case 'sem'
        vfcn = @(x)nanstd(x)/sqrt(length(x));
    case 'std'
        vfcn = @nanstd;
    otherwise
        error('Unkown variance function')
end
% Locate models for likelihood scaling (first check type, then name)
if isempty(Nceil) % noise ceiling model
    for m=1:length(M)
        if (strcmp(M{m}.name,'noiseceiling'));
            Nceil = m; 
        end; 
    end; 
end;
% if set to NaN, don't use noiseceiling 
if isnan(Nceil) 
    Nceil =[]; 
end; 

% Determine models being plotted
if isempty(mindx)
    mindx = 1:length(M);
    mindx([Nnull,Nceil]) = [];
end

% Check colors
if length(mindx) > length(colors)
    error('More models to be plotted than colors available')
end;

% - - - - - - -
% Scale Likelihoods
% - - - - - - -
% Scale likelihoods between null and ceiling fits
if ~isfield(T,'likelihood_norm') % do we need to scale likelihoods?
    % Make all data relative to the the null model 
    T.likelihood_norm = bsxfun(@minus,T.likelihood,T.likelihood(:,Nnull)); % null fit is 0
    if (~isempty(upperceil))
        upperceil = upperceil - T.likelihood(:,Nnull); 
    end; 
    if (normalize) 
        if (isempty(upperceil))
            error('Please provide upper noise ceiling for normalization'); 
        end; 
        T.likelihood_norm = bsxfun(@rdivide,T.likelihood_norm,upperceil);     % ceiling fit is 1
        upperceil=1; % Set the upper noise ceiling to 1 
    end; 
end;

% - - - - - - -
% Plot scaled fits
% - - - - - - -
i=1;
for m = mindx
    Y(i) = mean(T.likelihood_norm(Subj,m));
    U(i) = vfcn(T.likelihood_norm(Subj,m));
    if isfield(M{m},'name') && (~isempty(M{m}.name))
        labels{i} = M{m}.name;
    else
        labels{i} = sprintf('Model %d',m);
    end
    bar(i,Y(i),'FaceColor',colors{i},'EdgeColor',colors{i});
    hold on; 
    i=i+1;
end;

% Adjust figure settings
set(gca,'XTick',[1:i-1]);
set(gca,'XTickLabel',labels);
ylabel('Relative Likelihood');

% Plot grey noise ceiling if given 
lowerceil = T.likelihood_norm(:,Nceil); 
if (~isempty(lowerceil) && isempty(upperceil)) 
    line([0;i],[mean(lowerceil);mean(lowerceil)]); 
elseif ~isempty(lowerceil)
    v = [0,mean(lowerceil); 0,mean(upperceil); i,mean(upperceil); i,mean(lowerceil)];
    f = [1:4];
    patch('Vertices',v,'Faces',f,'EdgeColor',[0.9 0.9 0.9],...
        'FaceColor',[0.8 0.8 0.8],'FaceAlpha',.75);
end; 

% Plot errorbars (if group fits)
errorbar([1:i-1],Y,zeros(1,length(Y)),U,'LineWidth',1.25,'Color',[0 0 0],...
        'LineStyle','none');
hold off;

