function [T] = pcm_plotModelLikelihood(T,M,varargin)
% % function [T] = pcm_plotModelLikelihood(T,M,varargin)
% Plots the scaled log likelihood of the pattern component
% model(s) specified in M as a barplot. Errorbars reflect standard deviation (option
% to plot standard error of the mean- see OPTIONS). Noise ceiling is
% plotted as grey patch.
%
% Likelihoods are scaled to a null model (0) and the upper noiseceiling (1). If
% no model name in M matches 'null' or 'noiseceiling', the default is to
% scale likelihoods between the fits of the first (null) and last model fits.
%
% Upper bound of noise ceiling is the group fit of the noiseceiling model.
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
%     'Nnull':    Null model # (scale to 0). Default is 1
%     'Nceil':    Noiseceiling model # (scale to 1). Default is last model
%     'colors':   Cell of RGB cells for bar colors. Defaults included
%     'varfcn':   Error function used for errorbars. 'sem' or 'std' (default)
%     'mindx':    Models to plot as bars. Default is to plot all but those
%                 used for scaling likelihoods.
%     'Subj':     Specify single subjcet to plot fit of. NOTE: T will still
%                 be returned with scaled likelihoods for all subjects.
%
% OUTPUT:
%       T: Input structure returned with scaled likelihood fit fields:
%           'likelihood_norm' : Scaled crossvalidated fits
%           'likelihood_alln' : Scaled group (not crossvalidated) fits
%
% SArbuckle, 2016

% defaults
Nnull  = [];
Nceil  = [];
colors = {[.77 0 0],[0 0 .77],[.95 .6 0],[0 .8 0],[.8 0 .8],[.95 .95 0]};
varfcn = 'std';
mindx  = [];
Subj   = [];
vararginoptions(varargin,{'Nnull','Nceil','colors','varfcn','mindx','Subj'});
% - - - - - - -
% Check inputs
% - - - - - - -
if ~(isstruct(T)) || ~(isstruct(M))
    error('T or M are not pcm model structures. Check inputs.')
end
switch varfcn
    case 'sem'
        vfcn = @(x)nanstd(x)/sqrt(length(x));
    case 'std'
        vfcn = @nanstd;
    otherwise
        error('Unkown variance function')
end
% Locate models for likelihood scaling (first check type, then name)
if isempty(Nnull) % null model
    [~,Nnull] = find(strcmp({M.type},'null'));
    if (isempty(Nnull)) && (isfield(M,'name'))
        [~,Nnull] = find(strcmp({M.name},'null'));
        if isempty(Nnull)
            Nnull = 1;
        end
    else
        Nnull = 1;
    end
end;
if isempty(Nceil) % noise ceiling model
    [~,Nceil] = find(strcmp({M.type},'noiseceiling'));
    if (isempty(Nceil)) && (isfield(M,'name'))
        [~,Nceil] = find(strcmp({M.name},'noiseceiling'));
        if isempty(Nceil)
            Nceil = length(M);
        end
    else
        Nceil = length(M);
    end
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
    if (length(T.SN)>1) && (isfield(T,'likelihood_all')); % Plotting group level fit results
        T.likelihood_norm = bsxfun(@minus,T.likelihood,T.likelihood(:,Nnull)); % null fit is 0
        T.likelihood_alln = bsxfun(@minus,T.likelihood_all,T.likelihood(:,Nnull));
        noise_ceil        = T.likelihood_alln(:,Nceil);
        T.likelihood_norm = bsxfun(@rdivide,T.likelihood_norm,noise_ceil);     % ceiling fit is 1
        T.likelihood_alln = bsxfun(@rdivide,T.likelihood_alln,noise_ceil);
    elseif ~isfield(T,'likelihood_all') % Plotting fits for individual subject (only one subject in submitted structure)
        T.likelihood_norm = bsxfun(@minus,T.likelihood,T.likelihood(:,Nnull));
        noise_ceil        = T.likelihood_norm(:,Nceil);
        T.likelihood_norm = bsxfun(@rdivide,T.likelihood_norm,noise_ceil);
    end
end;

if (length(T.SN)>1) && (isfield(T,'likelihood_all') && isempty(Subj));
    lower_ceil = mean(T.likelihood_norm(:,Nceil));
elseif (~isfield(T,'likelihood_all') && (~isempty(Subj) || length(T.SN)==1))
    lower_ceil = 1;
else
    error('If plotting individual subject fits, please submit only one subject in T or specify ''Subj'' when calling pcm_plotModelLikelihood to plot models fit with pcm_fitModelIndivid. Note that you cannot plot single subject crossvalidated model fits.')
end;

% - - - - - - -
% Plot scaled fits
% - - - - - - -
figure; hold on;
i=1;
for m = mindx
    if isempty(Subj)
        Y(i) = mean(T.likelihood_norm(:,m));
        U(i) = vfcn(T.likelihood_norm(:,m));
    elseif Subj
        Y(i) = mean(T.likelihood_norm(Subj,m));
        U(i) = vfcn(T.likelihood_norm(Subj,m));
    end
    if isfield(M,'name') && (~isempty(M(m).name))
        labels{i} = M(m).name;
    else labels{i} = sprintf('Model %d',m);
    end
    bar(i,Y(i),'FaceColor',colors{i},'EdgeColor',colors{i});
    i=i+1;
end;

% Adjust figure settings
set(gca,'XTick',[1:i-1]);
set(gca,'XTickLabel',labels);
ylabel('Relative Likelihood');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) 90*(i-1) pos(4)]);
if (length(T.SN)>1 && isfield(T,'likelihood_all'))
    title('Crossval Group Model Fits');
elseif length(T.SN)==1
    title(sprintf('Subj %d Model Fits',T.SN));
elseif Subj
    title(sprintf('Subj %d Model Fits',Subj));
end;

% Plot grey noise ceiling
v = [0,lower_ceil; 0,1; i,1; i,lower_ceil];
f = [1:4];
patch('Vertices',v,'Faces',f,'EdgeColor',[0.9 0.9 0.9],...
    'FaceColor',[0.8 0.8 0.8],'FaceAlpha',.75);

% Plot errorbars (if group fits)
if isfield(T,'likelihood_all')
    errorbar([1:i-1],Y,zeros(1,length(Y)),U,'LineWidth',1.25,'Color',[0 0 0],...
        'LineStyle','none');
end
hold off;

