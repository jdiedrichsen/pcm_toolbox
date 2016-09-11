function pcm_plotFittedG(G_hat,T,M,varargin)
% % function pcm_fitModelGplot(G_hat,M)
% Imagescales the real and scaled model predicted second moment matrices (G) 
% of the pattern component model(s) in M. 
%
% Model G_preds are scaled by the scaling factor in T. This is done because
% the starting guess for theta params are scaled such that the first value
% becomes 1. This reduces the number of parameters to maximize the fit by
% 1, reducing computation time. To be comparable to the real G, the
% predicted matrices need to be "unscaled".
%
% INPUTS:
%   'G_hat': Real/observed G. Takes mean if G is 3 dimensions.
%   'T':     Output structure from pcm_fitModelCrossval or
%            pcm_fitModelIndivid functions. 
%   'M':     Output structure from pcm_fitModelCrossval or
%            pcm_fitModelIndivid functions.
%
%   *** If imagescaling individual fits, ensure option 'Subj' 
%       is also passed to specify subject data to plot.
%
% OPTIONS:
%   'cmap':     Specify colormap to use with imagesc (default = parula)
%   'Subj':     Subject number of individual subject data to plot. Include
%               only if outputs from pcm_fitModelIndivid are passed.
%   'mindx':    Specify which model predicted G-matries will be plotted.
%               Corresponds to model number(s) to plot. Default is to image
%               all but first and last models (both assumed to be used for
%               likelihood scaling purposes).
%   'clims':    Forced color scaling applied to all matrices. Default is 
%               forced between 0 to the max value across all G-matrices.
%
% SArbuckle, 2016

% defaults
cmap  = 'parula';
Subj  = [];
mindx = [];
clims = [];
vararginoptions(varargin,{'cmap','Subj','mindx','clims'});
% - - - - - - -
% Check inputs 
% - - - - - - -
% Check T and M
if ~(isstruct(T)) || ~(isstruct(M))
    error('T or M are not pcm model structures. Check inputs.')
end
% Check G_hat
if size(G_hat,3) > 1                % if G_hat for each subject
    if isempty(Subj)                % are we plotting group lvl G?
        G_real = mean(G_hat,3);     % if so, avg. G_hat across subjects
    else
        G_real = G_hat(:,:,Subj);   % otherwise take subject's G_hat
    end
end;
% Check Subj
if ~isempty(Subj)
    sf = @(x) x(Subj,:);
    T  = structfun(sf,T,'UniformOutput',false); % take only data for specified subject
end
% determine which models to imagescale
if isempty(mindx)
    mindx = [2:length(M)-1];
end
numG = length(mindx)+1;            
numC = size(G_real,1);


% - - - - - - -
% Scale the G-matrices to the same scale 
% - - - - - - -
G{1} = G_real;
max_G = max(max(G{1}));
for i = 1:length(mindx)
    m = mindx(i);
    if ~isempty(Subj)
        sc = mean(T.scale(:,m));
    else
        sc = mean(T.scale_all(:,m));    % get scaling param
    end
    G{i+1} = M(m).G_pred.*sc;       % apply scaling to model predicted G
    max_G(i+1) = max(max(G{i+1}));
end;


% - - - - - - -
% Imagescale G-matrices
% - - - - - - -    
% Check color scaling
if isempty(clims)
    clims = [0 max(max_G)];
end
% Now plot G-matrices
figure;
for i = 1:numG
    subplot(1,numG,i);
    imagesc(G{i},clims); 
    colormap(cmap);
    if i==1
        if ~isempty(Subj)
            title(sprintf('Subj %d  G',Subj));
        else
            title('Observed G');
        end
    else
        m = mindx(i-1);
        if isfield(M,'name') && (~isempty(M(m).name))
        	title(M(m).name);
        else
            title(sprintf('Model %d',m));
        end
    end
    axis equal
    yx = get(gca,'YLim');
    set(gca,'YLim',[(yx(1)+numC+(numC*.4)) (yx(2)-numC-(numC*.4))]);
    set(gca,'XTick',[],'YTick',[]);
end;    
% resize figure appropriately
p = get(gcf,'Position');
set(gcf,'Position',[p(1) p(2) p(3) p(4)/3]);
    
    
    
    
    
    
    
