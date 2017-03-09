function pcm_plotFittedG(G_hat,T,M,varargin)
% function pcm_fitModelGplot(G_hat,M)
%
% Plots the real and scaled model predicted second moment matrices (G) 
% of the pattern component model(s) in M. 
%
% Plotting using outputs from pcm_fitModelGroupCrossval:
% 
%   Model G_preds are scaled by the scaling factor in T. This is done because
%   the starting guess for theta params are scaled such that the first value
%   becomes 1. This reduces the number of parameters to maximize the fit by
%   1, reducing computation time. To be comparable to the real G, the
%   predicted matrices need to be "unscaled".
%
% Plotting using outputs from pcm_fitModelIndivid:
%
%   No scaling is applied for each subject. 
%   We also do not recommend plotting the group average of these Gs-
%   conceptually it does not make sense to do so.
%
% ----------------------------- Inputs ------------------------------------
%
%   'G_hat': Real/observed G. Takes mean if G is 3 dimensions.
%   'T':     Output structure from pcm_fitModelCrossval or
%             pcm_fitModelIndivid functions. 
%   'M':     Output structure from pcm_fitModelCrossval or
%             pcm_fitModelIndivid functions.
%
%   *** If imagescaling individual fits, ensure option 'Subj' 
%       is also passed to specify subject data to plot.
%
% ----------------------------- Options -----------------------------------
%
%   'mindx':    Specify which model predicted G-matries will be plotted.
%                Corresponds to model number(s) to plot. Default is to image
%                all but first and last models (both assumed to be used for
%                likelihood scaling purposes).
%   'style':    Specify if plotting 'imagesc' (default) or 'line' of model
%                Gs.
%   'cmap':     Specify colormap to use with imagesc (default = parula)
%   'subj':     Subject number of individual subject data to plot. Include
%                only if outputs from pcm_fitModelIndivid are passed.
%                Default plots group average.
%   'clims':    Forced color scaling applied to all matrices. Default is 
%                forced between 0 to the max value across all G-matrices.
%  'linewidth': Width of line plots if using style 'line'.
%    'colors':  Colors for lines if plotting using style 'line'. Includes
%                set of default colours (those in pcm_plotModelLikelihood). 
%

% - - - - - - -
% Defaults
% - - - - - - -
style     = 'imagesc';
cmap      = 'parula';
linewidth = 1.5;
subj      = [];
mindx     = [];
clims     = [];
colors    = {[0 0 0],...       % black (for observed G plot line)
             [.7 0 0],...      % red
             [0 0 .7],...      % blue
             [.9 .6 0],...     % orange
             [0 0.6 0.6],...   % cyan
             [0.5 0 0.5],...   % purple
             [0.2 0.6 0.2]};   % green
pcm_vararginoptions(varargin,{'cmap','subj','mindx','clims','style','linewidth','colors'});

% - - - - - - -
% Check inputs 
% - - - - - - -
% Check T and M.
if ~(isstruct(T)) || ~(isstruct(M{1}))
    error('T or M are not pcm model output structures. Check inputs.')
end

% Check subj (determine if plotting group average or individual's data).
if (~isempty(subj) || length(T.SN)==1)
    % User didn't specify specific subject, but T only has one subject.
    if isempty(subj)  
        subj = T.SN;
    end
    % Take only data for specified subject.
    sf = @(x) x(subj,:);
    T  = structfun(sf,T,'UniformOutput',false); 
end

% Check to see if we are plotting makes sense and inform user (if necessary).
if (isfield(M{1},'thetaIndiv') && (length(subj)>1 || isempty(subj)))
    warning('You are plotting the group average of fits from pcm_fitModelIndivid. This is not recommended as each subject will have different noise scaling.');
end

% Check G_hat.
if size(G_hat,3) > 1                % If we have G for >1 subject,...
    if isempty(subj)                % do we need to avg G (for group plot)?
        G_hat = mean(G_hat,3);      
    else
        G_hat = G_hat(:,:,subj);    % or take one subject's data?
    end
end;

% Determine which models to plot- don't plot null and noiseceiling models.
if isempty(mindx)
    mindx = [1:length(M)];
    for m = mindx
        if strcmp(M{m}.type,'noiseceiling');
            mindx(ismember(mindx,m)) = [];
        elseif strcmpi(M{m}.name,'null');
            mindx(ismember(mindx,m)) = [];
        end; 
    end; 
end
numG = length(mindx)+1;            
numC = size(G_hat,1);

% - - - - - - -
% Scale the G-matrices to the same scale. 
% (via subject scaling factors from fitting output)
% - - - - - - -
G{1} = G_hat;
max_G = max(max(G{1}));
for i = 1:length(mindx)
    m = mindx(i);
    if ~isempty(subj) % Plotting single subject data
        if isfield(M{m},'thetaCross')     % if using crossvalidated fitting output structure, need to calculate G
            G{i+1} = pcm_calculateG(M{m},M{m}.thetaCross(:,subj)).*T.scale(subj,m);
        elseif isfield(M{m},'thetaIndiv') % if using individual fitting output structure, G is already calcualted
            try
                G{i+1} = M{m}.G_pred(:,:,subj);
            catch
                G{i+1} = pcm_calculateG(M{m},M{m}.thetaIndiv(:,subj));
            end
        else
            error(sprintf(' Check the submitted structures for correct fields:\n\t- ''thetaCross'' (if submitting crossvalidated group structures)\n\t- ''thetaIndiv'' (if submitting individual (no cv) fitted structures)'))
        end
    else              % Plotting group Gs
        G_pred = zeros(numC,numC,length(unique(T.SN)));
        for s = T.SN'
            try
                G_pred(:,:,s) = pcm_calculateG(M{m},M{m}.thetaCross(:,s)).*T.scale(s,m);
            catch
                error(sprintf(' You are attempting to plot results from pcm_fitModelIndivid. \n You have submitted one subject''s G matrix but structure T has entries for multiple subjects. \n Please specify specific subject to plot using ''subj'' option.'));
            end
        end
        G{i+1} = mean(G_pred,3);
        clear G_pred
    end
    max_G(i+1) = max(max(G{i+1}));
end;

% - - - - - - -
% Plot scaled fitted G-matrices.
% - - - - - - -    
figure('Color',[1 1 1]);
switch style
    case 'imagesc'
        % Check color scaling
        if isempty(clims)
            clims = [0 max(max_G)];
        end
        % Imagesc G-matrices
        for i = 1:numG
            subplot(1,numG,i);
            imagesc(G{i},clims); 
            colormap(cmap);
            % Determine appropriate subplot title
            if i==1
                if ~isempty(subj)
                    title(sprintf('Subj %d  G',subj));
                else
                    title('Observed G');
                end
            else
                m = mindx(i-1);
                if isfield(M{m},'name') && (~isempty(M{m}.name))
                    title(M{m}.name);
                else
                    title(sprintf('Model %d',m));
                end
            end
            % reshape imagesc plot (clean some of surrounding whitespace)
            axis equal
            yx = get(gca,'YLim');
            set(gca,'YLim',[(yx(1)+numC+1) (yx(2)-numC-1)]);
            set(gca,'XTick',[],'YTick',[],'XColor',[1 1 1],'YColor',[1 1 1]);
        end;    
        % resize figure appropriately
        p = get(gcf,'Position');
        set(gcf,'Position',[p(1) p(2) p(3) p(4)/3]);
    case 'line'
        % Plot second moment matrices in 2 subplots: one for variances
        % (diag elements), the other for covariances (off-diag elements).
        subplot(1,3,1);     hold on; title('Variances');   xlabel('condition #');           ylabel('variance');
        subplot(1,3,[2:3]); hold on; title('Covariances'); xlabel('condition pair (a.u.)'); ylabel('covariance');
        leg = {};
        for i = 1:numG
            cv_indx = logical(tril(ones(size(G{i})),-1));
            subplot(1,3,1);     plot(diag(G{i}),    'LineWidth',linewidth,'Color',colors{i});
            subplot(1,3,[2:3]); plot(G{i}(cv_indx), 'LineWidth',linewidth,'Color',colors{i});
            % Determine appropriate subplot title.
            if i==1
                if ~isempty(subj)
                    leg{end+1} = sprintf('Subj %d  G',subj);
                else
                    leg{end+1} = 'Observed G';
                end
            else
                m = mindx(i-1);
                if isfield(M{m},'name') && (~isempty(M{m}.name))
                    leg{end+1} = M{m}.name;
                else
                    leg{end+1} = sprintf('Model %d',m);
                end
            end
        end
        % Apply the legend and scale both subplots to the same y-axis
        legend(leg,'location','northwest');
        legend boxoff;
        subplot(1,3,1);     yl(1,:) = get(gca,'YLim'); hold off;
        subplot(1,3,[2:3]); yl(2,:) = get(gca,'YLim'); hold off;
        subplot(1,3,1);     set(gca,'YLim',[0 max(max(yl))]); set(gca,'XLim',[1 numC]);
        subplot(1,3,[2:3]); set(gca,'YLim',[0 max(max(yl))]); set(gca,'XLim',[1 sum(cv_indx(:))]);
end
    
    
    
    
    
    
    
