function varargout = pcm_searchlight_run(Searchlights,ImagesIn,ImagesOut,pcmfcn,varargin)
%% function varargout = pcm_searchlight_run(Searchlights,Images,pcmfcn,varargin)
% Run pcm on the surface for multiple subjects from raw time series image.
% 
% INPUTS:
%   Searchlights: Searchlight definitions for multiple subjects. 
% 
%   ImagesIn:       Input images for multiple subjects.
% 
%   ImagesOut:      Output images for pcmfcn.
% 
%   pcmfcn:         user-difined pcm model fitting function (e.g., )
%                   1st input argument should be;
%                   Y: a cell array of Ndata x Nvoxel data
% 
% OUTPUTS:
% 
% 
% 
% 
% 
% JD & AY 2016
% 
% Recommended spec: memory larger than ** GB?

% note: good to be compatible for both one-step (prewhitening+fitting) and 
%       two-step (prewhitening+save, then fitting) style


%% Input check (cell/structure/char)
% Options
global verbose;
verbose = 0;       % Control output to commandline 
params  = {};      % Extra parameter pased to pcmfcn
tmpDir  = {};      % Directory for tmp*.mat files (for each subject) 

vararginoptions(varargin,{'params','tmpDir','verbose'});

% Number of subjects
sn = numel(Searchlights);
if sn~=numel(ImagesIn);
    error('Size mismatch between Searchlight and ImagesIn!\n');
end

% Filetype for Searchlights (cell array of file names? structure array?)
% Check searchlights
if iscell(Searchlights)
    SLnames = Searchlights;
    clear Searchlights
    for s=1:sn
        Searchlights(s) = load(SLnames{s});
    end
end
Nnode1      = numel(Searchlights(1).LI); % assuming the same node size across sbjs
NsubjNode   = zeros(Nnode1,1);
for s=1:sn    
    Searchlights(s).Good = ~isnan(Searchlights(s).vORr); % non-empty node
    NsubjNode = NsubjNode + double(Searchlights(s).Good);
end

% Output
Nout = numel(ImagesOut);
if (~isempty(ImagesOut))&&(~iscell(ImagesOut))
    error('outfiles must be given in a cell array'); 
end;
for i=1:Nout
    [d f t m]       = spm_fileparts(ImagesOut{i});
    metricNames{i}  = fullfile(d, [f t]);
    metricColumn{i} = m(2:end);
    if (~strcmp(t,'.metric'))
        error('file type must me a metric file');
    end;
end
%=== currently only allow .metric file type for output?

if isempty(tmpDir)
   % get directory from Searchlight 
end




%% Main part (how can we speedup things?)

%-------------------------------------------------------------------------% 
% Split searchlight definitions?
% local_splitSearchlights(Searchlights,length(ImagesIn{2}),tmpDir);

%-------------------------------------------------------------------------% 
% Read volumes into memory space 
% (this is memory-consuming. Any smarter way?)
for s=1:sn
    Vol{s} = spm_vol(ImagesIn{s});
end

%-------------------------------------------------------------------------% 
% loop over searchlight center nodes (this loop can be parallelized)
Result = zeros(Nnode1,Nout)*NaN;
for i_node=1:100;%Nnode1
    if NsubjNode(i_node) >= 0.8*sn % 80% of subjects
        Y = cell(NsubjNode(i_node));
        
        c = 1;tic;        
        for s=1:sn
            LI = Searchlights(s).LI{i_node};            
            if ~isempty(LI)
                % Sample images (time series data) for LI{i_centernode}
                linVox    = unique(cat(2,LI')); % unique voxel linear indices
                [I,J,K]   = ind2sub(Vol{s}(1).dim,linVox);
                
                X = sparse(length(Vol{s}),double(max(linVox))); % create sparse matrix for speed and memory?
                for i_image=1:length(Vol{s}) % sample all candidate voxels at once
                    X(i_image,linVox) = spm_sample_vol(Vol{s}(i_image),double(I),double(J),double(K),0);
                end;
                
                % stack data (either TxP or NxP) into Y{subj}
                Y{c} = full(X(:,LI));
                c = c+1;                
            end
        end
        toc;
        
        % run pcm for i_node
        Result(i_node,:) = feval(pcmfcn,Y,params{:});        
    else
        if verbose==0
            fprintf('Insufficient number of subjects. Ignoreing node# %d\n',i_node);
        end
    end    
end

%-------------------------------------------------------------------------% 
% Write result into .metric file
local_saveMetricFile(Result,metricNames,metricColumn,Nnode1);

end

%% Local functions

%-------------------------------------------------------------------------% 
% Split searchlight into optimal size to save memory
%-------------------------------------------------------------------------% 
function local_splitSearchlights(Searchlights,NInputFiles,tmpDir)
sn = length(Searchlights);
if length(NInputFiles)==1
    NInputFiles = repmat(NInputFiles,sn,1);
end

% Calculate memory size (=Ncenternodes*Nvoxel*Ndatapoints*Nsubjects*64bit)
for s=1:sn
    nVox                    = length([Searchlight(s).LI{Searchlights(s).Good}])/...
                                sum(Searchlights(s).Good);
    totalMemory(s)          = numel(Searchlights(s).LI)*nVox*NInputFiles(s);
end

b=1;

% cd(tmpDir);
% while (~isDone)
%     % Find the candidate search lights
%     % This starts in the x-direction and finds the slice of search lights
%     % that are within that slice, and similar for y and z.
%     % then it moves the block tightly to this location
%     isCandidate = ~isCalc;
%     for i=1:3
%         IJKc1(i)    = min(Searchlight.voxmin(isCandidate,i));
%         IJKc2(i)    = IJKc1(i)+blockSize(i);     
%         isCandidate = isCandidate & Searchlight.voxmin(:,i)>=IJKc1(i) & ...
%                                     Searchlight.voxmax(:,i)<=IJKc2(i);
%     end;
%     IJKc1  = min(Searchlight.voxmin(isCandidate,:),[],1);
%     IJKc2  = IJKc1+blockSize;
%     isCandidate=~isCalc & Searchlight.voxmin(:,1)>=IJKc1(1) & Searchlight.voxmin(:,2)>=IJKc1(2) &  ...
%                           Searchlight.voxmin(:,3)>=IJKc1(3) & Searchlight.voxmax(:,1)<=IJKc2(1) &  ...
%                           Searchlight.voxmax(:,2)<=IJKc2(2) & Searchlight.voxmax(:,3)<=IJKc2(3);
%     j=find(isCandidate);
%     if (isempty(j));
%         break;
%     end;
%     
%     % Save the substructure as a tempory file
%     T.LI    = {Searchlight.LI{j}}';
%     T.j     = j;
%     T.voxel = Searchlight.voxel(j,:);
%     if verbose==0; 
%         fprintf('block %d Corner: %d %d %d length:%d  \n',b,IJKc1,length(j));
%     end
%     save(sprintf('temp_%2.2d.mat',b),'-struct','T');
%     isCalc(j) = true;
%     isDone    = all(isCalc);
%     b = b+1;
% end
% numBlocks = b-1;
end


%-------------------------------------------------------------------------% 
% Save result as .metric file
%-------------------------------------------------------------------------% 
function local_saveMetricFile(R,metricNames,metricColumn,Nnode)
%make the struct
M = struct();
for i=1:size(R,2)     % Loop over all rows in the MVA-result
    j=1;
    endwhile = numel(M);
    while j<= endwhile                  % Loop over all existing metric files
        if ~isfield(M, 'save_names')    % If there is no one yet, initialize the first
            M.save_names        = metricNames{1};
            M.column_name       = metricColumn(1);
            M.data              = R(:,1);
            M.num_rows          = Nnode;
            M.num_cols          = 1;
            minmax_colormapping = [-1 1];
            M.column_color_mapping  = minmax_colormapping;
            M.encoding              = {'BINARY'};
            M.index                 = (0:(Nnode-1))';
            j = j+1;
            
        elseif strcmp(M(j).save_names, metricNames{i})      % If a metric file with this name is already existing: add column
            M(j).column_name            = [M(j).column_name, metricColumn(i)];
            M(j).data                   = [M(j).data R(:,i)];
            M(j).num_cols               = M(j).num_cols+1;
            M(j).column_color_mapping   = repmat(minmax_colormapping,M(j).num_cols,1);
            j=numel(M)+1;
            
        elseif j==numel(M)          % If not, generate additional metric file
            M(j+1).save_names       = metricNames{i};
            M(j+1).column_name      = metricColumn(i);
            M(j+1).data             = R(:,i);
            M(j+1).num_rows         = Nnode;
            M(j+1).num_cols         = 1;
            minmax_colormapping     = [-1 1];
            M(j+1).column_color_mapping = minmax_colormapping;
            M(j+1).encoding         = {'BINARY'};
            M(j+1).index            = (0:(Nnode-1))';
            j=j+1;
        else
            j=j+1;
        end
    end
end
for i=1:numel(M)
    caret_savemetric(M(i).save_names, rmfield(M(i),'save_names') );   %flexible
end

end